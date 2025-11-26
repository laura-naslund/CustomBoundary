# Reach Classification (Generic)
# Naslund.Laura@epa.gov, 20250116
# Ann.RoseberryLincoln@tetratech.com, 20240508
# Based on code written by Tom Barnum, USEPA, 20240229
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R v4.4.1

# User-specified variables:
#      state name (long form, not abbreviation)
#      qc_keep (a vector of which months' erom estimates to keep, along with mean annual erom estimate)
#      pct_var <- 60   # minimum percent of variation explained by components of PCA used for missing data imputation and hierarchical clustering
#      minCOMIDsCluster <- 0.1  # percent of COMIDs minimum per cluster, expressed as a fraction

# Input data required:
#      gadm41_USA_1 shapefile set of files containing state outline
#      StreamCat_clusterVars.csv containing StreamCat clustering variables

# Required input directory structure
#     - in.dir (contains csv file with potential StreamCat variables to include)
#     - in.dir/gadm41_USA_shp (contains USA state shapefiles)

# Creates output directory structure (starts at the parent directory of the
#     input directory), where brackets indicate the state abbreviation
#     - out.dir = ClusterOutput/[state]
#         * out.dir/Boundary
#         * out.dir/Ecoregions
#         * out.dir/HCPC
#         * out.dir/Histograms
#         * out.dir/NHDPlus
#         * out.dir/PCA
#         * out.dir/QC

#' @title Reach Classification (generic)
#'
#' @description Generates cluster identifiers for stream reaches within a selected state
#'
#' @details Identifies reach clusters based on stream flow and slope data obtained both from
#'          the nhdplusTools (using an API) and abiotic watershed summary data obtained using StreamCatTools.
#'          If needed, variables are Box-Cox transformed, missing values are imputed, and a PCA performed.
#'          Reaches are then clustered using kmeans on the principal components with a user specified number of clusters.
#'          A hierarchical clustering on the kmeans results is then performed and a final number
#'          of clusters is determined from a user-specified minimum proportion of stream reaches in each cluster.
#'          A final graphic with reach cluster assignments and PCA summary figures is then generated.
#'
#' Uses the libraries tidyverse, nhdplusTools, StreamCatTools, sf, moments, factoextra,
#' cowplot, FactoMineR, missMDA, tmap, viridis, ggrepel, readxl, tictoc
#'
#' @param state The desired state
#' @param qc_keep The EROM variables to be used in the PCA (QC_MM and QC_MA, where
#'                MM refers to the numeric monthly mean estimate and MA refers to
#'                the mean annual flow estimate, which is always required).
#' @param pct_var The minimum percentage variation that the PCA components must explain..
#' @param minCOMIDsCluster The minimum number of COMIDs required to be included in
#'                         the smallest cluster, to ensure that sufficient number
#'                         of sites are considered as comparators. This is expressed
#'                         as a proportion.
#' @param in.dir The path to the input data directory.
#'
#' @return Nothing is returned, but multiple files are written to the output directory.
#'
#' @keywords internal
#'
#' @export

boo.debug <- TRUE

clusterReachesCustom<- function(outputFolder, region_name, pct_var = 60, minCOMIDsCluster = 0.2, user_numclust = NULL, boo.debug = FALSE){

  if(boo.debug){
    outputFolder <- "C:/Users/lnaslund/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Test/Output"
    region_name <- "DEPied"
    pct_var <- 60
    minCOMIDsCluster <- 0.2
    user_numclust <- NULL
  }

  # STEP 1: Prep ----
  ## set vars ----
  tictoc::tic("Total")
  yyyymmdd <- format(lubridate::now(), "%Y%m%d")
  qc_keep <- c("QC_04", "QC_05", "QC_08","QC_MA")

  ## Declare functions ----
  `%>%` <- dplyr::`%>%`
  not_all_na <- function(x) {!all(is.na(x))}
  source("data-raw/clusterGraphic.R")
  source("data-raw/addClusterIDs.R")

  ## Load boundary from outputFolder ---
  load(file.path(outputFolder, "Boundary.rda"))

  ## Create output file paths ---
  if(dir.exists(file.path(outputFolder, "ClusterOutput")) == FALSE){dir.create(file.path(outputFolder, "ClusterOutput"))}

  if(dir.exists(file.path(outputFolder, "ClusterOutput", region_name)) == FALSE){dir.create(file.path(outputFolder, "ClusterOutput", region_name))}

  out.dir <- file.path(outputFolder, "ClusterOutput", region_name)

  out_folders <- c("Boundary", "HCPC", "Histograms", "NHDPlus", "PCA", "QC")

  for(i in 1:length(out_folders)){
    if(dir.exists(file.path(out.dir, out_folders[i]))==FALSE){
      dir.create(file.path(out.dir, out_folders[i]))
    }
  }

  ## Load/download required libraries
  libs <- c("tidyverse", "nhdplusTools", "StreamCatTools", "sf", "moments", "factoextra",
            "cowplot", "FactoMineR", "missMDA", "tmap", "viridis", "ggrepel", "readxl", "tictoc", "usethis", "caret")
  needed_libs <- setdiff(libs, .packages(all.available = TRUE))
  if(rlang::is_empty(needed_libs)==FALSE){
    install.packages(needed_libs)}
  lapply(libs, require, character.only = TRUE)
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # STEP 2: Get NHD+ data ----

  ## Get NHD+ data from API ----
  # Citation for NHDPlus data: McKay, L., Bondelid, T., Dewald, T., Johnston, J., Moore, R., and Rea, A., “NHDPlus Version 2: User Guide”, 2012 and U.S. Geological Survey, 2019, National Hydrography Dataset (ver. USGS National Hydrography Dataset Best Resolution (NHD) for Hydrologic Unit (HU) [specify number of HuC2s here - 2001 (published 20191002), accessed [date] at https://www.epa.gov/waterdata/get-nhdplus-national-hydrography-dataset-plus-data
  # Citation for nhdlusTools: Blodgett, D., Johnson, J.M., 2022, nhdplusTools: Tools for Accessing and Working with the NHDPlus, https://doi.org/10.5066/P97AS8JD
  # Desired NHD+ variables
   variables <- c("comid", tolower(qc_keep), "slope")
  #
  # # Check for existence of data
  if (file.exists(file.path(out.dir, "NHDPlus", paste0("NHD_", region_name, ".rda")))) {
    message("Previously saved NHDPlus data loaded")
    load(file.path(out.dir, "NHDPlus", paste0("NHD_", region_name, ".rda")))
  } else {
    message("Acquiring NHDPlus data")

    tictoc::tic("Get NHD+ data")
    reaches <- nhdplusTools::get_nhdplus(AOI = boundary) %>%
      dplyr::filter(ftype %in% c("Connector", "CanalDitch", "StreamRiver", "Drainageway", "ArtificialPath"))  %>%
      dplyr::select(all_of(variables))
    new.names <- c(toupper(variables), "geometry")
    colnames(reaches) <- paste(new.names)
    save(reaches, file = file.path(out.dir, "NHDPlus"
                                     , paste0("NHD_", region_name, ".rda")))
  }
  tictoc::toc(log = TRUE)

  file.copy(file.path(out.dir, "NHDPlus", paste0("NHD_", region_name, ".rda")), file.path(outputFolder, "Reaches.rda"))
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # STEP 3: Get StreamCat data ----
  ## Create a list of StreamCat variables & NHD variables used in the cluster analysis.
  ## Write a file of StreamCat variables used as stressors in the CASTool.
  tictoc::tic("Get StreamCat data")
  if (file.exists(file.path(out.dir, "NHDPlus", paste0("NHD_SC_", region_name, ".rda")))) {
    load(file.path(out.dir, "NHDPlus", paste0("NHD_SC_", region_name, ".rda")))
    message("Previously saved StreamCat data loaded")
  } else {
    message("Acquiring StreamCat data")

    ## Read clustering variables ----
    # Citation for StreamCat data: Hill, Ryan A., Marc H. Weber, Scott G. Leibowitz, Anthony R. Olsen, and Darren J. Thornbrugh, 2016. The Stream-Catchment (StreamCat) Dataset: A Database of Watershed Metrics for the Conterminous United States. Journal of the American Water Resources Association (JAWRA) 52:120-128. DOI: 10.1111/1752-1688.12372.
    # Citation for StreamCatTools:   Weber, Marc H, Hill, Ryan A., Brookes, Allen F. 2024, StreamCatTools: Tools to work with the StreamCat API within R and access the full suite of StreamCat and LakeCat metrics, https://usepa.github.io/StreamCatTools

    sc_ws_metrics <- read_csv(file.path("data-raw", "StreamCat_clusterVars.csv"))  %>%
      dplyr::filter(Type == "watershed") %>%
      pull(Variable)

    sc_ws_metrics_str <- paste(sc_ws_metrics, collapse = ",")

    desired_comids <- reaches %>% st_drop_geometry() %>% pull(COMID)
    n_desired_comids <- desired_comids %>% length()

    message("Reading in StreamCat variables for desired COMIDS")
    message("Requires ", ceiling(n_desired_comids/500), " requests")

    if(ceiling(n_desired_comids/500) == 0){
      message("n_stragglers = 0")
    } else {
      sc_comids <-NULL
      message("n_stragglers/500 > 0")
      for(q in 1:ceiling(n_desired_comids/500)){  # pulling in 500 COMIDs at a time to not overwhelm the server

        start_ind <- ((q-1)*500) + 1
        end_ind <- 500 * q
        print(paste0("getting ", start_ind, ":", end_ind))

        temp.comids <- desired_comids[start_ind:end_ind]

        temp.comids <- temp.comids[!is.na(temp.comids)] %>% paste(collapse = ",")

        tryCatch({temp_sc <- StreamCatTools::sc_get_data(metric = sc_ws_metrics_str,
                                                         aoi = 'watershed',
                                                         comid = temp.comids,
                                                         showAreaSqKm = TRUE) %>%
          dplyr::select(-catareasqkm, -catareasqkmrp100, -wsareasqkmrp100) %>%
          rename("COMID" = "comid")

        sc_comids <- sc_comids %>% bind_rows(temp_sc)
        }, error = function(msg){
          return(sc_comids)
        })
      }
    }

    WS.STATE <- sc_comids

    WS.STATE.SCvars <- dplyr::left_join(reaches, WS.STATE, by = "COMID")
    save(WS.STATE.SCvars, file = file.path(out.dir, "NHDPlus",
                                           paste0("NHD_SC_", region_name, ".rda")))

  }

  tictoc::toc(log = TRUE)
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # STEP 4: Data QC ----
  tictoc::tic("Perform QC of stream reach data")

  # check if there are any exact duplicates and remove
  n_dup <- WS.STATE.SCvars %>% group_by_all() %>% filter(n()>1) %>% nrow()
  if(n_dup == 0){
    message("No duplicated rows")
  }
  if(n_dup!= 0){
    message("Removing duplicate rows")
    WS.STATE.SCvars <- WS.STATE.SCvars %>% distinct()
  }

  # drop columns with all NAs
  na_cols <- WS.STATE.SCvars %>%
    st_drop_geometry() %>%
    select_if(function(x) all(is.na(x))) %>%
    names()

  if(length(na_cols) == 0){
    message("No completely empty columns")
  }
  if(length(na_cols) != 0){
    WS.STATE.SCvars <- WS.STATE.SCvars %>% dplyr::select(-all_of(na_cols))
  }

  # replace negative slopes with NAs
  n_neg_slope <- WS.STATE.SCvars %>% filter(SLOPE < 0) %>% nrow()
  WS.STATE.SCvars <- WS.STATE.SCvars %>% mutate(SLOPE = if_else(SLOPE < 0, NA, SLOPE))
  message(paste0("Replacing ", n_neg_slope, " negative slope values with NAs"))

  # remove COMIDs without any StreamCat data
  WS.STATE.FinalRaw <- WS.STATE.SCvars %>%
    st_drop_geometry() %>%
    rename("wsareasqkmws" = "wsareasqkm") %>%
    filter(!(if_all(ends_with("ws"), is.na)))

  n_na_row <- WS.STATE.SCvars %>% rename("wsareasqkmws" = "wsareasqkm") %>% filter(if_all(ends_with("ws"), is.na)) %>% nrow()
  message(paste0("Removing ", n_na_row, " rows without StreamCat data"))

  save(WS.STATE.FinalRaw, file = file.path(out.dir, "PCA",
                                           paste0("FinalRawData", region_name, ".rda")))

  tictoc::toc(log = TRUE)
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # STEP 5: Finalize dataset ----
  ## Get stats ----
  tictoc::tic("Get statistics & transform variables, if necessary")

  if (file.exists(file.path(out.dir, "Histograms", paste0(region_name, "_stats.tab")))) {
    WS.STATE.stats <- read.delim(file.path(out.dir, "Histograms", paste0(region_name, "_stats.tab")) , sep = "\t")
    message("Previously saved histogram stats")
  } else {

    WS.STATE.stats <- WS.STATE.FinalRaw %>%
      tidyr::pivot_longer(cols = !COMID, names_to = "Variable", values_to = "Value",
                          values_drop_na = FALSE) %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise(N = dplyr::n() - sum(is.na(Value))
                       , NumNAs = sum(is.na(Value))
                       , Min = min(Value, na.rm = TRUE)
                       , Max = max(Value, na.rm = TRUE)
                       , Mean = mean(Value, na.rm = TRUE)
                       , Median = median(Value, na.rm = TRUE)
                       , SD = sd(Value, na.rm = TRUE)
                       , Kurtosis = moments::kurtosis(Value, na.rm = TRUE)
                       , Skewness = moments::skewness(Value, na.rm = TRUE)
                       , SkewnessSq = Skewness * Skewness
                       , .groups = "drop_last")
    write.table(WS.STATE.stats
                , file.path(out.dir, "Histograms", paste0(region_name, "_stats.tab"))
                , sep = "\t", col.names = TRUE, row.names = FALSE, append = FALSE)
  }

  if (file.exists(file.path(out.dir, "QC", paste0(region_name, "_WsTransfData.tab"))) & file.exists(file.path(out.dir, "QC", paste0(region_name, "_lambdas.tab")))){
    WS.STATE.FinalTransf <- read.delim(file.path(out.dir, "QC", paste0(region_name, "_WsTransfData.tab")), sep = "\t")
    df.lambda <- read.delim(file.path(out.dir, "QC", paste0(region_name, "_lambdas.tab")), sep = "\t")
    message("Previously saved QC files")
  } else{

    ## Transform/scale variables ----
    ## Plot histograms before/after ----
    cols <- setdiff(colnames(WS.STATE.FinalRaw), "COMID")
    df.temp <- dplyr::select(WS.STATE.FinalRaw, COMID)
    df.temp_scaled <- dplyr::select(WS.STATE.FinalRaw, COMID)
    df.lambda <- data.frame(Variable = character(), Lambda = double())
    rm_vars <- c()

    for (c in seq_along(cols)) {
      col <- cols[c]
      print(col)
      fn <- paste0(col, ".png")
      fn2 <- paste0(col, "_transf.png")

      # plot histogram of data
      p <- ggplot2::ggplot(WS.STATE.FinalRaw, ggplot2::aes(x = .data[[col]])) +
        ggplot2::geom_histogram(bins = 500) +
        ggplot2::ggtitle(paste0("Histogram of ", col, " observations")) +
        ggplot2::xlab(col) +
        ggplot2::theme_bw()
      ggplot2::ggsave(file.path(out.dir, "Histograms", fn), p, width = 6, height = 4
                      , units = "in")

      sk_sq <- WS.STATE.stats %>%
        dplyr::filter(Variable == col) %>%
        dplyr::select(SkewnessSq)

      sd_dat <- WS.STATE.stats %>%
        dplyr::filter(Variable == col) %>%
        dplyr::select(SD)

      if (as.numeric(sd_dat) != 0) {

        if (grepl("pct", col)) {         # Do not transform PCT variables
          lambda <- NA_real_
          new_v <- WS.STATE.FinalRaw[[col]]
          df.temp <- df.temp %>% mutate(!!col := new_v)
          subtitle <- "Not transformed, centered and scaled"
        }

        else if (as.numeric(sk_sq)< 3){                                # Do not transform values with minimal skewness
          lambda <- NA_real_
          new_v <- WS.STATE.FinalRaw[[col]]
          df.temp <- df.temp %>% mutate(!!col := new_v)
          subtitle <- "Not transformed, centered and scaled"
        }

        else {     # Box-Cox transform highly skewed variables

          v_val <- WS.STATE.FinalRaw[[col]] + 1e-12

          if (min(v_val, na.rm = TRUE) < 0){
            v_val <- v_val - min(v_val, na.rm = TRUE)  + 1e-12 # add
          }

          bc <- caret::BoxCoxTrans(v_val, na.rm = TRUE)
          lambda <- bc$lambda

          # bc <- MASS::boxcox(lm(v_val ~ 1), plotit = FALSE) # for some reason, MASS::boxcox will not evaluate in a function
          # lambda <- bc$x[which.max(bc$y)]

          if (lambda != 0) {
            new_v <- ((v_val ^ lambda) - 1) / lambda
          } else {
            msg <- paste(col, " has lambda equal zero.")
            message(msg)
            new_v <- log10(v_val)
          }

          if (all(is.na(new_v))) {
            message("all new_v NA")
            new_v <- WS.STATE.FinalRaw[[col]]
            subtitle <- "Not transformed, centered and scaled"
          } else {
            subtitle <- paste0("lambda = ", round(lambda, 4), ", centered and scaled")
          }

          df.temp <- df.temp %>% mutate(!!col := new_v)
        }



        # Scale variables, this is just for visualization of the centered and scaled distribution.
        # The the unscaled and uncentered data are passed to the PCA functions which scale and center the data.
        new_v_scaled <- scale(new_v, center = TRUE, scale = TRUE)
        df.temp_scaled <- df.temp_scaled %>% mutate(!!col := new_v_scaled)

        p2 <- ggplot2::ggplot(df.temp_scaled, ggplot2::aes(x = .data[[col]])) +
          ggplot2::geom_histogram(bins = 500) +
          ggplot2::ggtitle(paste0("Histogram of ", col, " observations")) +
          ggplot2::labs(subtitle = subtitle) +
          ggplot2::xlab(col) +
          ggplot2::theme_bw()
        ggplot2::ggsave(file.path(out.dir, "Histograms", fn2), p2
                        , width = 6, height = 4, units = "in")

        df.lambda <- rbind(df.lambda, cbind(col, round(lambda, 4)))
      }

      if (as.numeric(sd_dat) == 0) {          # Do not include variables with zero variation
        rm_vars <- c(rm_vars, col)
      }

    }

    message(paste0("removed: ", paste(rm_vars, collapse = ", "), " for 0 variation"))

    WS.STATE.FinalTransf <- df.temp
    write.table(WS.STATE.FinalTransf, file.path(out.dir, "QC", paste0(region_name, "_WsTransfData.tab"))
                , sep = "\t", col.names = TRUE, row.names = FALSE, append = FALSE)

    write.table(df.lambda, file.path(out.dir, "QC", paste0(region_name, "_lambdas.tab"))
                , sep = "\t", col.names = TRUE, row.names = FALSE, append = FALSE)
  }

  if(file.exists(file.path(out.dir, "Histograms", paste0(region_name, "_TransfData_stats.tab")))){
    WS.STATE.transf.stats <- read.delim(file.path(out.dir, "Histograms", paste0(region_name, "_TransfData_stats.tab")), sep = "\t")
    message("Previous saved histogram transformed stats")
  } else{
    ## Get new stats ----
    WS.STATE.transf.stats <- WS.STATE.FinalTransf %>%
      tidyr::pivot_longer(cols = !COMID, names_to = "Variable"
                          , values_to = "Value", values_drop_na = FALSE) %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise(N = dplyr::n() - sum(is.na(Value))
                       , NumNAs = sum(is.na(Value))
                       , Min = min(Value, na.rm = TRUE)
                       , Max = max(Value, na.rm = TRUE)
                       , Mean = mean(Value, na.rm = TRUE)
                       , Median = median(Value, na.rm = TRUE)
                       , SD = sd(Value, na.rm = TRUE)
                       , Kurtosis = moments::kurtosis(Value, na.rm = TRUE)
                       , Skewness = moments::skewness(Value, na.rm = TRUE)
                       , SkewnessSq = Skewness * Skewness
                       , .groups = "drop_last")
    write.table(WS.STATE.transf.stats
                , file.path(out.dir, "Histograms", paste0(region_name, "_TransfData_stats.tab"))
                , sep = "\t", col.names = TRUE, row.names = FALSE, append = FALSE)
  }
  tictoc::toc(log = TRUE)
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # STEP 6: PCA ----

  if(file.exists(file.path(out.dir, "PCA", paste0(region_name, "_PCAresults.rda")))){
    load(file.path(out.dir, "PCA", paste0(region_name, "_PCAresults.rda")))
  } else {

    WS.STATE.FinalTransf.rownames <- WS.STATE.FinalTransf %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames(var = "COMID")

    WS.STATE.FinalTransf_complete <- WS.STATE.FinalTransf.rownames %>%
      filter(complete.cases(.))

    WS.STATE.PCA_complete <- FactoMineR::PCA(WS.STATE.FinalTransf_complete,
                                             scale.unit = TRUE,
                                             graph = FALSE)

    ncpGTpctvar <- min(which(factoextra::get_eigenvalue(WS.STATE.PCA_complete)[, 3] >= pct_var))

    ## Impute NA values ----
    tictoc::tic("Impute missing values")
    message("Number of missing values by variable")
    WS.STATE.FinalTransf %>% summarise(across(everything(), ~ sum(is.na(.)))) %>% pivot_longer(cols = everything()) %>% print(n=Inf) # print number of missing values imputed for each variable

    WS.STATE.impute.ALL.rownames <- missMDA::imputePCA(WS.STATE.FinalTransf.rownames,
                                                       ncp = ncpGTpctvar, scale = TRUE,
                                                       maxiter = 500)
    save(WS.STATE.impute.ALL.rownames,
         file = file.path(out.dir, "PCA", paste0(region_name, "_ImputedData.rda")))

    WS.STATE.imputedvals <- as.data.frame(WS.STATE.impute.ALL.rownames$completeObs)
    save(WS.STATE.imputedvals, file = file.path(out.dir, "PCA",
                                                paste0(region_name, "_Input4PCA.rda")))
    tictoc::toc(log = TRUE)

    ## Perform PCA ----
    tictoc::tic("Perform PCA")

    WS.STATE.PCA.rownames <- FactoMineR::PCA(WS.STATE.imputedvals,
                                             scale.unit = TRUE,
                                             ncp = ncpGTpctvar,
                                             graph = FALSE)

    ## Scree plot
    p <- factoextra::fviz_screeplot(WS.STATE.PCA.rownames, addlabels = TRUE, ylim = c(0, 50),
                                    linecolor = "black")
    ggplot2::ggsave(file = file.path(out.dir, "PCA"
                                     , paste0(region_name, "PCA_screeplot.png")),
                    p, dpi = 600, width = 5, height = 5, units = "in")
    ## Variable plot
    p <- factoextra::fviz_pca_var(WS.STATE.PCA.rownames, col.var = "contrib", add.labels = TRUE
                                  , gradient.cols = c("#85D54AFF", "#2D708EFF", "#440154FF")
                                  , labelsize = 3, repel = TRUE)
    ggplot2::ggsave(file = file.path(out.dir, "PCA", paste0(region_name, "_PCA_variables.png"))
                    , p, dpi = 600, width = 6, height = 6, units = "in")

    # Write outputs including eigenvalues, individuals, and variables
    write.table(WS.STATE.PCA.rownames$eig
                , file.path(out.dir, "PCA", paste0(region_name, "PCA_eigenvalues.tab"))
                , append = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
    write.table(WS.STATE.PCA.rownames$ind
                , file.path(out.dir, "PCA", paste0(region_name, "PCA_individuals.tab"))
                , append = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
    write.table(WS.STATE.PCA.rownames$var
                , file.path(out.dir, "PCA", paste0(region_name, "PCA_Variables.tab"))
                , append = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")

    save(WS.STATE.PCA.rownames, file = file.path(out.dir, "PCA", paste0(region_name, "_PCAresults.rda")))
  }
  tictoc::toc(log = TRUE)

  # STEP 7: HCPC ----
  tictoc::tic("Perform HCPC")
  time1 <- Sys.time()

  # Calculate the minimum number of COMIDs in a cluster
  numCOMIDs <- nrow(WS.STATE.PCA.rownames$ind$coord)
  minCOMIDs <- ceiling(numCOMIDs * minCOMIDsCluster)
  maxclust <- floor(1/minCOMIDsCluster)

  num_criteria <- data.frame(numclust = c(), meets_criteria = c())
  pick_list <- data.frame(region_name = c(), numclust = c(), fn = c())
  file_names <- c()

  ## Define region orientation
  boundary_bbox <- sf::st_bbox(boundary)
  bbox_length <- boundary_bbox["xmax"] - boundary_bbox["xmin"]
  bbox_width <- boundary_bbox["ymax"] - boundary_bbox["ymin"]

  if(bbox_length > bbox_width){
    orientation <- "horizontal"
  } else if(bbox_length < bbox_width){
    orientation <- "vertical"
  } else{
    orientation <- "horizontal"
  }

  # Perform clustering
  fc <- fastcluster::hclust.vector(as.matrix(WS.STATE.PCA.rownames$ind$coord), method = "ward")

  dated <- format(lubridate::now(), "%Y%m%d%H%M%S")

  if (is.null(user_numclust)) {

    message("User did not provide cluster number.")
    numclust <- seq(maxclust, 1, by = -1)

    for (c in seq_along(numclust)) {
      numtry <- numclust[c]
      print(numtry)

      tictoc::tic("Perform hierarchical clustering")

      temp_tree <- cutree(fc, k = numtry)
      temp_df <- data.frame(COMID = names(temp_tree), ClusterID = temp_tree %>% unname())

      temp_df_summary <- temp_df %>%
        dplyr::group_by(ClusterID) %>%
        dplyr::summarize(n = dplyr::n())

      fn_name <- paste0(region_name, "_ClusterGraphics_", dated, "_", numtry, ".png")
      fn = file.path(out.dir, fn_name)

      file_names <- c(file_names, fn)

      # doesn't generate the figure, ,I think it is because the region isn't in the landscape layout vector so you will need to add that as an option and provide as a parameter
      clusterGraphic(clusters = temp_df %>% mutate(COMID = as.integer(COMID)),
                     pca1 = WS.STATE.PCA.rownames,
                     flowlines = reaches,
                     sites = NA,
                     STATE.map = boundary,
                     map.title = region_name,
                     file.name = fn,
                     orient = orientation)


      if (min(temp_df_summary$n) > minCOMIDs) {
        msg <- paste0(numtry, " clusters meet the criteria specified.")
        message(msg)

        num_criteria <- num_criteria %>% bind_rows(data.frame(num_clust = numtry, meets_criteria = "yes"))

        write.csv(temp_df, file.path(outputFolder, "Clusters.csv"), row.names = FALSE)
        write.csv(temp_df, file.path(out.dir, "HCPC", paste0(region_name, "_Clusters_", numtry, ".csv")), row.names = FALSE)

      } else { # try then next lower number of clusters
        msg <- "Smallest cluster contains too few reaches."
        num_criteria <- num_criteria %>% bind_rows(data.frame(num_clust = numtry, meets_criteria = "no"))
        message(msg)

        write.csv(temp_df, file.path(out.dir, "HCPC", paste0(region_name, "_Clusters_", numtry, ".csv")), row.names = FALSE)
      }



      tictoc::toc(log = TRUE)
      } # end of for loop



    rm(msg, numtry)

  } else { # use the user-desired number of clusters

    numtry <- user_numclust
    message("User did provide cluster number.")


    temp_tree <- cutree(fc, k = numtry)
    temp_df <- data.frame(COMID = names(temp_tree), ClusterID = temp_tree %>% unname())


    fn_name <- paste0(region_name, "_ClusterGraphics_", dated, "_", numtry,"_fastclust", ".png")
    fn = file.path(out.dir, fn_name)


    clusterGraphic(clusters = temp_df %>% mutate(COMID = as.integer(COMID)),
                   pca1 = WS.STATE.PCA.rownames,
                   flowlines = reaches,
                   sites = NA,
                   STATE.map = boundary,
                   map.title = region_name,
                   file.name = fn)

    write.csv(temp_df, file.path(outputFolder, "Clusters.csv"), row.names = FALSE)
    write.csv(temp_df, file.path(out.dir, "HCPC", paste0(region_name, "_Clusters_", numtry, ".csv")), row.names = FALSE)

  } # End if user-desired number of clusters
  tictoc::toc(log = TRUE)

  num_criteria$file_names <- file_names

  default_file_name <- num_criteria %>% filter(meets_criteria == "yes") %>% arrange(desc(num_clust)) %>% slice(1) %>% pull(file_names)
  default_numclust <- num_criteria %>% filter(meets_criteria == "yes") %>% arrange(desc(num_clust)) %>% slice(1) %>% pull(num_clust)


  time2 <- Sys.time()

  time2 - time1
  beepr::beep(8)
  tictoc::toc(log = TRUE)

  ## print time spent ----
  print(region_name)
  tictoc::toc(log = TRUE)

  return(paste0("Default number of clusters: ", default_numclust))

}
