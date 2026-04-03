#' Get cluster file
#'
#' @param outputFolder file path containing output files
#' @param regionName region name
#' @param states vector of states in the analysis region
#'
#' @return writes a cluster files and graphics for a custom region to the outputFolder.
#' @return writes a folder with supplemental files like the variable stats file
#' @export
getClusterFile <- function(outputFolder, regionName, states, pct_var = 60, minCOMIDsCluster = 0.2, user_numclust = NULL){

  source("data-raw/clusterGraphic.R")

  if(dir.exists(file.path(outputFolder, paste0(regionName, "_ClusterSuppl"))) == FALSE){dir.create(file.path(outputFolder, paste0(regionName, "_ClusterSuppl")))}
  suppl_fp <- file.path(outputFolder, paste0(regionName, "_ClusterSuppl"))

  # Step 1: Get StreamCat data used for clustering  ----
  load(file.path(outputFolder, paste0(regionName, "_Reaches.rda")))

  sc_data <- NULL

  for(i in seq_len(length(states))){
    temp_state <- states[i]

    if(nchar(temp_state) == 2){
      temp_state <- state.name[which(state.abb == temp_state)]
    }

    message(paste0("Retrieving data from ", temp_state))

    temp_df <- CASToolHelperPckg::getSCClusterData(temp_state)
    sc_data <- sc_data |> dplyr::bind_rows(temp_df)
  }

  sc_data_region <- sc_data |> dplyr::filter(comid %in% reaches$COMID)

  # Step 2: Join to NHD data ----
  WS.STATE.SCvars <- reaches |> dplyr::full_join(sc_data_region, by = c("COMID" = "comid")) |>
    dplyr::rename("wsareasqkm"="TotDASqKM") |>
    dplyr::select(-FTYPE)

  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # Step 3: Data QC ----
  ## Check if there are any exact duplicates and remove ----
  n_dup <- WS.STATE.SCvars |>
    dplyr::group_by_all() |>
    dplyr::filter(dplyr::n()>1) |>
    nrow()

  if(n_dup == 0){
  }
  if(n_dup!= 0){
    message("Removing duplicate rows")
    WS.STATE.SCvars <- WS.STATE.SCvars |> dplyr::distinct()
  }

  ## Drop columns with all NAs ----
  na_cols <- WS.STATE.SCvars |>
    sf::st_drop_geometry() |>
    dplyr::select_if(function(x) all(is.na(x))) |>
    names()

  if(length(na_cols) == 0){
    message("No completely empty columns")
  }
  if(length(na_cols) != 0){
    WS.STATE.SCvars <- WS.STATE.SCvars |> dplyr::select(-all_of(na_cols))
  }

  ## Replace negative slopes with NAs ----
  n_neg_slope <- WS.STATE.SCvars |> dplyr::filter(SLOPE < 0) |> nrow()
  WS.STATE.SCvars <- WS.STATE.SCvars |> dplyr::mutate(SLOPE = dplyr::if_else(SLOPE < 0, NA, SLOPE))
  message(paste0("Replacing ", n_neg_slope, " negative slope values with NAs"))

  ## Remove COMIDs without any StreamCat data ----
  ## slightly different syntax from prior cluster code because WS area is from NHD not streamcat
  WS.STATE.FinalRaw <- WS.STATE.SCvars |>
    sf::st_drop_geometry() |>
    dplyr::filter(!(dplyr::if_all(dplyr::ends_with("ws"), is.na))) |>
    dplyr::rename("wsareasqkmws" = "wsareasqkm")

  n_na_row <- nrow(WS.STATE.SCvars) - nrow(WS.STATE.FinalRaw)
  message(paste0("Removing ", n_na_row, " rows without StreamCat data"))

  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # Step 4: Finalize dataset ----
  ## Get stats ----
  WS.STATE.stats <- WS.STATE.FinalRaw |>
    tidyr::pivot_longer(cols = !COMID, names_to = "Variable", values_to = "Value",
                        values_drop_na = FALSE) |>
    dplyr::group_by(Variable) |>
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

  write.csv(WS.STATE.stats, file.path(suppl_fp, paste0(regionName, "_stats.csv")), row.names = FALSE)

  ## Transform/scale variables ----
  ## Plot histograms before/after ----
  if(dir.exists(file.path(suppl_fp, "Histograms")) == FALSE){dir.create(file.path(suppl_fp, "Histograms"))}

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
    ggplot2::ggsave(file.path(suppl_fp, "Histograms", fn), p, width = 6, height = 4
                    , units = "in")

    sk_sq <- WS.STATE.stats |>
      dplyr::filter(Variable == col) |>
      dplyr::select(SkewnessSq)

    sd_dat <- WS.STATE.stats |>
      dplyr::filter(Variable == col) |>
      dplyr::select(SD)

    if (as.numeric(sd_dat) != 0) {

      if (grepl("pct", col)) {         # Do not transform PCT variables
        lambda <- NA_real_
        new_v <- WS.STATE.FinalRaw[[col]]
        df.temp <- df.temp |> dplyr::mutate(!!col := new_v)
        subtitle <- "Not transformed, centered and scaled"
      }

      else if (as.numeric(sk_sq)< 3){                                # Do not transform values with minimal skewness
        lambda <- NA_real_
        new_v <- WS.STATE.FinalRaw[[col]]
        df.temp <- df.temp |> dplyr::mutate(!!col := new_v)
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

        df.temp <- df.temp |> dplyr::mutate(!!col := new_v)
      }



      # Scale variables, this is just for visualization of the centered and scaled distribution.
      # The the unscaled and uncentered data are passed to the PCA functions which scale and center the data.
      new_v_scaled <- scale(new_v, center = TRUE, scale = TRUE)
      df.temp_scaled <- df.temp_scaled |> dplyr::mutate(!!col := new_v_scaled)

      p2 <- ggplot2::ggplot(df.temp_scaled, ggplot2::aes(x = .data[[col]])) +
        ggplot2::geom_histogram(bins = 500) +
        ggplot2::ggtitle(paste0("Histogram of ", col, " observations")) +
        ggplot2::labs(subtitle = subtitle) +
        ggplot2::xlab(col) +
        ggplot2::theme_bw()
      ggplot2::ggsave(file.path(suppl_fp, "Histograms", fn2), p2
                      , width = 6, height = 4, units = "in")

      df.lambda <- rbind(df.lambda, cbind(col, round(lambda, 4)))
    }

    if (as.numeric(sd_dat) == 0) {          # Do not include variables with zero variation
      rm_vars <- c(rm_vars, col)
    }

  }

  message(paste0("removed: ", paste(rm_vars, collapse = ", "), " for 0 variation"))

  WS.STATE.FinalTransf <- df.temp
  write.csv(df.lambda |> dplyr::rename("Variable" = "col", "Lambda" = "V2"), file.path(suppl_fp, paste0(regionName, "_lambdas.csv")), row.names = FALSE)

  ## Get new stats ----
  WS.STATE.transf.stats <- WS.STATE.FinalTransf |>
    tidyr::pivot_longer(cols = !COMID, names_to = "Variable"
                        , values_to = "Value", values_drop_na = FALSE) |>
    dplyr::group_by(Variable) |>
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

    write.csv(WS.STATE.transf.stats, file.path(suppl_fp, paste0(regionName, "_stats_transf.csv")), row.names = FALSE)

  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # Step 5: PCA ----

    WS.STATE.FinalTransf.rownames <- WS.STATE.FinalTransf |>
      tibble::remove_rownames() |>
      tibble::column_to_rownames(var = "COMID")

    WS.STATE.FinalTransf_complete <- WS.STATE.FinalTransf.rownames |>
      tidyr::drop_na()

    WS.STATE.PCA_complete <- FactoMineR::PCA(WS.STATE.FinalTransf_complete,
                                             scale.unit = TRUE,
                                             graph = FALSE)

    ncpGTpctvar <- min(which(factoextra::get_eigenvalue(WS.STATE.PCA_complete)[, 3] >= pct_var))

    ## Impute NA values ----
    message("Number of missing values by variable")
    WS.STATE.FinalTransf |>
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.)))) |>
      tidyr::pivot_longer(cols = dplyr::everything()) |>
      print(n=Inf) # print number of missing values imputed for each variable

    WS.STATE.impute.ALL.rownames <- missMDA::imputePCA(WS.STATE.FinalTransf.rownames,
                                                       ncp = ncpGTpctvar, scale = TRUE,
                                                       maxiter = 500)

    WS.STATE.imputedvals <- as.data.frame(WS.STATE.impute.ALL.rownames$completeObs)

    ## Perform PCA ----
    WS.STATE.PCA.rownames <- FactoMineR::PCA(WS.STATE.imputedvals,
                                             scale.unit = TRUE,
                                             ncp = ncpGTpctvar,
                                             graph = FALSE)

    ## Scree plot
    p <- factoextra::fviz_screeplot(WS.STATE.PCA.rownames, addlabels = TRUE, ylim = c(0, 50),
                                    linecolor = "black")
    ggplot2::ggsave(file = file.path(suppl_fp
                                     , paste0(regionName, "PCA_screeplot.png")),
                    p, dpi = 600, width = 5, height = 5, units = "in")
    ## Variable plot
    p <- factoextra::fviz_pca_var(WS.STATE.PCA.rownames, col.var = "contrib", add.labels = TRUE
                                  , gradient.cols = c("#85D54AFF", "#2D708EFF", "#440154FF")
                                  , labelsize = 3, repel = TRUE)
    ggplot2::ggsave(file = file.path(suppl_fp, paste0(regionName, "_PCA_variables.png"))
                    , p, dpi = 600, width = 6, height = 6, units = "in")

    # Write eigenvalues
    write.csv(WS.STATE.PCA.rownames$eig, file.path(suppl_fp, "PCA_eigenvalues.csv"))

    # Write outputs including eigenvalues, individuals, and variables
    # write.table(WS.STATE.PCA.rownames$eig
    #             , file.path(outputFolder, "PCA", paste0(regionName, "PCA_eigenvalues.tab"))
    #             , append = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
    # write.table(WS.STATE.PCA.rownames$ind
    #             , file.path(outputFolder, "PCA", paste0(regionName, "PCA_individuals.tab"))
    #             , append = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
    # write.table(WS.STATE.PCA.rownames$var
    #             , file.path(outputFolder, "PCA", paste0(regionName, "PCA_Variables.tab"))
    #             , append = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")

    # save(WS.STATE.PCA.rownames, file = file.path(outputFolder, "PCA", paste0(regionName, "_PCAresults.rda")))


  # STEP 7: HCPC ----
  tictoc::tic("Perform HCPC")

  # Calculate the minimum number of COMIDs in a cluster
  numCOMIDs <- nrow(WS.STATE.PCA.rownames$ind$coord)
  minCOMIDs <- ceiling(numCOMIDs * minCOMIDsCluster)
  maxclust <- floor(1/minCOMIDsCluster)

  num_criteria <- data.frame(numclust = c(), meets_criteria = c())

  ## Define region orientation
  load(file.path(outputFolder, paste0(regionName, "_Boundary.rda")))

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

  if (is.null(user_numclust)) {

    message("User did not provide cluster number.")
    numclust <- seq(maxclust, 1, by = -1)

    for (c in seq_along(numclust)) {
      numtry <- numclust[c]
      print(numtry)

      tictoc::tic("Perform hierarchical clustering")

      temp_tree <- cutree(fc, k = numtry)
      temp_df <- data.frame(COMID = names(temp_tree), ClusterID = temp_tree |> unname())

      temp_df_summary <- temp_df |>
        dplyr::group_by(ClusterID) |>
        dplyr::summarize(n = dplyr::n())

      fn_name <- paste0(regionName, "_ClusterGraphics_", numtry, ".png")
      fn = file.path(outputFolder, fn_name)

      # doesn't generate the figure, I think it is because the region isn't in the landscape layout vector so you will need to add that as an option and provide as a parameter
      suppressMessages(suppressWarnings(clusterGraphic(clusters = temp_df |> dplyr::mutate(COMID = as.integer(COMID)),
                                                       pca1 = WS.STATE.PCA.rownames,
                                                       flowlines = reaches,
                                                       sites = NA,
                                                       STATE.map = boundary,
                                                       map.title = regionName,
                                                       file.name = fn,
                                                       orient = orientation)))


      if (min(temp_df_summary$n) > minCOMIDs) {
        msg <- paste0(numtry, " clusters meet the criteria specified.")
        message(msg)

        num_criteria <- num_criteria |> dplyr::bind_rows(data.frame(num_clust = numtry, meets_criteria = "yes"))

        write.csv(temp_df, file.path(outputFolder, paste0(regionName, "_Clusters_", numtry, ".csv")), row.names = FALSE)

      } else { # try then next lower number of clusters
        msg <- "Smallest cluster contains too few reaches."
        num_criteria <- num_criteria |> dplyr::bind_rows(data.frame(num_clust = numtry, meets_criteria = "no"))
        message(msg)

        write.csv(temp_df, file.path(outputFolder, paste0(regionName, "_Clusters_", numtry, ".csv")), row.names = FALSE)
      }



      tictoc::toc(log = TRUE)
    } # end of for loop

    default_numclust <- num_criteria |>
      dplyr::filter(meets_criteria == "yes") |>
      dplyr::arrange(desc(num_clust)) |>
      dplyr::slice(1) |>
      dplyr::pull(num_clust)

    message(paste0("Default number of clusters: ", default_numclust))

    file.copy(file.path(outputFolder, paste0(regionName, "_ClusterGraphics_", default_numclust, ".png")), file.path(outputFolder, paste0(regionName, "_ClusterGraphic_default.png")))
    file.copy(file.path(outputFolder, paste0(regionName, "_Clusters_", default_numclust, ".csv")), file.path(outputFolder, paste0(regionName, "_Clusters_default.csv")))

    rm(msg, numtry)

  } else { # use the user-desired number of clusters

    numtry <- user_numclust
    message("User did provide cluster number.")


    temp_tree <- cutree(fc, k = numtry)
    temp_df <- data.frame(COMID = names(temp_tree), ClusterID = temp_tree |> unname())


    fn_name <- paste0(regionName, "_ClusterGraphics_", numtry, ".png")
    fn = file.path(outputFolder, fn_name)


    suppressMessages(suppressWarnings(clusterGraphic(clusters = temp_df |> dplyr::mutate(COMID = as.integer(COMID)),
                                                     pca1 = WS.STATE.PCA.rownames,
                                                     flowlines = reaches,
                                                     sites = NA,
                                                     STATE.map = boundary,
                                                     map.title = regionName,
                                                     file.name = fn,
                                                     orient = orientation)))

    write.csv(temp_df, file.path(outputFolder, paste0(regionName, "_Clusters_", numtry, ".csv")), row.names = FALSE)

  } # End if user-desired number of clusters
}
