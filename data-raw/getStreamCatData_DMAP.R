#' Get StreamCat Data
#'
#' @param outputFolder Complete file path to the output folder where the boundary rda file is stored
#' @param regionName Name of analysis region
#' @param states Vector of two letter abbreviations of states contained in the custom boundary
#'
#' @return
#' @export
#'
#' @examples
#'
getStreamCatData_custom <- function(outputFolder, regionName, states){

  # Load previously generated reaches file
  reaches_fp <- file.path(outputFolder, paste0(regionName, "_Reaches.rda"))
  load(reaches_fp)

  # Load parquet files for states requested by user
  states_df <- NULL
  for(i in 1:length(states)){
    message(paste0("Retrieving StreamCat data from ", states[i]))

    pq_fp <- paste0("s3://dmap-data-commons-ow/streamcat/", states[i], "_CASTool_StreamCatMetrics.parquet")
    temp_pq <- arrow::read_parquet(pq_fp, anonymous = TRUE)

    states_df <- states_df %>% dplyr::bind_rows(temp_pq)
  }

  # Get all SC parameter names
  data_bkginfo <- read.csv(file.path("data-raw", "SelectedStreamCatStressors.csv"))

  SCmetrics <- StreamCatTools::sc_get_params(param = 'metric_names')
  data_stressorinfoWS <- data.frame(SCmetrics)

  data_stressorinfoWS <- data_stressorinfoWS %>%
    dplyr::mutate(Year = suppressWarnings(
      dplyr::case_when(
        SCmetrics == "popden2010" ~ NA_integer_,
        grepl("^\\w*\\d{4}$", SCmetrics) ~ as.integer(sub("^\\w*(\\d{4})$",
                                                          "\\1", SCmetrics)),
        TRUE ~ NA_integer_))) %>%
    dplyr::mutate(StreamCatVar = suppressWarnings(
      dplyr::case_when(SCmetrics == "popden2010" ~ "popden2010",
                       grepl("\\d{4}$", SCmetrics) ~ sub("\\d{4}", "####", SCmetrics),
                       TRUE ~ SCmetrics))) %>%
    dplyr::select(SCmetrics, StreamCatVar, Year)

  data_stressorinfoWS <- dplyr::full_join(data_stressorinfoWS, data_bkginfo, by = c("StreamCatVar" = "variable")) %>%
    dplyr::filter(!is.na(description))%>%
    dplyr::mutate(StreamCatVar = sub("####", "", StreamCatVar),
                  Label = dplyr::case_when(grepl("^NABD", description) ~ description,
                                           grepl("^NPDES", description) ~ description,
                                           TRUE ~ stringr::str_to_sentence(description))) %>%
    dplyr::select(StreamCatVar, SCmetrics, Year, Label)

  sc_vars <- data_stressorinfoWS %>%
    dplyr::pull(SCmetrics)

  # Find COMIDs in reaches file not retrieved in user-requested state parquet files
  missing_coms <- setdiff(reaches$COMID, states_df$comid)

  message("Retrieve StreamCat data from missing COMIDS")
  tictoc::tic()

  missing_df <- NULL
  if(length(missing_coms)>0){
    tryCatch({
      missing_df <- StreamCatTools::sc_get_data(metric = sc_vars %>% paste(collapse = ","),
                                                aoi = 'watershed',
                                                comid = missing_coms %>% paste(collapse = ","))
    }, error = function(msg){
      print("Encountered error in downloading missing COMIDS. Possibly no additional COMIDS")
    })
  }

  tictoc::toc()

  # Select StreamCat vars from parquet
  sc_vars_pq <- c("comid", paste0(sc_vars, "ws"))
  states_sc_df <- states_df %>% dplyr::select(all_of(sc_vars_pq))

  # Bind user requested states and missing comids
  # Filter to comids in reaches file
  data_stressorWS <- states_sc_df %>%
    dplyr::bind_rows(missing_df)%>%
    dplyr::filter(comid %in% reaches$COMID)

  data_stressorWS <- data_stressorWS %>%
    dplyr::rename_with(~sub("ws$", "", .)) %>%
    tidyr::pivot_longer(cols = !c(comid),
                        names_to = "StreamCatVar",
                        values_to = "WatershedValue")

  data_stressorWS <- data_stressorWS %>%
    dplyr::mutate(Year = suppressWarnings(dplyr::case_when(
      StreamCatVar == "popden2010" ~ NA_integer_,
      grepl("^\\w*\\d{4}$", StreamCatVar) ~
        as.integer(sub("^\\w*(\\d{4})$", "\\1", StreamCatVar)),
      TRUE ~ NA_integer_)),
      StreamCatVar = suppressWarnings(dplyr::case_when(
        StreamCatVar == "popden2010" ~ "popden2010",
        grepl("\\d{4}$", StreamCatVar) ~ sub("\\d{4}", "", StreamCatVar),
        TRUE ~ StreamCatVar))) %>%
    dplyr::rename("COMID" = "comid")


  # write.csv(ret_df, file = file.path(outputFolder, paste0(regionName, "_WSStressor.csv")), row.names = FALSE)
  # write.csv(data_stressorinfoWS, file = file.path(outputFolder, paste0(regionName, "_WSStressorInfo.csv")), row.names = FALSE)

  save(data_stressorWS, file = file.path(outputFolder, paste0(regionName, "_WSStressor.rda")))
  save(data_stressorinfoWS, file = file.path(outputFolder, paste0(regionName, "_WSStressorInfo.rda")))

  message("Watershed stressor data and metadata saved to output folder")
  }
