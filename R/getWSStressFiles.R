#' Get watershed stressor data and metadata
#'
#' @param outputFolder file path containing output files
#' @param regionName region name
#' @param states vector of states in the analysis region
#'
#' @return Writes a csv with watershed stressor data within the analysis region formatted for the CASTool and watershed stressor metadata
#' @export
getWSStressFiles <- function(outputFolder, regionName, states){
  Sys.setenv(AWS_EC2_METADATA_DISABLED = "true")

  ws_ret <- data(list = "WSStressorInfo", package = "CustomBoundary", envir = environment())
  ws_info <- get(ws_ret, envir = environment())

  ws_stressors_vec <- ws_info |> dplyr::pull(SCmetrics)

  load(file.path(outputFolder, paste0(regionName, "_Reaches.rda")))

  reach_comids <- reaches |> sf::st_drop_geometry() |> dplyr::pull(COMID)

  ws_data <- NULL

  for(i in seq_len(length(states))){
    temp_state <- states[i]

    if(nchar(temp_state) != 2){
      temp_state <- state.abb[which(state.name == temp_state)]
    }

    message(paste0("Retrieving watershed stressor data from ", temp_state))

    state_fp <- paste0("s3://dmap-data-commons-ow/streamcat/CASTool_State_SC/", temp_state, "_CASTool_StreamCatMetrics.parquet")

    state_pq <- arrow::open_dataset(state_fp) |>
      dplyr::collect() |>
      dplyr::select(!dplyr::ends_with(".x")) |>
      dplyr::rename_with(~ stringr::str_remove(.x, "\\.y$"), .cols = dplyr::ends_with(".y")) |>
      dplyr::rename_all(~ stringr::str_remove(.x, "ws$")) |>
      dplyr::select(dplyr::all_of(c("comid", ws_stressors_vec))) |>
      dplyr::filter(comid %in% reach_comids)

    ws_data <- ws_data |> dplyr::bind_rows(state_pq)
  }

  ws_data_ret <- ws_data |>
    tidyr::pivot_longer(cols = !comid, names_to = "SCmetrics", values_to = "WatershedValue") |>
    dplyr::left_join(ws_info, by = "SCmetrics") |>
    dplyr::select(comid, StreamCatVar, WatershedValue, Year) |>
    dplyr::rename("COMID" = "comid")

  write.csv(ws_data_ret, file.path(outputFolder, paste0(regionName, "_WSStressor.csv")), row.names = FALSE)

  write.csv(ws_info, file.path(outputFolder, paste0(regionName, "_WSStressorInfo.csv")), row.names = FALSE)
}
