#' Get Reaches File
#' @param outputFolder file path to the folder for saving outputs
#' @param region_name name of the region
#'
#' @return writes an rda file containing the comids and geometry of flowines within the custom region to the outputFolder
#' @export

getReachesFile <- function(outputFolder, region_name, states){
  load(file.path(outputFolder, paste0(region_name, "_Boundary.rda")))

  states_df <- NULL

  for(i in seq_len(length(states))){
    temp_state <- states[i]

    message(paste0("Retrieving reaches from ", temp_state))

    if(nchar(temp_state) == 2){
      temp_state <- state.name[which(state.abb == temp_state)]
    }

    temp_df <- CASToolHelperPckg::getReaches(temp_state)

    states_df <- states_df |> dplyr::bind_rows(temp_df)

  }

  # Distinct because 300 m boundaries of states will overlap
  states_df <- states_df |>
    dplyr::distinct()

  state_index_flip <- sf::st_intersects(boundary, states_df)

  reaches <- states_df[state_index_flip[[1]],]

  save(reaches, file = file.path(outputFolder, paste0(region_name, "_Reaches.rda")))
}
