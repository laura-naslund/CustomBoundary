#' Get Boundary File
#'
#' @param inputFilePath file path to the geospatial file containing the region boundary
#' @param outputFolder file path to the folder for saving outputs
#' @param region_name name of the region
#'
#' @return writes an RDA file containing the region boundary transformed to EPSG 5070: NAD83/Conus Albers to the output folder
#' @export

getBoundaryFile <- function(inputFilePath, outputFolder, region_name){
  boundary <- sf::st_read(inputFilePath) |>
    sf::st_transform(crs = 5070)

  save(boundary, file = file.path(outputFolder, paste0(region_name, "_Boundary.rda")))
}
