writeBoundaryRDA <- function(inputFilePath, outputFolder, region_name){
  boundary <- sf::st_read(inputFilePath)

  save(boundary, file = file.path(outputFolder, paste0(region_name, "_Boundary.rda")))
}
