writeBoundaryRDA <- function(inputFilePath, outputFolder){
  boundary <- sf::st_read(inputFilePath)

  save(boundary, file = file.path(outputFolder, "Boundary.rda"))
}
