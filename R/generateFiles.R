#' generateFiles
#'
#' @return File paths of generated files
#' @export
generateFiles <- function(inputFilePath, outputFolder, regionName, states, clustering = TRUE, wsStressorData = TRUE, ...){

  message("Generating boundary file")
  suppressWarnings(getBoundaryFile(inputFilePath, outputFolder, regionName))

  message("Generating reaches file")
  suppressWarnings(getReachesFile(outputFolder, regionName, states))


  if(clustering){
    message("Generating clusters files")
    suppressWarnings(getClusterFile(outputFolder, regionName, states, ...))
  }


  if(wsStressorData){
    message("Generating watershed stressors file")
    suppressWarnings(getWSStressFiles(outputFolder, regionName, states))
  }

  message(paste0("Files written to ", outputFolder))
}
