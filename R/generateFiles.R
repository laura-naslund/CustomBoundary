#' generateFiles
#'
#' @return File paths of generated files
#' @export
#'
#' @examples
generateFiles <- function(inputFilePath, outputFolder, regionName, wsStressorData = FALSE, states, ...){
  source("data-raw/writeBoundaryRDA.R")
  source("data-raw/clusterReaches_customBoundary.R")
  source("data-raw/getStreamCatData_DMAP.R")

  message("Generating boundary file")
  suppressWarnings(writeBoundaryRDA(inputFilePath, outputFolder, regionName))

  message("Generating clusters file")
  suppressWarnings(clusterReachesCustom(outputFolder, regionName, ...))

  message("Generating watershed stressors file")
  if(wsStressorData){
    suppressWarnings(getStreamCatData_custom(outputFolder, regionName, states))
  }

  ws_stressor_out <- ""
  ws_stressorinfo_out <- ""

  if(wsStressorData == TRUE){
    ws_stressor_out <- paste0("Watershed stressor file: ", outputFolder, "/", regionName, "_WSStressor.rda\n")
    ws_stressorinfo_out <- paste0("Watershed stressor metadata file: ", outputFolder, "/", regionName, "_WSStressorInfo.rda\n")
  }

 cat("\n\nBoundary file: ", outputFolder, "/", regionName, "_Boundary.rda\n",
             "Cluster graphic file: ", outputFolder, "/", regionName, "_ClusterGraphic.png\n",
             "Cluster assignment file: ", outputFolder, "/", regionName, "_Clusters.csv\n",
             "Reaches file: ", outputFolder, "/", regionName, "_Reaches.rda\n",
             ws_stressor_out,
             ws_stressorinfo_out,
             sep = "")
}
