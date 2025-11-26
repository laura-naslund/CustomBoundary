#' generateFiles
#'
#' @return File paths of generated files
#' @export
#'
#' @examples
generateFiles <- function(inputFilePath, outputFolder, region_name, ...){
  source("data-raw/writeBoundaryRDA.R")
  source("data-raw/clusterReaches_customBoundary.R")

  writeBoundaryRDA(inputFilePath, outputFolder)

  suppressWarnings(clusterReachesCustom(outputFolder, region_name, ...))

  return(cat("\n\nBoundary file: ", outputFolder, "/Boundary.rda\n",
                "Cluster assignment file: ", outputFolder, "/Clusters.csv\n",
                "Reaches file: ", outputFolder, "/Reaches.rda", sep = ""))
}
