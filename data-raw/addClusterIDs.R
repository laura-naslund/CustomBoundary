# Add Cluster IDs to Sites
# Ann.RoseberryLincoln@tetratech.com, 20240508
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R v4.4.2
#
# library(devtools)
# install_github("ALincolnTt/ReachClustering")
# requires packages: dplyr, ggplot2, viridis

#' @title Add Cluster IDs to State Sites File
#'
#' @description Reads a state-specific sites file and adds a column "ClusterID".
#' Generates box plots for #reaches per cluster and #sites per cluster.
#' Currently accommodates between 2 and 6 clusters.
#'
#' This code is run when the clusterReaches.R code is run, but can be run again
#' at the command-line when new sites are added to the monitoring program.
#'
#' @details Both cluster and sites files are expected in csv format. If a
#' different clustering method was used (e.g., by ecoregion), then the sites
#' file should be subset into subregions according to the clustering method
#' and output files.
#'
#' @param clusters either a dataframe or a csv file containing COMID,
#' US_L3CODE, and ClusterID; output by clusterReaches.R
#' @param sites either a dataframe or a csv file containing COMID,
#' StationID, FinalLatitude, FinalLongitude, HorizontalDatum, RefSiteFlag,
#' and EcoL3Code (additional columns will neither be deleted nor used).
#' @param state the 2-character state abbreviation
#' @param outdir the directory into which to write the modified sites file,
#' with full path
#'
#' @return Writes a tab-delimited text file with two additional columns
#' compared to the input sites file: ClusterID and StationNote.
#
# no examples
#
#' @export


addClusterIDs <- function(clusters, sites, state, outdir) {

  boo.DEBUG <- FALSE
  if (boo.DEBUG) {
    out.dir <- "C:/Users/ann.lincoln/Documents/CASTool_DATA/WA/Clusters"
    clusters <- file.path(out.dir, "WashingtonClusterIDs_DEBUG_",
                          format(lubridate::now(), "%Y%m%d%h%m"), ".csv")
    sites <- file.path(out.dir, "WashingtonSites.csv")
    state <- "WA"
  }

  dated <- format(lubridate::now(), "%Y%m%d")

  if (!is.data.frame(clusters)) {
    if (file.exists(clusters)) {
      df.clusters <- read.csv(clusters, stringsAsFactors = FALSE
                              , header = TRUE, na.strings = c("", "NA", "na"))

    } else {
      message("Cluster file does not exist")
    }
  } else {
    df.clusters <- clusters
  }

  if (!is.data.frame(sites)) {
    if (file.exists(sites)) {
      STATE.sites <- read.csv(sites, stringsAsFactors = FALSE
                           , header = TRUE, na.strings = c("", "NA", "na"))

    } else {
      message("Sites file does not exist")
    }
  } else {
    STATE.sites <- sites
  }

  STATE.sites.clusterID <- merge(STATE.sites, df.clusters
                                 , by.x = c("COMID")
                                 , by.y = c("COMID")
                                 , all.x = TRUE)
  STATE.sites.clusterID.good <- STATE.sites.clusterID %>%
    dplyr::mutate(StationNote = "")
  STATE.sites.NOclusterID <- STATE.sites.clusterID %>%
    dplyr::filter(is.na(ClusterID)) %>%
    dplyr::mutate(StationNote = "COMID not in NHDPlus version used OR not provided")
  STATES.sites.clusterID.final <- rbind(STATE.sites.clusterID.good,
                                        STATE.sites.NOclusterID)

  fn <- paste0(state, "SitesWClusters_", dated, ".tab")
  write.table(STATES.sites.clusterID.final, file.path(out.dir, fn)
              , sep = "\t", col.names = TRUE, row.names = FALSE)


  # Create color vector (range 6 to 2)
  maxClusterID <- max(as.numeric(STATES.sites.clusterID.final$ClusterID),
                      na.rm = TRUE)
  if (maxClusterID == 6) {
    mag.vec <- viridis::viridis(23)[c(3,7,11,15,19,23)]
  } else if (maxClusterID == 5) {
    mag.vec <- viridis::viridis(19)[c(3,7,11,15,19)]
  } else if (maxClusterID == 4) {
    mag.vec <- viridis::viridis(15)[c(3,7,11,15)]
  } else if (maxClusterID == 3) {
    mag.vec <- viridis::viridis(11)[c(3,7,11)]
  } else { # max clusters = 2
    mag.vec <- viridis::viridis(7)[c(3,7)]
  }

  # Prepare count of reaches and sites in each cluster
  comid.count <- as.data.frame(table(df.clusters$ClusterID))
  names(comid.count)[1] <- "ClusterID"
  Total <- sum(comid.count$Freq)

  site.unique <- STATES.sites.clusterID.final[!duplicated(STATES.sites.clusterID.final$COMID), ]
  site.count <- as.data.frame(table(site.unique$ClusterID))
  names(site.count)[1] <- "ClusterID"
  site.Total <- sum(site.count$Freq)

  # Prepare bar graph of reaches in each cluster
  comid.title <- paste0("Number of COMIDs per cluster (N = ", Total, ")")
  all.comids.count <- ggplot2::ggplot(comid.count
                                      , ggplot2::aes(x = Freq, y = ClusterID
                                                     , fill = ClusterID)) +
  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()
                      , color = mag.vec, fill = mag.vec, alpha = 0.8) +
    ggplot2::labs(fill = "Cluster ID", x = "Frequency", y = "Cluster ID") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::geom_text(ggplot2::aes(label = Freq, hjust = ifelse(Freq >= 0.8 * max(Freq)
                                                                 , 1.5, -0.5))
                       , position = ggplot2::position_dodge(width = 0.9), size = 3.5
                       , fontface = "bold") +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(comid.title) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10, face = "bold"),
                   legend.text = ggplot2::element_text(size = 10, face = "bold"),
                   legend.position = "top",
                   axis.text = ggplot2::element_text(size = 11, face = "bold"),
                   axis.title = ggplot2::element_text(size = 11, face = "bold"),
                   strip.text = ggplot2::element_text(size = 12, face = "bold"),
                   axis.title.x = ggplot2::element_blank())
  ggplot2::ggsave(file.path(out.dir, paste0(state, "_NumCOMIDsPerCluster_", dated, ".png"))
                  , all.comids.count, dpi = 600, width = 6, height = 4, units = "in")

  # Prepare bar graph of number of sites per cluster (upper right)
  site.title = paste0("Number of sample sites per cluster (N = ", site.Total, ")")
  comid.site.count <- ggplot2::ggplot(site.count, ggplot2::aes(x = Freq, y = ClusterID
                                                     , fill = ClusterID)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()
                      , color = mag.vec, fill = mag.vec, alpha = 0.8) +
    ggplot2::labs(fill = "Cluster ID", x = "Frequency", y = "Cluster ID") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::geom_text(ggplot2::aes(label = Freq, hjust = ifelse(Freq >= 0.8 * max(Freq)
                                                                 , 1.5, -0.5))
                       , position = ggplot2::position_dodge(width = 0.9), size = 3.5
                       , fontface = "bold") +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(site.title) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10, face = "bold"),
                   legend.text = ggplot2::element_text(size = 10, face = "bold"),
                   legend.position = "top",
                   axis.text = ggplot2::element_text(size = 11, face = "bold"),
                   axis.title = ggplot2::element_text(size = 11, face = "bold"),
                   strip.text = ggplot2::element_text(size = 12, face = "bold"),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  ggplot2::ggsave(file.path(out.dir, paste0(state, "_NumSitesPerCluster_", dated, ".png"))
                  , comid.site.count, dpi = 600, width = 6, height = 4, units = "in")

}
