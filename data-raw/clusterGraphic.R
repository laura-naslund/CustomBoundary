# Cluster Graphics (Generic)
# Ann.RoseberryLincoln@tetratech.com, 20250103
# Based on code written by Tom Barnum, USEPA, 20240229
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R v4.4.2
#
# library(devtools)
# requires packages: cowplot, dplyr, ggplot2, sf, tibble, tmap, viridis

#' @title Reach Cluster Graphics
#'
#' @description Layout of 5 graphics color-coded by cluster identifier
#' (2-6 clusters are currently accommodated).
#'
#' @details Layout containing a state map with reaches colored by cluster
#' number, bar charts for number of reaches in each cluster and number of
#' sites in each cluster with the same color-coding, a dotplot showing the
#' cluster membership along PC1 and PC2, along with ellipses showing the
#' most dense portions of each cluster, and line graphics showing the 10
#' variables contributing most and least to each of the first two components.
#'
#' @param clusters a dataframe containing COMID and ClusterID
#' @param pca1 a PCA object for the desired state
#' @param flowlines a flowline shapefile for the state derived from NHDPlus
#' @param sites file that links sampling sites to stream reaches (COMIDs)
#' @param STATE.map a polygon shapefile for the desired state
#' @param map.title the name of the desired state or region
#' @param file.name the image filename to output, complete with path
#'
#' @return Writes a png map to the specified output directory.
#
# no examples
#
#' @export

# library(dplyr)
# library(ggplot2)
# library(sf)
# library(tibble)
# library(tmap)
# library(viridis)

### Create function to generate PCA plots using code from this website:
### https://tem11010.github.io/Plotting-PCAs/

### Inputs for ggplot.PCA function
#    df.clusters: clusters <- the exported csv file from the cluster analysis
#    WS.STATE.PCA: pca1 <- pca results in the WS.STATE.PCA.rownames object
#    WS.STATE.region: flowlines <- hydrology shape file of stream reaches
#    STATE.sites: sites <- file that links sampling sites to stream reaches (COMIDs)
#    STATE.shp: STATE.map <- state map to be plotted in the background
#    map.title <- the name of the state; check to make sure it's legible
#    file.name <- the name of the exported tiff

# require(MASS)
# require(cowplot)

clusterGraphic <- function(clusters, pca1, flowlines, sites, STATE.map, map.title, file.name, orient) {

  boo.debug <- FALSE
  if (boo.debug) {
    clusters = df.clusters
    pca1 = WS.STATE.PCA.rownames
    flowlines = WS.STATE.region
    sites = STATE.sites
    STATE.map = STATE.shp
    map.title = state
    file.name = fn
    rda.file.name = fn.rda
  }
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Required small functions:
  # the variable contribution plot has a circle around the variables that has a radius of 1. Here’s some code to make one.
  circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100) {
    r = diameter / 2
    tt <- seq(0, 2 * pi, length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }

  circ <- circleFun(c(0, 0), 2, npoints = 500)

  euc.dist <- function(x1, x2, y1, y2){
    euc <- sqrt(((abs(x2) - x1)^2) + ((abs(y2) - y1)^2))
    return(euc)
  }
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  clusters$pc1 <- pca1$ind$coord[, 1] # add PC1 score to clusters
  clusters$pc2 <- pca1$ind$coord[, 2] # add PC2 score to clusters

  ### extract the data for the variable contributions to each of the pc axes.
  pca.vars <- pca1$var$coord %>% data.frame
  pca.vars$vars <- gsub("completeObs.", "", rownames(pca.vars)) # convert rownames to vars column without removing rownames
  ### Calculate the length of the vector
  pca.vars$vec.dist <- euc.dist(0, pca.vars[, 1], 0, pca.vars[, 2])
  ### Put the dataframe in order based on length of vector
  pca.vars.desc <- pca.vars %>% dplyr::arrange(desc(vec.dist))
  pca.vars.asc <- pca.vars %>% dplyr::arrange(vec.dist)

  # Create color vector (range 1 to 6)
  maxClusterID <- max(as.numeric(clusters$ClusterID))
  if(maxClusterID == 10){
    mag.vec <- viridis::viridis(39)[c(3,7,11,15,19,23, 27,31, 35, 39)]
  } else if(maxClusterID == 9){
    mag.vec <- viridis::viridis(35)[c(3,7,11,15,19,23, 27, 31, 35)]
  } else if(maxClusterID == 8){
    mag.vec <- viridis::viridis(31)[c(3,7,11,15,19,23, 27, 31)]
  } else if(maxClusterID == 7){
    mag.vec <- viridis::viridis(27)[c(3,7,11,15,19,23, 27)]
  } else if (maxClusterID == 6) {
    mag.vec <- viridis::viridis(23)[c(3,7,11,15,19,23)]
  } else if (maxClusterID == 5) {
    mag.vec <- viridis::viridis(19)[c(3,7,11,15,19)]
  } else if (maxClusterID == 4) {
    mag.vec <- viridis::viridis(15)[c(3,7,11,15)]
  } else if (maxClusterID == 3) {
    mag.vec <- viridis::viridis(11)[c(3,7,11)]
  } else if (maxClusterID == 2) { # max clusters = 2
    mag.vec <- viridis::viridis(7)[c(3,7)]
  }
    else{
       mag.vec <- c("#31688EFF")
    }
  # mag.vec <- viridis::viridis(15)[c(3,7,11,15)]

  # Prepare PCA plot (individuals) ----
  pca.ind <- as.data.frame(pca1$ind$coord) %>%
    tibble::rownames_to_column(var = "COMID") %>%
    dplyr::mutate(COMID = as.integer(COMID))
  pca.ind <- dplyr::full_join(clusters, pca.ind, relationship = "one-to-one")
  pca.ind <- dplyr::mutate(pca.ind, ClusterID = as.character(ClusterID))

  pca12 <- ggplot2::ggplot(data = pca.ind, ggplot2::aes(x = Dim.1, y = Dim.2,
                                                        color = ClusterID)) +
    ggplot2::geom_hline(yintercept = 0, lty = 2, color = "gray45", alpha = 0.9) +
    ggplot2::geom_vline(xintercept = 0, lty = 2, color = "gray45", alpha = 0.9) +
    ggplot2::geom_point(alpha = 0.5, size = 0.7) +
    ggplot2::scale_shape_manual(values = mag.vec) +
    ggplot2::scale_color_manual(values = mag.vec) +
    ggplot2::stat_ellipse(geom = "polygon", ggplot2::aes(fill = ClusterID),
                          alpha = 0.2, show.legend = FALSE, level = 0.95) +
    ggplot2::scale_fill_manual(values = mag.vec) +
    ggplot2::xlab(paste("PC 1 (", round(pca1$eig[1, 2], 1), "%)", sep = "")) +
    ggplot2::ylab(paste("PC 2 (", round(pca1$eig[2, 2], 1), "%)", sep = "")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = "transparent",
                                                          color = "black",
                                                          linewidth = 0.5),
                   axis.text.y = ggplot2::element_text(size = 12, face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 12, face = "bold"
                                                         , vjust = 0.15, hjust = 0.1),
                   axis.title = ggplot2::element_text(size = 12, face = "bold"),
                   plot.title = ggplot2::element_text(size = 12, face = "bold"),
                   axis.ticks.length.y = ggplot2::unit(4, "pt"),
                   legend.direction = "horizontal",
                   legend.title = ggplot2::element_text(size = 10, face = "bold"),
                   legend.text = ggplot2::element_text(size = 10, face = "bold"),
                   legend.position = "none")

  # Prepare top 10 variable contribution ----
  vars.max <- ggplot2::ggplot() +
    ggplot2::geom_path(data = circ, ggplot2::aes(x, y), lty = 2, color = "grey",
                       alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    ggplot2::geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    ggplot2::geom_segment(data = pca.vars.asc,
                          ggplot2::aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.025, "npc")
                                                 , type = "open"),
                          lwd = 1, color = c(rep("grey95", nrow(pca.vars.asc) - 10),
                                             rep("black", 10))) +
    ggrepel::geom_text_repel(data = pca.vars.desc[1:10, ], size = 3,
                             ggplot2::aes(x = Dim.1, y = Dim.2, label = paste(vars)), max.overlaps = 15) +
    ggplot2::xlab("PC 1") +
    ggplot2::ylab("PC 2") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = "transparent",
                                                        color = "black",
                                                        linewidth = 0.5),
                   axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 10, face = "bold",
                                                       vjust = 0.15, hjust = 0.1),
                   axis.title = ggplot2::element_text(size = 12, face = "bold"),
                   plot.title = ggplot2::element_text(size = 12, face = "bold"),
                   axis.ticks.length.y = ggplot2::unit(4, "pt"))

  row.1 <- nrow(pca.vars.desc) - 10
  row.2 <- nrow(pca.vars.desc)

  # Prepare bottom 10 variable contribution ----
  vars.min <-  ggplot2::ggplot() +
    ggplot2::geom_path(data = circ, ggplot2::aes(x, y), lty = 2
                       , color = "grey", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    ggplot2::geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    ggplot2::geom_segment(data = pca.vars.desc
                          , ggplot2::aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2)
                          , arrow = ggplot2::arrow(length = ggplot2::unit(0.025, "npc")
                                                   , type = "open")
                          , lwd = 1, color = c(rep("grey95", nrow(pca.vars.desc) - 10)
                                               , rep("black", 10))) +
    ggrepel::geom_text_repel(data = pca.vars.desc[row.1:row.2, ], size = 3,
                             ggplot2::aes(x = Dim.1, y = Dim.2, label = paste(vars)), max.overlaps = 15) +
    ggplot2::xlab("PC 1") +
    ggplot2::ylab("PC 2") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = "transparent",
                                                        color = "black",
                                                        linewidth = 0.5),
                   axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 10, face = "bold",
                                                       vjust = 0.15, hjust = .1),
                   axis.title = ggplot2::element_text(size = 12, face = "bold"),
                   plot.title = ggplot2::element_text(size = 12, face = "bold"),
                   axis.ticks.length.y = ggplot2::unit(4, "pt"))

  # Calculate site/reach counts in each cluster ----
  comid.count <- as.data.frame(table(clusters$ClusterID))
  Total <- sum(comid.count$Freq)

  if(!is.na(sites)){ # LCN 1/31/25 added
    site.df <- dplyr::left_join(sites, clusters[c("COMID", "ClusterID")]
                                , relationship = "many-to-one")
    site.unique <- site.df[!duplicated(site.df$COMID), ]
    site.count <- as.data.frame(table(site.unique$ClusterID))
    site.Total <- sum(site.count$Freq)
    x.pt.1 <- max(comid.count[,2]) - (max(comid.count[,2]) * 0.08)
    x.pt.2 <- max(site.count[,2]) - (max(site.count[,2]) * 0.08)
  }



  # Prepare reaches/cluster bar graph ----
  comid.title <- paste0("Number of COMIDs per cluster (N = ", Total, ")")
  all.comids.count <- ggplot2::ggplot(comid.count
                                      , ggplot2::aes(x = Freq, y = as.factor(Var1)
                                                     , fill = as.factor(Var1))) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()
                      , color = mag.vec, fill = mag.vec, alpha = 0.8) +
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "Cluster ID", x = "Frequency", y = "Cluster ID") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::geom_text(ggplot2::aes(label = Freq, hjust = ifelse(Freq >= 0.8 * max(Freq)
                                                                 , 1.5, -0.5))
                       , position = ggplot2::position_dodge(width = 0.9), size = 3.5
                       , fontface = "bold") +
    ggplot2::theme_classic() +
    ggplot2::scale_y_discrete(breaks = c("1", "2", "3", "4"),
                              labels = c("1", "2", "3", "4")) +
    ggplot2::ggtitle(comid.title) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10, face = "bold"),
                   legend.text = ggplot2::element_text(size = 10, face = "bold"),
                   legend.position = "top",
                   axis.text = ggplot2::element_text(size = 11, face = "bold"),
                   axis.title = ggplot2::element_text(size = 11, face = "bold"),
                   strip.text = ggplot2::element_text(size = 12, face = "bold"),
                   axis.title.x = ggplot2::element_blank())

  # Prepare sites/cluster bar graph ----

  if(!is.na(sites)){
    site.title = paste0("Number of sample sites per cluster (N = ", site.Total, ")")
    comid.site.count <- ggplot2::ggplot(site.count
                                        , ggplot2::aes(x = Freq, y = as.factor(Var1)
                                                       , fill = as.factor(Var1))) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()
                        , color = mag.vec, fill = mag.vec, alpha = 0.8) +
      ggplot2::theme_classic() +
      ggplot2::labs(fill = "Cluster ID", x = "Frequency", y = "Cluster ID") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::geom_text(ggplot2::aes(label = Freq, hjust = ifelse(Freq >= 0.8 * max(Freq)
                                                                   , 1.5, -0.5))
                         , position = ggplot2::position_dodge(width = 0.9), size = 3.5
                         , fontface = "bold") +
      ggplot2::theme_classic() +
      ggplot2::scale_y_discrete(breaks = c("1", "2", "3", "4"),
                                labels = c("1", "2", "3", "4")) +
      ggplot2::ggtitle(site.title) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 10, face = "bold"),
                     legend.text = ggplot2::element_text(size = 10, face = "bold"),
                     legend.position = "top",
                     axis.text = ggplot2::element_text(size = 11, face = "bold"),
                     axis.title = ggplot2::element_text(size = 11, face = "bold"),
                     strip.text = ggplot2::element_text(size = 12, face = "bold"),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank())
  }


  # Prepare map ----
  # Get state bounding box
  bbox_new <- sf::st_bbox(STATE.map)
  xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
  yrange <- bbox_new$ymax - bbox_new$ymin # range of y values

  # State name will always appear above the state
  # bbox_new[1] <- bbox_new[1] - (0.25 * xrange) # xmin - left
  # bbox_new[3] <- bbox_new[3] + (0.25 * xrange) # xmax - right
  # bbox_new[2] <- bbox_new[2] - (0.25 * yrange) # ymin - bottom
  bbox_new[4] <- bbox_new[4] + (0.05 * yrange) # ymax - top

  bbox_new <- bbox_new %>%  # take the bounding box ...
    sf::st_as_sfc() # ... and make it a sf polygon

  NHD.clust <- dplyr::right_join(flowlines, clusters[c("COMID", "ClusterID")])
  names(NHD.clust)[names(NHD.clust) == 'ClusterID'] <- 'Cluster'
  state.map <- tmap::tm_shape(STATE.map, bbox = bbox_new) +
    tmap::tm_polygons(fill = "grey80") +
    tmap::tm_shape(NHD.clust) +
    tmap::tm_lines("Cluster", palette = mag.vec, legend.col.show = FALSE) +
    tmap::tm_shape(STATE.map) +
    tmap::tm_borders(col = "black", lwd = 1) +
    tmap::tm_layout(frame = FALSE, legend.show = FALSE, main.title = map.title,
                    main.title.size = 1, main.title.fontface = "bold",
                    # panel.show = TRUE, panel.labels = map.title,
                    # panel.label.bg.color = "white",
                    # panel.label.fontface = "bold",
                    # panel.label.height = 1.2,
                    inner.margins = c(0,0,0,0))+
    #LCN
    tmap::tm_options(component.autoscale = FALSE)
  state.map.grob <- tmap::tmap_grob(state.map)

  # Change layout according to to state orientation

  ## States not accommodated currently
  # Alaska -- not in NHD+
  # Florida
  # Hawaii -- not in NHD+
  # Maine
  # Michigan
  # Missouri
  # New York
  # Texas
  # Virginia - horizontal

  # Landscape layout ----
  # Sort of squarish or mostly horizontal states (tested with OR and WA)
  if (map.title %in% c("Arizona", "Arkansas", "Colorado", "Connecticut",
                       "Georgia", "Iowa", "Kansas", "Kentucky", "Louisiana",
                       "Maryland", "Massachusetts", "Montana", "Nebraska",
                       "New Mexico", "North Carolina", "North Dakota", "Ohio",
                       "Oklahoma", "Oregon", "Pennsylvania", "South Carolina",
                       "South Dakota", "Tennessee", "Utah", "Virginia",
                       "Washington", "West Virginia", "Wisconsin", "Wyoming") | orient == "horizontal") {

    # Prepare combination w/map upper left, 2 bar charts upper right,
    # and PCA plots along bottom
    if(!is.na(sites)){ # LCN added
      return.top.right <- cowplot::plot_grid(comid.site.count, all.comids.count,
                                             labels = c("B", "C"), nrow = 2)
    }
    else{
      return.top.right <- cowplot::plot_grid(all.comids.count,
                                             labels = c("B", "C"))
    }


    return.top <- cowplot::plot_grid(state.map.grob, return.top.right, ncol = 2,
                                     labels = c("A", ""), rel_heights = c(0.9, 1))

    return.bottom <- cowplot::plot_grid(pca12, vars.max, vars.min,
                                        rel_heights = c(0.6, 0.6, 0.6),
                                        ncol = 3, labels = c("D", "E", "F"))

    return.p <- cowplot::plot_grid(return.top, return.bottom, nrow = 2)

    png(file.name, width = 12, height = 8, units = "in", res = 600)

    return.p %>% print()

    dev.off()


  } else if (map.title %in% c("Alabama", "California", "Delaware", "Idaho",
                              "Illinois", "Indiana", "Minnesota", "Mississippi",
                              "Nevada", "New Hampshire", "New Jersey", "Utah",
                              "Vermont") | orient == "vertical") { # mostly vertical states (not tested with any)

    # Portrait layout ----

    # TODO: Figure out how to maximize map size to fill available space
    # Prepare combination w/map upper left, 2 bar charts bottom left (stacked),
    # and PCA plots along right side (stacked)
    if(!is.na(sites)){
      return.bottom.left <- cowplot::plot_grid(comid.site.count, all.comids.count,
                                               labels = c("B", "C"), nrow = 2)
    }
    else{
      return.bottom.left <- cowplot::plot_grid(all.comids.count,
                                               labels = c("B", "C"))
    }



    return.left <- cowplot::plot_grid(state.map.grob, return.bottom.left, nrow = 2,
                                      labels = c("A", ""), rel_heights = c(1, 0.6))

    return.right <- cowplot::plot_grid(pca12, vars.max, vars.min,
                                       rel_heights = c(0.6, 0.6, 0.6),
                                       nrow = 3, labels = c("D", "E", "F"))

    return.p <- cowplot::plot_grid(return.left, NULL, return.right, ncol = 3,
                                   rel_widths = c(0.45, 0.1, 0.45))

    png(file.name, width = 8, height = 12, units = "in", res = 600)

    return.p %>% print()

    dev.off()
  }

}
