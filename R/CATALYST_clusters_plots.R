#' Extract SCE information in preparation to plot the clusters.

#' Extract SCE information in preparation to plot the clusters.
#'
#' @param sce SingleCellExperiment.
#' @param metacluster Character, metacluster to visualize.
#' @param dr Character, dimensionality reduction to use (by default dr = UMAP").
#' @return Data.frame.
#' @export
extract_SCE_info <- function(sce, metacluster, dr = "UMAP"){
  # Extract DR info 
  dr_values <- reducedDim(sce, dr)
  ## Check colnames
  if (colnames(dr_values) != c("X1", "X2")) {
    colnames(dr_values) <- c("X1", "X2")
  }
  # Combine with experiment info
  ## Combine colData with dimension reduction values
  expinfo_dr_values <- data.frame(colData(sce), dr_values)
  # Get cluster IDs info
  sce_metacluster_ids <- as.data.frame(cluster_ids(sce, metacluster))
  colnames(sce_metacluster_ids) <- metacluster
  # Combine with experiment info and DR
  data.frame(expinfo_dr_values, sce_metacluster_ids)
}

#' Plot clusters.
#' 
#' Plot clusters.
#' @param df Data.frame, output of the extract_SCE_info function.
#' @param metacluster Character, metacluster to visualize.
#' @param point_size Numeric, (by default = 0.5)
#' @param point_alpha Numeric, (by default = 0.5)
#' @return ggplot2 plot.
#' @export
plot_metaclusters <- function(df, metacluster, point_size = 0.5, point_alpha = 0.5){
  ## CATALYST colors for meta-cluster 
  ## (keeping the cluster colors defined by the CATALYST package)
  catalyst_col <- CATALYST:::.cluster_cols
  # Plot
  ggplot(df, aes_string(x = "X1", y = "X2")) +
    geom_point(aes_string(col = metacluster), size = point_size, alpha = point_alpha, pch = 16) +
    scale_color_manual(values = catalyst_col, name = paste(metacluster, "meta-clusters")) +
    labs(x = paste(dr, "1"), y = paste(dr, "2")) + 
    theme_bw() + 
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1)))
}