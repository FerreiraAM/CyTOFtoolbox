# Modify the MDS plot from the CytoGLMM package: add a threshold about the markers
# (on the correlation with the MDS dimensions) and ggrepel to display the names of the
# markers

#' MDS on median marker expression
#' 
#' @param df_samples Data.frame, containing the cells intensities.
#' @param protein_names Vector, vector of markers.
#' @param sample_info_names Vector, vector of columns that are not marker's names.
#' @param color Character, indicates how to color the MDS dots.
#' @param sample_label Character, dot labels for MDS plot.
#' @param cor_threshold Numeric, threshold on which marker to display in the circle
#' plot (by default cor_threshold = 0.5).
#' 
#' @export
plot_MDS_withthreshold <- function (df_samples, 
                                    protein_names, 
                                    sample_info_names, 
                                    color, 
                                    sample_label = "",
                                    cor_threshold = 0.5) {
  # Compute the median of each marker
  expr_median <- df_samples %>% 
    dplyr::group_by(.dots = sample_info_names) %>% 
    dplyr::summarise_at(protein_names, stats::median) %>% as.data.frame
  # Compute the distance matrix based on the median
  dist_matrix <- stats::dist(expr_median[, -seq(sample_info_names)])
  # MDS of the distance matrix with the maximum dimension equals to 2
  mds_res <- stats::cmdscale(dist_matrix, eig = TRUE, k = 2)
  # Compute the explained variance
  explained_var <- (100 * mds_res$eig[1:2]/sum(mds_res$eig)) %>% round(digits = 1)
  # Combined the results
  tibble_MDS <- tibble::tibble("MDS1" = mds_res$points[,1], "MDS2" = mds_res$points[,2])
  expr_median <- dplyr::bind_cols(expr_median, tibble_MDS)
  # Compute the standard deviation
  protein_sd <- apply(expr_median[, protein_names], 2, stats::sd)
  # Select the marker with SD different from 0
  protein_selection <- protein_names[protein_sd != 0]
  
  # Correlation plot
  # How much each marker is correlated to MDS1/2?
  expr_cor <- stats::cor(expr_median[, protein_selection], expr_median[, c("MDS1", "MDS2")])
  expr_cor <- tibble::as_tibble(expr_cor)
  expr_cor <- tibble::add_column(expr_cor, protein_selection)
  expr_cor <- tibble::add_column(expr_cor, "x0" = rep(0, nrow(expr_cor)))
  expr_cor <- tibble::add_column(expr_cor, "y0" = rep(0, nrow(expr_cor)))
  # Do not display protein which are < cor_threshold of both MSD1/2
  expr_cor_filtered <- dplyr::filter(expr_cor, 
                                     abs(.data$MDS1) > cor_threshold & abs(.data$MDS2) > cor_threshold)
  # Draw the circle
  corcir <- circle(c(0, 0), npoints = 100)
  # Circle plot
  circle_plot <- ggplot() + 
    geom_path(data = corcir, aes_string(x = "x", y = "y"), colour = "gray65") + 
    geom_hline(yintercept = 0, colour = "gray65") +  
    geom_vline(xintercept = 0, colour = "gray65") + 
    xlim(-1.1, 1.1) + 
    ylim(-1.1, 1.1) + 
    geom_segment(data = expr_cor_filtered,
                 aes_string(x = "x0", y = "y0", xend = "MDS1", yend = "MDS2"), 
                 colour = "gray65") + 
    geom_text_repel(data = expr_cor_filtered, 
                    aes_string(x = "MDS1", y = "MDS2", label = "protein_selection")) + 
    labs(x = "axis 1") + labs(y = "axis 2") + 
    coord_fixed() 
  # MDS plot
  mds_plot <- ggplot(expr_median, aes_string(x = "MDS1", y = "MDS2", color = color)) + 
    geom_point(size = 2) + 
    coord_fixed(ratio = explained_var[2]/explained_var[1]) + 
    xlab(paste0("MDS1 (", explained_var[1], "%)")) + 
    ylab(paste0("MDS2 (", explained_var[2], "%)"))
  # If there is label specify for the MDS plot, add it to the MDS plot
  if (nchar(sample_label) > 1) {
    mds_plot <- mds_plot + geom_label(aes_string(label = "sample_label"))
  }
  cowplot::plot_grid(mds_plot, circle_plot, nrow = 1, rel_widths = c(0.6, 0.4))
}

# HELPER FUNCTIONS =================================================================================

# Draw the circle
# 
# @param center Vector, by default center = c(0,0).
# @param npoint Numeric, number of points.
# @return A tibble.
circle <- function(center = c(0, 0), npoints = 100) {
  r <- 1
  tt <- seq(0, 2 * pi, length = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[1] + r * sin(tt)
  return(tibble::tibble("x" = xx, "y" = yy))
}