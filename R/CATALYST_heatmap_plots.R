# Create function to plot the DA test heatmap

#' Plot the heatmap of the DA test
#'
#' @param x SingleCellExperiment.
#' @param y List, results of the diffcyt DA test.
#' @param comparison Character, on which comparison the DA test has been performed.
#' @param top_n Numeric, number of top cluster to display (by default top_n = 20).
#' @param all Logical, if TRUE all the clusters are displayed 
#' (if TRUE, top_n is ignored; by default all = FALSE).
#' @param th Numeric, p-value threshold (by default = 0.05).
#' @param order Logical, order the results by significance (by default order = TRUE).
#' @param normalize Logical, z-score normalization of the value. 
#' Relative population abundance are arcsine-square-root scaled prior to normalization.
#' (By default normalize = TRUE).
#' @param show_sample_ID Logical, specify is the sample IDs are displayed.
#' (By default show_sample_ID = TRUE).
#' @return HeatmapList-class object.
#' @export
myplotDAheatmap <- function(x, y, comparison, top_n = 20,
                            all = FALSE, th = 0.05, normalize = TRUE,
                            order = TRUE, show_sample_ID = TRUE){
  # Check if the input is an SCE object and contains cluster's information
  CATALYST:::.check_sce(x)
  # Extract expression values of the SCE
  es <- assay(x, "exprs")
  # Perform some checks on the input parameters
  stopifnot(is.numeric(top_n), length(top_n) == 1, is.logical(order), 
            length(order) == 1, is.numeric(th), length(th) == 1, 
            is.logical(normalize), 
            length(normalize) == 1)
  # # If show_sample_ID is true
  # ## The name of sample_ID should be provided
  # if (show_sample_ID & is.na(name_sample_ID)) {
  #   stop("The name of the sample ID column must be provided")
  # }
  # Check on which meta cluster the DA test has been performed
  stopifnot(!is.null(k <- metadata(y$res)$clustering_name))
  # Check that the meta cluster used for the DA test in the SCE object
  k <- CATALYST:::.check_validity_of_k(x, k)
  # Extract the cluster information about each cell
  x$cluster_id <- cluster_ids(x, k)
  # Extract the metadata information about each cell 
  factors <- dplyr::select(as.data.frame(colData(x)), -c("sample_id", 
                                                         "cluster_id"))
  # Extract the results of the DA test
  y <- rowData(y$res)
  # Identify the type of test
  type <- CATALYST:::.get_dt_type(y)
  # If order is true
  ## Order: specify if the results should be ordered by significance
  if (order) {
    # Reorder the results
    y <- y[order(y$p_adj), , drop = FALSE]
  }
  # If all is true or top_n superior than the number of lines in the test 
  ## all: specify if all the cluster should be displayed
  ## top_n: specify the number of top clusters
  if (all | top_n > nrow(y)) {
    # Assign the number of row of the result to top_n
    top_n <- nrow(y)
  }
  # Extract the results of the DA for the top_n clusters
  top <- as.data.frame(y[seq_len(top_n), ])
  ## Convert as character if any factor
  top <- dplyr::mutate_if(top, is.factor, as.character)
  
  # #### 1st heatmap: median type-marker expression by cluster ####
  # # If hm1 is true
  # ## hm1: left-hand side heatmap should be plotted
  # ## I did not change this part
  # if (hm1) {
  #   ms_by_k <- t(CATALYST:::.agg(x[type_markers(x)], "cluster_id"))[top$cluster_id,]
  #   qs <- quantile(ms_by_k, probs = c(0.01, 0.5, 0.99), na.rm = TRUE)
  #   hm_cols <- colorRamp2(qs, c("royalblue3", "white", "tomato2"))
  #   hm1 <- CATALYST:::.diff_hm(ms_by_k, hm_cols, "expression", cluster_rows = !order, 
  #                              xlab = "type_markers", row_title = "cluster_id"[!is.null(hm1)], 
  #                              row_names_side = "left")
  # } else {
  #   hm1 <- NULL
  # }
  
  #### TODO: probably a better way to do this ####
  # Code from CATALYST
  # Match the levels of the sample ID
  ## Match function returns the position of the first match
  m <- match(levels(x$sample_id), x$sample_id)
  # Get the patient IDs
  df <- data.frame(factors[m, ], row.names = NULL)
  ####
  
  #### Relative abundance by cluster ####
  # Count the numbers of cell in each cluster per sample
  cnts <- table(x$cluster_id, x$sample_id)
  # Convert frequency table into proportions
  frqs <- prop.table(cnts, 2)
  # Select the top clusters
  frqs <- frqs[top$cluster_id, ]
  # As matrix
  frqs <- as.matrix(unclass(frqs))
  
  #### Normalization ####
  # If normalize is true
  ## normalize: speficy if the values needs to be Z-score normalized
  ## in case it is a DA analysis, the relative population abundances are 
  ## arcsine-square-root prior to normalization
  if (normalize) {
    frqs <- CATALYST:::.z_normalize(asin(sqrt(frqs)))
    at <- seq(-2.5, 2.5, 0.5)
    labels <- at
    labels[-seq(2, length(at), 2)] <- ""
  }  else {
    min <- floor(min(frqs)/0.1) * 0.1
    max <- ceiling(max(frqs)/0.1) * 0.1
    at <- seq(min, max, 0.1)
    labels <- at
  }
  
  #### Heatmap ####
  # # Extract the values of the comparison
  # comparison_values <- unique(metadata(x)$experiment_info[,comparison])
  # ## Check if there are only 2
  # if (length(comparison_values) != 2) {
  #   stop("The number of values in the comparison is not 2")
  # }
  # Heatmap annotation
  # Get the different levels for the conditions
  lvls <- lapply(as.list(df), levels)
  # Get the number of levels per conditions
  nlvls <- vapply(lvls, length, numeric(1))
  # Get the colors for the different levels
  colors_for_levels <- sample(scales::hue_pal()(sum(nlvls)))
  # Assign a color to each levels
  uniq_combinaison_lvls <- unique(tidyr::gather(df, 
                                                key = "comparison", 
                                                value = "levels"))
  colors_withinfo <- cbind(colors_for_levels, uniq_combinaison_lvls)
  # As list
  cols <- sapply(names(nlvls), FUN = function(onecomp){
    onecomp_values <- dplyr::filter(colors_withinfo, comparison == onecomp)
    droplevels(tibble::deframe(onecomp_values[,c("levels", "colors_for_levels")]))
  }, USE.NAMES = TRUE)
  # Create a annotation for the Heatmap
  mycol_anno <- HeatmapAnnotation(which = "column", df = df, col = cols, 
                                  gp = gpar(col = "white"), 
                                  show_legend = c(TRUE, show_sample_ID, TRUE))
  # PLot the heatmap 
  myheatmap <- Heatmap(matrix = frqs,
                       col = rev(brewer.pal(9, "RdGy")),
                       name = paste0("normalized\n"[normalize], "frequency"),
                       cluster_columns = FALSE, 
                       column_title = "samples", 
                       column_title_side = "bottom",
                       clustering_distance_rows = "euclidean", 
                       clustering_method_rows = "median", 
                       column_names_gp = gpar(fontsize = 8), 
                       rect_gp = gpar(col = "white"),
                       row_title = "cluster ID", 
                       cluster_rows = !order,
                       show_row_names = TRUE,
                       row_names_side = "left",
                       show_column_names = show_sample_ID,
                       top_annotation = mycol_anno)
  
  #### Significance ####
  #### in CATALYST package: row_anno <- TRUE 
  # Which one are significant?
  s <- top$p_adj <= th
  # Check for NAs
  s[is.na(s)] <- FALSE
  # Add yes or no for the significant one or not
  s <- as.matrix(c("no", "yes")[as.numeric(s) + 1])
  # Reformat the p_values
  rownames(s) <- format(top$p_adj, scientific = TRUE, digits = 3)
  # Create the heatmap with 
  row_anno <- Heatmap(matrix = s, name = "significant", 
                      col = c(no = "lightgrey", yes = "limegreen"), 
                      width = unit(5, "mm"), rect_gp = gpar(col = "white"), show_row_names = TRUE, 
                      row_names_side = "right")
  # Combine
  main <- "Differential abundance tests between clusters"
  suppressWarnings(draw(myheatmap + row_anno, column_title = main, 
                        auto_adjust = FALSE, column_title_gp = gpar(fontface = "bold", 
                                                                    fontsize = 12)))
}
