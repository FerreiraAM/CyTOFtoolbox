# Goal create a volcano plot as a complementary plot for the output of the cytoGLMM package
# The cytoGLMM model(s) need to be run beforehand
# Packages: tydiverse, ggplot2.

#' Prepare data for the volcano plot
#'
#' @param data Tibble, raw CyTOF data.
#' @param protein_names Vector, vector of markers (by default = NULL, i.e. the 
#' protein names are extracted from the CytoGLMM object).
#' @param condition Character, column's name where the condition are stored.
#' @param CytoGLMM_fit Cytoglmm object, results from the CytoGLMM model.
#' @param alpha Numeric, threshold of p-value significance (by default alpha = 0.05)
#' @return Tibble.
#' @export
prepare_data_for_volcanoplot <- function(data, protein_names = NULL, condition, CytoGLMM_fit, alpha = 0.05){
  # If protein_names vector is not specify
  if (is.null(protein_names)){
    # Extract the information from the CytoGLMM fit object
    protein_names <- CytoGLMM_fit$protein_names
  } 
  # Check if the given names are in the protein_names of the object
  if (!all(protein_names %in% CytoGLMM_fit$protein_names)) {
    protein_not_in <- protein_names[which(!(protein_names %in% CytoGLMM_fit$protein_names))]
    stop(paste0(paste(protein_not_in, collapse = "; "), " not used in the CytoGLMM analysis.
                Please check the protein_names vector parameter."))
  }
  # Extract the comparison
  comparison <- function_extract_comparison(data = CytoGLMM_fit, condition = condition)
  # Filter the raw data
  data_conditions <- dplyr::pull(data, condition)
  data_filtered <- data[data_conditions %in% comparison, ]
  # Add 0.05 to the marker values (because many 1 in the data)
  data_05 <- dplyr::mutate_at(data_filtered, .vars = protein_names, .funs = function_add_05)
  
  # Compute log2 fold change
  data_05_log2foldchange <- function_compute_log2foldchange(data = data_05,
                                                            protein_names = protein_names,
                                                            condition = condition)
  # Compute MSI
  data_MSI <- function_compute_MSI(data, protein_names = protein_names)
  # Prepare CytoGLMM data
  formated_CytoGLMM_fit <- function_prepare_output_CytoGLMMmodel(data = CytoGLMM_fit, 
                                                                 alpha = alpha,
                                                                 protein_names = protein_names)
  # Combine data
  function_combine_datas(summary_CytoGLMM_fit = formated_CytoGLMM_fit,
                         data_log2foldchange = data_05_log2foldchange,
                         data_MSI = data_MSI)
  }

#' Volcano plot
#'
#' @param data Tibble, output of the \code{prepare_data_for_volcanoplot} function.
#' @return ggplot2 object.
#' @export
volcano_plot <- function(data){
  # Boundary for the MSI scale
  bound <- ceiling(log10(max(data$MSI)))
  # String of character for the x-axis
  exp_conditions <- paste0("log2 fold change (", unique(data$log2FC_ratio), ")")
  # Compute limit of the x-axis
  lim <- max(abs(data$log2foldchange))
  ## If the maximum is inferior than 1, set the limit at 1.1 to be able to draw the vertical lines
  if(lim < 1){
    lim <- 1.1
  }
  # Plot
  ggplot(data, aes_string(x = "log2foldchange", y = "log10_adjpval")) +
    geom_point(aes_string(colour = "adjpval_thres", size = "MSI")) +
    theme_bw() +
    ylab("-log10(adjusted p-value)") +
    xlab(exp_conditions) +
    xlim((0-lim), (0+lim)) +
    geom_vline(xintercept = 1, col = "gray41", linetype = "dotted", size = 1) + #fold-change less than 2 as log2(2) = 1
    geom_vline(xintercept = -1, col = "gray41", linetype = "dotted", size = 1) +
    geom_hline(yintercept = -log10(0.05), col = "gray41", linetype = "dotted", size = 1) +
    geom_vline(xintercept = 0, col = "red", size = 0.5) +
    geom_text_repel(aes_string(x = "log2foldchange", y = "log10_adjpval", label = "protein_name", colour = "adjpval_thres"), show.legend = FALSE) +
    scale_color_grey(start = 0.8, end = 0.2, name = "Adjusted p-value < 0.05") +
    scale_size_continuous(name = "MSI", breaks = c(0, 1:2 %o% 10^(0:bound)))
}

# HELPER FUNCTIONS =================================================================================

# Add 0.05 to the cyTOF values
# (in prevention of the log transformation that will be apply)
#
# @param x Numeric, vector.
# @return Numeric, vector.
function_add_05 <- function(x){
  x + 0.5
}

# Compute Mean Signal Intensity (MSI) per marker
#
# @param data Tibble, CyTOF data.
# @param protein_names Vector, vector of markers.
# @return Tibble.
function_compute_MSI <- function(data, protein_names){
  # Compute the mean for each marker
  MSI_data <- data %>%
    dplyr::select(protein_names) %>%
    dplyr::summarise_all("mean")
  # Return
  tidyr::gather(MSI_data, "protein_name", "MSI")
}

# Extract the condition used in the data (fit)
#
# @param data Cytoglmm object, results from the CytoGLMM model.
# @param condition Character, columns name where the condition are stored.
# @return Vector.
function_extract_comparison <- function(data, condition){
  two_conditions <- unique(dplyr::pull(data$df_samples_subset, condition))
  # Return
  as.character(two_conditions)
}

# Compute log2 fold change for each marker between 2 conditions
#
# @param data Tibble, CyTOF data.
# @param condition Character, columns name where the condition are stored.
# @param protein_names Vector, vector of markers.
# @return Tibble.
function_compute_log2foldchange <- function(data, condition, protein_names){
  # Compute mean per condition
  mean_per_condition <- data %>%
    dplyr::select(c(condition, protein_names)) %>%
    group_by(.dots=condition) %>%
    dplyr::summarise_all("mean")
  # log2 fold change
  mean_log2foldchange <- log2(mean_per_condition[1, protein_names]) - log2(mean_per_condition[2, protein_names])
  # Ratio
  ratio <- paste0(as.character(dplyr::pull(mean_per_condition[1, condition])), "/", 
                  as.character(dplyr::pull(mean_per_condition[2, condition])))
  # Return
  data.frame(tidyr::gather(mean_log2foldchange, "protein_name", "log2foldchange"),
             "log2FC_ratio" = ratio)
  
}

# Prepare the output from the CytoGLMM model
#
# @param data Cytoglmm object, results from the CytoGLMM model.
# @param alpha Numeric, threshold of p-value significance (by default alpha = 0.05).
# @param protein_names Vector, vector of markers.
# @return Tibble, formated CytoGLMM model results.
function_prepare_output_CytoGLMMmodel <- function(data, alpha, protein_names){
  # Extract the summary
  sum_data <- dplyr::filter(summary(data), .data$protein_name %in% protein_names)
  # Add a threshold column on the adjusted p-values
  sum_data <- dplyr::mutate(sum_data, adjpval_thres = ifelse(.data$pvalues_adj < 0.05, "significant", "non-significant"))
  # Transform the adjusted p-values with log10
  sum_data <- sum_data %>%
    dplyr::mutate(log10_adjpval = -log10(.data$pvalues_adj))
  # Return
  return(sum_data)
}


# Combine data from CytoGLMM, the log2 fold change and the MSI
#
# @param summary_CytoGLMM_fit Tibble, formated results from the CytoGLMM model.
# @param data_log2foldchange Tibble, log2 fold change per marker between 2 conditions.
# @param data_MSI Tibble, Mean Signal Intensity (MSI) per marker.
# @return Tibble.
function_combine_datas <- function(summary_CytoGLMM_fit, data_log2foldchange, data_MSI){
  # Combine log2foldchange with CytoGLMM p-values
  data_log2foldchange_CytoGLMMpvalues <- suppressMessages(left_join(summary_CytoGLMM_fit, data_log2foldchange))
  # Add MSI
  data_log2foldchange_CytoGLMMpvalues_MSI <- suppressMessages(left_join(data_log2foldchange_CytoGLMMpvalues, data_MSI))
  # Return
  return(data_log2foldchange_CytoGLMMpvalues_MSI)
}
