#' CyTOFtoolbox: A package regrouping tools to perform CyTOF data analysis.
#'
#' The CyTOFtoolbox package regroups functions to perfom:
#' clustering, differential expression analysis.
#'
#' @section Volcano plot functions:
#' The volcano plot functions prepare the data to plot an adapted volcano plot.
#' It used the CytoGLMM output as well as the raw data.
#'
#' @docType package
#' @name CyTOFtoolbox
#'
#' @import CytoGLMM
#' @import ggplot2
#'
#' @importFrom dplyr left_join group_by
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid gpar
#' @importFrom magrittr "%>%"
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom S4Vectors metadata
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation 
NULL
