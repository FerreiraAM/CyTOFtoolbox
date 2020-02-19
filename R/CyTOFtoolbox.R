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
#' @importFrom magrittr "%>%"
#' @importFrom ggrepel geom_text_repel
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
NULL
