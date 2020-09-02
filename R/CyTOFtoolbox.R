#' CyTOFtoolbox: A package regrouping tools to perform CyTOF data analysis.
#'
#' The CyTOFtoolbox package regroups functions to perfom:
#' differential expression analysis and clustering.
#'
#' @section CytoGLMM-extension:
#' The volcano plot functions prepare the data to plot a volcano plot adapted for
#' the CyTOF data. It uses the `CytoGLMM` output as well as the raw data.
#' A modified MDS plot function is also available, based on the original one from
#' the `CytoGLMM` package. It adds a threshold on the correlation marker/MDS axes.
#' 
#' @section Clustering:
#' A modified heatmap plot of the differential abundance (DA) test results from 
#' `CATALYST` package is available with the option of showing or not the sample
#' IDs as well as predefining the colors if needed.
#'
#' @docType package
#' @name CyTOFtoolbox
#' 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
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
NULL
