#' CyTOFtoolbox: A package regrouping visualization tools for CyTOF data analysis.
#'
#' The CyTOFtoolbox package regroups visualization tools for the analysis 
#' of CyTOF data.
#'
#' @section CytoGLMM-extension:
#' The volcano plot functions prepare the data to plot a volcano plot adapted for
#' the CyTOF data. It uses the `CytoGLMM` output as well as the raw data.
#' A modified MDS plot function is also available, based on the original one from
#' the `CytoGLMM` package. It adds a threshold on the correlation marker/MDS axes.
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
NULL
