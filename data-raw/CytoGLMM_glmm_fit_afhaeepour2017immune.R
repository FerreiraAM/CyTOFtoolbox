# Packages
library(CytoGLMM)
library(SummarizedExperiment)
library(magrittr)
library(tidyverse)

# Aghaeepour et al. (2017) dataset:
# 16 women from the Aghaeepour et al. (2017) dataset and samples from the first and third trimester of their pregnancy stimulated by IFN$\alpha$
# Path
zenodo_url <- "https://zenodo.org/record/2652578/files/"
cytof_data <- "se_aghaeepour2017immune.Rdata"
# Download data
## se_aghaeepour2017immune is an SummarizedExperiment object
load(url(paste0(zenodo_url, cytof_data)))

# Pre-processing of the data: the data are processed as in the example of the CytoGLMM workflow available online:
# Extract the values
exprs <- assay(se_aghaeepour2017immune)
# Extract the sample information
sample_info <- rowData(se_aghaeepour2017immune)
# Extract the sample names
sample_info_names <- names(sample_info)
# Create a data frame with the expression values and the sample information
df_samples <- cbind(as.data.frame(exprs), as.data.frame(sample_info))
df_samples %<>% as_tibble
# Extract protein names
protein_names <- colData(se_aghaeepour2017immune) %>%
  as.data.frame %>%
  dplyr::filter(type == "function") %>%
  .$protein_name

# Subset to NK cell population
df_samples_subset <- df_samples %>% dplyr::filter(celltype == "NK")
df_samples_subset %<>% dplyr::select(protein_names, sample_info_names)

# Transform the data
df_samples_subset_transformed <- df_samples_subset %>%
  dplyr::mutate_at(protein_names,
                   function(x) asinh(x/5))

# Fit a Generalized Linear Mixed Model (GLMM) with donor random effects.
CytoGLMM_glmm_fit_aghaeepour2017immune <- CytoGLMM::cytoglmm(df_samples_subset_transformed,
                                                    protein_names = protein_names,
                                                    condition = "term",
                                                    group = "donor",
                                                    cell_n_subsample = 1000, # Subsample is used ONLY for the sake of the data size
                                                    num_cores = 4)

# Store the data
usethis::use_data(CytoGLMM_glmm_fit_aghaeepour2017immune, overwrite = TRUE)
