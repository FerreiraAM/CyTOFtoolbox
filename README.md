# CyTOFtoolbox

An R package regrouping tools to perfom CyTOF data analysis.

## Installation

The R-package itself can be installed from source by using the devtools package. 
It also allows building vignettes.

```
devtools::install_github("FerreiraAM-stanford/CyTOFtoolbox",
  build_vignettes = TRUE)
```

Note: The installation of the CytoGLMM package is mandatory.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CytoGLMM")
```

## Summary

In this package, functions were created to complement the existing 
workflow:

### CytoGLMM extension

- Volcano plot as a complementary plot to the cytoglmm analysis.
- MDS plot with a threshold on the correlation marker/MDS axes.

To access the vignette: 

```
vignette("CytoGLMM-extension", package = "CyTOFtoolbox")
```

## Notes

Please be aware that some `CytoGLMM` workflows shared in the laboratory use the 
`ceiling` function while reading the FCS files.
