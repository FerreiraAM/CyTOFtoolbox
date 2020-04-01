# CyTOFtoolbox

An R package regrouping tools to perfom CyTOF data analysis.

## Installation

The R-package itself can be installed from source by using the devtools package. This also allows
building vignettes.
The repository is private, so be sure to use the token.

```
devtools::install_github("FerreiraAM-stanford/CyTOFtoolbox",
  auth_token = "1da45a4220297a9c55caf0d64794a87b3cd10c3a", 
  build_vignettes = TRUE)
```

Note: The installation of the CytoGLMM package is mandatory.
```
devtools::install_github("ChristofSeiler/CytoGLMM")
```

## Summary

In this package, multiple functions have been created to complement the existing workflow:

### CytoGLMM extension

- Volcano plot as a complementary plot to the cytoglmm analysis.
- MDS plot with a threshold on the correlation marker/MDS axes.
- To access the vignette: 
```
vignette("CytoGLMM-extension", package = "CyTOFtoolbox")
```

### Clustering

- Modified heatmap plot of the differential abundance (DA) test results from 
`CATALYST` package.
- To access the vignette:
```
vignette("clustering", package = "CyTOFtoolbox")
```