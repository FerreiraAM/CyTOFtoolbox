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

Note: The installation of the CytoGLMM package is needed.
```
devtools::install_github("ChristofSeiler/CytoGLMM")
```

## Summary

In this package, multiple functions have been created to complement the existing workflow:

### Vocano plot

- Volcano plot as a complementary plot to the cytoglmm analysis.
- To access the vignette: 
```
vignette("volcano-plot", package = "CyTOFtoolbox")
```

### Clustering

- Modified heatmap plot of the differential abundance (DA) test results from 
`CATALYST` package.
- To access the vignette:
```
vignette("clustering", package = "CyTOFtoolbox")
```