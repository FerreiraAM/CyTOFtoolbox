# CyTOFtoolbox

An R package regrouping tools to perfom CyTOF data analysis.

# Installation

The R-package itself can be installed from source by using the devtools package. This also allows
building vignettes.
```
# minimal install without vignettes
devtools::install_github("FerreiraAM-stanford/CyTOFtoolbox")

# install with building vignettes and also installing suggested dependencies
devtools::install_github("FerreiraAM-stanford/CyTOFtoolbox", build_vignettes = TRUE, dependencies = TRUE)
```

# Summary

In this package, there are the functions to create the volcano plot (complementary plot to the cytoglmm analysis).
