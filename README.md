# CyTOFtoolbox

An R package regrouping tools to perfom CyTOF data analysis.

# Installation

The R-package itself can be installed from source by using the devtools package. This also allows
building vignettes.
```
# The package is private, so be sure to use the token
devtools::install_github("FerreiraAM-stanford/CyTOFtoolbox",
  auth_token = "1da45a4220297a9c55caf0d64794a87b3cd10c3a", 
  build_vignettes = TRUE)
```

Note: Be sure to have install the CytoGLMM package.
```
devtools::install_github("ChristofSeiler/CytoGLMM")
```

# Summary

In this package, there are the functions to create the volcano plot (complementary plot to the cytoglmm analysis).
