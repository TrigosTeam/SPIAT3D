
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SPIAT3D

<!-- badges: start -->

<!-- badges: end -->

The goal of SPIAT3D (**SP**atial **I**mage **A**nalysis of **T**issues
**3D**) is to comprehensively analyse 3D spatial tissue data. SPIAT3D
does this through a diverse range of spatial metrics which fully exploit
the x, y, z and cell type identity. These include 3D cell colocalization
metrics, spatial heterogeneity metrics and clustering algorithms.

## Installation

You can install the development version of SPIAT3D from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("TrigosTeam/SPIAT3D")

if (!requireNamespace("BiocManager", quietly = TRUE))     
  install.packages("BiocManager")
BiocManager::install("SpatialExperiment")
```

## Vignette

The vignette with an overview of the package can be accessed from the
top Menu under About or clicking
[here](https://trigosteam.github.io/SPIAT3D/).
