
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- [![DOI](https://zenodo.org/badge/197059951.svg)](https://zenodo.org/badge/latestdoi/197059951) -->

# cdx2cea <img src='man/figures/logo.png' align="right" height="139" />

<!-- <img src="docs/figs/under_const.jpeg" align="center" alt="" width="360" /> -->

[`cdx2cea`](https://github.com/feralaes/cdx2cea) is an R package that
implements the cost-effectiveness analysis (CEA) of the of testin
avereag-risk Stage II colon cancer patients for the absence of CDX2
biomarker followed by adjuvant chemotherapy. The main website of
`cdx2cea` can be [found here](https://darth-git.github.io/cdx2cea/).

[`cdx2cea`](https://github.com/DARTH-git/cdx2cea) is part of the
following manuscript:

-   Alarid-Escudero F, Schrag D, Kuntz KM. “CDX2 biomarker testing and
    adjuvant therapy for stage II colon cancer: An exploratory
    cost-effectiveness analysis” (Under review)

<!-- The release that accompanies the published article has been archived in zenodo: https://zenodo.org/record/3445451. -->

# Preliminaries

-   Install
    [RStudio](https://www.rstudio.com/products/rstudio/download/)
-   Install `devtools` to install `darthpack` as a package and modify it
    to generate your own package

``` r
# Install release version from CRAN
install.packages("devtools")

# Or install development version from GitHub
# devtools::install_github("r-lib/devtools")
```

# Usage and installation

`cdx2cea` repository could be used in two different ways:

1.  [Regular coding
    template](#use-repository-as-a-regular-coding-template) for using it
    to generate a repository of your own model-based decision or
    cost-effectiveness analysis
2.  [R package](#use-as-an-r-package) for using it as a standalone
    package to run current functions of `darthpack`

The main website of the package could be found in:
<https://darth-git.github.io/darthpack/>

## Use repository as a regular coding template

1.  On the `cdx2cea` GitHub repository, navigate to the main page of the
    repository (<https://github.com/feralaes/cdx2cea>).
2.  Above the file list, click **Clone or download** and select either
    1.  **Open in desktop**, which requires the user to have a GitHub
        desktop installed, or
    2.  **Download zip** that will ask the user to download the whole
        repository as a .zip file.
3.  Open the RStudio project `cdx2cea.Rproj`.
4.  Install all the required and suggested packages listed in the
    [*DESCRIPTION*](https://github.com/feralaes/cdx2cea/blob/master/DESCRIPTION)
    file in the main folder of the repository
    -   To install `dampack`, please follow these instructions:

``` r
# Install development version from GitHub
devtools::install_github("feralaes/cdx2cea")
```

1.  In RStudio, load all the functions and data from the repository by
    typing `devtools::load_all(".")`
2.  Run all the decision modeling modules in the analysis folder.

## Use as an R package

1.  Install the development version of `cdx2cea` from
    [GitHub](https://github.com) with:

``` r
devtools::install_github("feralaes/cdx2cea")
```

1.  Load all the functions and data from the repository by typing

``` r
library(cdx2cea)
```
