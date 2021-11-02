
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://zenodo.org/badge/331070175.svg)](https://zenodo.org/badge/latestdoi/331070175)

# cdx2cea <img src='man/figures/logo.png' align="right" height="139" />

<!-- <img src="docs/figs/under_const.jpeg" align="center" alt="" width="360" /> -->

[`cdx2cea`](https://github.com/feralaes/cdx2cea) is an R package that
implements the cost-effectiveness analysis (CEA) of testing average-risk
Stage II colon cancer patients for the absence of CDX2 biomarker
expression followed by adjuvant chemotherapy.
<!-- The main website of `cdx2cea` can be [found here](https://feralaes.github.io/cdx2cea/). -->

[`cdx2cea`](https://github.com/feralaes/cdx2cea) is part of the
following manuscript:

-   Alarid-Escudero F, Schrag D, Kuntz KM. (2021) [“CDX2 biomarker
    testing and adjuvant therapy for stage II colon cancer: An
    exploratory cost-effectiveness
    analysis”](https://www.sciencedirect.com/science/article/pii/S1098301521017472)
    *Value in Health* (Online First).

The release that accompanies the published article has been archived in
zenodo: <https://zenodo.org/record/5093594#.YPYyDy1h1qs>

## How to cite this package in your article

You can cite this package like this “we based our analysis using the
cdx2cea R package (Alarid-Escudero F, Schrag D, and Kuntz KM 2021)”.
Here is the full bibliographic reference to include in your reference
list for the manuscript and the package (don’t forget to update the
‘last accessed’ date):

> Alarid-Escudero F, Schrag D, Kuntz KM (2021). [“CDX2 biomarker testing
> and adjuvant therapy for stage II colon cancer: An exploratory
> cost-effectiveness
> analysis”](https://www.sciencedirect.com/science/article/pii/S1098301521017472).
> *Value in Health* (In press).

> Alarid-Escudero F, Schrag D, Kuntz KM (2021). {cdx2cea}: A
> cost-efectiveness analysis of testing stage II colon cancer patients
> for the absence of CDX2 biomarker followed by adjuvant chemotherapy
> (Version v1.0.0). Zenodo.
> [10.5281/zenodo.5093594](https://www.doi.org/10.5281/zenodo.5093594).
> Last accessed 12 July 2021

## Preliminaries

-   Install
    [RStudio](https://www.rstudio.com/products/rstudio/download/)
-   Install `devtools` to install `cdx2cea` as a package and modify it
    to generate your own package

``` r
# Install release version from CRAN
install.packages("devtools")

# Or install development version from GitHub
# devtools::install_github("r-lib/devtools")
```

## Usage and installation

`cdx2cea` repository could be used in two different ways:

1.  [Regular coding
    template](#use-repository-as-a-regular-coding-template) for using it
    to generate a repository of your own model-based decision or
    cost-effectiveness analysis
2.  [R package](#use-as-an-r-package) for using it as a standalone
    package to run current functions of `cdx2cea`

<!-- The main website of the package could be found in: https://feralaes.github.io/feralaes/ -->

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
    -   To install `cdx2cea`, please follow these instructions:

``` r
# Install development version from GitHub
devtools::install_github("feralaes/cdx2cea")
```

5.  In RStudio, load all the functions and data from the repository by
    typing `devtools::load_all(".")`
6.  Run all the decision modeling modules in the analysis folder.

## Use as an R package

1.  Install the development version of `cdx2cea` from
    [GitHub](https://github.com) with:

``` r
devtools::install_github("feralaes/cdx2cea")
```

2.  Load all the functions and data from the repository by typing

``` r
library(cdx2cea)
```

## Citation

Alarid-Escudero F, Schrag D, Kuntz KM (2021). cdx2cea: A
cost-efectiveness analysis of testing stage II colon cancer patients for
the absence of CDX2 biomarker followed by adjuvant chemotherapy (Version
v1.0.0). Zenodo. <http://doi.org/10.5281/zenodo.5093594>

## Acknowledgements

This work was supported by a grant from Fulbright and the National
Council of Science and Technology of Mexico (CONACYT) and a Doctoral
Dissertation Fellowship from the Graduate School of the University of
Minnesota as part of Dr. Alarid-Escudero’s doctoral program. Drs. Kuntz
and Alarid-Escudero were supported by two grants from the National
Cancer Institute at the National Institutes of Health (grant numbers
U01-CA-199335 and U01-CA-253913) as part of the Cancer Intervention and
Surveillance Modeling Network (CISNET). The funding agencies had no role
in the design of the study, interpretation of results, or writing of the
manuscript. The content is solely the responsibility of the authors and
does not necessarily represent the official views of the National
Institutes of Health. The funding agreement ensured the authors’
independence in designing the study, interpreting the data, writing, and
publishing the report. No other funding noted.
