<!-- README.md is generated from README.Rmd. Please edit that file -->

# ortest <img src="man/figures/package-sticker.png" align="right" style="float:right; height:120px;"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/ortest)](https://CRAN.R-project.org/package=ortest)
[![R CMD
Check](https://github.com/chaoqiw0324/ortest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chaoqiw0324/ortest/actions/workflows/R-CMD-check.yaml)
[![Website](https://github.com/chaoqiw0324/ortest/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/chaoqiw0324/ortest/actions/workflows/pkgdown.yaml)
[![Test
coverage](https://github.com/chaoqiw0324/ortest/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/chaoqiw0324/ortest/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/chaoqiw0324/ortest/branch/master/graph/badge.svg)](https://codecov.io/gh/chaoqiw0324/ortest)
[![License: GPL (&gt;=
2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
<!-- badges: end -->

<p align="left">
• <a href="#overview">Overview</a><br> •
<a href="#features">Features</a><br> •
<a href="#installation">Installation</a><br> •
<a href="#get-started">Get started</a><br> •
<a href="#long-form-documentations">Long-form documentations</a><br> •
<a href="#citation">Citation</a><br> •
<a href="#contributing">Contributing</a><br> •
<a href="#acknowledgments">Acknowledgments</a><br> •
<a href="#references">References</a>
</p>

## Overview

The R package `ortest` provides an odds ratio based conditional independence test for mixed data types, using nonparametric methods for nuisance estimation, including random forests, XGboost, and SuperLearner

## Features

The main purpose of `ortest` is to test conditional independence test for mixed data types including numeric, binary and multi-level factor, and to implemented in pc algorithm. 

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

    ## Install < remotes > package (if not already installed) ----
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }

    ## Install < ortest > from GitHub ----
    remotes::install_github("chaoqiw0324/ortest")

Then you can attach the package `ortest`:

    library("ortest")

## Get started

For an overview of the main features of `ortest`, please read the [Using the Ortest package
](https://chaoqiw0324.github.io/ortest/docs/articles/ortest.html)
vignette.



## Citation

Please cite `ortest` as:

> Wu Chaoqi (2024) ortest: An R package to ortest. R package
> version 0.0.0.9000. <https://github.com/chaoqiw0324/ortest/>

## Contributing

All types of contributions are encouraged and valued. For more
information, check out our [Contributor
Guidelines](https://github.com/chaoqiw0324/ortest/blob/main/CONTRIBUTING.md).

Please note that the `ortest` project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.


