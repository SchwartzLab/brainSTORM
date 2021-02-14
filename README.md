
# brainSTORM

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.0.1-blue.svg)](https://github.com/SchwartzLab/brainSTORM)
<!-- badges: end -->

## Description

**brainSTORM** is a package that processes RNA-seq reads alignments into
transcriptomic-oriented objects, making use of the
[**txtools**](https://github.com/AngelCampos/txtools) package, to
process BRAINseq-based data.

## Quick example

## Installation

You can install the development version from
[GitHub](https://github.com/SchwartzLab/brainSTORM) typing in the
following commands in the R console and changing the string `YOUR_TOKEN`
for a Personal Access Token generated via GitHub:

  - Create an authentication token in Github following this
    [link](https://github.com/settings/tokens).
  - For “Scopes” select all under the “Repo: Full control of private
    repositories” category

<!-- end list -->

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("SchwartzLab/brainSTORM", auth_token = "YOUR_TOKEN")
```

## Further documentation

## Current limitations

## Licence

brainSTORM is developed under the [Artistic
Licence 2.0](https://opensource.org/licenses/Artistic-2.0).
