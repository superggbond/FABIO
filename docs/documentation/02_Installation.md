---
layout: page
title: Installation
description: ~
---

`FABIO` is implemented as a R package, which can be installed from GitHub.

### Dependencies 
* R libraries: Rcpp, RcppArmadillo, RcppProgress, data.table, snpStats, parallel, MASS
* PLINK is required, if you'd like to use the helper function to generate GReX input file using PLINK 1 binary files (.bed+.bim+.fam) of genotypes and pre-trained eQTL weights 

#### 1. Install `FABIO`
```r
install.packages("devtools")
devtools::install_github("superggbond/FABIO")
library(FABIO)
```
#### 2. Check the input options for the main function to perform TWAS fine-mapping on binary outcomes
```r
?fabio
```
#### 3. Check the input options for the helper function to prepare GReX input file in the format that FABIO requests
```r
?prepGReX
```
