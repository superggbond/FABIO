---
layout: page
title: Installation
description: ~
---

`FABIO` is implemented as a R package, which can be installed from GitHub.

### Dependencies 
* R libraries: Rcpp, RcppArmadillo, data.table

#### 1. Install `FABIO`
```
if(!require(devtools))
  install.packages(devtools)
devtools::install_github("superggbond/fabio")
library(fabio)
```
#### 2. Check the input options included in `FABIO`
```
?fabio
```
