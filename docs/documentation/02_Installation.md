---
layout: page
title: Installation
description: ~
---

`FABIO` is implemented as a command-line tool, which can be installed from GitHub.

### Dependencies 
* R libraries: Rcpp, RcppArmadillo, data.table, and optparse
* When running the main command-line tool 'FABIO.R', it will automatically check and install all required R packages.

#### 1. Install `FABIO`
```
git clone https://github.com/superggbond/FABIO-command-line-tool.git
```
#### 2. Check the input options included in `FABIO`
```
cd ./FABIO-command-line-tool/scripts
Rscript ./FABIO.R -h
```
