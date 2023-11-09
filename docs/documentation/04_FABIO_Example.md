---
layout: page
title: Tutorial
description: ~
---
This page provides a tutorial for TWAS fine-mapping using FABIO. Before runing the example code, make sure that the FABIO software is installed successfully. For instructions on installation, please see the [Installation section](https://superggbond.github.io/fabio/documentation/02_Installation.html).

## FABIO
The example data for FABIO tutorial can be accessed following this [page](https://superggbond.github.io/fabio/documentation/03_Data.html). Here are the details about the input data formats and how to run FABIO. 
### 1. Formats of input data for FABIO
* Predicted GReX matrix: We require the predicted GReX matrix of the TWAS cohort built up using standard softwares like [PredXican](https://github.com/hakyimlab/MetaXcan) or [BSLMM](https://github.com/genetics-statistics/GEMMA). The input GReX matrix will have gene names as the first column, with the following columns of preicted GReX at individual-level. Each following column represents GReX of genes for an individual.
* Binary phenotypes: We also require the observed binary phenotypes of the TWAS cohort. The input phenotypes should be formatted as a numeric vector, coded 1 for case and 0 for control. The order of the individuals here should be consistent with the order of columns in the GReX matrix input.

### 2. Running FABIO
The TWAS fine-mapping can be performed using the following scripts with our example data:
```r
library(fabio)
data(example)
grex <- obj[[1]]
y <- obj[[2]]
output <- fabio(grex, y, 100, 1000)
```
The essential inputs are:
- grex: the predicted GReX matrix dataframe
- y: the TWAS phenotype vector
- -w (--w-step): the number of warm-up steps in MCMC
- -s (--s-step): the number of sampling steps in MCMC
- -o (--out): specify the name of the output file

### 3. FABIO output
FABIO will output a summary table with two columns: 
- Gene: names of all input genes
- PIP: corresponding PIPs of all input genes
