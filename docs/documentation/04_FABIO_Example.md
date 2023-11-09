---
layout: page
title: Tutorial
description: ~
---
This page provides a tutorial for TWAS fine-mapping using FABIO. Before runing the example code, make sure that the FABIO software is installed successfully. For instructions on installation, please see the [Installation section](https://superggbond.github.io/FABIO-command-line-tool/documentation/02_Installation.html).

## FABIO
The example data for FABIO tutorial can be downloaded in this [page](https://superggbond.github.io/FABIO-command-line-tool/documentation/03_Data.html). Here are the details about the input data formats and FABIO commands. 
### 1. Formats of input data for FABIO
* Predicted GReX matrix: We require the predicted GReX matrix of the TWAS cohort built up using standard softwares like [PredXican](https://github.com/hakyimlab/MetaXcan) or [BSLMM](https://github.com/genetics-statistics/GEMMA). The input GReX matrix will have gene names as the first column, with the following columns of preicted GReX at individual-level. Each following column represents GReX of genes for an individual.
* Binary phenotypes: We also require the observed binary phenotypes of the TWAS cohort. The input phenotypes should be formatted as a plain text file, coded '1' for case and '0' for control under a column named 'pheno' (see [example](https://github.com/superggbond/FABIO-command-line-tool/blob/main/data/pheno.txt)). The order of the individuals here should be consistent with the order of columns in the GReX matrix file.

### 2. Running FABIO
Before running, make sure the utility script files [FABIO_utility.R](https://github.com/superggbond/FABIO-command-line-tool/blob/main/scripts/FABIO_utility.R) and [FABIO_utility.cpp](https://github.com/superggbond/FABIO-command-line-tool/blob/main/scripts/FABIO_utility.cpp) are under the same directory as the main tool [FABIO.R](https://github.com/superggbond/FABIO-command-line-tool/blob/main/scripts/FABIO.R), and set that directory as the working directory. Then, the TWAS fine-mapping can be performed using the following command:
```
$ Rscript FABIO.R -g ../data/grex.txt -p ../data/pheno.txt -w 100 -s 1000 -o ../FABIO_PIP.txt
```
The essential inputs are:
- -g (--gene): specify the path to the predicted GReX matrix file
- -p (--pheno): specify the path to the phenotype file
- -w (--w-step): the number of warm-up steps in MCMC
- -s (--s-step): the number of sampling steps in MCMC
- -o (--out): specify the name of the output file

### 3. FABIO output
FABIO will output a summary table with two columns: 
- Gene: names of all input genes
- PIP: corresponding PIPs of all input genes
