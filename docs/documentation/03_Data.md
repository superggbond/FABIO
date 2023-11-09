---
layout: page
title: Data Input
description: ~
---
There are two major input files for FABIO to perform TWAS fine-mapping on a binary trait of interest.
#### 1. Predicted GReX matrix of the TWAS cohort
  * An example R dataframe can be loaded from the FABIO package:
 ```
    library(fabio)
    data(example)
    grex = obj[[1]]
 ```
  
#### 2. Binary phenotype of the TWAS cohort
  * An example R vector coding case as 1 and control as 0 can be loaded from the FABIO package:
 ```
    library(fabio)
    data(example)
    y = obj[[2]]
 ``` 
