---
layout: page
title: Data Input
description: ~
---
There are two major input files for FABIO to perform TWAS fine-mapping on a binary trait of interest. Example data can be downloaded from: https://www.dropbox.com/scl/fi/hhdjl7ep881cv8fx2hbdk/example.rda?rlkey=zanf4tyeqar49lhpu0efbpipf&dl=0
#### 1. Predicted GReX matrix of the TWAS cohort
  * An example R dataframe can be loaded from the downloaded file:
 ```r
    load('/path/to/example.rda')
    grex = inputs[[1]]
 ```
  
#### 2. Binary phenotype of the TWAS cohort
  * An example R vector coding case as 1 and control as 0 can also be loaded from the downloaded file:
 ```r
    y = inputs[[2]]
 ``` 
