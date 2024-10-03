---
layout: page
title: Data Input
description: ~
---
There are two major input files for FABIO to perform TWAS fine-mapping on a binary trait of interest. Example data can be downloaded [here](https://www.dropbox.com/scl/fo/fxynm8uvedgvy7ni6hcbt/AAfTQVo89s78DsRNwpBH3lU?rlkey=nbqwrdi2r5y1bbojzf7z8ev7h&st=yz28n4nj&dl=0).
#### 1. Predicted GReX of the TWAS cohort
  * An example R dataframe can be loaded from the downloaded file:
 ```r
    grex <- data.table::fread('./example_grex.txt.gz')
 ```
  
#### 2. Binary phenotype of the TWAS cohort
  * An example R vector coding case as 1 and control as 0 can also be loaded from the downloaded file:
 ```r
    pheno <- scan('./example_pheno.txt',numeric())
 ``` 
