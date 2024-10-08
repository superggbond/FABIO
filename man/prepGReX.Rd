% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepGReX.R
\name{prepGReX}
\alias{prepGReX}
\title{The function to help prepare GReX input file in the format that FABIO requests}
\usage{
prepGReX(chr, geno_dir, weight_dir, out_dir = ".")
}
\arguments{
\item{chr}{Chromosome of the input plink files}

\item{geno_dir}{The directory containing input plink files}

\item{weight_dir}{The directory containing pre-trained eQTL effect sizes grouped by genes}

\item{save_dir}{The directory where the output file will be saved, set to be the working directory by default}
}
\value{
The results will be saved as a gzipped txt file
}
\description{
The function to help prepare GReX input file in the format that FABIO requests
}
\examples{
# pre-trained GEUVADIS eQTL effect sizes using BSLMM can be downloaded from: 
# https://www.dropbox.com/scl/fo/fxynm8uvedgvy7ni6hcbt/AAfTQVo89s78DsRNwpBH3lU?rlkey=nbqwrdi2r5y1bbojzf7z8ev7h&st=yz28n4nj&dl=0
# plink2 is required in the working environment

library(FABIO)

chr <- 22
geno_dir <- "/path/to/plink/files"
weight_dir <- "./BSLMM_weights/chr22"
prepGReX(chr, geno_dir, weight_dir)
# The results will be saved as a file named "grex_for_fabio.txt.gz"
}
