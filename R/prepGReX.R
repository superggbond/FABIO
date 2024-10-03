#'The function to help prepare GReX input file in the format that FABIO requests
#'@param chr Chromosome of the input plink files
#'@param geno_dir The directory containing input plink files
#'@param weight_dir The directory containing pre-trained eQTL effect sizes grouped by genes
#'@param save_dir The directory where the output file will be saved, set to be the working directory by default
#'
#'@return The results will be saved as a gzipped txt file 
#'@export
#'
#'@examples
#'# pre-trained GEUVADIS eQTL effect sizes using BSLMM can be downloaded from: 
#'# https://www.dropbox.com/scl/fo/fxynm8uvedgvy7ni6hcbt/AAfTQVo89s78DsRNwpBH3lU?rlkey=nbqwrdi2r5y1bbojzf7z8ev7h&st=yz28n4nj&dl=0
#'# plink2 is required in the working environment
#'
#'library(FABIO)
#'
#'chr <- 22
#'geno_dir <- "/path/to/plink/files"
#'weight_dir <- "./BSLMM_weights/chr22"
#'prepGReX(chr, geno_dir, weight_dir)
#'# The results will be saved as a file named "grex_for_fabio.txt.gz"

prepGReX = function(chr, geno_dir, weight_dir, out_dir='.') {
  output = list()
  weight_list = list.files(weight_dir)
  system('mkdir ./temp_files')
  
  numCores = ifelse(detectCores()>8,8,detectCores())
  results = mclapply(c(1:length(weight_list)), prepGReX_gene, mc.cores = numCores)
  
  output = do.call("rbind", results)
  output = as.data.frame(output)
  inds = 2:ncol(output)
  output[,inds] = apply(output[,inds], 2, function(x) as.numeric(x))
  fwrite(output,paste0(out_dir,'/grex_for_fabio.txt.gz'),col.names=F,sep=' ')
  
  system('rm -r temp_files')
}
