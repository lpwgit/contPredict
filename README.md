# R package for Sample Contamination Prediction

## Purpose
Identify sample contamination base on variant allele frequency (VAF) detected from genome sequencing on tumor samples. An example of VAF data and its configuration file are distributed with this package under inst/extdata/.

## Manual
[Refer manual here](/inst/extdata/sampleCont_vignettes.pdf)

## Installation
devtools::install_github("lpwgit/sampleCont")

## Example usage
library(sampleCont)  
data_path <-  system.file('extdata',package='sampleCont')  
output_path <- 'output_express'  
config_file <- system.file("extdata", 'config.txt', package = "sampleCont", mustWork = TRUE)  
run_sampleContamination(data_path = data_path, 
                        output_path = output_path, 
                        config_file = config_file)
