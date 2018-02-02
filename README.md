# R package for Sample Contamination Prediction

## Purpose
Identify sample contamination base on variant allele frequency (VAF) detected from genome sequencing on tumor samples. An example of VAF data and its configuration file are distributed with this package under inst/extdata/.

## Manual
<<<<<<< HEAD
[Refer manual here](/vignettes/sampleCont_manual.pdf)
=======
[Refer manual here](/inst/extdata/sampleCont_vignettes.pdf)
>>>>>>> 3b82f69fbb478a4e229c61fd993f97a6b55130cd

## Installation
devtools::install_github("lpwgit/sampleCont")

## Example usage
<<<<<<< HEAD
library(sampleCont)
data_path <- past0(getwd(),'/inst/extdata')
output_path <- past0(getwd() ,'/output')
config_file <- system.file("extdata", 'config.txt',
    package = "sampleCont", mustWork = TRUE)
run_sampleContamination(data_path = data_path,output_path = output_path, 
  config_file=config_file,rmcov_cutoff=0,manualsetPar=FALSE)
=======
library(sampleCont)  
data_path <-  system.file('extdata',package='sampleCont')  
output_path <- 'output_express'  
config_file <- system.file("extdata", 'config.txt', package = "sampleCont", mustWork = TRUE)  
run_sampleContamination(data_path = data_path, 
                        output_path = output_path, 
                        config_file = config_file)
>>>>>>> 3b82f69fbb478a4e229c61fd993f97a6b55130cd
