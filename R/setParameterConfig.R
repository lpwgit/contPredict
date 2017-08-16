#' Set parameters base on configuration file
#'
#'This function takes in data from configuration file to set parameters
#'
#' @param filename character, configuration file name and its path
#'
#' @usage setParameterConfig(filename)
#'
#' @references
#' {TBA}
#'
#' @export
#'
setParameterConfig <- function (configFile){
  if(missing(configFile)){
    stop("Empty file name!")
  }
  if(!file.exists(configFile)){
    stop(paste("Could not find",configFile))
  }else{
    config.data <- as.data.frame(read.csv(configFile,header=F,sep='\t'))
    ## set parameter
    VAF_cutoff1 <<- as.numeric(as.character(config.data[1,2])) # relation determination default: 0.05
    VAF_cutoff <<-  as.numeric(as.character(config.data[2,2])) # min VAF to be considered present in a particular sample, default: 0.002
    VAF_ignore <<-  as.numeric(as.character(config.data[3,2])) # variants inside the box, default: 0.2
    localPcomm_cutoff <<- as.numeric(as.character(config.data[4,2])) # classify high/low contamination cases, default: 10% as high otherwise low contamination pairs
    p.val_cutoff <<-  as.numeric(as.character(config.data[5,2])) #fisher test threshold for multisource contamination
    num_round_digit <<- as.numeric(as.character(config.data[6,2])) #rounding num digit

    ## relation determination
    center_cutoff <<-   as.numeric(as.character(config.data[7,2])) # globalPcomm > globalPcomm_cutoff, same patient when propCommC2 > center_cutoff (default: 50%)
    source_cutoff <<-  as.numeric(as.character(config.data[8,2]))  # globalPcomm > globalPcomm_cutoff, one way when propCommC1 > source_cutoff && propCommC1 < target_cutoff, x->y vise versa (default: 40%)
    target_cutoff <<-   as.numeric(as.character(config.data[9,2])) # globalPcomm > globalPcomm_cutoff, one way when propCommC1 > source_cutoff && propCommC1 < target_cutoff, x->y vise versa (default: 10%)
    region_cutoff <<-   as.numeric(as.character(config.data[10,2]))# globalPcomm < globalPcomm_cutoff, (default: 75%)
    n_sample <<-  as.numeric(as.character(config.data[11,2]))
  }
}
