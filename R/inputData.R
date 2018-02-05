#' Read in data from a file
#'
#' This function reads in file with header
#'
#' @param filename character, file name and its path
#'
#' @usage inputData(filename)
#'
#' @return data frame
#'
#' @references
#' {TBA}
#'
#' @export
#'
inputData <- function (filename){
  if(missing(filename)){
    stop("Empty filename!")
  }
  if(!file.exists(filename)){
    stop("Could not find file!")
  }else{
    stopifnot(is_scalar_character(filename))
    indata <- read.table(filename,header=T,check.names = F)
  }
indata
}
