#' Filter mutations by coverage
#'
#'This function filters mutation by coverage
#'
#' @param VAFdata data frame, VAF
#' @param VAFcov data frame, VAF cov
#'
#' @usage filterByCov(VAFdata,VAFcov,minCov)
#'
#' @return write to output directory
#'
#' @references
#' {TBA}
#'
#' @export
#'
filterByCov <- function(VAFdata,VAFcov,minCov=50){
  ncol <- dim(VAFdata)[2]
  vaf.used <- VAFdata[,-c((ncol-1),ncol)]
  id <- which(VAFcov<minCov, arr.ind=TRUE)
  vaf.used[id[which(id[,'col']!=1),]] <- 0.00000
  filter.vaf <- cbind(vaf.used,VAFdata[,c((ncol-1),ncol)])
filter.vaf

}
