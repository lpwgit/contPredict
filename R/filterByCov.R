#' Filter mutations base on coverage
#'
#' This function filters mutation base on coverage
#'
#' @param VAFdata data frame, mutation VAF
#' @param VAFcov data frame, mutation coverage
#' @param minCov numeric, minimum coverage (default: 50)
#'
#' @usage filterByCov(VAFdata,VAFcov,minCov)
#'
#' @return data frame, mutation VAF with coverage > minCov
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
