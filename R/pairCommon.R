#' Pairwise sample function
#'
#' @description
#' 4 different functions for pairwise sample
#'
#' @details
#' pairPCommon: calculates pcommon
#'
#' pairShare: counts number of common SNPs
#'
#' pairList: generates all pairwise sample list
#'
#' pairCommList: tabulates pcommon for all pairwise samples
#'
#' @param VAFdata data frame, mutation VAF
#' @param VAF_cutoff numeric, minimum VAF
#' @param n numeric, number of sample
#' @param num_round_digit numeric, number of rounding decimal digit (default: 3)
#' @param SNPcomm matrix, output from the calling of pairPCommon(VAFdata,VAF_cutoff,num_round_digit,n)
#' @param del character, file delimiter
#'
#' @usage
#' pairPCommon(VAFdata,VAF_cutoff,num_round_digit,n)
#' pairShare(VAFdata,VAF_cutoff,n)
#' pairList(VAFdata,n,del)
#' pairCommList(SNPcomm,n)
#'
#'
#' @return
#' pairPCommon: matrix, nxn pcommon, n: number of sample
#'
#' pairShare: matrix, nxn number of common snps
#'
#' pairList: vector, all pairwise samples
#'
#' pairCommList: list, all pairwise samples and their pcommon values
#'
#'
#' @references
#' {TBA}
#'
#' @export pairPCommon
#' @export pairShare
#' @export pairList
#' @export pairCommList
#'
pairPCommon<- function(VAFdata,VAF_cutoff,num_round_digit,n) {
  SNPcomm1 <- t(sapply(2:(n+1), function(i) sapply(2:(n+1), function(j) sum(VAFdata[,i]*VAFdata[,j]>0&(VAFdata[,i]>VAF_cutoff|VAFdata[,j]>VAF_cutoff))/sum(VAFdata[,i]>0&(VAFdata[,i]>VAF_cutoff|VAFdata[,j]>VAF_cutoff)))))
  colnames(SNPcomm1) <- names(VAFdata)[2:(n+1)]
  rownames(SNPcomm1) <- names(VAFdata)[2:(n+1)]
  SNPcomm_d <- round(SNPcomm1, digits = num_round_digit)
  SNPcomm_d
}

## pairwise common SNP count
pairShare<- function(VAFdata,VAF_cutoff,n) {
SNPshare <- t(sapply(2:(n+1), function(i) sapply(2:(n+1), function(j) sum(VAFdata[,i]*VAFdata[,j]>0&(VAFdata[,i]>VAF_cutoff|VAFdata[,j]>VAF_cutoff)))))
colnames(SNPshare) <- names(VAFdata)[2:(n+1)]
rownames(SNPshare) <- names(VAFdata)[2:(n+1)]
SNPshare
}

## all pairwise sample list
pairList <- function(VAFdata,n,del){
a <- sapply(1:(n-1), function(i) sapply((i+1):n, function(j) {
  paste(colnames(VAFdata)[i+1],colnames(VAFdata)[j+1], sep = del)
}))
unlist(a)
}

## all pairwise pcommon
pairCommList <- function (SNPcomm,n){
b <- sapply(1:(n-1), function(i) t(sapply((i+1):n, function(j) {
  c(SNPcomm[i,j],SNPcomm[j,i])
})))
b
}


