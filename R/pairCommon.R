#' Pairwise Sample Function
#'
#' @description
#' 4 different functions for pairwise sample
#'
#' @details
#' pairPCommon calculates pcommon of pairwise sample
#'
#' pairShare counts number of common SNPs for pairwise sample
#'
#' pairList generates all pairwise sample list
#'
#' pairCommList tabulates pairwise pcommon for each paired sample in a row
#'
#' @param VAFdata data frame, variant allele frequency (VAF)
#' @param VAF_cutoff numeric, minimum VAF
#' @param n_sample numeric, number of sample
#' @param num_round_digit numeric, number of rounding decimal digit (default:3)
#' @param SNPcomm matrix, output from the calling of pairPCommon(VAFdata,VAF_cutoff,num_round_digit,n_sample)
#' @param del character, file delimiter
#'
#' @usage
#' pairPCommon(VAFdata,VAF_cutoff,num_round_digit,n_sample)
#' pairShare(VAFdata,VAF_cutoff,n_sample)
#' pairList(VAFdata,n_sample,del)
#' pairCommList(SNPcomm,n_sample)
#'
#'
#' @return
#' pairPCommon: matrix, nxn pcommon of pairwise sample, n: number of sample
#'
#' pairShare: matrix, nxn number of common snps of pairwise sample
#'
#' pairList: vector, all possible pairwise samples
#'
#' pairCommList: list, all possible pairwise samples and the respective pcommon values
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
pairPCommon<- function(VAFdata,VAF_cutoff,num_round_digit,n_sample) {
  SNPcomm1 <- t(sapply(2:(n_sample+1), function(i) sapply(2:(n_sample+1), function(j) sum(VAFdata[,i]*VAFdata[,j]>0&(VAFdata[,i]>VAF_cutoff|VAFdata[,j]>VAF_cutoff))/sum(VAFdata[,i]>0&(VAFdata[,i]>VAF_cutoff|VAFdata[,j]>VAF_cutoff)))))
  colnames(SNPcomm1) <- names(VAFdata)[2:(n_sample+1)]
  rownames(SNPcomm1) <- names(VAFdata)[2:(n_sample+1)]
  SNPcomm_d <- round(SNPcomm1, digits = num_round_digit)
  SNPcomm_d
}

## pairwise common SNP count
pairShare<- function(VAFdata,VAF_cutoff,n_sample) {
SNPshare <- t(sapply(2:(n_sample+1), function(i) sapply(2:(n_sample+1), function(j) sum(VAFdata[,i]*VAFdata[,j]>0&(VAFdata[,i]>VAF_cutoff|VAFdata[,j]>VAF_cutoff)))))
colnames(SNPshare) <- names(VAFdata)[2:(n_sample+1)]
rownames(SNPshare) <- names(VAFdata)[2:(n_sample+1)]
SNPshare
}

## all pairwise sample list
pairList <- function(VAFdata,n_sample,del){
a <- sapply(1:(n_sample-1), function(i) sapply((i+1):n_sample, function(j) {
  paste(colnames(VAFdata)[i+1],colnames(VAFdata)[j+1], sep = del)
}))
unlist(a)
}

## all pairwise pcommon
pairCommList <- function (SNPcomm,n_sample){
b <- sapply(1:(n_sample-1), function(i) t(sapply((i+1):n_sample, function(j) {
  c(SNPcomm[i,j],SNPcomm[j,i])
})))
b
}


