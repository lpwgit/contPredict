#' Count number of SNPs in each sample
#'
#'This function counts number of SNPs in each sample
#'
#' @param VAFdata data frame, variant allele frequency (VAF)
#' @param VAF_cutoff numeric, minimum VAF
#' @param n numeric, number of sample
#'
#' @usage numMutationperSample (VAFdata,VAF_cutoff,n)
#'
#' @return list, SNP count per sample
#'
#' @references
#' {TBA}
#'
#' @export
#'
numMutationperSample <- function(VAFdata,VAF_cutoff,n) {
  SNPcount <- t(sapply(2:(n+1), function(i) sum(VAFdata[,i]>0&(VAFdata[,i]>VAF_cutoff))))
  colnames(SNPcount) <- names(VAFdata)[2:(n+1)]
SNPcount
}
