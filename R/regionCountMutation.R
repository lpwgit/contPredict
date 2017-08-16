#'Count number of mutation in 7 regions of variant allele frequency (VAF) scatter plot for a pair of sample (s1,s2)
#'
#'This function counts number of SNPs in 7 regions of pairwise sample VAF scatter plot
#'
#'@param sample_pairs data frame
#'@param VAFdata data frame, VAF
#'@param SNPcount integer matrix, number of SNPs for each sample
#'@param SNPshare integer matrix, number of common SNPs
#'@param VAF_cutoff numeric, minimum VAF (default: 0.002)
#'@param VAF_ignore numeric, ignore variant below this cutoff (default: 0.2)
#'@param n integer, number of sample
#'
#' @usage regionCountMutation(VAFdata,SNPcount,SNPshare,VAF_cutoff,VAF_ignore,n)
#'
#' @return matrix, n SNPs x 10 columns: number of SNPs in 7 regions, number of SNPs in s1, number of SNPs in s2, common SNPs between s1 and s2
#'
#' @references
#' {TBA}
#'
#' @export
#'
regionCountMutation <- function(sample_pairs,VAFdata,SNPcount,SNPshare,VAF_cutoff=0.002,VAF_ignore=0.2,n) {

  countPoint <- apply(sample_pairs, 1, function(i) {
    r1=r2=r3=r4=r5=r6=r7=0

    r1 <- sum(VAFdata[,i["pid1"]]==0 & VAFdata[,i["pid2"]]==0) #x=y=0
    r2 <- sum(VAFdata[,i["pid1"]]<=VAF_ignore & VAFdata[,i["pid2"]]<=VAF_ignore & !(VAFdata[,i["pid1"]]==0 & VAFdata[,i["pid2"]]==0)  )   #0<x&y<=0.2, corner>0
    r3 <- sum(VAFdata[,i["pid1"]]>VAF_ignore & VAFdata[,i["pid2"]]<=VAF_cutoff) ## X only
    r4 <- sum(VAFdata[,i["pid2"]]>VAF_ignore & VAFdata[,i["pid1"]]<=VAF_cutoff) ## Y only

    r5 <- sum(VAFdata[,i["pid1"]]>VAF_ignore &  VAFdata[,i["pid2"]] > VAF_cutoff & VAFdata[,i["pid2"]]<0.5*VAFdata[,i["pid1"]]) ## x>>y
    r6 <- sum(VAFdata[,i["pid2"]]>VAF_ignore &  VAFdata[,i["pid1"]] > VAF_cutoff & VAFdata[,i["pid2"]]>2*VAFdata[,i["pid1"]]) ## y>>x

    r7 <- sum(VAFdata[,i["pid2"]]<=2*VAFdata[,i["pid1"]] & VAFdata[,i["pid2"]]>=0.5*VAFdata[,i[["pid1"]]] & !(VAFdata[,i["pid2"]]<=VAF_ignore & VAFdata[,i["pid1"]]<=VAF_ignore)) ## center

    c1 <- SNPcount[,i["pid1"]]
    c2 <- SNPcount[,i["pid2"]]
    share <- SNPshare[i["pid1"],i["pid2"]] #all common SNP btw pid1&pid2

    c(r1,r2,r3,r4,r5,r6,r7,c1,c2,share)
  })
  countPoint <- t(countPoint)
  colnames(countPoint)=c("x0y0","xylt0d2","xonly","yonly","x2y","y2x","center","ttl_pid1","ttl_pid2","share")
  countPoint
}
