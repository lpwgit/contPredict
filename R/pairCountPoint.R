#' Count number of SNPs in 7 regions of a given pair of sample
#'
#'This function counts number of SNPs in 7 regions of a given pair of sample
#'
#' @param VAFdata data frame, Variant Allele Frequency
#' @param pid1 character, first sample id
#' @param pid2 character, second sample id
#' @param slope1 numeric, first slope that is closer to axis-X (default:0.5)
#' @param slope2 numeric, second slope that is closer to axis-Y (default:2)
#'
#' @usage pairCountPoint(VAFdata,pid1,pid2,slope1,slope2)
#'
#' @return matrix of 7 columns: number of SNPs in 7 regions
#'
#' @references
#' {TBA}
#'
#' @export
#'
#'
pairCountPoint <- function(VAFdata,pid1,pid2,slope1,slope2) {
  r1=r2=r3=r4=r5=r6=r7=0
  r1 <- sum(VAFdata[,pid1]==0 & VAFdata[,pid2]==0) #x=y=0
  r2 <- sum(VAFdata[,pid1]<=VAF_ignore & VAFdata[,pid2]<=VAF_ignore & !(VAFdata[,pid1]==0 & VAFdata[,pid2]==0)  )   #0<x&y<=0.2, corner>0
  r3 <- sum(VAFdata[,pid1]>VAF_ignore & VAFdata[,pid2]<=VAF_cutoff) ## X only
  r4 <- sum(VAFdata[,pid2]>VAF_ignore & VAFdata[,pid1]<=VAF_cutoff) ## Y only
  r5 <- sum(VAFdata[,pid1]>VAF_ignore &  VAFdata[,pid2] > VAF_cutoff & VAFdata[,pid2]<slope1*VAFdata[,pid1]) ## x>>y
  r6 <- sum(VAFdata[,pid2]>VAF_ignore &  VAFdata[,pid1] > VAF_cutoff & VAFdata[,pid2]>slope2*VAFdata[,pid1]) ## y>>x
  r7 <- sum(VAFdata[,pid2]<=slope2*VAFdata[,pid1] & VAFdata[,pid2]>=slope1*VAFdata[,pid1] & !(VAFdata[,pid2]<=VAF_ignore & VAFdata[,pid1]<=VAF_ignore)) ## center
  tab <- t(c(r1,r2,r3,r4,r5,r6,r7))
  colnames(tab)=c("x0y0","xylt0d2","xonly","yonly","x2y","y2x","center")
invisible(tab)
}
