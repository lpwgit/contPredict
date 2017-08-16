#' Filter mutations base on coefficient of variant
#'
#'This function filters common SNPs base on coefficient of variant (COV) at the specified cutoff
#'
#' @param VAFdata data frame, variant allele frequency (VAF)
#' @param percentage numeric, COV cutoff in percentage (default:5)
#' @param output_path character, tmp directory to store log file (SNPcv.txt)
#'
#' @usage filterCommonSNPs(VAFdata, percentage, output_path)
#'
#' @return integer vector, row index of VAF data to be used
#'
#' @references
#' {TBA}
#'
#' @export
#'
filterCommonSNPs<- function(VAFdata, percentage=5,output_path){
  lastCol = dim(VAFdata)[2]
  excludeCol =c(1,lastCol-1, lastCol)
  snp.mean <- as.matrix(unlist(apply(VAFdata[,-excludeCol],1,mean) ))
  snp.sd <-  as.matrix(unlist(apply(VAFdata[,-excludeCol],1,sd) ))
  df <- data.frame(mean = snp.mean,sd = snp.sd)
  df$cov <- df$sd/df$mean
  rownames(df) <- VAFdata$mutationID
  outfile <- paste0(output_path,"/tmp/SNPcv.txt")
  write.table(df,outfile,quote=F,sep='\t')
  cov.quantile <- quantile(df$cov,na.rm = T,seq(0,1,0.01))
  perc <- paste0(percentage,"%")
  which(df$cov!='NaN' & df$cov>cov.quantile[perc]) ## keep

}
