#' Main function to run sample contamination analysis
#'
#' This function detects sample contamination base on variant allele frequency (VAF) and variant coverage
#'
#' @param data_path character, VAF and COV data directory
#' @param output_path character, output directory
#' @param config_file character, configuration file
#' @param filterCOV numeric, below this cutoff
#' @param rmSNPcv_cutoff numeric, TRUE: filter SNPs with low covariance of coefficient (COV), FALSE: (default:FALSE)
#' @param manualsetPara logical, TRUE: parameter setting base on distribution of pairwise sample commonality, FALSE: default setting as in configuration file (default:FALSE)
#' @param manual_localPcomm numeric, localPcomm for classification low/high contamination (default: 10)
#' @param manual_center numeric, manual center cutoff for same subject determination (default: 50)
#'
#' @usage run_sampleContamination(data_path, output_path, config_file="config.txt",rmSNPcv_cutoff=0,manualsetPara=FALSE,manual_localPcomm=10,manual_center =50,filterCOV =0)
#'
#' @return contamination tab delimited text file with 3 columns: source sample, target sample, contamination percentage.Contamination and same subject circos plot figure in pdf
#'
#' @references
#' {TBA}
#'
#' @export
#'
run_sampleContamination<- function( data_path, output_path , config_file="config.txt",rmSNPcv_cutoff=0,manualsetPara=FALSE,manual_localPcomm=10,manual_center =50,filterCOV =0) {

ownLog<<-F
if(!file.exists(config_file) ){
  stop("Could not find config file!")
}

options(stringsAsFactors = FALSE)
options(warn=-1)
suppressPackageStartupMessages(library("sampleCont"))

# set parameters
setParameterConfig(config_file)

# set VAF and COV data file
vaf_file <<- paste0(data_path,"/VAF.out")
cov_file <<- paste0(data_path,"/VAF_cov.out")

## ensure data and output path exists
if (!file.exists(data_path)){
    dir.create(file.path(data_path))
}
if (!file.exists(output_path)){
  dir.create(file.path(output_path))
}
if (!file.exists(paste0(output_path,"/tmp"))){
  dir.create(file.path(paste0(output_path,"/tmp")))
}

## load data
VAFdata <<- inputData(vaf_file)
VAFcov <<- inputData(cov_file)

## get sample ids
sampleID <- colnames(VAFcov)[-1]

## remove snps with covariance of coefficient below cutoff
if(rmSNPcv_cutoff>0){
  	keep <- filterCommonSNPs(VAFdata,rmSNPcv_cutoff,output_path)
  	length(keep)
  	tmp <- VAFdata[keep,]
  	VAFdata <<- tmp
  	tmp <- VAFcov[keep,]
  	VAFcov <<- tmp
}

if(filterCOV>0){
  VAFdata <<- filterByCov(VAFdata,VAFcov,filterCOV)
}

## num snps per sample
SNPcount <- numMutationperSample(VAFdata,VAF_cutoff,n_sample)

process.msg <- "Processing SNPs..."
cat(process.msg,"\n")

## paiwise sample common SNPs
SNPshare<- pairShare(VAFdata,VAF_cutoff,n_sample)

## create all possible pairwise sample
delimeter <- ":"
pairs <- pairList(VAFdata,n_sample,delimeter)

## calculate pcomm
pcomm <- pairPCommon(VAFdata,VAF_cutoff,num_round_digit,n_sample)

## generate sample pairs info
pids <- t(sapply(as.matrix(pairs), function(i) unlist(strsplit(i, delimeter))))
colnames(pids) <- c("pid1", "pid2")
sample_pairs <- data.frame(pairID = pairs, pids)

## count mutation in 7 regions
countPoint <- regionCountMutation(sample_pairs,VAFdata,SNPcount,SNPshare,VAF_cutoff,VAF_ignore,n_sample)
sample_pairs <- cbind(sample_pairs, countPoint)

## manual set some parameters (localPcomm_cutoff,center_cutoff, target_cutoff, source_cutoff, region_cutoff) ignore those in configuration file
if(manualsetPara){
  ## manual localPcomm
  view.pcomm <- pcomm
  diag(view.pcomm) <-0
  pcomm.summ <- as.matrix(quantile(view.pcomm[view.pcomm>0],seq(0,1,0.01)))
  perc <- paste0(manual_localPcomm,"%")
  localPcomm_cutoff <- round(pcomm.summ[perc,],2)

  ## manual center_cutoff
  if(manual_center!=50){
    perc <- paste0(manual_center,"%")
    center_cutoff <- round(pcomm.summ[perc,],2)
    target_cutoff <-round((1-center_cutoff)/3,3)
    source_cutoff <- round((1-center_cutoff - target_cutoff),3)
    region_cutoff <- round((1-center_cutoff)/2,3)
  }
  rtn.out <-as.data.frame(cbind(pcomm=localPcomm_cutoff,center=center_cutoff,target=target_cutoff,source=source_cutoff,region=region_cutoff,n_sample =n_sample))

}else{
  rtn.out <-as.data.frame(cbind(pcomm=localPcomm_cutoff,center=center_cutoff,target=target_cutoff,source=source_cutoff,region=region_cutoff,n_sample =n_sample))
}

## identify relation btw pairwise sample
rel <- pairRelation(sample_pairs,center_cutoff,source_cutoff,target_cutoff,localPcomm_cutoff,region_cutoff,num_round_digit,output_path)

sample_pairs <- cbind(sample_pairs, rel)

if(ownLog==T){
  write.table(sample_pairs,quote=F,sep='\t',paste0(output_path,"/tmp/sample_pairs_info.txt"))
}

## check present of contamination pairs
continue.analysis <- TRUE
if(length(which(rel[,'rel'] !="00" & rel[,'rel'] !="-1")) ==0){
  cat("No contamination found!\n")
  if(length(which(rel[,'rel'] =="00")) >0){
  	cat("Only same subject circos\n")
  	same.subject <- sample_pairs[sample_pairs$rel=="00",]
  	same.subject.mr <- as.data.frame(mixingRatio(VAFdata,VAFcov,sample_pairs=sample_pairs,VAF_cutoff = VAF_cutoff,VAF_ignore = VAF_ignore,ALL_flag = TRUE,sameSubject=TRUE))
  	same.subject.out <- as.data.frame(cbind(source=same.subject.mr$source,target=same.subject.mr$target,link=round(10* as.numeric(same.subject.mr$lm_coeff),num_round_digit),rel=same.subject.mr$rel))
  	outfile <- paste0(output_path,"/tmp/circos_plotdata.txt")
  	write.table(same.subject.out,outfile,quote=F,sep='\t',row.names = F)
  	plot_circos_link (outfile,R=300,W=30,plotsize=800,titleStr="same subject only",seg.lab.size = 1.3,fig.file="sameSubjectOnly_circos_test.pdf",contaminatedOnly=FALSE,sameSubjectOnly=TRUE,allSamples=sampleID, output_path )
  	continue.analysis <-FALSE
  }
}

if(continue.analysis){
	## check multi-sources contamination
	final_rel <- multipleSource(sample_pairs,VAFdata ,VAFcov,VAF_cutoff,VAF_cutoff1,p.val_cutoff,output_path)
	write.table(as.matrix(final_rel),quote=F,sep='\t',paste(output_path,"/tmp/",n_sample,"sample_pairs_finalRel_multiSource.txt",sep = ''),row.names = F) #log

	## calculate mixing ratio
	process.msg <- "Calculating contamination level..."
	cat(process.msg,"\n")
	mr <- mixingRatio(VAFdata,VAFcov,sample_pairs,final_rel,VAF_cutoff,VAF_ignore,FALSE,output_path)
	write.table(as.matrix(mr),quote=F,sep='\t',paste(output_path,"/tmp/",n_sample,"sample_pairs_mixingRatio_b4elimination.txt",sep = ''),row.names = F) #own ref

	final.mr <- data.frame(cbind(mr[,'source'],mr[,'target'],as.numeric(mr[,'lm_coeff'])*100,mr[,'rel'],mr[,'flip']))
  colnames(final.mr) <- c('source','target','predicted_contamination_perc','rel','flip')

  ## write contamination output
  printmr <- data.frame(cbind(source=final.mr$source,target=final.mr$target,predicted_contamination_perc=final.mr$predicted_contamination_perc))
	write.table(unique(printmr),quote=F,sep='\t',paste(output_path,"/tmp/",center_cutoff,"center_",source_cutoff,"source_",target_cutoff,"target_",localPcomm_cutoff,"pcomm_",region_cutoff,"region_",n_sample,"sample_contaminationLevel.txt",sep = ''),row.names = F) ## final output
	result.file = paste(output_path,"/tmp/",center_cutoff,"center_",source_cutoff,"source_",target_cutoff,"target_",localPcomm_cutoff,"pcomm_",region_cutoff,"region_",n_sample,"sample_contaminationLevel.txt",sep = '')
	final.pred <- as.data.frame(filterMultiSources(result.file,output_path,VAFdata,VAFcov,uniq_both = 2))
} ## contamination exists

## clean up
tmp=paste0(output_path,"/tmp")
unlink(tmp, TRUE)

process.msg <-"Done!\nOutput written to:"
cat(process.msg, output_path,"\n")
rtn.out
} ## end function
