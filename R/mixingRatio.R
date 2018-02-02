#' Calculate contamination level
#'
#' This function calculates contamination level, ratio of VAF_target to VAF_source
#'
#' @param VAFdata data frame, variant allele frequency (VAF)
#' @param VAFcov data frame, variant coverage
#' @param sample_pairs data frame, sample pairs information
#' @param final_rel data frame, relation information
#' @param VAF_cutoff numeric, minimum VAF (default: 0.002)
#' @param VAF_ignore numeric, ignore variant less than this cutoff (default: 0.2)
#' @param ALL_flag logical, TRUE: calculate all pairwise mixing ratio regardless of relationship, FALSE: calculate pairwise sample with relationship of 10|01|11
#' @param output_path character, output directory
#' @param sameSubject logical, TRUE: calculate contamination level for same subject pairs (default: FALSE)
#'
#' @usage mixingRatio(VAFdata,VAFcov,sample_pairs,final_rel,VAF_cutoff,VAF_ignore,ALL_flag,output_path,sameSubject)
#'
#' @return matrix, source sample, target sample, relation, contamination level, flip_flag
#'
#' @references
#' {TBA}
#'
#' @export
#'
mixingRatio <- function(VAFdata,VAFcov,sample_pairs,final_rel,VAF_cutoff=0.002,VAF_ignore=0.2, ALL_flag=FALSE,output_path,sameSubject=FALSE){

if(!sameSubject){
  curated.file <- paste0(output_path,"/tmp/manual_curated.list")

  if(file.exists(curated.file)){
    colClasses=c("character","character","numeric", "numeric","character")
    curated.info <-read.table(curated.file,colClasses = colClasses)
    colnames(curated.info) <- c("source","target","slope1","slope2","rel")
    curation.flag<-TRUE
  }else{
    curation.flag<-FALSE
    slope1 <- 0.5 # default
    slope2 <- 2.0  # default
  }
  pid1.col <- which(colnames(sample_pairs) =="pid1")
  pid2.col <- which(colnames(sample_pairs) =="pid2")
  rel.col <- which(colnames(sample_pairs) =="rel")
  used_col <-c(pid1.col,pid2.col,rel.col) # pid1, pid2, rel

  num.bothway <- length(which(final_rel[,"rel"]=="11"))
  if(num.bothway==0){
    n <- length(which(match(rownames(sample_pairs),rownames(final_rel))!='NA'))
    flip_flag <-0
    n2 <- dim(final_rel)[1]
    if(n>0&n==n2){
      relation_for_mixingRatio <- cbind(sample_pairs[which(match(rownames(sample_pairs),rownames(final_rel))!='NA'),used_col],rep(flip_flag,n))
    }else if(n>0&n2>n){
        p1 <- cbind(sample_pairs[which(match(rownames(sample_pairs),rownames(final_rel))!='NA'),used_col],rep(flip_flag,n))
        p2 <- as.matrix(final_rel[!rownames(final_rel) %in% rownames(sample_pairs),])
        if(dim(p2)[2] ==1){
          p2 <- t(p2)
        }
        colnames(p1)<- c("pid1","pid2","rel","flip")
        colnames(p2)<- c("pid1","pid2","rel","flip")
        p2[,'rel'] <-"10"
        relation_for_mixingRatio <- rbind(p1,p2)
    }else{
      stop("Check! no candidate pairs for calculating contamination level")
    }
  }else{
    bothwayid <-which(final_rel[,"rel"]=="11")
    flip1 <-cbind(final_rel[bothwayid,'target'],final_rel[bothwayid,'source'],"10",1)
    flip2 <-cbind(final_rel[bothwayid,'source'],final_rel[bothwayid,'target'],"10",1)
    if(dim(flip1)[2] ==1){
      flip1=t(flip1)
      flip2=t(flip2)
    }
    colnames(flip1)<- c("pid1","pid2","rel","flip")
    colnames(flip2)<- c("pid1","pid2","rel","flip")
    n<- length(which(match(rownames(sample_pairs),rownames(final_rel))!='NA'))
    flip_flag <-0
    nonflip <- cbind(sample_pairs[which(match(rownames(sample_pairs),rownames(final_rel))!='NA'),used_col],rep(flip_flag,n))
    colnames(nonflip)<- c("pid1","pid2","rel","flip")
    relation_for_mixingRatio <- rbind(nonflip,flip1,flip2)
  }
}else{
  curation.flag<-FALSE
  same.subject <- sample_pairs[sample_pairs$rel=="00",]
  same.subject[,'pid1']<- as.character(same.subject[,'pid1'])
  same.subject[,'pid2']<-  as.character(same.subject[,'pid2'])
  same.subject[,'rel']<-  as.character(same.subject[,'rel'])
  relation_for_mixingRatio <- data.frame(cbind(pid1=same.subject$pid1,pid2=same.subject$pid2,rel=same.subject$rel,flip=rep(0,nrow(same.subject))),stringsAsFactors = FALSE)
}

  colnames(relation_for_mixingRatio) <- c("pid1","pid2","rel", "flip")
  ## calculate mixing ratio for 1-way & bothway
  mixingR<- NULL
  mixingR <- apply(relation_for_mixingRatio, 1, function(i) {
    v2_b <-0
    b<-0
    pval<-0
    v2_pval <-0
    flip <- as.character(unlist(i["flip"]))
    relation <-i["rel"]

    if((i["rel"] =="00" | i["rel"] =="-1") & !ALL_flag){
      source <- as.character(unlist(i["pid1"]))
      target <- as.character(unlist(i["pid2"]))
      c(source,target,rel=as.character(unlist(i["rel"])) , lm_coeff=0, mixingR_CovVaf = 0, mixingR_Cov =0,lmPval = 0,v2_lm_coeff=0, v2_mixingR_CovVaf = 0, v2_mixingR_Cov = 0,v2_lm_pval =0,flip=0)
    }
    else{
      if(i["rel"] =="01" | i["rel"] =="11" ){
        source <- as.character(unlist(i["pid2"]))
        target <- as.character(unlist(i["pid1"]))
      }else if(i["rel"] == "10"){
        source <- as.character(unlist(i["pid1"]))
        target <- as.character(unlist(i["pid2"]))
      }else if((i["rel"] == "00" | i["rel"] == "-1" ) & ALL_flag){
        source <- as.character(unlist(i["pid2"]))
        target <- as.character(unlist(i["pid1"]))
      }
      pid1 <-  as.character(unlist(i["pid1"]))
      pid2 <-  as.character(unlist(i["pid2"]))

      vaf_data<- cbind(VAFdata[,pid1],VAFdata[,pid2])
      colnames(vaf_data)<-c("VAF1","VAF2")

      cov_data<- cbind(VAFcov[,pid1],VAFcov[,pid2])
      colnames(cov_data)<-c("COV1","COV2")

      data<- as.data.frame(cbind(vaf_data,cov_data))

      if(curation.flag){
        curated <- which(curated.info$target == target & curated.info$source== source)

        if(length(curated)==1){
          slope1 <- as.numeric(curated.info[curated,'slope1'])
          slope2 <- as.numeric(curated.info[curated,'slope2'])
          manual_rel <-curated.info[curated,'rel']
          relation <- manual_rel
          if(manual_rel =="01" | manual_rel =="11" ){
            source <- as.character(unlist(i["pid2"]))
            target <- as.character(unlist(i["pid1"]))
          }else if(manual_rel == "10"){
            source <- as.character(unlist(i["pid1"]))
            target <- as.character(unlist(i["pid2"]))
          }else if((manual_rel == "00"| manual_rel == "-1" ) & ALL_flag){
            source <- as.character(unlist(i["pid2"]))
            target <- as.character(unlist(i["pid1"]))
          }
        }else{
          slope1 <- 0.5
          slope2 <- 2.0
        }
      }else{
        slope1 <- 0.5
        slope2 <- 2.0
      }
      if(target==pid1) {
        data1<-subset(data, VAF1>VAF_ignore|VAF2>VAF_ignore) # VAF1<VAF2 Y2X
        data1<-subset(data1, VAF1<=slope1*VAF2) #V1

        v2_data<-subset(data, VAF1>VAF_ignore|VAF2>VAF_ignore) # VAF1<VAF2 Y2X
        v2_data<-subset(v2_data, VAF1<VAF2) #V2

        if(nrow(data1)>=1){
          data1$minCov<-pmin(data1$COV1, data1$COV2)
          data1$R<-data1$VAF1/data1$VAF2
          data1$maxVAF <- pmax(data1$VAF1, data1$VAF2)
          fit1 <- lm(R~1, data = data1, weights = minCov)
          b <- as.numeric(coef(fit1)[1])
          pval <- summary(fit1)$coefficients[,4]

        }
        mixR1 <- sum(data1$R*data1$minCov)/sum(data1$minCov)
        mixR <- sum(data1$R*data1$minCov*data1$maxVAF)/sum(data1$minCov*data1$maxVAF)

        if(nrow(v2_data)>=1){
          v2_data$minCov<-pmin(v2_data$COV1, v2_data$COV2)
          v2_data$R<-v2_data$VAF1/v2_data$VAF2
          v2_data$maxVAF <- pmax(v2_data$VAF1, v2_data$VAF2)
          v2_fit1 <- lm(R~1, data = v2_data, weights = minCov)
          v2_b <- as.numeric(coef(v2_fit1)[1])
          v2_pval <- summary(v2_fit1)$coefficients[,4]

        }
        v2_mixR1 <- sum(v2_data$R*v2_data$minCov)/sum(v2_data$minCov)
        v2_mixR <- sum(v2_data$R*v2_data$minCov*v2_data$maxVAF)/sum(v2_data$minCov*v2_data$maxVAF)

      } else { # target == pid2
        data1<-subset(data, VAF1>VAF_ignore|VAF2>VAF_ignore) # VAF1>VAF2 X2Y
        data1<-subset(data1, VAF1>=slope2*VAF2) # & VAF2>0.02) #V1

        v2_data<-subset(data, VAF1>VAF_ignore|VAF2>VAF_ignore) # VAF1<VAF2 Y2X
        v2_data<-subset(v2_data, VAF1<VAF2) #V2

        if(nrow(data1)>=1){
          data1$minCov<-pmin(data1$COV1, data1$COV2)
          data1$R<-data1$VAF2/data1$VAF1
          data1$maxVAF <- pmax(data1$VAF1, data1$VAF2)
          fit2 <- lm(R~1, data = data1, weights = minCov)
          b <- as.numeric(coef(fit2)[1])
          pval <- summary(fit2)$coefficients[,4]
        }
        mixR1 <- sum(data1$R*data1$minCov)/sum(data1$minCov)
        mixR <- sum(data1$R*data1$minCov*data1$maxVAF)/sum(data1$minCov*data1$maxVAF)

        if(nrow(v2_data)>=1){
          v2_data$minCov<-pmin(v2_data$COV1, v2_data$COV2)
          v2_data$R<-v2_data$VAF1/v2_data$VAF2
          v2_data$maxVAF <- pmax(v2_data$VAF1, v2_data$VAF2)
          v2_fit1 <- lm(R~1, data = v2_data, weights = minCov)
          v2_b <- as.numeric(coef(v2_fit1)[1])
          v2_pval <- summary(v2_fit1)$coefficients[,4]
        }
        v2_mixR1 <- sum(v2_data$R*v2_data$minCov)/sum(v2_data$minCov)
        v2_mixR <- sum(v2_data$R*v2_data$minCov*v2_data$maxVAF)/sum(v2_data$minCov*v2_data$maxVAF)
      }
      num_round_digit = 5

      if(length(pval)==0){
        pval =1
      }
      if(length(v2_pval)==0){
        v2_pval =1
      }
      if(is.nan(pval)){
        pval =1
      }
      if(is.nan(v2_pval)){
        v2_pval =1
      }
      c(source,target,rel=as.character(unlist(i["rel"])) , lm_coeff=round(b,digits=num_round_digit), mixingR_CovVaf = round(mixR, digits = num_round_digit), mixingR_Cov = round(mixR1, digits = num_round_digit),lmPval = round(pval, digits = num_round_digit), v2_lm_coeff=round(v2_b,digits=num_round_digit), v2_mixingR_CovVaf = round(v2_mixR, digits = num_round_digit), v2_mixingR_Cov = round(v2_mixR1, digits = num_round_digit),v2_lmPval = round(v2_pval, digits = num_round_digit), flip=flip)
    }# !00
  })

  mixingR <- t(mixingR)
  colnames(mixingR)<-c ("source","target","rel","lm_coeff","mixingR_CovVaf","mixingR_Cov","lm_pval","v2_lm_coeff","v2_mixingR_CovVaf","v2_mixingR_Cov","v2_lm_pval","flip")
mixingR
}
