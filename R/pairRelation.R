#'Determine relationship between pairwise sample
#'
#'This function identifies relationship of pairwise sample [00: same patient| 10: X contaminate Y|01 Y contaminate X| 11: both way contamination]
#'
#' @param sample_pairs data frame, sample information
#' @param center_cutoff numeric, cutoff to be classified as same subject samples (default:0.5)
#' @param source_cutoff numeric, cutoff to be classififed as source sample (default:0.4)
#' @param target_cutoff numeric, cutoff to be classified as target sample (default:0.1)
#' @param localPcomm_cutoff numeric, cutoff to manage high or low SNPs sharing between a pair of sample
#' @param region_cutoff numeric, cutoff for case with low SNP sharing (default:0.75)
#' @param num_round_digit integer, rounding numeric up to this decimal point (default:3)
#'
#' @usage pairRelation(sample_pairs,center_cutoff,source_cutoff,target_cutoff,localPcomm_cutoff,region_cutoff,num_round_digit,output_path)
#'
#' @return matrix, 2 columns: pcommon and relation
#' @references
#' {TBA}
#'
#' @export
#'
pairRelation <- function(sample_pairs,center_cutoff=0.5,source_cutoff=0.4,target_cutoff=0.1,localPcomm_cutoff=0.1,region_cutoff=0.75,num_round_digit=3,output_path){
  ##

    #logfile <- paste0(output_path,"/tmp/pairRelation_log.txt")
    #out <- t(c('pid1','pid2',"ncommon","pcommon","pn1","pn2","pn3","qn1","qn3","rn1","rn3"))
    #write.table(out,quote=F,sep='\t',logfile,append=F,row.names = F,col.names = F)

  relation <- apply(sample_pairs, 1, function(i) {
    pn1=pn2=pn3=qn1=qn3=rn1=rn3=0
    rel <-"-1"

    ncommon <- as.numeric(i["x2y"])+as.numeric(i["y2x"])+as.numeric(i["center"])
    pcommon <- ncommon /(ncommon + as.numeric(i["xonly"])+as.numeric(i["yonly"]) )
    pcommon <- round(pcommon, digits = num_round_digit)

    pn1 <- as.numeric(i["x2y"])/ncommon
    pn1 <- round(pn1,digits =num_round_digit)

    pn2 <- as.numeric(i["center"])/ncommon
    pn2 <- round(pn2,digits =num_round_digit)

    pn3 <- as.numeric(i["y2x"])/ncommon
    pn3 <- round(pn3,digits =num_round_digit)

    qn1 <- as.numeric(i["x2y"])/(as.numeric(i["x2y"])+as.numeric(i["xonly"]))
    qn1 <- round(qn1,digits =num_round_digit)
    qn3 <- as.numeric(i["y2x"])/(as.numeric(i["y2x"])+as.numeric(i["yonly"]))
    qn3 <- round(qn3,digits =num_round_digit)

    rn1 <- as.numeric(i["x2y"])/(as.numeric(i["x2y"])+as.numeric(i["center"]))
    rn1 <- round(rn1,digits =num_round_digit)
    rn3 <- as.numeric(i["y2x"])/(as.numeric(i["y2x"])+as.numeric(i["center"]))
    rn3 <- round(rn3,digits =num_round_digit)

    if(is.na(qn1)){
      qn1=0
    }
    if(is.na(qn3)){
      qn3=0
    }
    if(is.na(rn1)){
      rn1=0
    }
    if(is.na(rn3)){
      rn3=0
    }

    if(is.na(pn1)){
      pn1=0
    }
    if(is.na(pn2)){
      pn2=0
    }
    if(is.na(pn3)){
      pn3=0
    }

    if(pcommon > localPcomm_cutoff){
      if(pn2>center_cutoff){
        rel <-"00"
      }else if(pn1>source_cutoff & pn3<target_cutoff){ #x2y
        rel<-"10"
      }else if(pn3>source_cutoff & pn1<target_cutoff){ #y2x
        rel<-"01"
      }else{ #bothway
        rel<-"11"
      }
    }else{## Mar 24 2017
      cutoff <- 20

      if((qn1>target_cutoff&rn1>region_cutoff) & (qn3<=target_cutoff&rn3<=region_cutoff)){
        rel <- "10" #x2y
      }else if((qn3>target_cutoff&rn3>region_cutoff) & (qn1<=target_cutoff&rn1<=region_cutoff)){
        rel <- "01" #y2x
      }else if((qn1>target_cutoff&rn1>region_cutoff) & (qn3>target_cutoff&rn3>region_cutoff)){
        rel <- "11" #x<->y
      }else if(rn1>region_cutoff|rn3>region_cutoff){ # little portion of center sharing
        if(qn1>qn3&qn3>0){
          ratio <- qn1/qn3
          if(ratio>=cutoff){
            rel <- "10"
          }
        }else if(qn3>qn1&qn1>0){
          ratio <- qn3/qn1
          if(ratio>=cutoff){
            rel <- "01"
          }
        }
      }else{ # center sharing larger portion
        if(qn1>qn3&qn3>0){
          ratio <- qn1/qn3
          if(ratio>=cutoff){
            rel <- "10"
          }
        }else if(qn3>qn1&qn1>0){
          ratio <- qn3/qn1
          if(ratio>=cutoff){
            rel <- "01"
          }
        }
      }
    }

    #out <- t(c(i['pid1'],i['pid2'],ncommon,pcommon,pn1,pn2,pn3,qn1,qn3,rn1,rn3))
    #write.table(out,quote=F,sep='\t',logfile,append=T,row.names = F,col.names = F)
    c(pcommon,rel)
  })
  relation = t(relation)
  colnames(relation) <- c("pcommon","rel")
  relation
}
