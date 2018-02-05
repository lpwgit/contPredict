#' Identify true source from multisource contamination
#'
#' This function uses Fisher exact test to filter multiple sources that contaminate a particular target sample
#'
#' @param sample_pairs data frame, sample infomation
#' @param VAFdata data frame, mutation VAF
#' @param VAFcov data frame, mutation coverage
#' @param VAF_cutoff numeric, minimum VAF (default: 0.002)
#' @param VAF_cutoff1 numeric, minimum VAF to be considered a source mutation (default: 0.05)
#' @param p.val_cutoff numeric, pval cutoff to be considered a significant contamination source (default: 0.05)
#' @param output_path character, output directory
#'
#' @usage multipleSource(sample_pairs,VAFdata,VAFcov,VAF_cutoff,VAF_cutoff1,p.val_cutoff,output_path)
#'
#' @return list, remaining paiwise sample relationship information after Fisher exact test that eliminates insignificant contamination pair
#'
#' @references
#' {TBA}
#'
#' @export
#'
multipleSource <- function(sample_pairs,VAFdata,VAFcov,VAF_cutoff=0.002,VAF_cutoff1=0.05,p.val_cutoff=0.05,output_path){
  pid1.col <- which(colnames(sample_pairs) =="pid1")
  pid2.col <- which(colnames(sample_pairs) =="pid2")
  rel.col <- which(colnames(sample_pairs) =="rel")
  used_col <-c(pid1.col,pid2.col,rel.col) # pid1, pid2, rel
  flip_flag <-0
  n <- length(which(sample_pairs[,"rel"]!="-1" & sample_pairs[,"rel"]!="00" & sample_pairs[,"rel"]!="11"))

  relation_for_multipleSource_part1 <- cbind(sample_pairs[which(sample_pairs[,"rel"]!="-1" & sample_pairs[,"rel"]!="00" & sample_pairs[,"rel"]!="11"),used_col],rep(flip_flag,n))
  colnames(relation_for_multipleSource_part1) <-c("pid1","pid2","rel","flip")

  num.bothway <-  length(which(sample_pairs[,"rel"]=="11"))
  relation_for_multipleSource_part2<-c()
  if(num.bothway>0){
    flip_flag <-0
    n<-length(which(sample_pairs[,"rel"]=="11"))
    bothway <-cbind(sample_pairs[which(sample_pairs[,"rel"]=="11"),used_col],rep(flip_flag,n))
    colnames(bothway) <- c("pid1","pid2","rel","flip")
    flip <- apply(bothway,1,function(i){
      c(i['pid2'],i['pid1'],"01",1) # 2017 Mar 28
    })
    flip<-t(flip)
    colnames(flip) <-c("pid1","pid2","rel","flip")
    rownames(flip) <- paste(flip[,'pid1'],flip[,'pid2'],sep=':') ##

    relation_for_multipleSource_part2 <- rbind(bothway,flip)
  }

  if(length(relation_for_multipleSource_part2)!=0){
    if(nrow(relation_for_multipleSource_part1) >0 & nrow(relation_for_multipleSource_part2) >0){
      relation_for_multipleSource <- rbind(relation_for_multipleSource_part1,relation_for_multipleSource_part2)
    }else if(nrow(relation_for_multipleSource_part1) >0 & nrow(relation_for_multipleSource_part2) ==0){
      relation_for_multipleSource <- relation_for_multipleSource_part1
    }else if(nrow(relation_for_multipleSource_part1) ==0 & nrow(relation_for_multipleSource_part2) >0){
      relation_for_multipleSource <- relation_for_multipleSource_part2
    }
  }else if(length(relation_for_multipleSource_part1)!=0){
      relation_for_multipleSource <- relation_for_multipleSource_part1
  }else{
    print ("Error!")
  }

  outfile <- paste0(output_path,"/tmp/eliminated.list")
  if(file.exists(outfile)){
    file.remove(outfile)
  }

  logfile<- paste0(output_path,"/tmp/fisherTest.list")
  if(file.exists(logfile)){
    file.remove(logfile)
  }
  realSourceEliminatefile <- paste0(output_path,"/tmp/realSource_eliminatedPairs.list")
  if(file.exists(realSourceEliminatefile)){
    file.remove(realSourceEliminatefile)
    cat(paste("remove","target","keep", sep=" "),file=realSourceEliminatefile,sep="\n")
  }
  colnames(relation_for_multipleSource) <- c("pid1","pid2","rel","flip")
  if(dim(relation_for_multipleSource)[1]>0){
    ## multiple source(tumor-Td,metastatis-Tm,normal-Nb) to a target, remove source that is insignificant
    source_target <-apply(relation_for_multipleSource,1,function (i){
      if(i['rel'] =="01" ){
        source<-i['pid2']
        target<-i['pid1']
      }else if(i['rel'] =="10") { ##10
        source<-i['pid1']
        target<-i['pid2']
      }else if(i['rel'] =="11") { ##11 bothway
        source<-i['pid2']
        target<-i['pid1']
      }else{
        print (paste(i['rel'] ,"check!"))
      }
      c(source,target,i['rel'],i['flip'])
    })
    source_target <-t(source_target)
    colnames(source_target)<-c("source","target","rel","flip")
    target_duplct <- unique(source_target[duplicated(source_target[, "target"]), "target"])
    rel_table <- source_target
  }else{
    target_duplct<-NULL
  }

  if(length(target_duplct)>0){
    for(j in 1:length(target_duplct)) {
      target <- target_duplct[j]
      source <- unique(source_target[source_target[, "target"]==target, "source"])
      relation <-as.matrix(source_target[source_target[, "target"]==target, "rel"])
      flip.info<-as.matrix(source_target[source_target[, "target"]==target, "flip"])
      num.source = length(source)

      source_pair <- combn(source, 2) ## regardless of patient
      sourcepair_chk_fishtest <- apply(source_pair, 2, function(k) {
          s1 <- k[1]
          s2 <- k[2]

            tmp.data <- as.data.frame(cbind(VAFdata[,target],VAFdata[,s1],VAFdata[,s2],as.matrix(unlist(VAFdata[,'publicDB']))))
            colnames(tmp.data) <-c(target,s1,s2,"publicDB")
            rownames(tmp.data) <- VAFdata[,'mutationID']
            cols = c(1,2,3)
            tmp.data[,cols] = apply(tmp.data[,cols], 2, function(x) as.numeric(as.character(x)))

            category <- apply(tmp.data,1,function(i){
              class(i[s1]) <- "numeric"
              class(i[s2]) <- "numeric"
              class(i[target]) <- "numeric"

              if(i[s1]>VAF_cutoff&i[s2]>VAF_cutoff&i[target]>VAF_cutoff&i[target]<i[s1]&i[target]<i[s2] ){ # 2017 Apr24
                # both present
                out <- c(unlist(i['publicDB']),"present","both")
              }else if(i[s1]>VAF_cutoff1&i[s2]<=VAF_cutoff&i[target]<i[s1]&i[target]>VAF_cutoff){
                out <- c(unlist(i['publicDB']),"present",s1)
              }else if(i[s2]>VAF_cutoff1&i[s1]<=VAF_cutoff&i[target]<i[s2]&i[target]>VAF_cutoff){
                out <- c(unlist(i['publicDB']),"present",s2)
              }else if(((i[s1]>VAF_cutoff1&i[s2]>VAF_cutoff) | (i[s1]>VAF_cutoff&i[s2]>VAF_cutoff1)) & !i[target]>VAF_cutoff ){
                out <- c(unlist(i['publicDB']),"absent","both")
              }else if((i[s2]>VAF_cutoff1&i[s1]<=VAF_cutoff) & !i[target]>VAF_cutoff ){
                out <- c(unlist(i['publicDB']),"absent",s2) # 2017Mar31
              }else if((i[s1]>VAF_cutoff1&as.numeric(i[s2])<=VAF_cutoff) & !i[target]>VAF_cutoff ){
                out <- c(unlist(i['publicDB']),"absent",s1) # 2017Mar31
              } else{
                out <- c(unlist(i['publicDB']),"dummy","dummy")
              }
              out
            })

            category <- t(category)
            colnames(category) <- c("publicDB","set","sampleID")
            data <- as.data.frame(category[which(category[,'sampleID']!="dummy"),])
            s1_tgt = sum(data$set=="present" & data$sampleID==s1)
            s2_tgt = sum(data$set=="present" & data$sampleID==s2)
            s1_not_in_tgt = sum(data$set=="absent" & data$sampleID==s1)
            s2_not_in_tgt = sum(data$set=="absent" & data$sampleID==s2)
            both_in_tgt = sum(data$set == "present" & data$sampleID =="both")
            both_not_in_tgt = sum(data$set == "absent" & data$sampleID =="both")

            data.outfile <- paste0(output_path,"/tmp/",paste(target,s1,s2,"present_absent_info.txt",sep="_"))
            write.table(data,file=data.outfile,sep="\t",col.names = T,row.names = T,quote=F)

          d <- matrix(c(s1_tgt,  s1_not_in_tgt, s2_tgt,s2_not_in_tgt), 2, 2) #multisource
          p.val <- fisher.test(d,alternative='greater')$'p.value'

          ## log file
          logprt <- cbind(target, s1, s2,p.val)
          write.table(logprt,file=logfile,append=TRUE,sep="\t",col.names = F,row.names = F,quote=F)

          d2 <- matrix(c(s2_tgt,  s2_not_in_tgt, s1_tgt,s1_not_in_tgt), 2, 2) #multisource
          p.val2 <- fisher.test(d2,alternative='greater')$'p.value'

          ## log file
          logprt <- cbind(target, s2, s1,p.val2)
          write.table(logprt,file=logfile,append=TRUE,sep="\t",col.names = F,row.names = F,quote=F)

          if(p.val >= p.val_cutoff & p.val2 <p.val_cutoff){
            remove.flag <- c(FALSE,TRUE)
          }else if(p.val2>=p.val_cutoff & p.val <p.val_cutoff){
            remove.flag <-c(TRUE,FALSE)
          }else{
            remove.flag <-c(TRUE,TRUE)
          }

          final.del <- remove.flag
          final.del
        })
        source_del <- unique(source_pair[!sourcepair_chk_fishtest])

        if(length(source_del)>0){
          rel_table <- rel_table [!((rel_table [, "source"] %in% source_del)&rel_table [, "target"]==target),]
          if(length(source_del)==1){
            prtout <- cbind(source_del,target)
          }else{
            prtout <- cbind(source_del,rep(target,length(source_del)))
          }
          write.table(prtout,file=outfile,append=TRUE,sep="\t",col.names = F,row.names = F,quote=F)
        }
    } #for(i in 1:length(target_duplct)) {
  } #if(length(target_duplct)>0){

  ## check if 1 record
  final_rel <- rel_table
  tmp <- as.matrix(rel_table)

  if(dim(tmp)[2] == 1){ # only 1 record
    tmp <-t(tmp)
    row1 <- paste0(tmp[,1],":",final_rel[2])
    row2 <- paste0(tmp[,2],":",final_rel[1])
    found1 <- which(rownames(sample_pairs) == row1)
    found2 <- which(rownames(sample_pairs) == row2)
    if(length(found1) ==1){
      rownames(tmp) <- rownames(sample_pairs)[found1]
    }else if(length(found2) ==1){
      rownames(tmp) <- rownames(sample_pairs)[found2]
    }
    rel_table<-tmp
  }
rel_table
}
