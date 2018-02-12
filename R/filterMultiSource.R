#' Filter multisource contamination
#'
#' This function filter multisource contamination
#'
#' @param filename character, contamination file name
#' @param output_path character, output directory
#' @param VAFdata data frame, mutation VAF
#' @param VAFcov data frame, mutation coverage
#' @param uniq_both integer, flag to control type of overlapping mutation in source and target samples
#' @usage filterMultiSources(filename,output_path,VAFdata,VAFcov,uniq_both)
#'
#' @return write contPredict.txt to output_path
#'
#' @references
#' {TBA}
#'
#' @export
#'
filterMultiSources <- function (filename,output_path,VAFdata,VAFcov,uniq_both=2){
  if(missing(filename)){
    stop("Empty filename!")
  }
  if(!file.exists(filename)){
    stop("Could not find file!")
  }

  outfile1 <- paste0(output_path,"/tmp/s1_vaf_cov.txt")
  outfile2 <- paste0(output_path,"/tmp/s2_vaf_cov.txt")
  outfile <- paste0(output_path,"/tmp/targetOnly_vaf_cov.txt")

  out <- cbind("target","s1","s2",'mutationID',"s1_vaf","target_vaf","vaf_diff","vaf_ratio","mutationID2","s1_cov","s1_cov","target_cov")
  write.table(out,quote = F,col.names=F, row.names = F,sep="\t",outfile1 )
  out <- cbind("target","s1","s2",'mutationID',"s1_vaf","target_vaf","vaf_diff","vaf_ratio","mutationID2","s1_cov","s1_cov","target_cov")
  write.table(out,quote = F,col.names=F, row.names = F,sep="\t",outfile2 )
  out <- cbind("target","s1","s2",'mutationID',"s1_vaf","s2_vaf","target_vaf","mutationID2","s1_cov","s1_cov","target_cov")
  write.table(out,quote = F,col.names=F, row.names = F,sep="\t",outfile )

  minCov <- 0
  predict.cont <- as.matrix(read.table(filename,header=T,check.names = F))
  process.predict.cont <- predict.cont
  dupl.target <- unique(predict.cont[duplicated(predict.cont[,'target']),'target'])
  n.dupl <- length(dupl.target)
  first = TRUE # append / write new file

  for(i in 1:n.dupl){
    current <- predict.cont[which(predict.cont[,"target"] == dupl.target[i]),]
    multi.source <-data.frame()
    if(nrow(current) ==2){

      s1 <- current[1,"source"]
      s2 <- current[2,"source"]

      f1 <- paste0(output_path, "/tmp/",dupl.target[i],"_",current[1,"source"],"_",current[2,"source"],"_present_absent_info.txt")
      f2 <- paste0(output_path, "/tmp/",dupl.target[i],"_",current[2,"source"],"_",current[1,"source"],"_present_absent_info.txt")
      if(file.exists(f1)){
        indata <- as.data.frame(read.table(f1,header=T))
      }else if(file.exists(f2)){
        indata <- as.data.frame(read.table(f2,header=T))
      }
      consider_3samples <- TRUE
      if(consider_3samples){
        # uniq in s1
        s1_tgt = as.matrix(rownames(indata)[indata$set=="present" & indata$sampleID==s1])

        # uniq in s2
        s2_tgt = as.matrix(rownames(indata)[indata$set=="present" & indata$sampleID==s2])

        # get both
        both_in_tgt = as.matrix(rownames(indata)[indata$set == "present" & indata$sampleID =="both"])

        s1_tgt_both <- rbind(s1_tgt,both_in_tgt)
        s2_tgt_both <- rbind(s2_tgt,both_in_tgt)
        s1_s2_tgt_both <- rbind(s1_tgt,s2_tgt,both_in_tgt)
      }else{
        # uniq in s1
        col.id <- c("mutationID",s1,s2,dupl.target[i])
        data<-VAFdata[, col.id]
        colnames(data) <- c("mutationID","s1", "s2","target")
        data1<-subset(data, s1>VAF_ignore)
        s1.data <- subset(data1,s1>VAF_cutoff1&target>VAF_cutoff&target<s1)

        data2<-subset(data, s2>VAF_ignore)
        s2.data <- subset(data2,s2>VAF_cutoff1&target>VAF_cutoff&target<s2)

        s1_tgt = as.matrix(s1.data$mutationID)

        # uniq in s2
        s2_tgt = as.matrix(s2.data$mutationID)

        # get both
        both_in_tgt = as.matrix(unique(intersect(s1_tgt,s2_tgt)))

        s1_tgt_both <- unique(rbind(s1_tgt,both_in_tgt))
        s2_tgt_both <- unique(rbind(s2_tgt,both_in_tgt))
        s1_s2_tgt_both <- unique(rbind(s1_tgt,s2_tgt,both_in_tgt))
      }

      if(uniq_both == 2){ ## uniq in s1/s2 and both
        col.id <- c("mutationID",s1)
        s1.vaf <- VAFdata[VAFdata$mutationID %in% s1_tgt_both, col.id] ## uniq in s1/s2 and both

        col.id <- c("mutationID",s2)
        s2.vaf <- VAFdata[VAFdata$mutationID %in% s2_tgt_both, col.id] ## uniq in s1/s2 and both

        col.id <- c("mutationID",dupl.target[i])
        s1.target.vaf <- VAFdata[VAFdata$mutationID %in% s1_tgt_both,  col.id]
        s2.target.vaf <- VAFdata[VAFdata$mutationID %in% s2_tgt_both,  col.id]

        col.id <- c("mutationID",s1,s2,dupl.target[i])
        s1.s2.target.cov <- VAFcov[VAFcov$mutationID %in% s1_tgt_both,  col.id]
        s2.s1.target.cov <- VAFcov[VAFcov$mutationID %in% s2_tgt_both,  col.id]

      }else if(uniq_both ==1){ ## both only
        col.id <- c("mutationID",s1)
        s1.vaf <- VAFdata[VAFdata$mutationID %in% both_in_tgt, col.id]
        col.id <- c("mutationID",s2)
        s2.vaf <- VAFdata[VAFdata$mutationID %in% both_in_tgt, col.id]
        col.id <- c("mutationID",dupl.target[i])
        s1.target.vaf <- VAFdata[VAFdata$mutationID %in% both_in_tgt,  col.id]
        s2.target.vaf <- VAFdata[VAFdata$mutationID %in% both_in_tgt,  col.id]

        col.id <- c("mutationID",s1,s2,dupl.target[i])
        s1.s2.target.cov <- VAFcov[VAFcov$mutationID %in% both_in_tgt,  col.id]
        s2.s1.target.cov <- VAFcov[VAFcov$mutationID %in% both_in_tgt,  col.id]

      }else if(uniq_both ==0){ ## uniq only
        col.id <- c("mutationID",s1)
        s1.vaf <- VAFdata[VAFdata$mutationID %in% s1_tgt, col.id]
        col.id <- c("mutationID",s2)
        s2.vaf <- VAFdata[VAFdata$mutationID %in% s2_tgt, col.id]
        col.id <- c("mutationID",dupl.target[i])
        s1.target.vaf <- VAFdata[VAFdata$mutationID %in% s1_tgt,  col.id]
        s2.target.vaf <- VAFdata[VAFdata$mutationID %in% s2_tgt,  col.id]

        col.id <- c("mutationID",s1,s2,dupl.target[i])
        s1.s2.target.cov <- VAFcov[VAFcov$mutationID %in% s1_tgt,  col.id]
        s2.s1.target.cov <- VAFcov[VAFcov$mutationID %in% s2_tgt,  col.id]
      }
      cont1 <- as.numeric(current[which(current[,'source'] == s1),"predicted_contamination_perc"])
      cont2 <- as.numeric(current[which(current[,'source'] == s2),"predicted_contamination_perc"])

      s1.target <- cbind(s1.vaf[order(s1.vaf$mutationID),], s1.target.vaf[order(s1.target.vaf$mutationID),2])
      colnames(s1.target)<-c("mutationID","source","target")
      s1.target$diff <- s1.target$source - s1.target$target
      s1.target$ratio <- s1.target$source / s1.target$target

      s1.diff.median <- median(s1.target$diff)
      s1.diff.sd <- sd(s1.target$diff)
      s1.diff.mean <- mean(s1.target$diff)
      s1.diff.max <- max(s1.target$diff)

      s2.target <- cbind(s2.vaf[order(s2.vaf$mutationID),], s2.target.vaf[order(s2.target.vaf$mutationID),2])
      colnames(s2.target)<-c("mutationID","source","target")
      s2.target$diff <- s2.target$source - s2.target$target
      s2.target$ratio <- s2.target$source / s2.target$target

      s2.diff.median <- median(s2.target$diff)
      s2.diff.max <- max(s2.target$diff)
      s2.diff.sd <- sd(s2.target$diff)
      s2.diff.mean <- mean(s2.target$diff)

      col.id <- s1
      s1.median.cov <- median(s1.s2.target.cov[,col.id])
      s1.min.cov <- min(s1.s2.target.cov[,col.id])

      col.id <- c(s1,dupl.target[i])
      s1.target.min.cov <- apply(s1.s2.target.cov[,col.id],1, function(each_row){
          min(each_row)
        })
      s1.target.min.cov <- as.matrix(unlist(s1.target.min.cov))

      col.id <- s2
      s2.median.cov  <- median(s2.s1.target.cov[,col.id])
      s2.min.cov  <- min(s2.s1.target.cov[,col.id])

      col.id <- c(s2,dupl.target[i])
      s2.target.min.cov <- apply(s2.s1.target.cov[,col.id],1, function(each_row){
        min(each_row)
      })
      s2.target.min.cov <- as.matrix(unlist(s2.target.min.cov))

      ## in target not in s1|s2
      col.id <- c("mutationID",s1,s2,dupl.target[i])
      current.vaf.data <- VAFdata[,col.id]
      current.cov.data <- VAFcov[,col.id]
      current.data <- cbind(current.vaf.data,current.cov.data)
      colnames(current.data)<- c("mutationID","s1","s2","target", "mutationID2","s1_cov","s2_cov","target_cov")

      current.data <- subset(current.data,target>VAF_cutoff&target>VAF_ignore)

      exclude.mut <- rbind(s1_tgt,s2_tgt,both_in_tgt)
      current.data <- current.data[!current.data$mutationID %in% exclude.mut,]

      target_only <- current.data$mutationID
      s1.s2.target.only.vaf <-VAFdata[VAFdata$mutationID %in% target_only,  col.id]
      s1.s2.target.only.cov <-VAFcov[VAFcov$mutationID %in% target_only,  col.id]

      min.cov.s1.target.only <- pmin(s1.s2.target.only.cov[,2],s1.s2.target.only.cov[,4]) ## snps in target not in s1/s2
      min.cov.s2.target.only <- pmin(s1.s2.target.only.cov[,3],s1.s2.target.only.cov[,4])

      ## log
      out <- cbind(dupl.target[i],s1,s2,s1.target[order(s1.target[,'mutationID']),],s1.s2.target.cov[order(s1.s2.target.cov[,'mutationID']),])
      write.table(out,quote = F,col.names=F, row.names = F,sep="\t",outfile1 ,append=T )
      out <- cbind(dupl.target[i],s1,s2,s2.target[order(s2.target[,'mutationID']),],s2.s1.target.cov[order(s2.s1.target.cov[,'mutationID']),])
      write.table(out,quote = F,col.names=F, row.names =F, sep="\t",outfile2 ,append=T  )

      out <- c(dupl.target[i],s1,s2, s1.s2.target.only.vaf, s1.s2.target.only.cov)
      write.table(out,quote = F, col.names=F,row.names = F,sep="\t",outfile ,append=T )

        #print("MY method")
        n.rep <- dim(s1.target)[1]
        s1.target$dist <-  sqrt( (s1.target$source/s1.target$target - rep(cont1/100,n.rep))^2 * log10(s1.target.min.cov))
        n.rep <- dim(s2.target)[1]
        s2.target$dist <-  sqrt( (s2.target$source/s2.target$target - rep(cont2/100,n.rep ))^2 * log10(s2.target.min.cov))

        row.min.cov <-current.data[which(current.data$target_cov > minCov & current.data$s1_cov >minCov & current.data$s2_cov> minCov),'mutationID']

        n.targetOnly = length(row.min.cov)
        n.rep <-n.targetOnly
        d1 <-  sqrt( (s1.s2.target.only.vaf[s1.s2.target.only.vaf$mutationID %in% row.min.cov,2]/s1.s2.target.only.vaf[s1.s2.target.only.vaf$mutationID %in% row.min.cov,4] - rep(cont2/100,n.rep ))^2 * log10(min.cov.s1.target.only[s1.s2.target.only.vaf$mutationID %in% row.min.cov]))
        d2 <-  sqrt( (s1.s2.target.only.vaf[s1.s2.target.only.vaf$mutationID %in% row.min.cov,3]/s1.s2.target.only.vaf[s1.s2.target.only.vaf$mutationID %in% row.min.cov,4] - rep(cont2/100,n.rep ))^2 * log10(min.cov.s2.target.only[s1.s2.target.only.vaf$mutationID %in% row.min.cov]))


      if(length(multi.source) ==0){
        weight.diff <- mean(s1.target$dist) + mean(d1)
        multi.source <- c(s1,s1.diff.median ,s1.diff.mean,s1.diff.sd,s2,dupl.target[i],dim(s1_tgt)[1],dim(both_in_tgt)[1],n.targetOnly,weight.diff,cont1,s1.median.cov,s1.min.cov,s1.diff.max,median(s1.target$ratio),mean(d1),mean(s1.target$dist))
        weight.diff <- mean(s2.target$dist) + mean(d2)
        multi.source <- rbind(multi.source,c(s2,s2.diff.median ,s2.diff.mean,s2.diff.sd,s1,dupl.target[i],dim(s2_tgt)[1],dim(both_in_tgt)[1],n.targetOnly,weight.diff,cont2,s2.median.cov,s2.min.cov,s2.diff.max,median(s2.target$ratio),mean(d2),mean(s2.target$dist)))
      }else{
        weight.diff <- mean(s1.target$dist) + mean(d1)
        multi.source <- rbind(multi.source,c(s1,s1.diff.median ,s1.diff.mean,s1.diff.sd,s2,dupl.target[i],dim(s1_tgt)[1],dim(both_in_tgt)[1],n.targetOnly,weight.diff,cont1,s1.median.cov,s1.min.cov,s1.diff.max,median(s1.target$ratio),mean(d1),mean(s1.target$dist)))
        weight.diff <- mean(s2.target$dist) + mean(d2)
        multi.source <- rbind(multi.source,c(s2,s2.diff.median ,s2.diff.mean,s2.diff.sd,s1,dupl.target[i],dim(s2_tgt)[1],dim(both_in_tgt)[1],n.targetOnly,weight.diff,cont2,s2.median.cov,s2.min.cov,s2.diff.max,median(s2.target$ratio),mean(d2),mean(s2.target$dist)))
      }
    }else if(nrow(current) >2){
      source.id <- current[,"source"]
      pair.file <- combn(source.id, 2)
      n.pair<- dim(pair.file)[2]
      discharge.source <- data.frame()
      store.source <- data.frame()

      for(j in 1:n.pair){
        s1 <- pair.file[1,j]
        s2 <- pair.file[2,j]

          f1 <- paste0(output_path, "/tmp/",dupl.target[i],"_",pair.file[1,j],"_",pair.file[2,j],"_present_absent_info.txt")
          f2 <- paste0(output_path, "/tmp/",dupl.target[i],"_",pair.file[2,j],"_",pair.file[1,j],"_present_absent_info.txt")
          if(file.exists(f1)){
            indata <- as.data.frame(read.table(f1))
          }else if(file.exists(f2)){
            indata <- as.data.frame(read.table(f2))
          }else{
            print ("Check file!")
          }

          # uniq in s1
          s1_tgt = as.matrix(rownames(indata)[indata$set=="present" & indata$sampleID==s1])
          # uniq in s2
          s2_tgt = as.matrix(rownames(indata)[indata$set=="present" & indata$sampleID==s2])
          # both
          both_in_tgt = as.matrix(rownames(indata)[indata$set == "present" & indata$sampleID =="both"])

          s1_tgt_both <- rbind(s1_tgt,both_in_tgt)
          s2_tgt_both <- rbind(s2_tgt,both_in_tgt)
          s1_s2_tgt_both <- rbind(s1_tgt,s2_tgt,both_in_tgt)

          col.id <- c("mutationID",s1)
          s1.vaf <- VAFdata[VAFdata$mutationID %in% s1_tgt_both, col.id]
          col.id <- c("mutationID",s2)
          s2.vaf <- VAFdata[VAFdata$mutationID %in% s2_tgt_both, col.id]
          col.id <- c("mutationID",dupl.target[i])
          s1.target.vaf <- VAFdata[VAFdata$mutationID %in% s1_tgt_both,  col.id]
          s2.target.vaf <- VAFdata[VAFdata$mutationID %in% s2_tgt_both,  col.id]

          s1.target <- cbind(s1.vaf[order(s1.vaf$mutationID),], s1.target.vaf[order(s1.target.vaf$mutationID),2])
          colnames(s1.target)<-c("mutationID","source","target")
          s1.target$diff <- s1.target$source - s1.target$target
          s1.target$ratio <- s1.target$source / s1.target$target

          s1.diff.median <- mean(s1.target$diff)
          s1.diff.sd <- sd(s1.target$diff)
          s1.diff.mean <- mean(s1.target$diff)
          s1.diff.max <- max(s1.target$diff)

          s2.target <- cbind(s2.vaf[order(s2.vaf$mutationID),], s2.target.vaf[order(s2.target.vaf$mutationID),2])
          colnames(s2.target)<-c("mutationID","source","target")
          s2.target$diff <- s2.target$source- s2.target$target
          s2.target$ratio <- s2.target$source/ s2.target$target

          s2.diff.median <- median(s2.target$diff)
          s2.diff.sd <- sd(s2.target$diff)
          s2.diff.mean <- mean(s2.target$diff)
          s2.diff.max <- max(s2.target$diff)
          cont1 <- as.numeric(current[which(current[,'source'] == s1),"predicted_contamination_perc"])
          cont2 <- as.numeric(current[which(current[,'source'] == s2),"predicted_contamination_perc"])

          col.id <- c("mutationID",s1,s2,dupl.target[i])
          s1.s2.target.cov <- VAFcov[VAFcov$mutationID %in% s1_tgt_both,  col.id]
          s2.s1.target.cov <- VAFcov[VAFcov$mutationID %in% s2_tgt_both,  col.id]

          col.id <- s1
          s1.median.cov <- median(s1.s2.target.cov[,col.id])
          s1.min.cov <- min(s1.s2.target.cov[,col.id])
          col.id <- s2
          s2.median.cov  <- median(s2.s1.target.cov[,col.id])
          s2.min.cov  <- min(s2.s1.target.cov[,col.id])

          col.id <- c(s1,dupl.target[i])
          s1.target.min.cov <- apply(s1.s2.target.cov[,col.id],1, function(each_row){
            min(each_row)
          })
          s1.target.min.cov <- as.matrix(unlist(s1.target.min.cov))

          col.id <- c(s2,dupl.target[i])
          s2.target.min.cov <- apply(s2.s1.target.cov[,col.id],1, function(each_row){
            min(each_row)
          })
          s2.target.min.cov <- as.matrix(unlist(s2.target.min.cov))

          ## in target not in s1|s2
          col.id <- c("mutationID",s1,s2,dupl.target[i])
          current.vaf.data <- VAFdata[,col.id]
          current.cov.data <- VAFcov[,col.id]
          current.data <- cbind(current.vaf.data,current.cov.data)
          colnames(current.data)<- c("mutationID","s1","s2","target", "mutationID2","s1_cov","s2_cov","target_cov")

          current.data <- subset(current.data,target>VAF_cutoff&target>VAF_ignore)

          exclude.mut <- rbind(s1_tgt,s2_tgt,both_in_tgt)
          current.data <- current.data[!current.data$mutationID %in% exclude.mut,]

          target_only <- current.data$mutationID
          s1.s2.target.only.vaf <-VAFdata[VAFdata$mutationID %in% target_only,  col.id]
          s1.s2.target.only.cov <-VAFcov[VAFcov$mutationID %in% target_only,  col.id]

          min.cov.s1.target.only <- pmin(s1.s2.target.only.cov[,2],s1.s2.target.only.cov[,4]) ## snps in target not in s1/s2
          min.cov.s2.target.only <- pmin(s1.s2.target.only.cov[,3],s1.s2.target.only.cov[,4])

          ## log
          out <- cbind(dupl.target[i],s1,s2,s1.target[order(s1.target[,'mutationID']),],s1.s2.target.cov[order(s1.s2.target.cov[,'mutationID']),])
          write.table(out,quote = F,col.names=F, row.names = F,sep="\t",outfile1 ,append=T )
          out <- cbind(dupl.target[i],s1,s2,s2.target[order(s2.target[,'mutationID']),],s2.s1.target.cov[order(s2.s1.target.cov[,'mutationID']),])
          write.table(out,quote = F,col.names=F, row.names =F, sep="\t",outfile2 ,append=T  )

          out <- c(dupl.target[i],s1,s2, s1.s2.target.only.vaf, s1.s2.target.only.cov)
          write.table(out,quote = F, col.names=F,row.names = F,sep="\t",outfile ,append=T )

          n.rep <- dim(s1.target)[1]
          s1.target$dist <-  sqrt( (s1.target$source/s1.target$target - rep(cont1/100,n.rep))^2 * log10(s1.target.min.cov))
          n.rep <- dim(s2.target)[1]
          s2.target$dist <-  sqrt( (s2.target$source/s2.target$target - rep(cont2/100,n.rep))^2 * log10(s2.target.min.cov))

          row.min.cov <-current.data[which(current.data$target_cov > minCov & current.data$s1_cov >minCov & current.data$s2_cov> minCov),'mutationID']
          n.targetOnly <- length(row.min.cov)
          n.rep <- n.targetOnly
          d1 <-  sqrt( (s1.s2.target.only.vaf[s1.s2.target.only.vaf$mutationID %in% row.min.cov,2]/s1.s2.target.only.vaf[,4] - rep(cont2/100,n.rep ))^2 * log10(min.cov.s1.target.only[s1.s2.target.only.vaf$mutationID %in% row.min.cov]))
          d2 <-  sqrt( (s1.s2.target.only.vaf[s1.s2.target.only.vaf$mutationID %in% row.min.cov,3]/s1.s2.target.only.vaf[,4] - rep(cont2/100,n.rep ))^2 * log10(min.cov.s2.target.only[s1.s2.target.only.vaf$mutationID %in% row.min.cov]))

          if(length(multi.source) ==0){
            weight.diff <- mean(s1.target$dist) + mean(d1)
            multi.source <- c(s1,s1.diff.median ,s1.diff.mean,s1.diff.sd,s2,dupl.target[i],dim(s1_tgt)[1],dim(both_in_tgt)[1],n.targetOnly,weight.diff,cont1,s1.median.cov,s1.min.cov,s1.diff.max,median(s1.target$ratio),mean(d1),mean(s1.target$dist))
            weight.diff <- mean(s2.target$dist) + mean(d2)
            multi.source <- rbind(multi.source,c(s2,s2.diff.median ,s2.diff.mean,s2.diff.sd,s1,dupl.target[i],dim(s2_tgt)[1],dim(both_in_tgt)[1],n.targetOnly,weight.diff,cont2,s2.median.cov,s2.min.cov,s1.diff.max,median(s2.target$ratio),mean(d2),mean(s2.target$dist)))
          }else{
            weight.diff <- mean(s1.target$dist) + mean(d1)
            multi.source <- rbind(multi.source,c(s1,s1.diff.median ,s1.diff.mean,s1.diff.sd,s2,dupl.target[i],dim(s1_tgt)[1],dim(both_in_tgt)[1],n.targetOnly,weight.diff,cont1,s1.median.cov,s1.min.cov,s1.diff.max,median(s1.target$ratio),mean(d1),mean(s1.target$dist)))
            weight.diff <- mean(s2.target$dist) + mean(d2)
            multi.source <- rbind(multi.source,c(s2,s2.diff.median ,s2.diff.mean,s2.diff.sd,s1,dupl.target[i],dim(s2_tgt)[1],dim(both_in_tgt)[1],n.targetOnly,weight.diff,cont2,s2.median.cov,s2.min.cov,s2.diff.max,median(s2.target$ratio),mean(d2),mean(s2.target$dist)))
          }
        } # each j pair
    } # > 2 source
    colnames(multi.source) <- c("source","median_diff","mean_diff","sd_diff","compare_source", "target","num_snps","num_both","num_targetOnly", "weight_diff","cont","cov_median","cov_min","vaf_diff_max","median_VAFratio", "dist_targetOnly","dist_source_target")

    col.id <- c("source","weight_diff")
    uniq.source <- as.matrix(unique(multi.source[,col.id]))
    remove.source <- uniq.source[-which.min(as.numeric(uniq.source[,'weight_diff'])),'source']

    ## remove source should not be the one with max weight diff
    if(any(remove.source %in% uniq.source[which.min(as.numeric(uniq.source[,'weight_diff'])),'source'])){
      remove.source <- remove.source[!remove.source %in% uniq.source[which.min(as.numeric(uniq.source[,'weight_diff'])),'source']]
    }

    rm.row <- which(match(process.predict.cont[,"source"], remove.source)!="NA" & process.predict.cont[,"target"] == dupl.target[i])
    if(length(rm.row)>0){
      process.predict.cont <- process.predict.cont[-rm.row,]

    }
  if(first){
    outfile <- paste0(output_path,"/tmp/log_dist.txt")
    write.table(multi.source,quote=F,outfile,append=F,row.names=F)
    first = FALSE
  }else{
    write.table(multi.source,quote=F,outfile,append=T,row.names=F)  }
  }# each dupl target

  outfile<- paste0(output_path,"/contPredict.txt")
  write.table(process.predict.cont, sep='\t',row.names = F, quote=F,outfile)
process.predict.cont
}
