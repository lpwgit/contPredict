#'Plot circos with links to show contaminated pairwise sample
#'
#'This function reads in circos plotdata and creates circos plot
#'
#' @param file character, circos data file name with its path
#' @param R integer, circle radius
#' @param W integer, segment width
#' @param plotsize integer, figure size
#' @param titleStr character, title of circos plot
#' @param seg.lab.size numeric, segment label font size
#' @param fig.file character, circos figure file name
#' @param contaminationOnly logical, only plot contaminated pairs? (default: TRUE)
#' @param sameSubjectOnly logical, only plot same subject pairs? (default: FALSE)
#' @param allSamples character vector, sample ids
#' @param output_path character, directory
#' @usage plot_circos_link (file,R=300,W=30,plotsize=800,titleStr="",seg.lab.size=1,fig.file="circos.png", contaminatedOnly=TRUE,sameSubjectOnly=FALSE,allSamples=sampleID)
#'
#' @return .tif figure at output directory
#'
#' @references
#' {TBA}
#'
#' @export
#'
plot_circos_link <- function (circos.indata,R=300,W=30,plotsize=800,titleStr="",seg.lab.size=1,fig.file="circos.pdf", contaminatedOnly=TRUE,sameSubjectOnly=FALSE,allSamples=sampleID,output_path){
  outfile<-paste0(output_path,"/allPairRelation.txt")
  write.table(circos.indata[,-3],outfile,sep='\t',quote=F,row.names = F)
  #circos.indata <- read.table(file,header=T,colClasses = c("character","character","numeric","character"))
  sample.id <- sort(unique(c(circos.indata$source,circos.indata$target,allSamples)))
  n.sample <- length(sample.id)

  ## segment data
  seg.pos <- c(0,10)
  dummy.link <-rep("NA",n.sample)
  dummy.link.pg <-rep("NA",n.sample)
  seg.data <- as.data.frame(cbind(seg.name=sample.id,seg.Start=seg.pos[1],seg.End=seg.pos[2],the.V=dummy.link,NO=dummy.link.pg),stringsAsFactors=FALSE)

  ## link data
  link.colors <- as.vector(circos.indata$rel)
  link.colors[which(link.colors=="00")] <- "grey"
  link.colors[which(link.colors=="-1")] <- "white"
  link.colors[which(link.colors=="10")] <- "blue"
  link.colors[which(link.colors=="11")] <- "red"
  link.colors[which(link.colors=="01")] <- "green"

  mid.point <- round(mean(seg.pos),digits=0)

  ## create 2 links for bothway contamination
  bothway.data <- circos.indata[which(circos.indata$rel =="11"),]
  if(nrow(bothway.data)==0){
    all.link.data <-as.data.frame(cbind(seg1=circos.indata$source,po1=mid.point,name1=circos.indata$source,seg2=circos.indata$target,po2=mid.point,name2=circos.indata$target,weight=circos.indata$link,rel=circos.indata$rel,linkColor = link.colors),stringsAsFactors=FALSE)
    double.line.flag <- rep(0,dim(all.link.data)[1])
  }else{
    tmp.p1 <- as.data.frame(cbind(seg1=circos.indata$source,po1=mid.point,name1=circos.indata$source,seg2=circos.indata$target,po2=mid.point,name2=circos.indata$target,weight=circos.indata$link,rel=circos.indata$rel,linkColor = link.colors),stringsAsFactors=FALSE)
    p1 <- tmp.p1[tmp.p1$rel !='11',]
    p3 <- as.data.frame(cbind(seg1=bothway.data$source,po1=mid.point,name1=bothway.data$source,seg2=bothway.data$target,po2=mid.point,name2=bothway.data$target,weight=bothway.data$link,rel="10",linkColor = "green"),stringsAsFactors=FALSE)
    p2 <- as.data.frame(cbind(seg1=bothway.data$source,po1=mid.point-1,name1=bothway.data$source,seg2=bothway.data$target,po2=mid.point-1,name2=bothway.data$target,weight=bothway.data$link,rel="01",linkColor = "green"),stringsAsFactors=FALSE)
    all.link.data <- rbind(p1,p2,p3)
    double.line.flag <- c(rep(0,dim(p1)[1]),rep(1,dim(p2)[1]),rep(1,dim(p3)[1]))
  }
  all.link.data$flag <- double.line.flag

  ## select samples
  seg.name = sample.id
  seg.num <- length(sample.id)
  plot.db <- segAnglePo(seg.data, seg=seg.name)
  seg.colors <- rainbow(seg.num, alpha=0.5)
  seg.info <- as.data.frame(cbind(id=sample.id, color=seg.colors),stringsAsFactors=FALSE)

  ## set link color base on contamination source
  source.color <- apply(all.link.data,1, function(i){
    if(i['rel'] =="01" ){
      seg.info[seg.info$id == i['name2'],'color']
    }else if (i['rel'] =="10"){
      seg.info[seg.info$id == i['name1'],'color']

    }else if (i['rel'] =="11"){
      seg.info[seg.info$id == i['name1'],'color']
    }else{
      i['linkColor']
    }
  })

  # check multisource
  only.contaminated <- all.link.data[all.link.data$rel!="00",]
  if(dim(only.contaminated)[1] >0){
  multisource <- apply(only.contaminated,1, function(i){
    if(i['rel'] =="01" ){
      c(i['name2'],i['name1'],i['flag'])
    }else if (i['rel'] =="10"){
      c(i['name1'],i['name2'],i['flag'])
    }else if (i['rel'] =="11"){
      c(i['name1'],i['name2'],i['flag'])
    }else{
      c(i['name2'],i['name1'],i['flag'])
    }
  })

  multisource <-t(multisource)
  colnames(multisource) <- c("source","target", "flag")
  multisource.target <- multisource[duplicated(multisource[,'target']),'target']
}
  all.link.data$linkColorBySource <- as.character((unlist(source.color)))
  link.type <- as.matrix(table(all.link.data$rel))
  sameSubject.flag =FALSE
  x2y.flag = FALSE
  y2x.flag = FALSE
  bothway.flag = FALSE

  #only same subjects
  if(length(which(rownames(link.type) =="00") ==1) & !contaminatedOnly){
    sameSubject.link.data <- all.link.data[all.link.data$rel=="00",]
    sameSubject.flag =TRUE
    sameSubject_line_type=rep("l",dim(sameSubject.link.data)[1])
  }

  ## same subject combine with contamination
  if(length(which(rownames(link.type) =="00") ==1) & contaminatedOnly){
    sameSubject.link.data <- all.link.data[all.link.data$rel=="00",]
    sameSubject.flag =TRUE
    sameSubject_line_type=rep("l",dim(sameSubject.link.data)[1])
  }

  if(length(which(rownames(link.type) =="01") ==1)){
    x2y.link.data <- all.link.data[all.link.data$rel=="01",]
    x2y.flag = TRUE

    # check for multisource
    x2y_line_type=rep("l",dim(x2y.link.data)[1])
    if(any(x2y.link.data$seg1 %in% multisource.target)){
      #x2y_line_type[which(x2y.link.data$seg1 %in% multisource.target ==TRUE)]= "c" # come back later
      x2y_line_type=rep("l",dim(x2y.link.data)[1])
    }else{
      x2y_line_type=rep("l",dim(x2y.link.data)[1])
    }
  }
  if(length(which(rownames(link.type) =="10")==1)){
    y2x.link.data <- all.link.data[all.link.data$rel=="10",]
    y2x.flag = TRUE

    # check for multisource
    y2x_line_type=rep("l",dim(y2x.link.data)[1])
    if(any(y2x.link.data$seg2 %in% multisource.target)){
      #y2x_line_type[which(y2x.link.data$seg2 %in% multisource.target ==TRUE)]= "c" ## come back later
      y2x_line_type=rep("l",dim(y2x.link.data)[1])
    }else{
      y2x_line_type=rep("l",dim(y2x.link.data)[1])
    }
  }
  if(length(which(rownames(link.type) =="11")==1)){
    bothway.link.data <- all.link.data[all.link.data$rel=="11",]
    bothway.flag = TRUE

    # check for multisource
    bothway_line_type=rep("l",dim(bothway.link.data)[1])
    if(any(bothway.link.data$seg1 %in% multisource.target)){
      #bothway_line_type[which(bothway.link.data$seg1 %in% multisource.target ==TRUE)]= "c" # come back later
      bothway_line_type=rep("l",dim(bothway.link.data)[1])
    }else{
      bothway_line_type=rep("l",dim(bothway.link.data)[1])
    }
  }

  if(any(c( sameSubject.flag, x2y.flag, y2x.flag,bothway.flag))){
    hsize <-6
    wsize <- 6
    pdf(paste0(output_path,"/",fig.file), height = hsize, width = wsize)
    par(mar=c(0, 0, 2, 0))
    plot(c(1,plotsize), c(1,plotsize), type="n", axes=FALSE, xlab="", ylab="", main=titleStr)
    circos_link(R=R, cir=plot.db, type="chr",  col=seg.colors, print.chr.lab=TRUE, W=W, seg.lab.size = seg.lab.size)
    R1 <- R-20
    if(sameSubjectOnly  & !contaminatedOnly){
      circos_link(R=R1, cir=plot.db, W=W, mapping=sameSubject.link.data, type="link",  col=sameSubject.link.data$linkColorBySource,link.wd=sameSubject.link.data$weight,lower_arch=2,line_type=sameSubject_line_type)
    }
    if(sameSubjectOnly & contaminatedOnly){
      circos_link(R=R1, cir=plot.db, W=W, mapping=sameSubject.link.data, type="link",  col=sameSubject.link.data$linkColorBySource,link.wd=sameSubject.link.data$weight,lower_arch=2,line_type=sameSubject_line_type)
      if(x2y.flag ){
        #print(x2y.link.data)
        circos_link(R=R1, cir=plot.db, W=W, mapping=x2y.link.data, type="link",  col=x2y.link.data$linkColorBySource,link.wd = x2y.link.data$weight,line_type= x2y_line_type)
      }
      if(y2x.flag ){
        #print(y2x.link.data)
        circos_link(R=R1, cir=plot.db, W=W, mapping=y2x.link.data, type="link",  col=y2x.link.data$linkColorBySource,link.wd = y2x.link.data$weight,line_type =y2x_line_type)
      }
      if(bothway.flag ){
        circos_link(R=R1, cir=plot.db, W=W, mapping=bothway.link.data, type="link", col=bothway.link.data$linkColorBySource,link.wd = bothway.link.data$weight,line_type=bothway_line_type)
      }
    }
    if(x2y.flag & !sameSubjectOnly){
      #print(x2y.link.data)
      circos_link(R=R1, cir=plot.db, W=W, mapping=x2y.link.data, type="link",  col=x2y.link.data$linkColorBySource,link.wd = x2y.link.data$weight,line_type= x2y_line_type)
    }
    if(y2x.flag & !sameSubjectOnly){
      #print(y2x.link.data)
      circos_link(R=R1, cir=plot.db, W=W, mapping=y2x.link.data, type="link",  col=y2x.link.data$linkColorBySource,link.wd = y2x.link.data$weight,line_type =y2x_line_type)
    }
    if(bothway.flag & !sameSubjectOnly){
      circos_link(R=R1, cir=plot.db, W=W, mapping=bothway.link.data, type="link", col=bothway.link.data$linkColorBySource,link.wd = bothway.link.data$weight,line_type=bothway_line_type)
    }
    invisible(dev.off())
  }

}
