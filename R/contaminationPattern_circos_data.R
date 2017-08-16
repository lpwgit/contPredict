#' Generate circos data base on contamination patterns
#'
#'This function generates data base on contamination patterns as input to circos plot
#'
#' @param contaminationPatternfile character, file name and its path
#'
#' @usage contaminationPattern_circos_data (contaminationPatternFile)
#'
#' @return NA
#'
#' @references
#' {TBA}
#'
#' @export
#'
contaminationPattern_circos_data <- function (contaminationPatternFile){
  if(!file.exists(contaminationPatternFile)){
    print(paste(contaminationPatternFile,"not found!"))
  }else{
    contamination.pattern <- read.table(contaminationPatternFile,header=TRUE)
    colnames(contamination.pattern) <-c("source","target","perc")
    n.pattern <- dim(contamination.pattern)[1]

    #print(n.pattern)
    #print(contamination.pattern)

    outfile<- paste0(contaminationPatternFile,".circosData")
    if(file.exists(outfile)){
      file.remove(outfile)
    }
    relation <-rep("10",n.pattern)

    ## bothways contamination
    circos.data <- data.frame(source=contamination.pattern$source,target=contamination.pattern$target,link=round(contamination.pattern$perc/10,3), rel=relation)
    write.table(circos.data,sep='\t',quote=F,outfile,row.names = F)
    #paste0("circos data:  ", outfile)
  }
}
