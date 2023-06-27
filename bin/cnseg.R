#! /usr/bin/env Rscript
suppressWarnings(suppressMessages(library(DNAcopy)))

args = R.utils::commandArgs(asValue = T)
datadir <- args$datadir
workdir <- args$workdir
seriesName <- args$seriesName
force <- as.logical(args$force)
undosd <- as.numeric(args$undosd)
genome <- args$genome

source(file.path(workdir,'bin','utils.R'))

if (!genome %in% c('hg19','hg38')) stop("genome version is invalid")
  
cids <- list.files(file.path(datadir, "processed",seriesName))

if (length(cids) == 0) stop("`seriesName` is invalid")

localProcessPath <- file.path(datadir,"processed",seriesName)

targetcids <- checkFileIncomplete(force, localProcessPath,'segments,cn.tsv', cids)

probefileName <- ifelse(genome == 'hg38','probes,cn,hg38.tsv','probes,cn.tsv')

if (length(targetcids) > 0) {
  for (cid in targetcids) {
    if (!file.exists(sprintf('%s/%s/%s',localProcessPath,cid,probefileName))) stop(paste('cn probe file not found for',seriesName, cid))
    
    fseg <- file.path(localProcessPath,cid,'segments,cn.tsv')
    cat("sample_id","chromosome","start", "end", "value", "probes\n",sep="\t",file=fseg,append = F)

    alldata <- read.table(sprintf('%s/%s/%s',
                                  localProcessPath,cid,probefileName),header=T,stringsAsFactors=F)
    
    # remove sex chromsome
    alldata <-  alldata[!alldata[,2] %in% c('X','Y',23,24),]
    
    for (chrname in 1:22){
      data <- alldata[alldata[,2]==chrname,]
      
      if ('BASEPOS' %in% colnames(data) & 'VALUE' %in% colnames(data)){
        posn <-as.numeric(data$BASEPOS)
        cnvalue <-as.numeric(data$VALUE)
      } else {
        posn <- as.numeric(data[,3])
        cnvalue<- as.numeric(data[,4])
      }
      
      cna1 <- CNA(genomdat=cnvalue,
                  chrom=rep(chrname, length(cnvalue)),
                  maploc=posn,
                  data.type="logratio",
                  sampleid=cid)
      smoo1 <- smooth.CNA(cna1, smooth.region=5,smooth.SD.scale=2)
      
      seg1 <- segment(smoo1, min.width=5, verbose=0, undo.splits='sdundo',undo.SD=undosd)
      ss1 <- segments.summary(seg1)[c(1:4,6,5)]
      write.table(ss1, file=fseg, sep="\t", quote=FALSE,
                  append = T, row.names=FALSE, col.names=F)
    }
  }
  
} 

log_path <- file.path(datadir,"processed",'logs')
dir.create(log_path, showWarnings = F, recursive = TRUE)

segfiles <- list.files(localProcessPath,recursive = T, full.names =T )
if (sum(grepl("segments,cn.tsv",segfiles)) == length(cids)){
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"CN-segmentation",seriesName,"complete") 
} else{
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"CN-segmentation",seriesName,"complete but file numbers aberrant") 
}

log_con <- file(file.path(log_path,paste0("process_",format(Sys.time(), "%y-%m-%d"),".log")),"a")  
writeLines(log,log_con)


  
  
  
  
  
  