#! /usr/bin/env Rscript
suppressWarnings(suppressMessages(library(DNAcopy)))
suppressWarnings(suppressMessages(library(pastecs)))


args = R.utils::commandArgs(asValue = T)
datadir <- args$datadir
workdir <- args$workdir
seriesName <- args$seriesName
force <- as.logical(args$force)
genome <- args$genome

source(file.path(workdir,'bin','utils.R'))

if (!genome %in% c('hg19','hg38')) stop("genome version is invalid")

cids <- list.files(file.path(datadir, "processed",seriesName))

if (length(cids) == 0) stop("`seriesName` is invalid")

localProcessPath <- file.path(datadir,"processed",seriesName)

targetcids <- checkFileIncomplete(force, localProcessPath,'segments,fracb.tsv', cids)

probefileName <- ifelse(genome == 'hg38','probes,fracb,hg38.tsv','probes,fracb.tsv')

if (length(targetcids) > 0) {
  for (cid in targetcids) {
    if (!file.exists(sprintf('%s/%s/%s',localProcessPath,cid,probefileName))) stop(paste('fracb probe file not found for',seriesName, cid))
    
    fseg <- file.path(localProcessPath,cid,'fracbseg.tsv')
    cat("ID","chrom","loc.start", "loc.end", "num.mark", "seg.mean","seg.sd","seg.median", 
        "seg.mad\n",sep="\t",file=fseg,append = F)
    
    ffinalseg <- file.path(localProcessPath,cid,"segments,fracb.tsv")
    cat("ID","chr","loc.start","loc.end","fracB\n",sep="\t",file=ffinalseg,append=F)
    
    alldata <- read.table(sprintf('%s/%s/%s',
                                  localProcessPath,cid,probefileName),header=T,stringsAsFactors=F)
    
    alldata <-  alldata[!alldata[,2] %in% c('X','Y',23,24),]
    
    for (chrname in 1:22){
      # take fracb probe data in single chr
      data <- alldata[alldata[,2]==chrname,]
      #################### segment in major allele frequency############## 
      if ('BASEPOS' %in% colnames(data) & 'VALUE' %in% colnames(data)){
        posn <-data$BASEPOS 
        baf1 <-data$VALUE
      } else {
        posn <- data[,3]
        baf1<- data[,4]
      }
      
      ### mBAF mirrored
      mbaf <- unlist(sapply(baf1, function(x) {
        if (is.na(x) == FALSE) {
          if (x < 0.5) x <- 1-x
          if (x >= 0.5) x<-x
        }
        else x<-x
      }))
      mbaf[abs(mbaf-0.5) >0.35] <- NA
      yy <- mbaf
      
      # segment to estimate both LOH and BAF
      cna1 <- CNA(as.matrix(yy),
                  chrom=rep(chrname, length(yy)),
                  maploc=posn,
                  data.type="logratio",
                  sampleid=cid)
      
      smoo1 <- smooth.CNA(cna1, smooth.region=10,smooth.SD.scale=4)
      
      
      seg1 <- segment(smoo1, min.width=5, verbose=0)
      ss1 <- segments.summary(seg1)
      write.table(ss1, file=fseg, sep="\t", quote=FALSE,
                  append = T, row.names=FALSE, col.names=F)
      #################### identify frac peaks in the corresponding segment region############## 
      ss1_sel <- ss1[ss1$num.mark > 0,]
      for (j in 1:nrow(ss1_sel)){
        range <- c(ss1_sel[j,"loc.start"], ss1_sel[j,"loc.end"])
        id <- which(data$BASEPOS>=range[1] & data$BASEPOS< range[2])
        subfracb <- data[id,]
        if (nrow(subfracb) <= 5) next
        dens <- density(subfracb$VALUE)
        y <- dens$y
        x <- dens$x
        tp<-turnpoints(y)
        peak <- x[tp$peaks]
        peak <- sort(c(peak,1-peak))
        peak <- rm.near.dy(peak)
        peak <- signif(peak, digits = 4)
        if (length(peak) == 0) next
        for (line in 1:length(peak)){
          write.table(data.frame(ID=cid,chr=chrname,loc.start=range[1],loc.end=range[2],fracB=peak[line]), 
                      file=ffinalseg, sep="\t", quote=FALSE,
                      append = T, row.names=FALSE, col.names=F)
        }
      }
    }
  }
    
} 

log_path <- file.path(datadir,"processed",'logs')
dir.create(log_path, showWarnings = F, recursive = TRUE)

segfiles <- list.files(localProcessPath,recursive = T, full.names =T )
if (sum(grepl("segments,fracb.tsv",segfiles)) == length(cids)){
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"BAF-segmentation",seriesName,"complete") 
} else{
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"BAF-segmentation",seriesName,"complete but file numbers aberrant") 
}

log_con <- file(file.path(log_path,paste0("process_",format(Sys.time(), "%y-%m-%d"),".log")),"a")  
writeLines(log,log_con)








