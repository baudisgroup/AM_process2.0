#! /usr/bin/env Rscript
suppressWarnings(suppressMessages(library(LabelSeg)))
args = R.utils::commandArgs(asValue = T)

datadir <- args$datadir
seriesName <- args$seriesName
workdir <- args$workdir
genome <- args$genome

files.sources <- list.files(path = file.path(workdir,'bin','visual_seg','scripts'),recursive = T,full.names = T)
invisible(sapply(files.sources,source))
load(file.path(workdir,'bin','visual_seg','sysdata.rda'))
source(file.path(workdir,'bin','visual_seg','visual.seg.R'))

source(file.path(workdir,'bin','utils.R'))

localProcessPath <- file.path(datadir,"processed",seriesName)
cids <- list.files(localProcessPath)
if (length(cids) == 0) stop("`seriesName` is invalid")

reportpath <- file.path(datadir,"processed",'data_quality_report','cnseg','calibration')
dir.create(reportpath,showWarnings=F,recursive = T)

report <- list()
for (cid in cids) {
  samplepath <- sprintf('%s/%s/segments,cn,%s.tsv',localProcessPath,cid,genome)
  if (!file.exists(samplepath)) stop(paste('cnseg file not found for',seriesName, cid))
    
  cidseg <- read.table(samplepath,sep = '\t',header=T)

  if (any(is.na(cidseg[,c(1:6)])) | dim(cidseg)[1] == 0 | any(cidseg[,4] < cidseg[,3])){
      stop(paste("cnseg is invalid in",seriesName, cid))
  }
    
    # remove sex chroms 
  cidseg <- cidseg[!cidseg[,2] %in%  c(23,'X',24,'Y'),]
    
    # label segment
  labelcidseg <- cidseg[,c(1,2,3,4,6,5)]
  labelcidseg <- labelseg(data=labelcidseg,genome = genome)
    
  cidsegprobe <- cidseg[,6]
    # check if files have abberant probe numbers
  if (any(cidsegprobe == 0)){stop(paste('probe number of cnseg is invalid in',seriesName, cid))}
    
    
  lrrsegsd <- round(sd(rep(cidseg[,5],abs(cidsegprobe))),3)
  segNum <- dim(cidseg)[1]
  cov <- seg_cov(labelcidseg)
  lowCov <- cov$lowCNA_cov
  highCov <-  cov$highCNA_cov
  normalCov <- cov$normal_cov
  lowCov_ratio <- cov$lowCNA_ratio
  
  plotfilename <- paste0('segments,cn,',genome,'.pdf')  
  tablefilename <- paste0('labelsegments,cn,',genome,'.tsv')
  if (!is.null(labelcidseg)){
      # plot and save image of labelled segments
    plot_segment_label(filepath = file.path(localProcessPath,cid),filename=plotfilename,data = labelcidseg,assembly = genome)
      # write labelled segment file
    write.table(labelcidseg, file=file.path(localProcessPath,cid,tablefilename),sep = '\t',row.names = F,quote=F)
    note <- "initial"
  } else{
    plot_segment_label(filepath = file.path(localProcessPath,cid),filename= plotfilename,data = cidseg[,c(1,2,3,4,6,5)],assembly = genome,no_label = T)
    note <- "failed-to-label"
  }
    
  report[[cid]] <- data.frame(series=seriesName, expeiment=cid,segment_num=segNum,
                              LLR_segsd = lrrsegsd, normal_cov= normalCov, lowCNA_cov=lowCov, lowCNA_ratio=lowCov_ratio , highCNA_cov= highCov, note=note)
}

report <- do.call(rbind, report)
  
# write report 
write.table(report, file=file.path(reportpath,paste0(seriesName,'.txt')),sep = '\t',quote = F,row.names = F)

# check and write log
log_path <- file.path(datadir,"processed",'logs')
dir.create(log_path, showWarnings = F, recursive = TRUE)

segfiles <- list.files(localProcessPath,recursive = T, full.names =T )
if (sum(grepl(plotfilename,segfiles)) == length(cids)){
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"labelseg",seriesName,"complete") 
} else{
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"labelseg",seriesName,"complete but file numbers aberrant") 
}

log_con <- file(file.path(log_path,paste0("process_",format(Sys.time(), "%y-%m-%d"),".log")),"a")  
writeLines(log,log_con)






