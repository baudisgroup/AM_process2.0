#! /usr/bin/env Rscript
suppressWarnings(suppressMessages(library(genlasso)))
suppressWarnings(suppressMessages(library(LabelSeg)))
args = R.utils::commandArgs(asValue = T)

datadir <- args$datadir
workdir <- args$workdir
seriesName <- args$seriesName
genome <- args$genome

files.sources <- list.files(path = file.path(workdir,'bin','visual_seg','scripts'),recursive = T,full.names = T)
invisible(sapply(files.sources,source))
load(file.path(workdir,'bin','visual_seg','sysdata.rda'))
source(file.path(workdir,'bin','visual_seg','visual.seg.R'))

source(file.path(workdir,'bin','utils.R'))

cids <- list.files(file.path(datadir, "processed",seriesName))

if (length(cids) == 0) stop("`arrayName` or `seriesName` is invalid")

localProcessPath <- file.path(datadir,"processed",seriesName)
reportpath <- file.path(datadir,"processed",'data_quality_report','cnseg','calibration',paste0(seriesName,'.txt'))
if (!file.exists(reportpath)) stop(paste('No data quality report for',seriesName))
report <- read.table(reportpath,sep = '\t',header=T)

noisy_idx <- which(report$highCNA_cov > 0.1 | report$segment_num > 1000)
fail_label_idx <- which(report$lowCNA_cov == 0 & report$highCNA_cov == 0 & report$normal_cov == 0)

targetcids <- cids[union(fail_label_idx,noisy_idx)]

for (cid in targetcids){
  cidseg <- stepFilter(datadir,seriesName,cid,lmd=1,genome)
  cidseg <- cidseg[!cidseg[,2] %in%  c(23,'X',24,'Y'),]
  # label new segment
  labelcidseg <- cidseg[,c(1,2,3,4,6,5)]
  labelcidseg <- labelseg(data=labelcidseg,genome = genome)
  cov <- seg_cov(labelcidseg)
  report <- update_seg(series_path=localProcessPath,experiment=cid,relabelSeg=labelcidseg,report=report,score=cov,problem='resegmentation',genome=genome,mergeseg=T)
}

write.table(report, file=reportpath,sep = '\t',quote = F,row.names = F)

# write log
log_path <- file.path(datadir,"processed",'logs')
dir.create(log_path, showWarnings = F, recursive = TRUE)
log_con <- file(file.path(log_path,paste0("process_",format(Sys.time(), "%y-%m-%d"),".log")),"a")  

log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"resegmentation complete","for",seriesName) 
writeLines(log,log_con)





