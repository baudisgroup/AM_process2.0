#! /usr/bin/env Rscript
suppressWarnings(suppressMessages(library(matrixStats)))
args = R.utils::commandArgs(asValue = T)

datadir <- args$datadir
seriesName <- args$seriesName
workdir <- args$workdir
adj_probe <- as.logical(args$adjust_probe)
genome <- args$genome

source(file.path(workdir,'bin','utils.R'))

cids <- list.files(file.path(datadir, "processed",seriesName))

log_path <- file.path(datadir,"processed",'logs')
dir.create(log_path, showWarnings = F, recursive = TRUE)
log_con <- file(file.path(log_path,paste0("process_",format(Sys.time(), "%y-%m-%d"),".log")),"a")  

for (cid in cids){
  adjusted_log <- adjustMedian(datadir,seriesName,cid,adj_probe,genome)
  writeLines(adjusted_log,log_con)
}

