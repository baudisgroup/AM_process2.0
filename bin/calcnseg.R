#! /usr/bin/env Rscript
suppressWarnings(suppressMessages(library(mongolite)))
suppressWarnings(suppressMessages(library(LabelSeg)))
args = R.utils::commandArgs(asValue = T)

datadir <- args$datadir
seriesName <- args$seriesName
workdir <- args$workdir
genome <- args$genome
docker <- as.logical(args$docker) 

source(file.path(workdir,'bin','utils.R'))
files.sources <- list.files(path = file.path(workdir,'bin','visual_seg','scripts'),recursive = T,full.names = T)
invisible(sapply(files.sources,source))
load(file.path(workdir,'bin','visual_seg','sysdata.rda'))
source(file.path(workdir,'bin','visual_seg','visual.seg.R'))

if (docker){
  connect_url <- "mongodb://host.docker.internal:27017"
} else{
  connect_url <- "mongodb://localhost:27017"
}

biosamples_collection <- mongo(collection="biosamples",
                                          db="progenetix",
                                          url =connect_url)

meta <- get_meta(collection=biosamples_collection,series=seriesName)

localProcessPath <- file.path(datadir,"processed",seriesName)
reportpath <- file.path(datadir,"processed",'data_quality_report','cnseg',paste0(seriesName,'.txt'))

if (!file.exists(reportpath)) stop(paste('No data quality report for',seriesName))

report <- read.table(reportpath,sep = '\t',header=T)

# find profiles failed to label
fail_label_idx <- which(report$lowCNA_cov == 0 & report$highCNA_cov == 0 & report$normal_cov == 0)

# find noisy profiles 
noisy_idx <- which(report$highCNA_cov > 0.1 | report$segment_num > 1000)
noisy_idx <-  setdiff(noisy_idx ,fail_label_idx)

# find profiles with possible abnormal baseline 
h_shift_idx <- which(report$LLR_segsd > 0.35 & report$lowCNA_cov > 0.5 & report$lowCNA_ratio >=3)
h_shift_idx <- unique(c(h_shift_idx, which(report$lowCNA_cov > 0.25 & report$lowCNA_ratio >= 5)))
h_shift_idx <- setdiff(setdiff(h_shift_idx,fail_label_idx),noisy_idx)

l_shift_idx <- which(report$LLR_segsd > 0.35  & report$lowCNA_cov > 0.5 & report$lowCNA_ratio <= 1/3)
l_shift_idx <- unique(c(l_shift_idx, which(report$lowCNA_cov > 0.25 & report$lowCNA_ratio <= 0.2)))
l_shift_idx <- setdiff(setdiff(l_shift_idx,fail_label_idx),noisy_idx)

bins.lst <- generate_genomic_intervals(workdir=workdir,genome=genome)
labelcidseg <- get_labelseg(datadir=datadir,series=seriesName,genome=genome)
if (is.null(labelcidseg)){return()}

fit <- compute_fit(labelcidseg)
whole_shift <- "n"
# Data likelihood is affected by number of samples to be calibrated 
half_sample_len <- 0.5*length(unique(labelcidseg[,1]))
# avoid whole shift inconsistency
if (length(h_shift_idx) >= half_sample_len){
  compare <- compare_whole_shift(datadir,workdir,seriesName,meta,labelcidseg,h_shift_idx,l_shift_idx,baseshift = 'h',report,fit,bins.lst,genome)
  whole_shift <- compare$shift
}

if (length(l_shift_idx) >= half_sample_len){
  compare_2 <- compare_whole_shift(datadir,workdir,seriesName,meta,labelcidseg,l_shift_idx,h_shift_idx,baseshift = 'l',report,fit,bins.lst,genome)
  whole_shift <- compare_2$shift
  if (length(h_shift_idx) >= half_sample_len){
    if (compare_2$posterior < compare$posterior){
      whole_shift <- compare$shift
    }
  }
}

if (length(h_shift_idx) > 0){
  for (experiment in report$expeiment[h_shift_idx]){
    prior_code <- get_prior_code(meta,experiment)
    if (whole_shift == 'n'){
      report <- shift_baseline(datadir,workdir,seriesName,prior_code,experiment,report,bins.lst,fit,whole_shift = F,genome,shift="higher")
    } else if (whole_shift == 'h'){
      report <- shift_baseline(datadir,workdir,seriesName,prior_code,experiment,report,bins.lst,fit,whole_shift = T,genome,shift="higher")
    } else{
      idx <- which(report$expeiment == experiment)
      report$note[idx] <- paste0(report$note[idx],";shift-baseline-higher_reverse")
    }
  }
}

if (length(l_shift_idx) > 0){
  for (experiment in report$expeiment[l_shift_idx]){
    prior_code <- get_prior_code(meta,experiment)
    if (whole_shift == 'n'){
      report <- shift_baseline(datadir,workdir,seriesName,prior_code,experiment,report,bins.lst,fit,whole_shift = F,genome,shift="lower")
    } else if (whole_shift == 'l'){
      report <- shift_baseline(datadir,workdir,seriesName,prior_code,experiment,report,bins.lst,fit,whole_shift = T,genome,shift="lower")
    } else{
      idx <- which(report$expeiment == experiment)
      report$note[idx] <- paste0(report$note[idx],";shift-baseline-lower_reverse")
    }
  }
}

if (length(noisy_idx) > 0){
  for (experiment in report$expeiment[noisy_idx]){
    prior_code <- get_prior_code(meta,experiment)
    report <- process_noise(datadir,workdir,seriesName,prior_code,experiment,report,bins.lst,fit,genome)
  }
}

write.table(report, file=reportpath ,sep = '\t',quote = F,row.names = F)

# write log
log_path <- file.path(datadir,"processed",'logs')
dir.create(log_path, showWarnings = F, recursive = TRUE)
log_con <- file(file.path(log_path,paste0("process_",format(Sys.time(), "%y-%m-%d"),".log")),"a")  
log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"calibration complete","for",seriesName) 
writeLines(log,log_con)



