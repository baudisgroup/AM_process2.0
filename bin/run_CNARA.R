#! /usr/bin/env Rscript
args = R.utils::commandArgs(asValue = T)
datadir <- args$datadir
seriesName <- args$seriesName
workdir <- args$workdir
genome <- args$genome

source(file.path(workdir,'bin','CNARA_source.R'))

assessmentpath <- file.path(datadir,"processed",'data_quality_report','cnseg','CNARA_assessment')
dir.create(assessmentpath,showWarnings=F,recursive = T)

samplePath <- file.path(datadir,"processed")
probeFileName <- paste0("probes,cn,",genome,".tsv")
segFileName <- paste0("segments,cn,",genome,".tsv")

trainingFile <- file.path(workdir,"data","CNARA","trainingSet.txt")
classifier <- trainSVM(trainingFile)


arrayDir <- file.path(samplePath, seriesName)
arrayList <- list.files(arrayDir)
noArray <- length(arrayList)
assessment <- list()
for (arr in seq_len(noArray)){
    probeFile <- file.path(arrayDir, arrayList[arr], probeFileName)
    segFile <- file.path(arrayDir, arrayList[arr], segFileName)
    newCNProbe <- readProbe(probeFile=probeFile, sampleID=arrayList[arr])
    newSpeakCNAno <- calSpeakCNAno(newCNProbe)
    segNumberCBS <- invisible(calCBSBreakpoints(newCNProbe))
    segSpread <- calSpread(newCNProbe, segFile=segFile)
    
    CNAno <- numberOfCNA(newSpeakCNAno)
    Speak <- speak(newSpeakCNAno)
    
    CNProfileMetrics <- createMetrics(sampleID=arrayList[arr], speak=Speak, numberOfCNA=CNAno, cbsBreakpoints=segNumberCBS, spread=segSpread)
    
    assessment[[arr]] <- cbind(experiment = arrayList[arr],assessQuality(CNProfileMetrics, svmClassifier=classifier))
  }
assessment <- do.call(rbind,assessment)
write.table(assessment, file=file.path(assessmentpath,paste0(seriesName,'.txt')),sep = '\t',quote = F,row.names = F)

# write log
log_path <- file.path(datadir,"processed",'logs')
dir.create(log_path, showWarnings = F, recursive = TRUE)
log_con <- file(file.path(log_path,paste0("process_",format(Sys.time(), "%y-%m-%d"),".log")),"a")  

log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"CNARA complete","for",seriesName) 
writeLines(log,log_con)

