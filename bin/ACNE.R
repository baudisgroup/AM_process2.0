#! /usr/bin/env Rscript
suppressWarnings(suppressMessages(library(aroma.affymetrix)))
suppressWarnings(suppressMessages(library(ACNE)))

args = R.utils::commandArgs(asValue = T)
datadir <- args$datadir
seriesName <- args$seriesName
cleanup <- as.logical(args$cleanup)
workdir <- args$workdir
memory <- as.numeric(args$memory)
force <- as.logical(args$force)

source(file.path(workdir,'bin','utils.R'))

dir.create(file.path(datadir,"processed"),showWarnings = F)
dir.create(file.path(datadir,"processed",'logs'),showWarnings=F)
localProcessPath <- file.path(datadir,"processed",seriesName)
dir.create(localProcessPath,showWarnings = F)

chipTypes <- list.files(file.path(datadir,"rawData",seriesName))

####### repeat for each chipType ######
for (chipType in chipTypes){
  if (!chipType %in% list.files(file.path(datadir,"annotationData","chipTypes"))){
    stop("chipType ", chipType," is not available in ",seriesName)
  }
  rawDatapath <- file.path(datadir,"rawData",seriesName, chipType)
  files <- list.files(rawDatapath)
  cids <- gsub(".CEL","",files)
  if (length(cids) == 0) stop("rawData is not found for ", seriesName)
  
  
  targetcids <- checkFileIncomplete(force, localProcessPath, 'probes,fracb.tsv', cids)
  
  ################## run ACNE########################
  if (length(targetcids) > 0){
    setOption(aromaSettings, "memory/ram", memory)  
    verbose <- Arguments$getVerbose(F, timestamp=TRUE)
    setwd(datadir)
    
    if (chipType %in% c("GenomeWideSNP_6","GenomeWideSNP_5")){
      cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full")
    }else {
      cdf <- AffymetrixCdfFile$byChipType(chipType)}
    
    gi <- getGenomeInformation(cdf)
    si <- getSnpInformation(cdf)
    acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE))
    cs <- AffymetrixCelSet$byName(seriesName, cdf=cdf)
    acc <- AllelicCrosstalkCalibration(cs, model="CRMAv2")
    csC <- process(acc, verbose=verbose)
    bpn <- BasePositionNormalization(csC, target="zero")
    csN <- process(bpn, verbose=verbose)
    plm <- NmfSnpPlm(csN, mergeStrands=TRUE)
    
    if (length(findUnitsTodo(plm)) > 0) {
      # Fit CN probes quickly (~5-10s/array + some overhead)
      units <- fitCnProbes(plm, verbose=verbose)
      
      # Fit remaining units, i.e. SNPs (~5-10min/array)
      units <- fit(plm, verbose=verbose)
    }
    
    ces <- getChipEffectSet(plm)
    
    options("scipen" = 9, digits=4)
    for (chromosome in 1:23) {
      units <- getUnitsOnChromosome(gi, chromosome=chromosome)
      pos <- getPositions(gi, units=units)
      unitnames <- getUnitNames(cdf,units=units)
      for (ii in 1:length(cs$Names)){
        cf <- ces[[ii]]
        data <- extractTotalAndFreqB(cf, units=units)
        beta <- data[,"freqB"]
        fracB <- RawAlleleBFractions(beta, pos, chromosome=chromosome)
        ID <- unitnames[which(!is.na(beta))]
        fracB <-extractSubset(fracB,which(!is.na(beta))) # to erase the CN probes
        fracB <-cbind(ID,fracB)
        colnames(fracB)[2:4] <- c("CHRO","BASEPOS","VALUE")
        fracB$VALUE <- round(fracB$VALUE,4)
        samplepath <- file.path(localProcessPath,cs$Names[ii])
        if (dir.exists(samplepath) == 0) dir.create(samplepath)
        write.table(fracB,sprintf("%s/probes,fracb,chr%d.tsv",samplepath,chromosome),quote=F,sep="\t",row.names = F)
      }
    }
    
    if (cleanup==T) {
      tmpfiledir <- paste0(datadir,c("/plmData/","/probeData/"))
      for (dir in tmpfiledir) {
        tmpfiles <- list.files(dir)
        tmpfiles <- tmpfiles[grep(seriesName,tmpfiles)]
        for (tmpfile in tmpfiles){
          unlink(file.path(tmpfiledir,tmpfile),recursive = TRUE)
        }
      }
    }
    
    system(sprintf('for i in %s/*/;
                    do cp $i/probes,fracb,chr1.tsv $i/probes,fracb.tsv;
                  for j in {2..23};
                    do tail -n +2 $i/probes,fracb,chr$j.tsv >> $i/probes,fracb.tsv;
                  done;
                  rm $i/probes,fracb,chr*;
                  done',localProcessPath))
    
  }
}

probefiles <- list.files(localProcessPath,recursive = T, full.names =T )  

if (sum(grepl("probes,fracb.tsv",probefiles)) == length(cids)){
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"ACNE",seriesName,"complete") 
} else{
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"ACNE",seriesName,"complete but file numbers aberrant")
}
log_con <- file(file.path(datadir,"processed",'logs',paste0("process_",format(Sys.time(), "%y-%m-%d"),".log")),"a")  
writeLines(log,log_con)


  