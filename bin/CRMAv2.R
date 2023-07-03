#! /usr/bin/env Rscript
suppressWarnings(suppressMessages(library(aroma.affymetrix)))

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
if (length(chipTypes) == 0) stop("rawData is not found for ", seriesName)
####### repeat for each chipType ######
cidcount = 0
for (chipType in chipTypes){
  if (!chipType %in% list.files(file.path(datadir,"annotationData","chipTypes"))){
    stop("chipType ", chipType," is not available in ",seriesName)
  } 
  
  rawDatapath <- file.path(datadir,"rawData",seriesName, chipType)
  files <- list.files(rawDatapath)
  
  cids <- gsub(".CEL","",files)
  if (length(cids) == 0) stop("rawData is not found for ", seriesName)
  cidcount = cidcount + length(cids)
  targetcids <- checkFileIncomplete(force, localProcessPath, 'probes,cn,hg19.tsv', cids)
  
  ################## run CRMAv2########################
  if (length(targetcids) > 0){
    
    setOption(aromaSettings, "memory/ram", memory)  
    verbose <- Arguments$getVerbose(F, timestamp=TRUE)
    setwd(datadir)

    if (chipType %in% c("GenomeWideSNP_6","GenomeWideSNP_5")){
      cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full")
    } else {
      cdf <- AffymetrixCdfFile$byChipType(chipType)
    }
    
    gi <- getGenomeInformation(cdf)
    si <- getSnpInformation(cdf)
    acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE))
    csR <- AffymetrixCelSet$byName(seriesName, cdf=cdf)
    acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2")
    csC <- process(acc, verbose=verbose)
    bpn <- BasePositionNormalization(csC, target="zero")
    csN <- process(bpn, verbose=verbose)
    plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE)
    if (length(findUnitsTodo(plm)) > 0) {
      # Fit CN probes quickly (~5-10s/array + some overhead)
      units <- fitCnProbes(plm, verbose=verbose)
      # Fit remaining units, i.e. SNPs (~5-10min/array)
      units <- fit(plm, verbose=verbose)
    }
    
    ces <- getChipEffectSet(plm)
    fln <- FragmentLengthNormalization(ces, target="zero")
    cesN <- process(fln, verbose=verbose) #26h for GSE19949
  
    # load external reference 
    load(file.path(datadir,"ReferenceFile",chipType,"ceR.RData"))
    ### change saveLoadReference path 
    saveLoadReferenceOrifilepath <- getPath(saveLoadReference)
    saveLoadReferenceNewfilepath <-gsub(".*plmData/",paste0(datadir,'/plmData/'),saveLoadReferenceOrifilepath)
    saveLoadReferenceNewfilepath <- Arguments$getReadablePathname(paste0(getFullName(saveLoadReference),'.CEL'), path=saveLoadReferenceNewfilepath, absolute=TRUE, mustExist=T)
    saveLoadReference$.pathname <- saveLoadReferenceNewfilepath
    
    saveLoadReferenceOricdfpath <- getPath(saveLoadReference$.cdf)
    saveLoadReferenceNewcdfpath <-gsub(".*annotationData/",paste0(datadir,'/annotationData/'),saveLoadReferenceOricdfpath)
    saveLoadReferenceNewcdfpath <- Arguments$getReadablePathname(paste0(getFullName(saveLoadReference$.cdf),'.CDF'), path=saveLoadReferenceNewcdfpath, absolute=TRUE, mustExist=T)
    saveLoadReference$.cdf$.pathname <- saveLoadReferenceNewcdfpath
    
    # Write out logR in a file
    for (chr in 1:23) {
      # Define the chromosome
      cdf <- getCdf(cesN)
      gi <- getGenomeInformation(cdf)
      units <- getUnitsOnChromosome(gi, chromosome=chr)
      pos <- getPositions(gi, units=units)
      pos <- as.numeric(pos)
      unitnames <- getUnitNames(cdf,units=units)
      # Retrieving CN estimates of the reference in the above region
      thetaR <- extractTheta(saveLoadReference,units=units)
      
      # Retrieving the corresponding estimates of samples
      sampleIDs <- gsub(".CEL","",list.files(paste(datadir,'rawData',seriesName,chipType,sep="/")))
      for (i in 1:length(sampleIDs)) {
        ce <- cesN[[indexOf(cesN, sampleIDs[i])]]
        theta <- extractTheta(ce, units=units)
        # Calculate the relative CNs
        R <- theta/thetaR
        logR <- round(log2(R),4)
        out <- data.frame(cbind(unitnames,rep(chr,length(pos)),pos,logR),stringsAsFactors = F)
        out <- out[order(as.numeric(out$pos)),]
        colnames(out) <- c("PROBEID","CHRO","BASEPOS","VALUE")
        samplepath <- file.path(localProcessPath,sampleIDs[i])
        if (!dir.exists(samplepath))  dir.create(samplepath)
        samplefile <- file.path(samplepath,'probes,cn,hg19.tsv')
        if (!file.exists(samplefile)){
          write.table(out,file = samplefile, sep = "\t",row.names = F,col.names = T,quote = F)
        } else{
          write.table(out,file = samplefile, sep = "\t",row.names = F,col.names = F,quote = F,append=T)
        }
      }
    }
    # If true, delete intermediate files under plmData, probeData
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
    
}
}

probefiles <- list.files(localProcessPath,recursive = T, full.names =T )  

if (sum(grepl("probes,cn,hg19.tsv",probefiles)) == cidcount){
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"CRMAv2",seriesName,"complete") 
} else{
  log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"CRMAv2",seriesName,"complete but file numbers aberrant") 
}

log_con <- file(file.path(datadir,"processed",'logs',paste0("process_",format(Sys.time(), "%y-%m-%d"),".log")),"a")  
writeLines(log,log_con)





