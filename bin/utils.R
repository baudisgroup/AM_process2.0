### return unprocessed files 
checkFileIncomplete <- function(force,path,filename,cids) {
  if (force){return(cids)}

  unprocessed_cids <- c()
  for(i in 1:length(cids)){
      checked_file <- file.path(path,cids[i],filename)

    if(!file.exists(checked_file)) {
      unprocessed_cids <- c(unprocessed_cids,cids[i])
    }
  }
  return(unprocessed_cids)
}

convertXYM <- function(aVector){
  aVector[aVector == 'X'] <- '23'
  aVector[aVector == 'Y'] <- '24'
  aVector[aVector == 'M'] <- '25'
  return(aVector)
}

# merge nearby signals by average
rm.near.dy <- function(x,distance=0.15) {
  count = 0
  i = 1
  
  while (i <= (length(x)+count-1)){
    c <- 0
    while (x[i+1-count]-x[i-c-count] <= distance) {
      i <- i+1
      c <- c+1
      if (i+1-count > length(x)) break
      #i+1-count > length(x)
    }
    x[i-c-count]<- mean(x[(i-c-count):(i-count)])
    if (c>0) x<- x[-c((i-c-count+1):(i-count))]
    count <- count+c
    i <- i+1
  }
  return(x)
}

adjustMedian <- function(datadir,seriesName,arrayName,adj_probe,genome){
  localProcessPath <- file.path(datadir,"processed",seriesName)
  currentsegpath <- file.path(localProcessPath,arrayName,'segments,cn.tsv')
  segfile <- read.table(currentsegpath,header = T,stringsAsFactors = F)
  newsegfile <- segfile
  # value column is 5; probe column is 6
  med <- round(weightedMedian(newsegfile[,5],newsegfile[,6]),5)
  newsegfile[,5] <- round(newsegfile[,5]-med,4)
  
  orisegpath <- file.path(localProcessPath,arrayName,'segments,cn,provenance.tsv')
  if (!file.exists(orisegpath)){
    write.table(segfile,file=orisegpath, sep="\t", quote=FALSE,row.names=FALSE)
  }
  write.table(newsegfile,file=currentsegpath, sep="\t", quote=FALSE,row.names=FALSE)
  
  if (adj_probe) {
    probefileName <- ifelse(genome == 'hg38','probes,cn,hg38.tsv','probes,cn.tsv')
    probefile <- read.table(file.path(localProcessPath,arrayName,probefileName),
                            header = T,stringsAsFactors = F)
    probefile[,4] <- as.numeric(probefile[,4])
    probefile <- probefile[!is.na(probefile[,4]),]
    probefile[,4] <- round(probefile[,4] -med,4)
    write.table(probefile,file.path(localProcessPath,arrayName,probefileName), 
                sep="\t", quote=FALSE,row.names=FALSE)
    log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"adjust segment and probe values with median:",med,"for",seriesName,arrayName)
  } else{
    log <- paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"),"adjust segment values with median:",med,"for",seriesName,arrayName)
  }
  
  return(log)
}


calculateLasso <- function(segmentData, lmd) {
  set.seed(1)
  x <- segmentData[,5]
  out = fusedlasso1d(x)
  minl <- min(out$lambda)
  lmd <- max(minl,lmd)
  beta1 = coef(out, lambda=lmd)$beta
  
  beta1 <- round(beta1,4)
  idx = which(abs(diff(beta1)) > 1e-4)
  return(list(beta1,idx,lmd))
}

calculateWeightedMean <- function(idx, segmentData) {
  lasti=1
  wm=rep(0,length(idx)+1)
  for(i in 1:(length(idx)+1)){
    wm[i] <- ifelse(i!=length(idx)+1,
                    weighted.mean(segmentData[,5][lasti:idx[i]],w=segmentData[,4][lasti:idx[i]]-segmentData[,3][lasti:idx[i]]+1),
                    weighted.mean(segmentData[,5][lasti:nrow(segmentData)],w=segmentData[,4][lasti:nrow(segmentData)]-segmentData[,3][lasti:nrow(segmentData)]+1))
    
    lasti <- idx[i]+1}
  wm <- round(wm, 4)
  return (wm)
}

stepFilter <- function(datadir,seriesName,arrayName,lmd){
  localProcessPath <- file.path(datadir,"processed",seriesName)
  segmentFile <- file.path(localProcessPath,arrayName,'segments,cn.tsv')
  seg <- read.table(segmentFile,stringsAsFactors = F,header = T,sep='\t')
  if (!is.numeric(seg[,2]))  seg[,2] <- as.numeric(convertXYM(seg[,2]))
  seg <- seg[order(seg[,2], seg[,3]),]
  #remove segments with probe < 10
  seg <- seg[seg[,6] >= 10,]
  if (dim(seg)[1] <= 1){return(seg)}
  
  #find point by abs(diff(beta)) > 1e-4 (the point differ from the next point)
  lassores <- calculateLasso(seg,lmd)
  idx <- lassores[[2]]
  ##chipType CHRO chr_start chr_end probes
  
  CHRO <- unique(seg[,2])
  chr_start <- sapply(CHRO, function(x) min(seg[which(seg[,2] == x),3]))
  chr_end <- sapply(CHRO, function(x) max(seg[which(seg[,2] == x),4]))
  probes <- sapply(CHRO, function(x) sum(seg[which(seg[,2] == x),6]))
  chrPos <- data.frame(CHRO = CHRO,
                       chr_start = chr_start,
                       chr_end = chr_end,
                       probes = probes)
  
  newseg <- data.frame()
  ### add all chromosome end index
  chro_end_idx <- NULL
  for (i in CHRO){
    chro_idx <- which(seg[,2] == i)
    chro_end_idx <- c(chro_end_idx, chro_idx[length(chro_idx)])
  }
  #change point, chr end, total end
  idx <- sort(union(c(idx,chro_end_idx),nrow(seg)))
  lmd <- lassores[[3]]
  #compute weighted means for interval between idx. The weight is seg len
  wm <- calculateWeightedMean(idx,seg)
  
  ##each idx is the end of segment and add the last row of all segments
  for (i in 1:length(idx)){
    if (i==1){
      countProbes <- 0 ## count probes that are in the hidden segments
      if(seg[idx[1],2]!=1){
        chrPassed <- 1 : seg[idx[1],2]
        for (j in 1:(length(chrPassed)-1)){
          ## j is the current chromosome
          if (!j %in% unique(chrPos[,1])) next
          tmpseg <- seg[idx[1],]
          tmpseg[,2] <- j
          tmpseg[,c(3,4)] <- chrPos[chrPos[,1]==j,c(2,3)]
          tmpseg[,5] <- wm[1]
          tmpseg[,6] <- chrPos[chrPos[,1]==j,4]
          newseg <- rbind(newseg, tmpseg)
          countProbes <- countProbes + tmpseg[,6]
        }
      }
      
      tmpseg <- seg[idx[1],]
      tmpseg[,3] <- chrPos[chrPos[,1]==tmpseg[,2],2]
      tmpseg[,5] <- wm[1]
      tmpseg[,6] <- sum(seg[c(1:idx[1]),6])-countProbes
      newseg<- rbind(newseg, tmpseg)
    }
    ## new segment and previous segment on the same chromosome
    else if (seg[idx[i],2] == seg[idx[i-1],2]){
      tmpseg <- vector()
      tmpseg <- seg[idx[i],]
      tmpseg[,3] <- seg[idx[i-1]+1,3] #start of last row +1's start
      tmpseg[,5] <- wm[i] # weighted mean
      tmpseg[,6] <- sum(seg[c((idx[i-1]+1):idx[i]),6]) # sum of probe numbers since last row +1
      newseg<- rbind(newseg, tmpseg)
    }
    ## new segment starts with a new chromosome
    else if (seg[idx[i],2] != seg[idx[i-1],2]){
      ## new segment starts with another chromosome instead of the same one
      countProbes <- 0 ## count probes that are in the hidden segments
      chrPassed <- seg[idx[i-1],2] : seg[idx[i],2]
      usedProbes <- sum(newseg[newseg[,2]==chrPassed[1],6]) ## probes used in the segments in that chromosome
      tmpseg <- seg[idx[i-1],]
      tmpseg[,3] <- tmpseg[,4]
      tmpseg[,4] <- chrPos[chrPos[,1]==chrPassed[1],3]
      tmpseg[,5] <- wm[i]
      tmpseg[,6] <- chrPos[chrPos[,1]==chrPassed[1],4] - usedProbes
      
      if (tmpseg[,6] > 0) newseg <- rbind(newseg, tmpseg)
      countProbes <- countProbes + tmpseg[,6]
      ## the chromosome is not the immediate next one
      if (length(chrPassed) > 2){
        for (j in 2:(length(chrPassed)-1)){
          ## 1. possibility: there's no record for the missing chromosome (no CNV)
          if (sum(chrPos[,1]==chrPassed[j]) == 0) next
          ## 2. possibility: CNV is at similar level as the previous and next chromosome (merging)
          tmpseg <- seg[idx[i-1],]
          tmpseg[,2] <- chrPassed[j]
          tmpseg[,c(3,4)] <- chrPos[chrPos[,1]==chrPassed[j],c(2,3)]
          tmpseg[,5] <- wm[i]
          tmpseg[,6] <- chrPos[chrPos[,1]==chrPassed[j],4]
          newseg <- rbind(newseg, tmpseg)
          countProbes <- countProbes + tmpseg[,6]
        }
      }
      tmpseg <- seg[idx[i],]
      tmpseg[,3] <- chrPos[chrPos[,1]==seg[idx[i],2],2]
      tmpseg[,5] <- wm[i]
      tmpseg[,6] <- sum(seg[(idx[i-1]+1):idx[i],6]) -countProbes
      newseg <- rbind(newseg, tmpseg)
    }
  }
  newseg[,5] <- round(newseg[,5],4)
  write.table(newseg, segmentFile, sep="\t", quote=FALSE,row.names=FALSE)
  return(newseg)
}


seg_cov <- function(seg){
  if (is.null(seg)){
    lowCov <- 0
    highCov <- 0
    normalCov <- 0
    lowCov_ratio <- 0
  } else{
    segNum <- dim(seg)[1]
    segLen <- seg[,4]-seg[,3]
    total.call.length <- sum(segLen)
    lowCov <- round(sum(segLen[seg$label %in% c('+1','-1')])/total.call.length,2)
    highCov <- round(sum(segLen[seg$label %in% c('+2','-2')])/total.call.length,2)
    normalCov <- round(sum(segLen[seg$label == '0'])/total.call.length,2)
    lowdupCov <- round(sum(segLen[seg$label %in% c('+1')])/total.call.length,2)
    lowdelCov <- round(sum(segLen[seg$label %in% c('-1')])/total.call.length,2)
    lowCov_ratio <- lowdupCov/lowdelCov 
    lowCov_ratio[is.na(lowCov_ratio)] <- 1
  }
  score <- data.frame('normal_cov'=normalCov,'lowCNA_cov'=lowCov,'lowCNA_ratio'=lowCov_ratio,'highCNA_cov'=highCov)
  return(score)
}

get_meta <- function(collection,series){
  if (!grepl('geo:',series)){
    series <- paste0("geo:",series)
  }
  # query data
  samples <- collection$find(paste0('{"external_references" : ', '{"$elemMatch" : ', '{"id" : ', '"',series,'"','}}}'))
  
  # check if samples are from geo database
  if ("analysis_info" %in% names(samples)){
    if ("id" %in% names(samples$histological_diagnosis)){
      colnames(samples$histological_diagnosis) <- c("histological_diagnosis_id","histological_diagnosis_label")
      samples_info <- cbind(samples$analysis_info,samples$histological_diagnosis)
      samples_info$experiment_id <- gsub('geo:','',samples_info$experiment_id)
      samples_info$series_id <- gsub('geo:','',samples_info$series_id)
      return(samples_info)
    }
  } 
  
  samples_info <- NULL
  return(samples_info)
  
}


generate_genomic_intervals <- function(workdir, genome='hg38', bin.size= 1000000, soft.expansion = 100000){
  cytoband.filepath <- file.path(workdir,'data','cytoband', paste0('UCSC_',genome,'_cytoband.txt'))
  cytoband.con <- file(cytoband.filepath,"r")
  cb <- read.csv(cytoband.con)
  close(cytoband.con)
  chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
  cb <- cb[cb$V1 %in% chroms, ]
  rownames(cb) <- seq(1:nrow(cb))
  arm.length <- numeric(48)
  
  for (i in seq(24)){
    temp <- cb[cb$V1 == chroms[i],]
    j <- 2*(i-1)
    arm.length[j+1] <- max(temp[temp$V4 %in% c('p11.1','p11','p11.11'),'V3'])
    arm.length[j+2] <- max(temp$V3)
  }
  
  bin.size = bin.size
  soft.expansion = soft.expansion
  count = 1
  result <- list()
  for ( i in seq(24)){
    p_max = arm.length[2*i-1]
    q_max = arm.length[2*i]
    arm = "p"
    start = 0
    p_first = p_max
    chr_result = list()
    chr_count = 1
    while (p_first >= bin.size + soft.expansion){
      p_first = p_first-bin.size
    }
    end = start + p_first
    while(start < q_max){
      interval.len = bin.size
      if (end > q_max){
        end = q_max
      }else if (q_max < end + soft.expansion){
        end = q_max
        interval.len = interval.len + soft.expansion
      }
      
      if (end >= p_max){
        arm = "q"
      }
      
      size = end-start
      chr_result[[chr_count]] <- data.frame(index= count, chromosome = gsub("chr","",chroms[i]), start=start, end=end, arm=arm, size= size )
      start = end
      end = end+interval.len
      count = count+ 1
      chr_count = chr_count+1
    }
    result[[chroms[i]]] <- do.call(rbind, chr_result)
  }
  return(result)
}


get_labelseg <- function(datadir,series,experiment = NULL){
  data <- list()
  series_path <- file.path(datadir,'processed',series)
  files <- list.files(series_path,recursive = F)
  if (!is.null(experiment)){
    files <- files[files %in% experiment]
  }
  for ( i in seq_len(length(files)) ){
    labelsegname <- file.path(series_path ,files[i],'labelsegments,cn.tsv')
    if(!file.exists(labelsegname)){ 
      next
      }
    data[[i]]  <- read.table(labelsegname, sep = '\t',header=T,colClasses = c("label"="character"))
  }
  data <- do.call(rbind,data)
  return(data)
}

fit_total <- function(seg,label){
  idx <- which(seg$label == label)
  if (length(idx) == 0){return()}
  sel.seg <- seg[idx,]
  sel.seg.len <- ceiling((sel.seg[,4]-sel.seg[,3]+1)/1000000)
  value <- rep(sel.seg[,6],sel.seg.len)
  fit <- MASS::fitdistr(value, "normal")
  return(fit)
}


compute_fit <- function(seg){
  low.dup.fit <- fit_total(seg,'+1')
  high.dup.fit <- fit_total(seg,'+2')
  low.del.fit <- fit_total(seg,'-1')
  high.del.fit <- fit_total(seg,'-2')
  normal.fit <- fit_total(seg,'0')
  
  fit.lst <- list()
  fit.lst[['+1']] <- low.dup.fit
  fit.lst[['+2']] <- high.dup.fit
  fit.lst[['-1']] <- low.del.fit
  fit.lst[['-2']] <- high.del.fit
  fit.lst[['0']] <- normal.fit
  return(fit.lst)
}

get_prior_code <- function(meta,experiment){
  if (!is.null(meta)){
    prior_code <- unique(meta$histological_diagnosis_id[meta$experiment_id %in% experiment])

    if (length(prior_code) > 1){ 
      stop(paste('metadata in Progenetix is invalid for',experiment))
    }
    }else{
    prior_code <- NULL
  }
  return(prior_code)
}


extract.bin.feature <- function(workdir,data,genome='hg38',bin.size= 1000000, soft.expansion = 100000,exclude.sex.chrom = NULL){
  bins <- generate_genomic_intervals(workdir = workdir, genome=genome, bin.size = bin.size,soft.expansion = soft.expansion)
  bins <- do.call(rbind, bins)
  if (!is.null(exclude.sex.chrom)){
    bins <- bins[!bins[,2] %in% exclude.sex.chrom,]
  }
  
  feature_dup_low <- list()
  feature_del_low <- list()
  feature_dup_high <- list()
  feature_del_high <- list()
  
  for (idx in c(1:length(unique(data[,1])))){
    ind.data <- data[data[,1] %in% unique(data[,1])[idx],]
    ind_dup_low <- rep(0,dim(bins)[1])
    ind_del_low <- rep(0,dim(bins)[1])
    ind_dup_high <- rep(0,dim(bins)[1])
    ind_del_high <- rep(0,dim(bins)[1])
    
    for (j in c(1:dim(ind.data)[1])){
      ind.seg.start <- ind.data[j,3]
      ind.seg.end <- ind.data[j,4]
      sel.bin <- which(bins$chromosome == ind.data[j,2] & bins[,4] > ind.seg.start & bins[,3] < ind.seg.end)
      if (length(sel.bin) == 0){next}
      ind_dup_high[sel.bin] <-   ind_dup_high[sel.bin] + as.numeric(ind.data[j,'label'] == '+2')
      ind_del_high[sel.bin] <- ind_del_high[sel.bin] + as.numeric(ind.data[j,'label'] == '-2')
      ind_dup_low[sel.bin] <- ind_dup_low[sel.bin] + as.numeric(ind.data[j,'label'] == '+1')
      ind_del_low[sel.bin] <-  ind_del_low[sel.bin] + as.numeric(ind.data[j,'label'] == '-1')
    }
    ind_dup_high[ind_dup_high > 1] <- 1
    ind_del_high[ind_del_high > 1] <- 1
    ind_dup_low[ind_dup_low > 1] <- 1
    ind_del_low[ind_del_low > 1] <- 1
    feature_dup_high[[idx]] <- ind_dup_high
    feature_del_high[[idx]] <- ind_del_high
    feature_dup_low[[idx]] <- ind_dup_low
    feature_del_low[[idx]] <- ind_del_low
  }
  feature_dup_high <- do.call(rbind,feature_dup_high)
  feature_del_high <- do.call(rbind,feature_del_high)
  feature_dup_low <- do.call(rbind,feature_dup_low)
  feature_del_low <- do.call(rbind,feature_del_low)
  
  if (length(unique(data[,1])) > 1){
    rownames(feature_dup_low) <- unique(data[,1])
    rownames(feature_dup_high) <- unique(data[,1])
    rownames(feature_del_low) <- unique(data[,1])
    rownames(feature_del_high) <- unique(data[,1])
  }
  
  
  feature.list <- list()
  feature.list[['low.dup']] <- feature_dup_low
  feature.list[['high.dup']] <- feature_dup_high
  feature.list[['low.del']] <- feature_del_low
  feature.list[['high.del']] <- feature_del_high
  feature.list[['samples']] <- unique(data[,1])
  return(feature.list)
}

compute_freq <- function(workdir,data,output,genome='hg38'){
  features <- extract.bin.feature(workdir,data,genome = genome,exclude.sex.chrom = c('X','Y'))
  
  freq_low_del <-  features$low.del
  freq_low_del <- colSums(freq_low_del)/nrow(freq_low_del)
  
  freq_low_dup <-  features$low.dup
  freq_low_dup <- colSums(freq_low_dup)/nrow(freq_low_dup)
  
  freq_high_del <-  features$high.del
  freq_high_del <- colSums(freq_high_del)/nrow(freq_high_del)
  
  freq_high_dup <- features$high.dup
  freq_high_dup <- colSums(freq_high_dup)/nrow(freq_high_dup)
  
  freq.lst <- list()
  freq.lst[['+1']] <- freq_low_dup
  freq.lst[['+2']] <- freq_high_dup
  freq.lst[['-1']] <- freq_low_del
  freq.lst[['-2']] <- freq_high_del
  freq.lst[['0']] <- 1-colSums(rbind(freq_low_dup,freq_low_del, freq_high_dup,freq_high_del))
    
  saveRDS(freq.lst,output)
}

compute_pointprob <- function(fit, point){
  prob <- pnorm(point, mean=fit$estimate[1], sd=fit$estimate[2], lower.tail = point <= fit$estimate[1]) 
  return(prob)
}

compute_posterior <- function(datadir,workdir,series,prior_code,seg,bins,fit,genome='hg38'){
  if (length(prior_code) > 0 & genome == 'hg38'){
    reference_freq <- readRDS(file.path(workdir,'data','calibration-rsc','hg38_precomputed_complete_freq_from_progenetix.rds'))
    reference_freq <- reference_freq$frequency[[prior_code]]
    prior <- list()
    prior[['+1']] <- reference_freq[c(1:2892)]
    prior[['-1']] <- reference_freq[c(2893:5784)]
    prior[['0']] <- 1- prior[['+1']] - prior[['-1']]
  } else{
    dir.create(file.path(datadir,'processed','seriesfreq'),showWarnings=F)
    freqpath <- file.path(datadir,'processed','seriesfreq',paste0(series,'-freq.rds'))
    if (!file.exists(freqpath)){
      data <- get_labelseg(datadir,series)
        # when a series contains 2 samples, average frequency doesn't make sense
      if (length(unique(data[,1])) == 2 ){
        data <- data[data[,1] != seg[1,1],]
      }
      compute_freq(workdir,data,freqpath,genome=genome)
    }
    prior <- readRDS(freqpath)
  } 
  
  # reference code only has normal, dup and del states
  if (length(prior_code) > 0){
    sel.seg <- seg[seg$label %in% c('0','+1','-1'),]
  } else{
    sel.seg <- seg
  }
  
  posterior <- vapply(seq_len(dim(sel.seg)[1]), function(i){
    x <- sel.seg[i,]
    ind.fit <- fit[[x$label]]
    lld <- compute_pointprob(ind.fit, x[,6]) * ceiling((x[,4]-x[,3])/1000000)
    ind.prior <- prior[[x$label]]
    ind.bins <- bins[[as.numeric(x[,2])]]
    ind.intersect.bins <- ind.bins[ind.bins[,3] < x[,4] & ind.bins[,4] > x[,3],'index']
    ind.avg.prior <- mean(ind.prior[ind.intersect.bins])
    return(ind.avg.prior * lld)
  },numeric(1))
  
    sel.seg$posterior <- posterior
    return(sel.seg)
}

compare_whole_shift <- function(datadir,workdir,series,meta,oriseriesSeg,shift_idx_type,opposite_shift_idx_type,baseshift,report,ori_fit,bins.lst,genome){
  seriesSeg <- oriseriesSeg
  seriesSeg[seriesSeg[,1] %in% report$expeiment[shift_idx_type],] <- labelseg(seriesSeg[seriesSeg[,1] %in% report$expeiment[shift_idx_type],],baseshift = baseshift,genome=genome)
  shift_fit <- compute_fit(seriesSeg)
  shiftSeg <- list()
  for (experiment in report$expeiment[shift_idx_type]){
    prior_code <- get_prior_code(meta,experiment)
    oriposteior <- compute_posterior(datadir,workdir,series,prior_code,oriseriesSeg[oriseriesSeg[,1] == experiment,],bins.lst,ori_fit,genome=genome)
    relabelposteior <- compute_posterior(datadir,workdir,series,prior_code,seriesSeg[seriesSeg[,1] == experiment,],bins.lst,shift_fit,genome=genome)
    shiftSeg[[experiment]] <- merge(oriposteior,relabelposteior,by=c(1,2,3,4),all=F)
  }
  
  for (experiment in report$expeiment[opposite_shift_idx_type]){
    prior_code <- get_prior_code(meta,experiment)
    if (length(prior_code) > 1){next}
    oriposteior <- compute_posterior(datadir,workdir,series,prior_code,oriseriesSeg[oriseriesSeg[,1] == experiment,],bins.lst,ori_fit,genome=genome)
    relabelposteior <- compute_posterior(datadir,workdir,series,prior_code,oriseriesSeg[oriseriesSeg[,1] == experiment,],bins.lst,shift_fit,genome=genome)
    shiftSeg[[experiment]] <- merge(oriposteior,relabelposteior,by=c(1,2,3,4),all=F)
  }
  
  shiftSeg <- do.call(rbind,shiftSeg)
  if (sum(shiftSeg$posterior.x) >= sum(shiftSeg$posterior.y)){
    shift = "n"
    posterior <- sum(shiftSeg$posterior.x)
  } else{
    shift = baseshift
    posterior <- sum(shiftSeg$posterior.y)
  }
  
  compare_list <- list()
  compare_list[['shift']] <- shift
  compare_list[['posterior']] <- posterior
  return(compare_list)
}

update_seg <- function(series_path,experiment,relabelSeg,report,score,problem,genome='hg38',mergeseg=F){
  #update local data
  if (is.null(relabelSeg)){
    relabelSeg <- read.table(file.path(series_path,experiment,'segments,cn.tsv'),sep = '\t',header = T)
    relabelSeg <- relabelSeg[,c(1,2,3,4,6,5)]
    plot_segment_label(filepath = file.path(series_path, experiment),filename='label_cnsegment.pdf',data = relabelSeg,assembly = genome,no_label = T)
    problem <- paste0(problem,';failed-to-label')
  } else{
    write.table(relabelSeg, file=file.path(series_path,experiment,'labelsegments,cn.tsv'),sep = '\t',quote=F,row.names = F)
    plot_segment_label(filepath = file.path(series_path,experiment),filename='label_cnsegment.pdf',data = relabelSeg,assembly = genome)
  }
    
  #update report
  idx <- which(report$expeiment == experiment)
  report[idx,c('normal_cov','lowCNA_cov','lowCNA_ratio','highCNA_cov')] <- c(score$normal_cov,score$lowCNA_cov,score$lowCNA_ratio,score$highCNA_cov)
  report$note[idx] <- paste0(report$note[idx],";", problem)
    
  #update other features for merged seg
  if (mergeseg){
    lrrsegsd <- round(sd(rep(relabelSeg[,6],abs(relabelSeg[,5]))),3)
    report[idx,c('segment_num',	'LLR_segsd')] <- c(dim(relabelSeg)[1],lrrsegsd)
  }
  
  return(report)
}

shift_baseline <- function(datadir,workdir,series,prior_code,experiment,report,bins.lst,fit,whole_shift,genome,shift){
  series_path <- file.path(datadir,'processed',series)

  labelSeg <- read.table(file.path(series_path ,experiment,'labelsegments,cn.tsv'),sep = '\t',header=T,colClasses = c("label"="character"))
  relabelSeg <- labelseg(data = labelSeg,baseshift = substr(shift,1,1),genome=genome)
  
  idx <- which(report$expeiment == experiment)
  if (all(relabelSeg$label == labelSeg$label)){
    report$note[idx] <- paste0(report$note[idx],";shift-baseline-",shift,"_failed")
    return(report)
  }
  
  if (whole_shift){
    cov <- seg_cov(relabelSeg)
    report <- update_seg(series_path=series_path,experiment=experiment,
                         relabelSeg=relabelSeg,report=report,score=cov,problem=paste0("shift-baseline-",shift),
                         genome=genome)
  } else{
    oriposteior <- compute_posterior(datadir,workdir,series,prior_code,labelSeg,bins.lst,fit,genome=genome)
    relabelposteior <- compute_posterior(datadir,workdir,series,prior_code,relabelSeg,bins.lst,fit,genome=genome)
    merge.seg <- merge(oriposteior,relabelposteior,by=c(1,2,3,4),all=F)
    if (sum(merge.seg$posterior.y) > sum(merge.seg$posterior.x)){
      cov <- seg_cov(relabelSeg)
      report <- update_seg(series_path=series_path,experiment=experiment,
                           relabelSeg=relabelSeg,report=report,score=cov,problem=paste0("shift-baseline-",shift),
                           genome=genome)
      
    } else{
      report$note[idx] <- paste0(report$note[idx],";shift-baseline-",shift,"_reverse")
    }
  }
  
  return(report)
}

cutoff_calling <- function(seg,low_thre,high_thre){
  lrr <- seg[,6]
  recall <- rep('0',length(lrr))
  recall[lrr >= low_thre] <- '+1'
  recall[lrr <= -low_thre] <- '-1'
  recall[lrr >= high_thre] <- '+2'
  recall[lrr <= -high_thre] <- '-2'
  seg$label <- recall
  return(seg)
}

process_noise <- function(datadir,workdir,series,prior_code,experiment,report,bins.lst,fit,genome){
  
  series_path <- file.path(datadir,'processed',series)

  labelSeg <- read.table(file.path(series_path ,experiment,'labelsegments,cn.tsv'),sep = '\t',header=T,colClasses = c("label"="character"))
  assess <- compute_posterior(datadir,workdir,series,prior_code,labelSeg,bins.lst,fit,genome=genome)
  target <- labelSeg
  method <- 'nochange'
  
  for ( l_thre in c(0.1,0.15,0.3)){
    for (h_thre in c(1,1.5,2)){
      relabelSeg <- cutoff_calling(labelSeg,l_thre, h_thre)
      reassess <-  compute_posterior(datadir,workdir,series,prior_code,relabelSeg,bins.lst,fit,genome=genome)
      merge.seg <- merge(assess,reassess,by=c(1,2,3,4),all=F)
      if (sum(merge.seg$posterior.y) > sum(merge.seg$posterior.x)){
        assess <- reassess
        target <- relabelSeg
        method <- paste0('lowthre_',l_thre,'_highthre_',h_thre)
      }
    }
  }
  if (method != 'nochange'){
    cov <- seg_cov(target)
    report <- update_seg(series_path=series_path,experiment=experiment,relabelSeg=target,report=report,score=cov,problem=paste0("noisy_",method),genome=genome)
  } else{
    idx <- which(report$expeiment == experiment)
    report$note[idx] <- paste0(report$note[idx],";noisy_",method)
  }
  
  return(report)
}





