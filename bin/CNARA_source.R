###source code for S4 package "CNAQA" by AI Ni (ainijulia AT gmail DOT com), 2014

##### S4 class definition
setClass("CNProbe", representation = representation(sampleID="character", probeID="character", chr="integer", pos="numeric", log2="numeric", lengthProbe="integer"))
setClass("SpeakCNAno", representation = representation(sampleID="character", x="numeric", xchr="integer", lengthx="integer", iterations="integer", quality="numeric", fcf="matrix"))
setClass("Metrics", representation = representation(sampleID="character", speak="numeric", numberOfCNA="integer", cbsBreakpoints="integer", spread="numeric"))

##### S4 methods

# generics
setGeneric("calSpeakCNAno", function(object, ...) standardGeneric("calSpeakCNAno"))
setGeneric("calCBSBreakpoints", function(object, ...) standardGeneric("calCBSBreakpoints"))
setGeneric("calCBS", function(object, ...) standardGeneric("calCBS"))
setGeneric("calSpread", function(object, ...) standardGeneric("calSpread"))
setGeneric("plotFCF", function(object) standardGeneric("plotFCF"))
setGeneric("numberOfCNA", function(object) standardGeneric("numberOfCNA"))
setGeneric("speak", function(object) standardGeneric("speak"))
setGeneric("quality", function(object) standardGeneric("quality"))
setGeneric("assessQuality", function(object, ...) standardGeneric("assessQuality"))

# methods
setMethod("initialize", "CNProbe",
	function(.Object, ..., sampleID="", probeID="", chr="", pos="", log2="", lengthProbe="") {
		.Object@sampleID <- sampleID
		.Object@probeID <- probeID
		.Object@chr <- chr
		.Object@pos <- pos
		.Object@log2 <- log2

		deleteXY <- which(.Object@chr==23 | .Object@chr==24 | .Object@chr=="X" | .Object@chr=="Y")
		if (length(deleteXY) > 0) {
			.Object@sampleID <- .Object@sampleID[-deleteXY]
			.Object@probeID <- .Object@probeID[-deleteXY]
			.Object@chr <- .Object@chr[-deleteXY]
			.Object@pos <- .Object@pos[-deleteXY]
			.Object@log2 <- .Object@log2[-deleteXY]
		}

		.Object@log2[which(.Object@log2 == Inf | .Object@log2 == -Inf)] <- NA
		.Object@lengthProbe <- length(.Object@log2)

		callNextMethod(.Object, ...) #call parent class initialize()
	})


setMethod("calSpeakCNAno", "CNProbe", 
	function(object, ..., windowSize=100, nskip=11, iterations=120){
		library(robfilter)

		y.rr <- med.filter(object@log2,width=windowSize+1,online=FALSE)
		x <- y.rr$level$MED

		sampleRate <- round(object@lengthProbe/10000)

		if (sampleRate < 1) {
			sampleRate <- 1
		}

		sampleID <- object@sampleID
		x <- x[sampleRate*(1:floor(object@lengthProbe/sampleRate))]
		xChr <- object@chr[sampleRate*(1:floor(object@lengthProbe/sampleRate))]
		lengthProbe <- length(x)

		changingPointIndex <- c(1,lengthProbe+1)
		k <- NULL
		cfPointIndex <- c(1,lengthProbe+1)
		kcf <- NULL

		diff <- NULL

		#creat the fit-counter-fit matrix, y0, y1, cf0, cf1
		FCF <- matrix(mean(x),4,lengthProbe)

		quality <- NULL

		for (j in 1:iterations) {

			FCF[1,] <- FCF[2,]
			tmpM <- x-FCF[1,]
		    resultH0 <- t(tmpM)%*%(tmpM)
			
			FCF[3,] <- FCF[4,]
			tmpM <- x-FCF[3,]
		    resultCF0 <- t(tmpM)%*%(tmpM)

		    quality[j] <- resultCF0/resultH0

		    #if (j > 40 && quality[j-30] == max(quality[(j-30):j])) 
			#	break

			if (j == 1) {
				resultH1 <- NULL
				for (i in 1:lengthProbe) {
					skipPoint <- changingPointIndex
					for (indexSkip in 1:nskip) {
						skipPoint <- c(skipPoint, changingPointIndex-indexSkip, changingPointIndex+indexSkip)
					}
					if (any(i==skipPoint)) {
						resultH1[i] <- Inf
						next
					}
							
					cpTmp <- sort(c(changingPointIndex, i)) 
					
					for (t in 1:(length(cpTmp)-1))
						FCF[2,cpTmp[t]:(cpTmp[t+1]-1)] <- mean(x[cpTmp[t]:(cpTmp[t+1]-1)])
					
					tmpM <- x-FCF[2,]
					resultH1[i] <- t(tmpM)%*%(tmpM)
				}
				k <- which(resultH1 == min(resultH1))

			} else {
				k <- kcf[which(diff == max(diff))]
			}

			cpTmp <- sort(c(changingPointIndex, k)) 
				
			for (t in 1:(length(cpTmp)-1))
				FCF[2,cpTmp[t]:(cpTmp[t+1]-1)] <- mean(x[cpTmp[t]:(cpTmp[t+1]-1)])

			tmpM <- x-FCF[2,]
			resultH1 <- t(tmpM)%*%(tmpM)

			if (resultH1 > resultH0) {
				break
			} else {
			
				changingPointIndex <- c(changingPointIndex, k)
				
				############################
				# Calculate counter fit
				############################

				tmpCPIndex <- sort(changingPointIndex)
				tmpCPIndex <- c(tmpCPIndex[which(tmpCPIndex==k)-1], tmpCPIndex[which(tmpCPIndex==k)], tmpCPIndex[which(tmpCPIndex==k)+1])
				if (j > 1) {
					deleteIndex <- which(kcf==k)
					kcf <- kcf[-deleteIndex]
					diff <- diff[-deleteIndex]
				}

				for (t in 1:2) {
					resultCF1 <- NULL
					diffM <- NULL
					for (i in (tmpCPIndex[t]+1):(tmpCPIndex[t+1]-1)) {
						if (any(i==tmpCPIndex[t]+1:floor(nskip/2) | i==tmpCPIndex[t+1]-1:floor(nskip/2))) {
							resultCF1[i] <- Inf
							next
						}
						
						FCF[4,tmpCPIndex[t]:(i-1)] <- mean(x[tmpCPIndex[t]:(i-1)])
						FCF[4,i:(tmpCPIndex[t+1]-1)] <- mean(x[i:(tmpCPIndex[t+1]-1)])
						tmpM <- x[tmpCPIndex[t]:(tmpCPIndex[t+1]-1)]-FCF[4,tmpCPIndex[t]:(tmpCPIndex[t+1]-1)]
						resultCF1[i] <- t(tmpM)%*%(tmpM)
					}

					minIndexTmp <- which(resultCF1 == min(resultCF1, na.rm = TRUE))
					tmpkcf <- minIndexTmp[ceiling(length(minIndexTmp)/2)]
					kcf <- c(kcf, tmpkcf)

					FCF[4,tmpCPIndex[t]:(tmpkcf-1)] <- mean(x[tmpCPIndex[t]:(tmpkcf-1)])
					FCF[4,tmpkcf:(tmpCPIndex[t+1]-1)] <- mean(x[tmpkcf:(tmpCPIndex[t+1]-1)])

					tmpM <- x[tmpCPIndex[t]:(tmpCPIndex[t+1]-1)]-FCF[4,tmpCPIndex[t]:(tmpCPIndex[t+1]-1)]
					diffM <- x[tmpCPIndex[t]:(tmpCPIndex[t+1]-1)]-mean(x[tmpCPIndex[t]:(tmpCPIndex[t+1]-1)])
					diff <- c(diff, t(diffM)%*%(diffM)-t(tmpM)%*%(tmpM))

				}
				
				cfPointIndex <- c(1,lengthProbe+1,kcf)
				cfTmp <- sort(cfPointIndex) 
				
				for (t in 1:(length(cfTmp)-1))
					FCF[4,cfTmp[t]:(cfTmp[t+1]-1)] <- mean(x[cfTmp[t]:(cfTmp[t+1]-1)])	

			}
		}

		newSpeakCNAno <- new("SpeakCNAno", sampleID=sampleID, x=x, xchr=xChr, lengthx=lengthProbe, iterations=j, quality=quality, fcf=FCF)
		return(newSpeakCNAno)

	})


setMethod("calCBSBreakpoints", "CNProbe", 
	function(object, ...) {
		library(DNAcopy)

		CNProbe <- object
		sampledIndex <- (1:floor(object@lengthProbe/ceiling(object@lengthProbe/100000)))*ceiling(object@lengthProbe/100000)

		CNProbe@chr <- CNProbe@chr[sampledIndex]
		CNProbe@pos <- CNProbe@pos[sampledIndex]
		CNProbe@log2 <- CNProbe@log2[sampledIndex]

		CNA.object <- CNA(CNProbe@log2, CNProbe@chr, CNProbe@pos, data.type="logratio", sampleid=object@sampleID)
		                  
		smoothed.CNA.object <- smooth.CNA(CNA.object)
		segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1, alpha=0)		

		segN <- length(segment.smoothed.CNA.object$output$ID)
		return(segN)

	})

setMethod("calCBS", "CNProbe", 
	function(object, ..., verbose=1, alpha=0) {
		library(DNAcopy)

		CNA.object <- CNA(object@log2, object@chr, object@pos, data.type="logratio", sampleid=object@sampleID)
		                  
		smoothed.CNA.object <- smooth.CNA(CNA.object)
		segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=verbose, alpha=alpha)		

		return(segment.smoothed.CNA.object)

	})

setMethod("calSpread", "CNProbe", 
	function(object, ..., segFile="", verbose=1, alpha=0) {
		if (file.exists(segFile)) {
			tmpSeg <- read.table(segFile, header = FALSE, skip = 1, sep = "\t")
		} else {
			print("No segment file exists. Call CBS...")
			tmpSeg <- calCBS(object, verbose=verbose, alpha=alpha)
			tmpSeg <- tmpSeg$output
			tmpSeg <- tmpSeg[replace(seq(tmpSeg), 5:6, 6:5)]
		}
		names(tmpSeg) <- c("sample", "chr",	"start", "stop", "mean", "probes")

		tmpStd <- NULL
		segL <- NULL

		for (s in 1 : length(tmpSeg$sample)) {
			tmpStart <- which(object@pos == tmpSeg$start[s])
			indexStart <- tmpStart[which(object@chr[tmpStart] == tmpSeg$chr[s])]
			#if (length(indexStart)>1) {
			#	print("Start")
			#	print(indexStart)
			#	print(object@sampleID)
			#}

			tmpStop <- which(object@pos == tmpSeg$stop[s])
			indexStop <- tmpStop[which(object@chr[tmpStop] == tmpSeg$chr[s])]
			#if (length(indexStop)>1) {
			#	print("Stop")
			#	print(indexStop)
			#	print(object@sampleID)
			#}

			segL[s] <- tmpSeg$probes[s]
			if (length(indexStart) > 0 && length(indexStop) > 0) {
				if (indexStart[1]+segL[s]-1 == indexStop[1]) {
					iStart <- indexStart[1]
					iStop <- indexStop[1]
				} else if (indexStart[1]+segL[s]-1 == indexStop[length(indexStop)]){
					iStart <- indexStart[1]
					iStop <- indexStop[length(indexStop)]
				} else if (indexStart[length(indexStart)]+segL[s]-1 == indexStop[1]) {
					iStart <- indexStart[length(indexStart)]
					iStop <- indexStop[1]
				} else {
					iStart <- indexStart[length(indexStart)]
					iStop <- indexStop[length(indexStop)]
				}

				tmpStd[s] <- sd(object@log2[iStart:iStop],na.rm = TRUE)
			}
		}

		profileStd <- sum(tmpStd[which(!is.na(tmpStd))]*segL[which(!is.na(tmpStd))])/sum(segL[which(!is.na(tmpStd))])
		return(profileStd)

	})

setMethod("plotFCF", "SpeakCNAno", 
	function(object) {
		index <- c(1:object@lengthx)
		plot(index, object@x, type = "l", xlab = "chromosome", ylab = "log2 ratio", col = colors()[170], lwd = 1, xaxt="n", ylim=c(-0.8,0.6), main=object@sampleID)

		labelTextPos <- vector(length=22)
		for (i in 1:22)
			labelTextPos[i] <- median(which(object@xchr==i))
		axis(1, at=labelTextPos, labels=seq(1:22), tick=FALSE, cex.axis=0.8)

		tickPos <- vector(length=23)
		for (i in 1:22)
			tickPos[i] <- min(which(object@xchr==i))
		tickPos[23] <- max(which(object@xchr==22))

		abline(v=tickPos, col="grey")
		lines(index, object@fcf[3,], col = colors()[123], lwd = 2)
		lines(index, object@fcf[1,], col = colors()[134], lwd = 2)

	})

setMethod("numberOfCNA", "SpeakCNAno",
	function(object) which(object@quality == max(object@quality)))

setMethod("speak", "SpeakCNAno",
	function(object) max(object@quality))

setMethod("quality", "SpeakCNAno",
	function(object) object@quality)


setMethod("assessQuality", "Metrics",
	function(object, ..., svmClassifier = svmClassifier, speakThdGood = 1.5, speakThdPoor = 2.5) {

		testset <- cbind(log(object@speak), log(object@numberOfCNA), log(object@cbsBreakpoints), object@spread)
		testset <- data.frame(testset)
		names(testset) <- names(svmClassifier[[2]])

		fit <- svmClassifier[[1]]
		meanSet <- svmClassifier[[2]]
		sdSet <- svmClassifier[[3]]

		#normalize testing set
		testset$Speak <- (testset$Speak - meanSet[1])/sdSet[1]
		testset$Breakpoint_CNA  <- (testset$Breakpoint_CNA - meanSet[2])/sdSet[2]
		testset$Breakpoint_CBS <- (testset$Breakpoint_CBS - meanSet[3])/sdSet[3]
		testset$Spread <- (testset$Spread - meanSet[4])/sdSet[4]

		#normalize SpeakThdGood and SpeakThdPoor
		thdGood <- (log(speakThdGood)-meanSet[1])/sdSet[1]
		thdPoor <- (log(speakThdPoor)-meanSet[1])/sdSet[1]

		predictQuality <- predict(fit, newdata=testset, decision.values=T)
		flag <- NULL
		caseDiag <- NULL
		if (predictQuality==1) {
			flag <- "good"
			if (testset[[1]] > thdGood) {
				caseDiag <- "2: good quality, discernible CNAs with few waves"
			} else {
				caseDiag <- "5: good quality, control or without many CNAs"
			}
		} else {
			flag <- "poor"
			test1 <- testset
			test1[[3]] <- svmClassifier[[4]]
			test2 <- testset
			test2[[4]] <- svmClassifier[[5]]

			if (predict(fit, newdata=test1, decision.values=T)==1){
				caseDiag <- "3: poor quality, indiscernible CNAs with many waves"
			} else if (predict(fit, newdata=test2, decision.values=T)==1){
				caseDiag <- "4: poor quality, few CNAs, few waves, high probe value dispersion"
			} else {
				caseDiag <- "3: poor quality, indiscernible CNAs with many waves"
			}

			if (testset[[1]] > thdPoor) 
				caseDiag <- "1: hypersegmented, discernible CNAs with some waves"
		}

		assessment <- data.frame(predictQuality, attr(predictQuality, "decision.values"), flag, caseDiag)
		names(assessment) <- c("label", "decision.values", "flag", "caseDiag")
		return(assessment)
		
	})



##### functions
readProbe <- function(probeFile, sampleID) {
	if (file.exists(probeFile)) {
		tmpProbe <- read.table(probeFile, header = FALSE, skip = 1, sep = "\t")
		names(tmpProbe) <- c("ID", "chr", "position", "log2")
		tmpProbe <- tmpProbe[order(tmpProbe[,2], tmpProbe[,3]), ]
		newCNProbe <- new("CNProbe", sampleID=sampleID, probeID=tmpProbe$ID, chr=tmpProbe$chr, pos=tmpProbe$position, log2=tmpProbe$log2)
	
	} else {
		stop("Cannot find CNA probe file.")
	}

	return(newCNProbe)
}

trainSVM <- function(trainingFile) {
	library(kernlab)
	library(e1071)
	
	if (!file.exists(trainingFile)) 
		stop("Cannot find training file.")

	trainingSet <- read.table(trainingFile, header = TRUE, sep ="\t")

	#log normalize speak, breakpoint_CNA and breakpoint_CBS
	set <- trainingSet[4:8]
	set[1:3] <- log(set[1:3])

	#find the mean and standard deviation of the training set and then normalize it
	meanSet <- rep(0,4)
	meanSet <- colMeans(set[1:4])

	sdSet <- rep(0,4)
	sdSet <- sapply(set[1:4], sd)

	set[[1]] <- (set[[1]]-meanSet[1])/sdSet[1]
	set[[2]] <- (set[[2]]-meanSet[2])/sdSet[2]
	set[[3]] <- (set[[3]]-meanSet[3])/sdSet[3]
	set[[4]] <- (set[[4]]-meanSet[4])/sdSet[4]

	fit <- svm(Label ~., data=set, kernel="radial",gamma=0.4, type="C-classification", cross=10)
	#summary(fit)
	svmClassifier <- list(fit, meanSet, sdSet, min(set[[3]]), min(set[[4]]))

	return(svmClassifier)
}

createMetrics <- function(sampleID="", speak=0, numberOfCNA=0, cbsBreakpoints=0, spread=0) {
	newMetrics <- new("Metrics", sampleID=sampleID, speak=speak, numberOfCNA=numberOfCNA, cbsBreakpoints=cbsBreakpoints, spread=spread)

	return(newMetrics)
}