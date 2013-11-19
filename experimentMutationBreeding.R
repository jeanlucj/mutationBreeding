##########################################################################################
# The experiment file has two kinds of objects:
# 1. Objects that all replications and settings will receive. These define the environment
#	of the experiment
# 2. Objects that will create the founder files for a particular replication. These define
#	the variation from one replication of the setting to the next.

##########################################################################################
# 1. Universal environment objects
# Minimally, it includes:
# experimentName -- a character variable that is a reference
# breedFuncList -- a character vector with the names of all the breeding simulation functions needed
# Note: the experiment file and the selection program file are sourced one after the
#	other so the selection program file can provide or supercede anything in the experiment file
#.libPaths("~/Rlibs")
#print(.libPaths())
library(rrBLUP)
# library(Rglpk)
# install.packages("regress", lib="/Library/Frameworks/R.framework/Versions/2.15/Resources/library", repos="http://lib.stat.cmu.edu/R/CRAN")
# library(regress)
#library(glmnet)
# Breeding simulation functions are all the functions needed to simulation gamete formation etc.
# The idea is to include only the functions needed, but for now it's easier to include them all...
source("~/Documents/Association/Functions/breedingSimulationFunctionsDocumented.R")

experimentName <- "GS_Update"

nReplications <- 360

####################################################################################
# Extra functions for all replications and settings
####################################################################################
###############################################################
oneLifecycle <- function(cycle, phenotype, select, crossAdvance, analyze){
	cat(cycle, "######################## ", Sys.getpid(), "\n")
	breedingData <<- phenotype(cycle, breedingData)
	breedingData <<- select(cycle, breedingData)
	breedingData <<- analyze(cycle, breedingData)
	breedingData <<- crossAdvance(cycle, breedingData)
	return(NULL)
}#END oneLifecycle

initializeReduceRecCentro <- function(speciesData){
# Reduce recombination in the middle of each chromosome.
# Specific values used are very ad hoc
# I'm going to
# 1. Project the first and last cut fraction of the physical map on to
#    the first and last 45% of the recombination map
# 2. Project the middle (1 - 2*cut) fraction of the physical map on to the
#    middle 10% of the recombination map using this sigmoidal function
# cut <- 0.25 reduces 50% of the physical map to 10% of the recombination map
	cut <- 0.25 # Has to be less than 0.5; smaller number compresses physical genome more
	theChr <- sort(unique(speciesData$map[,"Chr"]))
	for (chr in theChr){
		thePos <- speciesData$map[speciesData$map[,"Chr"] == chr, "Pos"]
		minPos <- min(thePos)
		maxPos <- max(thePos)
		newPos <- (thePos - minPos) / (maxPos - minPos)
		belCut <- newPos <= cut
		aboCut <- newPos > (1 - cut)
		newPos[belCut] <- newPos[belCut] * 0.45 / cut
		newPos[aboCut] <- 0.55 + (newPos[aboCut] - (1 - cut)) * 0.45 / cut
		midPos <- newPos[!(belCut | aboCut)]
		midPos <- 0.01 + (midPos - cut) * 0.98 / (1 - 2 * cut)
		invLogistic <- function(y) return(log(y / (1 - y)))
		midPos <- invLogistic(midPos) / invLogistic(0.99)
		midPos <- 0.45 + (midPos + 1) * 0.1 / 2
		newPos[!(belCut | aboCut)] <- midPos
		newPos <- minPos + newPos * (maxPos - minPos)		
		speciesData$map[speciesData$map[,"Chr"] == chr, "Pos"] <- newPos
	}
	# fix cumuPos; chrMax should still be fine
	speciesData$map <- makeChrMax(speciesData$map)$map
	
	return(speciesData)
}

####################################################################################
# Experiment-wide replacement functions
####################################################################################
# !!! No mutation version: if you want mutation, do it outside of makeGamete 
makeGamete <- function(genoVec, map, chrMax){
	nChr <- length(chrMax)
	genomeLength <- chrMax[nChr]
	nLoc <- nrow(map)
	# number of recombinations
	nCrossOver <- rpois(1, genomeLength/100)
	# positions of rec: two parts: cross overs and recombinations between chromosomes
	posRec <- sort(c(runif(nCrossOver, 0, genomeLength), chrMax[sample(nChr-1, rbinom(1, nChr-1, 0.5))]))
	nRec <- length(posRec)
	posRec <- c(-1, posRec, genomeLength)
	gamete <- numeric(nLoc)
	parentalStrand <- round(runif(1))
	for (rec in 1:(nRec+1)){
		lociThisSegment <- map[,"cumuPos"]>posRec[rec] & map[,"cumuPos"]<=posRec[rec+1]
		gamete[lociThisSegment] <- parentalStrand
		parentalStrand <- 1 - parentalStrand
	}
	newGamete <- genoVec[1:nLoc * 2 - gamete]
	return(newGamete)
}

# Take two parents, make a gamete from each one, and combine into one progeny
# markers has genotype vectors in rows and markers in columns
makeProgeny <- function(parents, markers, map, chrMax){
	nLoc <- nrow(map)
	progeny <- numeric(2*nLoc)
	progeny[1:nLoc*2-1] <- makeGamete(markers[parents[1],], map, chrMax)
	progeny[1:nLoc*2] <- makeGamete(markers[parents[2],], map, chrMax)
	return(progeny)
}

# Take mutation function out of makeGamete. First mate without mutation.
# Second add in mutation.  I did not find a flaw here (13 June 2013).
randomMateDiploid <- function(markers, map, popSize, chrMax, selfingAllowed=FALSE){
	nPar <- nrow(markers)
	parMat <- sapply(rep(nPar, popSize), sample, size=2, replace=selfingAllowed)
	genoMat <- t(apply(parMat, 2, makeProgeny, markers, map, chrMax))
	rownames(genoMat) <- paste("Progeny", 1:popSize, sep="")
	# Now add in newly mutated loci
	if (!is.null(speciesData$mutParms)){
		nChr <- 7 # WARNING 7 chr hardwired
		lengthChr <- 150 # WARNING 150 cM hardwired
		mutPerGam <- rpois(2*popSize, speciesData$mutParms$mutNum)
		nLoc <- nrow(map)
		nMut <- sum(mutPerGam)
		cumuPos <- runif(nMut, 0, nChr * lengthChr)
		Chr <- cumuPos %/% lengthChr + 1
		Pos <- cumuPos %% lengthChr
		Pos <- speciesData$unifToPos(Pos)
		cumuPos <- (Chr - 1) * lengthChr + Pos # Recalculate this after unifToPos
		# Make loci with one mutant individual
		nGam <- 2*popSize
		gamNumVec <- unlist(sapply(1:nGam, function(gamNum) rep(gamNum, mutPerGam[gamNum]))) - 1
		indNumVec <- gamNumVec %/% 2 + 1
		matOrPat <- gamNumVec %% 2
		coord <- cbind(matOrPat * popSize + indNumVec, 1:nMut)
		# ancAllele <- rbinom(nMut, 1, 0.5) * 2 - 1
		# newGeno <- matrix(rep(ancAllele, each=nGam), nrow=nGam)
		newGeno <- matrix(rep(-1, each=nGam), nrow=nGam)
		newGeno[coord] <- -newGeno[coord]
		# Integrate into map and genoMat: sort genoMat then reform it
		allCumuPos <- c(cumuPos , speciesData$map[,3])
		# Have had a bug where cumuPos contains the same value more than once (!)
		dupCP <- duplicated(allCumuPos)
		nSame <- sum(dupCP)
		if (nSame){
			allCumuPos[dupCP] <- allCumuPos[dupCP] + (rbinom(nSame, 1, 0.5) * 2 - 1) / 1e6
		}
		locOrd <- order(allCumuPos)
		locRnk <- rank(allCumuPos)
		speciesData$map <<- rbind(cbind(Chr, Pos, cumuPos), speciesData$map)[locOrd,]
		genoMat <- cbind(newGeno , rbind(genoMat[,1:nLoc * 2 - 1], genoMat[,1:nLoc * 2]))[,locOrd]
		genoMat <- matrix(genoMat, nrow=popSize)
		# mutInfo <- rbind(indNumVec, matOrPat, ancAllele, locRnk[1:nMut])
		# rownames(mutInfo) <- c("mutInd", "mutGam", "ancAllele", "locIdx")
		mutInfo <- rbind(indNumVec, matOrPat, locRnk[1:nMut])
		rownames(mutInfo) <- c("mutInd", "mutGam", "locIdx")
	} else mutInfo <- NULL
	return(list(pedigree=parMat, genoMat=genoMat, mutInfo=mutInfo))
}#END random mate diploid

# Take two inbred parents, make an F1, make a gamete from the F1, and double it
# markers has genotype vectors in rows and markers in columns
# Warning: this is good only if you know you are starting from complete inbreds
makeProgenyDH <- function(parents, markers, map, chrMax){
	nLoc <- nrow(map)
	F1 <- numeric(2*nLoc)
	F1[1:nLoc*2-1] <- markers[parents[1],1:nLoc*2]
	F1[1:nLoc*2] <- markers[parents[2],1:nLoc*2]
	return(doubleGamete(makeGamete(F1, map, chrMax)))
}

# Use the makeProgenyDH function to randomly-mate a set of individuals
# Warning: this is good only if you know you are starting from complete inbreds
# Take mutation function out of makeGamete. First mate without mutation.
# Second add in mutation
randomMateDH <- function(markers, map, popSize, chrMax){
	lengthChr <- 150 # WARNING 150 cm hardwired
	allMutLoc <<- NULL # NOTE: I still need to be able to recreate this
	nPar <- nrow(markers)
	parMat <- sapply(rep(nPar, popSize), sample, size=2)
	genoMat <- t(apply(parMat, 2, makeProgenyDH, markers, map, chrMax))
	rownames(genoMat) <- paste("DH", 1:popSize, sep="")
	# Now add in newly mutated loci
	if (!is.null(speciesData$mutParms)){
		mutPerInd <- rpois(popSize, speciesData$mutParms$mutNum)
		nMut <- sum(mutPerInd)
		locCumuPos <- runif(nMut, 0, nChr * lengthChr)
		locChr <- locCumuPos %/% (lengthChr) + 1
		locPos <- locCumuPos %% lengthChr
		locPos <- speciesData$unifToPos(locPos)
		# Make loci with one mutant individual
		indNumVec <- unlist(sapply(1:popSize, function(indNum) rep(indNum, mutPerInd[indNum])))
		coord <- cbind(indNumVec, 1:nMut)
		# ancAllele <- rbinom(nMut, 1, 0.5) * 2 - 1
		# newGeno <- matrix(rep(ancAllele, each=popSize), nrow=popSize)
		newGeno <- matrix(rep(-1, each=popSize), nrow=popSize) # -1 is the ancestral allele
		newGeno[coord] <- -newGeno[coord]
		# Integrate into map and genoMat
		allCumuPos <- c(locCumuPos , speciesData$map[,3])
		locOrd <- order(allCumuPos)
		locRnk <- rank(allCumuPos)
		speciesData$map <<- rbind(cbind(locChr, locPos, locCumuPos), speciesData$map)[locOrd,]
		genoMat <- cbind(newGeno, genoMat[,1:nLoc * 2])[,locOrd]
		genoMat <- matrix(rbind(genoMat, genoMat), nrow=popSize)
		# mutInfo <- rbind(indNumVec, ancAllele, locRnk[1:nMut])
		# rownames(mutInfo) <- c("mutInd", "ancAllele", "locIdx")
		mutInfo <- rbind(indNumVec, locRnk[1:nMut])
		rownames(mutInfo) <- c("mutInd", "locIdx")
	} else mutInfo <- NULL
	return(list(pedigree=parMat, genoMat=genoMat, mutInfo=mutInfo))
}#END random mate diploid

# Figure out what kind of progeny it is and randomly mate it
# Allowable progeny types are "DH"; "Outbred"; or "RIL,No." (where No. is the number of selfing generations, e.g., "RIL,4");
# WARNING! The pedigree is in two ROWS rather than the typical two columns
randomMate <- function(markers, map, popSize, chrMax, progType){
	# if you want to simulate RIL, figure out how many generations you have to inbreed
	if (substr(progType, 1, 3) == "RIL"){
		nGenInbreed <- as.numeric(strsplit(progType, ",")[[1]][2])
	}

	if (progType == "DH"){
		return(randomMateDH(markers, map, popSize, chrMax))
	} else if (progType == "Outbred"){
		return(randomMateDiploid(markers, map, popSize, chrMax))
	} else{
		return(randomMateRIL(markers, map, popSize, chrMax, nGenInbreed))
	}	
}

##########################################################################################
# 2. Objects to create variability from replication to replication
# Minimally, it includes:
# createSpeciesData -- a function to create species data (genome, genetic architecture, and breeding system)
# createBreedingData -- a function to create the base population for each rep
# In this iteration, the philosophy changes quite a bit: there is one category of locus
# with loci sometimes being assigned as QTL
# Only create as many loci as will be polymorphic
# Does not return a genotype matrix, only a matrix of founder haplotypes
# Percent four locus cat has to sum to 1 and be in this order
# causal_Obs, causal_notObs, notCausal_Obs, notCausal_notObs
# 4           3               2             1               : locusType numbers
createSpeciesData <- function(effPopSize, mutNum, perc4LocCat, seed=round(runif(1, 0, 1e9))){
	set.seed(seed)
		
	nChr <- 7
	lengthChr <- 150 # in cM
	progType <- "Outbred"
	
	# Calculation for the number of loci that will be polymorphic initially
	nLoc <- ceiling(4 * effPopSize * mutNum * sum(1 / 1:(effPopSize - 1)))

  # percQTL: percent QTL among all loci
	percQTL  <-  sum(perc4LocCat[1:2])
  
	# Start with a sample with reasonable polymorphism from a coalescent
	piecesPerM <- 2000
	nPiecesPerChr <- lengthChr / 100 * piecesPerM
	recBTpieces <- 1 / piecesPerM
	coalSim <- getCoalescentSim(effPopSize=effPopSize, nMrkOrMut=nLoc, nChr=nChr, nPiecesPerChr=nPiecesPerChr, recBTpieces=recBTpieces, minMAF=0.01, seed=seed)
	fHaps <- coalSim$markers
	# No longer: Scramble ancestral state
	# ancestralState <- rbinom(nLoc, 1, 0.5) * 2 - 1
	# fHaps[,ancestralState == 1] <- 1 - fHaps[,ancestralState == 1]
	fHaps <- fHaps * 2 - 1
	
  # Set up the map: set the first and last loci on each chr to pos 0 and lengthChr
  fMap <- coalSim$map
  for (chr in 1:nChr){
    firstChrLoc <- which(fMap[,"Chr"] == chr & fMap[,"Pos"] == min(fMap[fMap[,"Chr"] == chr,"Pos"]))
    fMap[firstChrLoc,"Pos"] <- 0
    lastChrLoc <- which(fMap[,"Chr"] == chr & fMap[,"Pos"] == max(fMap[fMap[,"Chr"] == chr,"Pos"]))
    fMap[lastChrLoc,"Pos"] <- lengthChr
  }
	fMap <- makeChrMax(fMap[,1:2])
	fChrMax <- fMap$chrMax
	fMap <- fMap$map

	# Function to sample the QTL effects
	sampleEffects <- function(geneActionList, locusList){
		nEffects <- length(locusList)
		effects <- rgamma(nEffects, 0.4)
		effTooBig <- effects > median(effects) # Simple censoring to get to equilibrium faster
		effects[effTooBig] <- rgamma(sum(effTooBig), 0.4)
		return(effects * (rbinom(nEffects, 1, 0.5) * 2 - 1))
	}
  
	# sample which loci are candidates to be QTL
	genArch  <- sampleGeneticArchitecture(nLoc, nLoc * percQTL, effectStdDev=sampleEffects, meanInteractionOrder=0, probDominance=0, inbredReference=FALSE, replace=FALSE)
	genArch$locusList <- sort(genArch$locusList)
	
	mutParms <- list(mutNum=mutNum, perc4LocCat=perc4LocCat)
	
	# Function to translate uniform positions to on with restricted rec in centromere
	unifToPos <- function(thePos){ # WARNING: assumes 0 and 150 are min & max of position
		cut <- 0.25 # 1 - 2*cut goes into the central 10%
		newPos <- thePos / 150
		belCut <- newPos <= cut
		aboCut <- newPos > (1 - cut)
		newPos[belCut] <- newPos[belCut] * 0.45 / cut
		newPos[aboCut] <- 0.55 + (newPos[aboCut] - (1 - cut)) * 0.45 / cut
		midPos <- newPos[!(belCut | aboCut)]
		midPos <- 0.01 + (midPos - cut) * 0.98 / (1 - 2 * cut)
		invLogistic <- function(y) return(log(y / (1 - y)))
		midPos <- invLogistic(midPos) / invLogistic(0.99)
		midPos <- 0.45 + (midPos + 1) * 0.1 / 2
		newPos[!(belCut | aboCut)] <- midPos
		return(newPos * 150)
	}
  # For now I just want this to be simple
  unifToPos <- function(thePos) return(thePos)
  
	# return(list(map=fMap, chrMax=fChrMax, genArch=genArch, founderHaps=fHaps, ancestralState=ancestralState, progType=progType, randSeed=seed, mutParms=mutParms, percQTL=percQTL, unifToPos=unifToPos, cycle=-1))
	return(list(map=fMap, chrMax=fChrMax, genArch=genArch, founderHaps=fHaps, progType=progType, randSeed=seed, mutParms=mutParms, percQTL=percQTL, unifToPos=unifToPos, cycle=-1))
}#END createSpeciesData

# Function to create the base population
createBreedingData <- function(nSelCan, nToSelect, nStoredGen, seed=round(runif(1, 0, 1e9))){
	set.seed(seed)

	# Heritability will emerge as a result of other parameters
	stdDevErr <- 1 # Error variance = 1 will be fixed across all simulations
		
	# Create the base population
	ftwoColMrk <- t(apply(speciesData$founderHaps, 1, doubleGamete))
	genoMat <- randomMate(ftwoColMrk, speciesData$map, nSelCan, speciesData$chrMax, speciesData$progType)
	mutInfo <- genoMat$mutInfo
	genoMat <- genoMat$genoMat
	nLoc <- nrow(speciesData$map) # Map is updated so this is the total # of loci
	# newAncAllele <- mutInfo["ancAllele",]
	nMut <- ncol(mutInfo)
	newLoc <- mutInfo["locIdx",]
	oldLoc <- (1:nLoc)[-newLoc]
	nNewQTL <- floor(nMut * speciesData$percQTL)
	# ancAllele <- c(speciesData$ancestralState, newAncAllele)[order(c(oldLoc, newLoc))]
	newQTLidx <- sample(newLoc, nNewQTL)
	newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.05) * 2 - 1)
	# newEff <- newEff * ancAllele[newQTLidx]
	oldQTLidx <- oldLoc[speciesData$genArch$locusList]
	allQTL <- c(oldQTLidx, newQTLidx)
	allQTLord <- order(allQTL)
	speciesData$genArch$effects <<- c(speciesData$genArch$effects, newEff)[allQTLord]
	speciesData$genArch$locusList <<- allQTL[allQTLord]
  
  # Figure out which loci are observed markers
	# Percent four locus cat has to sum to 1 and be in this order
	# causal_Obs, causal_notObs, notCausal_Obs, notCausal_notObs
	percObsQTL <- speciesData$mutParms$perc4LocCat
	percObsNonQTL <- percObsQTL[3] / sum(percObsQTL[3:4])
	percObsQTL <- percObsQTL[1] / sum(percObsQTL[1:2])
	obsQTL <- sample(allQTL, round(percObsQTL * length(allQTL)))
  obsNonQTL <- sample((1:nLoc)[-allQTL], round(percObsNonQTL * (nLoc - length(allQTL))))
  observedLoc <- sort(c(obsQTL, obsNonQTL))
  
	# How big can the training population become?
	maxRecordNum <- nStoredGen * nSelCan
	# These are the first selection candidates
	records <- data.frame(ID=1:nSelCan, MID=0, PID=0, cycle=1, genoVal=NA, phenoVal=NA, genoHat=NA)
	# return(list(genoMat=genoMat, observedLoc=observedLoc, records=records, ancAllele=ancAllele, allSelectSet=NULL, IDtoRow=1:nSelCan, stdDevErr=initStdDevErr, nSelCan=nSelCan, nToSelect=nToSelect, nStoredGen=nStoredGen, maxRecordNum=maxRecordNum, randSeed=seed))
	return(list(genoMat=genoMat, observedLoc=observedLoc, records=records, allSelectSet=NULL, IDtoRow=1:nSelCan, stdDevErr=initStdDevErr, nSelCan=nSelCan, nToSelect=nToSelect, nStoredGen=nStoredGen, maxRecordNum=maxRecordNum, randSeed=seed))
}#END createBreedingData

####################################################################################
# Now define the breedFuncList that has the name of all universal objects
####################################################################################
breedFuncList <- ls()
