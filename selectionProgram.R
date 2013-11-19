###############################################################
### Setting numbers for this selection program
###############################################################
settings <- c(2001, 2002)
allNcycles <- rep(400, length(settings))

settingReminder <- function(){
  cat("Using 4LocCat")
  cat("Setting 2001 Baseline: error std dev=1, mutation effect=1, prob. fav.=0.5, nSelCan=200, causal unobserved with mut. rate=1, non-causal observed mut. rate=8, non-causal non-observed mut. rate=1, markers not standardized, all progeny phenotyped", "\n")
  cat("Setting 2002: error std dev=1, mutation effect=1, prob. fav.=0.5, nSelCan=200, causal observed with mut. rate=1, non-causal observed mut. rate=0, non-causal non-observed mut. rate=1, markers not standardized, all progeny phenotyped", "\n")
}

###############################################################
#					Initialize
# Function to create or load speciesData or breedingData
###############################################################
initializeCreateData <- function(replication){
	switchCyc <<- nCycles / 2
  phenoOnlyParents <<- NULL # For settings with phenoOnlyParents, have nToSelect be 0.5 * nSelCan
	if (makeFounderFiles){
    # Percent four locus cat has to sum to 1 and be in this order
    # causal_Obs, causal_notObs, notCausal_Obs, notCausal_notObs
    if (setting == 2001){
      mutNum <- 10
      perc4LocCat <- c(0, 0.1, 0.8, 0.1)
    }
    if (setting == 2002){
      mutNum <- 2
      perc4LocCat <- c(0.5, 0, 0, 0.5)      
    }
    
		nSelCan <- 200
		nToSelect <- floor(runif(1, nSelCan*0.3, nSelCan*0.55))
		effPopSize <- nToSelect
		nStoredGen <- 8
		
		speciesData <<- createSpeciesData(effPopSize, mutNum, perc4LocCat, founderSeed)

		breedingData <<- createBreedingData(nSelCan, nToSelect, nStoredGen=nStoredGen, founderSeed)
    
    afterUpdate <- floor(runif(1, nSelCan*0.3, nSelCan*0.55))
    afterThat <- floor(runif(1, nSelCan*0.4, nSelCan*0.65))
    breedingData$nToSelectGS <<- c(afterUpdate, afterThat)
    breedingData$cumulativeIncidences <<- numeric(nrow(speciesData$map))
    
		save(replication, founderSeed, founderFile, speciesData, breedingData, file=founderFile)
	} else{ # Don't make founder file but load it
		load(founderFile, .GlobalEnv)
	}

	return(founderFile)
}

# WARNING throwing caution to the wind and implementing various short cuts
# Assumes simplest strictly additive model
###############################################################
#					PHENOTYPE
# Function to phenotype should expect to receive breeding records
# with all information in place except no genotypic or phenotypic
# values for the individuals of the current cycle.
# The role of the function is to fill in the records$phenoVal
# variable appropriately
###############################################################
# Phenotyping for genomic selection.
phenotypeGenomic <- function(cycle, breedingData){
  # Calculate genotypic values
	genRows <- which(breedingData$records$cycle == cycle)
	qtl <- speciesData$genArch$locusList
	qtlDosage <- (breedingData$genoMat[genRows, qtl*2 - 1] + breedingData$genoMat[genRows, qtl*2]) / 2
	breedingData$records$genoVal[genRows] <- speciesData$genArch$genoBase + c(qtlDosage %*% speciesData$genArch$effects)
	# No phenotyping at the transition or in odd cycles for GS
	if (cycle == switchCyc | (cycle > switchCyc & cycle %% 2 == 1)){ 
		phenRows <- 0
	} else{ # Here we are phenotyping
		if (cycle < switchCyc){
			phenRows <- which(breedingData$records$cycle == cycle)
			tmp <- calcPhenotype(breedingData$records$genoVal[phenRows], breedingData$stdDevErr, heritability=FALSE)
		} else{
      # Phenotype only latest parents with lower error: lower error
			if (setting %in% phenoOnlyParents){ 
				phenRows <- breedingData$selectSet
				tmp <- calcPhenotype(breedingData$records$genoVal[phenRows], breedingData$stdDevErr / sqrt(2), heritability=FALSE)
			} else{ # Phenotype all prior selection candidates
				phenRows <- which(breedingData$records$cycle %in% (cycle - (1:2)))
				tmp <- calcPhenotype(breedingData$records$genoVal[phenRows], breedingData$stdDevErr, heritability=FALSE)
			}
		}
		breedingData$records$phenoVal[phenRows] <- tmp$phenoVal
    # Keep track of how often each mutant has been phenotyped
		eo <- 1:nrow(speciesData$map) * 2 # WARNING: will need to refigure this out
		inst <- -breedingData$ancAllele * colSums(breedingData$genoMat[phenRows, eo-1] + breedingData$genoMat[phenRows, eo]) / 2 + length(phenRows)
		breedingData$cumulativeIncidences <- breedingData$cumulativeIncidences + inst
	}
	cat(cycle, "phenotypeGenomic", range(phenRows), "\n")
	return(breedingData)
}#END phenotypeGenomic

###############################################################
#					SELECT
# Function to select should expect to receive breeding records
# with all relevant information in place.
# The role of the function is to set the breedingData$selectSet
# variable
# Update marker effects cycle by cycle
# With no QTL in the model, I don't know what the correct expected
# effects and variances should be, so I use rrBLUP to estimate
# With QTL in the model, I cheat a little and use what I know from
# simulation to come up with those expectations:
# The true error variance is easy to come by.
# The true marker effect variance I think I can figure out.
# Variance of gamma is 0.4, then there is the variance from
# whether it's + or - which will be 0.05 * 0.95 * 0.8^2
# WARNING: assumes percQTL = 0.5 to calculate eMrkVar
# WARNING: I am dividing mrkDosage by 2 to fit the A.mat function
###############################################################
# On even cycles update the model. On odd cycles, use the marker effect estimates.
# Genomic selection where loci are split into four categories, all combinations of observed vs not and causal vs not.
selectGenomic4LocCat <- function(cycle, breedingData){
  thisCycle <- which(breedingData$records$cycle == cycle)
  if (cycle < switchCyc){ # Phenotypic selection
    cat(cycle, "selectPhenotypic_4LC", "\n")
    selInd <- order(breedingData$records$phenoVal[thisCycle], decreasing=TRUE)
    breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelect]]
  } else{
    mrk <- breedingData$observedLoc
    nMrk <- length(mrk)
    eo <- mrk * 2
    if (cycle %% 2 == 0){ # Calculate effect GEBVs with new phenotypes
      if (cycle == switchCyc){
        # Create marker effect and variance estimates coming out of phenotypic selection
        newPheno <- which(!is.na(breedingData$records$phenoVal))
        cat(cycle, "selectGenomic_InitialModel_4LC", length(newPheno), "\n")
        mrkDosage <- (breedingData$genoMat[newPheno,eo - 1] + breedingData$genoMat[newPheno,eo]) / 2
        mixedSolveOut <- mixed.solve(y=breedingData$records$phenoVal[newPheno], mrkDosage)
        breedingData$mrkEffEst <- mixedSolveOut$u
        breedingData$mrkEffVar <- rep(mixedSolveOut$Vu, nMrk)
        breedingData$meanPhen <- mixedSolveOut$beta
      } else{
        errVar <- breedingData$stdDevErr
        newPheno <- which(breedingData$records$cycle %in% (cycle - (1:2)))
        newPheno <- newPheno[!is.na(breedingData$records$phenoVal[newPheno])]
        cat(cycle, "selectGenomic_ModelUpdate_4LC", length(newPheno), "\n")
        mrkDosage <- (breedingData$genoMat[newPheno,eo - 1] + breedingData$genoMat[newPheno,eo]) / 2
        mrkEffEst <- breedingData$mrkEffEst
        # I am just doing diagonals on the variance even though I should do the whole matrix
        # How big of a problem is that?
        mrkEffVar <- breedingData$mrkEffVar
        covPhenMrk <- sapply(1:ncol(mrkDosage), FUN=function(col) mrkDosage[,col] * mrkEffVar[col]) # t(cov(beta, y))
        varPhen <- tcrossprod(covPhenMrk, mrkDosage) + diag(errVar, length(newPheno))
        regCoef <- t(solve(varPhen, covPhenMrk))
        breedingData$mrkEffEst <- mrkEffEst + regCoef %*% (breedingData$records$phenoVal[newPheno] - breedingData$meanPhen - mrkDosage %*% mrkEffEst)
        breedingData$mrkEffVar <- mrkEffVar - diag(regCoef %*% covPhenMrk)
      }
    } else{ # Calculate effect GEBVs with previously estimated effects
      cat(cycle, "selectGenomicRR_PreviousModel_4LC", "\n")
    }
    mrkDosage <- (breedingData$genoMat[thisCycle,eo - 1] + breedingData$genoMat[thisCycle,eo]) / 2
    breedingData$records$genoHat[thisCycle] <- mrkDosage %*% breedingData$mrkEffEst
    selInd <- order(breedingData$records$genoHat[thisCycle], decreasing=TRUE)
    breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelectGS[cycle %% 2 + 1]]]
  }#END doing GS
  breedingData$allSelectSet <- c(breedingData$allSelectSet, breedingData$selectSet)
  breedingData$selGenoVal <- mean(breedingData$records$genoVal[breedingData$selectSet])
  return(breedingData)
}#END selectGenomic4LocCat

# For setting xxxx where I want to use straight up ridge regression to compare with the updating method
selectGenomic4LocCatRR <- function(cycle, breedingData){
	thisCycle <- which(breedingData$records$cycle == cycle)
	if (cycle < switchCyc){ # Phenotypic selection
		cat(cycle, "selectPhenotypic_4LocCatRR", "\n")
		selInd <- order(breedingData$records$phenoVal[thisCycle], decreasing=TRUE)
		breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelect]]
	} else{
	  mrk <- breedingData$observedLoc
	  nMrk <- length(mrk)
	  eo <- mrk * 2
	  if (cycle %% 2 == 0){ # Calculate effect GEBVs with new phenotypes
		  # Create marker effect and variance estimates coming out of phenotypic selection
		  newPheno <- which(!is.na(breedingData$records$phenoVal))
		  cat(cycle, "selectGenomic4LocCatRR", length(newPheno), "\n")
		  mrkDosage <- (breedingData$genoMat[newPheno,eo - 1] + breedingData$genoMat[newPheno,eo]) / 2
		  mixedSolveOut <- mixed.solve(y=breedingData$records$phenoVal[newPheno], mrkDosage)
		  breedingData$mrkEffEst <- mixedSolveOut$u
		  breedingData$mrkEffVar <- rep(mixedSolveOut$Vu, nMrk)
		  breedingData$meanPhen <- mixedSolveOut$beta
		} else{ # Calculate effect GEBVs with previously estimated effects
			cat(cycle, "selectGenomic4LocCatRR_PreviousModel", "\n")
		}
	  mrkDosage <- (breedingData$genoMat[thisCycle,eo - 1] + breedingData$genoMat[thisCycle,eo]) / 2
	  breedingData$records$genoHat[thisCycle] <- mrkDosage %*% breedingData$mrkEffEst
	  selInd <- order(breedingData$records$genoHat[thisCycle], decreasing=TRUE)
	  breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelectGS[cycle %% 2 + 1]]]
	}
	breedingData$allSelectSet <- c(breedingData$allSelectSet, breedingData$selectSet)
	breedingData$selGenoVal <- mean(breedingData$records$genoVal[breedingData$selectSet])
	return(breedingData)
}#END selectGenomicObsQTLrr

###############################################################
#					ANALYZE
# This needs to go after selection but before cross & advance
# Analyze the results of one cycle
###############################################################
analysis4LocCat <- function(cycle, breedingData){
  cat(cycle, "analysis4LocCat", "\n")
  thisCycle <- which(breedingData$records$cycle == cycle)
  thin <- nCycles / 40
  # Calculate equiVar; segmentVar; chromosomeVar; genoVar
  nLoc <- nrow(speciesData$map)
  eo <- 1:nLoc*2
  locDosage <- (breedingData$genoMat[thisCycle, eo - 1] + breedingData$genoMat[thisCycle, eo]) / 2
  qtl <- speciesData$genArch$locusList
  qtlDosage <- locDosage[,qtl]
  qtlDosageVar <- apply(qtlDosage, 2, var)
  totEquiVar <- c(crossprod(speciesData$genArch$effects^2, qtlDosageVar))
  mrk <- breedingData$observedLoc
  
  # Calculate the variance of a segment (could be a whole chromosome)
  # NOTE: assumes that the segment does NOT overlap two chromosomes
  # Assume the simplest additive gene action for now
  # Only do this once equilibrium reached (see below)
  # WARNING depends on calculating equiVar before hand
  calcSegmentVar <- function(startPos, segSize, calcEquiVar=FALSE, calcSegAcc=FALSE){
    segGenoVals <- NA
    if (calcEquiVar) equiVar <- NA else equiVar <- NULL
    if (calcSegAcc) segAcc <- NA else segAcc <- NULL
    segLoci <- which(speciesData$map[,3] >= startPos & speciesData$map[,3] < startPos + segSize) # segLoci in locus order
    if (length(segLoci > 0)){
      segQTL <- which(qtl %in% segLoci) # segQTL in qtl order
      segGenoVals <- qtlDosage[,segQTL, drop=FALSE] %*% speciesData$genArch$effects[segQTL]
      if (calcEquiVar) equiVar <- crossprod(speciesData$genArch$effects[segQTL]^2, qtlDosageVar[segQTL])
      if (calcSegAcc){ # Here is where things need changing for markers
        segGenoHat <- locDosage[, intersect(mrk, segLoci), drop=FALSE] %*% breedingData$mrkEffEst[mrk %in% segLoci]
        segAcc <- cor(segGenoVals, segGenoHat)
      }
    }
    return(c(var(segGenoVals), equiVar, segAcc))
  }
  
  # Calculate segment variances every thin cycles on the whole genome
  chrVar <- NA
  segVar10 <- matrix(NA, 3, 105)
  if (cycle %% thin == 0){
    chrVar <- sapply(0:6 * 150, calcSegmentVar, segSize=150) # WARNING chr size of 150 hardwired
    seg10 <- rep(0:6 * 150, each=15) + 0:14 * 10
    segVar10 <- sapply(seg10, calcSegmentVar, segSize=10, calcEquiVar=TRUE, calcSegAcc=cycle >= switchCyc)
    breedingData$localStats <- c(breedingData$localStats, list(chrVar, segVar10))
  }
  # How many generations does it take to reach 5% allele frequency?
  # cycle; Chr and Pos; nGenerations polymorphic; effect; fixedAllele
  if (cycle == 1){
    breedingData$wasAbove05 <- logical(nLoc)
    breedingData$nGenPoly <- rep(1, nLoc)
  } 
  # mutFreq <- (-apply(locDosage, 2, mean) * breedingData$ancAllele + 1) / 2
  mutFreq <- (apply(locDosage, 2, mean) + 1) / 2
  newlyAbove05 <- which(!breedingData$wasAbove05 & mutFreq > 0.05)
  breedingData$wasAbove05[newlyAbove05] <- TRUE
  effects <- numeric(length(newlyAbove05)) # Effect of the mutation
  # effects[newlyAbove05 %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyAbove05] * -breedingData$ancAllele[intersect(newlyAbove05, qtl)]
  effects[newlyAbove05 %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyAbove05]
  nGenPoly <- breedingData$nGenPoly[newlyAbove05]
  locusType <- (newlyAbove05 %in% qtl)*2 + (newlyAbove05 %in% breedingData$observedLoc) + 1
  above05Info <- cbind(cycle, speciesData$map[newlyAbove05, c("Chr", "Pos")], nGenPoly, locusType)
  breedingData$above05Info <- rbind(breedingData$above05Info, above05Info)
  
  # Collect information on polymorphism
  # Number of loci polymorphic in the different categories. locDosage is for thisCycle
  locPoly <- which(apply(locDosage, 2, sd) > 0)
  mrkQtl <- intersect(locPoly, intersect(mrk, qtl))
  notMrkQtl <- intersect(locPoly, setdiff(qtl, mrk))
  mrkNotQtl <- intersect(locPoly, setdiff(mkr, qtl))
  notNot <- setdiff(locPoly, union(mrk, qtl))
  
  # Collect some information on the trait this cycle
  genoMean <- mean(breedingData$records$genoVal[thisCycle])
  genoVar <- var(breedingData$records$genoVal[thisCycle])
  thisCycleAccuracy <- cor(breedingData$records$genoVal[thisCycle], breedingData$records$genoHat[thisCycle])
  breedingData$history <- rbind(breedingData$history, data.frame(cycle=cycle, genoBase=speciesData$genArch$genoBase, nQTLpoly=length(c(mrkQtl, notMrkQtl)), nQTL=length(qtl), nMrkPoly=length(c(mrkQtl, mrkNotQtl)), nMrk=length(mrk), meanGenoVal=genoMean, accuracy=thisCycleAccuracy, genoVar=genoVar, totChrVar=sum(chrVar, na.rm=TRUE), totSegVar=sum(segVar10[1,], na.rm=TRUE), equiVar=totEquiVar, selGenoValMain=breedingData$selGenoVal), mrkQtlPoly=length(mrkQtl), notMrkQtlPoly=length(notMrkQtl), mrkNotQtlPoly=length(mrkNotQtl), neutralPoly=length(notNot))
  return(breedingData)
}#END analysis4LocCat

###############################################################
#					CROSS and ADVANCE
# pull out the select set
# create the next generation from the select set
# rbind to genoMat
# rbind to records
###############################################################
crossAdvanceMutUpdateUnobsQ <- function(cycle, breedingData){
	cat(cycle, "crossAdvanceMutate_UpdateUnobsQ", "\n")
	
	selectSet <- breedingData$selectSet
	nLoc <- nrow(speciesData$map)
	eo <- 1:nLoc * 2
	qtl <- speciesData$genArch$locusList
	mrk <- (1:nLoc)[-qtl]
	
	# Loci only become monomorphic and fall out if a cycle gets dropped
	dropCyc <- (cycle > 1 & cycle < switchCyc - breedingData$nStoredGen + 2) | (cycle >=switchCyc & cycle %% 2 == 0)
	if (dropCyc){
		cat(cycle, "Dropping past cycles", "\n")
		toDrop <- which(breedingData$records$cycle < cycle)
		breedingData$genoMat <- breedingData$genoMat[-toDrop,]
		selectSet <- selectSet - length(toDrop)
		breedingData$selectSet <- selectSet
		breedingData$records <- breedingData$records[-toDrop,]
		monomorphic <- which(apply(rbind(breedingData$genoMat[, eo], breedingData$genoMat[, eo-1]), 2, sd) == 0)
		# Take monomorphic out of the map
		speciesData$map <<- speciesData$map[-monomorphic,]
		# out of the QTL effects
		monoQTL <- which(qtl %in% monomorphic)
		fixQTLeff <- speciesData$genArch$effects[monoQTL]
		fixQTLal <- breedingData$genoMat[1, qtl[monoQTL]*2]
		shiftToGenoBase <- 2 * crossprod(fixQTLal, fixQTLeff)
		speciesData$genArch$effects <<- speciesData$genArch$effects[-monoQTL]
		# out of marker effect estimates
		monoMrk <- which(mrk %in% monomorphic)
		breedingData$mrkEffEst <- breedingData$mrkEffEst[-monoMrk]
		breedingData$mrkEffVar <- breedingData$mrkEffVar[-monoMrk]
		# out of the locus data
		breedingData$genoMat <- breedingData$genoMat[, -(rep(monomorphic*2, each=2)-0:1)]
		# out of nGenPoly
		breedingData$nGenPoly <- breedingData$nGenPoly[-monomorphic]
		# breedingData$ancAllele <- breedingData$ancAllele[-monomorphic]
		# Renumber qtl and mrk
		qtl <- which((1:nLoc)[-monomorphic] %in% qtl)
		mrk <- which((1:nLoc)[-monomorphic] %in% mrk)
	} else shiftToGenoBase <- 0
	
	# Make progeny with the remaining loci
	selGenoMat <- breedingData$genoMat[selectSet,]
	pedNgeno <- randomMate(selGenoMat, speciesData$map, breedingData$nSelCan, speciesData$chrMax, speciesData$progType)
	mutInfo <- pedNgeno$mutInfo
	newIDs <- max(breedingData$records$ID) + 1:breedingData$nSelCan
	rownames(pedNgeno$genoMat) <- newIDs
	newRecords <- data.frame(ID=newIDs, MID=selectSet[pedNgeno$pedigree[1,]], PID=selectSet[pedNgeno$pedigree[2,]], cycle=cycle + 1, genoVal=NA, phenoVal=NA, genoHat=NA)
	breedingData$records <- rbind(breedingData$records, newRecords)
	
	nLoc <- nrow(speciesData$map)
	eo <- 1:nLoc * 2
	speciesData$nLoc <<- c(speciesData$nLoc, nLoc)
	polyNow <- apply(rbind(pedNgeno$genoMat[, eo], pedNgeno$genoMat[, eo-1]), 2, sd) > 0
	
	# new genoMat has extra columns relative to old genoMat so add columns to the old
	if (!is.null(mutInfo)){
		nMut <- ncol(mutInfo)
		# mutAncAllele <- mutInfo["ancAllele",]
		newLoc <- mutInfo["locIdx",]
		oldLoc <- (1:nLoc)[-newLoc]
		allLocOrd <- order(c(oldLoc, newLoc))
		newGM <- breedingData$genoMat
		# newGM <- cbind(newGM, matrix(rep(mutAncAllele, each=nrow(newGM)*2), nrow=nrow(newGM)))
		newGM <- cbind(newGM, matrix(rep(-1, each=nrow(newGM)*2), nrow=nrow(newGM)))
		breedingData$genoMat <- rbind(newGM[,rep(allLocOrd, each=2)*2 - 1:0], pedNgeno$genoMat)
		# breedingData$ancAllele <- c(breedingData$ancAllele, mutAncAllele)[allLocOrd]
	
		# The map and the genoMat are taken care of; deal with QTL and effects
		nNewQTL <- floor(nMut * speciesData$percQTL)
		newQTLidx <- sample(nMut, nNewQTL)
		# qtlAncAl <- mutAncAllele[newQTLidx]
		# newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.95) * 2 - 1)
		# newEff <- newEff * qtlAncAl
		newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.05) * 2 - 1)
		allMutEff <- numeric(nMut)
		allMutEff[newQTLidx] <- newEff
		mutInfo <- rbind(mutInfo, allMutEff)
		newQTLidx <- newLoc[newQTLidx]
		oldQTLidx <- oldLoc[qtl]
		allQTL <- c(oldQTLidx, newQTLidx)
		allQTLord <- order(allQTL)
		speciesData$genArch$effects <<- c(speciesData$genArch$effects, newEff)[allQTLord]
		qtl <- allQTL[allQTLord]
		# shiftToGenoBase <- shiftToGenoBase - 2 * crossprod(qtlAncAl, newEff)
		shiftToGenoBase <- shiftToGenoBase + 2 * sum(newEff)
		
		# Do the same for markers
		nNewMrk <- nMut - nNewQTL
		newMrkIdx <- setdiff(newLoc, newQTLidx)
		oldMrkIdx <- oldLoc[mrk]
		allMrkOrd <- order(c(oldMrkIdx, newMrkIdx))
		breedingData$mrkEffEst <- c(breedingData$mrkEffEst, rep(0, nNewMrk))[allMrkOrd]
		eMrkVar <- 0.2476
		breedingData$mrkEffVar <- c(breedingData$mrkEffVar, rep(eMrkVar, nNewMrk))[allMrkOrd]
		
		# Take care of vector of nGenerations polymorphic
		if (cycle == 1) breedingData$nGenPoly <- rep(1, length(oldLoc))
		breedingData$nGenPoly <- c(breedingData$nGenPoly, numeric(nMut))[allLocOrd] + polyNow

	} else{# END mutations _did_ happen
		breedingData$genoMat <- rbind(breedingData$genoMat, pedNgeno$genoMat)
		breedingData$nGenPoly <- breedingData$nGenPoly + polyNow
	}
	
	speciesData$genArch$locusList <<- qtl
	speciesData$genArch$genoBase <<- speciesData$genArch$genoBase + shiftToGenoBase
	breedingData$mutHist <- c(breedingData$mutHist, list(mutInfo))
	# For loci that just became monomorphic, save:
	# cycle; Chr and Pos; nGenerations polymorphic; effect; fixedAllele
	# Which were poly last generation, but not now?  Store the values for those.
	lastCyc <- breedingData$records$cycle == cycle # (Just created cycle + 1)
	polyBefore <- apply(rbind(breedingData$genoMat[lastCyc, eo], breedingData$genoMat[lastCyc, eo-1]), 2, sd) > 0
	newlyMono <- which(polyBefore & !polyNow)
	effects <- numeric(length(newlyMono)) # Effect of the mutation
	# effects[newlyMono %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyMono] * -breedingData$ancAllele[intersect(newlyMono, qtl)]
	effects[newlyMono %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyMono]
	nGenPoly <- breedingData$nGenPoly[newlyMono]
	# ancAllele <- breedingData$ancAllele[newlyMono]
	fixedAllele <- pedNgeno$genoMat[1, newlyMono*2]
	# monoInfo <- cbind(cycle, speciesData$map[newlyMono, c("Chr", "Pos")], nGenPoly, effects, ancAllele, fixedAllele)
	monoInfo <- cbind(cycle, speciesData$map[newlyMono, c("Chr", "Pos")], nGenPoly, effects, fixedAllele)
	breedingData$monoInfo <- rbind(breedingData$monoInfo, monoInfo)
	
	# save(breedingData, speciesData, file=paste("bDsD",cycle,"RData", sep="."))	
	return(breedingData)
}#END cross advance _unObs_ QTL

# Update means that you cumulate only at the end of phenotypic selection and
# the two cycles that are phenotyped every other cycle
crossAdvanceMutObsQ <- function(cycle, breedingData){
  cat(cycle, "crossAdvanceMutate_UpdateObsQ", "\n")
  
  selectSet <- breedingData$selectSet
  nLoc <- nrow(speciesData$map)
  eo <- 1:nLoc * 2
  qtl <- speciesData$genArch$locusList
  #debugPlace <<- 0
  #cat(debugPlace, "\n"); debugPlace <<- debugPlace+1
  
  # Loci only become monomorphic and fall out if a cycle gets dropped
  dropCyc <- (cycle > 1 & cycle < switchCyc - breedingData$nStoredGen + 2) | (cycle >= switchCyc & cycle %% 2 == 0)
  if (dropCyc){
    if (setting %in% c(421) & cycle >= switchCyc){
      toDrop <- which(breedingData$records$cycle < cycle - breedingData$nStoredGen + 2)
    } else toDrop <- which(breedingData$records$cycle < cycle)
    cat(cycle, "Dropping cycles", range(breedingData$records$cycle[toDrop]), "\n")
    breedingData$genoMat <- breedingData$genoMat[-toDrop,]
    selectSet <- selectSet - length(toDrop)
    breedingData$selectSet <- selectSet
    breedingData$records <- breedingData$records[-toDrop,]
    monomorphic <- which(apply(rbind(breedingData$genoMat[, eo], breedingData$genoMat[, eo-1]), 2, sd) == 0)
    # Take monomorphic out of the map
    speciesData$map <<- speciesData$map[-monomorphic,]
    # out of the QTL effects
    monoQTL <- which(qtl %in% monomorphic)
    fixQTLeff <- speciesData$genArch$effects[monoQTL]
    fixQTLal <- breedingData$genoMat[1, qtl[monoQTL]*2]
    shiftToGenoBase <- crossprod(fixQTLal, fixQTLeff)
    speciesData$genArch$effects <<- speciesData$genArch$effects[-monoQTL]
    # out of marker effect estimates
    breedingData$mrkEffEst <- breedingData$mrkEffEst[-monomorphic]
    breedingData$mrkEffVar <- breedingData$mrkEffVar[-monomorphic]
    # out of the locus data
    breedingData$genoMat <- breedingData$genoMat[, -(rep(monomorphic*2, each=2)-0:1)]
    # out of nGenPoly
    breedingData$nGenPoly <- breedingData$nGenPoly[-monomorphic]
    # breedingData$ancAllele <- breedingData$ancAllele[-monomorphic]
    breedingData$cumulativeIncidences <- breedingData$cumulativeIncidences[-monomorphic]
    # Renumber qtl and mrk
    qtl <- which((1:nLoc)[-monomorphic] %in% qtl)
  } else shiftToGenoBase <- 0
  
  # Make progeny with the remaining loci
  selGenoMat <- breedingData$genoMat[selectSet,]
  pedNgeno <- randomMate(selGenoMat, speciesData$map, breedingData$nSelCan, speciesData$chrMax, speciesData$progType)
  mutInfo <- pedNgeno$mutInfo
  newIDs <- max(breedingData$records$ID) + 1:breedingData$nSelCan
  rownames(pedNgeno$genoMat) <- newIDs
  newRecords <- data.frame(ID=newIDs, MID=selectSet[pedNgeno$pedigree[1,]], PID=selectSet[pedNgeno$pedigree[2,]], cycle=cycle + 1, genoVal=NA, phenoVal=NA, genoHat=NA)
  breedingData$records <- rbind(breedingData$records, newRecords)
  
  nLoc <- nrow(speciesData$map)
  eo <- 1:nLoc * 2
  speciesData$nLoc <<- c(speciesData$nLoc, nLoc)
  polyNow <- apply(rbind(pedNgeno$genoMat[, eo], pedNgeno$genoMat[, eo-1]), 2, sd) > 0
  
  # new genoMat has extra columns relative to old genoMat so add columns to the old
  if (!is.null(mutInfo)){
    nMut <- ncol(mutInfo)
    # mutAncAllele <- mutInfo["ancAllele",]
    newLoc <- mutInfo["locIdx",]
    oldLoc <- (1:nLoc)[-newLoc]
    allLocOrd <- order(c(oldLoc, newLoc))
    newGM <- breedingData$genoMat
    # newGM <- cbind(newGM, matrix(rep(mutAncAllele, each=nrow(newGM)*2), nrow=nrow(newGM)))
    newGM <- cbind(newGM, matrix(rep(-1, each=nrow(newGM)*2), nrow=nrow(newGM)))
    breedingData$genoMat <- rbind(newGM[,rep(allLocOrd, each=2)*2 - 1:0], pedNgeno$genoMat)
    # breedingData$ancAllele <- c(breedingData$ancAllele, mutAncAllele)[allLocOrd]
    if (exists("cumulativeIncidences", breedingData)) breedingData$cumulativeIncidences <- c(breedingData$cumulativeIncidences, numeric(nMut))[allLocOrd]
    
    # The map and the genoMat are taken care of; deal with QTL and effects
    nNewQTL <- floor(nMut * speciesData$percQTL)
    newQTLidx <- sample(nMut, nNewQTL)
    # qtlAncAl <- mutAncAllele[newQTLidx]
    # newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.95) * 2 - 1)
    # newEff <- newEff * qtlAncAl
    newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.05) * 2 - 1)
    allMutEff <- numeric(nMut)
    allMutEff[newQTLidx] <- newEff
    mutInfo <- rbind(mutInfo, allMutEff)
    newQTLidx <- newLoc[newQTLidx]
    oldQTLidx <- oldLoc[qtl]
    allQTL <- c(oldQTLidx, newQTLidx)
    allQTLord <- order(allQTL)
    speciesData$genArch$effects <<- c(speciesData$genArch$effects, newEff)[allQTLord]
    qtl <- allQTL[allQTLord]
    # shiftToGenoBase <- shiftToGenoBase - crossprod(qtlAncAl, newEff)
    shiftToGenoBase <- shiftToGenoBase + sum(newEff)
    
    # Adjust mrkEffEst, which has an est for _all_ loci
    # E(var(mrkEff)) = rate / rLQ + rate^2 / rLQ - (rate * (2*pFav - 1)/2)^2
    breedingData$mrkEffVar <- c(breedingData$mrkEffVar, rep(0.2476, nMut))[allLocOrd]
    if (setting %in% c(421)) priorMrkEff <- 0 else priorMrkEff <- 0.18
    # breedingData$mrkEffEst <- c(breedingData$mrkEffEst, priorMrkEff * mutAncAllele)[allLocOrd]
    breedingData$mrkEffEst <- c(breedingData$mrkEffEst, -priorMrkEff)[allLocOrd]
    # Take care of vector of nGenerations polymorphic
    if (cycle == 1) breedingData$nGenPoly <- rep(1, length(oldLoc))
    breedingData$nGenPoly <- c(breedingData$nGenPoly, numeric(nMut))[allLocOrd] + polyNow
    
  } else{#END mutations _did_ happen above
    # Below, mutations did not happen
    breedingData$genoMat <- rbind(breedingData$genoMat, pedNgeno$genoMat)
    breedingData$nGenPoly <- breedingData$nGenPoly + polyNow
  }
  
  speciesData$genArch$locusList <<- qtl
  speciesData$genArch$genoBase <<- speciesData$genArch$genoBase + shiftToGenoBase
  # WARNING this is wrong, I'm pretty sure: should just take into account the _markers_ that have become monomorphic
  if (exists("meanPhen", breedingData)) breedingData$meanPhen <- breedingData$meanPhen + shiftToGenoBase
  breedingData$mutHist <- c(breedingData$mutHist, list(mutInfo))
  # For loci that just became monomorphic, save:
  # cycle; Chr and Pos; nGenerations polymorphic; effect; fixedAllele
  # Which were poly last generation, but not now?  Store the values for those.
  lastCyc <- breedingData$records$cycle == cycle # (Just created cycle + 1)
  polyBefore <- apply(rbind(breedingData$genoMat[lastCyc, eo], breedingData$genoMat[lastCyc, eo-1]), 2, sd) > 0
  newlyMono <- which(polyBefore & !polyNow)
  effects <- numeric(length(newlyMono)) # Effect of the mutation
  # effects[newlyMono %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyMono] * -breedingData$ancAllele[intersect(newlyMono, qtl)]
  effects[newlyMono %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyMono]
  nGenPoly <- breedingData$nGenPoly[newlyMono]
  # ancAllele <- breedingData$ancAllele[newlyMono]
  fixedAllele <- pedNgeno$genoMat[1, newlyMono*2]
  # monoInfo <- cbind(cycle, speciesData$map[newlyMono, c("Chr", "Pos")], nGenPoly, effects, ancAllele, fixedAllele)
  monoInfo <- cbind(cycle, speciesData$map[newlyMono, c("Chr", "Pos")], nGenPoly, effects, fixedAllele)
  breedingData$monoInfo <- rbind(breedingData$monoInfo, monoInfo)
  
  # save(breedingData, speciesData, file=paste("bDsD",cycle,"RData", sep="."))	
  return(breedingData)
}#END cross advance obs QTL

crossAdvance4LocCat <- function(cycle, breedingData){
  cat(cycle, "crossAdvanceMutate_4LocCat", "\n")
  
  selectSet <- breedingData$selectSet
  nLoc <- nrow(speciesData$map)
  eo <- 1:nLoc * 2
  qtl <- speciesData$genArch$locusList
  mrk <- breedingData$observedLoc
  
  # Loci only become monomorphic and fall out if a cycle gets dropped
  # During GS, even cycles get phenotyped, so don't drop them during odd cycles
  dropCyc <- (cycle > 1 & cycle < switchCyc - breedingData$nStoredGen + 2) | (cycle >= switchCyc & cycle %% 2 == 0)
  if (dropCyc){
    toDrop <- which(breedingData$records$cycle < cycle)
    cat(cycle, "Dropping cycles", range(breedingData$records$cycle[toDrop]), "\n")
    breedingData$genoMat <- breedingData$genoMat[-toDrop,]
    selectSet <- selectSet - length(toDrop)
    breedingData$selectSet <- selectSet
    breedingData$records <- breedingData$records[-toDrop,]
    monomorphic <- which(apply(rbind(breedingData$genoMat[, eo], breedingData$genoMat[, eo-1]), 2, sd) == 0)
    # Take monomorphic out of the map
    speciesData$map <<- speciesData$map[-monomorphic,]
    # out of the QTL effects
    monoQTL <- which(qtl %in% monomorphic)
    fixQTLeff <- speciesData$genArch$effects[monoQTL]
    fixQTLal <- breedingData$genoMat[1, qtl[monoQTL]*2]
    shiftToGenoBase <- crossprod(fixQTLal, fixQTLeff)
    speciesData$genArch$effects <<- speciesData$genArch$effects[-monoQTL]
    # out of marker effect estimates CURRENT CODING
    # 1. What do I do to the overall population mean?
    monoMrk <- which(mrk %in% monomorphic)
    breedingData$mrkEffEst <- breedingData$mrkEffEst[-monoMrk]
    breedingData$mrkEffVar <- breedingData$mrkEffVar[-monoMrk]
    # out of the locus data
    breedingData$genoMat <- breedingData$genoMat[, -(rep(monomorphic*2, each=2)-0:1)]
    # out of nGenPoly
    breedingData$nGenPoly <- breedingData$nGenPoly[-monomorphic]
    breedingData$wasAbove05 <- breedingData$wasAbove05[-monomorphic]
    # breedingData$ancAllele <- breedingData$ancAllele[-monomorphic]
    breedingData$cumulativeIncidences <- breedingData$cumulativeIncidences[-monomorphic]
    # Renumber qtl and mrk
    qtl <- which((1:nLoc)[-monomorphic] %in% qtl)
  } else shiftToGenoBase <- 0
  
  # Make progeny with the remaining loci
  selGenoMat <- breedingData$genoMat[selectSet,]
  pedNgeno <- randomMate(selGenoMat, speciesData$map, breedingData$nSelCan, speciesData$chrMax, speciesData$progType)
  mutInfo <- pedNgeno$mutInfo
  newIDs <- max(breedingData$records$ID) + 1:breedingData$nSelCan
  rownames(pedNgeno$genoMat) <- newIDs
  newRecords <- data.frame(ID=newIDs, MID=selectSet[pedNgeno$pedigree[1,]], PID=selectSet[pedNgeno$pedigree[2,]], cycle=cycle + 1, genoVal=NA, phenoVal=NA, genoHat=NA)
  breedingData$records <- rbind(breedingData$records, newRecords)
  
  nLoc <- nrow(speciesData$map)
  eo <- 1:nLoc * 2
  speciesData$nLoc <<- c(speciesData$nLoc, nLoc)
  polyNow <- apply(rbind(pedNgeno$genoMat[, eo], pedNgeno$genoMat[, eo-1]), 2, sd) > 0
  
  # new genoMat has extra columns relative to old genoMat so add columns to the old
  if (!is.null(mutInfo)){
    nMut <- ncol(mutInfo)
    # mutAncAllele <- mutInfo["ancAllele",]
    newLoc <- mutInfo["locIdx",]
    oldLoc <- (1:nLoc)[-newLoc]
    allLocOrd <- order(c(oldLoc, newLoc))
    newGM <- breedingData$genoMat
    # newGM <- cbind(newGM, matrix(rep(mutAncAllele, each=nrow(newGM)*2), nrow=nrow(newGM)))
    newGM <- cbind(newGM, matrix(rep(-1, each=nrow(newGM)*2), nrow=nrow(newGM)))
    breedingData$genoMat <- rbind(newGM[,rep(allLocOrd, each=2)*2 - 1:0], pedNgeno$genoMat)
    # breedingData$ancAllele <- c(breedingData$ancAllele, mutAncAllele)[allLocOrd]
    if (exists("cumulativeIncidences", breedingData)) breedingData$cumulativeIncidences <- c(breedingData$cumulativeIncidences, numeric(nMut))[allLocOrd]
    
    # The map and the genoMat are taken care of; deal with QTL and effects
    nNewQTL <- floor(nMut * speciesData$percQTL)
    newQTLidx <- sample(nMut, nNewQTL)
    # qtlAncAl <- mutAncAllele[newQTLidx]
    # newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.95) * 2 - 1)
    # newEff <- newEff * qtlAncAl
    newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.05) * 2 - 1)
    allMutEff <- numeric(nMut)
    allMutEff[newQTLidx] <- newEff
    mutInfo <- rbind(mutInfo, allMutEff)
    newQTLidx <- newLoc[newQTLidx]
    oldQTLidx <- oldLoc[qtl]
    allQTL <- c(oldQTLidx, newQTLidx)
    allQTLord <- order(allQTL)
    speciesData$genArch$effects <<- c(speciesData$genArch$effects, newEff)[allQTLord]
    qtl <- allQTL[allQTLord]
    # shiftToGenoBase <- shiftToGenoBase - crossprod(qtlAncAl, newEff)
    shiftToGenoBase <- shiftToGenoBase + sum(newEff)
    
    # Adjust mrkEffEst, which has an est for _all_ loci
    # E(var(mrkEff)) = rate / rLQ + rate^2 / rLQ - (rate * (2*pFav - 1)/2)^2
    breedingData$mrkEffVar <- c(breedingData$mrkEffVar, rep(0.2476, nMut))[allLocOrd]
    if (setting %in% c(421)) priorMrkEff <- 0 else priorMrkEff <- 0.18
    # breedingData$mrkEffEst <- c(breedingData$mrkEffEst, priorMrkEff * mutAncAllele)[allLocOrd]
    breedingData$mrkEffEst <- c(breedingData$mrkEffEst, -priorMrkEff)[allLocOrd]
    # Take care of vector of nGenerations polymorphic
    breedingData$nGenPoly <- c(breedingData$nGenPoly, numeric(nMut))[allLocOrd] + polyNow
    breedingData$wasAbove05 <- c(breedingData$wasAbove05, logical(nMut))[allLocOrd]
    
  } else{#END mutations _did_ happen above
    # Below, mutations did not happen
    breedingData$genoMat <- rbind(breedingData$genoMat, pedNgeno$genoMat)
    breedingData$nGenPoly <- breedingData$nGenPoly + polyNow
  }
  
  speciesData$genArch$locusList <<- qtl
  speciesData$genArch$genoBase <<- speciesData$genArch$genoBase + shiftToGenoBase
  if (exists("meanPhen", breedingData)) breedingData$meanPhen <- breedingData$meanPhen + shiftToGenoBase
  breedingData$mutHist <- c(breedingData$mutHist, list(mutInfo))
  # For loci that just became monomorphic, save:
  # cycle; Chr and Pos; nGenerations polymorphic; effect; fixedAllele
  # Which were poly last generation, but not now?  Store the values for those.
  lastCyc <- breedingData$records$cycle == cycle # (Just created cycle + 1)
  polyBefore <- apply(rbind(breedingData$genoMat[lastCyc, eo], breedingData$genoMat[lastCyc, eo-1]), 2, sd) > 0
  newlyMono <- which(polyBefore & !polyNow)
  effects <- numeric(length(newlyMono)) # Effect of the mutation
  # effects[newlyMono %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyMono] * -breedingData$ancAllele[intersect(newlyMono, qtl)]
  effects[newlyMono %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyMono]
  nGenPoly <- breedingData$nGenPoly[newlyMono]
  # ancAllele <- breedingData$ancAllele[newlyMono]
  fixedAllele <- pedNgeno$genoMat[1, newlyMono*2]
  # monoInfo <- cbind(cycle, speciesData$map[newlyMono, c("Chr", "Pos")], nGenPoly, effects, ancAllele, fixedAllele)
  monoInfo <- cbind(cycle, speciesData$map[newlyMono, c("Chr", "Pos")], nGenPoly, effects, fixedAllele)
  breedingData$monoInfo <- rbind(breedingData$monoInfo, monoInfo)
  
  # save(breedingData, speciesData, file=paste("bDsD",cycle,"RData", sep="."))  
  return(breedingData)
}#END crossAdvance4LocCat

###############################################################
#					PROCESS RESULTS
# Function to process the results once they are out
# Some possibilities are:
# 1. Just compare gains from selection
# 2. Look at the extent to which genetic relatedness increased
# 3. Compile information about favorable alleles that were / were not fixed
###############################################################
processResults <- function(breedingData, speciesData, setting, replication){
	cat("processResults", setting, replication, "\n")
	needGenoVal <- is.na(breedingData$records$genoVal)
	qtl <- speciesData$genArch$locusList
	qtlDosage <- breedingData$genoMat[needGenoVal, qtl*2 - 1] + breedingData$genoMat[needGenoVal, qtl*2]
	breedingData$records$genoVal[needGenoVal] <- speciesData$genArch$genoBase + c(qtlDosage %*% speciesData$genArch$effects)
	breedingData$cycleMeans <- tapply(breedingData$records$genoVal, factor(breedingData$records$cycle), mean, na.rm=TRUE)	
	return(breedingData)
}

###############################################################
#					Functions to use in sapply
###############################################################

# Setting 301, 302 and 303:
if (settings %in% c(301, 302, 303)){
settingList <- list(initialize=initializeCreateData, phenotype=phenotypeGenomic, select=selectGenomicUnobsQ, analyze=analysisUnobsQ, crossAdvance=crossAdvanceMutUpdateUnobsQ, processResults=processResults)
}
# Setting 305, 401, 402:
if (settings %in% c(305, 401, 402, 1004)){
settingList <- list(initialize=initializeCreateData, phenotype=phenotypeGenomic, select=selectGenomicObsQ, analyze=analysisObsQ, crossAdvance=crossAdvanceMutObsQ, processResults=processResults)
}
# Setting 421:
if (settings %in% c(421)){
settingList <- list(initialize=initializeCreateData, phenotype=phenotypeGenomic, select=selectGenomicObsQrr, analyze=analysisObsQ, crossAdvance=crossAdvanceMutUpdateObsQ, processResults=processResults)
}
# Setting 325, 326, 327, 411:
if (settings %in% c(325, 326, 327, 411)){
settingList <- list(initialize=initializeCreateData, phenotype=phenotypeGenomic, select=setselectGenomicObsQ, analyze=analysisObsQ, crossAdvance=crossAdvanceMutUpdateObsQ, processResults=processResults)
}
# Setting 314, 315, 317:
if (settings %in% c(314, 315, 317)){
settingList <- list(initialize=initializeCreateData, phenotype=phenotypeGenomic, select=selectGenomicObsQvarNGP, analyze=analysisObsQ, crossAdvance=crossAdvanceMutUpdateObsQ, processResults=processResults)
}

# Setting 1000:
if (settings %in% c(1000, 1001, 1002, 1003)){
  settingList <- list(initialize=initializeCreateData, phenotype=phenotypePheno, select=selectPheno, analyze=analysisPheno, crossAdvance=crossAdvanceMutPheno, processResults=processResults)
}

# allLists should be as long as the number of settings
allLists <- list(settingList)
