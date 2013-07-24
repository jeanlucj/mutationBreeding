###############################################################
### Setting numbers for this selection program
###############################################################
settings <- c(1004)
allNcycles <- rep(400, length(settings))

settingReminder <- function(){
  cat("Setting 1000 Baseline phenotypic selection, error std dev=8", "\n")
  cat("Setting 1001 Baseline phenotypic selection with free recombination between loci", "\n")
  cat("Setting 1002 Like 1000 but nSelCan=400", "\n")
  cat("Setting 1003 Like 1001 but nSelCan=400", "\n")
  cat("Setting 1004 Baseline genomic selection, seeking optimal nToSelect for GS", "\n")
  cat("Setting 401 initiates improved marker updating algorithm", "\n")
	cat("Setting 402 like 401 but takes expectation of markers into account", "\n")
	cat("Setting 411 try to use mixed integer programming", "\n")
	cat("Setting 421 compare to a slow mixed model analysis for every new update", "\n")
	
	cat("Setting 301 has unobserved QTL with equal numbers of markers", "\n")
	cat("Setting 302 has unobserved QTL with 4 times more markers", "\n")
	cat("Setting 303 has unobserved QTL with 7 times more markers", "\n")
	cat("Setting 305 has observed QTL, with nMrk =~ nQTL", "\n")
	cat("Setting 306 doubles the number of stored generations to 16", "\n")
	cat("Setting 307 like 305 but takes the ancAllele orientation into account", "\n")
	cat("Setting 308 like 307 but makes mutations _favorable_ initially", "\n")
	cat("Setting 309 like 305 but weights rare favorable alleles", "\n")
	cat("Setting 310 like 308 but makes mutations _favorable_ initially, with random boost", "\n")
	cat("Setting 311 combines 309 and 310", "\n")
	cat("Setting 312 stops mutation when you reach 200 cycles then does GS for 20 cycles", "\n")
	cat("Setting 314 like 305 but nGenPoly into account in eMrkVar", "\n")
	cat("Setting 315 like 314 but less severely", "\n")
	cat("Setting 316 like 309 but extent of weighting is random", "\n")
	cat("Setting 317 is a test of the machinery of 314", "\n")
	cat("Setting 318 is a more stringent test of the machinery of 314", "\n")
	cat("Setting 319 explores the prior effect, var, and weight space", "\n")
	cat("Setting 320 Garrick: large effects are overweighted: take sqrt", "\n")
	cat("Setting 321 & 322 Hill: explore different selection intensities", "\n")
	cat("Setting 323: common sense: only phenotype the parents, but at 2x replication")
	cat("Setting 325: first attempt to select the set of individuals to retain mutations")
	cat("Setting 326: second attempt: 325 penalized the loss of mutations or ancestral alleles")
	cat("Setting 327: third attempt: 327 increases penalty if mutation appears favorable")
	
	cat("Setting 206 doubles the number of stored generations to 16", "\n")
	cat("Setting 207 starts marker effects off at their expectations, -0.18", "\n")
	cat("Setting 208 doubles the number of candidates and selected", "\n")
	cat("Setting 209 & 210 assume some ability to decide a mutation is good if it is", "\n")
	cat("Setting 211 doubles the mutation rate", "\n")
	cat("Setting 212 stops mutation when you reach 200 cycles then does GS for 20 cycles", "\n")
	cat("Setting 213 attempts to use penalty.factor to shrink QTL effects less", "\n")
	cat("Setting 214 fits penalty.factor as a function of nGenPoly", "\n")
	cat("Setting 215 combines 214 with weighting of rare favorable alleles", "\n")
	cat("Setting 216 fits penalty.factor as a function of nGenPoly for _ridge_ penalty", "\n")
	cat("Setting 217 combines 216 with weighting of rare favorable alleles", "\n")
	cat("Setting 218 uses constant penalty but weights rare favorable alleles", "\n")
	cat("Setting 219 updates mrk eff by cycle starting at 0", "\n")
	cat("Setting 220 updates mrk eff by cycle starting at 10% correlation", "\n")
	cat("Setting 221 updates mrk eff by cycle starting at 50% correlation", "\n")
	
	cat("The following are all GBLUP based and were done with fixed loci (as opposed to locInOut)", "\n")
	cat("6: all loci are QTL; nSelCan=100, nToSelect=20; TP=8 generations", "\n")
	cat("7: like 6 but 50% of loci are QTL", "\n")
	cat("8: like 6 but 25% of loci are QTL", "\n")
	cat("9: like 6 but nSelCan=200, nToSelect=40", "\n")
	cat("10: like 9 but 50% of loci are QTL", "\n")
	cat("11: like 10 but rare favorable alleles are upweighted", "\n")
	cat("12: like 10 but TP=4 generations", "\n")
	cat("13: like 10 but separate relationship matrices are calculated for standing versus mutational variation", "\n")
	cat("14: like 13 but in cycles where the model is updated, selects from best of last cyc with pheno and next cyc without", "\n")
	cat("15: like 13 but in cycles where the model is updated, selects 30 from best of next cyc and 10 from last cyc with pheno that have high mutational BLUPs", "\n")
	cat("16-18: like 13-15 but nSelCan=400, nToSelect=80", "\n")
	cat("19: like 13 but rare favorable alleles are upweighted", "\n")
	cat("20: like 19 but nSelCan=400, nToSelect=80", "\n")
}

###############################################################
#					Initialize
# Function to create or load speciesData or breedingData
###############################################################
initializeCreateData <- function(replication){
	switchCyc <<- nCycles / 2
	if (setting == 312) switchCyc <<- 200
	if (makeFounderFiles){
		mutNum <- 2
		ratioLocToQTL <- 2
		if (setting == 302){
			mutNum <- 5
			ratioLocToQTL <- 5
		}
		if (setting == 303){
			mutNum <- 8
			ratioLocToQTL <- 8
		}
		nSelCan <- 200
		nToSelect <- 40
		effPopSize <- nToSelect
		
		speciesData <<- createSpeciesData(effPopSize, mutNum, ratioLocToQTL, founderSeed)
	
		if (setting %in% c(1001, 1003)) makeGamete <<- function(genoVec, map, chrMax){
		  nLoc <- nrow(map)
		  gamete <- rbinom(nLoc, 1, 0.5)
		  newGamete <- genoVec[1:nLoc * 2 - gamete]
		  return(newGamete)
		}
    
		nStoredGen <- 8
		if (setting == 306) nStoredGen <- 16
		if (setting %in% c(321, 322, 1000, 1001)) nToSelect <- floor(runif(1, 20, 121))
		if (setting %in% c(1002, 1003)){
      nToSelect <- floor(runif(1, 60, 221))
      nSelCan <- 400
		}
		if (setting %in% c(323, 325, 326, 327, 401, 402, 411, 421)) nToSelect <- 70
		breedingData <<- createBreedingData(nSelCan, nToSelect, nStoredGen=nStoredGen, founderSeed)
		breedingData$postPhenSelMrkEffEst <<- 0 # These will need to be tailored to the setting
		breedingData$postPhenSelMrkEffVar <<- 0.366
		if (setting %in% c(323, 325, 326, 327, 401, 402, 411, 421, 1004)) breedingData$nToSelectGS <<- c(80, 100)
		if (setting %in% c(310:311, 319)) breedingData$priorMutEff <<- runif(1, -0.2, 0.1)
		if (setting %in% c(319)) breedingData$priorMutVar <<- runif(1, 0.28, 1.4)
		if (setting %in% c(316, 319)) breedingData$rareFavWgt <<- runif(1)
		if (setting %in% c(320)) breedingData$effExponent <<- runif(1, 0.1, 0.9)
		if (setting %in% c(322, 1004)){
			afterUpD <- floor(runif(1, 50, 111))
			afterThat <- floor(runif(1, 80, 141))
			if (setting == 1004) breedingData$nToSelect <<- 85
			breedingData$nToSelectGS <<- c(afterUpD, afterThat)
			cat("Setting 322, 1004 ", nToSelect, afterUpD, afterThat, "\n")
		}
		breedingData$cumulativeIncidences <<- numeric(nrow(speciesData$map))
		if (setting %in% c(325, 326, 327, 411)){
			breedingData$penalty <<- replication %% 5 * 0.1 + runif(1, 0, 0.1)
			breedingData$penalty <<- replication %% 8 * 0.1 + runif(1, 0.1, 0.2) # for 411
			breedingData$maxMatters <<- round(exp(2 + replication %/% 8 %% 3 + runif(1)))
		}
		save(replication, founderSeed, founderFile, speciesData, breedingData, file=founderFile)
	} else{
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
	genRows <- which(breedingData$records$cycle == cycle)
	qtl <- speciesData$genArch$locusList
	qtlDosage <- (breedingData$genoMat[genRows, qtl*2 - 1] + breedingData$genoMat[genRows, qtl*2]) / 2
	breedingData$records$genoVal[genRows] <- speciesData$genArch$genoBase + c(qtlDosage %*% speciesData$genArch$effects)
	
	if (cycle == switchCyc | (cycle > switchCyc & cycle %% 2 == 1)){
		# No phenotyping at the transition or alternate cycles for GS
		phenRows <- 0
	} else{ # Here we are phenotyping
		if (cycle < switchCyc){
			phenRows <- which(breedingData$records$cycle == cycle)
			tmp <- calcPhenotype(breedingData$records$genoVal[phenRows], breedingData$stdDevErr, heritability=FALSE)
		} else{
			if (setting %in% c(323, 325, 326, 327, 401, 402, 411, 421, 1004)){ # Phenotype only latest parents with lower error
				phenRows <- breedingData$selectSet
				tmp <- calcPhenotype(breedingData$records$genoVal[phenRows], breedingData$stdDevErr / sqrt(2), heritability=FALSE)
			} else{ # Phenotype all prior selection candidates
				phenRows <- which(breedingData$records$cycle %in% (cycle - (1:2)))
				tmp <- calcPhenotype(breedingData$records$genoVal[phenRows], breedingData$stdDevErr, heritability=FALSE)
			}
		}
		eo <- 1:nrow(speciesData$map) * 2
		inst <- -breedingData$ancAllele * colSums(breedingData$genoMat[phenRows, eo-1] + breedingData$genoMat[phenRows, eo]) / 2 + length(phenRows)
		breedingData$cumulativeIncidences <- breedingData$cumulativeIncidences + inst
		breedingData$records$phenoVal[phenRows] <- tmp$phenoVal
	}
	cat(cycle, "phenotypeGenomic", range(phenRows), "\n")
	return(breedingData)
}#END phenotypeGenomic

# Phenotyping for genomic selection.
phenotypePheno <- function(cycle, breedingData){
  phenRows <- which(breedingData$records$cycle == cycle)
  qtl <- speciesData$genArch$locusList
  qtlDosage <- (breedingData$genoMat[phenRows, qtl*2 - 1] + breedingData$genoMat[phenRows, qtl*2]) / 2
  breedingData$records$genoVal[phenRows] <- speciesData$genArch$genoBase + c(qtlDosage %*% speciesData$genArch$effects)
  tmp <- calcPhenotype(breedingData$records$genoVal[phenRows], breedingData$stdDevErr, heritability=FALSE)
  breedingData$records$phenoVal[phenRows] <- tmp$phenoVal
  
  eo <- 1:nrow(speciesData$map) * 2
  inst <- -breedingData$ancAllele * colSums(breedingData$genoMat[phenRows, eo-1] + breedingData$genoMat[phenRows, eo]) / 2 + length(phenRows)
  breedingData$cumulativeIncidences <- breedingData$cumulativeIncidences + inst
  
  cat(cycle, "phenotypePheno", range(phenRows), "\n")
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
# WARNING: assumes ratioLocToQTL = 2 to calculate eMrkVar
# WARNING: I am dividing mrkDosage by 2 to fit the A.mat function
###############################################################
# On even cycles update the model. On odd cycles, use the marker effect estimates.
selectGenomicUnobsQ <- function(cycle, breedingData){
	thisCycle <- which(breedingData$records$cycle == cycle)
	if (cycle < switchCyc){ # Phenotypic selection
		cat(cycle, "selectPhenotypic_UnobsQupdME", "\n")
		selInd <- order(breedingData$records$phenoVal[thisCycle], decreasing=TRUE)
	} else{
		# I think the expected prior marker effect should be zero
		# This might be wrong: I will consider that we know nothing about the
		# marker effects, hence the prior variance will be the phenotypic variance
		mrk <- (1:nrow(speciesData$map))[-speciesData$genArch$locusList]
		nMrk <- length(mrk)
		eo <- mrk * 2
		if (cycle == switchCyc){ # Create prior marker effect estimates
			newPheno <- switchCyc - breedingData$nStoredGen + 1:(breedingData$nStoredGen - 1)
			newPheno <- breedingData$records$cycle %in% newPheno
			breedingData$mrkEffEst <- rep(0, nMrk)
			breedingData$mrkEffVar <- rep(var(breedingData$records$phenoVal[newPheno]), nMrk)
		}
		if (cycle %% 2 == 0){ # Calculate effect GEBVs with new phenotypes
			cat(cycle, "selectGenomicRR_ModelUpdate_UnobsQupdME", "\n")
			if (cycle > switchCyc) newPheno <- breedingData$records$cycle %in% (cycle - (1:2))
			mrkDosage <- breedingData$genoMat[newPheno,eo - 1] + breedingData$genoMat[newPheno,eo] # Just marker dosage not all loci
			# No point wasting breath for monomorphic markers
			polymorphic <- which(apply(mrkDosage, 2, sd) > 0)
			mrkDosage <- mrkDosage[,polymorphic]
			# Expectation and var of marker effects conditional on phenotypes Ober 2012 PLoS Gen
			relMat <- A.mat(mrkDosage, shrink=FALSE)
			mixedSolveOut <- mixed.solve(y=breedingData$records$phenoVal[newPheno], K=relMat)
			mrkFreq <- colMeans(mrkDosage) / 4 + 0.5
			eMrkVar <- mixedSolveOut$Vu / sum(mrkFreq * (1 - mrkFreq)) / 2
			mrkCenter <- scale(mrkDosage, scale=FALSE)
			varInv <- solve(mixedSolveOut$Vu * relMat + diag(mixedSolveOut$Ve, sum(newPheno)))
			regCoef <- eMrkVar * crossprod(mrkCenter, varInv)
			newEffEst <- regCoef %*% scale(breedingData$records$phenoVal[newPheno], scale=FALSE)
			newEffVar <- eMrkVar * (1 - diag(regCoef %*% mrkCenter))
			# Combine previous estimate with current estimate: wgt by inverse of variance
			newVar <- 1 / (1 / newEffVar + 1 / breedingData$mrkEffVar[polymorphic])
			newEst <- (newEffEst / newEffVar + breedingData$mrkEffEst[polymorphic] / breedingData$mrkEffVar[polymorphic]) * newVar
			breedingData$mrkEffEst[polymorphic] <- newEst
			breedingData$mrkEffVar[polymorphic] <- newVar
		} else{ # Calculate effect GEBVs with previously estimated effects
			cat(cycle, "selectGenomicRR_PreviousModel_ObsQupdME", "\n")
		}
		mrkDosage <- breedingData$genoMat[thisCycle,eo - 1] + breedingData$genoMat[thisCycle,eo]
		breedingData$records$genoHat[thisCycle] <- mrkDosage %*% breedingData$mrkEffEst
		selInd <- order(breedingData$records$genoHat[thisCycle], decreasing=TRUE)
	}
	breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelect]]
	breedingData$allSelectSet <- c(breedingData$allSelectSet, breedingData$selectSet)
	breedingData$selGenoVal <- c(mean(breedingData$records$genoVal[breedingData$selectSet]), NA)
	return(breedingData)
}#END selectGenomicUnobsQ

# 11 April 2013.  I _think_ I have figured out a correct function for updating
# the marker effects and variances.  I am retaining above the previous versions for historical
# reasons and in case I'm wrong.
selectGenomicObsQ <- function(cycle, breedingData){
	thisCycle <- which(breedingData$records$cycle == cycle)
	if (cycle < switchCyc){ # Phenotypic selection
		cat(cycle, "selectPhenotypic_ObsQ", "\n")
		selInd <- order(breedingData$records$phenoVal[thisCycle], decreasing=TRUE)
		breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelect]]
	} else{
		nLoc <- nrow(speciesData$map)
		eo <- 1:nLoc * 2
		if (cycle %% 2 == 0){ # Calculate effect GEBVs with new phenotypes
			if (cycle == switchCyc){ # Create prior marker effect estimates coming out of phenotypic selection
			  # Simulations in the baseline setting estimate that marker effects and variance thereof are:
			  # *Unconditional on the identity of the ancestral allele: mrkEffEst=0; mrkEffVar=0.366
			  # *Conditional on the identity of the ancestral allele: mrkEffEst=0.224 (ancestral is favorable); mrkEffVar=0.316
        # Those values need to be set in postPhenSelMrkEffEst and postPhenSelMrkEffVar
        # Currently set up for the unconditional version
			  newPheno <- which(!is.na(breedingData$records$phenoVal))
				cat(cycle, "selectGenomicRR_FirstModelRR", length(newPheno), "\n")
				errVar <- breedingData$stdDevErr^2
				breedingData$mrkEffEst <- breedingData$postPhenSelMrkEffEst * breedingData$ancAllele
				breedingData$mrkEffVar <- rep(breedingData$postPhenSelMrkEffVar, nLoc)
			  mrkDosage <- (breedingData$genoMat[newPheno,eo - 1] + breedingData$genoMat[newPheno,eo]) / 2
			  breedingData$meanPhen <- mean(breedingData$records$phenoVal[newPheno] - mrkDosage %*% breedingData$mrkEffEst)
			} else{
				errVar <- breedingData$stdDevErr^2 / 2 #WARNING! This assumes two-rep phenotyping
				newPheno <- cycle - (1:2)
				newPheno <- which(breedingData$records$cycle %in% newPheno)
				newPheno <- newPheno[!is.na(breedingData$records$phenoVal[newPheno])]
				cat(cycle, "selectGenomicRR_ModelUpdate_ObsQ", length(newPheno), "\n")
				mrkDosage <- (breedingData$genoMat[newPheno,eo - 1] + breedingData$genoMat[newPheno,eo]) / 2
			}
			mrkEffEst <- breedingData$mrkEffEst
			mrkEffVar <- breedingData$mrkEffVar
			covPhenMrk <- sapply(1:ncol(mrkDosage), FUN=function(col) mrkDosage[,col] * mrkEffVar[col]) # t(cov(\beta, y))
			varPhen <- tcrossprod(covPhenMrk, mrkDosage) + diag(errVar, length(newPheno))
			regCoef <- t(solve(varPhen, covPhenMrk))
			breedingData$mrkEffEst <- mrkEffEst + regCoef %*% (breedingData$records$phenoVal[newPheno] - breedingData$meanPhen - mrkDosage %*% mrkEffEst)
			breedingData$mrkEffVar <- mrkEffVar - diag(regCoef %*% covPhenMrk)
		} else{ # Calculate effect GEBVs with previously estimated effects
			cat(cycle, "selectGenomicRR_PreviousModel_ObsQ", "\n")
		}
		mrkDosage <- (breedingData$genoMat[thisCycle,eo - 1] + breedingData$genoMat[thisCycle,eo]) / 2
		breedingData$records$genoHat[thisCycle] <- mrkDosage %*% breedingData$mrkEffEst
		selInd <- order(breedingData$records$genoHat[thisCycle], decreasing=TRUE)
		breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelectGS[cycle %% 2 + 1]]]
	}#END doing GS
	breedingData$allSelectSet <- c(breedingData$allSelectSet, breedingData$selectSet)
	breedingData$selGenoVal <- c(mean(breedingData$records$genoVal[breedingData$selectSet]), NA)
	return(breedingData)
}#END selectGenomicObsQ

# Select based on phenotype
selectPheno <- function(cycle, breedingData){
  thisCycle <- which(breedingData$records$cycle == cycle)
  cat(cycle, "selectPheno", "\n")
  selInd <- order(breedingData$records$phenoVal[thisCycle], decreasing=TRUE)
  breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelect]]
  breedingData$allSelectSet <- c(breedingData$allSelectSet, breedingData$selectSet)
  breedingData$selGenoVal <- mean(breedingData$records$genoVal[breedingData$selectSet])
  return(breedingData)
}#END selectPheno

# For setting 421 where I want to use straight up rr to compare with the updating method
selectGenomicObsQrr <- function(cycle, breedingData){
	thisCycle <- which(breedingData$records$cycle == cycle)
	if (cycle < switchCyc){ # Phenotypic selection
		cat(cycle, "selectPhenotypic_ObsQ", "\n")
		selInd <- order(breedingData$records$phenoVal[thisCycle], decreasing=TRUE)
	} else{
		eo <- 1:nrow(speciesData$map) * 2
		if (cycle %% 2 == 0){ # Calculate effect GEBVs with new phenotypes
			hasPheno <- which(!is.na(breedingData$records$phenoVal))
			cat(cycle, "selectGenomicRR_ModelUpdate_ObsQrr", length(hasPheno), "\n")
			inRelMat <- c(thisCycle, hasPheno)
			mrkDosage <- breedingData$genoMat[inRelMat,eo - 1] + breedingData$genoMat[inRelMat,eo]
			# No point wasting breath for monomorphic markers
			polymorphic <- which(apply(mrkDosage, 2, sd) > 0)
			mrkDosage <- mrkDosage[,polymorphic]
			relMat <- A.mat(mrkDosage, shrink=FALSE)
			mixedSolveOut <- mixed.solve(y=breedingData$records$phenoVal[inRelMat], K=relMat)
			breedingData$records$genoHat[thisCycle] <- mixedSolveOut$u[1:breedingData$nSelCan]

			# Calculate relative marker effects for other uses from Ober 2012 PLoS Genetics
			mrkInc <- t(scale(mrkDosage[-(1:breedingData$nSelCan),], scale=FALSE))
			coefMat <- solve(mixedSolveOut$Vu * relMat[-(1:breedingData$nSelCan), -(1:breedingData$nSelCan)] + diag(mixedSolveOut$Ve, length(hasPheno)))
			yAdj <- breedingData$records$phenoVal[hasPheno] - mixedSolveOut$beta
			mrkFreq <- colMeans(mrkDosage) / 4 + 0.5
			constProp <- mixedSolveOut$Vu / sum(mrkFreq * (1 - mrkFreq)) / 2
			breedingData$mrkEffEst[polymorphic] <- constProp * mrkInc %*% coefMat %*% yAdj
		} else{ # Calculate effect GEBVs with previously estimated effects
			cat(cycle, "selectGenomicRR_PreviousModel_ObsQrr", "\n")
			mrkDosage <- breedingData$genoMat[thisCycle,eo - 1] + breedingData$genoMat[thisCycle,eo]
			breedingData$records$genoHat[thisCycle] <- mrkDosage %*% breedingData$mrkEffEst
		}
		selInd <- order(breedingData$records$genoHat[thisCycle], decreasing=TRUE)
	}
	if (cycle < switchCyc){
		breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelect]]
	} else{
		breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelectGS[cycle %% 2 + 1]]]
	}
	breedingData$allSelectSet <- c(breedingData$allSelectSet, breedingData$selectSet)
	breedingData$selGenoVal <- c(mean(breedingData$records$genoVal[breedingData$selectSet]), NA)
	return(breedingData)
}#END selectGenomicObsQTLrr

########################################################################################
# For now, this is just for setting 411
setSelectGenomicObsQ <- function(cycle, breedingData){
	thisCycle <- which(breedingData$records$cycle == cycle)
	if (cycle < switchCyc){ # Phenotypic selection
		cat(cycle, "selectPhenotypic_ObsQ", "\n")
		selInd <- order(breedingData$records$phenoVal[thisCycle], decreasing=TRUE)
		breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelect]]
	} else{
		nLoc <- nrow(speciesData$map)
		eo <- 1:nLoc * 2
		if (cycle %% 2 == 0){ # Calculate effect GEBVs with new phenotypes
			if (cycle == switchCyc){ # Create prior marker effect estimates coming out of phenotypic selection
				hasPheno <- which(!is.na(breedingData$records$phenoVal))
				cat(cycle, "selectGenomicRR_FirstModelRR", length(hasPheno), "\n")
				mrkDosage <- breedingData$genoMat[hasPheno,eo - 1] + breedingData$genoMat[hasPheno,eo]
				relMat <- A.mat(mrkDosage / 2, shrink=FALSE)
				mixedSolveOut <- mixed.solve(y=breedingData$records$phenoVal[hasPheno], K=relMat)
			
				# Calculate relative marker effects for other uses from Ober 2012 PLoS Genetics
				mrkFreq <- colMeans(mrkDosage) / 4 + 0.5
				mrkEffVar <- mixedSolveOut$Vu / sum(mrkFreq * (1 - mrkFreq)) / 2
				covPhenMrk <- mrkEffVar * scale(mrkDosage, scale=FALSE)
				regCoef <- t(solve(mixedSolveOut$Vu * relMat + diag(mixedSolveOut$Ve, length(hasPheno)), covPhenMrk))
				yAdj <- breedingData$records$phenoVal[hasPheno] - mixedSolveOut$beta
				breedingData$mrkEffEst <- regCoef %*% yAdj
				breedingData$mrkEffVar <- rep(mrkEffVar, nLoc)
			} else{
				cat(cycle, "selectGenomicRR_ModelUpdate_ObsQ", "\n")
				errVar <- breedingData$stdDevErr^2 / 2 #Warning! This assumes two-rep phenotyping
				newPheno <- cycle - (1:2)
				newPheno <- which(breedingData$records$cycle %in% newPheno)
				newPheno <- newPheno[!is.na(breedingData$records$phenoVal[newPheno])]
				mrkDosage <- breedingData$genoMat[newPheno,eo - 1] + breedingData$genoMat[newPheno,eo]
				# No point wasting breath for monomorphic markers
				polymorphic <- which(apply(mrkDosage, 2, sd) > 0)
				mrkDosage <- mrkDosage[,polymorphic]
				# Expectation and var of marker effects conditional on phenotypes Ober 2012 PLoS Gen
				mrkCenter <- scale(mrkDosage, scale=FALSE)
				mrkEffVarPoly <- breedingData$mrkEffVar[polymorphic]
				mrkTimesVar <- function(col) mrkCenter[,col] * mrkEffVarPoly[col]
				covPhenMrk <- sapply(1:ncol(mrkCenter), mrkTimesVar) # \Sigma_21
				regCoef <- t(solve(tcrossprod(covPhenMrk, mrkCenter) + diag(errVar, length(newPheno)), covPhenMrk))
				yAdj <- breedingData$records$phenoVal[newPheno]
				if (setting %in% c(402)) yAdj <- yAdj - mrkDosage %*% breedingData$mrkEffEst[polymorphic]
				breedingData$mrkEffEst[polymorphic] <- breedingData$mrkEffEst[polymorphic] + regCoef %*% scale(yAdj, scale=FALSE)
				breedingData$mrkEffVar[polymorphic] <- breedingData$mrkEffVar[polymorphic] - diag(regCoef %*% covPhenMrk)
			}
		} else{ # Calculate effect GEBVs with previously estimated effects
			cat(cycle, "selectGenomicRR_PreviousModel_ObsQ", "\n")
		}
		mrkDosage <- breedingData$genoMat[thisCycle,eo - 1] + breedingData$genoMat[thisCycle,eo]
		breedingData$records$genoHat[thisCycle] <- mrkDosage %*% breedingData$mrkEffEst

		# ancestral: -1 => mutation is +1 => positive effect is good => put the estimated effect in pnorm, then you want 1 - pnorm as the penalty: To simplify this, put ancestral * estimate in pnorm
		# ancestral: +1 => mutation is -1 => negative effect is good => put the estimated effect in pnorm, then you want pnorm as the penalty
		mutLossPenalty <- pnorm(0, breedingData$ancAllele * breedingData$mrkEffEst, sqrt(breedingData$mrkEffVar))
		maxMatters <- 1 - breedingData$cumulativeIncidences / breedingData$maxMatters
		maxMatters <- ifelse(maxMatters < 0, 0, maxMatters)
		calcMrkPenal <- function(col) (-breedingData$ancAllele[col] * mrkDosage[,col]/2 + 1) * mutLossPenalty[col] * maxMatters[col]
		penCoefVec <- rowSums(sapply(1:nLoc, calcMrkPenal))
		nToSelect <- breedingData$nToSelectGS[cycle %% 2 + 1]
		minPen <- sum(penCoefVec[order(penCoefVec, decreasing=FALSE)[1:nToSelect]])
		maxPen <- sum(penCoefVec[order(penCoefVec, decreasing=TRUE)[1:nToSelect]])
		obj <- breedingData$records$genoHat[thisCycle]
		mat <- rbind(1, penCoefVec)
		dir <- c(rep("<=", 2))
		rhs <- c(nToSelect, minPen + breedingData$penalty * (maxPen - minPen))
		save(breedingData, thisCycle, obj, mat, dir, rhs, file="./GSLocalBulmer/testLinProg.RData")
		selectSet <- thisCycle[Rglpk_solve_LP(obj, mat, dir, rhs, types="B", max=TRUE)$solution==1]	
		breedingData$selectSet <- selectSet
	}#END doing GS
	breedingData$allSelectSet <- c(breedingData$allSelectSet, breedingData$selectSet)
	breedingData$selGenoVal <- mean(breedingData$records$genoVal[breedingData$selectSet])
	return(breedingData)
}#END setSelectGenomicObsQ

# Model the marker variance as a function of the number of generations it has been polymorphic
# For settings 314 and 315 for the moment
selectGenomicObsQvarNGP <- function(cycle, breedingData){
	thisCycle <- which(breedingData$records$cycle == cycle)
	if (cycle < switchCyc){ # Phenotypic selection
		cat(cycle, "selectPhenotypic_ObsQupdME", "\n")
		selInd <- order(breedingData$records$phenoVal[thisCycle], decreasing=TRUE)
	} else{
		nLoc <- nrow(speciesData$map)
		eo <- 1:nLoc * 2
		# eMrkVar is going to be a _vector_ based on empirical observation from simulation
		# This is all very adhoc and should be guided by theory
		eMrkVar <- numeric(nLoc)
		nGP <- breedingData$nGenPoly
		nGP37 <- nGP[nGP < 37]
		eMrkVar[nGP < 37] <- cbind(1, nGP37, nGP37^2) %*% c(3.436, 0.0115, 0.04896)
		nGP37 <- nGP[nGP > 36]
		eMrkVar[nGP > 36] <- cbind(1, nGP37) %*% c(-43.733, 3.0288)
		eMrkVar <- 1 / eMrkVar
		if (setting == 315) eMrkVar <- (eMrkVar + 0.28) / 2 # Hedge in case too extreme
		if (setting == 317){
			qtl <- speciesData$genArch$locusList
			eMrkVar[qtl] <- 0.28
			eMrkVar[-qtl] <- 0.07
		}
		if (setting == 318){
			qtl <- speciesData$genArch$locusList
			eMrkVar[qtl] <- 0.28
			eMrkVar[-qtl] <- 0.01
		}
		if (cycle == switchCyc){ # Create prior marker effect estimates
			newPheno <- switchCyc - breedingData$nStoredGen + 1:(breedingData$nStoredGen - 1)
			breedingData$mrkEffVar <- eMrkVar
			breedingData$mrkEffEst <- rep(0, nLoc)
		}
		if (cycle %% 2 == 0){ # Calculate effect GEBVs with new phenotypes
			cat(cycle, "selectGenomicRR_ModelUpdate_ObsQnGP", "\n")
			if (cycle > switchCyc) newPheno <- cycle - (1:2)
			newPheno <- which(breedingData$records$cycle %in% newPheno)
			mrkDosage <- breedingData$genoMat[newPheno,eo - 1] + breedingData$genoMat[newPheno,eo]
			# No point wasting breath for monomorphic markers
			polymorphic <- which(apply(mrkDosage, 2, sd) > 0)
			mrkDosage <- mrkDosage[,polymorphic]
			# Expectation and var of marker effects conditional on phenotypes Ober 2012 PLoS Gen
			diagEMV <- diag(eMrkVar[polymorphic])
			mrkCenter <- scale(mrkDosage, scale=FALSE)
			covMrkY <- mrkCenter %*% diagEMV
			varYinv <- solve(tcrossprod(covMrkY,  mrkCenter) + diag(breedingData$stdDevErr^2, length(newPheno)))
			regCoef <- crossprod(covMrkY, varYinv)
			newEffEst <- regCoef %*% scale(breedingData$records$phenoVal[newPheno], scale=FALSE)
			newEffVar <- diag(diagEMV - regCoef %*% covMrkY)
			# Combine previous estimate with current estimate: wgt by inverse of variance
			newVar <- 1 / (1 / newEffVar + 1 / breedingData$mrkEffVar[polymorphic])
			newEst <- (breedingData$mrkEffEst[polymorphic] / breedingData$mrkEffVar[polymorphic] + newEffEst / newEffVar) * newVar
			breedingData$mrkEffEst[polymorphic] <- newEst
			breedingData$mrkEffVar[polymorphic] <- newVar
		} else{ # Calculate effect GEBVs with previously estimated effects
			cat(cycle, "selectGenomicRR_PreviousModel_ObsQnGP", "\n")
		}
		mrkDosage <- breedingData$genoMat[thisCycle,eo - 1] + breedingData$genoMat[thisCycle,eo]
		breedingData$records$genoHat[thisCycle] <- mrkDosage %*% breedingData$mrkEffEst
		selInd <- order(breedingData$records$genoHat[thisCycle], decreasing=TRUE)
	}
	breedingData$selectSet <- thisCycle[selInd[1:breedingData$nToSelect]]
	breedingData$allSelectSet <- c(breedingData$allSelectSet, breedingData$selectSet)
	breedingData$selGenoVal <- c(mean(breedingData$records$genoVal[breedingData$selectSet]), NA)
	return(breedingData)
}#END selectGenomicObsQvarNGP

###############################################################
#					ANALYZE
# This needs to go after selection but before cross & advance
# Analyze the results of one cycle
###############################################################
analysisUnobsQ <- function(cycle, breedingData){
	cat(cycle, "analysisLocalStats", "\n")
	thisCycle <- which(breedingData$records$cycle == cycle)
	
	thin <- nCycles / 40

	# Calculate equiVar; segmentVar; chromosomeVar; genoVar
	nLoc <- nrow(speciesData$map)
	eo <- 1:nLoc*2
	locDosage <- breedingData$genoMat[thisCycle, eo - 1] + breedingData$genoMat[thisCycle, eo]
	qtl <- speciesData$genArch$locusList
	mrk <- (1:nLoc)[-qtl]
	qtlDosage <- locDosage[,qtl]
	qtlDosageVar <- apply(qtlDosage, 2, var)
	totEquiVar <- c(crossprod(speciesData$genArch$effects^2, qtlDosageVar))
	mrkDosage <- locDosage[,mrk]
	
	# Calculate the variance of a segment (could be a whole chromosome)
	# NOTE: assumes that the segment does NOT overlap two chromosomes
	# Assume the simplest additive gene action for now
	# Only do this once equilibrium reached (see below)
	# WARNING depends on calculating equiVar before hand
	calcSegmentVar <- function(startPos, segSize, calcEquiVar=FALSE, calcSegAcc=FALSE, calcPoly=FALSE){
		segGenoVals <- NA
		if (calcEquiVar) equiVar <- NA else equiVar <- NULL
		if (calcSegAcc) segAcc <- NA else segAcc <- NULL
		if (calcPoly){
			segQTLpoly <- NA; segMrkPoly <- NA
		} else{
			segQTLpoly <- NULL; segMrkPoly <- NULL
		}
		segLoci <- which(speciesData$map[,3] >= startPos & speciesData$map[,3] < startPos + segSize) # segLoci in locus order
		if (length(segLoci > 0)){
			segQTL <- which(qtl %in% segLoci) # segQTL in qtl order
			segGenoVals <- qtlDosage[,segQTL, drop=FALSE] %*% speciesData$genArch$effects[segQTL]
			if (calcEquiVar) equiVar <- crossprod(speciesData$genArch$effects[segQTL]^2, qtlDosageVar[segQTL])
			if (calcSegAcc){
				segMrk <- which(mrk %in% segLoci)
				segGenoHat <- mrkDosage[,segMrk, drop=FALSE] %*% breedingData$mrkEffEst[segMrk]
				segAcc <- cor(segGenoVals, segGenoHat)
			}
			if (calcPoly){ # if all ind are heterozygous, it will look monomorphic
				segQTLpoly <- sum(apply(qtlDosage[,segQTL, drop=FALSE], 2, sd) > 0) / length(segQTL)
				segMrkPoly <- sum(apply(mrkDosage[,segMrk, drop=FALSE], 2, sd) > 0) / length(segMrk)
			}
		}
		return(c(var(segGenoVals), equiVar, segAcc, segQTLpoly, segMrkPoly))
	}
	
	# Calculate segment variances every thin cycles on the whole genome
	chrVar <- NA
	segVar10 <- matrix(NA, 5, 105)
	if (cycle >= switchCyc & cycle %% thin == 0){
		chrVar <- sapply(0:6 * 150, calcSegmentVar, segSize=150) # WARNING chr size of 150 hardwired
		seg10 <- rep(0:6 * 150, each=15) + 0:14 * 10
		segVar10 <- sapply(seg10, calcSegmentVar, segSize=10, calcEquiVar=TRUE, calcSegAcc=TRUE, calcPoly=TRUE)
		breedingData$localStats <- c(breedingData$localStats, list(chrVar, segVar10))
	}
		
	# Collect some information on this cycle
	genoMean <- mean(breedingData$records$genoVal[thisCycle])
	genoVar <- var(breedingData$records$genoVal[thisCycle])
	thisCycleAccuracy <- cor(breedingData$records$genoVal[thisCycle], breedingData$records$genoHat[thisCycle])
	breedingData$history <- rbind(breedingData$history, data.frame(cycle=cycle, genoBase=speciesData$genArch$genoBase, nQTLpoly=sum(apply(qtlDosage, 2, sd) > 0), nQTL=length(qtl), nMrkPoly=sum(apply(mrkDosage, 2, sd) > 0), nMrk=nLoc - length(qtl), meanGenoVal=genoMean, accuracy=thisCycleAccuracy, genoVar=genoVar, totChrVar=sum(chrVar, na.rm=TRUE), totSegVar=sum(segVar10[1,], na.rm=TRUE), equiVar=totEquiVar, selGenoValMain=breedingData$selGenoVal[1], selGenoValSpec=breedingData$selGenoVal[2]))
	
	# If there is mutation, recalibrate until equilibrium
	if (!(setting %in% c(321, 322, 323, 325, 326, 327, 401, 402, 411, 421))){
	if (exists("mutParms", where=speciesData)){
		reEquil <- switchCyc / 10
		if (cycle <= switchCyc & cycle %% reEquil == 0){
			meanGenoVar <- mean(breedingData$history$genoVar[(cycle - reEquil)+1:reEquil], na.rm=TRUE)
			h2 <- breedingData$heritability
			breedingData$stdDevErr <- sqrt(meanGenoVar * (1 - h2) / h2)
			cat(cycle, "Recalc stdDevErr:", breedingData$stdDevErr, "\n")
		}
	}
}
	return(breedingData)
}#END analysisUnobsQ

analysisObsQ <- function(cycle, breedingData){
	cat(cycle, "analysisLocalStats_ObsQ", "\n")
	thisCycle <- which(breedingData$records$cycle == cycle)
	thin <- nCycles / 40
	# Calculate equiVar; segmentVar; chromosomeVar; genoVar
	nLoc <- nrow(speciesData$map)
	eo <- 1:nLoc*2
	locDosage <- breedingData$genoMat[thisCycle, eo - 1] + breedingData$genoMat[thisCycle, eo]
	qtl <- speciesData$genArch$locusList
	qtlDosage <- locDosage[,qtl]
	qtlDosageVar <- apply(qtlDosage, 2, var)
	totEquiVar <- c(crossprod(speciesData$genArch$effects^2, qtlDosageVar))
	
	# Calculate the variance of a segment (could be a whole chromosome)
	# NOTE: assumes that the segment does NOT overlap two chromosomes
	# Assume the simplest additive gene action for now
	# Only do this once equilibrium reached (see below)
	# WARNING depends on calculating equiVar before hand
	calcSegmentVar <- function(startPos, segSize, calcEquiVar=FALSE, calcSegAcc=FALSE, calcPoly=FALSE){
		segGenoVals <- NA
		if (calcEquiVar) equiVar <- NA else equiVar <- NULL
		if (calcSegAcc) segAcc <- NA else segAcc <- NULL
		if (calcPoly) segLocPoly <- NA else segLocPoly <- NULL
		segLoci <- which(speciesData$map[,3] >= startPos & speciesData$map[,3] < startPos + segSize) # segLoci in locus order
		if (length(segLoci > 0)){
			segQTL <- which(qtl %in% segLoci) # segQTL in qtl order
			segGenoVals <- qtlDosage[,segQTL, drop=FALSE] %*% speciesData$genArch$effects[segQTL]
			if (calcEquiVar) equiVar <- crossprod(speciesData$genArch$effects[segQTL]^2, qtlDosageVar[segQTL])
			if (calcSegAcc){
				segGenoHat <- locDosage[, segLoci, drop=FALSE] %*% breedingData$mrkEffEst[segLoci]
				segAcc <- cor(segGenoVals, segGenoHat)
			}
			if (calcPoly){ # if all ind are heterozygous, it will look monomorphic
				segLocPoly <- sum(apply(locDosage[,segLoci, drop=FALSE], 2, sd) > 0) / length(segLoci)
			}
		}
		return(c(var(segGenoVals), equiVar, segAcc, segLocPoly))
	}
	
	# Calculate segment variances every thin cycles on the whole genome
	chrVar <- NA
	segVar10 <- matrix(NA, 4, 105)
	if (cycle >= switchCyc & cycle %% thin == 0){
		chrVar <- sapply(0:6 * 150, calcSegmentVar, segSize=150) # WARNING chr size of 150 hardwired
		seg10 <- rep(0:6 * 150, each=15) + 0:14 * 10
		segVar10 <- sapply(seg10, calcSegmentVar, segSize=10, calcEquiVar=TRUE, calcSegAcc=TRUE, calcPoly=TRUE)
		breedingData$localStats <- c(breedingData$localStats, list(chrVar, segVar10))
	}
		
	# Collect some information on this cycle
	genoMean <- mean(breedingData$records$genoVal[thisCycle])
	genoVar <- var(breedingData$records$genoVal[thisCycle])
	thisCycleAccuracy <- cor(breedingData$records$genoVal[thisCycle], breedingData$records$genoHat[thisCycle])
	breedingData$history <- rbind(breedingData$history, data.frame(cycle=cycle, genoBase=speciesData$genArch$genoBase, nQTLpoly=sum(apply(qtlDosage, 2, sd) > 0), nQTL=length(qtl), nMrkPoly=sum(apply(locDosage[,-qtl], 2, sd) > 0), nMrk=nLoc - length(qtl), meanGenoVal=genoMean, accuracy=thisCycleAccuracy, genoVar=genoVar, totChrVar=sum(chrVar, na.rm=TRUE), totSegVar=sum(segVar10[1,], na.rm=TRUE), equiVar=totEquiVar, selGenoValMain=breedingData$selGenoVal[1], selGenoValSpec=breedingData$selGenoVal[2]))
	
	# If there is mutation, recalibrate until equilibrium
	if (!(setting %in% c(321, 322, 323, 325, 326, 327, 401, 402, 411, 421))){
	if (exists("mutParms", where=speciesData)){
		reEquil <- switchCyc / 10
		if (cycle <= switchCyc & cycle %% reEquil == 0){
			meanGenoVar <- mean(breedingData$history$genoVar[(cycle - reEquil)+1:reEquil], na.rm=TRUE)
			h2 <- breedingData$heritability
			breedingData$stdDevErr <- sqrt(meanGenoVar * (1 - h2) / h2)
			cat(cycle, "Recalc stdDevErr:", breedingData$stdDevErr, "\n")
		}
		if (setting == 312 & cycle == switchCyc) speciesData$mutParms <<- NULL
	}
}
	return(breedingData)
}#END analysisObsQ

analysisPheno <- function(cycle, breedingData){
  cat(cycle, "analysisPheno", "\n")
  thisCycle <- which(breedingData$records$cycle == cycle)
  thin <- nCycles / 40
  # Calculate equiVar; segmentVar; chromosomeVar; genoVar
  nLoc <- nrow(speciesData$map)
  eo <- 1:nLoc*2
  locDosage <- (breedingData$genoMat[thisCycle, eo - 1] + breedingData$genoMat[thisCycle, eo]) / 2
  nInd <- nrow(locDosage)
  qtl <- speciesData$genArch$locusList
  qtlDosage <- locDosage[,qtl]
  qtlDosageVar <- apply(qtlDosage, 2, var)
  totEquiVar <- c(crossprod(speciesData$genArch$effects^2, qtlDosageVar))
  
  # Calculate the variance of a segment (could be a whole chromosome)
  # NOTE: assumes that the segment does NOT overlap two chromosomes
  # Assume the simplest additive gene action for now
  # Only do this once equilibrium reached (see below)
  # WARNING depends on calculating equiVar before hand
  calcSegmentVar <- function(startPos, segSize, calcEquiVar=FALSE){
    segGenoVals <- NA
    if (calcEquiVar) equiVar <- NA else equiVar <- NULL
    segLoci <- which(speciesData$map[,3] >= startPos & speciesData$map[,3] < startPos + segSize) # segLoci in locus order
    if (length(segLoci > 0)){
      segQTL <- which(qtl %in% segLoci) # segQTL in qtl order
      segGenoVals <- qtlDosage[,segQTL, drop=FALSE] %*% speciesData$genArch$effects[segQTL]
      if (calcEquiVar) equiVar <- crossprod(speciesData$genArch$effects[segQTL]^2, qtlDosageVar[segQTL])
    }
    return(c(var(segGenoVals), equiVar))
  }
  
  # Calculate segment variances every thin cycles on the whole genome
  chrVar <- NA
  segVar10 <- matrix(NA, 4, 105)
  if (cycle %% thin == 0){
    chrVar <- sapply(0:6 * 150, calcSegmentVar, segSize=150) # WARNING chr size of 150 hardwired
    seg10 <- rep(0:6 * 150, each=15) + 0:14 * 10
    segVar10 <- sapply(seg10, calcSegmentVar, segSize=10, calcEquiVar=TRUE)
    breedingData$localStats <- c(breedingData$localStats, list(chrVar, segVar10))
  }
  
  # Collect some information on this cycle
  genoMean <- mean(breedingData$records$genoVal[thisCycle])
  genoVar <- var(breedingData$records$genoVal[thisCycle])
  breedingData$history <- rbind(breedingData$history, data.frame(cycle=cycle, genoBase=speciesData$genArch$genoBase, nQTLNe=sum(apply(qtlDosage[sample(nInd, breedingData$nToSelect),], 2, sd) > 0), nQTL=length(qtl), nMrkNe=sum(apply(locDosage[sample(nInd, breedingData$nToSelect),-qtl], 2, sd) > 0), nMrk=nLoc - length(qtl), meanGenoVal=genoMean, genoVar=genoVar, totChrVar=sum(chrVar, na.rm=TRUE), totSegVar=sum(segVar10[1,], na.rm=TRUE), equiVar=totEquiVar, selGenoValMain=breedingData$selGenoVal))

  return(breedingData)
}#END analysisPheno

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
		breedingData$ancAllele <- breedingData$ancAllele[-monomorphic]
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
		mutAncAllele <- mutInfo["ancAllele",]
		newLoc <- mutInfo["locIdx",]
		oldLoc <- (1:nLoc)[-newLoc]
		allLocOrd <- order(c(oldLoc, newLoc))
		newGM <- breedingData$genoMat
		newGM <- cbind(newGM, matrix(rep(mutAncAllele, each=nrow(newGM)*2), nrow=nrow(newGM)))
		breedingData$genoMat <- rbind(newGM[,rep(allLocOrd, each=2)*2 - 1:0], pedNgeno$genoMat)
		breedingData$ancAllele <- c(breedingData$ancAllele, mutAncAllele)[allLocOrd]
	
		# The map and the genoMat are taken care of; deal with QTL and effects
		nNewQTL <- floor(nMut / speciesData$ratioLocToQTL)
		newQTLidx <- sample(nMut, nNewQTL)
		qtlAncAl <- mutAncAllele[newQTLidx]
		newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.95) * 2 - 1)
		newEff <- newEff * qtlAncAl
		allMutEff <- numeric(nMut)
		allMutEff[newQTLidx] <- newEff
		mutInfo <- rbind(mutInfo, allMutEff)
		newQTLidx <- newLoc[newQTLidx]
		oldQTLidx <- oldLoc[qtl]
		allQTL <- c(oldQTLidx, newQTLidx)
		allQTLord <- order(allQTL)
		speciesData$genArch$effects <<- c(speciesData$genArch$effects, newEff)[allQTLord]
		qtl <- allQTL[allQTLord]
		shiftToGenoBase <- shiftToGenoBase - 2 * crossprod(qtlAncAl, newEff)

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
	# cycle; Chr and Pos; nGenerations polymorphic; effect; ancAllele; fixedAllele
	# Which were poly last generation, but not now?  Store the values for those.
	lastCyc <- breedingData$records$cycle == cycle # (Just created cycle + 1)
	polyBefore <- apply(rbind(breedingData$genoMat[lastCyc, eo], breedingData$genoMat[lastCyc, eo-1]), 2, sd) > 0
	newlyMono <- which(polyBefore & !polyNow)
	effects <- numeric(length(newlyMono)) # Effect of the mutation
	effects[newlyMono %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyMono] * -breedingData$ancAllele[intersect(newlyMono, qtl)]
	nGenPoly <- breedingData$nGenPoly[newlyMono]
	ancAllele <- breedingData$ancAllele[newlyMono]
	fixedAllele <- pedNgeno$genoMat[1, newlyMono*2]
	monoInfo <- cbind(cycle, speciesData$map[newlyMono, c("Chr", "Pos")], nGenPoly, effects, ancAllele, fixedAllele)
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
		breedingData$ancAllele <- breedingData$ancAllele[-monomorphic]
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
		mutAncAllele <- mutInfo["ancAllele",]
		newLoc <- mutInfo["locIdx",]
		oldLoc <- (1:nLoc)[-newLoc]
		allLocOrd <- order(c(oldLoc, newLoc))
		newGM <- breedingData$genoMat
		newGM <- cbind(newGM, matrix(rep(mutAncAllele, each=nrow(newGM)*2), nrow=nrow(newGM)))
		breedingData$genoMat <- rbind(newGM[,rep(allLocOrd, each=2)*2 - 1:0], pedNgeno$genoMat)
		breedingData$ancAllele <- c(breedingData$ancAllele, mutAncAllele)[allLocOrd]
		if (exists("cumulativeIncidences", breedingData)) breedingData$cumulativeIncidences <- c(breedingData$cumulativeIncidences, numeric(nMut))[allLocOrd]
		
		# The map and the genoMat are taken care of; deal with QTL and effects
		nNewQTL <- floor(nMut / speciesData$ratioLocToQTL)
		newQTLidx <- sample(nMut, nNewQTL)
		qtlAncAl <- mutAncAllele[newQTLidx]
		newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.95) * 2 - 1)
		newEff <- newEff * qtlAncAl
		allMutEff <- numeric(nMut)
		allMutEff[newQTLidx] <- newEff
		mutInfo <- rbind(mutInfo, allMutEff)
		newQTLidx <- newLoc[newQTLidx]
		oldQTLidx <- oldLoc[qtl]
		allQTL <- c(oldQTLidx, newQTLidx)
		allQTLord <- order(allQTL)
		speciesData$genArch$effects <<- c(speciesData$genArch$effects, newEff)[allQTLord]
		qtl <- allQTL[allQTLord]
		shiftToGenoBase <- shiftToGenoBase - crossprod(qtlAncAl, newEff)

		# Adjust mrkEffEst, which has an est for _all_ loci
		# E(var(mrkEff)) = rate / rLQ + rate^2 / rLQ - (rate * (2*pFav - 1)/2)^2
		breedingData$mrkEffVar <- c(breedingData$mrkEffVar, rep(0.2476, nMut))[allLocOrd]
		if (setting %in% c(421)) priorMrkEff <- 0 else priorMrkEff <- 0.18
		breedingData$mrkEffEst <- c(breedingData$mrkEffEst, priorMrkEff * mutAncAllele)[allLocOrd]
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
  if (exists("meanPhen", breedingData)) breedingData$meanPhen <- breedingData$meanPhen + shiftToGenoBase
	breedingData$mutHist <- c(breedingData$mutHist, list(mutInfo))
	# For loci that just became monomorphic, save:
	# cycle; Chr and Pos; nGenerations polymorphic; effect; ancAllele; fixedAllele
	# Which were poly last generation, but not now?  Store the values for those.
	lastCyc <- breedingData$records$cycle == cycle # (Just created cycle + 1)
	polyBefore <- apply(rbind(breedingData$genoMat[lastCyc, eo], breedingData$genoMat[lastCyc, eo-1]), 2, sd) > 0
	newlyMono <- which(polyBefore & !polyNow)
	effects <- numeric(length(newlyMono)) # Effect of the mutation
	effects[newlyMono %in% qtl] <- speciesData$genArch$effects[qtl %in% newlyMono] * -breedingData$ancAllele[intersect(newlyMono, qtl)]
	nGenPoly <- breedingData$nGenPoly[newlyMono]
	ancAllele <- breedingData$ancAllele[newlyMono]
	fixedAllele <- pedNgeno$genoMat[1, newlyMono*2]
	monoInfo <- cbind(cycle, speciesData$map[newlyMono, c("Chr", "Pos")], nGenPoly, effects, ancAllele, fixedAllele)
	breedingData$monoInfo <- rbind(breedingData$monoInfo, monoInfo)
	
	# save(breedingData, speciesData, file=paste("bDsD",cycle,"RData", sep="."))	
	return(breedingData)
}#END cross advance obs QTL

crossAdvanceMutPheno <- function(cycle, breedingData){
  cat(cycle, "crossAdvanceMutatePheno", "\n")

  # Make progeny
  selectSet <- breedingData$selectSet
  selGenoMat <- breedingData$genoMat[selectSet,]
  pedNgeno <- randomMate(selGenoMat, speciesData$map, breedingData$nSelCan, speciesData$chrMax, speciesData$progType)
  mutInfo <- pedNgeno$mutInfo
  newIDs <- max(breedingData$records$ID) + 1:breedingData$nSelCan
  rownames(pedNgeno$genoMat) <- newIDs
  newRecords <- data.frame(ID=newIDs, MID=selectSet[pedNgeno$pedigree[1,]], PID=selectSet[pedNgeno$pedigree[2,]], cycle=cycle + 1, genoVal=NA, phenoVal=NA, genoHat=NA)
  breedingData$records <- newRecords
  
  qtl <- speciesData$genArch$locusList
  nLoc <- nrow(speciesData$map)
  eo <- 1:nLoc * 2
  speciesData$nLoc <<- c(speciesData$nLoc, nLoc)
  polyNow <- apply(rbind(pedNgeno$genoMat[, eo], pedNgeno$genoMat[, eo-1]), 2, sd) > 0
  
  # new genoMat has extra columns relative to old genoMat so add columns to the old
  nMut <- ncol(mutInfo)
  mutAncAllele <- mutInfo["ancAllele",]
  newLoc <- mutInfo["locIdx",]
  oldLoc <- (1:nLoc)[-newLoc]
  allLocOrd <- order(c(oldLoc, newLoc))
  breedingData$genoMat <- pedNgeno$genoMat
  breedingData$ancAllele <- c(breedingData$ancAllele, mutAncAllele)[allLocOrd]
  breedingData$cumulativeIncidences <- c(breedingData$cumulativeIncidences, numeric(nMut))[allLocOrd]
  if (cycle == 1) breedingData$nGenPoly <- rep(1, length(oldLoc))
  breedingData$nGenPoly <- c(breedingData$nGenPoly, numeric(nMut))[allLocOrd] + polyNow
  
  # Deal with QTL and effects
  nNewQTL <- floor(nMut / speciesData$ratioLocToQTL)
  newQTLidx <- sample(nMut, nNewQTL)
  qtlAncAl <- mutAncAllele[newQTLidx]
  newEff <- rgamma(nNewQTL, 0.4) * (rbinom(nNewQTL, 1, 0.95) * 2 - 1)
  newEff <- newEff * qtlAncAl
  allMutEff <- numeric(nMut)
  allMutEff[newQTLidx] <- newEff
  mutInfo <- rbind(mutInfo, allMutEff)
  newQTLidx <- newLoc[newQTLidx]
  oldQTLidx <- oldLoc[qtl]
  allQTL <- c(oldQTLidx, newQTLidx)
  allQTLord <- order(allQTL)
  speciesData$genArch$effects <<- c(speciesData$genArch$effects, newEff)[allQTLord]
  qtl <- allQTL[allQTLord]
  shiftToGenoBase <- -crossprod(qtlAncAl, newEff)
  
  monomorphic <- which(!polyNow)
  # Anything that is monomorphic now was polymorphic in the previous cycle
  # Save: cycle; Chr and Pos; nGenerations polymorphic; effect; ancAllele; fixedAllele
  # Which were poly last generation, but not now?  Store the values for those.
  effects <- numeric(length(monomorphic)) # Effect of the mutation
  effects[monomorphic %in% qtl] <- speciesData$genArch$effects[qtl %in% monomorphic] * -breedingData$ancAllele[intersect(monomorphic, qtl)]
  nGenPoly <- breedingData$nGenPoly[monomorphic]
  cumuIncid <- breedingData$cumulativeIncidences[monomorphic]
  ancAllele <- breedingData$ancAllele[monomorphic]
  fixedAllele <- pedNgeno$genoMat[1, monomorphic*2]
  monoInfo <- cbind(cycle, speciesData$map[monomorphic, c("Chr", "Pos")], nGenPoly, cumuIncid, effects, ancAllele, fixedAllele)
  breedingData$monoInfo <- rbind(breedingData$monoInfo, monoInfo)

  # Take monomorphic out of the map
  speciesData$map <<- speciesData$map[-monomorphic,]
  # out of the QTL effects
  monoQTL <- which(qtl %in% monomorphic)
  fixQTLeff <- speciesData$genArch$effects[monoQTL]
  fixQTLal <- breedingData$genoMat[1, qtl[monoQTL]*2]
  shiftToGenoBase <- shiftToGenoBase + crossprod(fixQTLal, fixQTLeff)
  speciesData$genArch$effects <<- speciesData$genArch$effects[-monoQTL]
  # out of the locus data
  breedingData$genoMat <- breedingData$genoMat[, -(rep(monomorphic*2, each=2)-0:1)]
  # out of nGenPoly
  breedingData$nGenPoly <- breedingData$nGenPoly[-monomorphic]
  breedingData$ancAllele <- breedingData$ancAllele[-monomorphic]
  breedingData$cumulativeIncidences <- breedingData$cumulativeIncidences[-monomorphic]
  # Renumber qtl and mrk
  qtl <- which((1:nLoc)[-monomorphic] %in% qtl)
  
  speciesData$genArch$locusList <<- qtl
  speciesData$genArch$genoBase <<- speciesData$genArch$genoBase + shiftToGenoBase
  breedingData$mutHist <- c(breedingData$mutHist, list(mutInfo))
  
  return(breedingData)
}#END cross advance pheno

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
