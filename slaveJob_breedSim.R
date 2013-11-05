# This is the function that should be given as the slaveJob function in the 
# MasterSlavePull generic parallel processing framework using Rmpi

# Some objects will come in through the jobParameters
# These are objects that change for every replication of every setting
	#  randSeed -- seed to initialize random number generator
	#  returnFileName -- name of the R file in which to store the results
slaveJob_breedSim <- function(jobParameters){
	if (DEBUG) options(warn=1) # Prints warnings when they happen
	
	for (parm in 1:length(jobParameters)) assign(names(jobParameters)[parm], jobParameters[[parm]])
	if (!exists("randSeed")) randSeed <- round(runif(1, 0, 1e9))
	if (!exists("returnFileName")) returnFileName <- "breedSimReturnFile.RData"
	if (exists("workingDirectory")) setwd(workingDirectory)
	
	# There are some objects that ALL settings in an experiment need.  Those should be
	# given in a baseline file. The breeding simulation functions are in this category.
	# The name of this baseline file can be given in the task list that is attached in
	# slaveBody else a generic name is assumed
	if (!exists("baselineFile")) baselineFile <- "baselineFile.RData"
	load(baselineFile, .GlobalEnv)

	# In experimental design terms, this is the complete block that is common
	# across all treatments. It may be that everyone has the same founders and
	# the only differentiation is the random seed.  That's "pseudo-replication"
	# The founder file contains
	#  replication -- an integer giving the block number
	#  breedingData -- a list with records and genoMat
	#  speciesData --  a list with attributes of the genome and breeding system 
	# that represent the population at the start of selection
	# It can contain replication, which would override (I hope) one given in the tasks list
	if (!exists("founderFile")) founderFile <- "founderFile.RData"
	load(founderFile, .GlobalEnv)
	
	# In experimental design terms, this is the treatment applied to an experimental
	# unit within the complete block.
	# The scheme file contains minimally:
	#  setting -- the setting name (or number) for reference
	#  nCycles -- the number of cycles of selection to be simulated
	# And the functions:
	#  initialize -- potentially modifies the founder population and prepares it to start breeding
	#  phenoList -- one or more functions that phenotype for selection or model updating
	#  selectList -- one or more functions that calculate the selection criterion
	#  crossAdvanceList -- one or more functions to determine crossing and generation advances
	#  processResults -- function that summarizes after each replication of the full breeding scheme
	if (!exists("schemeFile")) schemeFile <- "schemeFile.RData"
	print(schemeFile)
	load(schemeFile, .GlobalEnv)

	founderFile <- initialize(replication)

	set.seed(randSeed)

	dummy <- sapply(1:nCycles, oneLifecycle, phenotype, select, crossAdvance, analyze)

	breedingData <<- processResults(breedingData, speciesData, setting, replication)
	
	save(breedingData, speciesData, setting, replication, randSeed, file=returnFileName)

	if (!is.null(warnings())){
		theWarnings <- warnings()
		warningFile <- paste(substring(returnFileName, 1, nchar(returnFileName) - 6), "warnings.RData", sep="")
		print(warningFile)
		save(theWarnings, file=warningFile)
	}
	return(returnFileName)
}#END slaveJob_breedSim