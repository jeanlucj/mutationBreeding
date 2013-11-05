####################################################################################
####################################################################################
# Run simulations based on parameters in files specified by
# modificationPrograms.txt
#
# This code parses objects in the modification programs files to create a list
# of tasks. That list is then sent to an apply function that parallelizes code.
# At least three files need to be associated with each task
# baselineFile -- All settings in the experiment will use the functions in there
# founderFile -- Corresponds to the block, with breedingData and speciesData
# schemeFile -- Corresponds to the treatment, with setting, nCycles, initialize, 
#	phenoList, selectList, and crossAdvanceList
# install.packages("multicore", lib="/Library/Frameworks/R.framework/Versions/3.0/Resources/library", repos="http://lib.stat.cmu.edu/R/CRAN")
####################################################################################

# if you invoke from the command line, here are the arguments
# EXPERIMENT FILE NAME, -e, a file that has objects common to generating all settings
# MODIFICATION START, -m, where to start in modificationPrograms.txt
# the assumption is that if START == 1, then initial conditions need to be generated
# else they have been generated before and can be loaded
# OUTDIR, -o, a unique folder for this series of initial conditions
# OUTDIR needs to be a valid directory
# MAKE INPUT FILES, -i [TRUE / FALSE]. If TRUE generate all the input files. If FALSE generate only the scheme files.
# if START == 1 and MAKE FOUNDER FILES is TRUE, then a new series of founder files will be generated in OUTDIR, else
# a previously generated series will be found in OUTDIR
# MAXP, -p, maximum number of processors to use
# SEED, -s, a random number generator seed
# REPLICATION OFFSET, -r, a number to add to the replication for naming: multiple computers on the same setting
# DEBUG, -d, whether to put the analysis in debugging mode (use lapply instead of mclapply)
# example:
# nohup R --vanilla --slave --args -e ./experimentFile.R -m 1 -o ./testOut -i TRUE -p 2 -s 12345678 -r 20 -d FALSE < master_breedSim.R > output.txt &

####################################################################################
# The experiment file has objects useful to generating datasets for all settings
# Minimally, it includes:
# experimentName -- a character variable that is a reference
# breedFuncList -- a character vector with the names of all the breeding simulation functions needed
# Additionally it may contain:
# modificationsFile -- the name of the file pointing to modification programs
### Note: the experiment file and the selection program file are sourced one after the
### other so the selection program file can provide or supercede anything in the experiment file
# nReplications -- the number of replications for these settings
# processResults -- a function that summarizes after each replication of the full breeding scheme
# createSpeciesData -- a function to create species data (genome and breeding system)
# createBreedingData -- a function to create the base population for each rep

arguments <- commandArgs(trailingOnly = TRUE)
# Get a vector of argument strings and return a named list
extractArguments <- function(arguments){
	whereDash <- grep("-", arguments, fixed=TRUE)
	argNames <- substring(arguments[whereDash], 2)
	whereDash <- c(whereDash, length(arguments)+1)
	argList <- NULL
	for (argIdx in 1:(length(whereDash)-1)){
		argList <- c(argList, list(arguments[(whereDash[argIdx]+1):(whereDash[argIdx+1]-1)]))
	}
	names(argList) <- argNames
	return(argList)
}
argList <- extractArguments(arguments)
print(argList)
if ("o" %in% names(argList)) outFolder <- argList$o else outFolder <- paste("./simOut", as.integer(Sys.time()), round(runif(1, 0, 1e9)), sep=".")
	print(paste("Simulation output folder", outFolder))
if ("m" %in% names(argList)) modifStart <- as.integer(argList$m) else modifStart <- 1
if (!exists("modificationsFile")) modificationsFile <- "modificationPrograms.txt"
if ("i" %in% names(argList)) makeFounderFiles <- as.logical(argList$i) else makeFounderFiles <- TRUE
if ("p" %in% names(argList)) nCoresToUse <- as.integer(argList$p) else nCoresToUse <- 2
if ("s" %in% names(argList)) seed <- as.integer(argList$s) else seed <- as.integer(Sys.time())
	set.seed(seed)
if ("r" %in% names(argList)) repOffset <- as.integer(argList$r) else repOffset <- 0
if ("d" %in% names(argList)) DEBUG <- as.logical(argList$d) else DEBUG <- FALSE
if ("n" %in% names(argList)) system(paste("sleep", as.integer(argList$n)))

############################################################################################
### If needed, make space in the file system
############################################################################################
simInpFolder <- paste(outFolder, "simInp", sep="/")
simOutFolder <- paste(outFolder, "simOut", sep="/")
if (modifStart == 1 & makeFounderFiles){
	timeStamp <- Sys.time()
	try(system(paste("mkdir", outFolder)), silent = TRUE)
	try(system(paste("mkdir", simInpFolder)), silent = TRUE)
	try(system(paste("mkdir", simOutFolder)), silent = TRUE)
}
if ("e" %in% names(argList)) source(argList$e) else source("experimentFile.R")
	baselineFile <- paste(outFolder, "/baselineFile.", experimentName, ".RData", sep="")
	print(baselineFile)
	save(list=breedFuncList, file=baselineFile)
if (DEBUG) nReplications <- 2

############################################################################################
### The function that runs each replicate of each setting
############################################################################################
source("slaveJob_breedSim.R")

############################################################################################
### Loop through settings then through replications to create task list
############################################################################################
modificationPrograms <- readLines(modificationsFile)
for (modification in modifStart:length(modificationPrograms)){
	
	# get the name of the selection program script to modify the baseline
	# the selection programs should contain
	# settings -- a vector of the name of the settings (= breeding schemes = treatments) in the experiment
	# allNcycles -- a vector indicating the number of selection cycles for each setting
	# allInitialize -- a list of the initialization function for each setting
	# allLists -- a list of lists. Each sub-list contains
	#	phenoList -- a list of phenotyping functions over cycles
	#	selectList -- a list of selection functions over cycles
	#	analysisList -- a list of functions to calculate summaries and digest information over cycles
	#	crossAdvanceList -- a list of functions to cross and advance generations over cycles
	selectionProgram <- modificationPrograms[modification]
	source(selectionProgram)
	dummy <- settingReminder()
	
	# Run through the different settings
	schemeFileVec <- character(length(settings))
	
	makeSchemeFile <- function(setting, nCycles, functionList){
		initialize <- functionList$initialize
		phenotype <- functionList$phenotype
		select <- functionList$select
		analyze <- functionList$analyze
		crossAdvance <- functionList$crossAdvance
		processResults <- functionList$processResults
		schemeFile <- paste(simInpFolder, "/schemeFile", setting, ".RData", sep="")
		save(setting, nCycles, initialize, phenotype, select, analyze, crossAdvance, processResults, settingReminder, file=schemeFile)
		return(schemeFile)
	}#END makeSchemeFile
	schemeFileVec <- mapply(makeSchemeFile, settings, allNcycles, allLists)
	
	# Make the tasks
	taskFile <- paste(simInpFolder, "/taskFile.", experimentName, ".RData", sep="")
	tasks <- NULL
	for (setting in settings){
		for (rep in 1:nReplications){
			# Set up founder file
			founderFile <- paste(simInpFolder, "/founderFile", setting, ".", rep+repOffset, ".RData", sep="")
			if (makeFounderFiles){
				replication <- rep+repOffset
				founderSeed <- round(runif(1, 0, 1e9))
				save(replication, founderSeed, founderFile, file=founderFile)
			}
			
			# Set up tasks
			whichSet <- which(settings == setting)
			randSeed <- round(runif(1, 0, 1e9))
			returnFileName <- paste(simOutFolder, "/breedSimOut", setting, ".", rep+repOffset, ".RData", sep="")
			jobParameters <- list(randSeed=randSeed, returnFileName=returnFileName, baselineFile=baselineFile, founderFile=founderFile, schemeFile=schemeFileVec[whichSet])
			tasks <- c(tasks, list(list(slaveJob=slaveJob_breedSim, jobParameters=jobParameters)))
		}
	}
	save(tasks, file=taskFile)			
}#END loop through the modifications

############################################################################################
# Everything should be set up for the multitasking system to do it's work
# It needs
# nCoresToUse -- how many cores to use
# tasks -- a list of lists
# 	each list must have 
#	slaveJob -- a function that the slave executes 
#	jobParameters -- a list of parameters guiding the job
############################################################################################
objList <- c(ls(), "objList")
objList <- setdiff(objList, c("tasks", "nCoresToUse", "DEBUG", "makeFounderFiles"))
rm(list=objList)
doJobWithParameters <- function(jobAndParms){
	return(jobAndParms$slaveJob(jobAndParms$jobParameters))
}
if (DEBUG){ # Use lapply to get error messages
	debugPlace <- 0
	results <- lapply(tasks, doJobWithParameters)
	options(warn=1) # Prints warnings when they happen
	print(warnings())
} else {
	library("multicore")
	results <- mclapply(tasks, doJobWithParameters, mc.cores=nCoresToUse)
}