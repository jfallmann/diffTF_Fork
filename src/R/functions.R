
initFunctionsScript <- function(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE, verbose = TRUE) {
  
  checkAndLoadPackages("checkmate", verbose = verbose)
  assert(checkNull(packagesReq), checkCharacter(packagesReq, min.len = 1, min.chars = 1))
  assertCharacter(minRVersion, len = 1)
  assertInt(warningsLevel, lower = 0, upper = 2)
  assertFlag(disableScientificNotation)
  assertFlag(verbose)
  
  clearOpenDevices()
  
  
  # No annoying strings as factors by default
  options(stringsAsFactors = FALSE)
  
  # Print warnings as they occur
  options(warn = warningsLevel)
  
  # Just print 50 lines instead of 99999
  options(max.print = 200)
  
  # Disable scientific notation
  if (disableScientificNotation) options(scipen = 999)
  
  
  # We need at least R version 3.1.0 to continue
  stopifnot(getRversion() >= minRVersion)
  
  if (is.null(packagesReq)) {
    #packagesReq = .loadAllLibraries()
  }
  
  
  .detachAllPackages()
  
  checkAndLoadPackages(packagesReq, verbose = verbose)
  
}

checkAndLoadPackages <- function(packages, verbose = TRUE) {
  
  .checkAndInstallMissingPackages(packages, verbose = verbose)
  
  for (packageCur in packages) {
    library(packageCur, character.only = TRUE)
  }
  
}


startLogger <- function(logfile, level, removeOldLog = TRUE, appenderName = "consoleAndFile", verbose = TRUE) {
  
  checkAndLoadPackages(c("futile.logger"), verbose = verbose)
  
  assertSubset(level, c("TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
  assertFlag(removeOldLog)
  assertSubset(appenderName, c("console", "file", "consoleAndFile"))
  assertFlag(verbose)
  
  if (appenderName != "console") {
    assertDirectory(dirname(logfile), access = "w")
    if (file.exists(logfile)) {
      file.remove(logfile)
    }
  }
  
  # LEVELS: TRACE, DEBUG, INFO, WARN, ERROR, FATAL
  invisible(flog.threshold(level))
  
  
  if (appenderName == "console") {
    invisible(flog.appender(appender.console(file = logfile)))
  } else if (appenderName == "file")  {
    invisible(flog.appender(appender.file(file = logfile)))
  } else {
    invisible(flog.appender(appender.tee(file = logfile)))
  }
  
  
}

printParametersLog <- function(par.l, verbose = FALSE) {
  
  checkAndLoadPackages(c("futile.logger"), verbose = verbose)  
  assertList(par.l)
  flog.info(paste0("PARAMETERS:"))
  for (parCur in names(par.l)) {
    
    flog.info(paste0(" ", parCur, "=",  paste0(par.l[[parCur]], collapse = ",")))
    
  }
}




###########################################
# PACKAGE LOADING AND DETACHING FUNCTIONS #
###########################################

.checkAndInstallMissingPackages <- function(packages.vec, verbose = TRUE) {
  
  if (verbose) cat("Trying to automatically install missing packages. If this fails, install them manually...\n")
  
  packagesToInstall = setdiff(packages.vec, rownames(installed.packages()))
  
  
  if (length(packagesToInstall) > 0) {
    if (verbose) cat("Could not find the following packages: ", paste( packagesToInstall , collapse = ", "), "\n")
    install.packages(packagesToInstall, repos = "http://cran.rstudio.com/")  
    
    source("http://bioconductor.org/biocLite.R")
    for (packageCur in packagesToInstall) {
      biocLite(packageCur, suppressUpdates = TRUE)
    }
  } else {
    if (verbose) cat("All packages are already installed\n")
  }
  
}


.detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search())) == 1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list) > 0)  for (package in package.list) detach(package, character.only = TRUE)
  
}



clearOpenDevices <- function() {
  
  while (length(dev.list()) > 0) {
    dev.off()
  }
}

getComparisonFromDeSeqObject <- function(DeSeq.obj, designFormula) {
  
  assertClass(DeSeq.obj, "DESeqDataSet")
  
  mcol.df = mcols(results(DeSeq.obj))
  split1 = strsplit(mcol.df[which(mcol.df$type == "results"),][1,"description"], split = ":", fixed = TRUE)[[1]][2]
  
  formulaElements = gsub(pattern = "[~+*]", replacement = "", x = designFormula)
  elems = strsplit(formulaElements, split = "\\s+", perl = TRUE)[[1]]
  splitFinal = c()
  for (i in seq_len(length(elems))) {
    splitFinal = trimws(gsub(pattern = elems[i], replacement = "", split1))
  }
  
  return(splitFinal)
}


permuteSampleTable <- function(file_sampleTable, conditionComparison, factorVariableInFormula, stepsize, nRepetitionsPerStepMax) {
  
  checkAndLoadPackages(c("tidyverse", "checkmate"), verbose = TRUE)
  
  # Check if contrasts have been specified correctly
  conditionsContrast = strsplit(conditionComparison, split = ",", fixed = TRUE)[[1]]
  assertVector(conditionsContrast, len = 2)
  
  assertFileExists(file_sampleTable, access = "r")
  
  assertIntegerish(stepsize, lower = 1)
  assertIntegerish(nRepetitionsPerStepMax, lower = 1)
  
  sampleData.df = read_tsv(file_sampleTable, col_names = TRUE)
  
  assertSubset(conditionsContrast, sampleData.df$conditionSummary)
  
  # Get original ratio
  
  conditionCounter = table(sampleData.df$conditionSummary)
  
  nSamplesRareCondition     = min(conditionCounter)
  nSamplesFrequentCondition = max(conditionCounter)
  ratio = nSamplesRareCondition / (nSamplesFrequentCondition + nSamplesRareCondition)
  
  nameRareCondition      = names(conditionCounter)[conditionCounter == min(conditionCounter)]
  nameFrequentCondition  = names(conditionCounter)[conditionCounter == max(conditionCounter)]
  indexRareCondition     = which(sampleData.df$conditionSummary == nameRareCondition)
  indexFrequentCondition = which(sampleData.df$conditionSummary == nameFrequentCondition)
  
  stopifnot(length(indexRareCondition) == nSamplesRareCondition)
  
  if (nSamplesRareCondition > 10) {
    samplesRareCondition = c(2:5, seq(10, nSamplesRareCondition, stepsize))
    samplesRareCondition = c(3:9, seq(10, nSamplesRareCondition, stepsize))
    
    if (nSamplesRareCondition %% stepsize != 0) {
      samplesRareCondition = c(samplesRareCondition, nSamplesRareCondition)
    }
  } else {
    samplesRareCondition = c(2:nSamplesRareCondition)
  }
  
  nSamplesBase = length(samplesRareCondition)
  
  samplesFrequentCondition = ceiling(samplesRareCondition * (1/ratio - 1))
  
  # Correct rounding errors
  if (samplesFrequentCondition[nSamplesBase] > nSamplesFrequentCondition) samplesFrequentCondition[nSamplesBase] = nSamplesFrequentCondition
  
  # for each particular number of samples for the rare case
  
  subsamples.l = list()
  for (sampleBaseCur in 1:nSamplesBase) {
    
    nValidSamples = 0
    nPermutations = 0
    while (nValidSamples < nRepetitionsPerStepMax || nPermutations > 100) {
      
      nPermutations = nPermutations + 1
      # 1. Rare Condition 
      # How many different samples are actually possible?
      nCombinations = choose(nSamplesRareCondition, samplesRareCondition[sampleBaseCur])
      nSamplesCur = min(nCombinations, nRepetitionsPerStepMax)
      

      # Generate them
      table.l = list()
      sampleCombinations.l = list()
      for (i in 1:nSamplesCur) {
        
        table.l[[i]] = sample_n(sampleData.df[indexRareCondition,], samplesRareCondition[sampleBaseCur], replace = FALSE)
        
      }
      
      
      # 2. Frequent Condition 
      
      # How many different samples are actually possible?
      nCombinations = choose(nSamplesFrequentCondition, samplesFrequentCondition[sampleBaseCur])
      nSamplesCur = min(nCombinations, nRepetitionsPerStepMax)
      
      # Generate them
      table2.l = list()
      for (i in 1:nSamplesCur) {
        table2.l[[i]] = sample_n(sampleData.df[indexFrequentCondition,], samplesFrequentCondition[sampleBaseCur], replace = FALSE)
        
      }
      
      # Merge
      for (i in 1:min(length(table.l), length(table2.l))) {
        
        if (nValidSamples == nRepetitionsPerStepMax) {
          break
        }
        
        listname = paste0(samplesRareCondition[sampleBaseCur], "_", i)
        subsamples.l[[listname]]  = bind_rows(table.l[[i]], table2.l[[i]])
        
        
        # Check validity of sample
        currentPermutationSampleComb = sort(unique(subsamples.l[[listname]]$SampleID))
        nUnique = length(unique(unlist(subsamples.l[[listname]][,factorVariableInFormula])))
        
        if (any(sapply(sampleCombinations.l, function(x) {identical(currentPermutationSampleComb, x)} )) || nUnique == 1) {
          
          # sample invalid, redo
          sampleValid = FALSE
        } else {
          
          # sample unique and ok
          sampleCombinations.l[[i]] = currentPermutationSampleComb
          
          sampleValid = TRUE
          nValidSamples = nValidSamples + 1
        }
        
      } # end merge
      
    } # end while not enough valid samples
    
    
  }
  
  subsamples.l
}

testExistanceAndCreateDirectoriesRecursively <- function(directories) {
  
  for (dirname in unique(directories)) {
    
    if (!testDirectoryExists(dirname)) {
      dir.create(dirname, recursive = TRUE)
    } else {
      assertDirectoryExists(dirname, access = "w")
    }
    
  }
}


convertToFormula <- function(userFormula, validColnames = NULL) {
  
  # Create formula based on user-defined design
  designFormula = tryCatch({
    as.formula(userFormula)
  }, warning = function(w) {
    stop("Converting the design formula \"", userFormula, "\" created a warning, which should be checked carefully.")
  }, error = function(e) {
    stop("Design formula \"", userFormula, "\" not valid")
  })
  
  
  # Check colmn names
  if (!is.null(colnames)) {
    formulaVariables = attr(terms(designFormula), "term.labels")
    assertSubset(formulaVariables, validColnames)
  }

  
  designFormula
  
}
