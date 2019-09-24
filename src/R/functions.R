
initFunctionsScript <- function(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE, verbose = FALSE) {
  
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
  
  # Just print 200 lines instead of 99999
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

checkAndLoadPackages <- function(packages, verbose = FALSE) {
  
  .checkAndInstallMissingPackages(packages, verbose = verbose)
  
  for (packageCur in packages) {
    suppressMessages(library(packageCur, character.only = TRUE))
  }
  
}


startLogger <- function(logfile, level, removeOldLog = TRUE, appenderName = "consoleAndFile", verbose = TRUE) {
  
  checkAndLoadPackages(c("futile.logger"), verbose = FALSE)
  
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
    invisible(flog.appender(appender.console()))
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

.checkAndInstallMissingPackages <- function(packages.vec, verbose = FALSE) {
  
  if (verbose) cat("Trying to automatically install missing packages. If this fails, install them manually...\n")
  
  packagesToInstall = setdiff(packages.vec, rownames(installed.packages()))
  
  
  if (length(packagesToInstall) > 0) {
    if (verbose) cat("Could not find the following packages: ", paste( packagesToInstall , collapse = ", "), "\n")
    install.packages(packagesToInstall, repos = "http://cran.rstudio.com/")  
    
    if (getRversion() < "3.5.0") { 
        
        source("http://bioconductor.org/biocLite.R")
        for (packageCur in packagesToInstall) {
            biocLite(packageCur, suppressUpdates = TRUE)
        }
    } else {
        
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        
        for (packageCur in packagesToInstall) {
            BiocManager::install(packageCur)
        }
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


createFileList <- function(directory, pattern, recursive = FALSE, ignoreCase = FALSE, verbose = TRUE) {
  
  assertCharacter(directory, min.chars = 1, any.missing = FALSE, len = 1)
  assertCharacter(pattern, min.chars = 1, any.missing = FALSE, len = 1)
  assertFlag(recursive)  
  assertFlag(ignoreCase)  
  assertFlag(verbose)    
  
  assertDirectoryExists(directory)
  
  # Multiple patterns are now supported, integrate over them
  patternAll = strsplit(pattern, ",")[[1]]
  assertCharacter(patternAll, min.len = 1)
  
  if (verbose) cat("Found ", length(patternAll), " distinct pattern(s) in pattern string.\n")
  
  nFilesToProcessTotal = 0
  filesToProcess.vec = c()
  
  for (patternCur in patternAll) {
    
    # Replace wildcards by functioning patterns (such as .)
    patternMod = glob2rx(patternCur)
    
    # Remove anchoring at beginning and end
    patternMod = substr(patternMod, 2, nchar(patternMod) - 1)
    
    filesToProcessCur.vec = list.files(path = directory, pattern = patternMod, full.names = TRUE, recursive = recursive, ignore.case = ignoreCase)
    filesToProcess.vec = c(filesToProcess.vec, filesToProcessCur.vec)
    
    if (verbose) cat("Search for files with pattern \"", patternCur, "\" in directory ", directory, " (case insensitive:", ignoreCase, ")\n", sep ="")
    
    nFilesToProcessTotal = nFilesToProcessTotal + length(filesToProcessCur.vec)
  }
  
  
  
  if (nFilesToProcessTotal == 0) {
    stop(paste0("No files to process in folder ", directory, " that fulfill the desired criteria (", patternCur, ")."))
  } else {
    
    if (verbose) cat("The following", nFilesToProcessTotal, "files were found:\n", paste0(filesToProcess.vec, collapse = "\n "))
  }
  
  if (verbose) cat("\n")
  
  filesToProcess.vec
  
}

getComparisonFromDeSeqObject <- function(DeSeq.obj, designFormula, datatypeVariableToPermute = "factor") {
  
  assertClass(DeSeq.obj, "DESeqDataSet")
  
  mcol.df = mcols(results(DeSeq.obj))
  resultsRow = mcol.df[which(mcol.df$type == "results"),][1,"description"]
  split1 = strsplit(resultsRow, split = ":", fixed = TRUE)[[1]][2]
  
  if (datatypeVariableToPermute == "factor") {
    
    formulaElements = gsub(pattern = "[~+*]", replacement = "", x = designFormula)
    elems = strsplit(trimws(formulaElements), split = "\\s+", perl = TRUE)[[1]]
    
    # Remove names of variables from split1 variable
    split1New = split1
    for (varCur in elems) {
      split1New = trimws(gsub(pattern = varCur, replacement = "", split1New))
    }
    
    splitNew = trimws(gsub(pattern = "vs", replacement = "", split1New))
    splitFinal = strsplit(splitNew, split = "\\s+", perl = TRUE)[[1]]
    
    
  } else if (datatypeVariableToPermute == "integer") {
    
    # Here, we manually assign neg. and pos. change
    splitFinal = trimws(split1)
    splitFinal = paste0(splitFinal, " (", c("pos.", "neg."), " change)")
    
  }
  
  assertVector(splitFinal, len = 2)
  return(splitFinal)
}


permuteSampleTable <- function(file_sampleTable, conditionComparison, factorVariableInFormula, stepsize, nRepetitionsPerStepMax) {
  
  checkAndLoadPackages(c("tidyverse", "checkmate"), verbose = FALSE)
  
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

###########
# DESeq 2 #
###########

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
    
    indexInteractionTerms = which(grepl(":", formulaVariables))
    if (length(indexInteractionTerms) > 0) {
        
        interactionTerms = formulaVariables[indexInteractionTerms]
        # Check whether all individual items have been specified
        for (interactionTermCur in interactionTerms) {
            componentsCur = strsplit(interactionTermCur, ":")[[1]]
            if (!all(componentsCur %in% formulaVariables)) {
                message = paste0("Design formula is incorrect, not all terms from the interactions (", paste0(componentsCur, collapse = ","), ") have been specified individually.")
                checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
                
            }
            
        }
        formulaVariables = formulaVariables[-indexInteractionTerms]
        
        
    }
    
    assertSubset(formulaVariables, validColnames)
  }

  
  designFormula
  
}


glog2 <- function(x) ((asinh(x) - log(2))/log(2))

# Credits to Bernd Klaus
myMAPlot <- function(M, idx, main, minMean = 0) {
  M <- na.exclude(M[, idx])
  M <- (glog2(M))
  mean <- rowMeans(M)
  M <- M[mean > minMean, ]
  mean <- mean[mean > minMean]
  
  difference <- M[, 2] - M[, 1]
  
  main <- paste(main, ", Number of genes:", dim(M)[1]) 
  
  pl <- (qplot(mean, difference, main = main, ylim = c(-5,5), asp = 2, geom = "point", alpha = I(.5), color = I("grey30"), shape = I(16))  #  
         + geom_hline(aes(yintercept = 0), col = "#9850C3", show.legend = FALSE)
         + geom_smooth(method = "loess", se = FALSE, col = "#5D84C5", span = .4)
         + theme_bw()
  )
  
  return(pl)
}


plotDiagnosticPlots <- function(dd, differentialResults, conditionComparison, filename = NULL, maxPairwiseComparisons = 5, plotMA = FALSE, alpha = 0.05) {
  
  checkAndLoadPackages(c("tidyverse", "checkmate", "geneplotter", "DESeq2", "vsn", "RColorBrewer", "limma"), verbose = FALSE)
    flog.info(paste0("Plotting various diagnostic plots"))
  
  assertClass(dd, "DESeqDataSet")
  
  assert(testClass(differentialResults, "MArrayLM"), testClass(differentialResults, "DESeqDataSet"))
  assertVector(conditionComparison, len = 2)
  assert(checkNull(filename), checkDirectory(dirname(filename), access = "w"))
  
  if (!is.null(filename)) {
    pdf(filename)
  }
  
  if (plotMA) {
      
      
      
      if (testClass(differentialResults, "MArrayLM")) {
          title = paste0("limma results\n", conditionComparison[1], " vs. ", conditionComparison[2])
          isSign = ifelse(p.adjust(differentialResults$p.value[,ncol(differentialResults$p.value)], method = "BH") < alpha, paste0("sign. (BH, ", alpha, ")"), "not-significant")
          
          flog.info(paste0(" Plotting MA plots from limma..."))
          
          limma::plotMA(differentialResults, main = title, status = isSign)
          
      } else {
          
          # DESeq2 specific MA plots
          flog.info(paste0(" Plotting MA plot from DESeq (1)..."))
          # 1. Regular MA plot:
          # shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
          DESeq2::plotMA(dd, main = "Regular MA plot")
          
          # 2. MA plot based on shrunken log2 fold changes
          # Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the dds object to the function  lfcShrink. Below we specify to use the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.
          
          # It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
          flog.info(paste0(" Plotting MA plot from DESeq (2)..."))
          dd_LFC <- lfcShrink(dd, coef=resultsNames(dd)[length(resultsNames(dd))], type="apeglm")
          
          DESeq2::plotMA(dd_LFC, main = "MA plot based on shrunken log2 fold changes")
          
      }
  }
  
 

  
  # 2. Densities of counts for the different samples. 
  # Since most of the genes are (heavily) affected by the experimental conditions, a succesful normalization will lead to overlapping densities
  
  nSamples = nrow(colData(dd))
  colors = colorRampPalette(brewer.pal(9, "Set1"))(nSamples)
  
  xlabCur = "Mean log counts"
  flog.info(paste0(" Plotting multidensity plot..."))
   
  if (nSamples < 10) {

    multidensity(log(DESeq2::counts(dd, normalized = FALSE) + 0.5), xlab = xlabCur, main = "Non-normalized log counts", col = colors)
    multidensity(log(DESeq2::counts(dd, normalized = TRUE) + 0.5) , xlab = xlabCur, main = "Normalized log counts", col = colors) 
    
    multiecdf(log(DESeq2::counts(dd, normalized = FALSE) + 0.5), xlab = xlabCur, main = "Non-normalized log counts", col = colors)
    multiecdf(log(DESeq2::counts(dd, normalized = TRUE) + 0.5) , xlab = xlabCur, main = "Normalized log counts", col = colors) 
    
  } else {
    
    warning("Omitting legend due to large number of samples (threshold is 10 at the moment)")
    
    multidensity(log(DESeq2::counts(dd, normalized = FALSE) + 0.5), xlab = xlabCur, main = "Non-normalized log counts", legend = NULL, col = colors)
    multidensity(log(DESeq2::counts(dd, normalized = TRUE) + 0.5) , xlab = xlabCur, main = "Normalized log counts",     legend = NULL, col = colors)
    
    multiecdf(log(DESeq2::counts(dd, normalized = FALSE) + 0.5), xlab = xlabCur, main = "Non-normalized log counts", legend = NULL, col = colors)
    multiecdf(log(DESeq2::counts(dd, normalized = TRUE) + 0.5) , xlab = xlabCur, main = "Normalized log counts",     legend = NULL, col = colors)
  }
  
  if (maxPairwiseComparisons > 0) {
      
      flog.info(paste0(" Plotting pairwise sample comparisons. This amy take a while."))
      
      # 3. Pairwise sample comparisons.
      # To further assess systematic differences between the samples, we can also plot pairwise mean–average plots: We plot the average of the log–transformed counts vs the fold change per gene for each of the sample pairs.
      MA.idx = t(combn(seq_len(dim(colData(dd))[1]), 2))
      
      if (nrow(MA.idx) > maxPairwiseComparisons) {
        
        flog.info(paste0("The number of pairwise comparisons to plot exceeds the current maximum of ", maxPairwiseComparisons, ". Only ", maxPairwiseComparisons, " pairwise comparisons will be shown in the PDF."))
        MA.idx.filt = MA.idx[1:maxPairwiseComparisons,, drop = FALSE]
    
      } else {
        MA.idx.filt = MA.idx
      }
      
      for (i in seq_along(MA.idx.filt[,1])) { 
        
        flog.info(paste0(" Plotting pairwise comparison ", i, " out of ", nrow(MA.idx.filt)))
        label = paste0(colnames(dd)[MA.idx.filt[i,1]], " vs ", colnames(dd)[MA.idx.filt[i,2]])
        suppressWarnings(print(myMAPlot(DESeq2::counts(dd, normalized = TRUE), c(MA.idx[i,1], MA.idx.filt[i,2]), main =  label)))
      }
      
      # Show an empty page with a warning if plots have been omitted
      if (nrow(MA.idx) > nrow(MA.idx.filt)) {
        
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        message = paste0("All remaining pairwise comparisons plots\nbetween samples have been omitted\nfor time and memory reasons.\nThe current maximum is set to ", maxPairwiseComparisons, ".")
        text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
      }
  
  }
  
  # 4. Mean SD plot: Plot row standard deviations versus row means
  notAllZeroPeaks <- (rowSums(DESeq2::counts(dd)) > 0)
  
  suppressWarnings(meanSdPlot(assay(dd[notAllZeroPeaks,])))
  
  if (!is.null(filename)) {
    dev.off()
  }
}

createDebugFile <- function(snakemake) {
    
    checkAndLoadPackages(c("checkmate", "tools", "futile.logger"), verbose = FALSE)
    
    if (!testClass(snakemake, "Snakemake")) {
        flog.warn(paste0("Could not find snakemake object, therefore not saving anyting."))
        return()
    }
    
    logfile = snakemake@log[[1]][1]
    
    assertCharacter(logfile, any.missing = FALSE)
    
    filename = paste0(tools::file_path_sans_ext(logfile), ".rds")
    
    # flog.info(paste0("Saved Snakemake object for manually rerunning the R script to ", filename))
    
    saveRDS(snakemake, filename)
    
}

#' @importFrom BiocParallel multicoreWorkers MulticoreParam
#' @import checkmate 
.initBiocParallel <- function(nWorkers, verbose = FALSE) {
  
  checkAndLoadPackages(c("BiocParallel"), verbose = verbose)  
  assertInt(nWorkers, lower = 1)    
  assertInt(multicoreWorkers())
  
  if (nWorkers > multicoreWorkers()) {
    warning("Requested ", nWorkers, " CPUs, but only ", multicoreWorkers(), " are available and can be used.")
    nWorkers = multicoreWorkers()
  }
  
  MulticoreParam(workers = nWorkers, progressBar = TRUE, stop.on.error = TRUE)
  
}

#' @importFrom BiocParallel bplapply bpok
#' @import checkmate 
.execInParallelGen <- function(nCores, returnAsList = TRUE, listNames = NULL, iteration, abortIfErrorParallel = TRUE, verbose = TRUE, functionName, ...) {
  
  checkAndLoadPackages(c("BiocParallel"), verbose = verbose)  
  start.time  <-  Sys.time()
  
  assertInt(nCores, lower = 1)
  assertFlag(returnAsList)
  assertFunction(functionName)
  assertVector(iteration, any.missing = FALSE, min.len = 1)
  assert(checkNull(listNames), checkCharacter(listNames, len = length(iteration)))
  
  res.l = list()
  
  if (nCores > 1) {
    
    res.l = tryCatch( {
      bplapply(iteration, functionName, ..., BPPARAM = .initBiocParallel(nCores))
      
    }#, error = function(e) {
    #     warning("An error occured while executing the function with multiple CPUs. Trying again using only only one CPU...")
    #     lapply(iteration, functionName, ...)
    # }
    )
    
    failedTasks = which(!bpok(res.l))
    if (length(failedTasks) > 0) {
      warning("At least one task failed while executing in parallel, attempting to rerun those that failed: ",res.l[[failedTasks[1]]])
      if (abortIfErrorParallel) stop()
      
      res.l = tryCatch( {
        bplapply(iteration, functionName, ..., BPPARAM = .initBiocParallel(nCores), BPREDO = res.l)
        
      }, error = function(e) {
        warning("Once again, an error occured while executing the function with multiple CPUs. Trying again using only only one CPU...")
        if (abortIfErrorParallel) stop()
        lapply(iteration, functionName, ...)
      }
      )
    }
    
  } else {
    res.l = lapply(iteration, functionName, ...)
  }
  
  end.time  <-  Sys.time()
  
  if (nCores > multicoreWorkers()) {
    nCores = multicoreWorkers()
  }
  
  flog.info(paste0(" Finished execution using ",nCores," cores. TOTAL RUNNING TIME: ", round(end.time - start.time, 1), " ", units(end.time - start.time),"\n"))
  
  
  if (!returnAsList) {
    return(unlist(res.l))
  }
  
  if (!is.null(listNames)) {
    names(res.l) = listNames
  }
  
  res.l
  
}


.printExecutionTime <- function(startTime) {
  
  endTime  <-  Sys.time()
  flog.info(paste0("Script finished sucessfully. Execution time: ", round(endTime - startTime, 1), " ", units(endTime - startTime)))
}


checkAndLogWarningsAndErrors <- function(object, checkResult, isWarning = FALSE) {
  
  assert(checkCharacter(checkResult, len = 1), checkLogical(checkResult))
  
  if (checkResult != TRUE) {

    objectPart = ""
    if (!is.null(object)) {
      objectname = deparse(substitute(object))
      objectPart = paste0("Assertion on variable \"", objectname, "\" failed: ")
    } 
    
    lastPartError   = "# An error occurred. See details above. If you think this is a bug, please contact us. #\n"
    hashesStrError = paste0(paste0(rep("#", nchar(lastPartError) - 1), collapse = ""), "\n")
    messageError    = paste0(objectPart, checkResult, "\n\n", hashesStrError, lastPartError, hashesStrError)
    
    lastPartWarning = "# This warning may or may not be ignored. Carefully check the diagnostic files and downstream results. #\n"
    hashesStrWarning = paste0(paste0(rep("#", nchar(lastPartWarning) - 1), collapse = ""), "\n")
    messageWarning  = paste0(objectPart, checkResult, "\n\n", hashesStrWarning, lastPartWarning, hashesStrWarning)
    
    
    
    if (isWarning) {
      flog.warn(messageWarning)
      warning(messageWarning)
    } else {
      flog.error(messageError)
      stop(messageError)
    }
  }
}



#############
# FUNCTIONS #
#############

# Code from Armando Reyes
heatmap.act.rep <- function(df.tf.peak.matrix, HOCOMOCO_mapping.df.exp, cor.m, par.l, median.cor.tfs, median.cor.tfs.non, act.rep.thres.l){
    
    missingGenes = which(!HOCOMOCO_mapping.df.exp$ENSEMBL %in% colnames(cor.m))
    if (length(missingGenes) > 0) {
        HOCOMOCO_mapping.df.exp = dplyr::filter(HOCOMOCO_mapping.df.exp, ENSEMBL %in% colnames(cor.m))
    }
    
    cor.r.pearson.m <- cor.m[,as.character(HOCOMOCO_mapping.df.exp$ENSEMBL)]
    
    stopifnot(identical(colnames(df.tf.peak.matrix), as.character(HOCOMOCO_mapping.df.exp$HOCOID)))
    stopifnot(identical(colnames(cor.r.pearson.m), as.character(HOCOMOCO_mapping.df.exp$ENSEMBL)))
    colnames(cor.r.pearson.m) <- HOCOMOCO_mapping.df.exp$HOCOID
    BREAKS = seq(-1,1,0.05)
    diffDensityMat = matrix(NA, nrow = ncol(cor.r.pearson.m), ncol = length(BREAKS) - 1)
    rownames(diffDensityMat) = HOCOMOCO_mapping.df.exp$HOCOID
    
    TF_Peak_all.m <- df.tf.peak.matrix
    TF_Peak.m <- TF_Peak_all.m
    
    for (i in 1:ncol(cor.r.pearson.m)) {
        TF = colnames(cor.r.pearson.m)[i]
        TF_name = TF #as.character(HOCOMOCO_mapping.df.exp$HOCOID[HOCOMOCO_mapping.df.exp$HOCOID==TF])
        ## for the background, use all peaks
        h_noMotif = hist(cor.r.pearson.m[,TF][TF_Peak_all.m[,TF] == 0], breaks = BREAKS, plot = FALSE)
        ## for the foreground use only peaks with less than min_mot_n different TF motifs
        h_Motif = hist(cor.r.pearson.m[,TF][TF_Peak.m[,TF] != 0], breaks = BREAKS, plot = FALSE)
        diff_density = h_Motif$density - h_noMotif$density
        diffDensityMat[rownames(diffDensityMat) == TF_name[1], ] <- diff_density
    }
    diffDensityMat = diffDensityMat[!is.na(diffDensityMat[,1]),]
    colnames(diffDensityMat) = signif(h_Motif$mids,1)
    quantile(diffDensityMat)
    
    ## check to what extent the number of TF motifs affects the density values
    n_min = ifelse(colSums(TF_Peak.m) < nrow(TF_Peak.m),colSums(TF_Peak.m), nrow(TF_Peak.m)-colSums(TF_Peak.m))
    names(n_min) = HOCOMOCO_mapping.df.exp$HOCOID#[match(names(n_min), as.character(tf2ensg$ENSEMBL))]
    n_min <- sapply(split(n_min,names(n_min)),sum)
    quantile(n_min)
    remove_smallN = which(n_min < 100)
    cor(n_min[-remove_smallN],rowMax(diffDensityMat)[-remove_smallN], method = 'pearson')
    
    factorClassificationPlot <- sort(median.cor.tfs, decreasing = TRUE)
    diffDensityMat_Plot = diffDensityMat[match(names(factorClassificationPlot), rownames(diffDensityMat)), ]
    diffDensityMat_Plot = diffDensityMat_Plot[!is.na(rownames(diffDensityMat_Plot)),]
    annotation_rowDF = data.frame(median_diff = factorClassificationPlot[match(rownames(diffDensityMat_Plot), names(factorClassificationPlot))])
    

    for (thresCur in names(act.rep.thres.l)) {
        thresCur.v = act.rep.thres.l[[thresCur]]
        
        labelMain = paste0(as.numeric(thresCur)*100, " / ", (1 - as.numeric(thresCur))*100, " % quantiles")
        
        colBreaks = unique(c((-1),
                             thresCur.v[1], 
                             thresCur.v[2],
                             1))
        

        anno_rowDF = data.frame(threshold = cut(annotation_rowDF$median_diff, breaks = colBreaks))
        rownames(anno_rowDF) = rownames(diffDensityMat_Plot)
        colors = c(par.l$colorCategories["repressor"],par.l$colorCategories["not-expressed"], par.l$colorCategories["activator"])
        names(colors) = levels(anno_rowDF$threshold)
        
        pheatmap(diffDensityMat_Plot, cluster_rows = FALSE, cluster_cols = FALSE,
                 fontsize_row = 1.25, scale = 'row' , fontsize_col = 10, fontsize = 8, labels_col = c(-1, -0.5, 0, 0.5, 1),
                 annotation_row = anno_rowDF,annotation_legend = FALSE,
                 annotation_colors = list(threshold = colors), legend = TRUE, annotation_names_row = FALSE, main = labelMain)
    }
    
    
} # end function


# Which transformation of the y values to do?
transform_yValues <- function(values, addPseudoCount = TRUE, nPermutations, onlyForZero = TRUE) {
    
    # Should only happen with the permutation-based approach
    zeros = which(values == 0)
    if (length(zeros) > 0 & addPseudoCount) {
        values[zeros] = 1 / nPermutations
    }
    
    -log10(values)
}

transform_yValues_caption <- function() {
    "-log10"
}

my.median = function(x) median(x, na.rm = TRUE)
my.mean   = function(x) mean(x, na.rm = TRUE)


checkDesignIntegrity <- function(snakemake, par.l, sampleData.df, useRNA = FALSE) {
    
    
    par.l$designFormulaVariableTypes = snakemake@config$par_general$designVariableTypes
    checkAndLogWarningsAndErrors(par.l$designFormulaVariableTypes, checkCharacter(par.l$designFormulaVariableTypes, len = 1, min.chars = 3))
    par.l$designFormulaVariableTypes = gsub(" ", "", par.l$designFormulaVariableTypes)
    components = strsplit(par.l$designFormulaVariableTypes, ",")[[1]]
    
    par.l$conditionComparison  = snakemake@config$par_general$conditionComparison
    checkAndLogWarningsAndErrors(par.l$conditionComparison, checkCharacter(par.l$conditionComparison, len = 1))
    
    if (useRNA) {
        
        # The design formula for RNA-Seq is different from the one we used before for ATAC-Seq
        # Either take the one that the user provided or, if he did not, use a general one with only the condition
        par.l$designFormulaRNA = snakemake@config$par_general$designContrastRNA
        if (is.null(par.l$designFormulaRNA)) {
            par.l$designFormulaRNA = "~conditionSummary"
            flog.warn(paste0("Could not find the parameter designContrastRNA in the configuration file. The default of \"~conditionSummary\" will be taken as formula. If you know about confounding variables, rerun this step and add the parameter (see the Documentation for details)"))
        }
        
        designFormula = as.formula(par.l$designFormulaRNA)
    } else {
        
        par.l$designFormula= snakemake@config$par_general$designContrast
        designFormula = as.formula(par.l$designFormula)
    }
    
    formulaVariables = attr(terms(designFormula), "term.labels")
    
    indexInteractionTerms = which(grepl(":", formulaVariables))
    if (length(indexInteractionTerms) > 0) {
        
        interactionTerms = formulaVariables[indexInteractionTerms]
        # Check whether all individual items have been specified
        for (interactionTermCur in interactionTerms) {
            componentsCur = strsplit(interactionTermCur, ":")[[1]]
            if (!all(componentsCur %in% formulaVariables)) {
                message = paste0("Design formula is incorrect, not all terms from the interactions (", paste0(componentsCur, collapse = ","), ") have been specified individually.")
                checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
                
            }
            
        }
        formulaVariables = formulaVariables[-indexInteractionTerms]
        
        
    }
    
    checkAndLogWarningsAndErrors(components, checkVector(components, min.len = length(formulaVariables)))
    
    # Extract the variable that defines the contrast. Always the last element in the formula
    variableToPermute = formulaVariables[length(formulaVariables)]
    
    # Split further
    components2 = strsplit(components, ":")
    
    if (!all(sapply(components2,length) == 2)) {
        
        message = "The parameter \"designVariableTypes\" has not been specified correctly. It must contain all the variables that appear in the parameter \"designContrast\". See the documentation for details"
        checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        
    }
    
    components3 = unlist(lapply(components2, "[[", 1))
    components3types = tolower(unlist(lapply(components2, "[[", 2)))
    names(components3types) = components3
    checkAndLogWarningsAndErrors(formulaVariables, checkSubset(formulaVariables, components3))
    checkAndLogWarningsAndErrors(components3types, checkSubset(components3types, c("factor", "integer", "numeric", "logical")))
    
    # Check the sample table. Distinguish between the two modes: Factor and integer for conditionSummary
    if (components3types["conditionSummary"] == "logical" | components3types["conditionSummary"] == "factor") {
        
        nDistValues = length(unique(sampleData.df$conditionSummary))
        if (nDistValues != 2) {
            message = paste0("The column 'conditionSummary' must contain exactly 2 different values, but ", nDistValues, " were found.") 
            checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
        conditionsVec = strsplit(par.l$conditionComparison, ",")[[1]]
        if (!testSubset(as.character(unique(sampleData.df$conditionSummary)), conditionsVec)) {
            message = paste0("The specified elements for the parameter conditionComparison (", par.l$conditionComparison,   ") do not correspond to what is specified in the sample summary table (", paste0(unique(sampleData.df$conditionSummary), collapse = ","), "). All elements of conditionComparison must be present in the column conditionSummary.")
            checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
    } else {
        
        # Test whether at least two distinct values are present
        nDistinct = length(unique(sampleData.df$conditionSummary))
        if (nDistinct < 2) {
            message = paste0("At least 2 distinct values for the column 'conditionSummary' must be present in the sample file in the quantitative mode.") 
            checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
    }


    list(types = components3types, variableToPermute = variableToPermute)
    
}
