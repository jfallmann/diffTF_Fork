
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
         + geom_hline(aes(yintercept=0), col = "#9850C3", show.legend = FALSE)
         + geom_smooth(method = "loess", se = FALSE, col = "#5D84C5", span = .4)
         + theme_bw()
  )
  
  return(pl)
}


printMultipleGraphsPerPage <- function(plots.l, nCol = 1, nRow = 1, pdfFile = NULL, height = NULL, width = NULL, verbose = FALSE) {
  
  checkAndLoadPackages(c("grDevices", "gridExtra"), verbose = verbose)  
  
  assertList(plots.l, min.len = 1)
  assertInt(nCol, lower = 1)
  assertInt(nRow, lower = 1)
  assert(checkNull(pdfFile), checkCharacter(pdfFile))
  if (!testNull(pdfFile)) assertDirectory(dirname(pdfFile), access = "r")
  assert(checkNull(height), checkNumber(height, lower = 1))
  assert(checkNull(width), checkNumber(width, lower = 1))
  
  if (testNull(height)) {
    height = 7
  }
  
  if (testNull(width)) {
    width = 7
  }
  
  
  for (i in seq_len(length(plots.l))) {
    assertClass(plots.l[[i]], classes = c("ggplot", "gg"))
  }
  
  
  nPlotsPerPage =  nCol * nRow
  
  plotsNew.l = list()
  
  if (!testNull(pdfFile))  {
    clearOpenDevices()
    pdf(pdfFile, height = height, width = width)
  }
  
  index = 0
  
  for (indexAll in 1:length(plots.l)) {
    
    index = index + 1
    plotsNew.l[[index]] = plots.l[[indexAll]]
    
    # Print another page
    if (index %% nPlotsPerPage == 0) { ## print 8 plots on a page
      suppressMessages(print(do.call(grid.arrange,  c(plotsNew.l, list(ncol = nCol, nrow = nRow)))))
      plotsNew.l = list() # reset plot 
      index = 0 # reset index
    }
    
  }
  
  # Print the remainding plots in case they don't perfectly fit with the layout
  if (length(plotsNew.l) != 0) { 
    suppressMessages(print(do.call(grid.arrange,  c(plotsNew.l, list(ncol = nCol, nrow = nRow)))))
  }
  
  if (!testNull(pdfFile)) 
    dev.off()
  
  cat("Finished writing plots to file ", pdfFile, "\n")
}




computeDESeqDiagnosticPlots <- function(dd, filename = NULL, maxPairwiseComparisons = 5) {
  
  checkAndLoadPackages(c("tidyverse", "checkmate", "geneplotter", "DESeq2", "vsn", "RColorBrewer"), verbose = FALSE)
  
  
  assertClass(dd, "DESeqDataSet")
  assert(checkNull(filename), checkDirectory(dirname(filename), access = "w"))
  
  if (!is.null(filename)) {
    pdf(filename)
  }
  
  # 1. MA plots: An MA-plot (R. Dudoit et al. 2002) provides a useful overview for the distribution of the estimated coefficients in the model, e.g. the comparisons of interest, across all genes. On the y-axis, the “M” stands for “minus” – subtraction of log values is equivalent to the log of the ratio – and on the x-axis, the “A” stands for “average”. You may hear this plot also referred to as a mean-difference plot, or a Bland-Altman plot.
  
 #  Before making the MA-plot, we use the lfcShrink function to shrink the log2 fold changes for the comparison of dex treated vs untreated samples:
  
  # The log2 fold change for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor is shown on the x-axis. Each gene is represented with a dot. Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
  
  # The DESeq2 package uses a Bayesian procedure to moderate (or “shrink”) log2 fold changes from genes with very low counts and highly variable counts, as can be seen by the narrowing of the vertical spread of points on the left side of the MA-plot. As shown above, the lfcShrink function performs this operation. For a detailed explanation of the rationale of moderated fold changes, please see the DESeq2 paper (Love, Huber, and Anders 2014).
  
  # TODO: dd2 <- lfcShrink(dds, contrast=c("dex","trt","untrt"), res=res)
  # vary alpha from 0.01 to 0.2
  for (alphaCur in c(0.001, 0.01, 0.05, 0.1, 0.2)) {
    DESeq2::plotMA(dd, alpha = alphaCur, main = paste0("MA plot for alpha = ", alphaCur))
  }
  
  # Another useful diagnostic plot is the histogram of the p values (figure below). 
  # This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
  
  # TODO: hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
  
  # 2. Densities of counts for the different samples. 
  # Since most of the genes are (heavily) affected by the experimental conditions, a succesful normalization will lead to overlapping densities
  
  nSamples = nrow(colData(dd))
  colors = colorRampPalette(brewer.pal(9, "Set1"))(nSamples)
  
  xlabCur = "Mean log counts"
   
  if (nSamples < 10) {

    multidensity(log(counts(dd, normalized = FALSE) + 0.5), xlab = xlabCur, main = "Non-normalized log counts", col = colors)
    multidensity(log(counts(dd, normalized = TRUE) + 0.5) , xlab = xlabCur, main = "Normalized log counts", col = colors) 
    
    multiecdf(log(counts(dd, normalized = FALSE) + 0.5), xlab = xlabCur, main = "Non-normalized log counts", col = colors)
    multiecdf(log(counts(dd, normalized = TRUE) + 0.5) , xlab = xlabCur, main = "Normalized log counts", col = colors) 
    
  } else {
    
    warning("Omitting legend due to large number of samples (threshold is 10 at the moment)")
    
    multidensity(log(counts(dd, normalized = FALSE) + 0.5), xlab = xlabCur, main = "Non-normalized log counts", legend = NULL, col = colors)
    multidensity(log(counts(dd, normalized = TRUE) + 0.5) , xlab = xlabCur, main = "Normalized log counts",     legend = NULL, col = colors)
    
    multiecdf(log(counts(dd, normalized = FALSE) + 0.5), xlab = xlabCur, main = "Non-normalized log counts", legend = NULL, col = colors)
    multiecdf(log(counts(dd, normalized = TRUE) + 0.5) , xlab = xlabCur, main = "Normalized log counts",     legend = NULL, col = colors)
  }
  
  
  # 3. Pairwise sample comparisons.
  # To further assess systematic differences between the samples, we can also plot pairwise mean–average plots: We plot the average of the log–transformed counts vs the fold change per gene for each of the sample pairs.
  MA.idx = t(combn(seq_len(dim(colData(dd))[1]), 2))
  

  if (nrow(MA.idx) > maxPairwiseComparisons) {
    
    warning("The number of pairwise comparisons to plot exceeds the current maximum of ", maxPairwiseComparisons, ". Only ", maxPairwiseComparisons, " pairwise comparisons will be shown in the PDF.")
    MA.idx.filt = MA.idx[1:maxPairwiseComparisons,, drop = FALSE]

  } else {
    MA.idx.filt = MA.idx
  }
  
  for (i in seq_along(MA.idx.filt[,1])) { 
    
    
    label = paste0(colnames(dd)[MA.idx.filt[i,1]], " vs ", colnames(dd)[MA.idx.filt[i,2]])
    print(myMAPlot(counts(dd, normalized = TRUE), c(MA.idx[i,1], MA.idx.filt[i,2]), main =  label))
  }
  
  if (nrow(MA.idx) > nrow(MA.idx.filt)) {
    
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    message = paste0("All remaining pairwise comparisons plots\nbetween samples have been omitted\nfor time and memory reasons.\nThe current maximum is set to ", maxPairwiseComparisons, ".")
    text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
  }
  
  # Show an empty page with a warning if plots have been omitted
  
  # 4. Mean SD plot: Plot row standard deviations versus row means
  notAllZeroPeaks <- (rowSums(counts(dd)) > 0)
  
  meanSdPlot(assay(dd[notAllZeroPeaks,]))
  
  if (!is.null(filename)) {
    dev.off()
  }
}

createDebugFile <- function(snakemake) {
    
    checkAndLoadPackages(c("checkmate", "tools"), verbose = FALSE)
    
    if (!testClass(snakemake, "Snakemake")) {
        warning("Could not find snakemake object, therefore not saving anyting.")
        return()
    }
    
    logfile = snakemake@log[[1]][1]
    
    assertCharacter(logfile, any.missing = FALSE)
    
    filename = paste0(tools::file_path_sans_ext(logfile), ".rds")
    
    cat(paste0("Saved snakemake object in file ", filename, "\n"))
    
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
  
  if (verbose) message(" Finished execution using ",nCores," cores. TOTAL RUNNING TIME: ", round(end.time - start.time, 1), " ", units(end.time - start.time),"\n")
  
  
  if (!returnAsList) {
    return(unlist(res.l))
  }
  
  if (!is.null(listNames)) {
    names(res.l) = listNames
  }
  
  res.l
  
}


.printExecutionTime <- function(startTime, verbose = TRUE) {
  
  endTime  <-  Sys.time()
  if (verbose) message(" Execution time: ", round(endTime - startTime, 1), " ", units(endTime - startTime))
}


checkAndLogWarningsAndErrors <- function(object, checkResult, isWarning = FALSE) {
  
  assert(checkCharacter(checkResult, len = 1), checkLogical(checkResult))
  
  if (checkResult != TRUE) {

    objectPart = ""
    if (!is.null(object)) {
      objectname = deparse(substitute(object))
      objectPart = paste0("Assertion on variable \"", objectname, "\" failed: ")
    } 
    
    lastPartError   = "# This is an unexpected error and should not happen. Please contact us with error details. You may also run the R script manually and troubleshoot. #\n"
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

