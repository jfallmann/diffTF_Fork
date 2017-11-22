start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################

library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/2.DESeqPeaks.R.rds")
createDebugFile(snakemake)

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "vsn", "csaw", "checkmate", "limma", "tools", "EDASeq", "geneplotter", "RColorBrewer", "BiocParallel", "rlist"), verbose = FALSE)



###################
#### PARAMETERS ###
###################

par.l = list()

par.l$verbose = TRUE
par.l$log_minlevel = "INFO"

#####################
# VERIFY PARAMETERS #
#####################

checkAndLogWarningsAndErrors(snakemake, checkClass(snakemake, "Snakemake"))

## INPUT ##
checkAndLogWarningsAndErrors(snakemake@input, checkList(snakemake@input, min.len = 1))
checkAndLogWarningsAndErrors(snakemake@input, checkSubset(names(snakemake@input), c("", "sampleData", "peaks")))

par.l$file_input_sampleData = snakemake@input$sampleData
checkAndLogWarningsAndErrors(par.l$file_input_sampleData, checkFileExists(par.l$file_input_sampleData, access = "r"))

  
par.l$files_input_peaks = snakemake@input$peaks

for (fileCur in par.l$files_input_TF_summary) {
  checkAndLogWarningsAndErrors(fileCur, checkFileExists(fileCur, access = "r"))
}

## OUTPUT ##
checkAndLogWarningsAndErrors(snakemake@output, checkList(snakemake@output, min.len = 1))
checkAndLogWarningsAndErrors(names(snakemake@output), checkSubset(names(snakemake@output), c("", "sampleDataR", "peakFile", "peaks_tsv", "condComp", "normFacs", "plots", "plotsPerm", "DESeqObj")))

par.l$file_output_metadata  = snakemake@output$sampleDataR
par.l$file_output_peaks     = snakemake@output$peakFile
par.l$file_output_peaksTSV  = snakemake@output$peaks_tsv
par.l$file_output_condComp  = snakemake@output$condComp  
par.l$file_output_normFacs  = snakemake@output$normFacs
par.l$file_output_plots     = snakemake@output$plots
par.l$file_output_DESeq     = snakemake@output$DESeqObj

## CONFIG ##
checkAndLogWarningsAndErrors(snakemake@config, checkList(snakemake@config, min.len = 1))

par.l$designFormula = snakemake@config$par_general$designContrast
checkAndLogWarningsAndErrors(par.l$designFormula, checkCharacter(par.l$designFormula, len = 1, min.chars = 3))

par.l$designFormulaVariableTypes = snakemake@config$par_general$designVariableTypes
checkAndLogWarningsAndErrors(par.l$designFormulaVariableTypes, checkCharacter(par.l$designFormulaVariableTypes, len = 1, min.chars = 3))

par.l$nPermutations = snakemake@config$par_general$nPermutations
checkAndLogWarningsAndErrors(par.l$nPermutations, checkIntegerish(par.l$nPermutations, lower = 0))


## PARAMS ##
checkAndLogWarningsAndErrors(snakemake@params, checkList(snakemake@params, min.len = 1))
checkAndLogWarningsAndErrors(names(snakemake@params), checkSubset(names(snakemake@params), c("", "doCyclicLoess")))

par.l$doCyclicLoess = as.logical(snakemake@params$doCyclicLoess)
checkAndLogWarningsAndErrors(par.l$doCyclicLoess, checkFlag(par.l$doCyclicLoess))

## LOG ##
checkAndLogWarningsAndErrors(snakemake@log, checkList(snakemake@log, min.len = 1))
par.l$file_log = snakemake@log[[1]]


allDirs = c(dirname(par.l$file_output_metadata), 
            dirname(par.l$file_output_peaks), 
            dirname(par.l$file_output_normFacs), 
            dirname(par.l$file_output_plots), 
            dirname(par.l$file_output_peaksTSV),
            dirname(par.l$file_output_DESeq),
            dirname(par.l$file_log)
            )

testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
startLogger(par.l$file_log, par.l$log_minlevel,  removeOldLog = TRUE)
printParametersLog(par.l)



#################
# READ METADATA #
#################

sampleData.df = read_tsv(par.l$file_input_sampleData, col_names = TRUE, col_types = cols())

checkAndLogWarningsAndErrors(colnames(sampleData.df), checkSubset(c("bamReads"), colnames(sampleData.df)))

designFormula = as.formula(par.l$designFormula)
formulaVariables = attr(terms(designFormula), "term.labels")

par.l$designFormulaVariableTypes = gsub(" ", "", par.l$designFormulaVariableTypes)
components = strsplit(par.l$designFormulaVariableTypes, ",")[[1]]
checkAndLogWarningsAndErrors(components, checkVector(components, len = length(formulaVariables)))
# Split further
components2 = strsplit(components, ":")

if (!all(sapply(components2,length) == 2)) {
  
  message = "The parameter \"designVariableTypes\" has not been specified correctly. It must contain all the variables that appear in the parameter \"designContrast\". See the documentation for details"
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)

}

components3 = unlist(lapply(components2, "[[", 1))
variableToPermute = components3[length(components3)]
components3types = tolower(unlist(lapply(components2, "[[", 2)))
names(components3types) = components3
checkAndLogWarningsAndErrors(sort(formulaVariables), checkSetEqual(sort(formulaVariables), sort(components3)))
checkAndLogWarningsAndErrors(components3types, checkSubset(components3types, c("factor", "integer", "numeric", "logical")))

datatypeVariableToPermute = components3types[variableToPermute]

# Read and modify samples metadata
sampleData.df = mutate(sampleData.df, name = file_path_sans_ext(basename(sampleData.df$bamReads)))
  

# Check and change column types as specified in the design formula
for (colnameCur in names(components3types)) {
  
  coltype = components3types[colnameCur]
  if (coltype == "factor") {
    sampleData.df[,colnameCur] = as.factor(unlist(sampleData.df[,colnameCur]))
  } else if (coltype == "numeric") {
    sampleData.df[,colnameCur] = as.numeric(unlist(sampleData.df[,colnameCur]))
  } else if (coltype == "integer") {
    sampleData.df[,colnameCur] = as.integer(unlist(sampleData.df[,colnameCur]))
  } else if (coltype == "logical") {
    sampleData.df[,colnameCur] = as.logical(unlist(sampleData.df[,colnameCur]))
  } else {
    message = paste0("Unknown type: ", colnameCur)
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
}

# If variable to permute is a factor, check that is has 2 levels 
nLevels = length(unique(unlist(sampleData.df[,variableToPermute])))
if (datatypeVariableToPermute == "factor" & nLevels != 2) {
  message = paste0("The variable ", variableToPermute, " was specified as a factor, but it does not have two different levels but instead ", nLevels, ".")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}





##############################
# ITERATE THROUGH PEAK FILES #
##############################

peaks.df = NULL
coverageAll.m = NULL

# TODO: Peak annotation verifiation: Really unique? Number of columns correct? Must be between 3 and 6


flog.info(paste0("Iterating over ", length(par.l$files_input_peaks), " peak files "))

for (fileCur in par.l$files_input_peaks) {

  # TODO: Change accordingly: Between 4 and 7 here, set column names accordingly and
  peaks.df =  read_tsv(fileCur, col_names = c("chr", "PSS", "PES", "annotation", "ID", "coverage"), col_types = cols())
  flog.info(paste0("Parsed peak file ", fileCur, " with ", nrow(peaks.df)," rows"))
  
  if (nrow(peaks.df) == 0) {
    
    message = paste0("The file ", fileCur, " is empty, no overlaps have been identified. This file will be skipped.")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    
  } else {
    
    #TODO: Needed?
    peaks.df$identifier = paste0(peaks.df$chr,":", peaks.df$PSS,"-", peaks.df$PES)
    
    # Filter and retain only unique identifiers
    peaks.filtered.df = distinct(peaks.df, identifier, .keep_all = TRUE)
    
    nRowsFiltered = nrow(peaks.df) - nrow(peaks.filtered.df)
    if (par.l$verbose & nRowsFiltered  > 0) flog.info(paste0("Filtered ", nRowsFiltered, " non-unique positions out of ", nrow(peaks.df), " from peaks table."))
    
    peaks.df = peaks.filtered.df
    
    # concatenate results from COV from each iteration
    coverageAll.m = cbind(coverageAll.m, peaks.df$coverage)
    
  }
 
}


if (is.null(coverageAll.m)) {
  
  message = paste0("All overlap files are empty. Cannot continue. A different peak file may solve the issue.")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
} 


# Save the last, they are all identical anyway except for the count column
saveRDS(peaks.filtered.df, file = par.l$file_output_peaks)


## transform as matrix data frame with counts
coverageAll.m = as.matrix(coverageAll.m)
colnames(coverageAll.m) = sampleData.df$SampleID
rownames(coverageAll.m) = peaks.df$identifier # Take the first element as nameCur reprsentative, they are all identical anyway


#############################
# SAMPLE LABEL PERMUTATIONS #
#############################

sampleDataOrig.df = sampleData.df

samplesRare.l = NA


conditionCounter = table(sampleDataOrig.df[,variableToPermute])

nSamplesRareCondition     = min(conditionCounter)
nSamplesFrequentCondition = max(conditionCounter)

# Always take only the first element of the vectors, as conditions can have equal lengths
nameRareCondition      = names(conditionCounter)[conditionCounter == min(conditionCounter)][1]
nameFrequentCondition  = names(conditionCounter)[which(names(conditionCounter) != nameRareCondition)]
indexRareCondition     = which(sampleDataOrig.df[,variableToPermute] == nameRareCondition)
indexFrequentCondition = which(sampleDataOrig.df[,variableToPermute] == nameFrequentCondition)

nSamples = nrow(sampleData.df)

# Determine the number of permutations

if (datatypeVariableToPermute == "factor") {
  
  nSwaps = ceiling(nSamplesRareCondition / 2)
  nPermutationsTotal = choose(nSamplesRareCondition,nSwaps ) * choose(nSamplesFrequentCondition,nSwaps)

} else {
  
  nPermutationsTotal = factorial(nSamples)
  
}

if (nPermutationsTotal < par.l$nPermutations) {
  
  message = paste0("The total number of possible permutations is only ", nPermutationsTotal, ", but more have been requested. The value for the paramter nPermutations will be adjusted.")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  par.l$nPermutations = nPermutationsTotal
}


# Permute samples beforehand here so that each call to a permutation is unique
permutationsList.l = list()
samplesRare.l = list()
nPermutationsDone = 0
failsafeCounter   = 0
while (nPermutationsDone < (par.l$nPermutations + 1)) {

  #permutation 0 is always the original, non-permuted data
  if (nPermutationsDone == 0) {
    
    if (datatypeVariableToPermute == "factor") {
      indexRareConditionCur = indexRareCondition
    } else {
      
      # don't draw a sample just take the original vector
      indexRareConditionCur = unlist(sampleDataOrig.df[,variableToPermute])
    }

    
  } else {
    
    if (datatypeVariableToPermute == "factor") {
      
      indexRareConditionCur     = indexRareCondition
      indexFrequentConditionCur = indexFrequentCondition
      
      # Pick around 50% of the number of cases for the rare condition and swap
      sampleRareConditionSwap     = sort(sample(seq_len(nSamplesRareCondition), nSwaps))
      sampleFrequentConditionSwap = sort(sample(seq_len(nSamplesFrequentCondition), nSwaps))
      
      # Swap indexes
      temp = indexRareConditionCur[sampleRareConditionSwap]
      indexRareConditionCur[sampleRareConditionSwap] = indexFrequentConditionCur[sampleFrequentConditionSwap]
      indexFrequentConditionCur[sampleFrequentConditionSwap] = temp
      
      # Sort
      indexRareConditionCur     = sort(indexRareConditionCur)
      indexFrequentConditionCur = sort(indexFrequentConditionCur)
      
      stopifnot(length(unique(c(indexRareConditionCur, indexFrequentConditionCur))) == nSamples)
      
    } else {
      
      # Just draw a sample from the full vector
      indexRareConditionCur = sample(seq_len(nSamples), nSamples)
    }
    

  }

  
  # Check if this permutation has already been used. If yes, produce a different one
  permutationNameStr = paste0(indexRareConditionCur, collapse = ",")
  if (!permutationNameStr %in% names(samplesRare.l)) {
    failsafeCounter   = 0
    samplesRare.l[[permutationNameStr]] = indexRareConditionCur
    nPermutationsDone = nPermutationsDone + 1
  } else {
    failsafeCounter   =  failsafeCounter + 1
    if (failsafeCounter > 1000) {
      message = "Could not generate more permutations. This looks like a bug."
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    } 
  }
  
}

  
#Rename so it is easier to address in the following code
listNames = paste0("permutation", 0:par.l$nPermutations)
names(samplesRare.l) = listNames


runPermutation <- function(permutationCur, sampleDataOrig.df, variableToPermute, datatypeVariableToPermute, nameRareCondition, nameFrequentCondition, samplesRare.l, 
                           par.l, coverageAll.m) {
  
  sampleData.df = sampleDataOrig.df
  nSamples = nrow(sampleData.df)
  permutationName = paste0("permutation", permutationCur)
  
  ###################
  # PERMUTE SAMPLES #
  ###################
  
  if (permutationCur > 0) {
    flog.info(paste0("Running for permutation ", permutationCur))
    
    if (datatypeVariableToPermute == "factor") {
      
      rareConditionSamplesCur = samplesRare.l[[permutationName]]
      sampleData.df[rareConditionSamplesCur                            ,variableToPermute] = nameRareCondition
      sampleData.df[setdiff(seq_len(nSamples), rareConditionSamplesCur),variableToPermute] = nameFrequentCondition 
      
      
    } else if (datatypeVariableToPermute == "integer") {
    
      sampleData.df[,variableToPermute] = samplesRare.l[[permutationName]]
      
    }
   
    
  } else {
    flog.info(paste0("Running for original data"))
  }
 
  

  #############
  # RUN DESEQ #
  #############
  
  designFormula = convertToFormula(par.l$designFormula, colnames(sampleData.df))
  
  cds.peaks <- DESeqDataSetFromMatrix(countData = coverageAll.m,
                                      colData = sampleData.df,
                                      design = designFormula)
  
  
  
  # Do a regular size factor normalization
  if (!par.l$doCyclicLoess) {
    
    cds.peaks <- estimateSizeFactors(cds.peaks)
    
    normFacs = sizeFactors(cds.peaks)
    
  } else {
    
    # Perform a cyclic loess normalization
    # We use a slighlty more complicated setup to derive size factors for library normalization
    # Instead of just determining the size factors in DeSeq2 via cirtual samples, we use 
    # a normalization from the csaw package (see https://www.rdocumentation.org/packages/csaw/versions/1.6.1/topics/normOffsets)
    # and apply a non-linear normalization. 
    # For each sample, a lowess curve is fitted to the log-counts against the log-average count. 
    # The fitted value for each bin pair is used as the generalized linear model offset for that sample. 
    # The use of the average count provides more stability than the average log-count when low counts are present for differentially bound regions.
    
    
    # since counts returns,by default, non-normalized counts, the following code should be fine and there is no need to
    # also run estimateSizeFactors beforehand
    
    normFacs <- exp(normOffsets(counts(cds.peaks),
                                lib.sizes = colSums(counts(cds.peaks)),
                                type = "loess"))
    
    # sanity check: is the geometric mean across samples equal to one?
    
    #library("psych")
    #all.equal(geometric.mean(t(normFacs)), rep(1, dim(cds.peaks)[1]))
    
    rownames(normFacs) = rownames(coverageAll.m)
    colnames(normFacs) = colnames(coverageAll.m)
    
    # We now provide gene-specific normalization factors for each sample as a matrix, which will preempt sizeFactors
    normalizationFactors(cds.peaks) <- normFacs
    
  }
  
  compareNormalizations = FALSE
  
  if (compareNormalizations) {
    
    # Compare:
    # 2. Because size factors in DeSeq2 have a geometric mean of 1 and from normOffsets not, there needs to be another normalization with another geom_mean for the latter
    # 1. Geometic mean of all columns from normOffsets
    
    # dds <- makeExampleDESeqDataSet(n=1000, m=4)
    # TODO: An MA-plot here would be a good idea to test beforehand if a double normalization is needed
    
    # cds.peaks = dds
    sizeFactors1 = apply(normFacs, 2, geometric.mean)
    sizeFactors2 = sizeFactors1 / geometric.mean(sizeFactors1)
    
    # Compare with DeSeq2 size factors
    cds.peaks <- estimateSizeFactors(cds.peaks)
    sizeFactors(cds.peaks)
    sizeFactors(cds.peaks) <- NULL
  }
  
  
  # Filter peaks with zero counts
  cds.peaks.filt = cds.peaks[rowMeans(counts(cds.peaks)) > 0, ]
  
  
  # Deseq analysis
  cds.peaks.filt = tryCatch( {
    DESeq(cds.peaks.filt, fitType = 'local', quiet = TRUE)
    
  }, error = function(e) {
    message = "Warning: Could not run DESeq with local fitting, retry with default fitting type..."
    checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    DESeq(cds.peaks.filt, quiet = TRUE)
  }
  )
  
  cds.peaks.df <- as.data.frame(DESeq2::results(cds.peaks.filt))
  
  final.peaks.df = data_frame( 
                           "permutation" = permutationCur,
                           "position"    = rownames(cds.peaks.df), 
                           "D2_baseMean" = cds.peaks.df$baseMean,
                           "D2_l2FC"     = cds.peaks.df$log2FoldChange,
                           "D2_ldcSE"    = cds.peaks.df$lfcSE,
                           "D2_stat"     = cds.peaks.df$stat,
                           "D2_pval"     =  cds.peaks.df$pvalue, 
                           "D2_padj"     =  cds.peaks.df$padj
                  )
  
  
  ##################
  # PLOTS AND SAVE #
  ##################
  
  if (permutationCur == 0) {
    

    # Save the comparison that DeSeq made for later scripts
    comparisonDESeq = getComparisonFromDeSeqObject(cds.peaks.filt, par.l$designFormula, datatypeVariableToPermute)
    saveRDS(comparisonDESeq, file = par.l$file_output_condComp)
    saveRDS(cds.peaks.filt, file = par.l$file_output_DESeq)
    
    computeDESeqDiagnosticPlots(cds.peaks.filt, par.l$file_output_plots)
    # pdf(par.l$file_output_plots)
    # dev.off()
  
  } else {
   
    filenameCur = paste0(file_path_sans_ext(par.l$file_output_plots), "_permutation", permutationCur, ".pdf")
    computeDESeqDiagnosticPlots(cds.peaks.filt, filenameCur, maxPairwiseComparisons = 5)
    # pdf(filenameCur)
    # dev.off()
  }
  
  if (permutationCur > 0) {
    flog.info(paste0("Finished permutation ", permutationCur))
  } else {
    flog.info(paste0("Finished original data"))
  }
  
  
  return(list(peaks = final.peaks.df, normFacs = normFacs, sampleData = sampleData.df))
  
  
} # end of for each permutation

nCores = snakemake@threads
results.l = .execInParallelGen(nCores, 
                               returnAsList = TRUE, listNames = listNames, 
                               iteration = 0:par.l$nPermutations, 
                               abortIfErrorParallel = TRUE, 
                               verbose = TRUE,
                               runPermutation, sampleDataOrig.df, variableToPermute, datatypeVariableToPermute, nameRareCondition, nameFrequentCondition, samplesRare.l, 
                               par.l, coverageAll.m)

final.peaks.df = bind_rows(list.map(results.l, peaks))
sampleData.l   = list.map(results.l, sampleData)
normFacts.l    = list.map(results.l, normFacs)


################
# WRITE OUTPUT #
################

write_tsv(final.peaks.df, path = par.l$file_output_peaksTSV)

saveRDS(sampleData.l, par.l$file_output_metadata)
saveRDS(normFacts.l, par.l$file_output_normFacs)


.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)