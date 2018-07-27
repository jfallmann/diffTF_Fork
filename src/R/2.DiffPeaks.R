start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/2.DiffPeaks.R.rds")

library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################


createDebugFile(snakemake)

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "csaw", "checkmate", "limma", "tools"), verbose = FALSE)



###################
#### PARAMETERS ###
###################

par.l = list()

par.l$pseudocountAddition = 1
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"

#####################
# VERIFY PARAMETERS #
#####################

checkAndLogWarningsAndErrors(snakemake, checkmate::checkClass(snakemake, "Snakemake"))

## INPUT ##
checkAndLogWarningsAndErrors(snakemake@input, checkmate::checkList(snakemake@input, min.len = 1))
checkAndLogWarningsAndErrors(snakemake@input, checkmate::checkSubset(names(snakemake@input), c("", "sampleData", "BAMPeakoverlaps")))

par.l$file_input_sampleData = snakemake@input$sampleData
checkAndLogWarningsAndErrors(par.l$file_input_sampleData, checkmate::checkFileExists(par.l$file_input_sampleData, access = "r"))

  
par.l$file_input_peakOverlaps = snakemake@input$BAMPeakoverlaps

for (fileCur in par.l$files_input_TF_summary) {
  checkAndLogWarningsAndErrors(fileCur, checkmate::checkFileExists(fileCur, access = "r"))
}

## OUTPUT ##
checkAndLogWarningsAndErrors(snakemake@output, checkmate::checkList(snakemake@output, min.len = 1))
checkAndLogWarningsAndErrors(names(snakemake@output), checkmate::checkSubset(c("", "sampleDataR", "peakFile", "peaks_tsv", "condComp", "normFacs", "normCounts", "plots", "DESeqObj"), names(snakemake@output)))

par.l$file_output_metadata     = snakemake@output$sampleDataR
par.l$file_output_peaks        = snakemake@output$peakFile
par.l$file_output_peaksTSV     = snakemake@output$peaks_tsv
# par.l$file_output_peaksPermTSV = snakemake@output$peaksPerm_tsv
par.l$file_output_condComp     = snakemake@output$condComp  
par.l$file_output_normFacs     = snakemake@output$normFacs
par.l$file_output_normCounts   = snakemake@output$normCounts
par.l$file_output_plots        = snakemake@output$plots
par.l$file_output_DESeqObj     = snakemake@output$DESeqObj


## CONFIG ##
checkAndLogWarningsAndErrors(snakemake@config,checkmate::checkList(snakemake@config, min.len = 1))

par.l$designFormula = snakemake@config$par_general$designContrast
checkAndLogWarningsAndErrors(par.l$designFormula, checkCharacter(par.l$designFormula, len = 1, min.chars = 3))

par.l$designFormulaVariableTypes = snakemake@config$par_general$designVariableTypes
checkAndLogWarningsAndErrors(par.l$designFormulaVariableTypes, checkCharacter(par.l$designFormulaVariableTypes, len = 1, min.chars = 3))

par.l$nPermutations = snakemake@config$par_general$nPermutations
checkAndLogWarningsAndErrors(par.l$nPermutations, checkIntegerish(par.l$nPermutations, lower = 0))

par.l$conditionComparison  = snakemake@config$par_general$conditionComparison
checkAndLogWarningsAndErrors(par.l$conditionComparison, checkCharacter(par.l$conditionComparison, len = 1, min.chars = 3))

## PARAMS ##
checkAndLogWarningsAndErrors(snakemake@params,checkmate::checkList(snakemake@params, min.len = 1))
checkAndLogWarningsAndErrors(names(snakemake@params), checkSubset(names(snakemake@params), c("", "doCyclicLoess")))

par.l$doCyclicLoess = as.logical(snakemake@params$doCyclicLoess)
checkAndLogWarningsAndErrors(par.l$doCyclicLoess, checkFlag(par.l$doCyclicLoess))

## LOG ##
checkAndLogWarningsAndErrors(snakemake@log,checkmate::checkList(snakemake@log, min.len = 1))
par.l$file_log = snakemake@log[[1]]


allDirs = c(dirname(par.l$file_output_metadata), 
            dirname(par.l$file_output_peaks), 
            dirname(par.l$file_output_normFacs), 
            dirname(par.l$file_output_peaksTSV),
            dirname(par.l$file_log)
            )

testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel,  removeOldLog = TRUE)
printParametersLog(par.l)


#################
# READ METADATA #
#################

sampleData.df = read_tsv(par.l$file_input_sampleData, col_names = TRUE, col_types = cols())

checkAndLogWarningsAndErrors(colnames(sampleData.df), checkSubset(c("bamReads"), colnames(sampleData.df)))

conditionsVec = strsplit(par.l$conditionComparison, ",")[[1]]
if (!testSubset(sampleData.df$conditionSummary, conditionsVec)) {
  message = paste0("The specified elements for the parameter conditionComparison (", par.l$conditionComparison,   ") do not correspond to what is specified in the sample summary table (", paste0(unique(sampleData.df$conditionSummary), collapse = ","), "). All elements of conditionComparison must be present in the column conditionSummary.")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}

designFormula = as.formula(par.l$designFormula)
formulaVariables = attr(terms(designFormula), "term.labels")

# Extract the variable that defines the contrast. Always the last element in the formula
variableToPermute = formulaVariables[length(formulaVariables)]

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

# Change the conditionSummary specifically and enforce the direction as specified in the config file
sampleData.df$conditionSummary = factor(sampleData.df$conditionSummary, levels = conditionsVec)


# If variable to permute is a factor, check that is has 2 levels 
nLevels = length(unique(unlist(sampleData.df[,variableToPermute])))
if (datatypeVariableToPermute == "factor" & nLevels != 2) {
  message = paste0("The variable ", variableToPermute, " was specified as a factor, but it does not have two different levels but instead ", nLevels, ".")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}





##############################
# ITERATE THROUGH PEAK FILES #
##############################

coverageAll.df = read_tsv(par.l$file_input_peakOverlaps, col_names = TRUE, comment = "#", col_types = cols())

if (nrow(problems(coverageAll.df)) > 0) {
    flog.fatal(paste0("Parsing errors: "), problems(coverageAll.df), capture = TRUE)
    stop("Error when parsing the file ", fileCur, ", see errors above")
}

if (nrow(coverageAll.df) == 0) {
    
    message = paste0("Empty file ", par.l$file_input_peaks, ".")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
} 

## transform as matrix data frame with counts
coverageAll.m = as.matrix(dplyr::select(coverageAll.df, -one_of("Geneid", "Chr", "Start", "End", "Strand", "Length")))

# Take the basenames of the files, which have not been modified in the sorted versions of the BAMs as they are just located in different folders
sampleIDs = sampleData.df$SampleID[which(basename(sampleData.df$bamReads) %in% basename(colnames(coverageAll.m)))]

if (length(unique(sampleIDs)) != nrow(sampleData.df)) {
    message = paste0("Colnames mismatch.")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
} 

if (length(sampleIDs) != ncol(coverageAll.m)) {
  message = paste0("Mismatch between number of sample IDs (", length(sampleIDs), ") and number of columns in coverage file (", ncol(coverageAll.m), "). It appears that the number of samples been changed after running the pipeline the first time. Rerun the full pipeline from scratch.")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
} 

colnames(coverageAll.m) = sampleIDs

rownames(coverageAll.m) = coverageAll.df$Geneid


peaks.df = dplyr::select(coverageAll.df, one_of(c("Chr", "Start", "End", "Geneid")))
colnames(peaks.df) = c("chr", "PSS", "PES", "annotation")
# Filter and retain only unique identifiers
peaks.filtered.df = distinct(peaks.df, annotation, .keep_all = TRUE)
nRowsFiltered = nrow(peaks.df) - nrow(peaks.filtered.df)
if (par.l$verbose & nRowsFiltered  > 0) flog.info(paste0("Filtered ", nRowsFiltered, " non-unique positions out of ", nrow(peaks.df), " from peaks table."))
peaks.df = peaks.filtered.df

# Save the last, they are all identical anyway except for the count column
saveRDS(peaks.filtered.df, file = par.l$file_output_peaks)


#############################
# SAMPLE LABEL PERMUTATIONS #
#############################

sampleData.l = list()
nSamples = nrow(sampleData.df)
sampleData.l[["permutation0"]] = sampleData.df

conditionCounter = table(sampleData.df[,variableToPermute])

# Record the frequency of the conditions to determine how many permutations are possibler
nSamplesRareCondition     = min(conditionCounter)
nSamplesFrequentCondition = max(conditionCounter)
nameRareCondition         = names(conditionCounter)[conditionCounter == min(conditionCounter)][1]
nameFrequentCondition     = names(conditionCounter)[which(names(conditionCounter) != nameRareCondition)]
nPermutationsTotal        = choose(nSamples, nSamplesFrequentCondition) # same as choose(nSamples, nSamplesRareCondition)

if (nPermutationsTotal < par.l$nPermutations) {
  
  message = paste0("The total number of possible permutations is only ", nPermutationsTotal, ", but more have been requested. The value for the parameter nPermutations will be adjusted.")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  par.l$nPermutations = nPermutationsTotal
}


# Permute samples beforehand here so that each call to a permutation is unique
permutationsList.l = list()
nPermutationsDone = 0
failsafeCounter   = 0
while (nPermutationsDone < par.l$nPermutations) {

  sampleCur = sample.int(nSamples)
  samplesRareCondShuffled  = sampleData.df$SampleID[which(sampleData.df$conditionSummary[sampleCur] == nameRareCondition)]
  indexNameCur = paste0(sort(samplesRareCondShuffled), collapse = ",")
  
  # Check if this permutation has already been used. If yes, produce a different one

  if (!indexNameCur %in% names(permutationsList.l)) {
    failsafeCounter   = 0
    permutationsList.l[[indexNameCur]] = sampleCur
    nPermutationsDone = nPermutationsDone + 1
  } else {
    failsafeCounter   =  failsafeCounter + 1
    if (failsafeCounter > 5000) {
      message = "Could not generate more permutations. This looks like a bug."
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    } 
  }
  
}

##############################################
# RUN DESEQ TO OBTAIN NORMALIZED COUNTS ONLY #
##############################################

designFormula = convertToFormula(par.l$designFormula, colnames(sampleData.df))

cds.peaks <- DESeqDataSetFromMatrix(countData = coverageAll.m,
                                    colData = sampleData.df,
                                    design = designFormula)

# Recent versions of DeSeq seem to do this automatically, whereas older versions don't, so enforce it here
if (!identical(colnames(cds.peaks), colnames(coverageAll.m))) {
    colnames(cds.peaks) = colnames(coverageAll.m)
}
if (!identical(rownames(cds.peaks), rownames(coverageAll.m))) {
    rownames(cds.peaks) = rownames(coverageAll.m)
}

# Do a regular size factor normalization
if (!par.l$doCyclicLoess) {
    
    cds.peaks <- estimateSizeFactors(cds.peaks)
    
    normFacs = sizeFactors(cds.peaks)
    
    # TODO: NormFacts identical between different permutations? yes or?
    
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


# Filter peaks with zero counts
cds.peaks.filt   = cds.peaks[rowMeans(counts(cds.peaks)) > 0, ]


# DESeq log2fc are not used at all afterwards, as we currently only take the normalization factors to normalize the TFBS subsequently


# The levels have to be reversed because the first element is the one appearing at the right of the plot, with positive values as. compared to the reference
comparisonDESeq = rev(levels(sampleData.df$conditionSummary))

##############
# GET LOG2FC #
##############

countsNorm        = counts(cds.peaks.filt, norm = TRUE)
countsNorm.df     = as.data.frame(countsNorm) %>%
  dplyr::mutate(peakID = rownames(cds.peaks.filt))  %>%
  dplyr::select(one_of("peakID", colnames(countsNorm)))


if (par.l$nPermutations == 0) {
  
  # Generate normalized counts for limma analysis
  countsNorm.transf = log2(countsNorm + par.l$pseudocountAddition)
  rownames(countsNorm.transf) = rownames(cds.peaks.filt)
  
  sampleData.df$conditionSummary = factor(sampleData.df$conditionSummary)
  
  designMatrix = model.matrix(designFormula, data = sampleData.df)
  
  if (nrow(designMatrix) < nrow(sampleData.df)) {
    missingRows = setdiff(1:nrow(sampleData.df), as.integer(row.names(designMatrix)))
    message = paste0("There is a problem with the specified design formula (parameter designContrast): The corresponding design matrix has fewer rows. This usually means that there are missing values in one of the specified variables. The problem comes from the following lines in the summary file: ", paste0(missingRows, collapse = ","), ".") 
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  fit        <- eBayes(lmFit(countsNorm.transf, design = designMatrix))
  results.df <- topTable(fit, coef = colnames(fit$design)[ncol(fit$design)], number = Inf, sort.by = "none")
  
  final.peaks.df = data_frame(  
    "permutation" = 0,
    "peakID"      = rownames(results.df), 
    "limma_avgExpr"     = results.df$AveExpr,
    "l2FC"        = results.df$logFC,
    "limma_B"           = results.df$B,
    "limma_t_stat"      = results.df$t,
    "pval"        = results.df$P.Value, 
    "pval_adj"    = results.df$adj.P.Val
  )
  
  
  plotDiagnosticPlots(cds.peaks.filt, fit, comparisonDESeq, par.l$file_output_plots, maxPairwiseComparisons = 20)
  
  
} else {

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
  
  # TODO: "peakID"      = rownames(results.df), CHECK
  final.peaks.df = data_frame( 
    "permutation" = 0,
    "peakID"    = rownames(cds.peaks.df), 
    "DESeq_baseMean" = cds.peaks.df$baseMean,
    "l2FC"     = cds.peaks.df$log2FoldChange,
    "DESeq_ldcSE"    = cds.peaks.df$lfcSE,
    "DESeq_stat"     = cds.peaks.df$stat,
    "pval"     =  cds.peaks.df$pvalue, 
    "pval_adj" =  cds.peaks.df$padj
  )
  
  plotDiagnosticPlots(cds.peaks.filt, cds.peaks.filt, comparisonDESeq, par.l$file_output_plots, maxPairwiseComparisons = 20)
  
}


saveRDS(cds.peaks.filt, file = par.l$file_output_DESeqObj)

####################
# RUN PERMUTATIONS #
####################

if (par.l$nPermutations > 0) {
  #Rename so it is easier to address in the following code
  listNames = paste0("permutation", seq_len(par.l$nPermutations))
  names(permutationsList.l) = listNames
  
  # final.peaks.perm.df = tribble(~permutation, ~peakID, ~l2FC)
  sampleDataOrig.df = sampleData.df
  
  for (permutationCur in names(permutationsList.l)) {
    
    flog.info(paste0("Running for permutation ", permutationCur))
    
    sampleData.df = sampleDataOrig.df
    sampleData.df[,variableToPermute] = unlist(sampleData.df[,variableToPermute]) [permutationsList.l[[permutationCur]]]
    
    sampleData.l[[permutationCur]] = sampleData.df

    # We don't need permuted peaks l2fc actually so permutation-specific l2fc can be skipped
  
  }
}




################
# WRITE OUTPUT #
################

saveRDS(sampleData.l, par.l$file_output_metadata)

final.peaks.df = mutate_if(final.peaks.df, is.numeric, as.character)
write_tsv(final.peaks.df, path = par.l$file_output_peaksTSV)

#final.peaks.perm.df = mutate_if(final.peaks.perm.df, is.numeric, as.character)
#write_tsv(final.peaks.perm.df, path = par.l$file_output_peaksPermTSV)

saveRDS(comparisonDESeq, file = par.l$file_output_condComp)

saveRDS(normFacs, par.l$file_output_normFacs)

# Deactivated, because column names starting with numbers will cause problems and crash. Since the names come from the sample file, this cannot be excluded.
#countsNorm.df = mutate_if(countsNorm.df, is.numeric, as.character)
write_tsv(countsNorm.df, path = par.l$file_output_normCounts)



.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)