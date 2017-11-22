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
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/7.summaryFinal.R.rds")
createDebugFile(snakemake)

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "lsr", "dplyr", "ggrepel", "checkmate", "tools", "locfdr", "DESeq2"), verbose = FALSE)


###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$FDR_threshold = c(0.001, 0.01, 0.05,0.1,0.2)
par.l$probsThreshold = c(0.01, 0.99)
par.l$probsThreshold = c(0.05, 0.95)
par.l$expressionThreshold = 2
par.l$cohensDThreshold = 0.1
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.1, 0.5, 0.8)
par.l$log_minlevel = "INFO"
par.l$circularPlot_minDimensions  = 12
par.l$corMethod = "pearson"

par.l$extension_x_limits = 0.025 
par.l$extension_x_limits = 0.15 # 20 % x axis extension, regardless of the limits
par.l$extension_y_limits = 1
par.l$plot_grayColor = "grey50"

# diverging, modified
par.l$colorCategories = c("activator" = "#d7191c", "undetermined" = "black", "repressor" = "#2b83ba", "not-expressed" = "slategrey")
par.l$colorCategories = c("activator" = "#4daf4a", "undetermined" = "black", "repressor" = "#e41a1c", "not-expressed" = "Snow3")

par.l$colorConditions = c("#ef8a62", "#67a9cf")
par.l$circularPlot_height = 12
par.l$circularPlot_width = 16
par.l$size_TFAnnotation = 5.5
par.l$sizeLegend  = 20
par.l$rootFontSize = 8
par.l$sizeHelperLines = 0.3
par.l$legend_position = c(0.1, 0.9)

# Which transformation of the y values to do?
transform_yValues <- function(values) {
  log10(-log10(values))
}

#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(names(snakemake@input), c("", "allPermutationResults", "condComp", "DeSeqObj"))

par.l$files_input_permResults  = snakemake@input$allPermutationResults
for (fileCur in par.l$files_input_permResults) {
  assertFileExists(fileCur, access = "r")
}

par.l$file_input_condCompDeSeq = snakemake@input$condComp
assertFileExists(par.l$file_input_condCompDeSeq, access = "r")

par.l$file_input_DESeq = snakemake@input$DeSeqObj
assertFileExists(par.l$file_input_DESeq, access = "r")

## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "summary", "circularPlot", "diagnosticPlots"))

par.l$file_output_summary = snakemake@output$summary
par.l$file_plotCircular   = snakemake@output$circularPlot
par.l$file_plotDiagnostic = snakemake@output$diagnosticPlots


## CONFIG ##
assertList(snakemake@config, min.len = 1)

par.l$plotRNASeqClassification = as.logical(snakemake@config$par_general$RNASeqIntegration)
assertFlag(par.l$plotRNASeqClassification)



par.l$nPermutations = snakemake@config$par_general$nPermutations
assertIntegerish(par.l$nPermutations, lower = 0)

par.l$outdir = snakemake@config$par_general$outdir


if (par.l$plotRNASeqClassification) {

  par.l$file_input_HOCOMOCO_mapping    = snakemake@config$additionalInputFiles$HOCOMOCO_mapping
  par.l$file_input_geneCountsPerSample = snakemake@config$additionalInputFiles$RNASeqCounts
  assertFileExists(par.l$file_input_HOCOMOCO_mapping, access = "r")
  assertFileExists(par.l$file_input_geneCountsPerSample, access = "r")
}




## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]


allDirs = c(dirname(par.l$file_output_summary), 
            dirname(par.l$file_plotCircular),
            dirname(par.l$file_plotDiagnostic),
            dirname(par.l$file_log)
)

testExistanceAndCreateDirectoriesRecursively(allDirs)


assertCharacter(par.l$colorCategories, len = 4)
assertSubset(names(par.l$colorCategories), c("activator", "undetermined", "repressor", "not-expressed"))

assertCharacter(par.l$colorConditions, len = 2)

######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel,  removeOldLog = TRUE)
printParametersLog(par.l)

conditionComparison = readRDS(par.l$file_input_condCompDeSeq)
assertVector(conditionComparison, len = 2)

# Assemble the final table and collect permutation information from all TFs
output.global.TFs.orig = NULL


nTF = length(par.l$files_input_permResults)
for (fileCur in par.l$files_input_permResults) {
  
  resultsCur.df =  read_tsv(fileCur, col_names = TRUE, col_types = cols())
  assertIntegerish(nrow(resultsCur.df), lower = 1, upper = par.l$nPermutations + 1)
  
  if (is.null(output.global.TFs.orig)) {
    output.global.TFs.orig = resultsCur.df
  } else {
    output.global.TFs.orig = rbind(output.global.TFs.orig, resultsCur.df)
  }
  
}

# Remove TFs with NA values
output.global.TFs = filter(output.global.TFs.orig, permutation == 0, is.finite(weighted_meanDifference))


output.global.TFs$weighted_meanDifference_enrichment = NA

variableXAxis = "weighted_meanDifference"


# Remove rows with NA
TF_NA = which(is.na(output.global.TFs$weighted_meanDifference) | is.na(output.global.TFs$variance) | is.na(output.global.TFs$weighted_Tstat))

if (length(TF_NA) > 0) {
  output.global.TFs = output.global.TFs[-TF_NA,]
  
  TFs_NA = output.global.TFs$TF[TF_NA]
  message = paste0("The following TF have been removed from the data due to NA values in either weighted_meanDifference, variance or weighted_Tstat (insufficient data in previous steps): ", paste0(TFs_NA, collapse = ", "))
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
}


######################################
# FILTER BY PERMUTATIONS AND COMPARE #
######################################
# Compare the distributions from the real and random permutations
if (par.l$nPermutations > 0) {
  
  output.global.TFs.permutations = filter(output.global.TFs.orig, permutation > 0)
  
  
  thresholdsPermutations = quantile(output.global.TFs.permutations$weighted_meanDifference, probs = par.l$probsThreshold, na.rm = TRUE)
  flog.info(paste0("Thresholds for the permutations: ", paste0(thresholdsPermutations, collapse = " - ")))
  
  # Calculate a signal to noise
  TFs_positive = which(output.global.TFs$weighted_meanDifference >= 0)
  TFs_negative = which(output.global.TFs$weighted_meanDifference <  0)
  output.global.TFs$weighted_meanDifference_enrichment[TFs_negative] = output.global.TFs$weighted_meanDifference[TFs_negative] / thresholdsPermutations[1]
  output.global.TFs$weighted_meanDifference_enrichment[TFs_positive] = output.global.TFs$weighted_meanDifference[TFs_positive] / thresholdsPermutations[2]
  
  # Summarize the random permutations
  permutationDensities = ggplot(output.global.TFs.orig, aes(x = weighted_meanDifference, fill = factor(permutation))) + geom_density(alpha = 0.25) + theme_bw() + geom_vline(xintercept =  thresholdsPermutations, linetype = "dotted") 
  
  permutationDensities2 = ggplot(output.global.TFs.permutations, aes(x = weighted_meanDifference, fill = factor(permutation))) + geom_density(alpha = 0.25) + theme_bw() + geom_vline(xintercept =  thresholdsPermutations, linetype = "dotted") 
  
  permutationDensities3 = ggplot(output.global.TFs, aes(x = weighted_meanDifference)) + geom_density(alpha = 0.25) + theme_bw() + geom_vline(xintercept =  thresholdsPermutations, linetype = "dotted") 
  
  variableXAxis = "weighted_meanDifference_enrichment"
  
  
}


output.global.TFs$Cohend_factor = ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[1], par.l$classes_CohensD[1], 
                                         ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[2] , par.l$classes_CohensD[2], 
                                                ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[3], par.l$classes_CohensD[3], par.l$classes_CohensD[4])))

output.global.TFs$Cohend_factor = factor(output.global.TFs$Cohend_factor, levels = par.l$classes_CohensD, labels = seq_len(length(par.l$classes_CohensD)))


##################################################
# DIAGNOSTIC PLOTS AND SIGNIFICANCE CALCULATIONS #
##################################################

pdf(par.l$file_plotDiagnostic)

if (par.l$nPermutations > 0) {
  plot(permutationDensities)
  plot(permutationDensities2)
  plot(permutationDensities3)
}

estimates = tryCatch( {
  
  locfdrRes = locfdr(output.global.TFs$weighted_Tstat, plot = 4)
  # Currently taken as default
  MLE.delta  = locfdrRes$fp0["mlest", "delta"]
  
  # Not the current default
  CME.delta  = locfdrRes$fp0["cmest", "delta"]
  
  c(MLE.delta, CME.delta)
  
}, error = function(e) {
  message = "Could not run locfdr, use the median instead..."
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  
  # For the fallback, simply return the median for both estimates for now
  c(median(output.global.TFs$weighted_Tstat, na.rm = TRUE), median(output.global.TFs$weighted_Tstat, na.rm = TRUE))
}
)

MLE.delta = estimates[1]
CME.delta = estimates[2]

# Compute two different measures for the mean
# Which one to take? Depends, see page 101 in Bradley Efron-Large-Scale Inference_ Empirical Bayes Methods for Estimation, Testing, and Prediction
# The default option in locfdr is the MLE method, not central matching. Slight irregularities in the central histogram,
# as seen in Figure 6.1a, can derail central matching. The MLE method is more stable, but pays the price of possibly increased bias.


# 2. Should the variance associated with this weighted.mean T stat value be based on the variance of the T stat scores, the variance of the mean, or both?

# 3. Why not estimate the SD of the distribution directly, without the need to provide manual variance estimations? Because then each weigted mean value
# is treated as if it came from the same population? 

# 4. library(Hmisc) -> wtd.var(x, weights)

# 5 . Variance estimation  is based on before the centralization, is this correct? Yes, independent

# 6. Is centralization actually correct at all since the scores come from a different population? Yes, if assumed they come from a normal distribution

# NEW: Estimate variance of T score estimate with Welch corrected df (test$parameter)

# READ http://genomicsclass.github.io/book/pages/t-tests_in_practice.html

output.global.TFs$weighted_Tstat_centralized = output.global.TFs$weighted_Tstat - MLE.delta

output.global.TFs$variance = as.numeric(output.global.TFs$variance)
output.global.TFs$pvalue    = 2*pnorm(-abs(output.global.TFs$weighted_Tstat_centralized), sd = sqrt(output.global.TFs$variance))
output.global.TFs$pvalueAdj = p.adjust(output.global.TFs$pvalue, method = "BH")

# Handle extreme cases with p-values that are practically 0 and would cause subsequent issues
index0 = which(output.global.TFs$pvalueAdj < .Machine$double.xmin)
if (length(index0) > 0) {
  output.global.TFs$pvalueAdj[index0] = .Machine$double.xmin
}

output.global.TFs$yValue = transform_yValues(output.global.TFs$pvalueAdj)


colnamesToPlot = colnames(output.global.TFs)[-which(colnames(output.global.TFs) %in% c("TF", "sign"))]

for (FDRThresholdCur in c(par.l$FDR_threshold, 1)) {
  
  filtered.df = filter(output.global.TFs, pvalueAdj <= FDRThresholdCur)
  
  title = paste0("FDR: ", FDRThresholdCur, " (retaining ", nrow(filtered.df), " TF)")
  
  for (measureCur in colnamesToPlot) {
    
    if (measureCur == "weighted_meanDifference_enrichment" & par.l$nPermutations == 0) {
      next
    }
    
    if (measureCur %in%  c("Cohend_factor")) {
      plot(ggplot(filtered.df, aes_string(measureCur))  + stat_count() + theme_bw() + ggtitle(title))
      
    } else {
      plot(ggplot(filtered.df, aes_string(measureCur))  + geom_histogram(bins = 50) + theme_bw() + ggtitle(title))
    }
    
    
    
  }
}

dev.off()


##########################
# INTEGRATE RNA-Seq DATA #
##########################
if (par.l$plotRNASeqClassification) {
  
  classesList.l = list(c("activator","undetermined","repressor","not-expressed"),
                       c("activator","undetermined","repressor"),
                       c("activator","repressor")
  )
  
  extensionSize = as.integer(snakemake@config$par_general$regionExtension)
  assertIntegerish(extensionSize)
  
  rootOutdir = snakemake@config$par_general$outdir
  assertCharacter(rootOutdir)
  
  comparisonType = snakemake@config$par_general$comparisonType
  assertCharacter(comparisonType)
  
  if (nchar(comparisonType) > 0) {
    comparisonType = paste0(comparisonType, ".")
  }
  
  file_sampleSummary = snakemake@config$samples$summaryFile
  assertFileExists(file_sampleSummary)
  
  if (grepl(pattern = "/scratch/carnold/CLL/TF_act", rootOutdir)) {
    RNASeqSummary = "/scratch/carnold/CLL/TF_act_10samplesOnlyRNASeq/output/FINAL_OUTPUT/extension100/tf_activator_repressor_classification.RData"
    flog.info(paste0("Load object ", RNASeqSummary))
    load(RNASeqSummary)
  } else {
    
    sampleSummary.df = read_tsv(file_sampleSummary, col_types = cols())
    
    
    # Loading TF data
    HOCOMOCO_mapping.df = read.table(file = par.l$file_input_HOCOMOCO_mapping, header = TRUE)
    assertSubset(c("ENSEMBL", "HOCOID"), colnames(HOCOMOCO_mapping.df))
    
    # This is the gene read count matrix only for TFs coming from RNAseq
    TF.counts.df.all = read.table(par.l$file_input_geneCountsPerSample, header = TRUE)
    
    # Check sample names and set column names
    # Match the column names and do the intersections
    sharedColumns = intersect(colnames(TF.counts.df.all)[-1], sampleSummary.df$SampleID)
    
    if (length(sharedColumns) == 0) {
      message = paste0("No shared samples with RNA-Seq samples between sample table ", file_sampleSummary, " and RNA-Seq table ", par.l$file_input_geneCountsPerSample, ".")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }

    colnames(TF.counts.df.all)[1] = "ENSEMBL"
    
    # Clean ENSEMBL IDs
    TF.counts.df.all$ENSEMBL = gsub("\\..+", "", TF.counts.df.all$ENSEMBL, perl = TRUE)
    
    # TODO: Matching of names
    
    # Filter them by the IDs that correspond to the TFs
    TF.counts.df = filter(TF.counts.df.all, ENSEMBL %in% HOCOMOCO_mapping.df$ENSEMBL)
    
    if (nrow(TF.counts.df) == 0) {
      message = "No rows remaining after filtering against ENSEMBL IDs in HOCOMOCO. Check your ENSEMBL IDs for overlap with the HOCOMOCO translation table."
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    
    # Read counts per TFBS, coming from DESeqPeaks.R
    DESeq.obj = readRDS(par.l$file_input_DESeq)
    peak.counts = counts(DESeq.obj, normalized = TRUE)
    
    # Match the column names and do the intersections
    sharedColumns = intersect(colnames(peak.counts), colnames(TF.counts.df))
    
    if (length(sharedColumns) == 0) {
      message = "No shared samples with RNA-Seq samples."
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    flog.info(paste0(length(sharedColumns), " samples are shared between the input data and RNA-Seq data"))
    
    #peak.counts.orig = peak.counts
    peak.counts  = peak.counts [, which(colnames(peak.counts) %in% sharedColumns)]
    TF.counts.df = TF.counts.df[, which(colnames(TF.counts.df) %in% c(sharedColumns, "ENSEMBL"))]

    
    HOCOMOCO_mapping.df.exp <- filter(HOCOMOCO_mapping.df, ENSEMBL %in%  TF.counts.df$ENSEMBL, HOCOID %in% unique(output.global.TFs$TF))
    
    if (nrow(HOCOMOCO_mapping.df.exp) == 0) {
      message = paste0("Number of rows of HOCOMOCO_mapping.df.exp is 0. Something is wrong with the mapping table or the filtering")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    TF.peakMatrix.l = list()
    for (TFCur in HOCOMOCO_mapping.df.exp$HOCOID) {
      
      HOCOMOCO_mapping.subset.df = subset(HOCOMOCO_mapping.df.exp, HOCOID == TFCur)
      gene.sel = unique(HOCOMOCO_mapping.subset.df$ENSEMBL)
      if (length(gene.sel) > 1) {
        message = paste0("Mapping for ", TFCur, " not unique, take only the first mapping (", gene.sel[1], ") and discard the others (", paste0(gene.sel[-1], collapse = ","), ")")
        checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        
      }
      
      
      TF.output.df = read.table(file = paste0(rootOutdir, "/TF-SPECIFIC/",TFCur,"/extension", extensionSize, "/", comparisonType, TFCur,  ".output.tsv"), header = TRUE)
      TF.output.df = filter(TF.output.df, permutation == 0)
      assertSubset(c("chr", "PSS", "PES"), colnames(TF.output.df))
      TF.output.df$ID2 = paste0(TF.output.df$chr, ":", TF.output.df$PSS, "-", TF.output.df$PES)
      TF.peakMatrix.l[[TFCur]] = rownames(peak.counts) %in% TF.output.df$ID2
    }
    
    
    # cor.m = is the peaks names separated with # 
    
    # This is the peak (rows) and TF binding sites (columns)
    TF.peakMatrix.df = as.data.frame(TF.peakMatrix.l)
    
    # Sanity check
    if (all(rowSums(TF.peakMatrix.df) == 0)) {
      message = paste0("All counts are 0, something is wrong.")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    
    # Expressed TFs only 
    HOCOMOCO_mapping.subset.df = subset(HOCOMOCO_mapping.df.exp, HOCOID %in% colnames(TF.peakMatrix.df))
    
    # Remove first column, retain only counts
    expressed.TF.counts.df = t(TF.counts.df[,-c(1)])
    colnames(expressed.TF.counts.df) = TF.counts.df$ENSEMBL
    expressed.TF.counts.df = t(expressed.TF.counts.df)
    
    assertSubset(colnames(expressed.TF.counts.df), colnames(peak.counts))
    
    # Some rownames may be identical because of the mapping from HOCOMOCO
    
    peak.counts = peak.counts[,order(colnames(peak.counts))]
    expressed.TF.counts.df = expressed.TF.counts.df[,order(colnames(expressed.TF.counts.df))]
    
    if (!all(colnames(expressed.TF.counts.df) == colnames(peak.counts))) {
      message = "Colnames of expressed.TF.counts.df and peak.counts must be identical"
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    # Filter by rowMeans to eliminate rows with an sd of 0
    rowMeans1 = rowMeans(expressed.TF.counts.df)
    rowsToDelete = which(rowMeans1 < 1)
    if (length(rowsToDelete) > 0) {
      expressed.TF.counts.df = expressed.TF.counts.df[-rowsToDelete,]
      flog.info(paste0("Removed ", length(rowsToDelete), " TFs out of ", nrow(expressed.TF.counts.df), " because they had a row mean of < 1."))
    }
    rowMeans2 = rowMeans(peak.counts)
    rowsToDelete = which(rowMeans2 == 0)
    if (length(rowsToDelete) > 0) {
      peak.counts = peak.counts[-rowsToDelete,]
      flog.info(paste0("Removed ", length(rowsToDelete), "peaks out of ", nrow(peak.counts), " because they had a row mean of < 1."))
    }
    
    cor.m = t(cor(t(expressed.TF.counts.df), t(peak.counts), method = par.l$corMethod))
    
    # Mapping TFBS to TF 
    sort.cor.m = cor.m[,names(sort(colMeans(cor.m)))] 
    
    HOCOMOCO_mapping.df.exp = HOCOMOCO_mapping.df.exp[which(HOCOMOCO_mapping.df.exp$ENSEMBL %in% colnames(sort.cor.m)),]
    assertSubset(HOCOMOCO_mapping.df.exp$ENSEMBL, colnames(sort.cor.m))
    sort.cor.m = sort.cor.m[,HOCOMOCO_mapping.df.exp$ENSEMBL] 
    colnames(sort.cor.m) = as.character(HOCOMOCO_mapping.df.exp$HOCOID)
    sel.TF.peakMatrix.df = TF.peakMatrix.df[,colnames(sort.cor.m)]
    
    t.cor.sel.matrix = sort.cor.m
    t.cor.sel.matrix[(sel.TF.peakMatrix.df == 0)] = NA
    t.cor.sel.matrix = sel.TF.peakMatrix.df * t.cor.sel.matrix
    
    my.median = function(x) median(x, na.rm = TRUE)
    my.mean   = function(x) mean(x, na.rm = TRUE)
    median.cor.tfs = sort(apply(t.cor.sel.matrix, MARGIN = 2, FUN = my.median))
    t.cor.sel.matrix.non = sort.cor.m
    t.cor.sel.matrix.non[(sel.TF.peakMatrix.df == 1)] = NA
    t.cor.sel.matrix.non = sel.TF.peakMatrix.df + t.cor.sel.matrix.non
    act.rep.thres = quantile(sort(apply(t.cor.sel.matrix.non, MARGIN = 2, FUN = my.median)), probs = c(.05, .95))
    
    
    
    
  }  # end if special case for analysis
  
  AR.data = as.data.frame(median.cor.tfs)
  AR.data$TF = rownames(AR.data)
  
  output.global.TFs = merge(output.global.TFs, AR.data, by = "TF",all.x = TRUE)
  
  output.global.TFs$classification = ifelse(is.na(output.global.TFs$median.cor.tfs), "not-expressed",
                                            ifelse(output.global.TFs$median.cor.tfs <= act.rep.thres[1], "repressor",
                                                   ifelse(output.global.TFs$median.cor.tfs > act.rep.thres[2], "activator", "undetermined")))
  
  output.global.TFs$classification = factor(output.global.TFs$classification, levels = names(par.l$colorCategories))
  
} else {
  classesList.l = list(c())
}



output.global.TFs.origReal = output.global.TFs


###########################
# PLOT FOR DIFFERENT FDRs #
###########################

# Set the page dimensions to the maximum across all plotted variants
output.global.TFs.filteredSummary = filter(output.global.TFs, pvalueAdj < max(par.l$FDR_threshold))

if (par.l$nPermutations > 0) {
  output.global.TFs.filteredSummary = filter(output.global.TFs.filteredSummary, weighted_meanDifference < thresholdsPermutations[1] | weighted_meanDifference > thresholdsPermutations[2])
} else {
  output.global.TFs.filteredSummary = filter(output.global.TFs.filteredSummary, weighted_CD > par.l$cohensDThreshold)
}

nTF_label = nrow(output.global.TFs.filteredSummary)

height = width = max(nTF_label / 5, par.l$circularPlot_minDimensions)

# Increase a bit towards smaller heights by a factor, empirical observation
if (height < 20) {
  height = width = height * 1.2
}


permutationCategories = c(FALSE)
if (par.l$nPermutations > 0) permutationCategories = c(TRUE, FALSE)

pdf(file = par.l$file_plotCircular, height = height, width = width, useDingbats = FALSE)


for (includePermutationCur in permutationCategories) {
  

  
  for (FDRThresholdCur in par.l$FDR_threshold) {
    
    output.global.TFs.orig = output.global.TFs.origReal
    output.global.TFs.orig$sign = output.global.TFs.orig$pvalueAdj < FDRThresholdCur
    
    
    indexLoop = 0
    for (showClasses in classesList.l) {
      
      #showClasses = classesList.l[1]
      output.global.TFs = output.global.TFs.orig
      # First create the data and subset the TFs that should be plotted
      
      
      if (par.l$plotRNASeqClassification) {
        # Filter by classes to show
        output.global.TFs = filter(output.global.TFs, classification %in% showClasses)
        output.global.TFs$classification = factor(output.global.TFs$classification, levels = showClasses)
        
        
        colorCategoriesCur = par.l$colorCategories[which(names(par.l$colorCategories) %in% showClasses)]
      }
      
      
      # Limits for x-axis. Multiple by a factor > 1 to account for the white legend part in which no point should be located
      # For this, extend the x axis limit so that radial positions are smaller and do not cross the border of the legend
      limit_max = max(abs(output.global.TFs[,variableXAxis])) * (1 + par.l$extension_x_limits)
      
      
      limit_min = 0
      
      # Determine the value for 1 degree radial positions. 180 is the max
      radial_coeff = 180/limit_max
      
      # Calculate the radial positions
      output.global.TFs$radial = 0
      index_neg = which(output.global.TFs$weighted_meanDifference < 0)
      index_pos = which(output.global.TFs$weighted_meanDifference > 0)
      
      if (includePermutationCur) {
        index_neg_thresholdPermuations = which(thresholdsPermutations < 0)
        index_pos_thresholdPermuations = which(thresholdsPermutations > 0)
        threshold_radial_neg = NULL
        threshold_radial_pos = NULL
        if (length(index_neg_thresholdPermuations) > 0) {
          threshold_radial_neg = 180 + (radial_coeff*abs(thresholdsPermutations[index_neg_thresholdPermuations]))
          threshold_radial_neg = 180 + (radial_coeff*abs(1)) # at boundary, enrichment is by definition 1
        }
        if (length(index_pos_thresholdPermuations) > 0) {
          threshold_radial_pos = 180 - (radial_coeff*abs(thresholdsPermutations[index_pos_thresholdPermuations]))
          threshold_radial_pos = 180 - (radial_coeff*abs(1))
        }
        
        threshold_radial_all = c(threshold_radial_neg, threshold_radial_pos)
        
      }
      
      
      # extracting the labels from object in order to add them manually later
      
      
      output.global.TFs$radial[index_neg] = 180 + (radial_coeff*abs(unlist(output.global.TFs[, variableXAxis])[index_neg]))
      output.global.TFs$radial[index_pos] = 180 - (radial_coeff*abs(unlist(output.global.TFs[, variableXAxis])[index_pos]))
      
      # Which measure to use for the y value?
      signThresholdPlot = transform_yValues(FDRThresholdCur)
      
      # Dont display TFs that are deemed non-significant or that have a small Cohens D value
      
      ggrepel_df = filter(output.global.TFs, sign == TRUE)
      
      if (includePermutationCur) {
        ggrepel_df = filter(ggrepel_df, weighted_meanDifference < thresholdsPermutations[1] | weighted_meanDifference > thresholdsPermutations[2])
      } else {
        ggrepel_df = filter(ggrepel_df, weighted_CD > par.l$cohensDThreshold)
      }
      
      ## ggplot p objects
      ylimit = sign(max(output.global.TFs$yValue)) * ceiling(abs(max(output.global.TFs$yValue))) + par.l$extension_y_limits
      
      labels        = list()
      #labels$range  = layer_scales(p1)$y$range$range
      
      startNo = sign(min(output.global.TFs$yValue)) * ceiling(abs(min(output.global.TFs$yValue)))
      endNo   = ylimit
      stepsize = 1
      labels$breaks = round(seq(startNo, endNo, stepsize),0)
      
      
      # 
      # fillVar  = NULL
      # if (par.l$plotRNASeqClassification) {
      #   fillVar = ggrepel_df$classification
      # }
      
      width_whiteArea = 15
      
      #p1 = ggplot(environment = localenv) + 
      p1 = ggplot() +   
        geom_rect(data = NULL,
                  aes(xmin = 0,
                      xmax = 180,
                      ymin = -Inf, 
                      ymax = ylimit),
                  alpha = .45, 
                  fill = par.l$colorConditions[1]) + 
        geom_rect(data = NULL,
                  aes(xmin = 180,
                      xmax = 360,
                      ymin = -Inf, 
                      ymax = ylimit),
                  alpha = .45, 
                  fill = par.l$colorConditions[2]) + 
        geom_rect(data = NULL, aes(xmin = 0, xmax = 360 , ymin = -Inf, ymax = signThresholdPlot), alpha = 0.5, fill = "white") +
        coord_polar() + 
        # limits of the angular plot 
        scale_x_continuous(limits = c(0,360)) +
        # sizes of the points for significance
        scale_size_manual(values = c(0.5,1), guide = FALSE) 
      
      
      if (includePermutationCur) {
        
        p1 = p1 + geom_segment(aes(x = threshold_radial_all, xend = threshold_radial_all, y=signThresholdPlot, yend = max(labels$breaks)), size = 1.2, linetype = "solid", alpha = 0.5, color = "red") 
        
        p1 = p1 + geom_rect(data = NULL, aes(xmin = threshold_radial_all[1], xmax = threshold_radial_all[2] , ymin = signThresholdPlot, ymax = max(labels$breaks)), alpha = 0.5, fill = "white") 
      }
      
      
      
      p1 = p1 + geom_hline(yintercept = ylimit,
                           size = 0.8,
                           linetype = "solid", 
                           color = par.l$plot_grayColor,
                           alpha = .9) +
        # Make the white area separate into two to avoid strange artefacts
        geom_rect(data=NULL,aes(xmin = 360 - width_whiteArea, xmax = 360 , ymin = -Inf, ymax = ylimit + 0.1), alpha = 1, fill = "white") +
        geom_rect(data=NULL,aes(xmin = 0 , xmax = width_whiteArea , ymin = -Inf, ymax = ylimit + 0.1), alpha = 1, fill = "white")
      
      
      
      if (par.l$plotRNASeqClassification) {
        p1 = p1 +  geom_point(data = output.global.TFs, aes(x = radial, 
                                                            y = yValue,
                                                            #label = output.global.TFs$TF,
                                                            size = sign,
                                                            fill = classification,
                                                            color = classification), 
                              alpha = 1, shape = 21)
      } else {
        p1 = p1 +  geom_point(data = output.global.TFs, aes(x = radial, 
                                                            y = yValue,
                                                            #label = TF,
                                                            size = sign), 
                              alpha = 1, shape = 21)
      }
      
      
      # Draw p value significance threshold line: 1 or 2 lines, depending on the permutations
      if (includePermutationCur) {
        
        p1 = p1 + geom_segment(aes(x = 0, xend = threshold_radial_all[2]  , y = signThresholdPlot, yend = signThresholdPlot),  
                               size = 1.2,
                               linetype = "solid",  # before: "longdash"
                               color = "red",
                               alpha = .5) + 
          geom_segment(aes(x = threshold_radial_all[1], xend = 360, y = signThresholdPlot, yend = signThresholdPlot),  
                       size = 1.2,
                       linetype = "solid",  # before: "longdash"
                       color = "red",
                       alpha = .5)
      } else {
        
        p1 = p1 + geom_hline( yintercept = signThresholdPlot,
                              size = 1.2,
                              linetype = "longdash",
                              color = "red",
                              alpha = .5)
        
      }
      
      
      fontSize = par.l$rootFontSize + nrow(ggrepel_df) * 0.01
      
      # y axis labels
      for (i in 2:(length(labels$breaks) - 1)) {
        # p3 = p3 + annotate(geom = "text", x = 15, y = -labels$breaks[i], label = paste0(labels$breaks[i]), vjust = 0.5, hjust = 1.5, angle = -15, size = fontSize)
        p1 = p1 + annotate(geom = "text", x = 0, y = labels$breaks[i], label = paste0(as.numeric(labels$breaks[i])), vjust = 0.5, hjust = 0.5, angle = 0, size = fontSize)
        
      }
      
      p2 = p1
      
      
      # Increase the size of the poitns in the legend, see https://stackoverflow.com/questions/20415963/how-to-increase-the-size-of-points-in-legend-of-ggplot2
      #p2 = p2 + guides(fill = guide_legend(override.aes = list(size = par.l$sizeLegend), nrow = 2))
      p2 = p2 + guides(fill = guide_legend(override.aes = list(size = par.l$sizeLegend), nrow = 1))
      
      # Add the annotation outside of the plot (y axis) and draw helper lines
      
      ## scaling to the angles automatically
      
      limitMaxAxis = limit_max
      
      
      coeff_angle = (limitMaxAxis * 2) / 360
      
      anglesPos = c(seq(180, 15, -30))
      anglesNeg = c(seq(180 + 30, 345, 30))
      
      df.axis = data.frame(angles = c(anglesPos, anglesNeg), 
                           annotation = c(coeff_angle  * abs(anglesPos - 180),
                                          -coeff_angle * abs(anglesNeg - 180)))
      
      
      
      p2 = p2 + geom_segment(aes(x=df.axis$angles, xend = df.axis$angles , y=min(labels$breaks), yend = max(labels$breaks)), 
                             size = par.l$sizeHelperLines, linetype = "dotted", color = par.l$plot_grayColor, alpha = .9) 
      
      if (includePermutationCur) {
        
        df.axis$annotation = round(abs(df.axis$annotation), 1)
        
      } else {
        df.axis$annotation = round(df.axis$annotation, 3)
      }
      
      
      for (i in anglesPos) {
        
        p2 = p2 + annotate(geom = "text", x = i, y = ylimit + 0.55, label = paste0(df.axis[which(df.axis$angles == i),]$annotation), vjust = 0, hjust = 0, size = fontSize)
      }
      for (i in anglesNeg) {
        
        p2 = p2 + annotate(geom = "text", x = i, y = ylimit + 0.55, label = paste0(df.axis[which(df.axis$angles == i),]$annotation), vjust = 0, hjust = 1, size = fontSize)
      }
      
      
      
      p3 = p2 
      
      labelY = paste0("log10(-log10(adj. p-value) (FDR: ", round(FDRThresholdCur*100,0), "%)")
      labelY = paste0("Significance (FDR: ", round(FDRThresholdCur*100,2), "%)")
      
      #p3 = p3 + annotate(geom = "text", x = 0, y = 0, label = "log10 T stat.", vjust = 0.5, hjust = 0.1, angle = 90, size = fontSize, fontface = 'italic') + 
      p3 = p3 + annotate(geom = "text", x = 0, y =  ylimit + 1, label = labelY, vjust = 0, hjust = 0.5, angle = 0, size = fontSize, fontface = 'italic')
      
      angleLegAct = 180
      length_arrow = 80
      yOffset = 2.25
      yOffsetLabels = yOffset - 0.25
      
      
      
      # Draw the legend in the lower part (arrows)
      p3 = p3 + 
        annotate(geom = "text", x = angleLegAct + 0, y = ylimit + yOffset, label = " TF activity ", vjust = 0, angle = 0, size = fontSize) + 
        # First right side for positive values, corresponding to positive weighted mean difference values (first element of DESeq comparison)
        annotate(geom = "text", x = angleLegAct - 30, y = ylimit + yOffsetLabels, label = conditionComparison[1], vjust = 0, angle = 30, size = fontSize) + 
        # Now left side
        annotate(geom = "text", x = (angleLegAct + 30), y = ylimit + yOffsetLabels, label = conditionComparison[2], vjust = 0, angle = -30, size = fontSize) + 
        geom_segment(aes(x = angleLegAct + 10, xend = angleLegAct + 10 + length_arrow , y = ylimit + yOffset, yend = ylimit + yOffset), size=0.3, arrow = arrow(length = unit(0.6,"cm")))  +
        geom_segment(aes(x=angleLegAct - 10, xend = angleLegAct - 10 - length_arrow , y = ylimit + yOffset, yend = ylimit + yOffset), size=0.3, arrow = arrow(length = unit(0.6,"cm"))) 
      
      # Add minor y circular lines 
      for (x in 1:length(labels$breaks)) {
        p3 = p3 + geom_hline(yintercept = labels$breaks[x], size = par.l$sizeHelperLines, linetype = "dotted",  color = par.l$plot_grayColor, alpha = .9) 
      }
      
      
      # ggrepel function
      
      if (par.l$plotRNASeqClassification) {
        p3 = p3 + scale_color_manual(values = colorCategoriesCur, guide = FALSE) + 
          geom_label_repel(data = ggrepel_df, aes(x = radial,
                                                  y = yValue,
                                                  label = TF,
                                                  fill = classification),
                           # size of the label
                           size = par.l$size_TFAnnotation,
                           fontface = 'bold', color = 'white',
                           segment.size = 0.15,
                           # how thick is connectin line
                           label.padding = unit(0.2, "lines"),
                           # how far from center points
                           nudge_y = 0.15,
                           nudge_x = 0, 
                           segment.alpha = .8,
                           segment.color = par.l$plot_grayColor, show.legend = FALSE) + 
          scale_fill_manual(values = colorCategoriesCur, name = "TF class")
      } else {
        
        p3 = p3 + geom_label_repel(data = ggrepel_df, aes(x = radial,
                                                          y = yValue,
                                                          label = TF),
                                   # size of the label
                                   fill = par.l$plot_grayColor,
                                   size = par.l$size_TFAnnotation,
                                   fontface = 'bold', color = 'white',
                                   segment.size = 0.15,
                                   # how thick is connectin line
                                   label.padding = unit(0.1, "lines"),
                                   # how far from center points
                                   nudge_y = 0.15,
                                   nudge_x = 0, 
                                   segment.alpha = .8,
                                   segment.color = par.l$plot_grayColor, show.legend = FALSE)
      }
      
      
      
      p3 = p3 +  
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_text(),
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid = element_blank(),
              legend.position = "top",
              legend.text = element_text(size = par.l$sizeLegend),
              legend.title = element_text(size = par.l$sizeLegend, face = "bold"),
              legend.margin = margin(t = 0, unit = "cm"),
              legend.key = element_blank(),
              legend.justification = "center",
              #plot.margin = grid::unit(c(0, 0, 0, 0), "mm")) + 
              #plot.margin = unit(c(0,-7,-4,-7), units = "cm")) + 
              plot.margin = unit(c(0,0,0,0), units = "cm")) +
        xlab("") + ylab("") 
      
      print(p3)
    } # end for all showClasses
    
    
    
  } # end for different FDRs
  
} # end for with and without permutations

dev.off()

write_tsv(output.global.TFs.origReal, path = par.l$file_output_summary, col_names = TRUE)

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)
