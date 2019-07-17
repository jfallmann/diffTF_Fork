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
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/summaryFinal.R.rds")
# Note that one currently cannot overwrite valeus within the snakemake object; instead, assign them to a variable as done below and change if necessary
createDebugFile(snakemake)

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "ggrepel", "checkmate", "tools", "grDevices", "locfdr"), verbose = FALSE)


###################
#### PARAMETERS ###
###################

# Currently hard-coded parameters
par.l = list()

# 1. Misc
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"

# 2. Statistical thresholds and values
par.l$significanceThresholds  = c(0.001, 0.01, 0.05,0.1,0.2) # p-value thresholds
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.1, 0.5, 0.8)

# 3. RNA-Seq specific
par.l$corMethod = "pearson" # Expression-peak count correlation method. As we quantile normalize now, should be pearson
par.l$regressionMethod = "glm" # for correlating RNA-Seq classification with TF activity
par.l$filter_minCountsPerCondition = 5 # For filtering RNA-seq genes, see Documentation

# 4. Volcano plot settings
par.l$maxTFsToLabel = 150 # Maximum Tfs to label in the Volcano plot
par.l$volcanoPlot_minDimensions  = 12
par.l$minPointSize = 0.3
par.l$plot_grayColor = "grey50"
par.l$colorCategories = c("activator" = "#4daf4a", "undetermined" = "black", "repressor" = "#e41a1c", "not-expressed" = "Snow3") # diverging, modified
par.l$colorConditions = c("#ef8a62", "#67a9cf") # Colors of the two conditions for the background

#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(c("", "allPermutationResults", "condComp", "normCounts"), names(snakemake@input))


par.l$files_input_permResults  = snakemake@input$allPermutationResults
for (fileCur in par.l$files_input_permResults) {
  assertFileExists(fileCur, access = "r")
}

par.l$file_input_condCompDeSeq = snakemake@input$condComp
assertFileExists(par.l$file_input_condCompDeSeq, access = "r")

par.l$file_input_countsNorm = snakemake@input$normCounts
assertFileExists(par.l$file_input_countsNorm, access = "r")

par.l$file_input_metadata = snakemake@input$sampleDataR
assertFileExists(par.l$file_input_metadata, access = "r")

## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(c("", "summary", "volcanoPlot", "diagnosticPlots", "plotsRDS"), names(snakemake@output))

par.l$file_output_summary  = snakemake@output$summary
par.l$file_plotVolcano     = snakemake@output$volcanoPlot
par.l$files_plotDiagnostic = snakemake@output$diagnosticPlots
par.l$file_output_plots    = snakemake@output$plotsRDS

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

allDirs = c(dirname(par.l$file_output_summary), dirname(par.l$files_plotDiagnostic),dirname(par.l$file_log))
testExistanceAndCreateDirectoriesRecursively(allDirs)

assertCharacter(par.l$colorCategories, len = 4)
assertSubset(names(par.l$colorCategories), c("activator", "undetermined", "repressor", "not-expressed"))
assertCharacter(par.l$colorConditions, len = 2)

######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel,  removeOldLog = TRUE)
printParametersLog(par.l)


################
# COLLECT DATA #
################

conditionComparison = readRDS(par.l$file_input_condCompDeSeq)
assertVector(conditionComparison, len = 2)

# Assemble the final table and collect permutation information from all TFs
output.global.TFs.orig = NULL

nTF = length(par.l$files_input_permResults)
for (fileCur in par.l$files_input_permResults) {
  
  resultsCur.df =  read_tsv(fileCur, col_names = TRUE)
  # resultsCur.df =  read_tsv(fileCur, col_names = TRUE, col_types = list(
  #   col_integer(), # "permutation"
  #   col_character(), # "TF",
  #   col_double(), # "weighted_meanDifference
  #   col_double(), # weighted_CD
  #   col_double(), # TFBS
  #   col_double(), # weighted_Tstat
  #   col_double() # variance
  # ))
  assertIntegerish(nrow(resultsCur.df), lower = 1, upper = par.l$nPermutations + 1)
  
  if (is.null(output.global.TFs.orig)) {
    output.global.TFs.orig = resultsCur.df
  } else {
    output.global.TFs.orig = rbind(output.global.TFs.orig, resultsCur.df)
  }
  
}

# Convert columns to numeric if they are not already
output.global.TFs.orig = mutate(output.global.TFs.orig,
                                weighted_meanDifference = as.numeric(weighted_meanDifference),
                                variance                = as.numeric(variance),
                                weighted_CD             = as.numeric(weighted_CD),
                                weighted_Tstat          = as.numeric(weighted_Tstat))

# Remove rows with NA
TF_NA = which(is.na(output.global.TFs.orig$weighted_meanDifference))

if (length(TF_NA) > 0) {
    output.global.TFs.orig = output.global.TFs.orig[-TF_NA,]
  
  TFs_NA = output.global.TFs.orig$TF[TF_NA]
  message = paste0("The following TF have been removed from the data due to NA values in weighted_meanDifference (insufficient data in previous steps): ", paste0(unique(TFs_NA), collapse = ", "))
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
}


########################################################
# FILTER BY PERMUTATIONS AND COMPARE, DIAGNOSTIC PLOTS #
########################################################
# Compare the distributions from the real and random permutations
diagPlots.l = list()

pdf(par.l$files_plotDiagnostic[1])
output.global.TFs.orig$pvalue = NA
if (par.l$nPermutations > 0) {

    dataReal.df = dplyr::filter(output.global.TFs.orig, permutation == 0)
    dataPerm.df = dplyr::filter(output.global.TFs.orig, permutation > 0)
  

    # ( (#perm >TH)/#perm ) / ( (#perm >TH)/#perm + (#real > TH)/n_real )
    
    xrange = range(output.global.TFs.orig$weighted_meanDifference, na.rm = TRUE)
    
    plot(density(dataPerm.df$weighted_meanDifference, na.rm = TRUE), col = "black", main = "Weighted mean difference values (black = permuted)", xlim = xrange)
    lines(density(dataReal.df$weighted_meanDifference, na.rm = TRUE), col = "red")

    for (TFCur in unique(output.global.TFs.orig$TF)) {
        
        dataCur.df = dplyr::filter(output.global.TFs.orig, TF == TFCur)
        dataReal.df = dplyr::filter(dataCur.df, permutation == 0)
        dataPerm.df = dplyr::filter(dataCur.df, permutation > 0)
        
        
        if (nrow(dataPerm.df) == 0) {
            message = paste0("For TF : ", TFCur, ", no permutation data could be found. This is not supposed to happen. Rerun the binning step for this TF.")
            checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
            next
        }
        
        rangeX = range(dataCur.df$weighted_meanDifference)
        rowCur = which(output.global.TFs.orig$TF == TFCur)
        
        g = ggplot(dataPerm.df, aes(weighted_meanDifference)) + geom_density() + geom_vline(xintercept = dataReal.df$weighted_meanDifference[1], color = "red") + ggtitle(TFCur) + xlim(c(rangeX * 1.5)) # + scale_x_continuous(limits = c(min(dataCur.df$weighted_meanDifference) - 0.5, max(dataCur.df$weighted_meanDifference) + 0.5))
        
        diagPlots.l[[TFCur]] = g
        plot(g)

        nPermThreshold = length(which(abs(dataPerm.df$weighted_meanDifference) > abs(dataReal.df$weighted_meanDifference[1])))
        pvalueCur = nPermThreshold / par.l$nPermutations
        output.global.TFs.orig$pvalue[rowCur] = pvalueCur 

    }
    
    # TODO: Diagnostic plot for pvalue
    plot(ggplot(output.global.TFs.orig, aes(pvalue)) + geom_density() + ggtitle("Local fdr density across all TF"))

} else {

  # # Calculate adjusted p values out of variance 
  # estimates = tryCatch( {
  #   
  #   locfdrRes = locfdr(output.global.TFs.orig$weighted_Tstat, plot = 4)
  #   # Currently taken as default
  #   MLE.delta  = locfdrRes$fp0["mlest", "delta"]
  #   
  #   # Not the current default
  #   CME.delta  = locfdrRes$fp0["cmest", "delta"]
  #   
  #   c(MLE.delta, CME.delta)
  #   
  # }, error = function(e) {
  #   message = "Could not run locfdr, use the median instead..."
  #   checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  #   
  #   # For the fallback, simply return the median for both estimates for now
  #   c(median(output.global.TFs.orig$weighted_Tstat, na.rm = TRUE), median(output.global.TFs.orig$weighted_Tstat, na.rm = TRUE))
  # }
  # )
  # 
  # MLE.delta = estimates[1]
  # CME.delta = estimates[2]
  # 
  # # Compute two different measures for the mean
  # # Which one to take? Depends, see page 101 in Bradley Efron-Large-Scale Inference_ Empirical Bayes Methods for Estimation, Testing, and Prediction
  # # The default option in locfdr is the MLE method, not central matching. Slight irregularities in the central histogram,
  # # as seen in Figure 6.1a, can derail central matching. The MLE method is more stable, but pays the price of possibly increased bias.

  # output.global.TFs.orig$weighted_Tstat_centralized = output.global.TFs.orig$weighted_Tstat - MLE.delta

  populationMean = 0
  zScore = (output.global.TFs.orig$weighted_Tstat - populationMean) / sqrt(output.global.TFs.orig$variance)
  
  # 2-sided test
  output.global.TFs.orig$pvalue   = 2*pnorm(-abs(zScore))
  
  # Handle extreme cases with p-values that are practically 0 and would cause subsequent issues
  index0 = which(output.global.TFs.orig$pvalue < .Machine$double.xmin)
  if (length(index0) > 0) {
    output.global.TFs.orig$pvalue[index0] = .Machine$double.xmin
  }
}

output.global.TFs.permutations = dplyr::filter(output.global.TFs.orig, permutation > 0)
output.global.TFs              = dplyr::filter(output.global.TFs.orig, permutation == 0)



output.global.TFs = mutate(output.global.TFs, 
                           Cohend_factor = ifelse(weighted_CD < par.l$thresholds_CohensD[1], par.l$classes_CohensD[1], 
                                                  ifelse(weighted_CD < par.l$thresholds_CohensD[2] , par.l$classes_CohensD[2], 
                                                         ifelse(weighted_CD < par.l$thresholds_CohensD[3], par.l$classes_CohensD[3], par.l$classes_CohensD[4]))),
                           Cohend_factor = factor(Cohend_factor, levels = par.l$classes_CohensD, labels = seq_len(length(par.l$classes_CohensD))),
                           pvalueAdj     = p.adjust(pvalue, method = "BH"))


colnamesToPlot = c("weighted_meanDifference", "weighted_CD", "TFBS", "weighted_Tstat", "variance", "pvalue", "pvalueAdj")

for (pValueCur in c(par.l$significanceThresholds , 1)) {
  
  filtered.df = dplyr::filter(output.global.TFs, pvalueAdj <= pValueCur)
  
  title = paste0("p-value: ", pValueCur, " (retaining ", nrow(filtered.df), " TF)")
  
  for (measureCur in colnamesToPlot) {
      
      if (all(!is.finite(unlist(filtered.df[,measureCur])))) {
          next
      }
      
      if (! measureCur %in% colnames(output.global.TFs)) {
          next
      }
    
    if (measureCur %in%  c("Cohend_factor")) {
      plot(ggplot(filtered.df, aes_string(measureCur))  + stat_count() + theme_bw() + ggtitle(title))
      
    } else {
      plot(ggplot(filtered.df, aes_string(measureCur))  + geom_histogram(bins = 50) + theme_bw() + ggtitle(title))
    }
    
  }
}

stats.df = group_by(output.global.TFs.orig, permutation) %>% summarise(max = max(weighted_meanDifference), min = min(weighted_meanDifference))
ggplot(stats.df, aes(min)) + geom_density()
ggplot(stats.df, aes(max)) + geom_density()
dev.off()

#################################
# mode of change quantification #
#################################

sampleData.l = readRDS(par.l$file_input_metadata)
sampleData.df = sampleData.l[["permutation0"]]
designComponents.l = checkDesignIntegrity(snakemake, par.l, sampleData.df, useRNA = TRUE)

components3types   = designComponents.l$types

# Which of the two modes should be done, pairwise or quantitative?
comparisonMode = "quantitative"
if (components3types["conditionSummary"] == "logical" | components3types["conditionSummary"] == "factor") {
    comparisonMode = "pairwise"
}


##########################
# INTEGRATE RNA-Seq DATA #
##########################
if (par.l$plotRNASeqClassification) {
  
    # Require some more packages here
    checkAndLoadPackages(c( "lsr", "DESeq2",  "matrixStats",  "pheatmap", "preprocessCore"), verbose = FALSE)
    
    
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
    
    # Read RNAseq counts
    TF.counts.df.all = read_tsv(par.l$file_input_geneCountsPerSample, col_names = TRUE)
    if (nrow(problems(TF.counts.df.all)) > 0) {
        flog.fatal(paste0("Parsing errors: "), problems(TF.counts.df.all), capture = TRUE)
        stop("Error when parsing the file ", par.l$file_input_geneCountsPerSample, ", see errors above")
    }
    
    # Add row names
    TF.counts.df.all = as.data.frame(TF.counts.df.all)
    rownames(TF.counts.df.all) = TF.counts.df.all$ENSEMBL
    TF.counts.df.all = TF.counts.df.all[,-1]
    
    # Clean ENSEMBL IDs
    rownames(TF.counts.df.all) = gsub("\\..+", "", rownames(TF.counts.df.all), perl = TRUE)

    
    sampleData.df = dplyr::filter(sampleData.df, SampleID %in% colnames(TF.counts.df.all))
    
    nFiltRows = nrow(sampleData.l[["permutation0"]]) - nrow(sampleData.df)
    if (nFiltRows > 0) {
        flog.warn(paste0("Filtered ", nFiltRows, " sample IDs after comparising sample names with RNA-Seq table"))
    }
    
    # The design formula for RNA-Seq is different from the one we used before for ATAC-Seq
    # Either take the one that the user provided or, if he did not, use a general one with only the condition
    par.l$designFormulaRNA = snakemake@config$par_general$designContrastRNA
    if (is.null(par.l$designFormulaRNA)) {
        par.l$designFormulaRNA = "~conditionSummary"
        flog.warn(paste0("Could not find the parameter designContrastRNA in the configuration file. The default of \"~conditionSummary\" will be taken as formula. If you know about confounding variables, rerun this step and add the parameter (see the Documentation for details)"))
    }
    designFormulaRNA = convertToFormula(par.l$designFormulaRNA, colnames(sampleData.df))
    
    file_sampleSummary = snakemake@config$samples$summaryFile
    assertFileExists(file_sampleSummary)

    sampleSummary.df = read_tsv(file_sampleSummary, col_types = cols())
    
    # Loading TF data
    HOCOMOCO_mapping.df = read.table(file = par.l$file_input_HOCOMOCO_mapping, header = TRUE)
    assertSubset(c("ENSEMBL", "HOCOID"), colnames(HOCOMOCO_mapping.df))
    
    # Clean ENSEMBL IDs
    HOCOMOCO_mapping.df$ENSEMBL = gsub("\\..+", "", HOCOMOCO_mapping.df$ENSEMBL, perl = TRUE)
    

    ####################################
    # Run DeSeq2 on raw RNA-Seq counts #
    ####################################
    dd <- DESeqDataSetFromMatrix(countData = TF.counts.df.all[,sampleData.df$SampleID],
                                 colData = sampleData.df,
                                 design = designFormulaRNA)
    
    dd = estimateSizeFactors(dd)
    # dd = DESeq(dd)
    dd_counts =  counts(dd, normalized=TRUE)

    ######################################
    # Filtering of lowly expressed genes #
    ######################################
    nMin = par.l$filter_minCountsPerCondition
    if (comparisonMode == "pairwise") {
        
        # If in either of the two conditions counts are below a minimum OR the median is 0, we discard the gene and mark them as non-expressed
        samples_cond1 = colData(dd)$SampleID[which(colData(dd)$conditionSummary == levels(colData(dd)$conditionSummary)[1])]
        samples_cond2 = colData(dd)$SampleID[which(colData(dd)$conditionSummary == levels(colData(dd)$conditionSummary)[2])]
        
        idx <- (rowMeans(dd_counts[,samples_cond1]) > nMin | rowMeans(dd_counts[,samples_cond2]) > nMin) & rowMedians(dd_counts) > 0
        
    } else {
        
        # Here, we have to employ a different approach for when to filter genes given that we may have more than two conditions
        idx <- rowMedians(dd_counts) > 0
        idx <- rowMeans(dd_counts) > 0
    }

    nFiltered = length(which(idx == FALSE))
    
    if (nFiltered > 0) {
        flog.info(paste0("Filtered ", nFiltered, " genes from RNA-Seq table because of low counts."))
        dd.filt = dd[idx,]
    } else {
        dd.filt = dd
    }

    dd.filt <- DESeq(dd.filt)
    dd_counts.filt =  counts(dd.filt, normalized=TRUE)
    RNA.counts.filt.df = dd_counts.filt %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% as.tibble()
    
    # Raw counts, used for other types of normalization thereafter
    dd_counts.raw.filt =  counts(dd.filt, normalized=FALSE)
   
    dd_counts.filt.quantile = normalize.quantiles(as.matrix(dd_counts.raw.filt))
    RNA.counts.quantile.df.all = dd_counts.filt.quantile %>% as.data.frame()  %>% as.tibble()
    
    # Fix row and column names
    colnames(RNA.counts.quantile.df.all) = colnames(dd_counts.raw.filt)
    RNA.counts.quantile.df.all = RNA.counts.quantile.df.all %>%
        mutate(ENSEMBL = RNA.counts.filt.df$ENSEMBL) %>%
        dplyr::select(ENSEMBL, colnames(dd_counts.raw.filt))
    
    # TODO: Change and make permanent
    RNA.counts.filt.df = RNA.counts.quantile.df.all
    
    # Check sample names and set column names
    # Match the column names and do the intersections
    sharedColumns = intersect(colnames(RNA.counts.filt.df)[-1], sampleSummary.df$SampleID)
    
    if (length(sharedColumns) == 0) {
      message = paste0("No shared samples with RNA-Seq samples between sample table ", file_sampleSummary, " and RNA-Seq table ", par.l$file_input_geneCountsPerSample, ".")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    } 
    
    if (length(sharedColumns) < nrow(sampleSummary.df)) {
        message = paste0("Only ", length(sharedColumns), " out of ", ncol(RNA.counts.filt.df) - 1, " columns are shared between RNA-Seq counts and the sample table. Make sure the columns names in the RNA-seq counts file are identical to the names in the sample table.")
        checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    } else {
      flog.info(paste0(length(sharedColumns), " out of ", ncol(RNA.counts.filt.df) - 1, " columns are shared between RNA-Seq counts and the sample table"))
    }
    
    colnames(RNA.counts.filt.df)[1] = "ENSEMBL"


    # Filter them by the IDs that correspond to the TFs
    TF.counts.filt.df = dplyr::filter(RNA.counts.filt.df, ENSEMBL %in% HOCOMOCO_mapping.df$ENSEMBL)
    
    if (nrow(TF.counts.filt.df) == 0) {
      message = "No rows remaining after filtering against ENSEMBL IDs in HOCOMOCO. Check your ENSEMBL IDs for overlap with the HOCOMOCO translation table."
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    
    # Read counts per TFBS, coming from previous steps
    peak.counts = read_tsv(par.l$file_input_countsNorm, col_types = cols())
    if (nrow(problems(peak.counts)) > 0) {
        flog.fatal(paste0("Parsing errors: "), problems(peak.counts), capture = TRUE)
        stop("Error when parsing the file ", par.l$file_input_countsNorm, ", see errors above")
    }
    
    # Match the column names and do the intersections
    sharedColumns = intersect(colnames(peak.counts), colnames(TF.counts.filt.df))
    
    if (length(sharedColumns) == 0) {
      message = "No shared samples with RNA-Seq samples."
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    flog.info(paste0(length(sharedColumns), " samples are shared between the input data and RNA-Seq data"))
    
    #peak.counts.orig = peak.counts
    peak.counts  = dplyr::select(peak.counts, one_of("peakID", sharedColumns)) 
    TF.counts.filt.df = TF.counts.filt.df[, which(colnames(TF.counts.filt.df) %in% c(sharedColumns, "ENSEMBL"))]
    
    
    HOCOMOCO_mapping.df.exp <- dplyr::filter(HOCOMOCO_mapping.df, ENSEMBL %in%  TF.counts.filt.df$ENSEMBL, HOCOID %in% output.global.TFs$TF)
    
    # Filter genes that might not be in the output.global.TFs$TF list 
    # Dont include, cuases errors downstream. Fitler this in the heatmap function
    # TF.counts.filt.df = filter(TF.counts.filt.df, ENSEMBL %in% HOCOMOCO_mapping.df.exp$ENSEMBL)
    
    if (nrow(HOCOMOCO_mapping.df.exp) == 0) {
      message = paste0("Number of rows of HOCOMOCO_mapping.df.exp is 0. Something is wrong with the mapping table or the filtering")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    # Loop through all TFs and determine whether or not a peak has a binding site for this TF, resulting in a binary matrix
    TF.peakMatrix.l = list()
    for (TFCur in HOCOMOCO_mapping.df.exp$HOCOID) {
      
      HOCOMOCO_mapping.subset.df = subset(HOCOMOCO_mapping.df.exp, HOCOID == TFCur)
      gene.sel = unique(HOCOMOCO_mapping.subset.df$ENSEMBL)
      if (length(gene.sel) > 1) {
        message = paste0("Mapping for ", TFCur, " not unique, take only the first mapping (", gene.sel[1], ") and discard the others (", paste0(gene.sel[-1], collapse = ","), ")")
        checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        
      }
      
      TF.output.df = read.table(file = paste0(rootOutdir, "/TF-SPECIFIC/",TFCur,"/extension", extensionSize, "/", comparisonType, TFCur,  ".output.tsv.gz"), header = TRUE)
      TF.peakMatrix.l[[TFCur]] = peak.counts$peakID %in% TF.output.df$peakID
    }
    
    # cor.m = peaks names separated with # 
    
    # This is the peak (rows) and TF binding sites (columns)
    TF.peakMatrix.df = as.data.frame(TF.peakMatrix.l)
    
    # Sanity check
    if (all(rowSums(TF.peakMatrix.df) == 0)) {
      message = paste0("All counts are 0, something is wrong.")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }

    # Expressed TFs only 
    # TODO: Filter for expressed TFs, how is this done and WHY
    HOCOMOCO_mapping.subset.df = subset(HOCOMOCO_mapping.df.exp, HOCOID %in% colnames(TF.peakMatrix.df))
    
    # Remove first column, retain only counts
    expressed.TF.counts.df = t(TF.counts.filt.df[,-c(1)])
    colnames(expressed.TF.counts.df) = TF.counts.filt.df$ENSEMBL
    expressed.TF.counts.df = t(expressed.TF.counts.df)
    
    assertSubset(colnames(expressed.TF.counts.df), colnames(peak.counts))
    
    # Some rownames may be identical because of the mapping from HOCOMOCO
    
    peak.counts = peak.counts[,order(colnames(peak.counts))]
    expressed.TF.counts.df = expressed.TF.counts.df[,order(colnames(expressed.TF.counts.df))]
    
    index_lastColumn = which(colnames(peak.counts) == "peakID")
    
    if (!all(colnames(expressed.TF.counts.df) == colnames(peak.counts)[-index_lastColumn])) {
      message = "Colnames of expressed.TF.counts.df and peak.counts must be identical"
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }

    rowMeans2 = rowMeans(peak.counts[,-index_lastColumn])
    rowsToDelete = which(rowMeans2 < 1)
    if (length(rowsToDelete) > 0) {
        flog.info(paste0("Removed ", length(rowsToDelete), " peaks out of ", nrow(peak.counts), " because they had a row mean of < 1."))
    
        # Filter these peaks also from the peakCount matrix  
        stopifnot(nrow(TF.peakMatrix.df) == nrow(peak.counts))
        TF.peakMatrix.df = TF.peakMatrix.df[-rowsToDelete,]
        peak.counts = peak.counts[-rowsToDelete,]
    }
    
    peak.counts = dplyr::select(peak.counts, -one_of("peakID"))
    
    # Same dimensions as sel.TF.peakMatrix.df:  peaks as rows and TFs as columns. Stores correlations of peak counts vs RNA-Seq counts for each TF-Gene / peak
    cor.m = t(cor(t(expressed.TF.counts.df), t(peak.counts), method = par.l$corMethod))
    
    # Mapping TFBS to TF 
    sort.cor.m = cor.m[,names(sort(colMeans(cor.m)))] 
    
    HOCOMOCO_mapping.df.exp = HOCOMOCO_mapping.df.exp[which(HOCOMOCO_mapping.df.exp$ENSEMBL %in% colnames(sort.cor.m)),]
    
    # Some entries in the HOCOMOCO mapping can be repeated (i.e., the same ID for two different TFs, such as ZBTB4.S and ZBTB4.D)
    # Originally, we deleted these rows from the mapping and took the first entry only
    # However, since TFs with the same ENSEMBL ID can still be different with respect to their TFBS, we now duplicate such genes also in the correlation table
    #HOCOMOCO_mapping.df.exp = HOCOMOCO_mapping.df.exp[!duplicated(HOCOMOCO_mapping.df.exp[, c("ENSEMBL")]),]
    #assertSubset(as.character(HOCOMOCO_mapping.df.exp$ENSEMBL), colnames(sort.cor.m))
    
    
    # Change the column names from ENSEMBL ID to TF names. Reorder the columns first to make sure the order is the same. Due to the duplication ID issue, the number of columns may increase after the column selection below
    sort.cor.m = sort.cor.m[,as.character(HOCOMOCO_mapping.df.exp$ENSEMBL)] 
    colnames(sort.cor.m) = as.character(HOCOMOCO_mapping.df.exp$HOCOID)
    
    # This binary matrix has peaks as rows and TFs as columns of whether or not a particular peak has a TFBS from this TF or not
    sel.TF.peakMatrix.df = TF.peakMatrix.df[,colnames(sort.cor.m)]

    
    # 1. Focus on peaks with TFBS overlaps
    t.cor.sel.matrix = sort.cor.m
    # Transform 0 values into NA for the matrix to speed up subsequent analysis
    t.cor.sel.matrix[(sel.TF.peakMatrix.df == 0)] = NA
    # Goal: Eliminate all correlation values for cases in which a peak has no TFBS for the particular TF
    # Same dimensions as the two matrices used for input.
    # Matrix multiplication here essentially means we only multiply each individual entry with either 1 or NA
    # The result is a matrix full of NAs and the remaining entries are the correlation values for peaks that have a TFBS for the particular TF
    t.cor.sel.matrix = sel.TF.peakMatrix.df * t.cor.sel.matrix
    # Gives one value per TF, designating the median correlation per TF across all peaks
    median.cor.tfs = sort(apply(t.cor.sel.matrix, MARGIN = 2, FUN = my.median))
    
    # 2. Background
    # Start with the same correlation matrix
    t.cor.sel.matrix.non = sort.cor.m
    # Transform 1 values into NA for the matrix to speed up subsequent analysis
    t.cor.sel.matrix.non[(sel.TF.peakMatrix.df == 1)] = NA
    t.cor.sel.matrix.non = sel.TF.peakMatrix.df + t.cor.sel.matrix.non
    # Gives one value per TF, designating the median correlation per TF
    median.cor.tfs.non <- sort(apply(t.cor.sel.matrix.non, MARGIN=2, FUN = my.median))
    # Not used thereafter
    median.cor.tfs.rest <- sort(median.cor.tfs - median.cor.tfs.non[names(median.cor.tfs)])
    
    # 3. Final classification: Calculate thresholds by calculating the quantiles of the background anhd compare the real values to the background
    act.rep.thres.l = list()
    for (thresholdCur in c(0.1, 0.05, 0.01, 0.001)) {
        
        act.rep.cur = quantile(sort(apply(t.cor.sel.matrix.non, MARGIN = 2, FUN = my.median)), probs = c(thresholdCur, 1-thresholdCur))
        # Enforce the thresholds to be at least 0, so we never have an activator despite a negative median correlation 
        # and a repressor despite a positive one
        act.rep.thres.l[[as.character(thresholdCur)]][1] = min(0, act.rep.cur[1])
        act.rep.thres.l[[as.character(thresholdCur)]][2] = max(0, act.rep.cur[2])

        flog.info(paste0("Thresholds for repressor/activator for threshold ", thresholdCur, ": ", act.rep.thres.l[[as.character(thresholdCur)]][1], " and ", act.rep.thres.l[[as.character(thresholdCur)]][2]))
        
    }
   
    # AR.data = as.data.frame(median.cor.tfs)
    # AR.data$TF = rownames(AR.data)
    # 
    # output.global.TFs = merge(output.global.TFs, AR.data, by = "TF",all.x = TRUE)
    # 
    # output.global.TFs$classification = ifelse(is.na(output.global.TFs$median.cor.tfs), "not-expressed",
    #                                         ifelse(output.global.TFs$median.cor.tfs <= act.rep.thres[1], "repressor",
    #                                                ifelse(output.global.TFs$median.cor.tfs > act.rep.thres[2], "activator", "undetermined")))
    # 
    # output.global.TFs$classification = factor(output.global.TFs$classification, levels = names(par.l$colorCategories))
    # 
    # 
    # 
    # 
    
    # TODO: for each TFBS, a p-value and a correlation value
    colnameClassification    = paste0("classification")
    colnameMedianCor         = paste0("median.cor.tfs")
    colnameClassificationPVal= paste0("classification_distr_rawP")
    
    AR.data = as.data.frame(median.cor.tfs)
    AR.data$TF = rownames(AR.data)
    colnames(AR.data)[1] = colnameMedianCor
    
    output.global.TFs[,colnameMedianCor] = NULL
    output.global.TFs = merge(output.global.TFs, AR.data, by = "TF",all.x = TRUE)
    
    
    # Merge new ones into old
    for (thresCur in names(act.rep.thres.l)) {
        thresCur.v = act.rep.thres.l[[thresCur]]
        colnameClassificationCur = paste0(colnameClassification, "_q", thresCur)
        output.global.TFs[, colnameClassificationCur] = ifelse(is.na(output.global.TFs[,colnameMedianCor]), "not-expressed",
                                                               ifelse(output.global.TFs[,colnameMedianCor] <= thresCur.v[1], "repressor",
                                                                      ifelse(output.global.TFs[,colnameMedianCor] >= thresCur.v[2], "activator", "undetermined")))
        output.global.TFs[, colnameClassificationCur] = factor(output.global.TFs[, colnameClassificationCur], levels = names(par.l$colorCategories))
        
    }
    
    
    output.global.TFs[,colnameClassificationPVal] = NULL
    # Do a Wilcoxon test for each TF as a 2nd filtering criterion
    for (TFCur in colnames(t.cor.sel.matrix)) {
        
        rowNo = which(output.global.TFs$TF == TFCur)
        
        # Removing NAs actually makes a difference, as these are "artifical" anyway here due to the two matrices let's remove them
        dataMotif      = na.omit(t.cor.sel.matrix[,TFCur])
        dataBackground = na.omit(t.cor.sel.matrix.non[,TFCur])
        
        # Test the distributions
        if (output.global.TFs[rowNo, colnameMedianCor] > 0) {
            alternativeTest = "greater"
        } else {
            alternativeTest = "less"
        }
        
        testResults = wilcox.test(dataMotif, dataBackground, alternative = alternativeTest)
        
        stopifnot(length(rowNo) == 1)
        output.global.TFs[rowNo,colnameClassificationPVal] = testResults$p.value
        
    }
    
    ################################################
    # POST-FILTER: CHANGE SOME TFs TO UNDETERMINED #
    ################################################
    
    par.l$thresholds_pvalue_Wilcoxon = 0.05
    # Change the classification with the p-value from the distribution test
    
    colnameClassificationPVal = paste0("classification_distr_rawP")
    
    for (thresholdCur in c(0.1, 0.05, 0.01, 0.001)) {
        
        flog.info(paste0("Post filter: Doing Wilcoxon test for threshold ", thresholdCur))
        
        colnameClassification      = paste0("classification_q", thresholdCur)
        colnameClassificationFinal = paste0("classification_q", thresholdCur, "_final")
        
        output.global.TFs[,colnameClassificationFinal] = output.global.TFs[,colnameClassification]
        
        TFs_to_change = dplyr::filter(output.global.TFs, (!!as.name(colnameClassification) == "activator" | !!as.name(colnameClassification) == "repressor") & 
                                   !!as.name( colnameClassificationPVal) > !!par.l$thresholds_pvalue_Wilcoxon)$TF
        
        # Filter some TFs to be undetermined
        if (length(TFs_to_change) > 0) {
            flog.info(paste0(" Changing the following TFs to 'undetermined' because they were classified as either activator or repressor before but the Wilcoxon test was not significant: ", paste0(TFs_to_change, collapse = ",")))

            output.global.TFs[which(output.global.TFs$TF  %in% TFs_to_change), colnameClassificationFinal] = "undetermined"
        }

        
    }
        
    
    
    
    
    
    ####################
    ####################
    # DIAGNOSTIC PLOTS #
    ####################
    ####################

    pdf(file = par.l$files_plotDiagnostic[2], width = 3, height = 8)
    xlab="median pearson correlation (r)"
    ylab=""
    xlim= c(-max(abs(range(median.cor.tfs))) - 0.05, max(abs(range(median.cor.tfs))) + 0.05)
    ylim=c(1,length(median.cor.tfs.non))
    
    par(mfrow=c(1,1))
    
    for (thresCur in names(act.rep.thres.l)) {
        thresCur.v = act.rep.thres.l[[thresCur]]
        
        thresCur_upper = (1 - as.numeric(thresCur)) * 100
        thresCur_lower = as.numeric(thresCur) * 100
        
        
        plot(median.cor.tfs.non[names(median.cor.tfs)], 1:length(median.cor.tfs.non),
             xlim=xlim, ylim=ylim, main="", xlab=xlab, ylab=ylab,
             col=adjustcolor("darkgrey",alpha=1), pch = 16, cex = 0.5,axes = FALSE)
        points(median.cor.tfs, 1:length(median.cor.tfs),
               pch=16,  cex=0.5, 
               col=ifelse(median.cor.tfs>thresCur.v[2], par.l$colorCategories["activator"] ,ifelse(median.cor.tfs<thresCur.v[1], par.l$colorCategories["repressor"], par.l$colorCategories["undetermined"]))
        ) 
        text(x =c((thresCur.v[1]-0.01),(thresCur.v[2]+0.01)),
             y=c((length(median.cor.tfs.non)+5),(length(median.cor.tfs.non)+5)), pos=c(2,4),
             labels =c(paste0(thresCur_lower, "\npercentile"), paste0(thresCur_upper, "\npercentile")),cex=0.7, col=c("black","black"))
        abline(v=thresCur.v[1], col=par.l$colorCategories["repressor"])
        abline(v=thresCur.v[2], col=par.l$colorCategories["activator"])
        axis(side = 1, lwd = 1, line = 0, at = c(-0.2,0,0.2), cex=1)
        
    }
    
    
    
    heatmap.act.rep(TF.peakMatrix.df, HOCOMOCO_mapping.df.exp, cor.m, par.l, median.cor.tfs, median.cor.tfs.non, act.rep.thres.l)

    dev.off()


    # Process results
    res.peaks.filt  = results(dd.filt) %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% as.tibble()

  
    # TODO
    # In ivans version, he used a non-filtered HOCOMOCO table. I however filter it before, so that some ENSEMBL IDs might already be filtered
    
    TF.specific = left_join(HOCOMOCO_mapping.subset.df, res.peaks.filt, by = "ENSEMBL") %>% dplyr::filter(!is.na(baseMean))
    
    if (nrow(TF.specific) == 0) {
        message = "The Ensembl IDs from the translation table do not match with the IDs from the RNA-seq counts table."
        checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }

    # TODO:_ deal with NAs, where do they come from?

    output.global.TFs$weighted_meanDifference = as.numeric(output.global.TFs$weighted_meanDifference)
    
    pdf(par.l$files_plotDiagnostic[3])
    thresholds  = c(0.1, 0.05, 0.01, 0.001)
    
    #######################################
    # Correlation plots for the 3 classes #
    #######################################
    for (thresholdCur in thresholds) {
        
        colnameClassificationCur = paste0("classification_q", thresholdCur, "_final")
        
        output.global.TFs.merged = output.global.TFs %>%
          dplyr::filter(!!as.name(colnameClassificationCur) != "not-expressed")  %>%
          full_join(TF.specific, by = c( "TF" = "HOCOID"))  %>%
          mutate(baseMeanNorm = (baseMean - min(baseMean, na.rm = TRUE)) / (max(baseMean, na.rm = TRUE) - min(baseMean, na.rm = TRUE)) + par.l$minPointSize)   %>%
          dplyr::filter(!is.na(!!as.name(colnameClassificationCur))) 

       
        for (classificationCur in unique(output.global.TFs.merged[,colnameClassificationCur])) {
          
          output.global.TFs.cur = dplyr::filter(output.global.TFs.merged, !!as.name(colnameClassificationCur) == classificationCur)
          
          cor.res.l = list()
          for (corMethodCur in c("pearson", "spearman")) {
            cor.res.l[[corMethodCur]] = cor.test(output.global.TFs.cur$weighted_meanDifference, output.global.TFs.cur$log2FoldChange, method = corMethodCur)
          }
          
          titleCur = paste0(classificationCur, ": R=", 
                            signif(cor.res.l[["pearson"]]$estimate, 2), "/", 
                            signif(cor.res.l[["spearman"]]$estimate, 2), ", p-value ", 
                            signif(cor.res.l[["pearson"]]$p.value,2),  "/", 
                            signif(cor.res.l[["spearman"]]$p.value,2), "\n(Pearson/Spearman, stringency: ", thresholdCur, ")")
          
          g = ggplot(output.global.TFs.cur, aes(weighted_meanDifference, log2FoldChange)) + geom_point(aes(size = baseMeanNorm)) + 
            geom_smooth(method = par.l$regressionMethod, color = par.l$colorCategories[classificationCur]) + 
            ggtitle(titleCur) + 
            ylab("log2 fold-change RNA-seq") + 
            theme_bw() + theme(plot.title = element_text(hjust = 0.5))
          plot(g)
          
        }
    }

    #############################
    # Density plots for each TF #
    #############################
    
    stopifnot(identical(colnames(t.cor.sel.matrix), colnames(t.cor.sel.matrix.non)))
    
    for (colCur in seq_len(ncol(t.cor.sel.matrix))) {
      
      TFCur = colnames(t.cor.sel.matrix)[colCur]
      dataMotif      = t.cor.sel.matrix[,colCur]
      dataBackground = t.cor.sel.matrix.non[,colCur]
      mainLabel = paste0(TFCur," (#TFBS = ",length(which(!is.na(dataMotif)))," )")
      
      plot(density(dataMotif, bw=0.1, na.rm=TRUE), xlim=c(-1,1), ylim=c(0,2),
           main=, mainLabel, lwd=2.5, col="red", axes = FALSE, xlab = "Pearson correlation")
      abline(v=0, col="black", lty=2)
      legend("topleft",box.col = adjustcolor("white",alpha.f = 0), legend = c("Motif","Non-motif"), lwd=c(2,2),cex = 0.8, col=c("red","darkgrey"), lty=c(1,1) )
      axis(side = 1, lwd = 1, line = 0)
      axis(side = 2, lwd = 1, line = 0, las = 1)
      lines(density(dataBackground, bw=0.1, na.rm=T), lwd=2.5, col="darkgrey")
      
    } 
    dev.off()

} else {
  classesList.l = list(c())
}

output.global.TFs$yValue = transform_yValues(output.global.TFs$pvalueAdj, addPseudoCount = TRUE, nPermutations = par.l$nPermutations)
output.global.TFs.origReal = output.global.TFs

#########################################
# PLOT FOR DIFFERENT P VALUE THRESHOLDS #
#########################################

# Set the page dimensions to the maximum across all plotted variants
output.global.TFs.filteredSummary = dplyr::filter(output.global.TFs, pvalue <= max(par.l$significanceThresholds))

nTF_label = min(par.l$maxTFsToLabel, nrow(output.global.TFs.filteredSummary))

TFLabelSize = ifelse(nTF_label < 20, 8,
                     ifelse(nTF_label < 40, 7,
                      ifelse(nTF_label < 60, 6,
                        ifelse(nTF_label < 60, 5,
                          ifelse(nTF_label < 60, 4, 3)))))



variableXAxis = "weighted_meanDifference"

thresholds = "NA" # not applicable here
if (par.l$plotRNASeqClassification) {
    thresholds  = c(0.1, 0.05, 0.01, 0.001)
}

allPlots.l = list()

for (thresholdCur in thresholds) {
    
    flog.info(paste0("Stringency threshold for classification: ", thresholdCur))
    
    allPlots.l[[as.character(thresholdCur)]] = list()
    colnameClassificationFinal = paste0("classification_q", thresholdCur, "_final")

    for (significanceThresholdCur in par.l$significanceThresholds) {
      
      pValThrStr = as.character(significanceThresholdCur)
    
      for (showClasses in classesList.l) {
    
        output.global.TFs = output.global.TFs.origReal %>%
            mutate( pValueAdj_log10 = transform_yValues(pvalueAdj, addPseudoCount = TRUE, nPermutations = par.l$nPermutations),
                    pValue_log10 = transform_yValues(pvalue, addPseudoCount = TRUE, nPermutations = par.l$nPermutations),
                    pValueAdj_sig = pvalueAdj <= significanceThresholdCur,
                    pValue_sig = pvalue <= significanceThresholdCur) 
           
        if (par.l$plotRNASeqClassification) {
            #output.global.TFs = filter(output.global.TFs, classification %in% showClasses)
             output.global.TFs = dplyr::filter(output.global.TFs, !!as.name(colnameClassificationFinal) %in% showClasses)
    
        }
        
        for (pValueStrCur in c("pvalue", "pvalueAdj")) {
            
            if (pValueStrCur == "pvalue") {
                
                pValueScoreCur = "pValue_log10"
                pValueSigCur = "pValue_sig"
                pValueStrLabel = "raw p-value"
                
                ggrepel_df = dplyr::filter(output.global.TFs, pValue_sig == TRUE)
                maxPValue = max(output.global.TFs$pValue_log10, na.rm = TRUE)
                
            } else {
                
                pValueScoreCur = "pValueAdj_log10"
                pValueSigCur = "pValueAdj_sig"
                pValueStrLabel = "adj. p-value"
                
                ggrepel_df = dplyr::filter(output.global.TFs, pValueAdj_sig == TRUE)
                maxPValue = max(output.global.TFs$pValueAdj_log10, na.rm = TRUE)
            }
        
            # Increase the ymax a bit more
            ymax = max(transform_yValues(significanceThresholdCur, addPseudoCount = TRUE, nPermutations = par.l$nPermutations), maxPValue, na.rm = TRUE) * 1.1
            alphaValueNonSign = 0.3
            
            # Reverse here because negative values at left mean that the condition that has been specified in the beginning is higher. 
            # Reverse the rev() that was done before for this plot therefore to restore the original order
            labelsConditionsNew = rev(conditionComparison)
        
            g = ggplot()
            
            label_nameChange = 'TF activity higher in'
            if (comparisonMode != "pairwise") {
                label_nameChange = 'Change with increasing\nvalues of conditionSummary'
            }
            
            if (par.l$plotRNASeqClassification) {
            # TODO: finalize
              # g = g + geom_point(data = output.global.TFs, aes_string("weighted_meanDifference", pValueScoreCur, 
              #        alpha = pValueSigCur, size = "TFBS", fill = "classification"),  shape=21, stroke = 0.5, color = "black") +  
              #     scale_fill_manual("TF class", values = par.l$colorCategories)
             labelLegendCur = paste0("TF class (stringency: ", thresholdCur, ")")
              g = g + geom_point(data = output.global.TFs, aes_string("weighted_meanDifference", pValueScoreCur, 
                     alpha = pValueSigCur, size = "TFBS", fill = colnameClassificationFinal),  shape=21, stroke = 0.5, color = "black") +  
                  scale_fill_manual(labelLegendCur, values = par.l$colorCategories)
              
              g = g + geom_rect(aes(xmin = -Inf,xmax = 0,ymin = -Inf, ymax = Inf, color = par.l$colorConditions[2]),
                            alpha = .3, fill = par.l$colorConditions[2], size = 0) +
                      geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf,ymax = Inf, color = par.l$colorConditions[1]), alpha = .3, fill = par.l$colorConditions[1], size = 0) + 
                      scale_color_manual(name = label_nameChange, values = par.l$colorConditions, labels = conditionComparison)
              
            } else {
                
              g = g + geom_point(data = output.global.TFs, aes_string("weighted_meanDifference", pValueScoreCur, alpha = pValueSigCur, size = "TFBS"), shape=21, stroke = 0.5, color = "black")
              g = g + geom_rect(aes(xmin = -Inf, xmax = 0,   ymin = -Inf, ymax = Inf, fill = par.l$colorConditions[2]), alpha = .3) + 
                      geom_rect(aes(xmin = 0,    xmax = Inf, ymin = -Inf, ymax = Inf, fill = par.l$colorConditions[1]), alpha = .3)
              g = g + scale_fill_manual(name = 'TF activity higher in', values = rev(par.l$colorConditions), labels = labelsConditionsNew)
            }
         
            g = g + ylim(-0.1,ymax) + 
                ylab(paste0(transform_yValues_caption(), " (", pValueStrLabel, ")")) + 
                xlab("weighted mean difference") + 
                scale_alpha_manual(paste0(pValueStrLabel, " < ", significanceThresholdCur), values = c(alphaValueNonSign, 1), labels = c("no", "yes")) + 
                geom_hline(yintercept = transform_yValues(significanceThresholdCur, addPseudoCount = TRUE, nPermutations = par.l$nPermutations), linetype = "dotted") 
            
            if (nrow(ggrepel_df) <= par.l$maxTFsToLabel) {
                
                if (par.l$plotRNASeqClassification) {
                    # TODO: Originally "classification" 
                    g = g +  geom_label_repel(data = ggrepel_df, aes_string("weighted_meanDifference", pValueScoreCur, label = "TF", fill = colnameClassificationFinal),
                                              size = TFLabelSize, fontface = 'bold', color = 'white',
                                              segment.size = 0.3, box.padding = unit(0.2, "lines"), max.iter = 5000,
                                              label.padding = unit(0.2, "lines"), # how thick is connectin line
                                              nudge_y = 0.05, nudge_x = 0,  # how far from center points
                                              segment.alpha = .8, segment.color = par.l$plot_grayColor, show.legend = FALSE)
                } else {
                    g = g +  geom_label_repel(data = ggrepel_df, aes_string("weighted_meanDifference", pValueScoreCur, label = "TF"),
                                              size = TFLabelSize, fontface = 'bold', color = 'black',
                                              segment.size = 0.3, box.padding = unit(0.2, "lines"), max.iter = 5000,
                                              label.padding = unit(0.2, "lines"), # how thick is connectin line
                                              nudge_y = 0.05, nudge_x = 0,  # how far from center points
                                              segment.alpha = .8, segment.color = par.l$plot_grayColor, show.legend = FALSE)
                }
            } else {
                
                flog.warn(paste0(" Not labeling significant TFs, maximum of ", par.l$maxTFsToLabel, " exceeded for ", pValThrStr, " and ", pValueStrCur))
                
                labelPlot = paste0("*TF labeling skipped because number of significant TFs\nexceeds the maximum of ", par.l$maxTFsToLabel, " (", nrow(ggrepel_df), ")")
    
                g = g + annotate("text", label = labelPlot, x = 0, y = ymax, size = 3)
                
            } # end else
    
              g = g + theme_bw() + 
                  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)), 
                        axis.title.x = element_text(size=rel(1.5)), axis.title.y = element_text(size=rel(1.5)),
                        legend.title=element_text(size=rel(1.5)), legend.text=element_text(size=rel(1.5))) 
              
              if (par.l$plotRNASeqClassification) {
                g = g + guides(alpha = guide_legend(override.aes = list(size=5), order = 2),
                               fill = guide_legend(override.aes = list(size=5), order = 3),
                               color = guide_legend(override.aes = list(size=5), order = 1))
                
                allPlots.l[[as.character(thresholdCur)]] [[pValThrStr]] [[paste0(showClasses,collapse = "-")]] [[pValueStrCur]] = g
              
                } else {
                g = g + guides(alpha = guide_legend(override.aes = list(size=5), order = 2),
                               fill = guide_legend(override.aes = list(size=5), order = 3))
         
                allPlots.l[[as.character(thresholdCur)]] [[pValThrStr]] [[pValueStrCur]] = g
              }
              
          } # end separately for raw and adjusted p-values
            
        } # end for all showClasses
    
    } # end for different significance thresholds
    
    
    ####################
    # VOLCANO PLOT PDF #
    ####################
    height = width = max(nTF_label / 15 , par.l$volcanoPlot_minDimensions)
    
    if (par.l$plotRNASeqClassification) {
        
        pdf(file = par.l$file_plotVolcano, height = height, width = width, useDingbats = FALSE)
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        message = paste0("This file is intentionally empty.\n",
                         "\nSee the other similarly named files with additional \".q\" for results.\n",
                        "\nIncreasing q-values indicate higher stringency\n and therefore more TFs are classified as undetermined.\n")
        text(x = 0.5, y = 0.5, message, cex = 0.9, col = "black")
        dev.off()
        

        fileCur = gsub(".pdf", paste0(".q", thresholdCur, ".pdf"), par.l$file_plotVolcano)
        pdf(file = fileCur, height = height, width = width, useDingbats = FALSE)
        
    } else {
        
        # Only one variant here
        pdf(file = par.l$file_plotVolcano, height = height, width = width, useDingbats = FALSE)

    }
    
    
    
    for (pValueStrCur in c("pvalueAdj", "pvalue")) {
        
        for (significanceThresholdCur in par.l$significanceThresholds) {
          
          for (showClasses in classesList.l) {
    
                  if (par.l$plotRNASeqClassification) {
                      plot(allPlots.l[[as.character(thresholdCur)]] [[as.character(significanceThresholdCur)]] [[paste0(showClasses,collapse = "-")]] [[pValueStrCur]])
                  } else {
                      plot(allPlots.l[[as.character(thresholdCur)]] [[as.character(significanceThresholdCur)]] [[pValueStrCur]])
                  }
          } # end for each class
        } # end for each threshold
    } # end for both raw and adjusted p-values
    dev.off()

    
    
} # end for all stringency thresholds for the classification

#########################
# FINALLY, SAVE TO DISK #
#########################
output.global.TFs.origReal = dplyr::select(output.global.TFs.origReal, -one_of("permutation", "yValue"))
output.global.TFs.origReal.transf = dplyr::mutate_if(output.global.TFs.origReal, is.numeric, as.character)
write_tsv(output.global.TFs.origReal.transf, path = par.l$file_output_summary, col_names = TRUE)
saveRDS(allPlots.l, file = par.l$file_output_plots)

.printExecutionTime(start.time)
flog.info("Session info: ", sessionInfo(), capture = TRUE)
