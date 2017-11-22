start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################

library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "vsn", "modeest", "checkmate", "limma", "geneplotter", "RColorBrewer", "tools"), verbose = FALSE)

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} and {TF} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/3.analyzeTF.{TF}.R.rds")
# snakemake = readRDS("/scratch/carnold/CLL/27ac_TF/output/Logs_and_Benchmarks/3.analyzeTF.R_TF=MAFK.S.rds")
createDebugFile(snakemake)

###################
#### PARAMETERS ###
###################

par.l = list()

par.l$verbose = TRUE
par.l$log_minlevel = "INFO"
par.l$maxPairwiseComparisonsDiagnosticPermutations = 2

#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(names(snakemake@input), c("", "overlapFile", "sampleDataR", "peakFile", "peakFile2", "normFacs", "plotsPerm"))

par.l$file_input_peakTFOverlaps  = snakemake@input$overlapFile
assertFileExists(par.l$file_input_peakTFOverlaps, access = "r")

par.l$file_input_metadata = snakemake@input$sampleDataR
assertFileExists(par.l$file_input_metadata, access = "r")

par.l$file_input_peaks = snakemake@input$peakFile
assertFileExists(par.l$file_input_peaks, access = "r")

par.l$file_input_peak2 = snakemake@input$peakFile2
assertFileExists(par.l$file_input_peak2, access = "r")

par.l$file_input_normFacs = snakemake@input$normFacs
assertFileExists(par.l$file_input_normFacs, access = "r")


## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "outputTSV", "outputRDS", "plot_diagnostic", "plot_diagnosticPerm", "plot_TFSummary", "plot_TFSummaryPerm", "DESeqObj"))

par.l$file_output_summaryAll      = snakemake@output$outputTSV
par.l$file_output_summaryStats    = snakemake@output$outputRDS
par.l$file_output_plot_diagnostic = snakemake@output$plot_diagnostic
par.l$file_output_plot_TFSummary  = snakemake@output$plot_TFSummary
par.l$file_output_DESeq           = snakemake@output$DESeqObj

## WILDCARDS ##
assertList(snakemake@wildcards, min.len = 1)
assertSubset(names(snakemake@wildcards), c("", "TF"))

par.l$TF = snakemake@wildcards$TF
assertCharacter(par.l$TF, len = 1, min.chars = 1)

## CONFIG ##
assertList(snakemake@config, min.len = 1)

par.l$designFormula = snakemake@config$par_general$designContrast
assertCharacter(par.l$designFormula, len = 1, min.chars = 3)

par.l$designFormula = snakemake@config$par_general$designContrast
checkAndLogWarningsAndErrors(par.l$designFormula, checkCharacter(par.l$designFormula, len = 1, min.chars = 3))

par.l$designFormulaVariableTypes = snakemake@config$par_general$designVariableTypes
checkAndLogWarningsAndErrors(par.l$designFormulaVariableTypes, checkCharacter(par.l$designFormulaVariableTypes, len = 1, min.chars = 3))

par.l$nPermutations = snakemake@config$par_general$nPermutations
assertIntegerish(par.l$nPermutations, lower = 0)

## PARAMS ##
assertList(snakemake@params, min.len = 1)
assertSubset(names(snakemake@params), c("", "doCyclicLoess", "allBAMS"))

par.l$doCyclicLoess = as.logical(snakemake@params$doCyclicLoess)
assertFlag(par.l$doCyclicLoess)

par.l$allBAMS = snakemake@params$allBAMS
for (fileCur in par.l$allBAMS) {
  checkAndLogWarningsAndErrors(fileCur, checkFileExists(fileCur))
}

## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]


allDirs = c(dirname(par.l$file_output_summaryAll), 
            dirname(par.l$file_output_summaryStats), 
            dirname(par.l$file_output_plot_diagnostic), 
            dirname(par.l$file_output_plot_TFSummary),
            dirname(par.l$file_output_DESeq),
            dirname(par.l$file_log)
)


testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel, removeOldLog = TRUE)
printParametersLog(par.l)


#################
# READ METADATA #
#################

sampleData.l = readRDS(par.l$file_input_metadata)

# Initiate data structures that are populated hereafter
res_DESeq.l = list()

TF_output.df = tribble(~permutation, ~TF, ~chr, ~MSS, ~MES, ~strand, ~PSS, ~PES, ~annotation, ~ID, ~identifier, ~baseMean, ~log2FoldChange, ~lfcSE, ~stat, ~pvalue, ~padj)

outputSummary.df = tribble(~permutation, ~TF, ~Pos_l2FC, ~Mean_l2FC, ~Median_l2FC, ~Mode_l2FC, ~sd, ~Ttest_pval, ~Modeskewness, ~T_statistic, ~TFBS_num)

# TODO
coverageAll.df = NULL

if (length(sampleData.l) == 0) {
  message = "Length of sampleData.l list is 0 but is has to be at least 1. Rerun the prepareData rule."
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}

sampleData.df = sampleData.l[["permutation0"]]

#####################
# READ OVERLAP FILE #
#####################

colnamesNew = sampleData.df$SampleID[which(sampleData.df$bamReads %in% par.l$allBAMS)]
if (length(colnamesNew) != nrow(sampleData.df)) {
    message = "Could not grep sampleIDs from filenames."
    checkAndLogWarningsAndErrors(NULL,  message, isWarning = FALSE)
}
overlapsAll.df = read_tsv(par.l$file_input_peakTFOverlaps, col_names = FALSE, col_types = cols())

if (nrow(problems(overlapsAll.df)) > 0) {
  flog.fatal(paste0("Parsing errors: "), problems(overlapsAll.df), capture = TRUE)
  stop("Error when parsing the file ", fileCur, ", see warnings")
}

colnames(overlapsAll.df) = c("chr","MSS","MES","annotation","ID","strand","fileOrigin", colnamesNew)


overlapsAll.df = overlapsAll.df %>%
                  dplyr::mutate(identifier = paste0(chr,":", MSS, "-",MES)) %>%
                  dplyr::mutate(mean = apply(dplyr::select(overlapsAll.df, one_of(colnamesNew)), 1, mean))  %>%
                  dplyr::distinct(identifier, .keep_all = TRUE) %>%
                  dplyr::select(-one_of("fileOrigin"))


if (nrow(overlapsAll.df) > 0) {
  
  # Group by ID
  # take only the maximum row mean of all samples, sample with biggest coverage
  coverageAll_grouped.df = overlapsAll.df %>%
    dplyr::group_by(ID) %>%
    dplyr::slice(which.max(mean))
  
  TF.table.m = as.matrix(coverageAll_grouped.df[,sampleData.df$SampleID])
  colnames(TF.table.m) = sampleData.df$SampleID
  rownames(TF.table.m) = coverageAll_grouped.df$identifier
  
  
  
  # Create formula based on user-defined design
  designFormula = convertToFormula(par.l$designFormula, colnames(sampleData.df))
  
  
} else {
  message <<- "Could not find any overlaps with peaks. This TF will be ignored in subsequent steps"
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
}


################################
# ITERATE THROUGH PERMUTATIONS #
################################

normFacs.l = readRDS(par.l$file_input_normFacs)
peaksFiltered.df = readRDS(par.l$file_input_peaks)

# Adjust the number of permutations in case less have been computed
 if (par.l$nPermutations + 1 < length(sampleData.l)) {
   message = paste0("In the output objects, more permutations seem to be stored. They will be ignored and the uoriginal value of nPermutations will be used")
   checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
 } else if (par.l$nPermutations + 1 > length(sampleData.l)) {
   valueNew = length(sampleData.l) - 1
   message = paste0("The value of the parameter nPermutations differs from what is saved in the output objects. The value of nPermutations will be adjusted to ", valueNew)
   checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
   par.l$nPermutations = valueNew
 }



errorOccured = FALSE
errorMessage = ""



for (permutationCur in 0:par.l$nPermutations) {
  
  flog.info(paste0("Running permutation ", permutationCur))
  permutationName = paste0("permutation", permutationCur)
  
  assertSubset(permutationName, names(sampleData.l))
  sampleData.df = sampleData.l[[permutationName]]

  
  
  TF.cds = tryCatch( {
    
    # create Deseq object from the TF specific data
    TF.cds <- DESeqDataSetFromMatrix(countData = TF.table.m,
                                     colData = sampleData.df,
                                     design = designFormula)
    
   
    
  }, error = function(e) {
    errorMessage <<- "Could not initiate DESeq."
    checkAndLogWarningsAndErrors(NULL,  errorMessage, isWarning = TRUE)
    errorOccured <<- TRUE
    
  }
  )
  
  # We here again take the gene-specific normalization factors from the previously calculated ones based on the peaks,
  # but subset this to only those corresponding to the TF of interest
  
  normFacs = normFacs.l[[permutationName]]
  
  # Sanity check
  if (length(which(!coverageAll_grouped.df$annotation %in% rownames(normFacs))) > 0) {
    errorMessage <<- "Inconsistency detected between the normalization factor rownames and the object coverageAll_grouped.df"
    checkAndLogWarningsAndErrors(NULL,  errorMessage, isWarning = TRUE)
  }
  
  if (!identical(colnames(TF.cds), colnames(normFacs))) {
    errorMessage <<- "Column names differenht between TF.cds and normFacs"
    checkAndLogWarningsAndErrors(NULL,  errorMessage, isWarning = TRUE)
  }
  
  # assertSubset(rownames(normFacs))
  if (par.l$doCyclicLoess) {

    
    rownamesTFs = which(rownames(normFacs) %in% coverageAll_grouped.df$annotation)
    normalizationFactors(TF.cds) <- normFacs[rownamesTFs,, drop = FALSE]
  } else {
    sizeFactors(TF.cds) = normFacs
  }
  
  # low RC, check by rowMean
  TF.cds.filt = TF.cds[rowMeans(counts(TF.cds)) > 0, ]
  
  if (!errorOccured) {
    
    # Deseq 2 functions
    # off/on shrinkage of the log2FC - betaPrior
    # with the simulations of negative binomial distribution increase the sample size
    
    # Run the local fit first, if that throws an error try the default fit type
    
    res_DESeq = tryCatch( {
        suppressMessages(DESeq(TF.cds.filt,fitType = 'local'))
      
    }, error = function(e) {
      message = "Could not run DESeq with local fitting, retry with default fitting type..."
      checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      
      res_DESeq = tryCatch( {
        DESeq(TF.cds.filt)
        
      }, error = function(e) {
        errorMessage <<- "Could not run DESeq with regular fitting either, set all values to NA."
        checkAndLogWarningsAndErrors(NULL, errorMessage, isWarning = TRUE)
        errorOccured <<- TRUE
      }
      )
      
    }
    )
    
  }
  
  if (errorOccured) {
    
    # TODO: check with new NA row if this still works downstream
    
    TF_output.df = add_row(TF_output.df, 
                             permutation = permutationCur, 
                             TF          = par.l$TF
                            # Rest is set to NA
                           )
    
    outputSummary.df = add_row(outputSummary.df, 
                               permutation = permutationCur,
                               TF  = par.l$TF
                               # Rest is set to NA
                              )
    
    
    # Plot dummy graphs that clearly state that something went wrong
    filenameCurSummary = par.l$file_output_plot_TFSummary
    filenameCurDiag    = par.l$file_output_plot_diagnostic
    
    if (permutationCur > 0) {
      
      filenameCurSummary = paste0(file_path_sans_ext(par.l$file_output_plot_TFSummary) , "_permutation", permutationCur, ".pdf")
      filenameCurDiag    = paste0(file_path_sans_ext(par.l$file_output_plot_diagnostic), "_permutation", permutationCur, ".pdf")
    } 
    
    par(mar = c(0,0,0,0))
    
    for (filenameCur in c(filenameCurDiag, filenameCurSummary)) {
      pdf(filenameCur)
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      message = paste0("Insufficient data to run analysis.\n", errorMessage, "\nThis TF will be ignored in subsequent steps.")
      text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
      dev.off()
    }
    
    # Save an empty list
    res_DESeq = list()
  
    # Restore mar to the default
    par(mar = c(5, 4, 4, 2) + 0.1)
    
    
  } else {
    
  
    res_DESeq.df <- as.data.frame(DESeq2::results(res_DESeq))
    
    
    # addition  02.06
    final.TF.df = data_frame("position"    = rownames(res_DESeq.df), 
                             "D2_baseMean" = res_DESeq.df$baseMean,
                             "D2_l2FC"     = res_DESeq.df$log2FoldChange,
                             "D2_ldcSE"    = res_DESeq.df$lfcSE,
                             "D2_stat"     = res_DESeq.df$stat,
                             "D2_pval"     = res_DESeq.df$pvalue, 
                             "D2_padj"     = res_DESeq.df$padj
    )
    
    # assign final.peaks.df to the peaks.df and filter away NAs at the p.adjust
    
    
    
    
    assertSubset(rownames(res_DESeq.df), overlapsAll.df$identifier)
    
    # todo: check the coverage columns
    rm_col = c("chr.y","annotation.y","identifier.y" )
    order = c("permutation", "TF", "chr","MSS","MES","strand", "PSS","PES","annotation","ID", "identifier","baseMean", "log2FoldChange","lfcSE","stat", "pvalue","padj")
    
    # dplyr::mutate(TF = par.l$TF) gives the following weird error message: Error: Unsupported type NILSXP for column "TF"
    TFCur = par.l$TF
    
    
    TF_outputCur.df = res_DESeq.df %>%
      rownames_to_column(var = "identifier") %>%
      dplyr::full_join(overlapsAll.df,by = c("identifier")) %>%
      dplyr::full_join(peaksFiltered.df, by = "ID") %>%
      dplyr::filter(!is.na(baseMean)) %>%
      dplyr::rename(annotation = annotation.x, chr = chr.x, identifier = identifier.x) %>%
      dplyr::select(-one_of(rm_col)) %>%
      dplyr::mutate(TF = TFCur, permutation = permutationCur) %>%
      dplyr::select(one_of(order)) %>%
      dplyr::arrange(chr)
    
    TF_output.df = rbind(TF_output.df, TF_outputCur.df)
    
    # d) Comparisons between peaks and binding sites
    
    
    # TODO: Not needed
    #peaks_C = nrow(peaks.df[peaks.df$log2FoldChange > 0,])/nrow(peaks.df)
    
    peaks.df = read_tsv(par.l$file_input_peak2, col_types = cols())
    
    modeNum     = mlv(round(final.TF.df$D2_l2FC, 2), method = "mfv", na.rm = TRUE)
    Ttest       = t.test(final.TF.df$D2_l2FC, peaks.df$D2_l2FC)
    
    outputSummary.df = add_row(outputSummary.df,
                          permutation     = permutationCur,
                          TF              = par.l$TF,
                          Pos_l2FC        = nrow(final.TF.df[final.TF.df$D2_l2FC > 0,]) / nrow(final.TF.df),
                          Mean_l2FC       = mean(final.TF.df$D2_l2FC, na.rm = TRUE),
                          Median_l2FC     = median(final.TF.df$D2_l2FC, na.rm = TRUE),
                          Mode_l2FC       = modeNum[[1]],
                          sd              = sd(final.TF.df$D2_l2FC, na.rm = TRUE),
                          Ttest_pval      = Ttest$p.value,
                          Modeskewness    = modeNum[[2]], 
                          T_statistic     = Ttest$statistic[[1]], 
                          TFBS_num        = nrow(final.TF.df)
                        )
  
    
    ############
    ############
    ## GRAPHS ##
    ############
    ############
    
    filenameCurSummary = par.l$file_output_plot_TFSummary
    filenameCurDiag    = par.l$file_output_plot_diagnostic
    
    if (permutationCur > 0) {
      
      filenameCurSummary = paste0(file_path_sans_ext(par.l$file_output_plot_TFSummary) , "_permutation", permutationCur, ".pdf")

      # Disabled for now
      # filenameCurDiag    = paste0(file_path_sans_ext(par.l$file_output_plot_diagnostic), "_permutation", permutationCur, ".pdf")
      # computeDESeqDiagnosticPlots(res_DESeq, filenameCurDiag, maxPairwiseComparisons = par.l$maxPairwiseComparisonsDiagnosticPermutations)
      
    } else {
      computeDESeqDiagnosticPlots(res_DESeq, filenameCurDiag)
    }
    
    # Determine the type of the variable that is part of the DESeq comparison:
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
    
    
    pdf(filenameCurSummary)
    
    comparisonDESeq = getComparisonFromDeSeqObject(res_DESeq, par.l$designFormula, datatypeVariableToPermute = datatypeVariableToPermute)
    xlabLabel = paste0(" log2 FC ", comparisonDESeq)
    
    # Density plot
    TF_dens = ggplot() + geom_density(aes(x = peaks.df$D2_l2FC,fill = "A" ),
                                      alpha = .5, color = "black") +
      geom_density(aes(x = final.TF.df$D2_l2FC,fill = "B"),size = 1, alpha = .7) +
      xlab(xlabLabel) +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", colour = "black", size = 10, margin = margin(25,0,0,0)),
            axis.title.y = element_text(face = "bold", colour = "black", size = 10, margin = margin(0,25,0,0)),
            axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = c(0.9,0.9),
            legend.justification = "center",
            legend.title = element_blank()) +
      scale_fill_manual(values = c("A" = "grey50" , "B" = "blue"), labels = c("PEAKS", par.l$TF))
    
    plot(TF_dens)
    #ggsave(plot = TF_dens, filename = par.l$file_output_plot_density, width = 4, height = 4, useDingbats = FALSE, dpi = 600)
    
    ########
    # ECDF #
    ########
    
    # create ecdf plots for each TF
    ECDF_TF = ggplot() + 
      stat_ecdf(aes(x = final.TF.df$D2_l2FC,colour = paste0("", par.l$TF))) +
      stat_ecdf(aes(x = peaks.df$D2_l2FC,colour = "Peaks" )) +
      xlab(xlabLabel) + 
      guides(colour = guide_legend(title = "ORIGIN"))
    
    plot(ECDF_TF)
    #ggsave(plot = ECDF_TF, filename = par.l$file_output_plot_ecdf, width = 6, height = 4, useDingbats = FALSE, dpi = 600)
    
    dev.off()
    
    
    
  } # end if (!errorOccured)
  
  res_DESeq.l[[permutationName]] = res_DESeq
  
} # end for each permutation

write_tsv(TF_output.df, path = par.l$file_output_summaryAll)
saveRDS(outputSummary.df, file = par.l$file_output_summaryStats)
saveRDS(res_DESeq.l, file = par.l$file_output_DESeq)

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)
