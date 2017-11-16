start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################

assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "vsn", "modeest", "checkmate", "limma", "geneplotter", "RColorBrewer", "tools"), verbose = TRUE)

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# snakemake = readRDS("/scratch/carnold/CLL/TF_actNew/output/Logs_and_Benchmarks/3.analyzeTF.R_TF=VDR.B.rds")
createDebugFile(snakemake, "3.analyzeTF.R")

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
assertSubset(names(snakemake@input), c("", "overlapFiles", "sampleDataR", "peakFile", "peakFile2", "normFacs", "plotsPerm"))

par.l$file_input_peakTFOverlaps  = snakemake@input$overlapFiles
for (fileCur in par.l$file_input_peakTFOverlaps) {
  assertFileExists(fileCur, access = "r")
}

par.l$file_input_metadata = snakemake@input$sampleDataR
assertFileExists(par.l$file_input_metadata, access = "r")

par.l$file_input_peaks = snakemake@input$peakFile
assertFileExists(par.l$file_input_peaks, access = "r")

par.l$file_input_peak2 = snakemake@input$peakFile2
assertFileExists(par.l$file_input_peak2, access = "r")

par.l$file_input_normFacs = snakemake@input$normFacs
assertFileExists(par.l$file_input_normFacs, access = "r")

# TODO: par.l$file_input_metadata


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

par.l$nPermutations = snakemake@config$par_general$nPermutations
assertIntegerish(par.l$nPermutations, lower = 0)

## PARAMS ##
assertList(snakemake@params, min.len = 1)
assertSubset(names(snakemake@params), c("", "doCyclicLoess"))

par.l$doCyclicLoess = as.logical(snakemake@params$doCyclicLoess)
assertFlag(par.l$doCyclicLoess)

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
startLogger(par.l$file_log, par.l$log_minlevel, appenderName = "file", removeOldLog = TRUE)
printParametersLog(par.l)


#################
# READ METADATA #
#################

sampleData.l = readRDS(par.l$file_input_metadata)

# Initiate data structures that are populated hereafter
res_DESeq.l = list()

TF_output.df = tribble(~permutation, ~TF, ~chr, ~MSS, ~MES, ~strand, ~PSS, ~PES, ~annotation, ~ID, ~identifier, ~baseMean, ~log2FoldChange, ~lfcSE, ~stat, ~pvalue, ~padj)

outputSummary.df = tribble(~permutation, ~TF, ~Pos_l2FC, ~Mean_l2FC, ~Median_l2FC, ~Mode_l2FC, ~sd, ~Ttest_pval, ~Modeskewness, ~T_statistic, ~TFBS_num)

coverageAll.df = NULL

if (length(sampleData.l) == 0) {
  message = "Length of sampleData.l list is 0 but is has to be at least 1. Rerun the prepareData rule."
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}

sampleData.df = sampleData.l[["permutation0"]]

flog.info(paste0("Iterating over ", length(par.l$file_input_peakTFOverlaps), " TF-peak overlap files"))

errorOccured = FALSE
errorMessage = ""

for (fileCur in par.l$file_input_peakTFOverlaps) {
  
  flog.info(paste0("Read file ", fileCur))
  info = file.info(fileCur)
  if (info$size > 0) {
    
    peaksCur.df = read_tsv(fileCur, col_names = c("chr","MSS","MES","annotation","ID","strand","coverage"), 
                           col_types = cols(chr = col_character(), 
                                            MSS = col_integer(),
                                            MES = col_integer(),
                                            annotation = col_character(),
                                            ID = col_integer(),
                                            #ID = col_character(),
                                            strand = col_character(),
                                            coverage = col_integer()))
    
    if (nrow(problems(peaksCur.df)) > 0) {
      flog.fatal(paste0("Parsing errors: "), problems(peaksCur.df), capture = TRUE)
      stop("Error when parsing the file ", fileCur, ", see warnings")
    }
    
    
    assertDataFrame(peaksCur.df, ncols = 7)
    #assertDataFrame(peaksCur.df, ncols = 6)
    
    #colnames(peaksCur.df) = c("chr","MSS","MES","annotation","ID","coverage")
    
    peaksCur.df = mutate(peaksCur.df, identifier = paste0(chr,":", MSS, "-",MES))
    
    # get only unique identifiers
    peaksCur.df =  distinct(peaksCur.df, identifier, .keep_all = TRUE)
    
    if (!is.null(coverageAll.df) && nrow(peaksCur.df) != nrow(coverageAll.df)){
      message = paste0("Wrong number of rows in file ", fileCur, ". Expected: ", nrow(coverageAll.df), ", found: ", nrow(peaksCur.df), ". Some overlap files might be corrupt, delete the following files and reproduce them: ", paste0(par.l$file_input_peakTFOverlaps, collapse = ","))
      flog.fatal(message)
      stop(message)
    }
    
    # Only use the "coverage column and append to the final list
    coverageAll.df = dplyr::bind_cols(coverageAll.df, peaksCur.df[,"coverage"])
    
  } else {
      message = paste0("Empty file ", fileCur, ", skip...")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  }
  
  
}

if (!is.null(coverageAll.df)) {
  
  # Add row mean
  coverageAll.df$mean = apply(coverageAll.df, 1, mean)
  
  # Add metadata (all except coverage)
  columns = c("chr","MSS","MES","annotation","ID","strand","identifier")
  #columns = c("chr","MSS","MES","annotation","ID","identifier")
  assertSubset(columns, colnames(peaksCur.df))
  coverageAll.df = dplyr::bind_cols(coverageAll.df, peaksCur.df[, columns])
  colnames(coverageAll.df) = c(sampleData.df$name, "mean", columns)
  
  # Group by ID
  # take only the maximum row mean of all samples, sample with biggest coverage
  coverageAll_grouped.df = coverageAll.df %>%
    dplyr::group_by(ID) %>%
    dplyr::slice(which.max(mean))
  
  TF.table.m = as.matrix(coverageAll_grouped.df[,sampleData.df$name])
  colnames(TF.table.m) = sampleData.df$name
  rownames(TF.table.m) = coverageAll_grouped.df$identifier
  
  
  
  # Create formula based on user-defined design
  designFormula = convertToFormula(par.l$designFormula, colnames(sampleData.df))
  
  
} else {
  errorOccured <<-TRUE
  errorMessage <<- "Could not find any overlaps with peaks. This TF will be ignored in subsequent steps"
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


for (permutationCur in 0:par.l$nPermutations) {
  
  flog.info(paste0("Running permutation ", permutationCur))
  permutationName = paste0("permutation", permutationCur)
  
  assertSubset(permutationName, names(sampleData.l))
  sampleData.df = sampleData.l[[permutationName]]

  
  
  TF.cds.filt = tryCatch( {
    
    # create Deseq object from the TF specific data
    TF.cds <- DESeqDataSetFromMatrix(countData = TF.table.m,
                                     colData = sampleData.df,
                                     design = designFormula)
    
    # We here again take the gene-specific normalization factors from the previously calculated ones based on the peaks,
    # but subset this to only those corresponding to the TF of interest

    normFacs = normFacs.l[[permutationName]]
    
    # assertSubset(rownames(normFacs))
    if (par.l$doCyclicLoess) {
      # TODO: make this less fault tolerant...
      normalizationFactors(TF.cds) <- normFacs[coverageAll_grouped.df$ID,, drop = FALSE]
    } else {
      sizeFactors(TF.cds) = normFacs
    }
    
    # low RC, check by rowMean
    TF.cds[rowMeans(counts(TF.cds)) > 0, ]
    
  }, error = function(e) {
    errorMessage <<- "Could not initiate DESeq."
    checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    errorOccured <<- TRUE
    
  }
  )
  
  if (!errorOccured) {
    
    # Deseq 2 functions
    # off/on shrinkage of the log2FC - betaPrior
    # with the simulations of negative binomial distribution increase the sample size
    
    # Run the local fit first, if that throws an error try the default fit type
    
    res_DESeq = tryCatch( {
      DESeq(TF.cds.filt,fitType = 'local')
      
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
    
    
    
    
    assertSubset(rownames(res_DESeq.df), peaksCur.df$identifier)
    
    
    rm_col = c("coverage.x","chr.y","annotation.y","coverage.y","identifier.y" )
    order = c("permutation", "TF", "chr","MSS","MES","strand", "PSS","PES","annotation","ID", "identifier","baseMean", "log2FoldChange","lfcSE","stat", "pvalue","padj")
    
    # dplyr::mutate(TF = par.l$TF) gives the following weird error message: Error: Unsupported type NILSXP for column "TF"
    TFCur = par.l$TF
    
    
    TF_outputCur.df = res_DESeq.df %>%
      rownames_to_column(var = "identifier") %>%
      dplyr::full_join(peaksCur.df,by = c("identifier")) %>%
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
      filenameCurDiag    = paste0(file_path_sans_ext(par.l$file_output_plot_diagnostic), "_permutation", permutationCur, ".pdf")
      
      computeDESeqDiagnosticPlots(res_DESeq, filenameCurDiag, maxPairwiseComparisons = par.l$maxPairwiseComparisonsDiagnosticPermutations)
      
    } else {
      computeDESeqDiagnosticPlots(res_DESeq, filenameCurDiag)
    }
    
    
    pdf(filenameCurSummary)
    
    comparisonDESeq = getComparisonFromDeSeqObject(res_DESeq, par.l$designFormula)
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
