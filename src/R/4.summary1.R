start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################
library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "modeest", "checkmate", "ggrepel"), verbose = FALSE)

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/4.summary1.R.rds")
createDebugFile(snakemake)

###################
#### PARAMETERS ###
###################
par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$min_pValue   = .Machine$double.xmin
par.l$FDR_threshold = 0.05
par.l$plot_min_diffMean = 0.02
par.l$log_minlevel = "INFO"


#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(names(snakemake@input), c("", "peaks", "TF", "analysesTFOutput"))

par.l$file_input_peaks = snakemake@input$peaks
assertFileExists(par.l$file_input_peaks, access = "r")

par.l$files_input_TF_summary = snakemake@input$TF
for (fileCur in par.l$files_input_TF_summary) {
  assertFileExists(fileCur, access = "r")
}

## OUTPUT ##
assertList(snakemake@output, min.len = 1)
#assertSubset(names(snakemake@output), c("", "volcanoPlot", "outputTable"))
assertSubset(names(snakemake@output), c("", "outputTable"))

#par.l$file_output_volcanoPlot = snakemake@output$volcanoPlot
par.l$file_output_table       = snakemake@output$outputTable

## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]



allDirs = c(#dirname(par.l$file_output_volcanoPlot), 
            dirname(par.l$file_output_table),
            dirname(par.l$file_log)
          )

testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
######################

startLogger(par.l$file_log, par.l$log_minlevel, removeOldLog = TRUE)
printParametersLog(par.l)


################
# Collect data #
################

par.l$filesTransl.l = list()
nTF = length(par.l$files_input_TF_summary)

# Create translation table
for (fileCur in par.l$files_input_TF_summary) {
  elements = strsplit(fileCur, split = "/",  fixed = TRUE)[[1]]
  hit =  which(grepl(pattern = "^extension", elements))
  stopifnot(length(hit) == 1 & hit != 1)
  par.l$filesTransl.l[[fileCur]] = elements[hit - 1]  
}



peaks.df = read_tsv(par.l$file_input_peaks, col_types = cols())

assertSubset(colnames(peaks.df), c("permutation", "position", "D2_baseMean", "D2_l2FC", "D2_ldcSE", "D2_stat", "D2_pval", "D2_padj")) #, "vst_diff"))

summary.df = NULL

nTFs = length(par.l$filesTransl.l)

flog.info(paste0("Using ", nTFs, " TFs"))

for (fileCur in par.l$files_input_TF_summary) {
  
  TF = par.l$filesTransl.l[[fileCur]]
  stopifnot(!is.null(TF))
  if (file.exists(fileCur)) {
    
    stats.df = readRDS(fileCur)
    
    if (is.null(summary.df)) {
      summary.df = stats.df
    } else {
      summary.df = rbind(summary.df, stats.df)
    }
    

  } else {
    message = paste0("File missing: ", fileCur)
    flog.warn(message)
    warning(message)
  }
  
}

# TODO: WHy same values for both TF

nRowsSummary = nrow(summary.df)

flog.info(paste0(" Imported ", nRowsSummary, " TFs out of a list of ", nTFs, ". Missing: ", nTFs - nRowsSummary))

nTFMissing = length(which(is.na(summary.df$Pos_l2FC)))
if (nTFMissing == nrow(summary.df)) {
  error = "All TF have missing data. Cannot continue. Add more samples or change the peaks."
  flog.fatal(error)
  stop(error)
}
  

# Replace p-values of 0 with the smallest p-value on the system
summary.df$Ttest_pval[summary.df$Ttest_pval == 0] = .Machine$double.xmin

# TODO: was previously assigned to modeNum, but why?
mode_peaks = mlv(round(peaks.df$D2_l2FC, 2), method = "mfv", na.rm = TRUE)

summary.df = summary.df %>%
              dplyr::mutate(
                  adj_pvalue = p.adjust(Ttest_pval, method = "fdr"),
                  Diff_mean  = Mean_l2FC   - mean  (peaks.df$D2_l2FC, na.rm = TRUE), 
                  DiffMedian = Median_l2FC - median(peaks.df$D2_l2FC, na.rm = TRUE),
                  Diff_mode  = Mode_l2FC - mode_peaks[[1]],    
                  Diff_skew  = Modeskewness - mode_peaks[[2]])  %>%
              na.omit(summary.df)



write_tsv(summary.df, par.l$file_output_table) 

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)