start.time  <-  Sys.time()

#########################
# LIBRARY AND FUNCTIONS #
#########################

source("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/functions.R")

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "modeest", "checkmate", "ggrepel"), verbose = TRUE)


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


############################################
# READ AND VALIDATE COMMAND LINE ARGUMENTS #
############################################
args <- commandArgs(trailingOnly = TRUE)


if (length(args) != 6) {
  stop("Expecting 6 arguments but found ", length(args),". Exiting.")
} else {
  par.l$file_input_peaks       = args[1]
  par.l$files_input_TF_summary = args[2]
  par.l$file_output_volcanoPlot= args[3]
  par.l$file_output_table      = args[4]
  par.l$TFs                    = args[5]
  par.l$file_log               = args[6]

}

# args = c(
#   
#   "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/PEAKS/MC.peaks.tsv"                 ,
#   "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/TGIF1.D/extension50/TGIF1.D.MC.summary.rds,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.summary.rds"                 ,
#   "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/FINAL_OUTPUT/extension50/MC.TF_vs_peak_distribution.pdf"   ,              
#   "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/FINAL_OUTPUT/extension50/MC.TF_vs_peak_distribution.tsv"  ,             
#   "TGIF1.D,ZN350"             ,    
#   "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/Logs_and_Benchmarks/createVolcanoPlot.R.log"
#   
# )

#####################
# VERIFY PARAMETERS #
#####################

assertFileExists(par.l$file_input_peaks, access = "r")
assertCharacter(par.l$files_input_TF_summary, len = 1)
assertDirectoryExists(dirname(par.l$file_log), access = "w")

allDirs = c(dirname(par.l$file_output_volcanoPlot), 
            dirname(par.l$file_output_table)
          )


for (dirname in unique(allDirs)) {
  
  if (!testDirectoryExists(dirname)) {
    dir.create(dirname, recursive = TRUE)
  } else {
    assertDirectoryExists(dirname, access = "w")
  }
}

for (dirname in unique(allDirs)) {
  
  if (!testDirectoryExists(dirname)) {
    dir.create(dirname, recursive = TRUE)
  } else {
    assertDirectoryExists(dirname, access = "w")
  }
}

fileList = strsplit(par.l$files_input_TF_summary, split = ",", fixed = TRUE)[[1]]
fileListTransl.l = list()
nTF = length(fileList)

for (fileCur in fileList) {
  assertFileExists(fileCur, access = "r")
  TFCur = strsplit(basename(fileCur), split = ".",  fixed = TRUE)[[1]][1] # Start of basename is TF name
  fileListTransl.l[[fileCur]] = TFCur
}


######################
# FINAL PREPARATIONS #
######################

startLogger(par.l$file_log, par.l$log_minlevel, appenderName = "file", removeOldLog = TRUE)
printParametersLog(par.l)


################
# Collect data #
################

peaks.df = read_tsv(par.l$file_input_peaks, col_types = cols())
assertSubset(colnames(peaks.df), c("position", "D2_baseMean", "D2_l2FC", "D2_ldcSE", "D2_stat", "D2_pval", "D2_padj")) #, "vst_diff"))

summary.df = tribble(~TF_name,
                     ~Pos_l2FC,
                     ~Mean_l2FC,
                     ~Median_l2FC,
                     ~Mode_l2FC,
                     ~Ttest_pval,
                     ~Mode_skewness,
                     ~T_statistic,
                     ~TFBS_num
)


nTFs = length(fileListTransl.l)

flog.info(paste0("Using ", nTFs, " TFs"))

for (fileCur in fileList) {
  
  TF = fileListTransl.l[[fileCur]]
  stopifnot(!is.null(TF))
  if (file.exists(fileCur)) {
    stat = readRDS(fileCur)
    
    summary.df[TF,"TF_name"]         = stat$TF
    summary.df[TF,"Pos_l2FC"]    = stat$Pos_l2FC
    summary.df[TF,"Mean_l2FC"]   = stat$Mean_l2FC
    summary.df[TF,"Median_l2FC"] = stat$Median_l2FC
    summary.df[TF,"Mode_l2FC"]   = stat$Mode_l2FC
    summary.df[TF,"Ttest_pval"]  = stat$Ttest_pval
    summary.df[TF,"Mode_skewness"]   = stat$Modeskewness
    summary.df[TF,"T_statistic"]     = stat$T_statistic
    summary.df[TF,"TFBS_num"]        = stat$TFBS_num
    
  }
  
}

# TODO: WHy same values for both TF

nRowsSummary = nrow(summary.df)

flog.info(paste0(" Imported ", nRowsSummary, " TFs out of a list of ", nTFs, ". Missing: ", nTFs - nRowsSummary))

# Replace p-values of 0 with the smallest p-value on the system
summary.df$Ttest_pval[summary.df$Ttest_pval == 0] = par.l$min_pValue

# TODO: was previously assigned to modeNum, but why?
mode_peaks = mlv(round(peaks.df$D2_l2FC, 2), method = "mfv", na.rm = TRUE)

summary.df = summary.df %>%
              dplyr::mutate(
                  adj_pvalue = p.adjust(Ttest_pval, method = "BH"),
                  Diff_mean  = Mean_l2FC   - mean  (peaks.df$D2_l2FC, na.rm = TRUE), 
                  DiffMedian = Median_l2FC - median(peaks.df$D2_l2FC, na.rm = TRUE),
                  Diff_mode  = Mode_l2FC - mode_peaks[[1]],    
                  Diff_skew  = Mode_skewness - mode_peaks[[2]])  %>%
              na.omit(summary.df)


# Loop through summary files and use the TFBS_num column




# Automatically calculate the significance thresholds
# Reverse Ivans heuristic approach earlier
threshold1 = par.l$FDR_threshold / nTF
min_T_stat = qt(threshold1/2, median(summary.df$TFBS_num), lower.tail = FALSE)



# Filter rows
plot_thr.df = summary.df %>%
                filter(Diff_mean > par.l$plot_min_diffMean | Diff_mean < -par.l$plot_min_diffMean)  %>%
                filter(abs(T_statistic) > min_T_stat)



TF_volcano = ggplot() +
  geom_point(aes(x = summary.df$Diff_mean,
                 y = abs(summary.df$T_statistic),
                 label = summary.df$TF_name),size = 1)   +
  geom_vline(xintercept = 0, size = 0.7,
             linetype = "longdash", color = "blue") +
  geom_hline(yintercept = min_T_stat, size = 0.7,
             linetype = "longdash", color = "red") +
  geom_text_repel(aes(x = plot_thr.df$Diff_mean,
                      y = abs(plot_thr.df$T_statistic),
                      label = plot_thr.df$TF_name),size = 2.5,
                  segment.size = 0.5,box.padding = unit(0.05,"lines")) +
  ylab("Absolute T-statistic") +
  xlab(paste0("Mean(TF distr) - mean(peaks)")) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(face = "bold", color = "black", size = 20),
        axis.title.x = element_text(face = "bold", colour = "black", size = 24,margin = margin(25,0,0,0)),
        axis.title.y = element_text(face = "bold", colour = "black", size = 24,margin = margin(0,25,0,0)),
        axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.1,0.9),
        legend.justification = "center",
        legend.title = element_blank())

ggsave(plot = TF_volcano, filename = par.l$file_output_volcanoPlot,width = 12, height = 8, useDingbats = FALSE, dpi = 600)


#write.table(d, file = file_output_table, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
write_tsv(summary.df, par.l$file_output_table) # TODO: check the dec = "." parameter

end.time  <-  Sys.time()
message(" Finished execution. TOTAL RUNNING TIME: ", round(end.time - start.time, 1), " ", units(end.time - start.time),"\n")

