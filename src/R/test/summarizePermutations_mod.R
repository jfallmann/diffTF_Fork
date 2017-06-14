## README: this script is to remove CG bias by comparing permutations from the same CG bin 

start.time  <-  Sys.time()

#########################
# LIBRARY AND FUNCTIONS #
#########################

source("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/functions.R")

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "lsr", "ggrepel", "checkmate", "tools"), verbose = TRUE)


###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$FDR_threshold = c(0.01, 0.05,0.1,0.2, 0.5, 0.7)
par.l$volcanoPlot_height = 9
par.l$volcanoPlot_width  = 12
par.l$volcanoPlot_dpi  = 600
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.2, 0.5, 0.8)
par.l$excludeSmallCohensPrinting = FALSE
par.l$log_minlevel = "INFO"

############################################
# READ AND VALIDATE COMMAND LINE ARGUMENTS #
############################################
args <- commandArgs(trailingOnly = TRUE)

HOCOMOCO_pwm_dir = "/g/scb/zaugg/berest/ATACSeq/Analysis/Ackermann/HOCOMOCO.human/PWMscan/"
command = paste0("ls  ",HOCOMOCO_pwm_dir,"* | awk -F\"_\" \'{print $1}\' | awk -F\"/\" \'{print $NF}\'")
TFs <- system(command,intern = TRUE)
samples_TF_act = c()
for ( x in TFs ) {
  iter = paste0("/scratch/carnold/CLL/TF_act/output/TF-SPECIFIC/",x,"/extension100/",x,".CLL.permutationSummary.tsv" )
  samples_TF_act = append(samples_TF_act, iter)
}
samples_TF_act = paste(samples_TF_act, collapse = ",")
TFs = paste(TFs, collapse = ",")


args=c(
  samples_TF_act,
"/scratch/carnold/CLL/TF_act/output/TEMP/extension100/conditionComparison.rds",
"/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/CLL.all.volcano.l2FC.removedCGbias.AR.tsv",
"/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/CLL.all.volcano.l2FC.removedCGbias.AR.pdf",
"/scratch/carnold/CLL/TF_act/output/Logs_and_Benchmarks/summarizePermutations.R.log"
)


if (length(args) != 5) {
  stop("Expecting 5 arguments but found ", length(args),". Exiting.")
} else {
  par.l$files_input_permResults   = args[1]
  par.l$file_input_condCompDeSeq  = args[2]
  par.l$file_output_summary       = args[3]
  par.l$file_plotVolcano          = args[4]
  par.l$file_log                  = args[5]
}

#####################
# VERIFY PARAMETERS #
#####################

assertDirectoryExists(dirname(par.l$file_log), access = "w")
assertCharacter(par.l$files_input_permResults, len = 1)

for (dirname in unique(c(dirname(par.l$file_output_summary), dirname(par.l$file_plotVolcano)))) {
  
  if (!testDirectoryExists(dirname)) {
    dir.create(dirname, recursive = TRUE)
  } else {
    assertDirectoryExists(dirname, access = "w")
  }
  
}

assertFileExists(par.l$file_input_condCompDeSeq, access = "r")

######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel, appenderName = "file", removeOldLog = TRUE)
printParametersLog(par.l)


conditionComparison = readRDS(par.l$file_input_condCompDeSeq)

# Assemble the final table and collect permutation information from all TFs
output.global.TFs = data_frame(TF = character(), 
                               weighted_mean = numeric(),
                               weighted_Tstat = numeric(), 
                               weighted_CD = numeric(), 
                               weighted_median = numeric(), 
                               TFBS = numeric(),
                               Cohend_factor = character()
) 

fileList = strsplit(par.l$files_input_permResults, split = ",", fixed = TRUE)[[1]]
nTF = length(fileList)
for (fileCur in fileList) {
  
  assertFileExists(fileCur, access = "r")
  resultsCur.df =  read_tsv(fileCur, col_names = TRUE, col_types = cols())
  assertIntegerish(nrow(resultsCur.df), lower = 1, upper = 1)
  output.global.TFs = rbind(output.global.TFs, resultsCur.df[1,])
}


output.global.TFs$Cohend_factor = ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[1], par.l$classes_CohensD[1], 
                                         ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[2] , par.l$classes_CohensD[2], 
                                                ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[3], par.l$classes_CohensD[3], par.l$classes_CohensD[4])))

output.global.TFs$Cohend_factor = factor(output.global.TFs$Cohend_factor, levels = par.l$classes_CohensD, labels = seq_len(length(par.l$classes_CohensD)))



write_tsv(output.global.TFs, path = par.l$file_output_summary, col_names = TRUE)



########
# PLOT #
########

# Automatically calculate the significance thresholds
# Reverse Ivans heuristic approach earlier

par.l$FDR_threshold = seq(0.05, 1, 0.05)

T_stat = c()
for (FDRCur in par.l$FDR_threshold) {
    threshold1 = FDRCur /  nTF
  threshold1 = FDRCur
  T_stat = c(T_stat, qt(threshold1/2, median(output.global.TFs$TFBS,  na.rm = TRUE), lower.tail = FALSE))
  
  output.global.TFs[, paste0("FDR_", FDRCur)] = qt(threshold1/2, output.global.TFs$TFBS, lower.tail = FALSE)
}

output.global.TFs$pvalues = 2 * pt(abs(output.global.TFs$weighted_Tstat), output.global.TFs$TFBS, lower.tail = FALSE)
res.l = performIHW(output.global.TFs$pvalues, output.global.TFs$weighted_mean, alpha = 0.2, covariate_type = "ordinal",
                       nbins = 3, m_groups = NULL, quiet = TRUE, nfolds = 5L,
                       nfolds_internal = 5L, nsplits_internal = 1L, lambdas = "auto",
                       seed = 1L, distrib_estimator = "grenander", lp_solver = "lpsymphony",
                       adjustment_type = "BH", return_internal = FALSE, doDiagnostics = TRUE, pdfFile = NULL, verbose = TRUE)

res.df = as.data.frame(res.l$ihwResults)
output.global.TFs$FDR.0.2.adjpval = res.df$adj_pvalue 

min_T_stat = min(T_stat)
min_T_stat = 1.96
min_T_stat = 2.7

#min_pval = 0.05

# if (par.l$excludeSmallCohensPrinting) {
#   plot_thr_df = output.global.TFs[output.global.TFs$Cohend_factor != "1", ]
# } else {
#   plot_thr_df = output.global.TFs
# }

plot_thr_df = output.global.TFs[ which((output.global.TFs$Cohend_factor != "1"  | abs(output.global.TFs$weighted_Tstat) > min_T_stat ) & 
                                   abs(output.global.TFs$weighted_mean) > 0.02), ] 

## plot a distribtuion 
axis_colour  = "black"
panel_colour = "white"


labelsCohensD = c()
for (i in seq_len(length(par.l$classes_CohensD))) {
  if (i != length(par.l$classes_CohensD)) {
    threshold = par.l$thresholds_CohensD[i]
    signCur = "<"
  } else {
    threshold = par.l$thresholds_CohensD[i - 1]
    signCur = ">="
  }
  
  labelCur = paste0(par.l$classes_CohensD[i], " (", signCur, " ", threshold, ")")
  labelsCohensD = c(labelsCohensD, labelCur)
}

classes_CohenD = sort(as.numeric(unique(output.global.TFs$Cohend_factor)))


TF_volcano1 = ggplot() +
  geom_point(aes(x = output.global.TFs$weighted_mean,
                 y = abs(output.global.TFs$weighted_Tstat),
                 size = output.global.TFs$Cohend_factor))  +
  # geom_hline(yintercept = T_stat, size = 0.5,
  #            linetype = "longdash", color = "red") +
  # annotate(geom="text", label= par.l$FDR_threshold, x=Inf, y = T_stat, vjust = 0, color = "red", size = 3, hjust = 1) +

  geom_hline(yintercept = min_T_stat, size = 0.5, linetype = "longdash", color = "darkgreen") + 
  geom_label_repel(aes(x = plot_thr_df$weighted_mean,
                       y = abs(plot_thr_df$weighted_Tstat),
                       size = plot_thr_dfs$Cohend_factor,
                       label = plot_thr_df$TF), 
                   size = 1.75, segment.size = 0.25, box.padding = unit(0.05,"lines")) +
  ylab("Absolute T-statistic\n(with various FDR thresholds)") +
  xlab(paste0("TF 'activity' ", conditionComparison)) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_text(face = "bold", colour = "black", size = 15, margin = margin(25,0,0,0)),
        axis.title.y = element_text(face = "bold", colour = "black", size = 15, margin = margin(0,25,0,0)),
        axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.grid.major = element_line(color = panel_colour),
        panel.grid.minor = element_line(color = panel_colour),
        panel.background =  element_rect(fill = panel_colour, colour = panel_colour),
        legend.justification = "center", 
        legend.background = element_rect(fill = panel_colour, colour = panel_colour),
        legend.key = element_rect(fill = panel_colour, colour = panel_colour))  +
  scale_size_manual(values = classes_CohenD, labels = labelsCohensD[classes_CohenD]) +
  #  scale_colour_gradient2(low = "blue", mid = "green" , high = "red", midpoint = 0.5) + 
  labs(size = "Cohen's D")  
TF_volcano1



ggsave(plot = TF_volcano1, height = par.l$volcanoPlot_height, width = par.l$volcanoPlot_width, dpi = par.l$volcanoPlot_dpi, filename = par.l$file_plotVolcano, useDingbats = FALSE)

