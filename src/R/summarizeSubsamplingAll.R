start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################

source("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/functions.R")

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
packages = c("tidyverse", "futile.logger", "checkmate", "tools", "methods", "ggrepel", "purrr", "plyr", "pheatmap", "reshape2", "viridis", "jsonlite", "scales")
checkAndLoadPackages(packages, verbose = TRUE)


###################
#### PARAMETERS ###
###################

par.l = list()


par.l$rootFolder = "/scratch/carnold/CLL/TF_act_downsampling0.25"
par.l$stepsize = 5
par.l$nRepetitionsPerStepMax = 20
par.l$file_output_log = "summarizeSubsamplingAll.R.log"

par.l$file_TF_list = paste0("/scratch/carnold/CLL/TF_act/input/names.TF.all.tsv")

assertFileExists(par.l$file_input_config)

outputFolder = paste0(par.l$rootFolder, "/output/FINAL_OUTPUT/subsampleSamples")




par.l$file_input_config  = paste0(par.l$rootFolder, "/input/config.json")
config.l = fromJSON(file(par.l$file_input_config), simplifyVector = FALSE)

assertCharacter(config.l$samples$summaryFile, len = 1, min.chars = 1)
file_sampleData = config.l$samples$summaryFile
assertFileExists(file_sampleData)

assertIntegerish(config.l$par_general$regionExtension, lower = 1, len = 1)
extensionSize = config.l$par_general$regionExtension

assertCharacter(config.l$par_general$conditionComparison, len = 1, min.chars = 1)
conditionComparison = config.l$par_general$conditionComparison
assertCharacter(conditionComparison, len = 1, min.chars = 3)


assertCharacter(config.l$par_general$analysisName, len = 1, min.chars = 1)
analysisName = config.l$par_general$analysisName
  
#################
### SUMMARIZE ###
#################

startLogger(par.l$file_output_log, par.l$log_minlevel,  removeOldLog = TRUE)

allTF.df = read_tsv(par.l$file_TF_list, col_names = FALSE)
allTF = allTF.df$X1
nTF = length(allTF)

files = createFileList(paste0(par.l$rootFolder, "/output/FINAL_OUTPUT/extension", extensionSize), ".StatsAllTF.tsv", recursive = TRUE, ignoreCase = FALSE, verbose = TRUE)

assertCharacter(files, min.len = 109)


permSummaryAll.l = list()

# get min and max TF activity values
activityMin = 0
activityMax = 0

for (fileCur in files) {
  
  #flog.info(paste0("Read file ", fileCur))
  nameSubsampleCurVec = strsplit(basename(fileCur), split = ".", fixed = TRUE)[[1]]
  nameSubsampleCur = gsub(pattern = "subsample", "", nameSubsampleCurVec[2])
  stats.subsample.df = read_tsv(fileCur, col_names = TRUE)
  permSummaryAll.l[[nameSubsampleCur]] = stats.subsample.df
  
  activityMinCur = min(stats.subsample.df$weighted_mean)
  activityMaxCur = max(stats.subsample.df$weighted_mean)
  
  if (activityMinCur < activityMin) {
    activityMin = activityMinCur
  }
  
  if (activityMaxCur > activityMax) {
    activityMax = activityMaxCur
  }
  
}

sampleSizeStr = strsplit(names(permSummaryAll.l), split = "_")

sampleSizes = as.numeric(unlist(map(sampleSizeStr, 1)))
sampleSizesSortedAsc = sort(unique(sampleSizes))

repetitions = as.numeric(unlist(map(sampleSizeStr, 2)))

###############################
# 1. TF-SPECIFIC SUMMARY PLOT #
###############################



# Construct label for x axis sto replace default label
xAxisLabel = c()
subsamples.l = permuteSampleTable(file_sampleData, conditionComparison, "Treatment", par.l$stepsize, par.l$nRepetitionsPerStepMax)
for (sampleSizeCur in sampleSizesSortedAsc) {
  distCur = table(subsamples.l[[paste0(sampleSizeCur, "_1")]]$conditionSummary)
  stopifnot(length(distCur) > 0)
  xAxisLabel = c(xAxisLabel, paste0(distCur[2], "+", distCur[1], " (", round(sum(distCur)/84 * 100,0), "%)"))
}


pdf(paste0(outputFolder, "/subsamplingResults.pdf"))


resultsPerTF.l = list()
for (TFCur in allTF) {
  
  rowTF = which(permSummaryAll.l[[1]]$TF == TFCur)
  stopifnot(length(rowTF) == 1)
  
  
  meanValues = map(permSummaryAll.l,"weighted_mean")
  
  result.df = tibble(sampleSize = sampleSizes,
                     weightedMean = unlist(map(meanValues, rowTF)))
  
  resultsPerTF.l[[TFCur]] = result.df
  
  p <- ggplot(result.df , aes(x = sampleSize, y = weightedMean, group = sampleSize))
  p <- p + geom_boxplot() + geom_jitter(alpha = 0.1, height = 0, color = "black")
  p <- p + .getThemeForGGPlot()
  p <- p + geom_hline(yintercept = 0, color = "gray")
  p <- p + scale_y_continuous(limits = c(activityMin,activityMax))
  p <- p + ggtitle(TFCur)
  p <- p + scale_x_continuous(breaks=sampleSizesSortedAsc, labels = xAxisLabel)
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p <- p + ylab("TF activity: Weighted mean")
  p <- p + xlab("Sample size total (unmutated + mutated)")
  plot(p)
}

dev.off()

saveRDS(resultsPerTF.l, paste0(outputFolder, "/subsamplingResults.rds"))

###########################
# 2. OVERALL SUMMARY PLOT #
###########################



# Summarize variability
variabilitySummaryPerTF.l = list()
for (TFCur in allTF) {
  
  summaryCur.l = list()
  
  for (sampleSizeCur in sampleSizesSortedAsc) {
    
    data                                          = dplyr::filter(resultsPerTF.l[[TFCur]], sampleSize == sampleSizeCur)$weightedMean
    summaryCur.l[[paste0(sampleSizeCur)]]$data    = data
    
    summaryCur.l[[paste0(sampleSizeCur)]]$median  = median(data)
    summaryCur.l[[paste0(sampleSizeCur)]]$mean    = mean(data)
    summaryCur.l[[paste0(sampleSizeCur)]]$sd      = sd(data)
    summaryCur.l[[paste0(sampleSizeCur)]]$mad     = mad(data)
    summaryCur.l[[paste0(sampleSizeCur)]]$madm    = summaryCur.l[[paste0(sampleSizeCur)]]$mad / summaryCur.l[[paste0(sampleSizeCur)]]$median
    summaryCur.l[[paste0(sampleSizeCur)]]$cv      = summaryCur.l[[paste0(sampleSizeCur)]]$sd  / summaryCur.l[[paste0(sampleSizeCur)]]$mean
    summaryCur.l[[paste0(sampleSizeCur)]]$Q95_Q05 = abs(quantile(data, probs = c(.95), na.rm = TRUE) - quantile(data, probs = c(.05), na.rm = TRUE))
    
    # Count how often the activity is > 0
    summaryCur.l[[paste0(sampleSizeCur)]]$activityPos     = data > 0
    
  }
  
  # TODO: does mean of sd make more sense ?
  
  summaryCur.l[["all"]] = list()
  summaryCur.l[["all"]]$sd     = sd(unlist(map(summaryCur.l, "sd")), na.rm  = TRUE)
  summaryCur.l[["all"]]$mean   = mean(unlist(map(summaryCur.l, "mean")), na.rm  = TRUE)
  summaryCur.l[["all"]]$median = median(unlist(map(summaryCur.l, "median")), na.rm  = TRUE)
  summaryCur.l[["all"]]$mad    = mad(unlist(map(summaryCur.l, "mad")), na.rm  = TRUE)
  summaryCur.l[["all"]]$madm   = summaryCur.l[["all"]]$mad / summaryCur.l[["all"]]$median
  summaryCur.l[["all"]]$cv     = summaryCur.l[["all"]]$sd  / summaryCur.l[["all"]]$mean 
  
  activityPosGroundTruth = summaryCur.l[[paste0(max(sampleSizesSortedAsc))]]$mean > 0
  
  for (sampleSizeCur in sampleSizesSortedAsc) {
    
    summaryCur.l[[paste0(sampleSizeCur)]]$fractionActivityDirWrong = length(which(summaryCur.l[[paste0(sampleSizeCur)]]$activityPos != activityPosGroundTruth)) / length(summaryCur.l[[paste0(sampleSizeCur)]]$activityPos)
    summaryCur.l[[paste0(sampleSizeCur)]]$nActivityDirWrong = length(which(summaryCur.l[[paste0(sampleSizeCur)]]$activityPos != activityPosGroundTruth))
  }
  
  
  #map(summaryCur.l, ~length(which(.$data > 0)))
  #map(summaryCur.l, ~length(.$data))
  
  variabilitySummaryPerTF.l[[TFCur]] = summaryCur.l
  
}

saveRDS(variabilitySummaryPerTF.l, paste0(outputFolder, "/subsamplingResultsSummaryPerTF.rds"))


pdf(paste0(outputFolder, "/subsamplingResultsClassification.pdf"), height = 8, width = 8)

# Plot classification accurary whether or not TF activity is positive or negative
for (TFCur in allTF) {
  
  data = unlist(map(variabilitySummaryPerTF.l[[TFCur]], ~.$fractionActivityDirWrong))
  
  result.df = tibble(sampleSize = as.integer(names(data)),
                     fractionFalse = data)
  
  p <- ggplot(result.df , aes(x = sampleSize, y = fractionFalse))
  p <- p + geom_point() 
  p <- p + .getThemeForGGPlot()
  p <- p + scale_y_continuous(limits = c(0,1))
  p <- p + ggtitle(TFCur)
  p <- p + geom_smooth() 
  p <- p + scale_x_continuous(breaks=sampleSizesSortedAsc, labels = xAxisLabel)
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p <- p + ylab("Fraction of activities in the opposite direction")
  p <- p + xlab("Sample size total (unmutated + mutated)")
  plot(p)
}

dev.off()



pdf(paste0(outputFolder, "/subsamplingResultsClassificationSummary.pdf"), height = 10, width = 12)

activityValuesFull = unlist(map(variabilitySummaryPerTF.l, c(paste0(max(sampleSizesSortedAsc)), "mean")))
resultsAll.df = data_frame(activityTrue = activityValuesFull)

for (sampleSizeCur in sampleSizesSortedAsc) {
  
  for (takeAbsolute in c(TRUE, FALSE)) {
    
    result.df = tibble(activityTrue    = activityValuesFull,
                       fractionFalse   = unlist(map(variabilitySummaryPerTF.l, c(paste0(sampleSizeCur), "fractionActivityDirWrong"))),
                       nFalse          = unlist(map(variabilitySummaryPerTF.l, c(paste0(sampleSizeCur), "nActivityDirWrong"))))
    
    resultsAll.df[,paste0("fractionFalse_", sampleSizeCur)] = result.df$fractionFalse
    resultsAll.df[,paste0("nFalse_", sampleSizeCur)] = result.df$nFalse
    
    
    if (takeAbsolute) {
      p <- ggplot(result.df , aes(x = abs(activityTrue), y = fractionFalse))
    } else {
      p <- ggplot(result.df , aes(x = activityTrue, y = fractionFalse))
    }
    
    
    p <- p + geom_point() 
    p <- p + .getThemeForGGPlot()
    p <- p + scale_y_continuous(limits = c(0,1))
    p <- p + ggtitle(paste0("Sample size total (unmutated + mutated): ", xAxisLabel[which(sampleSizesSortedAsc == sampleSizeCur)]))
    p <- p + geom_smooth() 
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    p <- p + ylab("Fraction of times for which TF activity is in opposite direction")
    
    p <- p + xlab(ifelse(takeAbsolute, "TF activity (absolute value)", "TF activity"))
    plot(p)
    
    
  }
  
  
}



for (binsCur in c(5,10,20,30,40)) {
  
  nBins = binsCur
  resultsAll.df = transform(resultsAll.df, activityBin = cut(activityTrue, nBins, include.lowest = TRUE))
  resultsAll.df = transform(resultsAll.df, activityBinAbs = cut(abs(activityTrue), nBins, include.lowest = TRUE))
  
  binDistAll = table(resultsAll.df$activityBinAbs)
  binDist = binDistAll[which(binDistAll > 0)]
  
  nBinsReal    = length(unique(resultsAll.df$activityBin))
  nBinsAbsReal = length(unique(resultsAll.df$activityBinAbs))
  
  res_melted.df = melt(select(resultsAll.df, -activityTrue,-activityBin, -starts_with("fraction"), -ends_with("34")))
  
  res_melted2.df = ddply(res_melted.df, "activityBinAbs", mutate, percent = value/sum(value) * 100)
  
  nClasses = length(unique(res_melted2.df$variable))
  
  xAxisName   = paste0("Binned TF activity from small to large (absolute values, ", nBinsAbsReal, " non-empty bins)")
  xAxisLabels = paste0("Bin ", seq(1: nBinsAbsReal), "(", binDist, " TF)\n", sort(unique(res_melted2.df$activityBinAbs)))
  
  p <- ggplot(res_melted2.df, aes(activityBinAbs, y = percent, fill = variable))
  p <- p + geom_bar(stat='identity') 
  p <- p + .getThemeForGGPlot()
  p <- p + scale_y_continuous(name="Relative of cases for which TF activity has an inverted sign")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p <- p + scale_x_discrete(name = xAxisName, labels = xAxisLabels)
  #p <- p + scale_fill_manual(name="Subsampling size", values = rev(viridis(nClasses)), labels = xAxisLabel)
  #p <- p + scale_fill_manual(name="Subsampling size", values = rev(magma(nClasses)), labels = xAxisLabel)
  #p <- p + scale_fill_manual(name="Subsampling size", values = rev(plasma(nClasses)), labels = xAxisLabel)
  p <- p + scale_fill_manual(name="Subsampling size", values = rev(inferno(nClasses)), labels = xAxisLabel)
  plot(p)
  
  p <- ggplot(res_melted2.df, aes(activityBinAbs, y = value, fill = variable))
  p <- p + geom_bar(stat='identity') 
  p <- p + .getThemeForGGPlot()
  p <- p + scale_y_continuous(name="Absolute number of cases for which TF activity has an inverted sign")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p <- p + scale_x_discrete(name = xAxisName, labels = xAxisLabels)
  #p <- p + scale_fill_manual(name="Subsampling size", values = rev(viridis(nClasses)), labels = xAxisLabel)
  #p <- p + scale_fill_manual(name="Subsampling size", values = rev(magma(nClasses)), labels = xAxisLabel)
  #p <- p + scale_fill_manual(name="Subsampling size", values = rev(plasma(nClasses)), labels = xAxisLabel)
  p <- p + scale_fill_manual(name="Subsampling size", values = rev(inferno(nClasses)), labels = xAxisLabel)
  plot(p)
  
}


dev.off()


# Plot variability summary
# Specific for each, which TF has the least amount of variation across the different subsampling data?
TFList = names(variabilitySummaryPerTF.l)
medianVal = map(variabilitySummaryPerTF.l, "median")

cv   = unlist(map(variabilitySummaryPerTF.l, c("all", "cv")))
sd   = unlist(map(variabilitySummaryPerTF.l, c("all", "sd")))
mad  = unlist(map(variabilitySummaryPerTF.l, c("all", "mad")))
madm = unlist(map(variabilitySummaryPerTF.l, c("all", "madm")))

varSummary.df = data_frame(TF = names(cv), cv = cv, sd = sd, mad = mad, madm = madm)

pdf(paste0(outputFolder, "/subsamplingResultsSummary.pdf"), height = 20, width = 8)


p1 <- ggplot(varSummary.df, aes(x = reorder(TF, -cv), y = cv))
p1 <- p1 + geom_point(size = 0.5) 
p1 <- p1 + coord_flip()
p1 <- p1 + .getThemeForGGPlot()
p1 <- p1 + theme(axis.text.y=element_text(size=2))
p1 <- p1 + ylab("Coefficient of variation of weighted TF activity across subsamples")
p1 <- p1 + xlab("TF (sorted)")
p1 <- p1 + geom_hline(yintercept = 0, color = "gray")
plot(p1)

p2 <- ggplot(varSummary.df, aes(x = reorder(TF, -sd), y = sd))
p2 <- p2 + geom_point(size = 0.5) 
p2 <- p2 + coord_flip()
p2 <- p2 + .getThemeForGGPlot()
p2 <- p2 + theme(axis.text.y=element_text(size=2))
p2 <- p2 + ylab("Standard deviation of weighted TF activity across subsamples")
p2 <- p2 + xlab("TF (sorted)")
p2 <- p2 + geom_hline(yintercept = 0, color = "gray")
plot(p2)

p3 <- ggplot(varSummary.df, aes(x = reorder(TF, -mad), y = mad))
p3 <- p3 + geom_point(size = 0.5) 
p3 <- p3 + coord_flip()
p3 <- p3 + .getThemeForGGPlot()
p3 <- p3 + theme(axis.text.y=element_text(size=2))
p3 <- p3 + ylab("MAD of weighted TF activity across subsamples")
p3 <- p3 + xlab("TF (sorted)")
p3 <- p3 + geom_hline(yintercept = 0, color = "gray")
plot(p3)

saveRDS(varSummary.df, paste0(outputFolder, "/subsamplingResultsSummary.rds"))


# Summarize across subsample sizes across all TF
summary2.l = list()



for (varMeasureCur in c("cv", "sd", "mad", "madm")) {
  
  summary2.df = data_frame(subsampleSize = numeric(), data = numeric())
  
  for (sampleSizeCur in sampleSizesSortedAsc) {
    
    summary2.df = full_join(summary2.df, data_frame(subsampleSize = sampleSizeCur, 
                                                    data = unlist(map(variabilitySummaryPerTF.l, c(paste0(sampleSizeCur), varMeasureCur)))))
    
  }
  
  summary2.l[[varMeasureCur]] = summary2.df
  
  p <- ggplot(summary2.df , aes(x = subsampleSize, y = data, group = subsampleSize))
  p <- p + geom_boxplot() 
  p <- p + .getThemeForGGPlot()
  p <- p + scale_x_continuous(breaks=sampleSizesSortedAsc, labels = xAxisLabel)
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p <- p + ylab(paste0(varMeasureCur, " of weighted TF activity"))
  p <- p + xlab("Sample size total (unmutated + mutated)")
  plot(p)
}

saveRDS(summary2.l, paste0(outputFolder, "/subsamplingResultsSummary2.rds"))


dev.off()  

###################
# 3. VOLCANO PLOT #
###################

final_TF_Act_data =  paste0(par.l$rootFolder, "/output/FINAL_OUTPUT/extension100/CLL.all.volcano.l2FC.removedCGbias.AR.FDR5.tsv")
data = read_tsv(final_TF_Act_data, col_names = TRUE)

data.df = data_frame(weightedMean = data$weighted_mean, sd = sd, TF = data$TF)
plot_thr_df = data.df[which(data.df$sd > 0.02),]

TF_volcano1 = ggplot() +
  geom_point(aes(x = data.df$weightedMean, y = data.df$sd))  +
  geom_hline(yintercept = 0.02, size = 0.5, linetype = "longdash", color = "darkgreen") + 
  geom_label_repel(aes(x = plot_thr_df$weightedMean,
                       y = plot_thr_df$sd,
                       label = plot_thr_df$TF), 
                   size = 1.75, segment.size = 0.25, box.padding = unit(0.05,"lines")) +
  ylab("Variability (sd) of weighted TF activity estimates across all subsamples") +
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
        legend.key = element_rect(fill = panel_colour, colour = panel_colour))



end.time  <-  Sys.time()
message(" Finished execution of script. TOTAL RUNNING TIME: ", round(end.time - start.time, 1), " ", units(end.time - start.time),"\n")






# Check the subsampling results and where to find the permutation results

subsamples.l = permuteSampleTable(file_sampleData, conditionComparison, "Treatment", par.l$stepsize, par.l$nRepetitionsPerStepMax)
allFractions = c("1", "0.75", "0.50",  "0.25", "0.125", "0.06", "0.02", "0.01")



nameFull = "34_1"

allMeasures  = c("weighted_mean", "weighted_Tstat", "weighted_CD", "weighted_median")

# TODO: Code for 3 bins, equally bin all TF activities into three bins
nBins = 3
binsChar = as.character(1:nBins)

finalMatrices.l = list()
finalMatricesClassification1Bin.l = list()
finalMatricesClassification3Bin.l = list()
subsampleSizesSorted = sort(as.numeric(unique(unlist(map(strsplit(names(subsamples.l), split = "_", fixed = TRUE), 1)))))
for (measureCur in allMeasures) {
  finalMatrices.l[[measureCur]] = finalMatricesClassification1Bin.l[[measureCur]] = matrix(nrow = length(allFractions), ncol = length(subsampleSizesSorted))
  rownames(finalMatrices.l[[measureCur]]) =  rownames(finalMatricesClassification1Bin.l[[measureCur]]) = allFractions
  colnames(finalMatrices.l[[measureCur]]) =  colnames(finalMatricesClassification1Bin.l[[measureCur]]) = subsampleSizesSorted 
  
  finalMatricesClassification3Bin.l[[measureCur]] = list()
  for (binCur in binsChar) {
    finalMatricesClassification3Bin.l[[measureCur]][[binCur]] = matrix(nrow = length(allFractions), ncol = length(subsampleSizesSorted))
    rownames(finalMatricesClassification3Bin.l[[measureCur]][[binCur]]) = allFractions
    colnames(finalMatricesClassification3Bin.l[[measureCur]][[binCur]]) = subsampleSizesSorted 
  }
}

# Stores the results of the full dataset to compare against.
truth.l = list()



for (fractionCur in allFractions) {
  
  for (measureCur in allMeasures) {
    
    root = "/scratch/carnold/CLL"
    rootNew = ifelse(fractionCur == "1", paste0(root, "/TF_act"), paste0(root, "/TF_act_downsampling", fractionCur))
    
    remainder = paste0("/output/FINAL_OUTPUT/extension", extensionSize, "/subsamplingSamples/CLL.permResultsAllSummary.rds")
   
    filename = paste0(rootNew, remainder)
    # Create the summary rds file if it does not exist yet
    if (!file.exists(filename)  | file.exists(filename)) {
      
      permSummaryAll.l = list()
      
      for (nameSubsampleCur in names(subsamples.l)) {
        stats_filename = paste0(rootNew, "/output/FINAL_OUTPUT/extension", extensionSize, "/subsamplingSamples/", analysisName, ".subsample", nameSubsampleCur , ".StatsAllTF.tsv")
        
        if (!file.exists(stats_filename)) {
          
          stats.df.subsample = list()
          if (fractionCur != "1") warning(paste0("Could not find results for ", nameSubsampleCur, " and ", fractionCur))
          
        } else {
          stats.df.subsample = read_tsv(stats_filename, col_names = TRUE, col_types = cols())
        }

        permSummaryAll.l[[nameSubsampleCur]] = stats.df.subsample
      }
      
      
      names(permSummaryAll.l) = names(subsamples.l)
      
      saveRDS(permSummaryAll.l, file = filename)
      
      
    }
    
  
   results.l = readRDS(filename)
    
    weighted_mean.l = map(results.l, measureCur)
    subsampleSizes = sort(unique(unlist(map(strsplit(names(weighted_mean.l), split = "_", fixed = TRUE), 1))))
    
    if (fractionCur == "1") {
      truth.l[[measureCur]] = weighted_mean.l 
    } 
    
    
    # Calculate the median
    valuesPerSubsampleSize.l = list()
    valuesPerSubsampleSizeClassification1Bin.l = list()
    valuesPerSubsampleSizeClassification3Bin.l = list()
    
    for (binCur in binsChar) {
      valuesPerSubsampleSizeClassification3Bin.l[[binCur]] = list()
    }
    
    for (nameCur in names(weighted_mean.l)) {
      
      
      subsampleSizeCur = strsplit(nameCur, split = "_", fixed = TRUE)[[1]][1]
      
      # Set all values to NA
      if (is.null(weighted_mean.l[[nameCur]])) {
        
        corCur = NA
        correct = NA
        
        for (binCur in binsChar) {
          
          valuesPerSubsampleSizeClassification3Bin.l[[binCur]][[subsampleSizeCur]] = c(valuesPerSubsampleSizeClassification3Bin.l[[binCur]][[subsampleSizeCur]], NA)
        }
        
      } else {
        
        corCur = cor(weighted_mean.l[[nameCur]], truth.l[[measureCur]][[nameFull]])
        

        # Should always be equal due to the constant nameFull, even though it is recomputed here all the time
        quantileThresholds = quantile(abs(weighted_mean.l[[nameCur]]), seq(1/nBins, 1, 1/nBins)) 
        TF_quantile_Indexes = list()
        for (i in 1:nBins) {
          
          if (i == 1) {
            TF_quantile_Indexes[[as.character(i)]] =  which(abs(weighted_mean.l[[nameCur]]) <=  quantileThresholds[i])
          } else {
            TF_quantile_Indexes[[as.character(i)]] =  which(abs(weighted_mean.l[[nameCur]]) >  quantileThresholds[i-1] & abs(weighted_mean.l[[nameCur]]) <=  quantileThresholds[i])
            
          }
        }
  
        # Classify how often we predict the correct sign of the TF activity
        correct1Bin = length(which((weighted_mean.l[[nameCur]] > 0 & truth.l[[measureCur]][[nameFull]] > 0) | (weighted_mean.l[[nameCur]] < 0 & truth.l[[measureCur]][[nameFull]] < 0))) / nTF
        
        
        for (binCur in binsChar) {
          
          indexesCur = TF_quantile_Indexes[[binCur]]
          valuesCur = weighted_mean.l[[nameCur]] [indexesCur]
          correct = length(which((valuesCur > 0 & truth.l[[measureCur]][[nameFull]][indexesCur] > 0)  | (valuesCur < 0 & truth.l[[measureCur]][[nameFull]] [indexesCur] < 0))) / length(TF_quantile_Indexes[[binCur]])
          
          stopifnot(is.finite(correct))
          
          
          valuesPerSubsampleSizeClassification3Bin.l[[binCur]][[subsampleSizeCur]] = c(valuesPerSubsampleSizeClassification3Bin.l[[binCur]][[subsampleSizeCur]], correct)
        }
        
        
      }
     
      
      valuesPerSubsampleSize.l[[subsampleSizeCur]] = c(valuesPerSubsampleSize.l[[subsampleSizeCur]], corCur)
      
      valuesPerSubsampleSizeClassification1Bin.l[[subsampleSizeCur]] = c(valuesPerSubsampleSizeClassification1Bin.l[[subsampleSizeCur]], correct1Bin)
      

      
    }
    
    medianValues = unlist(map(valuesPerSubsampleSize.l, median, na.rm = TRUE))
    medianValuesClassification = unlist(map(valuesPerSubsampleSizeClassification1Bin.l, median, na.rm = TRUE))
   
 
    finalMatrices.l[[measureCur]][fractionCur, ] = medianValues
    finalMatricesClassification1Bin.l[[measureCur]][fractionCur, ] = medianValuesClassification
    
    
    for (binCur in binsChar) {
      
      finalMatricesClassification3Bin.l[[measureCur]][[binCur]][fractionCur, ] = unlist(map(valuesPerSubsampleSizeClassification3Bin.l[[binCur]], median, na.rm = TRUE))

    }
  }
  

}


# Determine median number of reads

median_nReads = c()

for (fractionCur in allFractions) {
  
  if (fractionCur == "1") {
    
    dirFull = "/scratch/carnold/CLL/TF_act/input/bam"
    files = createFileList(dirFull, paste0("*.bam$"), recursive = FALSE, ignoreCase = FALSE, verbose = FALSE)
    
    stopifnot(length(files) == 88)
    # We have four extra fiels we did not analyze actually, but the difference in numbers will be very minor so ok for now...
    
  } else {
    
    downsampleDir = "/scratch/carnold/CLL/Downsampling/output"
    files = createFileList(downsampleDir, paste0("*", fractionCur, "*.bam$"), recursive = FALSE, ignoreCase = FALSE, verbose = FALSE)
    
    stopifnot(length(files) == 84)
  }
  
  
  
  
  nReads = c()
  for (fileCur in files) {
    command =  paste0("samtools idxstats ", fileCur, " | cut -f3 | awk '{s+=$1} END {print s}'")
    nReads[fileCur] = as.numeric(system(command, intern = TRUE))
  }
  
  median_nReads[fractionCur] = median(as.numeric(nReads))
  
}

rowLabel = paste0(round(as.numeric(allFractions) * 100, 1), "%: ", .prettyNum(round(median_nReads/1000000, 1)), " m")


resultsAll.l = list()
resultsAll.l$medianNoReads = median_nReads
resultsAll.l$finalMatricesClassification1Bin = finalMatricesClassification1Bin.l
resultsAll.l$finalMatricesClassification3Bin = finalMatricesClassification3Bin.l
resultsAll.l$finalMatrices = finalMatrices.l
saveRDS(resultsAll.l, file = "/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/subsamplingSamples/resultsALL.rds")

pdf("/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/subsamplingSamples/heatmapSummaryAll.pdf")
pdf("heatmapSummaryAll.pdf")



for (classification in c(0,1,nBins)) {
  
  if (classification == 1) {
    data.l = finalMatricesClassification1Bin.l
  } else if (classification == nBins) {
    data.l = finalMatricesClassification3Bin.l
  } else {
    data.l = finalMatrices.l
  }
  
  
  for (measureCur in allMeasures) {
    
    fontsize = 10
    
 
    if (classification == 1 | classification == nBins) {
      breaks_num = seq(0.5, 1, 0.01)
      color = viridis(length(breaks_num) + 1, option = "D")
      color = colorRampPalette( c("red", "yellow", "darkgreen"), space = "rgb")(length(breaks_num) + 1)
      label = paste0("% of correctly classified (positive or \nnegative, based on full data) across all TF\n for ", measureCur)
    } else {
      breaks_num = seq(0, 1, 0.01)
      color = viridis(length(breaks_num) + 1, option = "D")
      color = colorRampPalette( c("red", "yellow", "darkgreen"), space = "rgb")(length(breaks_num) + 1)
      label = paste0("Correlation with full data\nacross all TF for ", measureCur, "\n")
    }
    
    dataToPlot = list()
    if (classification == nBins) {
      
      for (binCur in binsChar) {
        dataToPlot[[binCur]] = data.l[[measureCur]][[binCur]]
      }
    } else {
      dataToPlot[["all"]] = data.l[[measureCur]]
    }
    
    
    
    for (binCur in names(dataToPlot))  {
      
      labelReal = label
      if (classification == nBins) {
        labelReal = paste0(label, " (bin ", binCur, ")")
      }
      
      pheatmap(dataToPlot[[binCur]], kmeans_k = NA, breaks = breaks_num, border_color = "grey60", color = color,
                cellwidth = 25, cellheight = 25, scale = "none", cluster_rows = FALSE,
                cluster_cols = FALSE, clustering_distance_rows = "euclidean",
                cutree_rows = NA, cutree_cols = NA,
                legend = TRUE, legend_breaks = NA,
                legend_labels = NA, annotation_row = NA, 
                #annotation_col = sampleData.df[,c("gender","mutation","batch_proc","time","age_coll")],
                annotation = NA, annotation_colors = NA, annotation_legend = TRUE,
                annotation_names_row = TRUE, annotation_names_col = TRUE,
                drop_levels = TRUE, show_rownames = T, show_colnames = T, main = labelReal,
                fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize,
                display_numbers = T, number_format = "%.2f", number_color = "black",
                fontsize_number = 0.8 * fontsize, gaps_row = NULL, gaps_col = NULL,
                labels_row = rowLabel, labels_col = xAxisLabel, 
                width = 5,
                filename = NA)
    }
    
   
  }
}


dev.off()
