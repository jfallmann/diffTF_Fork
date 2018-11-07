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
createDebugFile(snakemake)

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "lsr", "DESeq2",  "matrixStats", "ggrepel", "checkmate", "tools", "grDevices", "locfdr", "pheatmap"), verbose = FALSE)

# TODO: grdevices realyl needed? Only for adjustcolor so far

###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$significanceThresholds  = c(0.001, 0.01, 0.05,0.1,0.2)

par.l$expressionThreshold = 2
par.l$cohensDThreshold = 0.1
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.1, 0.5, 0.8)
par.l$regressionMethod = "glm"
par.l$filter_minCountsPerCondition = 5
par.l$log_minlevel = "INFO"
par.l$circularPlot_minDimensions  = 12
par.l$corMethod = "pearson"

par.l$extension_x_limits = 0.025 
par.l$extension_x_limits = 0.15 # 15 % x axis extension, regardless of the limits
par.l$extension_y_limits = 1.1
par.l$plot_grayColor = "grey50"

par.l$maxTFsToDraw = 150

par.l$pseudocountLogTransform = 0.0001

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

transform_yValues <- function(values, addPseudoCount = TRUE, onlyForZero = TRUE) {
    
    # Should only happen with the permutation-based approach
    zeros = which(values == 0)
    if (length(zeros) > 0 & addPseudoCount) {
        values[zeros] = 1 / par.l$nPermutations
    }
    
    -log10(values)
}

transform_yValues_caption <- function() {
  "-log10"
}

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
assertSubset(c("", "summary", "circularPlot", "diagnosticPlots", "plotsRDS"), names(snakemake@output))

par.l$file_output_summary  = snakemake@output$summary
par.l$file_plotCircular    = snakemake@output$circularPlot
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


allDirs = c(dirname(par.l$file_output_summary), 
            dirname(par.l$file_plotCircular),
            dirname(par.l$files_plotDiagnostic),
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


#############
# FUNCTIONS #
#############
plotCircular <- function(dataCur, par.l, showClasses, significanceThresholdCur, conditionComparison, variableXAxis) {
  
  if (par.l$plotRNASeqClassification) {
    # Filter by classes to show
    dataCur = filter(dataCur, classification %in% showClasses)
    dataCur$classification = factor(dataCur$classification, levels = showClasses)
    
    colorCategoriesCur = par.l$colorCategories[which(names(par.l$colorCategories) %in% showClasses)]
  }
  

  # TODO: When adjust p-value? Before or after filtering?
  dataCur = mutate(dataCur, sign = pvalueAdj <= significanceThresholdCur)
  
  
  # Limits for x-axis. Multiple by a factor > 1 to account for the white legend part in which no point should be located
  # For this, extend the x axis limit so that radial positions are smaller and do not cross the border of the legend
  # Enforce a minimum of above 1 for enrichment analysis
  limit_max = max(abs(dataCur[,variableXAxis])) * (1 + par.l$extension_x_limits)
  limit_min = 0
  
  # Determine the value for 1 degree radial positions. 180 is the max
  radial_coeff = 180/limit_max
  
  # Calculate the radial positions
  dataCur$radial = 0
  index_neg = which(dataCur$weighted_meanDifference < 0)
  index_pos = which(dataCur$weighted_meanDifference > 0)
  
  # extracting the labels from object in order to add them manually later
  dataCur$radial[index_neg] = 180 + (radial_coeff*abs(unlist(dataCur[, variableXAxis])[index_neg]))
  dataCur$radial[index_pos] = 180 - (radial_coeff*abs(unlist(dataCur[, variableXAxis])[index_pos]))
  
  # Which measure to use for the y value?
  signThresholdPlot = transform_yValues(significanceThresholdCur, addPseudoCount = FALSE)
  
  # Dont display TFs that are deemed non-significant
  ggrepel_df = filter(dataCur, sign == TRUE)
  # ggrepel_df = filter(ggrepel_df, weighted_CD > par.l$cohensDThreshold)
  
  
  ylimit = max(signThresholdPlot, sign(max(dataCur$yValue)) * ceiling(abs(max(dataCur$yValue)))) * par.l$extension_y_limits
  
  labels        = list()
  startNo = sign(min(dataCur$yValue)) * ceiling(abs(min(dataCur$yValue)))
  endNo   = ylimit
  stepsize = 1
  labels$breaks = round(seq(startNo, endNo, stepsize),0)
  
  
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

  p1 = p1 + geom_hline(yintercept = ylimit, size = 0.8, linetype = "solid", color = par.l$plot_grayColor, alpha = .9) +
    # Make the white area separate into two to avoid strange artefacts
    geom_rect(data = NULL,aes(xmin = 360 - width_whiteArea, xmax = 360 , ymin = -Inf, ymax = ylimit + 0.1), alpha = 1, fill = "white") +
    geom_rect(data = NULL,aes(xmin = 0 , xmax = width_whiteArea , ymin = -Inf, ymax = ylimit + 0.1), alpha = 1, fill = "white")
  
  
  if (par.l$plotRNASeqClassification) {
    p1 = p1 +  geom_point(data = dataCur, aes(x = radial, y = yValue, size = sign, fill = classification, color = classification), alpha = 1, shape = 21)
  } else {
    p1 = p1 +  geom_point(data = dataCur, aes(x = radial, y = yValue, size = sign), alpha = 1, shape = 21)
  }
  
  p1 = p1 + geom_rect(data = NULL, aes(xmin = 180, xmax = 360 - width_whiteArea , ymin = signThresholdPlot, ymax =  ylimit), alpha = 0.5, fill = "white") 
    
  p1 = p1 + geom_rect(data = NULL, aes(xmin = 0 + width_whiteArea, xmax = 180, ymin = signThresholdPlot, ymax =  ylimit), alpha = 0.5, fill = "white")
    
  
  fontSize = par.l$rootFontSize + nrow(ggrepel_df) * 0.01
  
  # y axis labels
  for (i in 1:(length(labels$breaks))) {
    p1 = p1 + annotate(geom = "text", x = 0, y = labels$breaks[i], label = paste0(as.numeric(labels$breaks[i])), vjust = 0.5, hjust = 0.5, angle = 0, size = fontSize)
  }
  
  p2 = p1

  p2 = p2 + geom_segment(aes(x = 0 + width_whiteArea, xend = 360 - width_whiteArea, y = signThresholdPlot, yend = signThresholdPlot), size = 0.5, linetype = "solid", color = "red", alpha = .5) 
 
  
  
  # Increase the size of the poitns in the legend, see https://stackoverflow.com/questions/20415963/how-to-increase-the-size-of-points-in-legend-of-ggplot2
  #p2 = p2 + guides(fill = guide_legend(override.aes = list(size = par.l$sizeLegend), nrow = 2))
  p2 = p2 + guides(fill = guide_legend(override.aes = list(size = par.l$sizeLegend), nrow = 1))
  
  # Add the annotation outside of the plot (y axis) and draw helper lines
  ## scaling to the angles automatically
  
  coeff_angle = (limit_max * 2) / 360
  
  anglesPos = c(seq(180, 15, -30))
  anglesNeg = c(seq(180 + 30, 345, 30))
  
  df.axis = data.frame(angles = c(anglesPos, anglesNeg), 
                       annotation = c(coeff_angle  * abs(anglesPos - 180),
                                      -coeff_angle * abs(anglesNeg - 180)))
  
  
  p2 = p2 + geom_segment(aes(x = df.axis$angles, xend = df.axis$angles , y = min(labels$breaks), yend = ylimit), 
                         size = par.l$sizeHelperLines, linetype = "dotted", color = par.l$plot_grayColor, alpha = .9) 
  
  df.axis$annotation = signif(df.axis$annotation, 3)
  
  
  for (i in anglesPos) {
    
    p2 = p2 + annotate(geom = "text", x = i, y = ylimit + 0.55, label = paste0(df.axis[which(df.axis$angles == i),]$annotation), vjust = 0, hjust = 0, size = fontSize)
  }
  for (i in anglesNeg) {
    
    p2 = p2 + annotate(geom = "text", x = i, y = ylimit + 0.55, label = paste0(df.axis[which(df.axis$angles == i),]$annotation), vjust = 0, hjust = 1, size = fontSize)
  }
  
  
  
  p3 = p2 
  
  labelY = paste0("Sign. threshold: ", signif(significanceThresholdCur*100,2), "%\n", transform_yValues_caption(), " (adj. p-value)")


  
  #p3 = p3 + annotate(geom = "text", x = 0, y = 0, label = "log10 T stat.", vjust = 0.5, hjust = 0.1, angle = 90, size = fontSize, fontface = 'italic') + 
  p3 = p3 + annotate(geom = "text", x = 0, y =  ylimit + 0.5, label = labelY, vjust = 0, hjust = 0.5, angle = 0, size = fontSize, fontface = 'italic')
  
  angleLegAct = 180
  length_arrow = 80
  yPos = ylimit * 1.4
  yPosLabels = yPos - 0.1
  
  
  
  # Draw the legend in the lower part (arrows)
  p3 = p3 + 
    annotate(geom = "text", x = angleLegAct + 0, y = yPos, label = " TF activity ", vjust = 0, angle = 0, size = fontSize) + 
    # First right side for positive values, corresponding to positive weighted mean difference values (first element of DESeq comparison)
    annotate(geom = "text", x = angleLegAct - 30, y =  yPosLabels, label = conditionComparison[1], vjust = 0, angle = 30, size = fontSize) + 
    # Now left side
    annotate(geom = "text", x = (angleLegAct + 30), y = yPosLabels, label = conditionComparison[2], vjust = 0, angle = -30, size = fontSize) + 
    geom_segment(aes(x = angleLegAct + 10, xend = angleLegAct + 10 + length_arrow , y = yPos, yend = yPos), size = 0.3, arrow = arrow(length = unit(0.6,"cm")))  +
    geom_segment(aes(x = angleLegAct - 10, xend = angleLegAct - 10 - length_arrow , y = yPos, yend = yPos), size = 0.3, arrow = arrow(length = unit(0.6,"cm"))) 
  
  # Add minor y circular lines 
  for (x in 1:length(labels$breaks)) {
    p3 = p3 + geom_hline(yintercept = labels$breaks[x], size = par.l$sizeHelperLines, linetype = "dotted",  color = par.l$plot_grayColor, alpha = .9) 
  }
  
  
  # ggrepel function
  
  if (par.l$plotRNASeqClassification) {
    p3 = p3 + scale_color_manual(values = colorCategoriesCur, guide = FALSE) + 
      geom_label_repel(data = ggrepel_df, aes(x = radial, y = yValue, label = TF, fill = classification),
                       size = par.l$size_TFAnnotation, fontface = 'bold', color = 'white',
                       segment.size = 0.15,
                       label.padding = unit(0.2, "lines"), # how thick is connectin line
                       nudge_y = 0.15, nudge_x = 0,  # how far from center points
                       segment.alpha = .8, segment.color = par.l$plot_grayColor, show.legend = FALSE) + 
      scale_fill_manual(values = colorCategoriesCur, name = "TF class")
  } else {
    
    p3 = p3 + geom_label_repel(data = ggrepel_df, aes(x = radial,
                                                      y = yValue,
                                                      label = TF),
                               fill = par.l$plot_grayColor, size = par.l$size_TFAnnotation, fontface = 'bold', color = 'white',
                               segment.size = 0.15, label.padding = unit(0.1, "lines"), nudge_y = 0.15, nudge_x = 0, 
                               segment.alpha = .8, segment.color = par.l$plot_grayColor, show.legend = FALSE)
  }
  
  
  
  p3 = p3 + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(), axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(), panel.background = element_blank(), panel.grid = element_blank(), legend.position = "top", legend.text = element_text(size = par.l$sizeLegend), legend.title = element_text(size = par.l$sizeLegend, face = "bold"),
                  legend.margin = margin(t = 0, unit = "cm"), legend.key = element_blank(), legend.justification = "center",
                  #plot.margin = grid::unit(c(0, 0, 0, 0), "mm")) + 
                  #plot.margin = unit(c(0,-7,-4,-7), units = "cm")) + 
                  plot.margin = unit(c(0,0,0,0), units = "cm")) + xlab("") + ylab("") 
  
  p3
  
} # end function plotCircular


heatmap.act.rep <- function(df.tf.peak.matrix, tf2ensg.exp){
  
  
  cor.r.pearson.m <- cor.m[,tf2ensg.exp$ENSEMBL]
  
  identical(colnames(df.tf.peak.matrix), tf2ensg.exp$HOCOID)
  identical(colnames(cor.r.pearson.m), tf2ensg.exp$ENSEMBL)
  colnames(cor.r.pearson.m) <- tf2ensg.exp$HOCOID
  BREAKS = seq(-1,1,0.05)
  diffDensityMat = matrix(NA, nrow = ncol(cor.r.pearson.m), ncol = length(BREAKS)-1)
  rownames(diffDensityMat) = tf2ensg.exp$HOCOID
  
  TF_Peak_all.m <- df.tf.peak.matrix
  TF_Peak.m <- TF_Peak_all.m
  
  for (i in 1:ncol(cor.r.pearson.m)){
    TF = colnames(cor.r.pearson.m)[i]
    TF_name = TF #as.character(tf2ensg.exp$HOCOID[tf2ensg.exp$HOCOID==TF])
    ## for the background, use all peaks
    h_noMotif = hist(cor.r.pearson.m[,TF][TF_Peak_all.m[,TF]==0], breaks = BREAKS, plot = FALSE)
    ## for the foreground use only peaks with less than min_mot_n different TF motifs
    h_Motif = hist(cor.r.pearson.m[,TF][TF_Peak.m[,TF]!=0], breaks = BREAKS, plot = FALSE)
    diff_density = h_Motif$density - h_noMotif$density
    diffDensityMat[rownames(diffDensityMat)==TF_name[1], ] <- diff_density
  }
  diffDensityMat = diffDensityMat[!is.na(diffDensityMat[,1]),]
  colnames(diffDensityMat) = signif(h_Motif$mids,1)
  quantile(diffDensityMat)
  
  ## check to what extent the number of TF motifs affects the density values
  n_min = ifelse(colSums(TF_Peak.m)<nrow(TF_Peak.m),colSums(TF_Peak.m), nrow(TF_Peak.m)-colSums(TF_Peak.m))
  names(n_min) = tf2ensg.exp$HOCOID#[match(names(n_min), as.character(tf2ensg$ENSEMBL))]
  n_min <- sapply(split(n_min,names(n_min)),sum)
  quantile(n_min)
  remove_smallN = which(n_min<100)
  cor(n_min[-remove_smallN],rowMax(diffDensityMat)[-remove_smallN], method = 'pearson')
  
  factorClassificationPlot <- sort(median.cor.tfs, decreasing = TRUE)
  diffDensityMat_Plot = diffDensityMat[match(names(factorClassificationPlot), rownames(diffDensityMat)), ]
  diffDensityMat_Plot = diffDensityMat_Plot[!is.na(rownames(diffDensityMat_Plot)),]
  annotation_rowDF = data.frame(median_diff = 
                                  factorClassificationPlot[match(rownames(diffDensityMat_Plot), names(factorClassificationPlot))])
  colBreaks = unique(c((-1),round(quantile(median.cor.tfs.non, probs = c(.05, .95) )[1], digits = 3), round(quantile(median.cor.tfs.non, probs = c(.05, .95) )[2], digits=3),1))
  anno_rowDF = data.frame(
    threshold = cut(annotation_rowDF$median_diff, 
                    breaks = colBreaks))
  rownames(anno_rowDF) = rownames(diffDensityMat_Plot)
  colors = c(par.l$colorCategories["repressor"],par.l$colorCategories["not-expressed"], par.l$colorCategories["activator"])
  names(colors) = levels(anno_rowDF$threshold)
  
  
  pheatmap(diffDensityMat_Plot, cluster_rows = FALSE, cluster_cols = FALSE,
           fontsize_row = 1.25, scale = 'row' , fontsize_col = 10, fontsize=8, labels_col=c(-1, -0.5, 0, 0.5, 1),
           annotation_row = anno_rowDF,annotation_legend=F,
           annotation_colors = list(threshold = colors), legend = T, annotation_names_row = FALSE)
  
} # end function


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
output.global.TFs.orig$weighted_meanDifference = as.numeric(output.global.TFs.orig$weighted_meanDifference)
output.global.TFs.orig$variance = as.numeric(output.global.TFs.orig$variance)
output.global.TFs.orig$weighted_CD = as.numeric(output.global.TFs.orig$weighted_CD)
output.global.TFs.orig$weighted_Tstat = as.numeric(output.global.TFs.orig$weighted_Tstat)



# Remove rows with NA
TF_NA = which(is.na(output.global.TFs.orig$weighted_meanDifference))

if (length(TF_NA) > 0) {
    output.global.TFs.orig = output.global.TFs.orig[-TF_NA,]
  
  TFs_NA = output.global.TFs.orig$TF[TF_NA]
  message = paste0("The following TF have been removed from the data due to NA values in weighted_meanDifference (insufficient data in previous steps): ", paste0(unique(TFs_NA), collapse = ", "))
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
}


######################################
# FILTER BY PERMUTATIONS AND COMPARE, DIAGNOSTIC PLOTS #
######################################
# Compare the distributions from the real and random permutations
diagPlots.l = list()


pdf(par.l$files_plotDiagnostic[1])
output.global.TFs.orig$pvalue = 0
if (par.l$nPermutations > 0) {

    dataReal.df = filter(output.global.TFs.orig, permutation == 0)
    dataPerm.df = filter(output.global.TFs.orig, permutation > 0)

    # ( (#perm >TH)/#perm ) / ( (#perm >TH)/#perm + (#real > TH)/n_real )
    
    xrange = range(output.global.TFs.orig$weighted_meanDifference, na.rm = TRUE)
    
    plot(density(dataPerm.df$weighted_meanDifference, na.rm = TRUE), col = "black", main = "Weighted mean difference values (black = permuted)", xlim = xrange)
    lines(density(dataReal.df$weighted_meanDifference, na.rm = TRUE), col = "red")
    

  
    
    for (TFCur in unique(output.global.TFs.orig$TF)) {
        
        dataCur.df = filter(output.global.TFs.orig, TF == TFCur)
        dataReal.df = filter(dataCur.df, permutation == 0)
        dataPerm.df = filter(dataCur.df, permutation > 0)
        
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
  # 
  # 
  # output.global.TFs.orig$weighted_Tstat_centralized = output.global.TFs.orig$weighted_Tstat - MLE.delta
  # 
  
  
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


#output.global.TFs.orig$percentileAbs =  pmin(output.global.TFs.orig$percentile, 1 - output.global.TFs.orig$percentile)
output.global.TFs.permutations = filter(output.global.TFs.orig, permutation > 0)
output.global.TFs              = filter(output.global.TFs.orig, permutation == 0)

output.global.TFs$Cohend_factor = ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[1], par.l$classes_CohensD[1], 
                                         ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[2] , par.l$classes_CohensD[2], 
                                                ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[3], par.l$classes_CohensD[3], par.l$classes_CohensD[4])))

output.global.TFs$Cohend_factor = factor(output.global.TFs$Cohend_factor, levels = par.l$classes_CohensD, labels = seq_len(length(par.l$classes_CohensD)))


output.global.TFs$pvalueAdj = p.adjust(output.global.TFs$pvalue, method = "BH")


colnamesToPlot = c("weighted_meanDifference", "weighted_CD", "TFBS", "weighted_Tstat", "variance", "pvalue", "pvalueAdj")


for (pValueCur in c(par.l$significanceThresholds , 1)) {
  
  filtered.df = filter(output.global.TFs, pvalueAdj <= pValueCur)
  
  title = paste0("p-value: ", pValueCur, " (retaining ", nrow(filtered.df), " TF)")
  
  for (measureCur in colnamesToPlot) {
      
      if (all(!is.finite(unlist(filtered.df[,measureCur])))) {
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

    
    sampleSummary.df = read_tsv(file_sampleSummary, col_types = cols())
    
    
    # Loading TF data
    HOCOMOCO_mapping.df = read.table(file = par.l$file_input_HOCOMOCO_mapping, header = TRUE)
    assertSubset(c("ENSEMBL", "HOCOID"), colnames(HOCOMOCO_mapping.df))
    
    # Read RNAseq counts
    # par.l$file_input_geneCountsPerSample = "/g/scb2/zaugg/berest/Projects/CLL/PREPARE/Armando.AR/52samplesAR/RNA/rawCounts_8samples.tsv"
    
    TF.counts.df.all = read_tsv(par.l$file_input_geneCountsPerSample, col_names = TRUE)
    if (nrow(problems(TF.counts.df.all)) > 0) {
      flog.fatal(paste0("Parsing errors: "), problems(TF.counts.df.all), capture = TRUE)
      stop("Error when parsing the file ", par.l$file_input_geneCountsPerSample, ", see errors above")
    }
    
    # Add row names
    TF.counts.df.all = as.data.frame(TF.counts.df.all)
    rownames(TF.counts.df.all) = TF.counts.df.all$ENSEMBL
    TF.counts.df.all = TF.counts.df.all[,-1]

    
    sampleData.l = readRDS(par.l$file_input_metadata)
    sampleData.df = sampleData.l[["permutation0"]]
    sampleData.df = filter(sampleData.df, SampleID %in% colnames(TF.counts.df.all))
    
    nFiltRows = nrow(sampleData.l[["permutation0"]]) - nrow(sampleData.df)
    if (nFiltRows > 0) {
      flog.warn(paste0("Filtered ", nFiltRows, " sample IDs afteer comparising sample names with RNA-Seq table"))
    }
    
    par.l$designFormula = snakemake@config$par_general$designContrast
    designFormula = convertToFormula(par.l$designFormula, colnames(sampleData.df))
 
    ####################################
    # Run DeSeq2 on raw RNA-Seq counts #
    ####################################
    dd <- DESeqDataSetFromMatrix(countData = TF.counts.df.all,
                                 colData = sampleData.df,
                                 design = designFormula)
    
    dd = estimateSizeFactors(dd)
    dd = DESeq(dd)
    dd_counts =  counts(dd, normalized=TRUE)
    TF.counts.df.all = dd_counts %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% as.tibble()
    
    # Check sample names and set column names
    # Match the column names and do the intersections
    sharedColumns = intersect(colnames(TF.counts.df.all)[-1], sampleSummary.df$SampleID)
    
    if (length(sharedColumns) == 0) {
      message = paste0("No shared samples with RNA-Seq samples between sample table ", file_sampleSummary, " and RNA-Seq table ", par.l$file_input_geneCountsPerSample, ".")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    } else {
      flog.info(paste0(length(sharedColumns), " out of ", ncol(TF.counts.df.all) - 1, " columns are shared between RNA-Seq counts and the sample table"))
    }
    
    colnames(TF.counts.df.all)[1] = "ENSEMBL"
    
    # Clean ENSEMBL IDs
    TF.counts.df.all$ENSEMBL = gsub("\\..+", "", TF.counts.df.all$ENSEMBL, perl = TRUE)
    

    
    # Filter them by the IDs that correspond to the TFs
    TF.counts.df = filter(TF.counts.df.all, ENSEMBL %in% HOCOMOCO_mapping.df$ENSEMBL)
    
    if (nrow(TF.counts.df) == 0) {
      message = "No rows remaining after filtering against ENSEMBL IDs in HOCOMOCO. Check your ENSEMBL IDs for overlap with the HOCOMOCO translation table."
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    
    # Read counts per TFBS, coming from DESeqPeaks.R
    # TODO: Provide counts matrix
    peak.counts = read_tsv(par.l$file_input_countsNorm, col_types = cols())
    if (nrow(problems(peak.counts)) > 0) {
        flog.fatal(paste0("Parsing errors: "), problems(peak.counts), capture = TRUE)
        stop("Error when parsing the file ", par.l$file_input_countsNorm, ", see errors above")
    }
    
    # Match the column names and do the intersections
    sharedColumns = intersect(colnames(peak.counts), colnames(TF.counts.df))
    
    if (length(sharedColumns) == 0) {
      message = "No shared samples with RNA-Seq samples."
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    flog.info(paste0(length(sharedColumns), " samples are shared between the input data and RNA-Seq data"))
    
    #peak.counts.orig = peak.counts
    peak.counts  = dplyr::select(peak.counts, one_of("peakID", sharedColumns)) 
    TF.counts.df = TF.counts.df[, which(colnames(TF.counts.df) %in% c(sharedColumns, "ENSEMBL"))]
    
    
    HOCOMOCO_mapping.df.exp <- filter(HOCOMOCO_mapping.df, ENSEMBL %in%  TF.counts.df$ENSEMBL, HOCOID %in% output.global.TFs$TF)
    
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
      
      TF.output.df = read.table(file = paste0(rootOutdir, "/TF-SPECIFIC/",TFCur,"/extension", extensionSize, "/", comparisonType, TFCur,  ".output.tsv.gz"), header = TRUE)
      TF.peakMatrix.l[[TFCur]] = peak.counts$peakID %in% TF.output.df$peakID
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
    
    index_lastColumn = which(colnames(peak.counts) == "peakID")
    
    if (!all(colnames(expressed.TF.counts.df) == colnames(peak.counts)[-index_lastColumn])) {
      message = "Colnames of expressed.TF.counts.df and peak.counts must be identical"
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    # Filter by rowMeans to eliminate rows with an sd of 0
    rowMeans1 = rowMeans(expressed.TF.counts.df)
    rowsToDelete = which(rowMeans1 == 0)
    rowsToDelete = which(rowMeans1 < 1)
    if (length(rowsToDelete) > 0) {
      flog.info(paste0("Removed ", length(rowsToDelete), " TFs out of ", nrow(expressed.TF.counts.df), " because they had a row mean of 0."))
       expressed.TF.counts.df = expressed.TF.counts.df[-rowsToDelete,]
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
    median.cor.tfs.non <- sort(apply(t.cor.sel.matrix.non, MARGIN=2, FUN = my.median))
    
    
    median.cor.tfs.rest <- sort(median.cor.tfs - median.cor.tfs.non[names(median.cor.tfs)])
    
    act.rep.thres = quantile(sort(apply(t.cor.sel.matrix.non, MARGIN = 2, FUN = my.median)), probs = c(.05, .95))

    
    AR.data = as.data.frame(median.cor.tfs)
    AR.data$TF = rownames(AR.data)
    
    output.global.TFs = merge(output.global.TFs, AR.data, by = "TF",all.x = TRUE)
    
    output.global.TFs$classification = ifelse(is.na(output.global.TFs$median.cor.tfs), "not-expressed",
                                            ifelse(output.global.TFs$median.cor.tfs <= act.rep.thres[1], "repressor",
                                                   ifelse(output.global.TFs$median.cor.tfs > act.rep.thres[2], "activator", "undetermined")))
    
    output.global.TFs$classification = factor(output.global.TFs$classification, levels = names(par.l$colorCategories))
    
    ####################
    ####################
    # DIAGNOSTIC PLOTS #
    ####################
    ####################
    
    
    pdf(file = par.l$files_plotDiagnostic[2], width = 3, height = 8)
    xlab="median pearson correlation (r)"
    ylab=""
    xlim= c(-0.2,0.2)
    ylim=c(1,length(median.cor.tfs.non))
    
    par(mfrow=c(1,1))
    
    plot(median.cor.tfs.non[names(median.cor.tfs)], 1:length(median.cor.tfs.non),
         xlim=xlim,
         ylim=ylim,
         main="",
         xlab=xlab,
         ylab=ylab,
         col=adjustcolor("darkgrey",alpha=1),
         pch=16,
         cex=0.5,
         axes = FALSE)
    points(median.cor.tfs, 1:length(median.cor.tfs),
           pch=16,
           col=ifelse(median.cor.tfs>act.rep.thres[2], par.l$colorCategories["activator"] ,ifelse(median.cor.tfs<act.rep.thres[1], par.l$colorCategories["repressor"], par.l$colorCategories["undetermined"])), 
           cex=0.5
    ) 
    text(x =c((act.rep.thres[1]-0.01),(act.rep.thres[2]+0.01)),
         y=c((length(median.cor.tfs.non)+5),(length(median.cor.tfs.non)+5)), pos=c(2,4),
         labels =c("5th\nPercentile","95th\nPercentile"),cex=0.7, col=c("black","black"))
    abline(v=act.rep.thres[1], col=par.l$colorCategories["repressor"])
    abline(v=act.rep.thres[2], col=par.l$colorCategories["activator"])
    axis(side = 1, lwd = 1, line = 0, at = c(-0.2,0,0.2), cex=1)
    

    heatmap.act.rep(TF.peakMatrix.df, HOCOMOCO_mapping.df.exp)
    dev.off()
    
    
    # TODO: TF_names_update = read_delim("/g/scb2/zaugg/berest/Projects/CLL/PREPARE/Armando.AR/52samplesAR/HOCOTFID2ENSEMBL.txt",delim = " ")
    
    # Filter genes
    samples_cond1 = colData(dd)$SampleID[which(colData(dd)$conditionSummary == levels(colData(dd)$conditionSummary)[1])]
    samples_cond2 = colData(dd)$SampleID[which(colData(dd)$conditionSummary == levels(colData(dd)$conditionSummary)[2])]
    
    nMin = par.l$filter_minCountsPerCondition
    idx <- (rowMeans(dd_counts[,samples_cond1]) > nMin | rowMeans(dd_counts[,samples_cond2]) > nMin) & rowMedians(dd_counts) > 0
    dd.filt = dd[idx,]
    dd.filt <- DESeq(dd.filt)
    
    # Process results
    res.peaks.filt     = results(dd.filt) %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% as.tibble()
    expression.df.filt = counts(dd.filt , normalized=TRUE) %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% as.tibble()
    expresssion.df.all = full_join(expression.df.filt, res.peaks.filt, by = 'ENSEMBL')
    
    # expresssion.df.all = plyr::join_all(list(expression.df.filt,res.peaks.filt), by = 'ENSEMBL', type = 'inner')
    
    # TODO
    # In ivans version, he used a non-filtered HOCOMOCO table. I however filter it before, so that some ENSEMBL IDs might already be filtered
    
    TF.specific = left_join(HOCOMOCO_mapping.subset.df, res.peaks.filt, by = "ENSEMBL") %>% filter(!is.na(baseMean))
    
    
    
    # TODO:_ deal with NAs, where do they come from?
    

    
    output.global.TFs$weighted_meanDifference = as.numeric(output.global.TFs$weighted_meanDifference)
    
    par.l$minPointSize = 0.3
    
    output.global.TFs.merged = output.global.TFs %>%
      filter(classification != "not-expressed")  %>%
      full_join(TF.specific, by = c( "TF" = "HOCOID"))  %>%
      # TODO check and implemenmt differently
      mutate(log2FoldChange = log2FoldChange * -1) %>%
      mutate(baseMeanNorm = (baseMean - min(baseMean, na.rm = TRUE)) / (max(baseMean, na.rm = TRUE) - min(baseMean, na.rm = TRUE)) + par.l$minPointSize)  %>%
      filter(!is.na(classification)) 
    

    #######################################
    # Correlation plots for the 3 classes #
    #######################################
    
    pdf(par.l$files_plotDiagnostic[3])
    for (classificationCur in unique(output.global.TFs.merged$classification)) {
      
      output.global.TFs.cur = filter(output.global.TFs.merged, classification == classificationCur)
      
      cor.res.l = list()
      for (corMethodCur in c("pearson", "spearman")) {
        cor.res.l[[corMethodCur]] = cor.test(output.global.TFs.cur$weighted_meanDifference, output.global.TFs.cur$log2FoldChange, method = corMethodCur)
      }
      
      titleCur = paste0(classificationCur, ": R=", 
                        signif(cor.res.l[["pearson"]]$estimate, 2), "/", 
                        signif(cor.res.l[["spearman"]]$estimate, 2), ", p-value ", 
                        signif(cor.res.l[["pearson"]]$p.value,2),  "/", 
                        signif(cor.res.l[["spearman"]]$p.value,2), "\n(Pearson/Spearman)")
      
      g = ggplot(output.global.TFs.cur, aes(weighted_meanDifference, log2FoldChange)) + geom_point(aes(size = baseMeanNorm)) + 
        geom_smooth(method = par.l$regressionMethod, color = par.l$colorCategories[classificationCur]) + 
        ggtitle(titleCur) + 
        ylab("log2 fold-change RNA-seq") + 
        theme_bw() + theme(plot.title = element_text(hjust = 0.5))
      plot(g)
      
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
      legend("topleft",box.col = adjustcolor("white",alpha.f = 0),
             legend = c("Motif","Non-motif"),
             lwd=c(2,2),cex = 0.8,
             col=c("red","darkgrey"), lty=c(1,1) )
      axis(side = 1, lwd = 1, line = 0)
      axis(side = 2, lwd = 1, line = 0, las = 1)
      
      lines(density(dataBackground, bw=0.1, na.rm=T), lwd=2.5, col="darkgrey")
      
    } 
    dev.off()
    
    

} else {
  classesList.l = list(c())
}

output.global.TFs$yValue = transform_yValues(output.global.TFs$pvalueAdj)
output.global.TFs.origReal = output.global.TFs

#########################################
# PLOT FOR DIFFERENT P VALUE THRESHOLDS #
#########################################

# TODO: change pvalue to pvalueAdj ?
# Set the page dimensions to the maximum across all plotted variants
output.global.TFs.filteredSummary = filter(output.global.TFs, pvalue <= max(par.l$significanceThresholds))

nTF_label = min(par.l$maxTFsToDraw, nrow(output.global.TFs.filteredSummary))

TFLabelSize = ifelse(nTF_label < 20, 8,
                     ifelse(nTF_label < 40, 7,
                      ifelse(nTF_label < 60, 6,
                        ifelse(nTF_label < 60, 5,
                          ifelse(nTF_label < 60, 4, 3)))))





allPlots.l = list("circular" = list(), "volcano" = list())



variableXAxis = "weighted_meanDifference"



for (significanceThresholdCur in par.l$significanceThresholds) {
  
  pValThrStr = as.character(significanceThresholdCur)

  for (showClasses in classesList.l) {

    #################
    # CIRCULAR PLOT #
    #################
    plotFinal = plotCircular(output.global.TFs.origReal, par.l, showClasses, significanceThresholdCur, conditionComparison, variableXAxis)
    #plotFinal = NA
    if (par.l$plotRNASeqClassification) {
      allPlots.l[["circular"]] [[pValThrStr]] [[paste0(showClasses,collapse = "-")]] = plotFinal
    } else {
      allPlots.l[["circular"]] [[pValThrStr]]  = plotFinal
    }
    
    ################
    # VOLCANO PLOT #
    ################
    
    output.global.TFs = output.global.TFs.origReal %>%
        mutate( pValueAdj_log10 = transform_yValues(pvalueAdj),
                pValue_log10 = transform_yValues(pvalue),
                pValueAdj_sig = pvalueAdj <= significanceThresholdCur,
                pValue_sig = pvalue <= significanceThresholdCur) 
       
    if (par.l$plotRNASeqClassification) {
        output.global.TFs = filter(output.global.TFs, classification %in% showClasses)
    }
    
    for (pValueStrCur in c("pvalue", "pvalueAdj")) {
        
        
        if (pValueStrCur == "pvalue") {
            
            pValueScoreCur = "pValue_log10"
            pValueSigCur = "pValue_sig"
            pValueStrLabel = "raw p-value"
            
            ggrepel_df = filter(output.global.TFs, pValue_sig == TRUE)
            maxPValue = max(output.global.TFs$pValue_log10, na.rm = TRUE)
            
        } else {
            
            pValueScoreCur = "pValueAdj_log10"
            pValueSigCur = "pValueAdj_sig"
            pValueStrLabel = "adj. p-value"
            
            ggrepel_df = filter(output.global.TFs, pValueAdj_sig == TRUE)
            maxPValue = max(output.global.TFs$pValueAdj_log10, na.rm = TRUE)
        }
    
        
        # Increase the ymax a bit more
        ymax = max(transform_yValues(significanceThresholdCur), maxPValue, na.rm = TRUE) * 1.1
        alphaValueNonSign = 0.3
        
            
        # Reverse here because negative values at left mean that the condition that has been specified in the beginning is higher. 
        # Reverse the rev() that was done before for this plot therefore to restore the original order
        labelsConditionsNew = rev(conditionComparison)
    
        g = ggplot()
        
        if (par.l$plotRNASeqClassification) {
          g = g + geom_point(data = output.global.TFs, aes_string("weighted_meanDifference", pValueScoreCur, alpha = pValueSigCur, size = "TFBS", fill = "classification"), shape=21, stroke = 0.5, color = "black") +  scale_fill_manual("TF class", values = par.l$colorCategories)
          
          g = g +
              geom_rect(aes(xmin = -Inf,xmax = 0,ymin = -Inf, ymax = Inf, color = par.l$colorConditions[2]),
                        alpha = .3, fill = par.l$colorConditions[2], size = 0) +
              geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf,ymax = Inf, color = par.l$colorConditions[1]),                                                                                   alpha = .3, fill = par.l$colorConditions[1], size = 0) + 
              scale_color_manual(name = 'TF activity higher in', values = par.l$colorConditions, labels = conditionComparison)
          
        } else {
            
          g = g + geom_point(data = output.global.TFs, aes_string("weighted_meanDifference", pValueScoreCur, alpha = pValueSigCur, size = "TFBS"), shape=21, stroke = 0.5, color = "black")
          g = g + geom_rect(aes(xmin = -Inf,
                                xmax = 0,
                                ymin = -Inf, 
                                ymax = Inf, fill = par.l$colorConditions[2]),
                            alpha = .3) + 
              geom_rect(aes(xmin = 0,
                            xmax = Inf,
                            ymin = -Inf, 
                            ymax = Inf, fill = par.l$colorConditions[1]),
                        alpha = .3)
          g = g + scale_fill_manual(name = 'TF activity higher in', values = rev(par.l$colorConditions), labels = labelsConditionsNew)
        }
     
    
        g = g + ylim(-0.1,ymax) + 
            ylab(paste0(transform_yValues_caption(), " (", pValueStrLabel, ")")) + 
            xlab("weighted mean difference") + 
            scale_alpha_manual(paste0(pValueStrLabel, " < ", significanceThresholdCur), values = c(alphaValueNonSign, 1), labels = c("no", "yes")) + 
            geom_hline(yintercept = transform_yValues(significanceThresholdCur), linetype = "dotted") 
        
        if (nrow(ggrepel_df) <= par.l$maxTFsToDraw) {
            
            if (par.l$plotRNASeqClassification) {
                g = g +  geom_label_repel(data = ggrepel_df, aes_string("weighted_meanDifference", pValueScoreCur, label = "TF", fill = "classification"),
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
            
            flog.warn(paste0("Not labeling significant TFs, maximum of ", par.l$maxTFsToDraw, " exceeded for ", pValThrStr, " and ", pValueStrCur))
            
            
            if (nrow(ggrepel_df) > par.l$maxTFsToDraw) {
                
                labelPlot = paste0("*TF labeling skipped because number of significant TFs\nexceeds the maximum of ", par.l$maxTFsToDraw, " (", nrow(ggrepel_df), ")")
                flog.warn(labelPlot)
                
                g = g + annotate("text", label = labelPlot, x = 0, y = ymax, size = 3)
                
            }
            
            
            
        }
        
      
          g = g + theme_bw() + 
              theme(axis.text.x = element_text(size=rel(1.5)),
                    axis.text.y = element_text(size=rel(1.5)), 
                    axis.title.x = element_text(size=rel(1.5)),
                    axis.title.y = element_text(size=rel(1.5)),
                    legend.title=element_text(size=rel(1.5)), 
                    legend.text=element_text(size=rel(1.5))) 
          
          if (par.l$plotRNASeqClassification) {
            g = g + guides(alpha = guide_legend(override.aes = list(size=5), order = 2),
                           fill = guide_legend(override.aes = list(size=5), order = 3),
                           color = guide_legend(override.aes = list(size=5), order = 1))
            
            allPlots.l[["volcano"]] [[pValThrStr]] [[paste0(showClasses,collapse = "-")]] [[pValueStrCur]] = g
          } else {
            g = g + guides(alpha = guide_legend(override.aes = list(size=5), order = 2),
                           fill = guide_legend(override.aes = list(size=5), order = 3))
     
            
            allPlots.l[["volcano"]] [[pValThrStr]] [[pValueStrCur]] = g
          }
          
      } # end separately for raw and adjusted p-values
        
    } # end for all showClasses



} # end for different significance thresholds


#####################
# CIRCULAR PLOT PDF #
#####################
height = width = max(nTF_label / 5, par.l$circularPlot_minDimensions)

# Increase a bit towards smaller heights by a factor, empirical observation
if (height < 20) {
  height = width = height * 1.2
}
pdf(file = par.l$file_plotCircular, height = height, width = width, useDingbats = FALSE)
for (significanceThresholdCur in par.l$significanceThresholds) {
  
  for (showClasses in classesList.l) {
    if (par.l$plotRNASeqClassification) {
      plot(allPlots.l[["circular"]] [[as.character(significanceThresholdCur)]] [[paste0(showClasses,collapse = "-")]])
    } else {
      plot(allPlots.l[["circular"]] [[as.character(significanceThresholdCur)]])
    }
    
  }
}
dev.off()


####################
# VOLCANO PLOT PDF #
####################
height = width = max(nTF_label / 15 , par.l$circularPlot_minDimensions)
pdf(file = par.l$file_plotVolcano, height = height, width = width, useDingbats = FALSE)


for (pValueStrCur in c("pvalueAdj", "pvalue")) {
    
    for (significanceThresholdCur in par.l$significanceThresholds) {
      
      for (showClasses in classesList.l) {

              if (par.l$plotRNASeqClassification) {
                  plot(allPlots.l[["volcano"]] [[as.character(significanceThresholdCur)]] [[paste0(showClasses,collapse = "-")]] [[pValueStrCur]])
              } else {
                  plot(allPlots.l[["volcano"]] [[as.character(significanceThresholdCur)]] [[pValueStrCur]])
              }
          
          
        
        
      }
    }
}
dev.off()


output.global.TFs.origReal = dplyr::select(output.global.TFs.origReal, -one_of("permutation", "yValue"))
output.global.TFs.origReal.transf = dplyr::mutate_if(output.global.TFs.origReal, is.numeric, as.character)
write_tsv(output.global.TFs.origReal.transf, path = par.l$file_output_summary, col_names = TRUE)
saveRDS(allPlots.l, file = par.l$file_output_plots)

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)
