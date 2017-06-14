## README: this script is to remove CG bias by comparing permutations from the same CG bin 

start.time  <-  Sys.time()

#########################
# LIBRARY AND FUNCTIONS #
#########################

source("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/functions.R")

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "lsr", "dplyr", "ggrepel", "checkmate", "tools"), verbose = TRUE)



###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$FDR_threshold = c(0.01, 0.05,0.1,0.2, 0.5, 0.7)
par.l$volcanoPlot_height = 16
par.l$volcanoPlot_width  = 24
par.l$volcanoPlot_dpi  = 600
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.2, 0.5, 0.8)
par.l$excludeSmallCohensPrinting = FALSE
par.l$log_minlevel = "INFO"
par.l$fixedMinTStatValue  = 1.96
par.l$circularPlot_minDimensions  = 12

par.l$extension_x_limits = 0.025 # TODO
par.l$extension_x_limits = 0.2 # 20 % x axis extension, regardless of the limits
par.l$extension_y_limits = 1
par.l$color_points = c("darkgreen","tomato1","darkred","slategrey")
par.l$circularPlot_height = 12
par.l$circularPlot_width = 16
par.l$legend_position = c(0.1, 0.9)


############################################
# READ AND VALIDATE COMMAND LINE ARGUMENTS #
############################################

args <- commandArgs(trailingOnly = TRUE)

args = c(
  condComp \
  "/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/CLL.all.volcano.l2FC.removedCGbias.AR.csv",
  "/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/CLL.all.volcano.l2FC.removedCGbias.AR.pdf",
  ,
  TRUE
  
  
)

par.l$file_plotCircular = "/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/CLL.all.circular.l2FC.removedCGbias.AR.pdf"
output.global.TFs = read_tsv("/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/CLL.all.volcano.l2FC.removedCGbias.AR.FDR5.tsv")
  


if (length(args) != 7) {
  stop("Expecting 7 arguments but found ", length(args),". Exiting.")
} else {
  par.l$files_input_permResults   = args[1]
  par.l$file_input_condCompDeSeq  = args[2]
  par.l$file_output_summary       = args[3]
  par.l$file_plotVolcano          = args[4]
  par.l$file_plotCircular         = args[5]
  par.l$plotTFClassification      = as.logical(args[6])
  par.l$file_log                  = args[7]
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

# TODO: Right now this is fixed and no check is performed etc if the object is valid. It only works for CLl right now as it is dataset-specific
if (par.l$plotTFClassification) {
  par.l$activactorRepressorObject = "/scratch/carnold/CLL/TF_act_10samplesOnlyRNASeq/output/FINAL_OUTPUT/extension100/tf_activator_repressor_classification.RData"
  assertFileExists(par.l$activactorRepressorObject)
} else {
  par.l$activactorRepressorObject = NULL
}


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

T_stat = c()
for (FDRCur in par.l$FDR_threshold) {
  threshold1 = FDRCur /  nTF
  threshold1 = FDRCur
  T_stat = c(T_stat, qt(threshold1/2, median(output.global.TFs$TFBS,  na.rm = TRUE), lower.tail = FALSE))
  
  output.global.TFs[, paste0("FDR_", FDRCur)] = qt(threshold1/2, output.global.TFs$TFBS, lower.tail = FALSE)
}



min_T_stat = min(T_stat)
min_T_stat = 1.96

if (par.l$excludeSmallCohensPrinting) {
  plot_thr_df = output.global.TFs[output.global.TFs$Cohend_factor != "1", ]
} else {
  plot_thr_df = output.global.TFs
}

plot_thr_df = plot_thr_df[abs(plot_thr_df$weighted_Tstat) > min_T_stat, ]


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
  geom_hline(yintercept = T_stat, size = 0.5,
             linetype = "longdash", color = "red") +
  annotate(geom="text", label= par.l$FDR_threshold, x=Inf, y = T_stat, vjust = 0, color = "red", size = 3, hjust = 1) + 
  
  geom_hline(yintercept = min_T_stat, size = 0.5, linetype = "longdash", color = "darkgreen") + 

  geom_label_repel(aes(x = plot_thr_df$weighted_mean,
                       y = abs(plot_thr_df$weighted_Tstat),
                       size = output.global.TFs$Cohend_factor,
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


ggsave(plot = TF_volcano1, height = par.l$volcanoPlot_height, width = par.l$volcanoPlot_width, dpi = par.l$volcanoPlot_dpi, filename = par.l$file_plotVolcano, useDingbats = FALSE)


########################
#### LETS GET FANCY ####
########################


## circular need 
limit_max = max(abs(output.global.TFs$weighted_mean)) * (1 + par.l$extension_x_limits)
limit_min = 0

radial_coeff = (180*180)/limit_max

output.global.TFs.neg = output.global.TFs[which(output.global.TFs$weighted_mean < 0),]
output.global.TFs.pos = output.global.TFs[which(output.global.TFs$weighted_mean > 0),]
output.global.TFs.neg$radial = -((radial_coeff*abs(output.global.TFs.neg$weighted_mean))/180) + 360
output.global.TFs.pos$radial = (radial_coeff*abs(output.global.TFs.pos$weighted_mean))/180

output.global.TFs.comb = rbind(output.global.TFs.neg,output.global.TFs.pos)


output.global.TFs.comb$sign = as.factor(ifelse(abs(output.global.TFs.comb$weighted_Tstat) > min_T_stat, TRUE, FALSE))


if (par.l$plotTFClassification) {
  
  ## this is AR reporessors part 
  load(par.l$activactorRepressorObject)
  
  stopifnot(!is.null(median.cor.tfs))
  stopifnot(!is.null(act.rep.thres))
  
  AR.data = as.data.frame(median.cor.tfs)
  AR.data$TF = rownames(AR.data)
  output.global.TFs.comb2 = merge(output.global.TFs.comb, AR.data, by = "TF",all.x = T)
  
  output.global.TFs.comb2$AR2 = ifelse(is.na(output.global.TFs.comb2$median.cor.tfs), "not-expressed",
                                       ifelse(output.global.TFs.comb2$median.cor.tfs < act.rep.thres[1], "repressor",
                                              ifelse(output.global.TFs.comb2$median.cor.tfs > act.rep.thres[2], "activator", "undetermined")))
  
  output.global.TFs.comb3 = merge(output.global.TFs.comb,output.global.TFs.comb2, by = "TF" )
  output.global.TFs.comb = output.global.TFs.comb3
  output.global.TFs.comb$AR2 = factor(output.global.TFs.comb$AR2, levels = c("activator","undetermined","repressor","not-expressed"))
  
  ggrepel_df = output.global.TFs.comb[which(output.global.TFs.comb$sign.x == TRUE ),]
  
} else {
  
  output.global.TFs.comb$weighted_Tstat.x = output.global.TFs.comb$weighted_Tstat
  output.global.TFs.comb$radial.x = output.global.TFs.comb$radial
  output.global.TFs.comb$sign.x = output.global.TFs.comb$sign
  
  ggrepel_df = output.global.TFs.comb[which(output.global.TFs.comb$sign == TRUE ),]
  
}

## scaling to the angles automatically

coeff_angle = (limit_max * 2) / 360

anglesPart1 = seq(0,150,30)
anglesPart2 = seq(210,330,30)

df.axis = data.frame(angles = c(anglesPart1,165,195,anglesPart2), 
                     annotation = c(anglesPart1 * coeff_angle,
                                    coeff_angle*175,
                                    -coeff_angle*175,
                                    -(rev(anglesPart1) * coeff_angle)[-length(anglesPart1)]))
df.axis$annotation = round(df.axis$annotation,2)



## ggplot p objects
ylimit = max(log10(abs(output.global.TFs.comb$weighted_Tstat.x))) + par.l$extension_y_limits

fillVar  = NULL
fillVarSign = NULL
if (par.l$plotTFClassification) {
  fillVar = output.global.TFs.comb$AR2
  fillVarSign = output.global.TFs.comb$AR2[which(output.global.TFs.comb$sign.x == TRUE )]
}

p1 = ggplot( ) + 
  geom_rect(data = NULL,
            aes(xmin = 0,
                xmax = 180,
                ymin = -Inf, 
                ymax = ylimit),
            alpha = .6, 
            fill = "pink") + 
  geom_rect(data = NULL,
            aes(xmin = 180,
                xmax = 360,
                ymin = -Inf, 
                ymax = ylimit),
            alpha = .6, 
            fill ="lightblue") +
  geom_rect(data = NULL,
            aes(xmin = -Inf,
                xmax = Inf,
                ymin = -Inf, 
                ymax = log10(min_T_stat)),
            alpha = .3, 
            fill = "white") + 
  geom_hline(yintercept = ylimit,
             size = 0.7,
             linetype = "solid", 
             color = "grey50",
             alpha = .9) +
  geom_rect(data=NULL,aes(xmin = 165, xmax = 195, ymin = -Inf, ymax = ylimit), alpha = 1, fill = "white") 

if (par.l$plotTFClassification) {
  p1 = p1 +  geom_point(aes(x = output.global.TFs.comb$radial.x, 
                            y = log10(abs(output.global.TFs.comb$weighted_Tstat.x)),
                            label = output.global.TFs.comb$TF,
                            size = output.global.TFs.comb$sign.x,
                            fill = fillVar,
                            color = fillVar), 
                        alpha = 1, shape = 21)
} else {
  p1 = p1 +  geom_point(aes(x = output.global.TFs.comb$radial.x, 
                            y = log10(abs(output.global.TFs.comb$weighted_Tstat.x)),
                            label = output.global.TFs.comb$TF,
                            size = output.global.TFs.comb$sign.x), 
                        alpha = 1, shape = 21)
}

p1 = p1 + coord_polar() + 
  # limits of the angular plot 
  scale_x_continuous(limits = c(0,360)) +
  # sizes of the points for significance
  scale_size_manual(values = c(0.5,1), guide = FALSE) 

if (par.l$plotTFClassification) {
  p1 = p1 + scale_color_manual(values = par.l$color_points, guide = FALSE) + 
    geom_label_repel( aes(x = ggrepel_df$radial.x,
                          y = log10(abs(ggrepel_df$weighted_Tstat.x)),
                          label = ggrepel_df$TF,
                          fill = fillVarSign),
                      # size of the label
                      size = 3.5,
                      fontface = 'bold', color = 'white',
                      segment.size = 0.25,
                      # how thick is connectin line
                      label.padding = unit(0.1, "lines"),
                      # how far from center points
                      nudge_y = 0.15,
                      nudge_x = 0, 
                      segment.alpha = .8,
                      segment.color = 'grey50',show.legend = FALSE)
} else {
  
  p1 = p1 + geom_label_repel( aes(x = ggrepel_df$radial.x,
                                  y = log10(abs(ggrepel_df$weighted_Tstat.x)),
                                  label = ggrepel_df$TF),
                              # size of the label
                              fill = 'grey50',
                              size = 3.5,
                              fontface = 'bold', color = 'white',
                              segment.size = 0.25,
                              # how thick is connectin line
                              label.padding = unit(0.1, "lines"),
                              # how far from center points
                              nudge_y = 0.15,
                              nudge_x = 0, 
                              segment.alpha = .8,
                              segment.color = 'grey50',show.legend = FALSE)
}


# ggrepel function

# significance threshold line 
p1 = p1 + geom_hline(yintercept = log10(min_T_stat), 
                     size = 0.7,
                     linetype = "longdash", 
                     color = "red",
                     alpha = .2)

if (par.l$plotTFClassification) {
  p1 = p1 + 
    scale_fill_manual(values = par.l$color_points) + 
    labs(size = "Significance", fill = "Role:") 
}

p1 = p1 + guides(fill = guide_legend( override.aes = list(size = 8)))



# extracting the labels from object in order to add them manually later
labels        = list()
labels$range  = layer_scales(p1)$y$range$range
labels$sum    = sum(abs(labels$range[1]) + abs(labels$range[2]))
labels$breaks = round(seq(labels$range[1],labels$range[2],1),0)



p2 = p1 +  
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = par.l$legend_position,
        legend.text=element_text(size=16),
        legend.justification = "center",
        plot.margin = grid::unit(c(0, 0, 0, 0), "mm")) + 
  xlab("") + ylab("") 


# Set some size-dependent parameters


fontSize = 4.5 + nrow(ggrepel_df)* 0.01


height = width = max(nrow(ggrepel_df) / 5, par.l$circularPlot_minDimensions)

# Increase a bit towards smaller heights by a factor, empirical observation
if (height < 20) {
  height = width = height * 1.2
}

p2 = p2 + geom_segment(aes(x=df.axis$angles, xend = df.axis$angles , y=min(labels$breaks), yend = max(labels$breaks)), 
                       size = 0.7,linetype = "dotted", color = "grey50",alpha = .4) 
for (i in df.axis$angles) {
  p2 = p2 + annotate(geom = "text", x = i, y = ylimit + 0.25, label = paste0(df.axis[which(df.axis$angles == i),]$annotation), vjust = 0.5, size = fontSize)
}


p3 = p2 

# y axis labels
for (i in 2:length(labels$breaks)) {
  p3 = p3 + annotate(geom = "text", x = 167, y = labels$breaks[i], label = paste0(labels$breaks[i]), vjust = 0.5,hjust = 1, angle = 15, size = fontSize) 
}

p3 = p3 + annotate(geom = "text", x = 183, y = 0, label = "log10 of T-statistic", vjust = -1.5, angle = 80, size = fontSize) + 
  geom_segment(aes(x=10, xend = 90 , y=ylimit + 1, yend = ylimit + 1), size=0.3,
               arrow = arrow(length = unit(0.6,"cm"))) + 
  annotate(geom = "text", x = 0, y = ylimit + 1, label = " TF activity ", vjust = 0, angle = 0, size = fontSize) + 
  
  annotate(geom = "text", x = 30, y = ylimit + 1.25, label = strsplit(conditionComparison, " ")[[1]][[1]], 
           vjust = 0, angle = -26, size = fontSize) + 
  annotate(geom = "text", x = 330, y = ylimit + 1.25, label =  strsplit(conditionComparison, " ")[[1]][[3]], 
           vjust = 0, angle = 26, size = fontSize) + 
  geom_segment(aes(x=350, xend = 270 , y=ylimit + 1, yend = ylimit + 1), size=0.3,
               arrow = arrow(length = unit(0.6,"cm"))) 

# Minor y circular lines 
for (x in 1:length(labels$breaks)){
  p3 = p3 + geom_hline(yintercept = labels$breaks[x],
                       size = 1,
                       linetype = "dotted", 
                       color = "grey50",
                       alpha = .4) 
}



ggsave(plot = p3, height = height, width = width, dpi = par.l$volcanoPlot_dpi, filename = par.l$file_plotCircular, useDingbats = FALSE, limitsize = FALSE)

