library(readr)
library(dplyr)
library(lsr)
library(ggplot2)
library(ggrepel)

par.l = list()

par.l$extension_x_limits = 0.025 # TODO
par.l$extension_x_limits = 0.2 # 20 % x axis extension, regardless of the limits
par.l$extension_y_limits = 1
par.l$color_points = c("darkgreen","tomato1","darkred","slategrey")
par.l$volcanoPlot_height = 12
par.l$volcanoPlot_width = 16
par.l$legend_position = c(0.1, 0.9)
par.l$activactorRepressorObject = "/scratch/carnold/CLL/TF_act_10samplesOnlyRNASeq/output/FINAL_OUTPUT/extension100/tf_activator_repressor_classification.RData"
par.l$activactorRepressorObject = NULL


conditionComparison = readRDS("/scratch/carnold/CLL/TF_act/output/TEMP/extension100/conditionComparison.rds")



volc.new = read_tsv("/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/CLL.all.volcano.l2FC.removedCGbias.AR.FDR5.tsv")
#volc.new = permSummaryAll.df

## circular need 
limit_max = max(abs(volc.new$weighted_mean)) * (1 + par.l$extension_x_limits)
limit_min = 0

radial_coeff = (180*180)/limit_max

volc.new.neg = volc.new[which(volc.new$weighted_mean < 0),]
volc.new.pos = volc.new[which(volc.new$weighted_mean > 0),]
volc.new.neg$radial = -((radial_coeff*abs(volc.new.neg$weighted_mean))/180) + 360
volc.new.pos$radial = (radial_coeff*abs(volc.new.pos$weighted_mean))/180

volc.new.comb = rbind(volc.new.neg,volc.new.pos)

# TODO: min_T_stat
min_T_stat = 2.5

volc.new.comb$sign = as.factor(ifelse(abs(volc.new.comb$weighted_Tstat) > min_T_stat, TRUE, FALSE))

if (!is.null(par.l$activactorRepressorObject)) {
  
  ## this is AR reporessors part 
  load(par.l$activactorRepressorObject)
  
  stopifnot(!is.null(median.cor.tfs))
  stopifnot(!is.null(act.rep.thres))
  
  AR.data = as.data.frame(median.cor.tfs)
  AR.data$TF = rownames(AR.data)
  volc.new.comb2 = merge(volc.new.comb, AR.data, by = "TF",all.x = T)
  
  volc.new.comb2$AR2 = ifelse(is.na(volc.new.comb2$median.cor.tfs), "not-expressed",
                              ifelse(volc.new.comb2$median.cor.tfs < act.rep.thres[1], "repressor",
                                     ifelse(volc.new.comb2$median.cor.tfs > act.rep.thres[2], "activator", "undetermined")))
  
  volc.new.comb3 = merge(volc.new.comb,volc.new.comb2, by = "TF" )
  volc.new.comb = volc.new.comb3
  volc.new.comb$AR2 = factor(volc.new.comb$AR2, levels = c("activator","undetermined","repressor","not-expressed"))
  
  ggrepel_df = volc.new.comb[which(volc.new.comb$sign.x == TRUE ),]
  
} else {
  
  volc.new.comb$weighted_Tstat.x = volc.new.comb$weighted_Tstat
  volc.new.comb$radial.x = volc.new.comb$radial
  volc.new.comb$sign.x = volc.new.comb$sign
  
}

############### new plot 
## scaling to the angles automatic 

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
ylimit = max(log10(abs(volc.new.comb$weighted_Tstat.x))) + par.l$extension_y_limits

fillVar  = NULL
if (!is.null(par.l$activactorRepressorObject)) {
  fillVar = volc.new.comb$AR2
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

if (!is.null(par.l$activactorRepressorObject)) {
  p1 = p1 +  geom_point(aes(x = volc.new.comb$radial.x, 
                            y = log10(abs(volc.new.comb$weighted_Tstat.x)),
                            label = volc.new.comb$TF,
                            size = volc.new.comb$sign.x,
                            fill = fillVar,
                            color = fillVar), 
                        alpha = 1, shape = 21)
} else {
  p1 = p1 +  geom_point(aes(x = volc.new.comb$radial.x, 
                            y = log10(abs(volc.new.comb$weighted_Tstat.x)),
                            label = volc.new.comb$TF,
                            size = volc.new.comb$sign.x), 
                        alpha = 1, shape = 21)
}

p1 = p1 + coord_polar() + 
  # limits of the angular plot 
  scale_x_continuous(limits = c(0,360)) +
  # sizes of the points for significance
  scale_size_manual(values = c(0.5,1), guide = FALSE) 

if (!is.null(par.l$activactorRepressorObject)) {
  p1 = p1 + scale_color_manual(values = par.l$color_points, guide = FALSE) + 
            geom_label_repel( aes(x = ggrepel_df$radial.x,
                          y = log10(abs(ggrepel_df$weighted_Tstat.x)),
                          label = ggrepel_df$TF,
                          fill = fillVar),
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
  
if (!is.null(par.l$activactorRepressorObject)) {
  p1 = p1 + 
       scale_fill_manual(values = par.l$color_points) + 
       labs(size = "Significance", fill = "Role:") 
}
  
p1 = p1 + guides(fill = guide_legend( override.aes = list(size = 8)))



# extracting the labels from object in order to add them manually later
labels = list()
labels$range = layer_scales(p1)$y$range$range
labels$sum = sum(abs(labels$range[1]) + abs(labels$range[2]))
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



p2 = p2 + geom_segment(aes(x=df.axis$angles, xend = df.axis$angles , y=min(labels$breaks), yend = max(labels$breaks)), 
                       size = 0.7,linetype = "dotted", color = "grey50",alpha = .4) 
for (i in df.axis$angles) {
  p2 = p2 + annotate(geom = "text", x = i, y = 2.35, label = paste0(df.axis[which(df.axis$angles == i),]$annotation), vjust = 0.5, size = 4.5)
}

p3 = p2 

for (i in 2:length(labels$breaks)) {
  p3 = p3 + annotate(geom = "text", x = 167, y = labels$breaks[i], label = paste0(labels$breaks[i]), vjust = 0.5,hjust = 1, angle = 15) 
}

p3 = p3 + annotate(geom = "text", x = 183, y = 0, label = "log10 of T-statistic", vjust = -1.5, angle = 80) + 
  geom_segment(aes(x=10, xend = 90 , y=ylimit + 1, yend = ylimit + 1), size=0.3,
               arrow = arrow(length = unit(0.6,"cm"))) + 
  annotate(geom = "text", x = 0, y = ylimit + 1, label = " TF activity ", vjust = 0, angle = 0, size = 4.5) + 
  
  annotate(geom = "text", x = 30, y = ylimit + 1.25, label = strsplit(conditionComparison, " ")[[1]][[1]], 
           vjust = 0, angle = -26, size = 4.5)+ 
  annotate(geom = "text", x = 330, y = ylimit + 1.25, label =  strsplit(conditionComparison, " ")[[1]][[3]], 
           vjust = 0, angle = 26, size = 4.5) + 
  geom_segment(aes(x=350, xend = 270 , y=ylimit + 1, yend = ylimit + 1), size=0.3,
               arrow = arrow(length = unit(0.6,"cm"))) 

# for loop for minor y circular lines 
for (x in 1:length(labels$breaks)){
  p3 = p3 + geom_hline(yintercept = labels$breaks[x],
                       size = 1,
                       linetype = "dotted", 
                       color = "grey50",
                       alpha = .4) 
  
}
p3
