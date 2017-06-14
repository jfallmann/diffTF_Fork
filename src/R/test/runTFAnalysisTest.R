start.time  <-  Sys.time()

#########################
# LIBRARY AND FUNCTIONS #
#########################

source("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/functions.R")

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "vsn", "modeest", "checkmate", "limma"), verbose = TRUE)


###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"
par.l$doVSTTransformation = FALSE
par.l$signThresholdMAPlot = 0.1


args <- commandArgs(trailingOnly = TRUE)
# 
# #"RUNX3", "JUN"
# args = c("/scratch/leyva/PWM/output/Age/PEAKS/Age.sampleMetadata.rds"        ,
# "/scratch/leyva/PWM/output/Age/PEAKS/Age.normFacs.rds"     ,
# "/scratch/leyva/PWM/output/Age/PEAKS/Age.peaks.rds"         ,
# "/scratch/leyva/PWM/output/Age/PEAKS/Age.peaks.tsv"       ,
# "/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.187_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.187_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.275_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.275_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.385_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.385_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.387_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.387_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.391_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.391_rep3.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.395_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.395_rep3.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.400_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.400_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.403_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.403_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.409_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.409_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.414_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.414_rep3.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.415_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.415_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.418_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.418_rep3.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.419_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.419_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.419_rep3.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.419_rep4.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.421_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.421_rep2.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.421_rep3.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.421_rep4.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.436_rep1.final.overlapPeaks.bed,/scratch/leyva/PWM/output/Age/TF-SPECIFIC/RUNX3/extension100/RUNX3.Age.436_rep2.final.overlapPeaks.bed"       ,
# "RUNX3.Age.output.tsv"    ,
# "RUNX3.Age.summary.rds"  ,
# "RUNX3.Age.MA.realcounts.pdf"   ,
# "RUNX3.Age.meanSD.pdf"    ,
# "RUNX3.Age.log2.dens.pdf"   ,
# "RUNX3.Age.ECDF.pdf"      ,
# "RUNX3"                ,
# "~Treatment + conditionSummary"                  ,
# "testOut.log"
# )
# 

args = c(
  "/scratch/bunina/TFpipeline/wt129/output0vs4/PEAKS/day0vs4.sampleMetadata.rds",                     
  "/scratch/bunina/TFpipeline/wt129/output0vs4/PEAKS/day0vs4.normFacs.rds",                     
  "/scratch/bunina/TFpipeline/wt129/output0vs4/PEAKS/day0vs4.peaks.rds" ,                    
  "/scratch/bunina/TFpipeline/wt129/output0vs4/PEAKS/day0vs4.peaks.tsv" ,                    
  "/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.day0_rep1_run1.final.s.overlapPeaks.bed,/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.day0_rep2_run1.final.s.overlapPeaks.bed,/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.day0_rep1.final.s.overlapPeaks.bed,/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.day0_rep2.final.s.overlapPeaks.bed,/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.day4_rep1_run1.final.s.overlapPeaks.bed,/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.day4_rep2_run1.final.s.overlapPeaks.bed,/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.day4_rep1.final.s.overlapPeaks.bed,/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.day4_rep2.final.s.overlapPeaks.bed"                    ,
  "/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.output.tsv" ,                   
  "/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.summary.rds",                  
  "/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.MA.realcounts.pdf" ,         
  "/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.meanSD.pdf" ,      
  "/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.log2.dens.pdf",   
  "/scratch/bunina/TFpipeline/wt129/output0vs4/TF-SPECIFIC/TCF7/extension50/TCF7.day0vs4.ECDF.pdf"   ,                 
  "TCF7"                ,    
  "~ Treatment + Condition"                     ,
  "/scratch/bunina/TFpipeline/wt129/output0vs4/Logs_and_Benchmarks/analyzeTF.TCF7.R.log"
  
  
  
)



# args = c("/scratch/leyva/PWM3/output/Age/PEAKS/Age.sampleMetadata.rds"    ,               
# "/scratch/leyva/PWM3/output/Age/PEAKS/Age.normFacs.rds"     ,                
# "/scratch/leyva/PWM3/output/Age/PEAKS/Age.peaks.rds"      ,               
# "/scratch/leyva/PWM3/output/Age/PEAKS/Age.peaks.tsv"     ,                
# "/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.382_rep1.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.382_rep2.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.387_rep1.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.387_rep2.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.409_rep1.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.409_rep2.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.385_rep1.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.385_rep2.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.400_rep1.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.400_rep2.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.436_rep1.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.436_rep2.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.414_rep2.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.414_rep3.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.391_rep1.final.s.overlapPeaks.bed,/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.391_rep3.final.s.overlapPeaks.bed"                     ,
# "/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.output.tsv"  ,                  
# "/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.summary.rds"  ,                  
# "/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.MA.realcounts.pdf"   ,           
# "/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.meanSD.pdf"      ,             
# "/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.log2.dens.pdf"  ,             
# "/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/EPAS1/extension50/EPAS1.Age.ECDF.pdf"    ,                
# "EPAS1"                   ,
# "~ Treatment + conditionSummary"                    ,
# "/scratch/leyva/PWM3/output/Age/Logs_and_Benchmarks/analyzeTF.EPAS1.R.log"
# )


if (length(args) != 11) {
  stop("Expecting 11 arguments but found ", length(args),". Exiting.")
} else {
  par.l$file_input_metadata        = args[1]
  par.l$file_input_normFacs        = args[2]
  par.l$file_input_peaks           = args[3]
  par.l$file_input_peak2           = args[4]
  par.l$file_input_peakTFOverlaps  = args[5]
  par.l$file_output_summaryAll     = args[6]
  par.l$file_output_summaryStats   = args[7]
  par.l$file_output_plotsAll       = args[8]
  par.l$TF                         = args[9]
  par.l$designFormula              = args[10]
  par.l$file_log                   = args[11]
}




#####################
# VERIFY PARAMETERS #
#####################


assertFileExists(par.l$file_input_peak2)
assertFileExists(par.l$file_input_peaks)
assertFileExists(par.l$file_input_metadata)
assertFileExists(par.l$file_input_normFacs)

assertDirectoryExists(dirname(par.l$file_log), access = "w")


allDirs = c(dirname(par.l$file_output_summaryAll), 
            dirname(par.l$file_output_summaryStats), 
            dirname(par.l$file_output_plotsAll)
)

for (dirname in unique(allDirs)) {
  
  if (!testDirectoryExists(dirname)) {
    dir.create(dirname, recursive = TRUE)
  } else {
    assertDirectoryExists(dirname, access = "w")
  }
  
}

fileList = strsplit(par.l$file_input_peakTFOverlaps, split = ",", fixed = TRUE)[[1]]

for (fileCur in fileList) {
  assertFileExists(fileCur, access = "r")
}

######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel, appenderName = "file", removeOldLog = TRUE)
printParametersLog(par.l)


#################
# READ METADATA #
#################

sampleData.df = readRDS(par.l$file_input_metadata)


coverageAll.df = NULL

flog.info(paste0("Iterating over ", length(fileList), " TF-peak overlap files"))


for (fileCur in fileList) {
  
  peaksCur.df = read_tsv(fileCur, col_names = FALSE, col_types = cols())
  assertDataFrame(peaksCur.df, ncols = 7)
  #assertDataFrame(peaksCur.df, ncols = 6)
  
  colnames(peaksCur.df) = c("chr","MSS","MES","annotation","ID","strand","coverage")
  #colnames(peaksCur.df) = c("chr","MSS","MES","annotation","ID","coverage")
  
  peaksCur.df = mutate(peaksCur.df, identifier = paste0(chr,":", MSS, "-",MES))
  
  # get only unique identifiers
  peaksCur.df =  distinct(peaksCur.df, identifier, .keep_all = TRUE)
  
  # Only use the "coverage column and append to the final list
  coverageAll.df = dplyr::bind_cols(coverageAll.df, peaksCur.df[,"coverage"])
  
}

# Add row mean
coverageAll.df$mean = apply(coverageAll.df, 1, mean)

# Add metadata (all except coverage)
columns = c("chr","MSS","MES","annotation","ID","strand","identifier")
#columns = c("chr","MSS","MES","annotation","ID","identifier")
assertSubset(columns, colnames(peaksCur.df))
coverageAll.df = dplyr::bind_cols(coverageAll.df, peaksCur.df[, columns])
colnames(coverageAll.df) = c(sampleData.df$name, "mean", columns)


# add step to add ID MARIANA
original.peaks.df = read_tsv("/scratch/leyva/PWM3/input/cons_ATACpeaks_LinearModel_5col.bed", col_names = F)
colnames(original.peaks.df) = c("chr","PSS","PES","annotation","ID")
original.peaks.df$NR = rownames(original.peaks.df)

coverageAll.corr.df = coverageAll.df %>% 
  dplyr::full_join(original.peaks.df, by = "annotation") %>% 
  dplyr::mutate(ID = NR, chr = chr.x) %>% 
  dplyr::select(-one_of(c("chr.x","chr.y","PSS","PES","NR","ID.x","ID.y"))) %>%
  na.omit()

coverageAll.df = coverageAll.corr.df 
# Group by ID
coverageAll_grouped.df = coverageAll.df %>%
  dplyr::group_by(ID) %>%
  dplyr::slice(which.max(mean))

TF.table.m = as.matrix(coverageAll_grouped.df[,sampleData.df$name])
colnames(TF.table.m) = sampleData.df$name
rownames(TF.table.m) = coverageAll_grouped.df$identifier


# Create formula based on user-defined design
designFormula = tryCatch({
  as.formula(par.l$designFormula)
}, warning = function(w) {
  stop("Converting the design formula \"", par.l$designFormula, "\" created a warning, which should be checked carefully.")
}, error = function(e) {
  stop("Design formula \"", par.l$designFormula, "\" not valid")
})

# Check colmn names
formulaVariables = attr(terms(designFormula), "term.labels")
assertSubset(formulaVariables, colnames(sampleData.df))

# create Deseq object from the TF specific data
TF.cds <- DESeqDataSetFromMatrix(countData = TF.table.m,
                                 colData = sampleData.df,
                                 design = designFormula)

# normalize with normFacs
normFacs = readRDS(par.l$file_input_normFacs)
normalizationFactors(TF.cds) <- normFacs[as.numeric(coverageAll_grouped.df$ID),]
# low RC, check by rowMean
TF.cds.filt = TF.cds[rowMeans(counts(TF.cds)) > 0, ]



# Deseq 2 functions
# off/on shrinkage of the log2FC - betaPrior
# with the simulations of negative binomial distribution increase the sample size
res_DESeq <- DESeq(TF.cds.filt,fitType = 'local')

res_DESeq = tryCatch( {
  DESeq(TF.cds.filt,fitType = 'local')
  
}, error = function(e) {
  warning("Warning: Could not run DESeq with local fitting, retry with default fitting type...")
  DESeq(TF.cds.filt)
}
)



res_DESeq.df <- as.data.frame(DESeq2::results(res_DESeq))

# TODO: Extract the correct order of the comparison via resultsNames(res_DESeq). The last element is the one we want to capture


if (par.l$doVSTTransformation) {
  
  #vst transformation
  vsd.TF <- DESeq2::varianceStabilizingTransformation(res_DESeq,fitType = 'local', blind = FALSE)
  # extract log2Fc from vsd object
  TF_row.df = as.data.frame(assay(vsd.TF))
  # remove batch effects
  TF_row.df = as.data.frame(removeBatchEffect(TF_row.df, batch = as.character(sampleData.df$Treatment)))
  # TODO: Does not work anymore, conditionsContrast has been removed
  colname_condition1 = paste0("mean_",conditionsContrast[1])
  colname_condition2 = paste0("mean_",conditionsContrast[2])
  colnames1 = sampleData.df$name[which(sampleData.df$conditionSummary == conditionsContrast[1])]
  colnames2 = sampleData.df$name[which(sampleData.df$conditionSummary == conditionsContrast[2])]
  TF_row.df[,colname_condition1] = apply(TF_row.df[, colnames(TF_row.df) %in% colnames1],1, mean)
  TF_row.df[,colname_condition2] = apply(TF_row.df[, colnames(TF_row.df) %in% colnames2],1, mean)
  TF_row.df$diff = round(TF_row.df[,colname_condition1] - TF_row.df[,colname_condition2], 3)
  
  
}
# addition  02.06
final.TF.df = data_frame("position"    = rownames(res_DESeq.df), 
                         "D2_baseMean" = res_DESeq.df$baseMean,
                         "D2_l2FC"     = res_DESeq.df$log2FoldChange,
                         "D2_ldcSE"    = res_DESeq.df$lfcSE,
                         "D2_stat"     = res_DESeq.df$stat,
                         "D2_pval"     = res_DESeq.df$pvalue, 
                         "D2_padj"     = res_DESeq.df$padj#, 
                         #"vst_diff"    = TF_row.df$diff
)

# assign final.peaks.df to the peaks.df and filter away NAs at the p.adjust

peaksFiltered.df = readRDS(par.l$file_input_peaks)


assertSubset(rownames(res_DESeq.df), peaksCur.df$identifier)


rm_col = c("coverage.x","chr.y","coverage.y","identifier.y" )
order = c("TF","chr","MSS","MES","strand", "PSS","PES","annotation","ID", "identifier","baseMean", "log2FoldChange","lfcSE","stat", "pvalue","padj")

# dplyr::mutate(TF = par.l$TF) gives the following weird error message: Error: Unsupported type NILSXP for column "TF"
TFCur = par.l$TF


TF_output.df = res_DESeq.df %>%
  rownames_to_column(var = "identifier") %>%
  dplyr::full_join(peaksCur.df,by = c("identifier")) %>%
  dplyr::full_join(peaksFiltered.df, by = "annotation") %>%
  dplyr::filter(!is.na(baseMean)) %>%
  dplyr::rename(ID = ID.y, chr = chr.x, identifier = identifier.x) %>%
  dplyr::select(-one_of(rm_col)) %>%
  dplyr::mutate(TF = TFCur) %>%
  dplyr::select(one_of(order)) %>%
  dplyr::arrange(chr)


write_tsv(TF_output.df, path = par.l$file_output_summaryAll)



# d) Comparisons between peaks and binding sites


# TODO: Not needed
#peaks_C = nrow(peaks.df[peaks.df$log2FoldChange > 0,])/nrow(peaks.df)

peaks.df = read_tsv(par.l$file_input_peak2, col_types = cols())

modeNum     = mlv(round(final.TF.df$D2_l2FC, 2), method = "mfv", na.rm = TRUE)
Ttest       = t.test(final.TF.df$D2_l2FC, peaks.df$D2_l2FC)

output.df = data_frame(TF              = par.l$TF,
                       Pos_l2FC    = nrow(final.TF.df[final.TF.df$D2_l2FC > 0,]) / nrow(final.TF.df),
                       Mean_l2FC   = mean(final.TF.df$D2_l2FC, na.rm = TRUE),
                       Median_l2FC = median(final.TF.df$D2_l2FC, na.rm = TRUE),
                       Mode_l2FC   = modeNum[[1]],
                       Ttest_pval  = Ttest$p.value,
                       Modeskewness    = modeNum[[2]], 
                       T_statistic     = Ttest$statistic[[1]], 
                       TFBS_num        = nrow(final.TF.df)
)


saveRDS(output.df,file = par.l$file_output_summaryStats)

############
############
## GRAPHS ##
############
############

pdf(par.l$file_output_plotsAll)

######
# MA #
######

# MA plot for TF
#pdf(par.l$file_output_plot_MA)
DESeq2::plotMA(res_DESeq, main = "MA plot", alpha = par.l$signThresholdMAPlot)

#dev.off()

#############
# VSTvsReal #
#############

# transformation for the TF
notAllZeroTF <- (rowSums(counts(res_DESeq)) > 0)
#pdf(par.l$file_output_plotvsReal)
title = "not supported for meanSdPlot"
meanSdPlot(assay(res_DESeq[notAllZeroTF,]))
#meanSdPlot(assay(vsd.TF[notAllZeroTF,]))
#dev.off()


comparisonDESeq = getComparisonFromDeSeqObject(res_DESeq, par.l$designFormula)
xlabLabel = paste0(" log2 FC ", comparisonDESeq)

# density plot ## nice addition 27.04
TF_dens = ggplot() + geom_density(aes(x = peaks.df$D2_l2FC,fill = "A" ),
                                  alpha = .5, color = "black") +
  geom_density(aes(x = final.TF.df$D2_l2FC,fill = "B"),size = 1, alpha = .7) +
  xlab(xlabLabel) +
  ggtitle("TF density") + 
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
#ggsave(plot = TF_dens, filename = par.l$file_output_plot, width = 4, height = 4, useDingbats = FALSE, dpi = 600)

########
# ECDF #
########

# create ecdf plots for each TF
ECDF_TF = ggplot() + 
  stat_ecdf(aes(x = final.TF.df$D2_l2FC,colour = paste0("", par.l$TF))) +
  stat_ecdf(aes(x = peaks.df$D2_l2FC,colour = "Peaks" )) +
  ggtitle("ECDF") + 
  xlab("Log2FC WT/KO") + 
  guides(colour = guide_legend(title = "ORIGIN"))

plot(ECDF_TF)
#ggsave(plot = ECDF_TF, filename = par.l$file_output_plot_ecdf, width = 6, height = 4, useDingbats = FALSE, dpi = 600)

dev.off()

