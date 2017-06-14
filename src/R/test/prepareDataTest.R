start.time  <-  Sys.time()

#########################
# LIBRARY AND FUNCTIONS #
#########################

source("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/functions.R")

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "vsn", "csaw", "checkmate", "limma", "tools"), verbose = TRUE)


###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$doVSTTransformation = FALSE


par.l$log_minlevel = "INFO"

############################################
# READ AND VALIDATE COMMAND LINE ARGUMENTS #
############################################
args <- commandArgs(trailingOnly = TRUE)

# args = c(
# 
#   "/scratch/leyva/PWM/input/sampleTable_all.tab",
#   "/scratch/carnold/Age/PEAKS/Age.sampleMetadata.rds"  ,
#   "/scratch/carnold/Age/PEAKS/Age.peaks.rds" ,
#   "/scratch/carnold/Age/PEAKS/Age.peaks.tsv"   ,
#   "/scratch/carnold/Age/TEMP/extension100/conditionComparison.rds" ,
#   "/scratch/carnold/Age/PEAKS/Age.normFacs.rds" ,
#   "/scratch/carnold/Age/PEAKS/Age.MAplot.peaks.pdf" ,
#   "/scratch/carnold/Age/PEAKS/Age.meanSD.peaks.pdf",
#   "/scratch/leyva/PWM/output/Age/PEAKS",
#   ".overlapPeaks.bed" ,
#   "old,young" ,
#   "~ Treatment + Condition",
#   "/scratch/carnold/prepareData.R.log"
# 
# )

# args = c("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/data/IvanNew/sampleTable_TET2.GMP.all.csv"   ,             
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/outputReal/PEAKS/sampleMetadata.rds"             ,   
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/outputReal/PEAKS/peaks.rds"                 ,
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/outputReal/PEAKS/peaks.tsv"                 ,
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/outputReal/PEAKS/normFacs.rds"                 ,
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/outputReal/PEAKS/MAplot.peaks.pdf"            ,    
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/outputReal/PEAKS/VSTvsReal.peaks.pdf"           ,     
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/outputReal/PEAKS"              ,   
# ".overlapPeaks.bed"                ,
# "KO,WT"                 ,
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/outputReal/Logs_and_Benchmarks/prepareData.R.log"
# )

args = c("/scratch/leyva/PWM3/input/sampleTable_wo_275_421.tab"   ,
"/scratch/leyva/PWM3/output/Age/PEAKS/Age.sampleMetadata.rds"             ,
"/scratch/leyva/PWM3/output/Age/PEAKS/Age.peaks.rds"                 ,
"/scratch/leyva/PWM3/output/Age/PEAKS/Age.peaks.tsv"                 ,
"/scratch/leyva/PWM3/output/Age/TEMP/extension50/conditionComparison.rds"                 ,
"/scratch/leyva/PWM3/output/Age/PEAKS/Age.normFacs.rds"            ,
"/scratch/leyva/PWM3/output/Age/PEAKS/Age.MAplot.peaks.pdf"           ,
"/scratch/leyva/PWM3/output/Age/PEAKS/Age.meanSD.peaks.pdf"              ,
"/scratch/leyva/PWM3/output/Age/PEAKS",
".overlapPeaks.bed"                ,
"old,young"                 ,
"~ Treatment + conditionSummary",
"/scratch/leyva/PWM3/output/Age/Logs_and_Benchmarks/prepareData.R.log2"
)


if (length(args) != 13) {
  stop("Expecting 13 arguments but found ", length(args),". Exiting.")
} else {
  par.l$file_input_sampleData      = args[1]
  par.l$file_output_metadata       = args[2]
  par.l$file_output_peaks          = args[3]
  par.l$file_output_peaksTSV       = args[4]
  par.l$file_output_condComp       = args[5]  
  par.l$file_output_normFacs       = args[6]
  par.l$file_output_plot_MA        = args[7]
  par.l$file_output_plot_VSTvsReal = args[8]
  par.l$dir_peaks                  = args[9]
  par.l$pattern_peaks              = args[10]
  par.l$conditionComparison        = args[11]
  par.l$designFormula              = args[12]
  par.l$file_output_log            = args[13]
}



#####################
# VERIFY PARAMETERS #
#####################

assertDirectoryExists(par.l$dir_peaks, access = "r")
assertFileExists(par.l$file_input_sampleData)

assertCharacter(par.l$pattern_peaks, len = 1)

assertDirectoryExists(dirname(par.l$file_output_log), access = "w")
assertCharacter(par.l$conditionComparison, len = 1, min.chars = 3)


allDirs = c(dirname(par.l$file_output_metadata), 
            dirname(par.l$file_output_peaks), 
            dirname(par.l$file_output_normFacs), 
            dirname(par.l$file_output_plot_MA), 
            dirname(par.l$file_output_plot_VSTvsReal),
            dirname(par.l$file_output_peaksTSV)
            )

for (dirname in unique(allDirs)) {
  
  if (!testDirectoryExists(dirname)) {
    dir.create(dirname, recursive = TRUE)
  } else {
    assertDirectoryExists(dirname, access = "w")
  }
  
}

######################
# FINAL PREPARATIONS #
startLogger(par.l$file_output_log, par.l$log_minlevel,  removeOldLog = TRUE)
printParametersLog(par.l)




#################
# READ METADATA #
#################

# sampleData.df = read_tsv(par.l$file_input_sampleData, col_types =
#                            cols(
#                              SampleID = col_character(),
#                              bamReads = col_character(),
#                              Peaks = col_character(),
#                              Treatment = col_character(),
#                              Condition = col_factor(levels = c("PAT", "CONT")),
#                              Gender = col_factor(levels = c("M", "F")),
#                              Mutation = col_factor(levels = c("MUT", "WT")),
#                              conditionSummary = col_character()
#                            ))

sampleData.df = read_tsv(par.l$file_input_sampleData, col_names = TRUE)

assertSubset(c("bamReads", "conditionSummary"), colnames(sampleData.df))

#sampleData.df = read_tsv("/scratch/leyva/PWM3/input/sampleTable_wo_275_421.tab", col_names = TRUE)

# Check if contrasts have been specified correctly
conditionsContrast = strsplit(par.l$conditionComparison, split = ",", fixed = TRUE)[[1]]
assertVector(conditionsContrast, len = 2)
assertSubset(conditionsContrast, sampleData.df$conditionSummary)

# EDIT; Commented out to not enforce particular column names
# mutate(Treatment = as.factor(Treatment)) %>%

# Read and modify samples metadata
sampleData.df = sampleData.df %>%
                    mutate(name = file_path_sans_ext(basename(sampleData.df$bamReads))) %>%
                    filter(conditionSummary %in% conditionsContrast) %>%
                    mutate(conditionSummary = as.factor(conditionSummary))

saveRDS(sampleData.df, par.l$file_output_metadata)

##############################
# ITERATE THROUGH PEAK FILES #
##############################

#On this step it depends what waas the initial peaks file(how many columns)
peaks.l = list()
coverageAll.m = NULL

flog.info(paste0("Iterating over ", length(sampleData.df$name), " peak files in ", par.l$dir_peaks))

for (nameCur in sampleData.df$name) {

  file_output_PEAKS = paste0(par.l$dir_peaks, "/", nameCur, par.l$pattern_peaks)
  flog.info(paste0(" Parsing file ", file_output_PEAKS))
  assertFileExists(file_output_PEAKS)

  peaks.l[[nameCur]] =  read_tsv(file_output_PEAKS, col_names = c("chr", "PSS", "PES", "annotation", "ID", "coverage"), col_types = cols())
  peaks.l[[nameCur]]$identifier = paste0(peaks.l[[nameCur]]$chr,":", peaks.l[[nameCur]]$PSS,"-", peaks.l[[nameCur]]$PES)

  # Filter and retain only unique identifiers
  peaks.filtered.df = distinct(peaks.l[[nameCur]], identifier, .keep_all = TRUE)
  
  nRowsFiltered = nrow(peaks.l[[nameCur]]) - nrow(peaks.filtered.df)
  if (par.l$verbose & nRowsFiltered  > 0) flog.info(paste0("Filtered ", nRowsFiltered, " non-unique positions out of ", nrow(peaks.l[[nameCur]]), " from peaks table."))
  
  peaks.l[[nameCur]] = peaks.filtered.df

  # concatenate results from COV from each iteration
  coverageAll.m = cbind(coverageAll.m, peaks.l[[nameCur]]$coverage)
}

# TODO: Should contain all of peaks.l[[nameCur]] and not only the last one

# add step to add ID MARIANA
original.peaks.df = read_tsv("/scratch/leyva/PWM3/input/cons_ATACpeaks_LinearModel_5col.bed", col_names = F)
colnames(original.peaks.df) = c("chr","PSS","PES","annotation","ID")
original.peaks.df$NR = rownames(original.peaks.df)

peaks.filtered.corr.df = peaks.filtered.df %>% 
  dplyr::full_join(original.peaks.df, by = "annotation") %>% 
  dplyr::mutate(ID = NR, chr = chr.x, PSS = PSS.x, PES = PES.x ) %>% 
  dplyr::select(-one_of(c("chr.x","chr.y","PSS.x","PSS.y","PES.x","PES.y","NR","ID.x","ID.y"))) %>%
  na.omit() %>%
  select(one_of(c("chr","PSS","PES","annotation","ID","identifier","coverage")))

peaks.filtered.df = peaks.filtered.corr.df 
### 
saveRDS(peaks.filtered.df, file = par.l$file_output_peaks)

## transform as matrix data frame with counts
coverageAll.m = as.matrix(coverageAll.m)
colnames(coverageAll.m) = sampleData.df$name # TODO
rownames(coverageAll.m) = peaks.l[[1]]$identifier # Take the first element as nameCurresentative, they are all identical anyway

#############
# RUN DESEQ #
#############

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

# get Deseq object


cds.peaks <- DESeqDataSetFromMatrix(countData = coverageAll.m,
                                    colData = sampleData.df,
                                    design = designFormula)

counts(cds.peaks) = counts(cds.peaks)
# vst data dont use data about modelling of the linear model in our case batch
# Normalize with LOESS
normFacs <- exp(normOffsets(counts(cds.peaks),
                            lib.sizes = colSums(counts(cds.peaks)),
                            type = "loess"))
rownames(normFacs) = rownames(coverageAll.m)
# add normalization factor
normalizationFactors(cds.peaks) <- normFacs

saveRDS(normFacs, par.l$file_output_normFacs)




# filter
# 02.12.16 add rowMeans > par.l$filter_rowMeansThreshold to check fot the peaks (Armando mail)

cds.peaks.filt = cds.peaks[rowMeans(counts(cds.peaks)) > 0, ]
#cds.peaks.filt$conditionSummary = factor(cds.peaks.filt$conditionSummary, levels = unique(sampleData.df$conditionSummary))


# Deseq analysis
cds.peaks <- DESeq(cds.peaks.filt, fitType = 'local')

# Save the comparison that DeSeq made for later scripts
comparisonDESeq = getComparisonFromDeSeqObject(cds.peaks, par.l$designFormula)
saveRDS(comparisonDESeq, file = par.l$file_output_condComp)

cds.peaks.df <- as.data.frame(DESeq2::results(cds.peaks))

###############################################
# VST TRANSFORMATION AND BATCH EFFECT REMOVAL #
###############################################
if (par.l$doVSTTransformation) {
  # TODO: Maybe remove, not needed after
  # do the vst transfromation
  vsd.peaks.raw = DESeq2::varianceStabilizingTransformation(cds.peaks, fitType = 'local', blind = FALSE)
  # extract vst log2FC
  peaks.df_row = as.data.frame(assay(vsd.peaks.raw))

  # TODO: generalize
  if (length(formulaVariables) > 1) {
    
  }
  stop("Not working")
  peaks.df_row = as.data.frame(removeBatchEffect(peaks.df_row, batch = as.character(sampleData.df$Treatment)))
 
  colname_condition1 = paste0("mean_",conditionsContrast[1])
  colname_condition2 = paste0("mean_",conditionsContrast[2])
  colnames1 = sampleData.df$name[which(sampleData.df$conditionSummary == conditionsContrast[1])]
  colnames2 = sampleData.df$name[which(sampleData.df$conditionSummary == conditionsContrast[2])]
  peaks.df_row[,colname_condition1] = apply(peaks.df_row[, colnames(peaks.df_row) %in% colnames1],1, mean)
  peaks.df_row[,colname_condition2] = apply(peaks.df_row[, colnames(peaks.df_row) %in% colnames2],1, mean)
  peaks.df_row$diff = round(peaks.df_row[,colname_condition1] - peaks.df_row[,colname_condition2], 3)
}


final.peaks.df = data_frame("position"    = rownames(cds.peaks.df), 
                            "D2_baseMean" = cds.peaks.df$baseMean,
                            "D2_l2FC"     = cds.peaks.df$log2FoldChange,
                            "D2_ldcSE"    = cds.peaks.df$lfcSE,
                            "D2_stat"     = cds.peaks.df$stat,
                            "D2_pval"     =  cds.peaks.df$pvalue, 
                            "D2_padj"     =  cds.peaks.df$padj #, 
                            #"vst_diff"    =  peaks.df_row$diff
                            )

# assign final.peaks.df to the peaks.df and filter away NAs at the p.adjust
#peaks.df = final.peaks.df[complete.cases(final.peaks.df[,7]),]
# 23.02 talk with Armando and Judith decided to bring NAs back

################
# WRITE OUTPUT #
################

write_tsv(final.peaks.df, path = par.l$file_output_peaksTSV)

#########
# PLOTS #
#########

pdf(par.l$file_output_plot_MA)
DESeq2::plotMA(cds.peaks)
dev.off()

# see how change the distr after transformation
notAllZeroPeaks <- (rowSums(counts(cds.peaks)) > 0)


pdf(par.l$file_output_plot_VSTvsReal)
meanSdPlot(assay(cds.peaks[notAllZeroPeaks,]))
# meanSdPlot(assay(vsd.peaks.raw[notAllZeroPeaks,]))
dev.off()


end.time  <-  Sys.time()
message(" Finished execution. TOTAL RUNNING TIME: ", round(end.time - start.time, 1), " ", units(end.time - start.time),"\n")
