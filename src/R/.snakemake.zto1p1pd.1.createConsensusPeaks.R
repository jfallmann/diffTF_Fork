
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('data/peaks/GMP.WT.1.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/GMP.WT.2.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/GMP.WT.3.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/GMP.WT.4.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/MPP.WT.1.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/MPP.WT.2.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/MPP.WT.3.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/MPP.WT.4.final.sort.CGcorr_peaks.narrowPeak', "peaks" = c('data/peaks/GMP.WT.1.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/GMP.WT.2.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/GMP.WT.3.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/GMP.WT.4.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/MPP.WT.1.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/MPP.WT.2.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/MPP.WT.3.final.sort.CGcorr_peaks.narrowPeak', 'data/peaks/MPP.WT.4.final.sort.CGcorr_peaks.narrowPeak')),
    output = list('../output/PEAKS/GMPvsMPP.consensusPeaks_lengthDistribution.pdf', '../output/PEAKS/GMPvsMPP.consensusPeaks.bed', "summaryPlot" = '../output/PEAKS/GMPvsMPP.consensusPeaks_lengthDistribution.pdf', "consensusPeaks_bed" = '../output/PEAKS/GMPvsMPP.consensusPeaks.bed'),
    params = list(),
    wildcards = list(),
    threads = 1,
    log = list('../output/LOGS_AND_BENCHMARKS/1.produceConsensusPeaks.R.log'),
    resources = list(),
    config = list("par_general" = list("outdir" = '../output', "regionExtension" = 100, "comparisonType" = 'GMPvsMPP', "designContrast" = '~ conditionSummary', "designVariableTypes" = 'conditionSummary:factor', "nPermutations" = 3, "nBootstraps" = 0, "nCGBins" = 10, "TFs" = 'CTCF,CEBPB,SNAI2,CEBPA,UBIP1,CEBPG,CEBPD,ZFX,AP2D,PAX5.S,SNAI1,ZEB1,SP4,MBD2,IRF1,MECP2,PAX5.D,SP3,NFIA.C,SP1.A,IRF7,MYF6, NRF1,DBP,MAZ,NKX28,DLX2,GATA1,P53,ZN143,AIRE,NR2C2,HMGA1,FUBP1,TEAD3,OVOL1,HXD4,KLF1,RXRG,HNF1B,ZIC3,HNF1A,NANOG.S,GFI1,PO3F1,NR2C1,ELF5,TF65.C,NFAC3,TEAD1', "dir_scripts" = '../../src/R', "RNASeqIntegration" = TRUE), "samples" = list("summaryFile" = 'sampleData.tsv'), "additionalInputFiles" = list("refGenome_fasta" = 'referenceGenome/mm10.fa', "dir_PWMScan" = 'PWMscan', "RNASeqCounts" = 'data/RNA-Seq/RNA.counts.tsv', "HOCOMOCO_mapping" = '../../src/HOCOMOCO/translationTable_mouse.csv'), "peaks" = list("consensusPeaks" = '', "peakType" = 'narrow', "minOverlap" = 2)),
    rule = 'produceConsensusPeaks'
)
######## Original script #########
start.time  <-  Sys.time()



#########################
# LIBRARY AND FUNCTIONS #
#########################

library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DiffBind", "checkmate", "stats"), verbose = FALSE)

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/1.createConsensusPeaks.R.rds")
createDebugFile(snakemake)

###################
#### PARAMETERS ###
###################

par.l = list()

par.l$verbose = TRUE
par.l$log_minlevel = "INFO"

#####################
# VERIFY PARAMETERS #
#####################


assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(names(snakemake@input), c("", "peaks"))

par.l$file_input_peaks = snakemake@input$peaks


## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "consensusPeaks_bed", "summaryPlot"))

par.l$file_output_consensusPeaks = snakemake@output$consensusPeaks_bed
par.l$file_output_plot           = snakemake@output$summaryPlot


## CONFIG ##
assertList(snakemake@config, min.len = 1)

file_sampleData = snakemake@config$samples$summaryFile
assertFileExists(file_sampleData)

peakType = snakemake@config$peaks$peakType
assertSubset(peakType, c("raw", "bed", "narrow", "macs", "swembl", "bayes", "peakset", "fp4"), empty.ok = FALSE)

minOverlap = snakemake@config$peaks$minOverlap
assertIntegerish(minOverlap, lower = 0)

## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]


allDirs = c(dirname(par.l$file_output_consensusPeaks), 
            dirname(par.l$file_output_plot),
            dirname(par.l$file_log)
)

testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel,  removeOldLog = TRUE)
printParametersLog(par.l)

##########################
# Create consensus peaks #
##########################

# Provide the metadata file and parse the CSV here
sampleMetaData.df = read_tsv(file_sampleData, col_types = cols())
assertSubset(c("SampleID", "bamReads", "conditionSummary", "Peaks"), colnames(sampleMetaData.df))
assertIntegerish(minOverlap, upper = nrow(sampleMetaData.df))

sampleMetaData.df$Condition  = sampleMetaData.df$conditionSummary
sampleMetaData.df$PeakCaller = peakType
sampleMetaData.df$conditionSummary = sampleMetaData.df$Condition

# sampleMetaData.df = sampleMetaData.df[complete.cases(sampleMetaData.df),]

if (nrow(sampleMetaData.df) == 0) {
  message = paste0("No rows left to analyze, the data frame for DiffBind contains no rows after filtering for complete cases")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}
## DiffBind assigning. Convert to regular data frame first so DiffBind does not complain
dba <- dba(sampleSheet = as.data.frame(sampleMetaData.df))
consensus = dba.peakset(dba, minOverlap = minOverlap)

# Get the consensus peaks (last row)
consensusPeaks.df = as.data.frame(consensus$peaks[[nrow(sampleMetaData.df) + 1]])

colnames(consensusPeaks.df) = c("chr", "start", "end", "pvalue")
consensusPeaks.df$length = consensusPeaks.df$end - consensusPeaks.df$start

# Graph
# TODO: Make a bit nicer
g = ggplot(consensusPeaks.df) + geom_density(aes(x = length)) + scale_x_continuous(labels = seq(0,3000,250), breaks = seq(0,3000,250))
ggsave(plot = g, filename = par.l$file_output_plot)


consensusPeaks.df$name = paste0(consensusPeaks.df$chr, ":", consensusPeaks.df$start, "-", consensusPeaks.df$end)
consensusPeaks.df$id = rownames(consensusPeaks.df)

# Delete columns we don't need
consensusPeaks.df$pvalue = NULL
consensusPeaks.df$length = NULL

write_tsv(consensusPeaks.df, path = par.l$file_output_consensusPeaks, col_names = FALSE)

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)
