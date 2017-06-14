## README: this script is to remove CG bias by comparing permutations from the same CG bin 

start.time  <-  Sys.time()

#########################
# LIBRARY AND FUNCTIONS #
#########################

source("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/functions.R")

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "lsr", "ggrepel", "checkmate", "tools", "methods"), verbose = TRUE)

# methods needed here because Rscript does not loads this package automatically, see http://stackoverflow.com/questions/19468506/rscript-could-not-find-function


###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$FDR_threshold = 0.05
par.l$volcanoPlot_height = 16
par.l$volcanoPlot_width  = 24
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.2, 0.5, 0.8)
par.l$log_minlevel = "INFO"



############################################
# READ AND VALIDATE COMMAND LINE ARGUMENTS #
############################################
args <- commandArgs(trailingOnly = TRUE)

                 
# args = c("/media/carnold/s3/PAH/APC/OUTPUT/FINAL_OUTPUT/extension100/ABC.allMotifs.tsv"            ,     
# "/media/carnold/s3/PAH/APC/OUTPUT/TEMP/ABC.motifs.coord.nucContent.bed"                 ,
# "GATA6.B"                 ,
# "/media/carnold/s3/PAH/APC/OUTPUT/TF-SPECIFIC/GATA6.B/extension100/GATA6.B.ABC.permutationResults.rds" ,                
# "/media/carnold/s3/PAH/APC/OUTPUT/TF-SPECIFIC/GATA6.B/extension100/GATA6.B.ABC.permutationSummary.tsv"  ,               
# "/media/carnold/s3/PAH/APC/OUTPUT/Logs_and_Benchmarks/doPermutations.GATA6.B.log"
# )

# args=c("/scratch/leyva/PWM3/output/Age/TEMP/extension50/Age.allTFData_processedForPermutations.rds" ,
# "/scratch/leyva/PWM3/output/Age/TEMP/extension50/Age.allTFUniqueData_processedForPermutations.rds" ,
# "SP2"                 ,
# "/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/SP2/extension50/SP2.Age.permutationResults.rds"  ,
# "/scratch/leyva/PWM3/output/Age/TF-SPECIFIC/SP2/extension50/SP2.Age.permutationSummary.tsv" ,
# "/scratch/leyva/PWM3/output/Age/Logs_and_Benchmarks/doPermutations.SP2.log")


if (length(args) != 6) {
  stop("Expecting 6 arguments but found ", length(args),". Exiting.")
} else {
  #par.l$file_input_TF_orig           = args[1]
  #par.l$file_input_nucContentGenome  = args[2]
  par.l$file_input_TF_all            = args[1]
  par.l$file_input_TF_unique_all     = args[2] 
  par.l$TF                           = args[3]
  par.l$file_output_permResults      = args[4]
  par.l$file_output_summary          = args[5]
  par.l$file_log                     = args[6]
}

#####################
# VERIFY PARAMETERS #
#####################

assertFileExists(par.l$file_input_TF_all, access = "r")
assertFileExists(par.l$file_input_TF_all, access = "r")
assertDirectoryExists(dirname(par.l$file_log), access = "w")

assertCharacter(par.l$TF, len = 1, min.chars = 1)

for (dirname in unique(c(dirname(par.l$file_output_permResults), dirname(par.l$file_output_summary)))) {
  
  if (!testDirectoryExists(dirname)) {
    dir.create(dirname, recursive = TRUE)
  } else {
    assertDirectoryExists(dirname, access = "w")
  }
  
}

stopifnot(length(par.l$classes_CohensD) == length(par.l$thresholds_CohensD) + 1)


######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel, removeOldLog = TRUE)
printParametersLog(par.l)


TF.motifs.all        = readRDS(par.l$file_input_TF_all)
TF.motifs.all.unique = readRDS(par.l$file_input_TF_unique_all)


TFCur = par.l$TF
perm.l = list()

allTF = TFCur

output.global.TFs = data_frame(TF              = TFCur, 
                               weighted_mean   = numeric(length(allTF)),
                               weighted_Tstat  = numeric(length(allTF)),
                               weighted_CD     = numeric(length(allTF)),
                               weighted_median = numeric(length(allTF)), 
                               TFBS            = numeric(length(allTF))
                               )

uniqueBins = unique(TF.motifs.all$CG.bins)
TFCounter = 1


perm.l[[TFCur]] = list()
nameCur = paste0(TFCur,"_summary.df") # delete
perm.l[[TFCur]][[nameCur]] = data_frame(bin        = character(), 
                                        mean       = numeric(),
                                        pval       = numeric(), 
                                        ratio_TFBS = numeric(),
                                        cohensD    = numeric(),  
                                        median     = numeric()
                                                    )

nCol = ncol(perm.l[[TFCur]][[nameCur]])

nRowsCurTF = nrow(TF.motifs.all[which(TF.motifs.all$TF == TFCur),])

TF.subsetCur.df = TF.motifs.all[TF.motifs.all$TF == TFCur,]

for (bin in uniqueBins) {
  
  bin = as.character(bin)
  flog.info(paste0("Bin ", bin))
  rowsCur = which(TF.subsetCur.df$CG.bins == bin)
  TF.subset.df = TF.subsetCur.df[rowsCur,]
  TFBS.quota.bin = nrow(TF.subset.df)
  # be careful in the binned_TF_all i use motifs without duplocated regions
  binned_TF_all = TF.motifs.all.unique[which(TF.motifs.all.unique$CG.bins == bin & TF.motifs.all.unique$TF != TFCur),] # delete
  
  nRowsCur = nrow(binned_TF_all)

  # add <= 1 because of the t-test
  if (TFBS.quota.bin  <= 1 |  nRowsCur == 0) {
    
    #flog.warn(paste0("0 rows for TF ", TFCur, " and bin ", bin))

    perm.l[[TFCur]][[nameCur]][bin, "bin"] = bin
    perm.l[[TFCur]][[nameCur]][bin, 2:nCol] = NA

  } else {
    
    l2FC_data = as.vector(binned_TF_all$l2FC)
    statistical.test = t.test(l2FC_data, TF.subset.df$l2FC)
    
    
    perm.l[[TFCur]][[nameCur]][bin,1] = bin
    
    perm.l[[TFCur]][[nameCur]][bin,"mean"] = mean(TF.subset.df$l2FC) - mean(l2FC_data)
    
    perm.l[[TFCur]][[nameCur]][bin, "pval"] = statistical.test$statistic[[1]]
    
    perm.l[[TFCur]][[nameCur]][bin, "ratio_TFBS"] = TFBS.quota.bin/nRowsCurTF
    
    perm.l[[TFCur]][[nameCur]][bin, "cohensD"] = cohensD(l2FC_data, TF.subset.df$l2FC)
    
    perm.l[[TFCur]][[nameCur]][bin, "median"] = median(TF.subset.df$l2FC, na.rm = TRUE) - median(l2FC_data, na.rm = TRUE)
    
  }
  
  
  
}
  


output.global.TFs[TFCounter, "TF"]   = TFCur

output.global.TFs[TFCounter, "TFBS"] = nRowsCurTF

output.global.TFs[TFCounter, "weighted_mean"] = weighted.mean(perm.l[[TFCur]][[nameCur]]$mean ,
                                                          perm.l[[TFCur]][[nameCur]]$ratio_TFBS, na.rm = TRUE)
output.global.TFs[TFCounter, "weighted_Tstat"] = weighted.mean(perm.l[[TFCur]][[nameCur]]$pval,
                                                           perm.l[[TFCur]][[nameCur]]$ratio_TFBS, na.rm = TRUE)
output.global.TFs[TFCounter, "weighted_CD"] = weighted.mean(perm.l[[TFCur]][[nameCur]]$cohensD,
                                                        perm.l[[TFCur]][[nameCur]]$ratio_TFBS, na.rm = TRUE)
output.global.TFs[TFCounter, "weighted_median"] = weighted.mean(perm.l[[TFCur]][[nameCur]]$median,
                                                            perm.l[[TFCur]][[nameCur]]$ratio_TFBS, na.rm = TRUE)





output.global.TFs$Cohend_factor = ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[1], par.l$classes_CohensD[1], 
                                           ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[2] , par.l$classes_CohensD[2], 
                                                  ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[3], par.l$classes_CohensD[3], par.l$classes_CohensD[4])))

output.global.TFs$Cohend_factor = factor(output.global.TFs$Cohend_factor, levels = par.l$classes_CohensD, labels = seq_len(length(par.l$classes_CohensD)))

# Save objects

saveRDS(perm.l, file = par.l$file_output_permResults)

write_tsv(output.global.TFs, path = par.l$file_output_summary, col_names = TRUE)


