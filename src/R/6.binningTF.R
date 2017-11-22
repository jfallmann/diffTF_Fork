start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################
library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "lsr", "ggrepel", "checkmate", "tools", "methods", "boot"), verbose = TRUE)

# methods needed here because Rscript does not loads this package automatically, see http://stackoverflow.com/questions/19468506/rscript-could-not-find-function

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} and {TF} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/6.binningTF.{TF}.R.rds")
createDebugFile(snakemake)

###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$volcanoPlot_height = 16
par.l$volcanoPlot_width  = 24
par.l$minNoDatapoints = 5
par.l$log_minlevel = "INFO"


#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(names(snakemake@input), c("", "allTFData", "allTFUniqueData"))

par.l$file_input_TF_all = snakemake@input$allTFData
assertFileExists(par.l$file_input_TF_all, access = "r")

par.l$file_input_TF_unique_all = snakemake@input$allTFUniqueData
assertFileExists(par.l$file_input_TF_unique_all, access = "r")

## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "permResults", "summary", "covResults"))

par.l$file_output_permResults       = snakemake@output$permResults
par.l$file_output_summary           = snakemake@output$summary
par.l$file_output_covarianceResults = snakemake@output$covResults


## CONFIG ##
assertList(snakemake@config, min.len = 1)

par.l$nPermutations = snakemake@config$par_general$nPermutations
assertIntegerish(par.l$nPermutations, lower = 0)

## PARAMS ##
assertList(snakemake@params, min.len = 1)
assertSubset(names(snakemake@params), c("", "nBootstraps"))

par.l$nBootstraps = as.integer(snakemake@params$nBootstraps)
assertIntegerish(par.l$nBootstraps)



## WILDCARDS ##
assertList(snakemake@wildcards, min.len = 1)
assertSubset(names(snakemake@wildcards), c("", "TF"))

par.l$TF = snakemake@wildcards$TF
assertCharacter(par.l$TF, len = 1, min.chars = 1)

## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]

assertDirectoryExists(dirname(par.l$file_log), access = "w")

allDirs = c(dirname(par.l$file_output_permResults), 
            dirname(par.l$file_output_summary),
            dirname(par.l$file_output_covarianceResults)
)

testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel, removeOldLog = TRUE)
printParametersLog(par.l)

# Function for the bootstrap
ttest <- function(x, d, all) {
  
  statistical.test = t.test(all, x[d])
  return(statistical.test$statistic[[1]])
}


#############
# LOAD DATA #
#############
TF.motifs.all.orig         = readRDS(par.l$file_input_TF_all)
TF.motifs.all.unique.orig  = readRDS(par.l$file_input_TF_unique_all)

TFCur                      = par.l$TF
perm.l                     = list()
boostrapResults.l          = list()
boostrapResults.l[[TFCur]] = list()

output.global.TFs = tribble(~permutation, ~TF, ~weighted_meanDifference, ~weighted_Tstat, ~weighted_CD, ~weighted_median, ~weighted_sd, ~TFBS, ~variance)
perm.l[[TFCur]]   = tribble(~permutation, ~bin, ~meanDifference, ~nDataAll, ~nDataBin, ~pvalue, ~Tstat, ~df, ~sd, ~ratio_TFBS, ~cohensD, ~median, ~variance)
summaryCov.df     = tribble(~permutation, ~bin1, ~bin2, ~weight1, ~weight2, ~cov)

uniqueBins = unique(TF.motifs.all.orig$CG.bins)
nBins      = length(uniqueBins)
nCol       = ncol(perm.l[[TFCur]])

################
# PERMUTATIONS #
################

nPermData = length(unique(TF.motifs.all.orig$permutation))
# Adjust the number of permutations in case less have been computed
if (par.l$nPermutations + 1 < nPermData) {
  message = "In the output objects, more permutations seem to be stored. They will be ignored and the original value of nPermutations will be used"
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
} else if (par.l$nPermutations + 1 > nPermData) {
  valueNew = nPermData - 1
  message = paste0("The value of the parameter nPermutations differs from what is saved in the output objects. The value of nPermutations will be adjusted to ", valueNew)
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  par.l$nPermutations = valueNew
}


for (permutationCur in 0:par.l$nPermutations) {
  
  permutationName = paste0("permutation", permutationCur)
  
  if (permutationCur > 0) {
    flog.info(paste0("Running permutation ", permutationCur))
  } else {
    flog.info(paste0("Running for real data ", permutationCur))
  }
  
  TF.motifs.all         = filter(TF.motifs.all.orig       , permutation == permutationCur)
  TF.motifs.all.unique  = filter(TF.motifs.all.unique.orig, permutation == permutationCur)
  
  nRowsTF = nrow(TF.motifs.all[which(TF.motifs.all$TF == TFCur),])
  
  TF.subsetCur.df = TF.motifs.all[TF.motifs.all$TF == TFCur,]
  
  boostrapResults.l[[TFCur]][[permutationName]] = list()
  
  nBinsWithData = 0
  
  for (bin in uniqueBins) {
    
    boostrapResults.l[[TFCur]][[permutationName]][[bin]] = list()

    bin = as.character(bin)
    flog.info(paste0("Bin ", bin))
    rowsCur = which(TF.subsetCur.df$CG.bins == bin)
    binned.curTF.df = TF.subsetCur.df[rowsCur,]
    nRowsBinCurTF = nrow(binned.curTF.df)
    # be careful in the binned.allTF.df i use motifs without duplicated regions
    binned.allTF.df = TF.motifs.all.unique[which(TF.motifs.all.unique$CG.bins == bin & TF.motifs.all.unique$TF != TFCur),] # delete
    
    nRowsAllTF = nrow(binned.allTF.df)
  
  
    if (nRowsBinCurTF  <= par.l$minNoDatapoints |  nRowsAllTF <= par.l$minNoDatapoints) {
      
      flog.info(paste0("  Skip bin, not enough data: ", nRowsBinCurTF, " and ", nRowsAllTF))
      # Remaining columns automaitcally set to NA
      perm.l[[TFCur]] = add_row(perm.l[[TFCur]], permutation = permutationCur, bin = bin)
                                        
      boostrapResults.l[[TFCur]][[permutationName]][[bin]] = NA
  
    } else {
      
      nBinsWithData = nBinsWithData + 1
      l2fc_allBins = binned.allTF.df$l2FC
      flog.info(paste0("  Calculating data for bin based on ", nRowsBinCurTF, " and ", nRowsAllTF, " rows"))
      
      statistical.test = t.test(binned.allTF.df$l2FC, binned.curTF.df$l2FC)
      
      flog.info(paste0("  Running bootstrap... "))
      
      if (permutationCur == 0 && par.l$nBootstraps > 1) {
        boostrapResults.l[[TFCur]][[permutationName]][[bin]] = boot(data = binned.curTF.df$l2FC, statistic = ttest, R = par.l$nBootstraps, all = binned.allTF.df$l2FC)
     
        # Get the estimated boostrap variance
        varianceCur = var(boostrapResults.l[[TFCur]][[permutationName]][[bin]]$t[,1])
        
      } else {
        boostrapResults.l[[TFCur]][[permutationName]][[bin]] = NA
        varianceCur = NA
      }
      
      perm.l[[TFCur]] = add_row(perm.l[[TFCur]],
                                permutation    = permutationCur,
                                bin            = bin,
                                meanDifference = mean(binned.curTF.df$l2FC) - mean(binned.allTF.df$l2FC),
                                nDataAll       = length(binned.allTF.df$l2FC),
                                nDataBin       = length(binned.curTF.df$l2FC),
                                pvalue         = statistical.test$p.value, 
                                Tstat          = statistical.test$statistic[[1]],
                                df             = statistical.test$parameter,
                                sd             = sd(binned.curTF.df$l2FC, na.rm = TRUE), 
                                ratio_TFBS     = nRowsBinCurTF/nRowsTF,
                                cohensD        = cohensD(binned.allTF.df$l2FC, binned.curTF.df$l2FC),  
                                median         = median(binned.curTF.df$l2FC, na.rm = TRUE) - median(binned.allTF.df$l2FC, na.rm = TRUE),
                                variance       = varianceCur
      )
    
      
    }
    
  }
  
  if (nBinsWithData == 0) {
    
    message = paste0(" Not enough data for any of the ", nBins, " bins, this TF will be skipped in subsequent steps")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  } else {
    flog.info(paste0("\nFinished calculation across bins successfully for ", nBinsWithData, " out of ", nBins, " bins"))
  }
    
  ###################################################################
  # Summarize bootstrap results and estimate covariance across bins #
  ###################################################################
  
  if (permutationCur == 0) {
    
    for (bin1 in seq_len(nBins - 1)) {
      
      for (bin2 in (bin1 + 1):nBins) {
        
        bin1C = as.character(uniqueBins[bin1])
        bin2C = as.character(uniqueBins[bin2])
        weight1 = filter(perm.l[[TFCur]], bin == bin1C, permutation == permutationCur)$ratio_TFBS
        weight2 = filter(perm.l[[TFCur]], bin == bin2C, permutation == permutationCur)$ratio_TFBS
        
        if (is.na(boostrapResults.l[[TFCur]][[permutationName]][[bin1C]]) || is.na(boostrapResults.l[[TFCur]][[permutationName]][[bin2C]])) {
          
          summaryCov.df = add_row(summaryCov.df, permutation = permutationCur, bin1 = bin1, bin2 = bin2)
          
        } else {
          cov = cov(boostrapResults.l[[TFCur]][[permutationName]][[bin1C]]$t[,1], boostrapResults.l[[TFCur]][[permutationName]][[bin2C]]$t[,1], method = "pearson")
          
          summaryCov.df = add_row(summaryCov.df, permutation = permutationCur, bin1 = bin1, bin2 = bin2, weight1 = weight1, weight2 = weight2, cov = cov)
        }
        
        
      }
      
    }
    
    # Estimate the variance
    
    # Filter for bins for which we actually have data for
    perm.filtered.df   = filter(perm.l[[TFCur]], permutation == permutationCur, !is.na(ratio_TFBS), !is.na(df))
    
    if (nrow(perm.filtered.df) > 0) {
      
      weights = perm.filtered.df$ratio_TFBS
      
      if (par.l$nBootstraps > 1) {
        varianceIndividual = perm.filtered.df$variance
      } else {
        df = perm.filtered.df$df
        varianceIndividual = df / (df - 2)
      }
      
      summaryCov.filt.df = filter(summaryCov.df, permutation == permutationCur, !is.na(weight1), !is.na(weight2), !is.na(cov))

      # See https://en.wikipedia.org/wiki/Variance#Weighted_sum_of_variables
      # Function to estimate the variance of the weighted mean
      # Original proposal by Bernd
      # see the paper for a derivation of the formula
      varianceFinal     = sum(weights^2 * varianceIndividual)     + (2 * sum(summaryCov.filt.df$weight1 * summaryCov.filt.df$weight2 * summaryCov.filt.df$cov))
      
      
    } else {
      
      message = paste0("Could not calculate variance due to missing values. Set variance to NA")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      varianceFinal = NA

    }
    
    
  
  } else {
    
    varianceFinal = NA
    perm.filtered.df   = filter(perm.l[[TFCur]], permutation == permutationCur, !is.na(ratio_TFBS), !is.na(df))
  }
  
  
 
  output.global.TFs = add_row(output.global.TFs,
                              permutation             = permutationCur,
                              TF                      = TFCur,
                              weighted_meanDifference = weighted.mean(perm.filtered.df$meanDifference, perm.filtered.df$ratio_TFBS, na.rm = TRUE),
                              weighted_Tstat          = weighted.mean(perm.filtered.df$Tstat         , perm.filtered.df$ratio_TFBS, na.rm = TRUE),
                              weighted_CD             = weighted.mean(perm.filtered.df$cohensD       , perm.filtered.df$ratio_TFBS, na.rm = TRUE),
                              weighted_median         = weighted.mean(perm.filtered.df$median        , perm.filtered.df$ratio_TFBS, na.rm = TRUE),
                              weighted_sd             = weighted.mean(perm.filtered.df$sd            , perm.filtered.df$ratio_TFBS, na.rm = TRUE),
                              TFBS                    = nRowsTF,
                              variance                = varianceFinal
  )
  
  
} # end for each permutation 

# Save objects

saveRDS(perm.l, file = par.l$file_output_permResults)
saveRDS(summaryCov.df, file = par.l$file_output_covarianceResults)


write_tsv(output.global.TFs, path = par.l$file_output_summary, col_names = TRUE)

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)
