## README: prepare Permutations, to make it more efficient
start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################
library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "lsr", "ggrepel", "checkmate", "tools", "methods", "boot"), verbose = FALSE)

# methods needed here because Rscript does not loads this package automatically, see http://stackoverflow.com/questions/19468506/rscript-could-not-find-function

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} and {TF} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/binningTF.{TF}.R.rds")

createDebugFile(snakemake)

###################
#### PARAMETERS ###
###################

par.l = list()
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"
par.l$minNoDatapoints = 5

# Used for plotting
par.l$includePlots = FALSE

#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(names(snakemake@input), c("", "sampleDataR", "nucContent", "motifes"))

par.l$file_input_metadata = snakemake@input$sampleDataR
assertFileExists(par.l$file_input_metadata, access = "r")

par.l$file_input_nucContentGenome  = snakemake@input$nucContent
assertFileExists(par.l$file_input_nucContentGenome , access = "r")

par.l$files_input_TF_allMotives = snakemake@input$motifes
for (fileCur in par.l$files_input_TF_allMotives) {
  assertFileExists(fileCur, access = "r")
}


## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "permResults", "summary"))

par.l$file_output_permResults       = snakemake@output$permResults
par.l$file_output_summary           = snakemake@output$summary


## CONFIG ##
assertList(snakemake@config, min.len = 1)

par.l$nPermutations = snakemake@config$par_general$nPermutations
assertIntegerish(par.l$nPermutations, lower = 0, len = 1)

par.l$nBins = snakemake@config$par_general$nCGBins
assertIntegerish(par.l$nBins, lower = 1, upper = 100, len = 1)

par.l$nBootstraps = as.integer(snakemake@config$par_general$nBootstraps)
assertIntegerish(par.l$nBootstraps, len = 1)

## WILDCARDS ##
assertList(snakemake@wildcards, min.len = 1)
assertSubset(names(snakemake@wildcards), c("", "TF"))

par.l$TF = snakemake@wildcards$TF
assertCharacter(par.l$TF, len = 1, min.chars = 1)


## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]

allDirs = c(dirname(par.l$file_output_permResults), 
            dirname(par.l$file_output_summary),
            dirname(par.l$file_log)
)

testExistanceAndCreateDirectoriesRecursively(allDirs)


TFCur                      = par.l$TF
# TODO: perm.l with all TF instead of TFCur
perm.l                     = list()
calculateVariance = par.l$nPermutations == 0

if (calculateVariance) {
  boostrapResults.l          = list()
  boostrapResults.l[[TFCur]] = list() 
}

output.global.TFs = tribble(~permutation, ~TF, ~weighted_meanDifference, ~weighted_CD, ~TFBS, ~weighted_Tstat, ~variance)
perm.l[[TFCur]]   = tribble(~permutation, ~bin, ~meanDifference, ~nDataAll, ~nDataBin, ~ratio_TFBS, ~cohensD, ~variance, ~df, ~pvalue, ~Tstat)
summaryCov.df     = tribble(~permutation, ~bin1, ~bin2, ~weight1, ~weight2, ~cov)


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


if (calculateVariance && par.l$nBootstraps < 1000) {
  flog.warn(paste0("The value for nBootstraps is < 1000. We strongly recommend using a higher value in order to reliably estimate the statistical variance."))
}

sampleData.l = readRDS(par.l$file_input_metadata)

if (length(sampleData.l) == 0) {
  message = "Length of sampleData.l list is 0 but is has to be at least 1. Rerun the rule DiffPeaks."
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}


# Adjust the number of permutations in case less have been computed
if (par.l$nPermutations + 1 < length(sampleData.l)) {
  message = paste0("In the output objects, more permutations seem to be stored. They will be ignored and the currently specified value of nPermutations will be used")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
} else if (par.l$nPermutations + 1 > length(sampleData.l)) {
  valueNew = length(sampleData.l) - 1
  message = paste0("The value of the parameter nPermutations differs from what is saved in the output objects. The value of nPermutations will be adjusted to ", valueNew)
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  par.l$nPermutations = valueNew
  par.l$files_input_TF_allMotives = par.l$files_input_TF_allMotives[1:(par.l$nPermutations+1)]
}

####################
# READ NUC CG FILE #
####################

flog.info(paste0("Reading and processing file ", fileCur))

TF.motifs.CG  = read_tsv(par.l$file_input_nucContentGenome, col_names = TRUE, 
                         col_types = list(
                           col_skip(), # "chr"
                           col_skip(), # "MSS"
                           col_skip(), # "MES"
                           col_character(), # "TFBSID"
                           col_character(), # "TF",
                           col_skip(), # "AT"
                           col_double(), # "CG"
                           col_skip(), # "A"
                           col_skip(), # "C"
                           col_skip(), # "G"
                           col_skip(), # "T"
                           col_skip(), # "N"
                           col_skip(), # "other_nucl"
                           col_skip() # "length"
                         )
)

if (nrow(problems(TF.motifs.CG)) > 0) {
  flog.fatal(paste0("Parsing errors: "), problems(TF.motifs.CG), capture = TRUE)
  stop("Error when parsing the file ", fileCur, ", see errors above")
}


colnames(TF.motifs.CG) = c("TFBSID","TF","CG")

nRowsNA = length(which(is.na(TF.motifs.CG$CG)))
if (nRowsNA > 0) {
  message <- paste0("The file ", fileCur, " contains ", nRowsNA, " rows out of ", nrow(TF.motifs.CG), " with a missing value for the CG content, which most likely results from an assembly discordance between the BAM files and the specified fasta file. These regions will be removed in subsequent steps. The first 10 are printed with their first 3 columns here for debugging purposes:") 
  message = paste0(message, paste0(unlist(TF.motifs.CG[1:10,"chr"]), ":", unlist(TF.motifs.CG[1:10,"MSS"]), "-", unlist(TF.motifs.CG[1:10,"MES"]), collapse = ", "))
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  TF.motifs.CG = TF.motifs.CG[-nRowsNA,]
}


TF.motifs.CG = TF.motifs.CG %>%
  mutate(CG.identifier = paste0(TF,":" ,TFBSID))  %>%
 dplyr:: select(-one_of("TFBSID"))


CGBins = seq(0,1, 1/par.l$nBins)


#########################
# READ ALL MOTIVES FILE #
#########################

################
# PERMUTATIONS #
################

for (fileCur in par.l$files_input_TF_allMotives) {
    
  # Log 2 fold-changes from the particular permutation
  TF.motifs.ori  = read_tsv(fileCur, col_names = TRUE, 
                            col_types = list(
                              col_character(), # "TF",
                              col_character(), # "TFBSID"
                              col_double() # "log2FoldChange", important to have double here as col_number would not parse numbers in scientific notation correctly
                            )
  )

 
  if (nrow(problems(TF.motifs.ori)) > 0) {
    flog.fatal(paste0("Parsing errors: "), problems(TF.motifs.ori), capture = TRUE)
    stop("Error when parsing the file ", fileCur, ", see errors above")
  }
  
  if (nrow(TF.motifs.ori) == 0) {
    message = paste0("The file ", fileCur, " is empty. Something went wrong before. Make sure the previous steps succeeded.")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  if (ncol(TF.motifs.ori) != 3) {
    message = paste0("The file ", fileCur, " does not have 3 columns. Something is wrong with the number of permutations. We recommend restarting the pipeline from the DiffPeaks step.")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
 
  colnames(TF.motifs.ori) = c("TF", "TFBSID", "log2FoldChange")
  
  permutationCur = as.numeric(gsub(".*perm([0-9]+).tsv.gz", '\\1', fileCur))
  TF.motifs.ori$permutation = permutationCur
  
  if (permutationCur > 0) {
    flog.info(paste0("Running permutation ", permutationCur))
  } else {
    flog.info(paste0("Running for real data ", permutationCur))
  }
  
  
  #Filter permutations in the original files that the user does not want anymore
  TF.motifs.ori = TF.motifs.ori %>% 
    dplyr::filter(permutation <= par.l$nPermutations) %>%
    mutate(CG.identifier = paste0(TF,":",TFBSID)) %>%
    dplyr::select(-one_of("TF"))
  

  #########
  # MERGE #
  #########
  
  # TODO: full join necessary?
  TF.motifs.all =  TF.motifs.ori %>% 
    full_join(TF.motifs.CG, by = c("CG.identifier"))  %>% 
    mutate(CG.bins = cut(CG, breaks = CGBins, labels = paste0(round(CGBins[-1] * 100,0),"%"), include.lowest = TRUE))  %>%  
    dplyr::select(-one_of("CG.identifier", "CG"))
  
  
  # Not needed anymore, delete
  rm(TF.motifs.ori)
  
  # remove duplicated TFBS from different TFs to use in the permuations 
  TF.motifs.all.unique = TF.motifs.all[!duplicated(TF.motifs.all[,c("permutation", "TFBSID")]),]
  

  uniqueBins = unique(TF.motifs.all$CG.bins)
  nBins      = length(uniqueBins)
  nCol       = ncol(perm.l[[TFCur]])
  
  # TODO: TFBSID needed?
  
  flog.info(paste0(" Found ", nrow(TF.motifs.all) - nrow(TF.motifs.all.unique), " duplicated TFBS across all TF."))
  

  nRowsTF = nrow(TF.motifs.all[which(TF.motifs.all$TF == TFCur),])
  
  TF.subsetCur.df = TF.motifs.all[TF.motifs.all$TF == TFCur,]
  
  nBinsWithData = 0
  if (calculateVariance)
    boostrapResults.l[[TFCur]][[as.character(permutationCur)]] = list()
  
  if (par.l$includePlots) binnedCombined.df  = tribble(~bin, ~type, ~value)
  
  for (bin in uniqueBins) {
    
    bin = as.character(bin)
    if (calculateVariance) boostrapResults.l[[TFCur]][[as.character(permutationCur)]][[bin]] = list()
    #flog.info(paste0("Bin ", bin))
    rowsCur = which(TF.subsetCur.df$CG.bins == bin)
    binned.curTF.df = TF.subsetCur.df[rowsCur,]
    nRowsBinCurTF = nrow(binned.curTF.df)
    # be careful in the binned.allTF.df i use motifs without duplicated regions
    binned.allTF.df = dplyr::filter(TF.motifs.all.unique, CG.bins == bin, TF != TFCur, is.finite(log2FoldChange))
    
    nRowsAllTF = nrow(binned.allTF.df)

   
    if (nRowsBinCurTF  <= par.l$minNoDatapoints |  nRowsAllTF <= par.l$minNoDatapoints) {
      
      #flog.info(paste0("  Skip bin, not enough data: ", nRowsBinCurTF, " and ", nRowsAllTF))
      # Remaining columns automatically set to NA
      perm.l[[TFCur]] = add_row(perm.l[[TFCur]], permutation = permutationCur, bin = bin)
      if (calculateVariance) boostrapResults.l[[TFCur]][[as.character(permutationCur)]][[bin]] = NA

      #binnedCombined.df = add_row(binnedCombined.df, bin = bin, type = paste0(TFCur, "-only"), value = NA)
      #binnedCombined.df = add_row(binnedCombined.df, bin = bin, type = paste0("all_other"), value = NA)
    } else {
      
      nBinsWithData = nBinsWithData + 1

      if (calculateVariance) {
        
        flog.info(paste0("  Running bootstrap... "))
        
        boostrapResults.l[[TFCur]][[as.character(permutationCur)]][[bin]] = boot(data = binned.curTF.df$log2FoldChange, statistic = ttest, R = par.l$nBootstraps, all = binned.allTF.df$log2FoldChange)

        # Get the estimated boostrap variance
        varianceCur = var(boostrapResults.l[[TFCur]][[as.character(permutationCur)]][[bin]]$t[,1])
        
        statistical.test = t.test(binned.allTF.df$log2FoldChange, binned.curTF.df$log2FoldChange)
        df             = statistical.test$parameter
        pvalue         = statistical.test$p.value 
        Tstat          = statistical.test$statistic[[1]]
        
      } else {
        
        if (calculateVariance) boostrapResults.l[[TFCur]][[as.character(permutationCur)]][[bin]] = NA
        df =  pvalue = Tstat = varianceCur = NA
       
      }
      
      perm.l[[TFCur]] = add_row(perm.l[[TFCur]],
                                permutation    = permutationCur,
                                bin            = bin,
                                meanDifference = mean(binned.curTF.df$log2FoldChange) - mean(binned.allTF.df$log2FoldChange),
                                nDataAll       = length(binned.allTF.df$log2FoldChange),
                                nDataBin       = length(binned.curTF.df$log2FoldChange),
                                ratio_TFBS     = nRowsBinCurTF/nRowsTF,
                                cohensD        = cohensD(binned.allTF.df$log2FoldChange, binned.curTF.df$log2FoldChange),
                                variance       = varianceCur,
                                df             = df,
                                pvalue         = pvalue, 
                                Tstat          = Tstat
      )
      
      if (par.l$includePlots) {
          binnedCombined.df = add_row(binnedCombined.df, bin = bin, type = paste0(TFCur, "-only"), value = binned.curTF.df$log2FoldChange)
          binnedCombined.df = add_row(binnedCombined.df, bin = bin, type = paste0("all_other"), value =  binned.allTF.df$log2FoldChange)
      }
      
      
      
  
      
      
    }
    
  } # end for each bin
  
 
  if (nBinsWithData == 0) {
    
    message = paste0(" Not enough data for any of the ", nBins, " bins, this TF will be skipped in subsequent steps")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  } else {
    flog.info(paste0(" Finished calculation across bins successfully for ", nBinsWithData, " out of ", nBins, " bins"))
  }
  
  ###################################################################
  # Summarize bootstrap results and estimate covariance across bins #
  ###################################################################
  
  varianceFinal = NA
  

  # if par.l$nPermutations == 0  && par.l$nBootstraps > 1
  
  if (calculateVariance) {
    
    # Estimate variance
    
    for (bin1 in seq_len(nBins - 1)) {
      
      for (bin2 in (bin1 + 1):nBins) {
        
        bin1C = as.character(uniqueBins[bin1])
        bin2C = as.character(uniqueBins[bin2])
        weight1 = filter(perm.l[[TFCur]], bin == bin1C, permutation == permutationCur)$ratio_TFBS
        weight2 = filter(perm.l[[TFCur]], bin == bin2C, permutation == permutationCur)$ratio_TFBS
        
        if (is.na(boostrapResults.l[[TFCur]][[as.character(permutationCur)]][[bin1C]]) || is.na(boostrapResults.l[[TFCur]][[as.character(permutationCur)]][[bin2C]])) {
          
          summaryCov.df = add_row(summaryCov.df, permutation = permutationCur, bin1 = bin1, bin2 = bin2)
          
        } else {
          cov = cov(boostrapResults.l[[TFCur]][[as.character(permutationCur)]][[bin1C]]$t[,1], boostrapResults.l[[TFCur]][[as.character(permutationCur)]][[bin2C]]$t[,1], method = "pearson")
          
          summaryCov.df = add_row(summaryCov.df, permutation = permutationCur, bin1 = bin1, bin2 = bin2, weight1 = weight1, weight2 = weight2, cov = cov)
        }
        
        
      }
      
    }
    
    # Estimate the variance
    
    # Filter for bins for which we actually have data for
    perm.filtered.df   = filter(perm.l[[TFCur]], permutation == permutationCur, !is.na(ratio_TFBS))
    
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
      varianceFinal = sum(weights^2 * varianceIndividual) + (2 * sum(summaryCov.filt.df$weight1 * summaryCov.filt.df$weight2 * summaryCov.filt.df$cov))
      
      
    } else {
      
      message = paste0("Could not calculate variance due to missing values. Set variance to NA")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  
    }
  } 
  
  #############
  # SUMMARIZE #
  #############
  
  perm.filtered.df   = filter(perm.l[[TFCur]], permutation == permutationCur, !is.na(ratio_TFBS))
  
  if (calculateVariance) {
    perm.filtered.df   = filter(perm.filtered.df, !is.na(df))
  }
  
  wmd = weighted.mean(perm.filtered.df$meanDifference, perm.filtered.df$ratio_TFBS, na.rm = TRUE)
  
  output.global.TFs = add_row(output.global.TFs,
                              permutation             = permutationCur,
                              TF                      = TFCur,
                              weighted_meanDifference = wmd,
                              weighted_CD             = weighted.mean(perm.filtered.df$cohensD, perm.filtered.df$ratio_TFBS, na.rm = TRUE),
                              weighted_Tstat          = weighted.mean(perm.filtered.df$Tstat  , perm.filtered.df$ratio_TFBS, na.rm = TRUE),
                              TFBS                    = nRowsTF,
                              variance                = varianceFinal)
  

  if (par.l$includePlots) {
      xlabStr = paste0("log2 fold-change of TFBS")
      binnedCombined.df$bin = factor(binnedCombined.df$bin, levels = levels(uniqueBins))
      g1 = ggplot(binnedCombined.df, aes(value, fill = type)) + geom_density(alpha = 0.5) + 
          xlab(xlabStr) + ylab("Density") + xlim(c(-1,1)) + scale_fill_manual("TF", values=c("gray", "darkblue")) +
          geom_vline(xintercept = mean(binned.curTF.df$log2FoldChange) - mean(binned.allTF.df$log2FoldChange), linetype = "dotted") + 
          theme_bw() + facet_wrap(~ bin, nrow = 6) + theme(legend.position="top") 
      
      plot.df = filter(perm.l[[TFCur]], !is.na(meanDifference)) %>%
          mutate(binNo = as.numeric(gsub("%", "", bin)))  %>%
          dplyr::select(one_of("binNo", "meanDifference", "ratio_TFBS"))  %>%
          dplyr::rename(weight = ratio_TFBS)
      plot.new.df = reshape2::melt(plot.df, id = "binNo")
      
      g2 = ggplot(plot.new.df, aes(binNo)) + geom_line(aes(y = weight, color = "weight")) + geom_point(aes(y = weight, color = "weight")) +
                  scale_y_continuous(sec.axis = sec_axis(~./5, name = "Mean difference")) + 
                  xlab("Bin (in % of GC-content)") + ylab("Weight") + geom_hline(yintercept = wmd, linetype = "dotted") + theme_bw()
                  
                  gridExtra::grid.arrange(g1,g2,ncol=2)
  }
 
  
} # end for each permutation 


# Save objects

saveRDS(perm.l, file = par.l$file_output_permResults)

# Convert all numeric data types to character in order to prevent any scientific notation
output.global.TFs = mutate_if(output.global.TFs, is.numeric, as.character)
write_tsv(output.global.TFs, path = par.l$file_output_summary, col_names = TRUE)

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)
