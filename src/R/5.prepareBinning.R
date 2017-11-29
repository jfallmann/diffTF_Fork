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
checkAndLoadPackages(c("tidyverse", "futile.logger", "checkmate", "tools", "methods", "BiocParallel", "rlist"), verbose = FALSE)

# methods needed here because Rscript does not loads this package automatically, see http://stackoverflow.com/questions/19468506/rscript-could-not-find-function

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/5.prepareBinning.R.rds")

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
assertSubset(names(snakemake@input), c("", "nucContent", "motifes"))

par.l$files_input_nucContentGenome  = snakemake@input$nucContent
for (fileCur in par.l$files_input_nucContentGenome) {
  assertFileExists(fileCur , access = "r")
}


par.l$file_input_TF_allMotives = snakemake@input$motifes
assertFileExists(par.l$file_input_TF_allMotives, access = "r")

## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "allTFData", "allTFUniqueData"))

par.l$file_output_allTF       = snakemake@output$allTFData
par.l$file_output_allTFUnique = snakemake@output$allTFUniqueData


## CONFIG ##
assertList(snakemake@config, min.len = 1)

par.l$nPermutations = snakemake@config$par_general$nPermutations
assertIntegerish(par.l$nPermutations, lower = 0)


## PARAMS ##
assertList(snakemake@params, min.len = 1)
assertSubset(names(snakemake@params), c("", "nCGBins"))

par.l$nBins = as.integer(snakemake@params$nCGBins)
assertIntegerish(par.l$nBins, lower = 1, upper = 100)

## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]

allDirs = c(dirname(par.l$file_output_allTF), 
            dirname(par.l$file_output_allTFUnique),
            dirname(par.l$file_log)
)

testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel, removeOldLog = TRUE)
printParametersLog(par.l)


#########################
# READ ALL MOTIVES FILE #
#########################
# import original dataframe with motifs. Has to be read only once
TF.motifs.ori = read_tsv(par.l$file_input_TF_allMotives, col_names = TRUE, col_types = cols())

if (nrow(TF.motifs.ori) == 0) {
  message = paste0("The file ", par.l$file_input_TF_allMotives, " is empty. Something went wrong before. Make sure the previous steps succeeded.")
  flog.fatal(message)
  stop(message)
}

# TODO: eventuell löschen: chr, MSS, MES, PSS, PES

#Filter permutations in the original files that the user does not want anymore
TF.motifs.ori = filter(TF.motifs.ori, permutation <= par.l$nPermutations)

TF.motifs.ori$CG.identifier = paste0(TF.motifs.ori$TF,":",TF.motifs.ori$TFBSID)


#########################
# READ ALL NUC CG FILES #
#########################


readFile <- function(permutationCur, par.l) {
    
    fileCur = par.l$files_input_nucContentGenome[permutationCur + 1]
    
    if (permutationCur > 0) {
        flog.info(paste0("Reading and processing file ", fileCur, " from permutation ", permutationCur))
    } else {
        flog.info(paste0("Reading and processing file ", fileCur, " from original data ", permutationCur))
    }
    
    TF.motifs.CG.cur  = read_tsv(fileCur, col_names = TRUE, 
                                 col_types = list(
                                     col_skip(), # "chr"
                                     col_skip(), # "MSS"
                                     col_skip(), # "MES"
                                     col_character(), # "strand"
                                     col_character(), # "TF",
                                     col_skip(), # "AT"
                                     col_number(), # "CG"
                                     col_skip(), # "A"
                                     col_skip(), # "C"
                                     col_skip(), # "G"
                                     col_skip(), # "T"
                                     col_skip(), # "N"
                                     col_skip(), # "other_nucl"
                                     col_skip() # "length"
                                 )
    )
    if (nrow(problems(TF.motifs.CG.cur)) > 0) {
        flog.fatal(paste0("Parsing errors: "), problems(TF.motifs.CG.cur), capture = TRUE)
        stop("Error when parsing the file ", fileCur, ", see errors above")
    }
    
    
    colnames(TF.motifs.CG.cur) = c("TFBSID","TF","CG")
    
    nRowsNA = length(which(is.na(TF.motifs.CG.cur$CG)))
    if (nRowsNA > 0) {
        message <- paste0("The file ", fileCur, " contains ", nRowsNA, " rows out of ", nrow(TF.motifs.CG.cur), " with a missing value for the CG content, which most likely results from an assembly discordance between the BAM files and the specified fasta file. These regions will be removed in subsequent steps. The first 10 are printed with their first 3 columns here for debugging purposes:") 
        message = paste0(message, paste0(unlist(TF.motifs.CG.cur[1:10,"chr"]), ":", unlist(TF.motifs.CG.cur[1:10,"MSS"]), "-", unlist(TF.motifs.CG.cur[1:10,"MES"]), collapse = ", "))
        checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        TF.motifs.CG.cur = TF.motifs.CG.cur[-nRowsNA,]
    }
    
    
    TF.motifs.CG.cur = TF.motifs.CG.cur %>%
        mutate(CG.identifier = paste0(TF,":" ,TFBSID), permutation = permutationCur)  %>%
        select(-one_of("TFBSID", "TF"))
    
    return(TF.motifs.CG.cur)
    
}


nCores = snakemake@threads
results.l = .execInParallelGen(nCores, 
                               returnAsList = TRUE, listNames = paste0("permutation", 0:par.l$nPermutations), 
                               iteration = 0:par.l$nPermutations, 
                               abortIfErrorParallel = TRUE, 
                               verbose = TRUE,
                               readFile, par.l)

TF.motifs.CG = do.call("rbind", results.l)



#########
# MERGE #
#########

CGBins = seq(0,1, 1/par.l$nBins)

TF.motifs.all =  TF.motifs.ori %>% 
  full_join(TF.motifs.CG, by = c("CG.identifier", "permutation"))  %>% 
  select(-one_of("CG.identifier")) %>% 
  mutate(CG.bins = cut(CG, breaks = CGBins, labels = paste0(round(CGBins[-1] * 100,0),"%"), include.lowest = TRUE)) 


# Not needed anymore, delete
rm(TF.motifs.CG)
rm(TF.motifs.ori)

# remove duplicated TFBS from different TFs to use in the permuations 
TF.motifs.all.unique = TF.motifs.all[!duplicated(TF.motifs.all[,c("permutation", "TFBSID")]),]

# TODO: TFBSID and peakID needed?

flog.info(paste0("Found ", nrow(TF.motifs.all) - nrow(TF.motifs.all.unique), " duplicated TFBS across all TF."))

saveRDS(TF.motifs.all, file = par.l$file_output_allTF)
saveRDS(TF.motifs.all.unique, file = par.l$file_output_allTFUnique)