## README: prepare Permutations, to make it more efficient
start.time  <-  Sys.time()

#########################
# LIBRARY AND FUNCTIONS #
#########################

source("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/functions.R")

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "checkmate", "tools", "methods"), verbose = TRUE)

# methods needed here because Rscript does not loads this package automatically, see http://stackoverflow.com/questions/19468506/rscript-could-not-find-function


###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"



############################################
# READ AND VALIDATE COMMAND LINE ARGUMENTS #
############################################
args <- commandArgs(trailingOnly = TRUE)

# args=c("/scratch/leyva/PWM3/output/Age/FINAL_OUTPUT/extension50/Age.allMotifs.tsv",
# "/scratch/leyva/PWM3/output/Age/TEMP/extension50/Age.motifs.coord.nucContent.bed",
# "/scratch/leyva/PWM3/output/Age/TEMP/extension50/Age.allTFData_processedForPermutations.rds",
# "/scratch/leyva/PWM3/output/Age/TEMP/extension50/Age.allTFUniqueData_processedForPermutations.rds",
# "/scratch/leyva/PWM3/output/Age/Logs_and_Benchmarks/preparePermutations.log"
# )

if (length(args) != 5) {
  stop("Expecting 5 arguments but found ", length(args),". Exiting.")
} else {
  par.l$file_input_TF_allMotives          = args[1]
  par.l$file_input_nucContentGenome  = args[2]
  par.l$file_output_allTF            = args[3]
  par.l$file_output_allTFUnique      = args[4]
  par.l$file_log                     = args[5]
}

#####################
# VERIFY PARAMETERS #
#####################

assertFileExists(par.l$file_input_TF_allMotives, access = "r")
assertFileExists(par.l$file_input_nucContentGenome, access = "r")
assertDirectoryExists(dirname(par.l$file_log), access = "w")


for (dirname in unique(c(dirname(par.l$file_output_allTF), dirname(par.l$file_output_allTFUnique)))) {
  
  if (!testDirectoryExists(dirname)) {
    dir.create(dirname, recursive = TRUE)
  } else {
    assertDirectoryExists(dirname, access = "w")
  }
  
}


######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel, removeOldLog = TRUE)
printParametersLog(par.l)



# import original dataframe with motifs 
TF.motifs.ori = read_tsv(par.l$file_input_TF_allMotives, col_names = TRUE, col_types = cols())
TF.motifs.CG  = read_tsv(par.l$file_input_nucContentGenome, col_names = TRUE, col_types = cols())

# Rename some 
colnames(TF.motifs.ori) = c("TF","chr","MSS","MES", "strand", "PSS","PES","annotation","ID","identifier","baseMean",
                            "l2FC","lfcSE","stat","pval","padj") # ,"VST_diff")

colnames(TF.motifs.CG) = c("chr","MSS","MES","strand","TF","AT","CG","A","C","G","T","N","other_nucl","length")

# create the identifier as column to merge CG and ori 
TF.motifs.CG$CG.identifier  = paste0(TF.motifs.CG$TF,":" ,TF.motifs.CG$chr ,":", TF.motifs.CG$MSS,  "-", TF.motifs.CG$MES)
TF.motifs.ori$CG.identifier = paste0(TF.motifs.ori$TF,":",TF.motifs.ori$chr,":", TF.motifs.ori$MSS, "-", TF.motifs.ori$MES)


# concatenate the data in one df 

drop.cols = c( "A","C","G","T","N","other_nucl","length","chr.y","MSS.y","MES.y","strand.y","AT","CG.identifier","TF.y")

TF.motifs.all =  TF.motifs.ori %>% 
  full_join(TF.motifs.CG, by = "CG.identifier")  %>% 
  select(-one_of(drop.cols)) %>% 
  mutate(CG.bins = cut(CG, breaks = seq(0,1,0.1), 
                       labels = paste0(seq(10,100,10),"%"), include.lowest = TRUE)) %>%
  #rename(TF = TF.x, chr = chr.x, MSS = MSS.x, MES = MES.x, strand = strand.x)
  dplyr::rename(TF = TF.x, chr = chr.x, MSS = MSS.x, MES = MES.x)

# Not needed anymore, delete
rm(TF.motifs.CG)
rm(TF.motifs.ori)

# remove duplicated TFBS from different TFs to use in the permuations 
TF.motifs.all.unique = TF.motifs.all[!duplicated(TF.motifs.all[,"identifier"]),]

saveRDS(TF.motifs.all, file = par.l$file_output_allTF)
saveRDS(TF.motifs.all.unique, file = par.l$file_output_allTFUnique)