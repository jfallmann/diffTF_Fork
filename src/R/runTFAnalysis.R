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
# 
#  args = c(
#   "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/PEAKS/MC.sampleMetadata.rds"        ,
#    "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/PEAKS/MC.normFacs.rds"               ,
#    "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/PEAKS/MC.peaks.rds"                   ,
#    "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/PEAKS/MC.peaks.tsv"                     ,
#    "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_131106_MONK_0322_BC2YWVACXX_L4_ATTCCT_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_GTGGCC_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_TTAGGC_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_150408_BRISCOE_0218_AC69K1ACXX_L1_CGTACG_output.filtered.subsample.overlapPeaks.bed",
#    "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN740/extension50/ZN350.MC.output.tsv"  ,
#    "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN740/extension50/ZN350.MC.summary.rds"    ,
#    "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN740/extension50/ZN350.MC.MA.realcounts.pdf"   ,
#   "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN740/extension50/ZN350.MC.VSTvsReal.peaks.pdf"   ,
#    "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN740/extension50/ZN350.MC.log2.dens.VST.pdf"     ,
#    "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/output/TF-SPECIFIC/ZN740/extension50/ZN350.MC.ECDF.pdf",
#    "ZN350"                    ,
#    "CONT_WT,PAT_WT" ,
#    "~ Treatment + conditionSummary",
#   "test.log"
# 
#  )

# args = c("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/PEAKS/MC.sampleMetadata.rds"        ,             
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/PEAKS/MC.normFacs.rds"     ,              
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/PEAKS/MC.peaks.rds"         ,          
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/PEAKS/MC.peaks.tsv"       ,             
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.h27ac_131106_MONK_0322_BC2YWVACXX_L4_ATTCCT_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_GTGGCC_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_TTAGGC_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.h27ac_150408_BRISCOE_0218_AC69K1ACXX_L1_CGTACG_output.filtered.subsample.overlapPeaks.bed"        ,            
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.output.tsv"    ,  
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.summary.rds"  ,  
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.MA.realcounts.pdf"   , 
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.meanSD.pdf"    ,
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.log2.dens.pdf"   ,  
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/AHR/extension50/AHR.MC.ECDF.pdf"      ,     
# "AHR"                ,
# "~Treatment + conditionSummary"                  ,
# "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/Logs_and_Benchmarks/analyzeTF.AHR.R.log"
# )
# 
# args = c("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/PEAKS/MC.sampleMetadata.rds"        ,             
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/PEAKS/MC.normFacs.rds"     ,              
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/PEAKS/MC.peaks.rds"         ,          
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/PEAKS/MC.peaks.tsv"       ,             
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_131106_MONK_0322_BC2YWVACXX_L4_ATTCCT_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_GTGGCC_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_TTAGGC_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_150408_BRISCOE_0218_AC69K1ACXX_L1_CGTACG_output.filtered.subsample.overlapPeaks.bed"        ,            
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.output.tsv"    ,  
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.summary.rds"  ,  
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.MA.realcounts.pdf"   , 
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.meanSD.pdf"    ,
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.log2.dens.pdf"   ,  
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.ECDF.pdf"      ,     
#          "ZN350"                ,
#          "~Treatment + conditionSummary"                  ,
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/output/Logs_and_Benchmarks/analyzeTF.ZN350.R.log")
# 
# args = c("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/PEAKS/MC.sampleMetadata.rds"        ,             
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/PEAKS/MC.normFacs.rds"     ,              
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/PEAKS/MC.peaks.rds"         ,          
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/PEAKS/MC.peaks.tsv"       ,             
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_131106_MONK_0322_BC2YWVACXX_L4_ATTCCT_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_GTGGCC_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_TTAGGC_output.filtered.subsample.overlapPeaks.bed,/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.h27ac_150408_BRISCOE_0218_AC69K1ACXX_L1_CGTACG_output.filtered.subsample.overlapPeaks.bed"        ,            
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.output.tsv"    ,  
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.summary.rds"  ,  
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.MA.realcounts.pdf"   , 
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.meanSD.pdf"    ,
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.log2.dens.pdf"   ,  
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/TF-SPECIFIC/ZN350/extension50/ZN350.MC.new.ECDF.pdf"      ,     
#          "ZN350"                ,
#          "~Treatment + conditionSummary"                  ,
#          "/g/scb2/zaugg/carnold/Projects/PWM_Ivan/example/outputOld/Logs_and_Benchmarks/analyzeTF.ZN350.R.log")
# 
# 
# part1 = "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_50.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_244_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_244_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_552_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_552_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_653_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_653_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_680_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_680_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_981_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_981_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_1125_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_1125_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_1303.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_1781.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2132_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2132_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2459_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2459_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2483.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2613_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2613_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2886_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2886_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2938_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2938_2.final.s.downsample0.001.overlapPeaks.bed,"
# 
# part2 = "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2938_3.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2938_4.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2938_5.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2977_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_2977_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3069_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3069_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3142.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3156_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3156_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3215_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3215_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3215_3.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3215_4.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3240.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3263_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3263_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3386.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3439_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3439_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3439_3.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3492_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3492_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3756.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3811.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3823_1.final.s.downsample0.001.overlapPeaks.bed,"
# 
# part3 = "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3823_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3873.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3943.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_3980.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4034.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4078_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4078_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4080.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4102_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4102_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4102_3.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4189.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4251_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4251_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4333.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4621_1.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4621_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4621_3.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4668_1.final.s.downsample0.001.overlapPeaks.bed,"
# 
# 
# part4 = "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4668_2.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4747.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4784.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4963.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_4989.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_5044.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_5048.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_5129.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_5147.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_5199.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_5204.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_5229.final.s.downsample0.001.overlapPeaks.bed,/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ATAC_5263.final.s.downsample0.001.overlapPeaks.bed"

# args = c(
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/PEAKS/CLL.sampleMetadata.rds",
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/PEAKS/CLL.normFacs.rds",
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/PEAKS/CLL.peaks.rds" ,
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/PEAKS/CLL.peaks.tsv",
#   paste0(part1,part2,part3,part4),
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.output.tsv" ,
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.summary.rds" ,
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.MA.realcounts.pdf",
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.meanSD.pdf" ,
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.log2.dens.pdf",
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/TF-SPECIFIC/HSFY1/extension100/HSFY1.CLL.ECDF.pdf",
#   "HSFY1",
#   "~ Treatment + Condition" ,
#   "/scratch/carnold/CLL/TF_act_downsampling0.001/output/Logs_and_Benchmarks/analyzeTF.HSFY1.R.log"
# )

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 14) {
  stop("Expecting 14 arguments but found ", length(args),". Exiting.")
} else {
  par.l$file_input_metadata        = args[1]
  par.l$file_input_normFacs        = args[2]
  par.l$file_input_peaks           = args[3]
  par.l$file_input_peak2           = args[4]
  par.l$file_input_peakTFOverlaps  = args[5]
  par.l$file_output_summaryAll     = args[6]
  par.l$file_output_summaryStats   = args[7]
  par.l$file_output_plot_MA        = args[8]
  par.l$file_output_plotvsReal     = args[9]
  par.l$file_output_plot           = args[10]
  par.l$file_output_plot_ecdf      = args[11]
  par.l$TF                         = args[12]
  par.l$designFormula              = args[13]
  par.l$file_log                   = args[14]
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
            dirname(par.l$file_output_plot_MA), 
            dirname(par.l$file_output_plotvsReal), 
            dirname(par.l$file_output_plot),
            dirname(par.l$file_output_plot_ecdf)
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

# Group by ID
# take only the maximum row mean of all samples, sample with biggest coverage
coverageAll_grouped.df = coverageAll.df %>%
  dplyr::group_by(ID) %>%
  dplyr::slice(which.max(mean))

TF.table.m = as.matrix(coverageAll_grouped.df[,sampleData.df$name])
colnames(TF.table.m) = sampleData.df$name
rownames(TF.table.m) = coverageAll_grouped.df$identifier


# Create formula based on user-defined design
designFormula = convertToFormula(par.l$designFormula, colnames(sampleData.df))

# create Deseq object from the TF specific data
TF.cds <- DESeqDataSetFromMatrix(countData = TF.table.m,
                                 colData = sampleData.df,
                                 design = designFormula)

# normalize with normFacs
normFacs = readRDS(par.l$file_input_normFacs)
normalizationFactors(TF.cds) <- normFacs[coverageAll_grouped.df$ID,]
# low RC, check by rowMean
TF.cds.filt = TF.cds[rowMeans(counts(TF.cds)) > 0, ]



# Deseq 2 functions
# off/on shrinkage of the log2FC - betaPrior
# with the simulations of negative binomial distribution increase the sample size

# Run the local fit first, if that throws an error try the default fit type

res_DESeq = tryCatch( {
  DESeq(TF.cds.filt,fitType = 'local')
  
}, error = function(e) {
  warning("Warning: Could not run DESeq with local fitting, retry with default fitting type...")
  DESeq(TF.cds.filt)
}
)

res_DESeq.df <- as.data.frame(DESeq2::results(res_DESeq))


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


rm_col = c("coverage.x","chr.y","annotation.y","coverage.y","identifier.y" )
order = c("TF","chr","MSS","MES","strand", "PSS","PES","annotation","ID", "identifier","baseMean", "log2FoldChange","lfcSE","stat", "pvalue","padj")

# dplyr::mutate(TF = par.l$TF) gives the following weird error message: Error: Unsupported type NILSXP for column "TF"
TFCur = par.l$TF


TF_output.df = res_DESeq.df %>%
  rownames_to_column(var = "identifier") %>%
  dplyr::full_join(peaksCur.df,by = c("identifier")) %>%
  dplyr::full_join(peaksFiltered.df, by = "ID") %>%
  dplyr::filter(!is.na(baseMean)) %>%
  dplyr::rename(annotation = annotation.x, chr = chr.x, identifier = identifier.x) %>%
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

######
# MA #
######

#pdf(par.l$file_output_plotsAll)

# MA plot for TF
pdf(par.l$file_output_plot_MA)
DESeq2::plotMA(res_DESeq)
dev.off()

#############
# VSTvsReal #
#############

# transformation for the TF
notAllZeroTF <- (rowSums(counts(res_DESeq)) > 0)
pdf(par.l$file_output_plotvsReal)
meanSdPlot(assay(res_DESeq[notAllZeroTF,]))
dev.off()


comparisonDESeq = getComparisonFromDeSeqObject(res_DESeq, par.l$designFormula)
xlabLabel = paste0(" log2 FC ", comparisonDESeq)

# density plot ## nice addition 27.04
TF_dens = ggplot() + geom_density(aes(x = peaks.df$D2_l2FC,fill = "A" ),
                                      alpha = .5, color = "black") +
  geom_density(aes(x = final.TF.df$D2_l2FC,fill = "B"),size = 1, alpha = .7) +
  xlab(xlabLabel) +
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

#plot(TF_dens)
ggsave(plot = TF_dens, filename = par.l$file_output_plot, width = 4, height = 4, useDingbats = FALSE, dpi = 600)

########
# ECDF #
########

# create ecdf plots for each TF
ECDF_TF = ggplot() + 
  stat_ecdf(aes(x = final.TF.df$D2_l2FC,colour = paste0("", par.l$TF))) +
  stat_ecdf(aes(x = peaks.df$D2_l2FC,colour = "Peaks" )) +
  xlab("Log2FC WT/KO") + 
  guides(colour = guide_legend(title = "ORIGIN"))

#plot(ECDF_TF)
ggsave(plot = ECDF_TF, filename = par.l$file_output_plot_ecdf, width = 6, height = 4, useDingbats = FALSE, dpi = 600)

#dev.off()
  

