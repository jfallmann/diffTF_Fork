# give the start for the timer to calculate how long the script take
start.time  <-  Sys.time()

###############
## PACKAGES ###
# TODO. Clean up, what is really needed?
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggrepel)
  library(gridExtra)
  library(grid)
  
  library(IDPmisc)

  library(foreach) # REMOVE MAYBE
  library(doParallel) # REMOVE MAYBE
  library(weights)  # TODO
  library(raster) # TODO
  library(csaw)
  library(modeest)  # TODO
  library(vsn)
  library(checkmate)
  library(limma)

})

###################
#### PARAMETERS ###
###################

rootDir = "/g/scb2/zaugg/berest/Projects/PAH"

wd = paste0(rootDir, "/PAH/data/peaks/")
wd.HOC = paste0(rootDir, "/PAH/analysis/")

par.l$analysisName = "PAH.H3K27ac.162k"
# TODO: Specify the pairing one wants to do. One particular pairing corresponds to one sample data file. This can be done nicer
par.l$comparisonType = "MC"


file_sampleData = "/g/scb/zaugg/reyespal/Projects/PAH/output/old/shared/sampleTable_H3K27ac.txt"
file_cleanPeaks = paste0(wd, "H3K27ac.162k.peaks.bed")

file_MA_peaks = paste0(dir_results,par.l$analysisName,"/MAplot.peaks.",par.l$analysisName,".",par.l$comparisonType,".pdf")



file_MA = paste0("MAplot.",TF,".",par.l$comparisonType,".realcounts.pdf")
file_vctVsReal = paste0("VSTvsReal.",TF,".",par.l$comparisonType,".peaks.pdf")
file_vst = paste0("DENS.",TF,".",par.l$comparisonType,".log2.dens.VST.pdf")
file_ecdf = paste0(TF,".",par.l$analysisName,".",EXT,".",par.l$comparisonType,".ECDF.pdf")
file_output1 = paste0(dir_results,par.l$analysisName,"/peaks.",EXT,".",par.l$comparisonType,".csv")

# TODO: Description (Armando)
rowMeansThreshold = 0.1

# TODO: delete
# bamwd = "/g/scb2/zaugg/reyespal/Projects/PAH/data/CHIPseq/bam/filtered/"
# bam_ext = "_output.filtered.bam"

# TODO: Not needed anymore, will be converted to Snakemake
EXT = 100

HOCOMOCO = "/g/scb/zaugg/berest/ATACSeq/Analysis/Ackermann/HOCOMOCO.human/PWMscan/"
nCores = 50 # TODO. Remove maybe
batch1 = "1batch"
batch2 = "2batch"

dir_results = paste0(rootDir, "/PAH/result/")



##########################
#### VERIFY PARAMETERS ###
##########################

assertFileExists(file_sampleData)
assertFileExists(file_cleanPeaks)
assertDirectoryExists(wd)
assertDirectoryExists(wd.HOC)


### PEAKS section #####



## Normalisation step


sampleData.df = read_delim(file_sampleData, delim = " ", col_types = 
                           cols(
                             SampleID = col_character(),
                             bamReads = col_character(),
                             Peaks = col_character(),
                             Treatment = col_character(),
                             Condition = col_factor(levels = c("PAT", "CONT")),
                             Gender = col_factor(levels = c("M", "F")),
                             Mutation = col_factor(levels = c("MUT", "WT"))
                           ))


sampleData.df$Treatment = as.factor(sampleData.df$Treatment)



#extract the data from the files and transform it to data frames
#On this step it depends what waas the initial peaks file(how many columns)
peaks.l = list()
coverageAll.m = NULL


for (rep in names) {

  file_PEAKS = paste0(wd,rep,".",par.l$analysisName,".PEAKS.bed")
  file_PEAKS =  "/g/scb2/zaugg/berest/Projects/PAH/PAH/data/peaks/h4me1_131122_TENNISON_0270_BC34YGACXX_L2_ATTCCT.PAH.ME1.PEAKS.bed"
  rep = "h4me1_131122_TENNISON_0270_BC34YGACXX_L2_ATTCCT"
  
  peaks.l[[rep]] =  read_tsv(file_PEAKS, col_names = c("chr", "start", "end", "annotation", "id", "coverage"))
  peaks.l[[rep]]$identifier = paste0(peaks.l[[rep]]$chr,":", peaks.l[[rep]]$start,"-", peaks.l[[rep]]$end)

  # Filter and retain only unique identifiers
  peaks.filtered.df = distinct(peaks.l[[rep]], identifier, .keep_all = TRUE)
  
  if (par.l$verbose) flog.info(paste0("Filtered ", nrow(peaks.l[[rep]]) - nrow(peaks.filtered.df), " non-unique positions out of ", nrow(peaks.l[[rep]]), " from peaks table."))
  
  peaks.l[[rep]] = peaks.filtered.df

    
  # concatenate results from COV from each iteration
  coverageAll.m = cbind(coverageAll.m,peaks.l[[rep]]$coverage)
}

## transform as matrix data frame with counts
coverageAll.m = as.matrix(coverageAll.m)
colnames(coverageAll.m) = names
rownames(coverageAll.m) = peaks.l[[Cond1Rep1]]$identifier # TODO
#

### create colData annotation file

# get Deseq object
cds.peaks <- DESeqDataSetFromMatrix(countData = coverageAll.m,
                                    colData = sampleData.df,
                                    design = ~ Batch + Condition)
# vst data dont use data about modelling of the linear model in our case batch
# Normalize with LOESS
normFacs <- exp(normOffsets(counts(cds.peaks),
                            lib.sizes = colSums(counts(cds.peaks)),
                            type = "loess"))
rownames(normFacs) = rownames(coverageAll.m)
# add normalization factor
normalizationFactors(cds.peaks) <- normFacs
# filter
# 02.12.16 add rowMeans > rowMeansThreshold to check fot the peaks (Armando mail)

cds.peaks.filt = cds.peaks[rowMeans(counts(cds.peaks)) >= rowMeansThreshold, ]
cds.peaks.filt$condition = factor(cds.peaks.filt$condition,
                                  levels = c(nameOfCond1,nameOfCond2))
cds.peaks.filt$batch = factor(cds.peaks.filt$batch,
                              levels = c(batch1,batch2))

# Deseq analysis
cds.peaks <- DESeq(cds.peaks.filt,fitType = 'local')
# should
res.peaks <- as.data.frame(DESeq2::results(cds.peaks))


# TODO: Maybe remove, not needed after
# do the vst transfromation
vsd.peaks.raw = DESeq2::varianceStabilizingTransformation(cds.peaks, fitType = 'local',blind = FALSE)
# extract vst log2FC
df.peaks_row = as.data.frame(assay(vsd.peaks.raw))
#remove batch factors

# TODO: Check if order is ok
df.peaks_row = as.data.frame(removeBatchEffect(df.peaks_row, batch = as.character(sampleData.df$Treatment)))
#
df.peaks_row[,paste0("mean_",nameOfCond1)] = apply(df.peaks_row[,cond1],1,mean)
df.peaks_row[,paste0("mean_",nameOfCond2)] = apply(df.peaks_row[,cond2],1,mean)
df.peaks_row[,"diff"] = round(df.peaks_row[,paste0("mean_",nameOfCond2)] - df.peaks_row[,paste0("mean_",nameOfCond1)],3)
# addition
final.df.peaks = data.frame(rownames(res.peaks), res.peaks$baseMean,res.peaks$log2FoldChange,res.peaks$lfcSE,res.peaks$stat,
                            res.peaks$pvalue, res.peaks$padj, df.peaks_row$diff)
colnames(final.df.peaks) = c("position","D2_baseMean", "D2_l2FC", "D2_ldcSE","D2_stat",
                                                         "D2_pval", "D2_padj", "vst_diff")
# assign final.df.peaks to the df.peaks and filter away NAs at the p.adjust
#df.peaks = final.df.peaks[complete.cases(final.df.peaks[,7]),]
# 23.02 talk with Armando and Judith decided to bring NAs back
df.peaks = final.df.peaks



## add everything to one data frame

## GRAPHICS
# create results folder
system(paste0("mkdir -p ", dir_results,par.l$analysisName), intern = T)
## need to have df peaks after pipeline to do graphs
# we need to have values from Deseq2 as wells as vst values


write.table(df.peaks, file = file_output1, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
# TODO: replace by write_tsv

# get MA plot for real counts

###########
# MA plot #
###########

pdf(file_MA_peaks)
DESeq2::plotMA(cds.peaks)
dev.off()

##################
# VST after plot #
##################
# see how change the distr after transformation
notAllZeroPeaks <- (rowSums(counts( cds.peaks)) > 0)
file_VSTvsReal = paste0(dir_results,par.l$analysisName,"/VSTvsReal.peaks.",par.l$analysisName,".",par.l$comparisonType,".pdf")
pdf(file_VSTvsReal)
meanSdPlot(assay( cds.peaks[notAllZeroPeaks,]))
meanSdPlot(assay( vsd.peaks.raw[notAllZeroPeaks,]))
dev.off()


runTFAnalysis <- function(par.l$analysisName, TF, extensionSize) {
  
  
  
  ######################################
  ###### Working for each TF ###########
  
  # TODO: Dependence on underscore. Change to basenamne, make better
  
  setwd(HOCOMOCO)
  # # get the array with all TF of interest
  command = paste0("ls  ",HOCOMOCO,"*_pwmscan.bed | awk -F\"_\" \'{print $1}\' | awk -F\"/\" \'{print $NF}\'")
  ret <- system(command,intern = TRUE)
  TF = ret[[var]]
   

  
    ## trhis is for only one site for TF
    setwd(dir)
    TF_df = list()
    TF.table.m = NULL
    TF_df1 = NULL
    TF_df2 = NULL
    for (rep in names) {
      #rep = "KO1"
      CovTFBSfile = paste0(dir,rep,".",TF,".",par.l$analysisName,".PEAKS.bed")
  
      TF_df[[rep]] = read.table(CovTFBSfile) # ? changed for ENH
      colnames(TF_df[[rep]]) = c("chr","MSS","MES","annotation","ID","strand","coverage")
      TF_df[[rep]]$identifier = factor(paste0(TF_df[[rep]]$chr,":",TF_df[[rep]]$MSS,"-",TF_df[[rep]]$MES),
                                       levels = unique(paste0(TF_df[[rep]]$chr,":",TF_df[[rep]]$MSS,"-",TF_df[[rep]]$MES)))
      # get only unique identifiers
  
      TF_df[[rep]] = TF_df[[rep]] %>% distinct(identifier, .keep_all = TRUE)
  
      TF_df2 = dplyr::bind_cols(TF_df2, TF_df[[rep]][,"coverage",drop=FALSE])
  
    }
  
    TF_df2 = TF_df2 %>%
      dplyr::bind_cols(TF_df[[Cond1Rep1]][,c("chr","MSS","MES","annotation","ID","strand","identifier")])
  
    TF_df2$mean = apply(TF_df2[,1:length(names)],1,mean)
    colnames(TF_df2) = c(names,"chr","MSS","MES","annotation","ID","strand","identifier","mean")
  
    TF_df1 = TF_df2 %>%
      dplyr::group_by(ID) %>%
      dplyr::slice(which.max(mean))
  
    TF.table.m = as.matrix(TF_df1[,names])
  
  
    colnames(TF.table.m) = names
    rownames(TF.table.m) = TF_df1$identifier
  
  
    # create Deseq object from the TF specific data
    TF.cds <- DESeqDataSetFromMatrix(countData = TF.table.m,
                                     colData = sampleData.df,
                                     design = ~ Batch + Condition)
  
    # normalize with normFacs
    normalizationFactors(TF.cds) <- normFacs[TF_df1$ID,]
    # low RC
    # check by rowMean
    TF.cds.filt = TF.cds[rowMeans(counts(TF.cds)) >= rowMeansThreshold, ]
  
    TF.cds.filt$condition = factor(TF.cds.filt$condition, levels = c(nameOfCond1,nameOfCond2))
    TF.cds.filt$batch = factor(TF.cds.filt$batch, levels = c(batch1,batch2))
    # Deseq 2 functions
    # off/on shrinkage of the log2FC - betaPrior
    # with the simulations of negative binomial distribution increase the sample size
    TF.cds.De <- DESeq( TF.cds.filt,fitType = 'local')
    res.TF <- as.data.frame(DESeq2::results(TF.cds.De))
    #vst transformation
    vsd.TF <- DESeq2::varianceStabilizingTransformation(TF.cds.De,fitType = 'local', blind = F)
    # extract log2Fc from vsd object
    df.TF_row = as.data.frame(assay(vsd.TF))
    # remove batch effects
    df.TF_row = as.data.frame(removeBatchEffect(df.TF_row, batch = as.character(coldata$batch)))
    #
    df.TF_row[,paste0("mean_",nameOfCond1)] = apply(df.TF_row[,cond1],1,mean)
    df.TF_row[,paste0("mean_",nameOfCond2)] = apply(df.TF_row[,cond2],1,mean)
    df.TF_row[,"diff"] = round(df.TF_row[,paste0("mean_",nameOfCond2)] - df.TF_row[,paste0("mean_",nameOfCond1)],3)
    # addition  02.06
    final.df.TF = data.frame(rownames(res.TF), res.TF$baseMean,res.TF$log2FoldChange,res.TF$lfcSE,res.TF$stat,
                                res.TF$pvalue, res.TF$padj, df.TF_row$diff)
    colnames(final.df.TF) = c("position","D2_baseMean", "D2_l2FC", "D2_ldcSE","D2_stat",
                                                             "D2_pval", "D2_padj", "vst_diff")
    # assign final.df.peaks to the df.peaks and filter away NAs at the p.adjust
    df.TF = final.df.TF
    df.TF_row = df.TF_row
  
    ################### Grapical output #######################################
    # output material
  
    TF_output = tbl_df(merge(res.TF, df.TF_row,  by = 0, all = FALSE))
    save_col = c("identifier","baseMean","log2FoldChange","lfcSE","stat",
                 "pvalue","padj","VST_diff")
    rm_col = c("COV.x","chr.y","ORI.y","COV.y","identifier.y" )
    order = c("TF","chr","MSS","MES","strand","start","end","annotation","ID",
              "identifier","baseMean", "log2FoldChange","lfcSE","stat",
              "pvalue","padj","VST_diff")
    
    # TODO check
  
    TF_output = TF_output %>%
      dplyr::rename(VST_diff = diff, identifier = Row.names) %>%
      dplyr::select(one_of(save_col)) %>%
      dplyr::arrange(identifier) %>%
      dplyr::full_join( TF_df[[1]],by = c("identifier")) %>%
      dplyr::full_join(peaks.l[[1]], by = "ID") %>%
      dplyr::filter(!is.na(baseMean)) %>%
      dplyr::rename(ORI = ORI.x, chr = chr.x, identifier = identifier.x) %>%
      dplyr::select(-one_of(rm_col)) %>%
      dplyr::mutate(TF = TF) %>%
      dplyr::select(one_of(order)) %>%
      dplyr::arrange(chr)
    
    file_output = paste0(TF,".output.",par.l$analysisName,".",extensionSize,".",par.l$comparisonType,".csv")
  
    write.table(TF_output, file = file_output, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
  
    
    # d) Comparisons between peaks and binding sites
  
    options(scipen = 111)
  
    peaks_C_VST = nrow(df.peaks[df.peaks$log2FoldChange > 0,])/nrow(df.peaks)
    df.output = data.frame(TF = character(),Pos_l2FC_VST = double(),
                           Mean_l2FC_VST = double(),Median_l2FC_VST = double(),
                           Mode_l2FC_VST = double(),Ttest_pval_VST = double(),
                           Modeskewness = double(), T_statistic = double(), TFBS_num = double(),
                           stringsAsFactors=FALSE)
    df.output[TF,1] = TF
    # VST transformation
    df.output[TF,2] = nrow(df.TF[df.TF$D2_l2FC > 0,])/nrow(df.TF)
    df.output[TF,3] = mean(df.TF$D2_l2FC, na.rm = T)
    df.output[TF,4] = median(df.TF$D2_l2FC, na.rm = T)
  
    modeNum_VST = mlv(round(df.TF$D2_l2FC, 2), method = "mfv", na.rm = T)
    df.output[TF,5] = modeNum_VST[[1]]
  
    Ttest_VST = t.test(df.TF$D2_l2FC, df.peaks$D2_l2FC)
    df.output[TF,6] = Ttest_VST$p.value
  
    df.output[TF,7] = modeNum_VST[[2]]
  
    df.output[TF,8] = Ttest_VST$statistic[[1]]
  
    df.output[TF,9] = nrow(df.TF)
    # add the median to hte output
    saveRDS(df.output,file = paste0(TF,".df.T.",par.l$comparisonType,".output"))
    
    ############
    ############
    ## GRAPHS ##
    ############
    ############
    
    ######
    # MA #
    ######

    # MA plot for TF
    pdf(file_MA)
    DESeq2::plotMA(TF.cds.De)
    dev.off()
    
    #############
    # VSTvsReal #
    #############

    # transformation for the TF
    notAllZeroTF <- (rowSums(counts(TF.cds.De)) > 0)
    pdf(file_vctVsReal)
    meanSdPlot(assay( TF.cds.De[notAllZeroTF,]))
    meanSdPlot(assay( vsd.TF[notAllZeroTF,]))
    dev.off()
    
    #######
    # VST #
    #######

    # density plot ## nice addition 27.04
    TF_dens_VST = ggplot() + geom_density(aes(x = df.peaks$D2_l2FC,fill = "A" ),
                                          alpha = .5, color = "black") +
      geom_density(aes(x = df.TF$D2_l2FC,fill = "B"),size = 1,alpha = .7) +
      xlab(paste0(" log2 of Normalised (",nameOfCond1,"_RC/",nameOfCond2,"_RC)")) +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", colour = "black", size = 10,margin = margin(25,0,0,0)),
            axis.title.y = element_text(face = "bold", colour = "black", size = 10,margin = margin(0,25,0,0)),
            axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = c(0.9,0.9),
            legend.justification = "center",
            legend.title = element_blank()) +
      scale_fill_manual(values = c("A" = "grey50" , "B" = "blue"),labels = c("PEAKS", TF))
    
    pdf(file_vst, width = 4, height = 4)
    print(TF_dens_VST)
    dev.off()
    
    
    ########
    # ECDF #
    ########
    

    pdf(file_ecdf, width = 6, height = 4)
    
    # create ecdf plots for each TF
    ECDF_TF = ggplot() + 
              stat_ecdf(aes(x = df.TF$D2_l2FC,colour = paste0("VST ", TF))) +
              stat_ecdf(aes(x = df.peaks$D2_l2FC,colour = "VST peaks" )) +
              xlab("Log2FC WT/KO") + 
              guides(colour = guide_legend(title = "ORIGIN"))
    
    print(ECDF_TF)
    dev.off()
    

  
}


