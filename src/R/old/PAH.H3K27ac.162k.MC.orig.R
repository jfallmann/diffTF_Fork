################
## PACKAGES ####
suppressPackageStartupMessages({
  library(DESeq)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(IDPmisc)
  library(ggrepel)
  library(foreach)
  library(doParallel)
  library(weights)
  library(raster)
  library(csaw)
  library(modeest)
  library(vsn)
  library(checkmate)
  library(limma)
})

###################
#### PARAMETERS ####

wd = "/g/scb2/zaugg/berest/Projects/PAH/PAH/data/peaks/"
wd.HOC = "/g/scb2/zaugg/berest/Projects/PAH/PAH/analysis/"
expirement = "PAH.H3K27ac.162k"
clean_peaks = "H3K27ac.162k.peaks.bed"

# Cond1 = PAT
Cond1Rep1 = "h27ac_131106_MONK_0322_BC2YWVACXX_L4_ATTCCT"
Cond1Rep2 = "h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_TTAGGC"
Cond1Rep3 = "h27ac_150408_BRISCOE_0218_AC69K1ACXX_L1_GAGTGG"
Cond1Rep4 = "h27ac_131106_MONK_0322_BC2YWVACXX_L4_TTAGGC"
Cond1Rep5 = "h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_GATCAG"
Cond1Rep6 = "h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_ATTCCT"
Cond1Rep7 = "h27ac_150408_BRISCOE_0218_AC69K1ACXX_L1_TTAGGC"
Cond1Rep8 = "h27ac_131106_MONK_0322_BC2YWVACXX_L4_CGTACG"

# Cond2 = CONT
Cond2Rep1 = "h27ac_131106_MONK_0322_BC2YWVACXX_L4_GATCAG"
Cond2Rep2 = "h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_CGTACG"
Cond2Rep3 = "h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_GTGGCC"
Cond2Rep4 = "h27ac_131114_PINKERTON_0285_BC34WKACXX_L1_GAGTGG"
Cond2Rep5 = "h27ac_150407_MARPLE_0300_BC69VUACXX_L2_GATCAG"
Cond2Rep6 = "h27ac_150408_BRISCOE_0218_AC69K1ACXX_L1_CGTACG"
Cond2Rep7 = "h27ac_150408_BRISCOE_0218_AC69K1ACXX_L1_GTGAAA"

nameOfCond1 = "PATwoM"
nameOfCond2 = "CONT"



bamwd = "/g/scb2/zaugg/reyespal/Projects/PAH/data/CHIPseq/bam/filtered/"
bam_ext = "_output.filtered.bam"

EXT = 100

HOCOMOCO = "/g/scb/zaugg/berest/ATACSeq/Analysis/Ackermann/HOCOMOCO.human/PWMscan/"
nCores = 50
batch1 = "1batch"
batch2 = "2batch"

result = "/g/scb2/zaugg/berest/Projects/PAH/PAH/result/"

extenstion = "PC"
#####################
### PEAKS section ###


# give the start for the timer to calculate how long the script take
start.time  <-  Sys.time()
## Normalisation step
setwd(wd)
# assign names and condtitions to the vector! This is where u change when u want to do controls!

cond1 = c(Cond1Rep1,Cond1Rep2,Cond1Rep3,Cond1Rep4,Cond1Rep5,Cond1Rep6,Cond1Rep7,Cond1Rep8)
cond2 = c(Cond2Rep1,Cond2Rep2,Cond2Rep3,Cond2Rep4,Cond2Rep5,Cond2Rep6,Cond2Rep7)
names = c(cond1,cond2)
batches_seq = c(rep(batch1,2),rep(batch2,1),rep(batch1,3),rep(batch2,1),rep(batch1,5),rep(batch2,3))
# calculate amount of reads at the peaks! Need this step only at first time!

# for ( file in names ) {
#
#   system(paste0("mkdir -p ",wd, "H3K27ac.140k/", file), intern = T)
#
# }
#
# temp_com = paste0("ls  ",wd,"H3K27ac.140k/bins/27ac* | awk -F\"/\" \'{print $NF}\'")
# temp <- system(temp_com,intern = TRUE)
# # ## parallelize for loop !!
# # #setup parallel backend to use n processors
#
# cl <- makeCluster(nCores)
# registerDoParallel(cl)
# # # # main for lopo
# foreach_packages = c('dplyr','ggplot2','grid','gridExtra','IDPmisc','ggrepel','raster','weights','modeest',
#                    'vsn','DESeq','DESeq2', 'checkmate', 'limma')
# foreach(sam = 1:length(names), .combine = c) %:% foreach(var = 1:length(temp),.packages = foreach_packages, .combine = c) %dopar% {
#
#     #for (var in ret){
#   tempfile = temp[[var]]
#   rep = names[[sam]]
#
#   temp_bed = paste0(wd,"/H3K27ac.140k/",rep,"/",tempfile,".PEAKS.bed")
#   temp_command = paste0("/g/software/bin/multiBamCov -bams ",bamwd,rep,bam_ext," -bed ",
#                   wd,"H3K27ac.140k/bins/",tempfile," > ", temp_bed)
#
#   system(temp_command, intern = T)
#
#
# }

#  for ( rep in names) {




#  PEAKSbed = paste0(wd,rep,".",expirement,".PEAKS.bed")
#  if (file.exists(PEAKSbed)) {
#    next
#  } else {
#    system(paste0("find ",wd, "H3K27ac.140k/",rep, "/ -type f -name \'27ac*\' -exec cat {} + >> ",
#                  wd,"pre.bed && sort -k1,1 -k2,2n ",wd,"pre.bed > ",
#                  PEAKSbed, "&& rm ",wd,"pre.bed" ),intern = TRUE)
#  }
#  assertFileExists(PEAKSbed)
# }




for (rep in names) {
  PEAKSbed = paste0(wd,rep,".",expirement,".PEAKS.bed")
  CovPEAKS = paste0("/g/software/bin/multiBamCov -bams ",bamwd,rep,bam_ext," -bed ",wd,clean_peaks," > ",PEAKSbed)
  if (file.exists(PEAKSbed)) {
    next
  } else {
    system(CovPEAKS,intern = TRUE)
  }
  assertFileExists(PEAKSbed)
}

#extract the data from the files and transform it to data frames
#On this step it depends what waas the initial peaks file(how many columns)
peaks_df = list()
peaks.table = NULL
for (rep in names) {

  PEAKSbed = paste0(wd,rep,".",expirement,".PEAKS.bed")
  peaks_df[[rep]] = read.table(PEAKSbed)
  colnames(peaks_df[[rep]]) = c("chr","PSS","PES","ORI","ID","COV")
  peaks_df[[rep]]$identifier = factor(paste0(peaks_df[[rep]]$chr,":",
                                             peaks_df[[rep]]$PSS,"-",
                                             peaks_df[[rep]]$PES),
                                      levels = unique(paste0(peaks_df[[rep]]$chr,
                                                             ":",peaks_df[[rep]]$PSS,
                                                             "-",peaks_df[[rep]]$PES)))
  # get only unique identifiers
  peaks_df[[rep]] = peaks_df[[rep]] %>%
    distinct(identifier, .keep_all = TRUE)
  # concatenate results from COV from each iteration
  peaks.table = cbind(peaks.table,peaks_df[[rep]]$COV)
}

## transform as matrix data frame with counts
peaks.table = as.matrix(peaks.table)
colnames(peaks.table) = names
rownames(peaks.table) = peaks_df[[Cond1Rep1]]$identifier
#

### create colData annotation file
coldata = data.frame(row.names = names,
                     condition = c(rep(nameOfCond1,length(cond1)),
                                   rep(nameOfCond2,length(cond2)))
)
coldata$condition = factor(x = coldata$condition,levels = c(nameOfCond1,nameOfCond2))
# for batch
coldata$batch = batches_seq
coldata$batch = factor(x = coldata$batch,levels = c(batch1,batch2))
# get Deseq object
cds.peaks <- DESeqDataSetFromMatrix(countData = peaks.table,
                                    colData = coldata,
                                    design = ~ batch + condition)
# vst data dont use data about modelling of the linear model in our case batch
# Normalize with LOESS
normFacs <- exp(normOffsets(counts(cds.peaks),
                            lib.sizes = colSums(counts(cds.peaks)),
                            type = "loess"))
rownames(normFacs) = rownames(peaks.table)
# add normalization factor
normalizationFactors(cds.peaks) <- normFacs
# filter
# 02.12.16 add rowMeans > 0.1 to check fot the peaks (Armando mail)
cds.peaks.filt = cds.peaks[rowMeans(counts(cds.peaks)) >= 0.1, ]
cds.peaks.filt$condition = factor(cds.peaks.filt$condition,
                                  levels = c(nameOfCond1,nameOfCond2))
cds.peaks.filt$batch = factor(cds.peaks.filt$batch,
                              levels = c(batch1,batch2))

# Deseq analysis
cds.peaks <- DESeq(cds.peaks.filt,fitType = 'local')
# should
res.peaks <- as.data.frame(DESeq2::results(cds.peaks))



# do the vst transfromation
vsd.peaks.raw = DESeq2::varianceStabilizingTransformation(cds.peaks, fitType = 'local',blind = FALSE)
# extract vst log2FC
df.peaks_row = as.data.frame(assay(vsd.peaks.raw))
#remove batch factors
df.peaks_row = as.data.frame(removeBatchEffect(df.peaks_row, batch = as.character(coldata$batch)))
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
system(paste0("mkdir -p ", result,expirement), intern = T)
## need to have df peaks after pipeline to do graphs
# we need to have values from Deseq2 as wells as vst values
write.table(df.peaks, file = paste0(result,expirement,"/peaks.",EXT,".",extenstion,".csv"),quote = FALSE, sep = "\t",
            dec = ".", row.names = FALSE, col.names = TRUE)
# get MA plot for real counts
pdf(paste0(result,expirement,"/MAplot.peaks.",expirement,".",extenstion,".pdf") )
DESeq2::plotMA(cds.peaks)
dev.off()
# see how change the distr after transformation
notAllZeroPeaks <- (rowSums(counts( cds.peaks)) > 0)
pdf(paste0(result,expirement,"/VSTvsReal.peaks.",expirement,".",extenstion,".pdf"))
meanSdPlot(assay( cds.peaks[notAllZeroPeaks,]))
meanSdPlot(assay( vsd.peaks.raw[notAllZeroPeaks,]))
dev.off()

######################################
###### Working for each TF ###########

setwd(HOCOMOCO)
# # get the array with all TF of interest
command = paste0("ls  ",HOCOMOCO,"*_pwmscan.bed | awk -F\"_\" \'{print $1}\' | awk -F\"/\" \'{print $NF}\'")
ret <- system(command,intern = TRUE)
# ## parallelize for loop !!
# #setup parallel backend to use n processors
cl <- makeCluster(nCores)
registerDoParallel(cl)
# # # main for lopo
foreach_packages = c('dplyr','ggplot2','grid','gridExtra','IDPmisc','ggrepel','raster','weights','modeest',
                     'vsn','DESeq','DESeq2', 'checkmate', 'limma')
foreach(var = 1:length(ret),.packages = foreach_packages) %dopar% {
  #for (var in ret){
  TF = ret[[var]]
  # working in shell to create pdf to the specific folder
  dir = paste0(wd.HOC,expirement,"/",TF,"/")
  ## counts the RC at the TFBS in the peaks
  CreateDir = paste0("mkdir -p ",dir)

  if (file.exists(dir)) {
    print(paste0(dir, " is already created!"))
  } else {
    system(CreateDir,intern = TRUE)
  }

  setwd(dir)
  # get the TFBS from the peaks
  TFBSinPeaks_raw = paste0(dir,TF,".",expirement,".peaks.bed")
  GetTFBSFromPeaks = paste0("/g/software/bin/intersectBed -a ",wd,clean_peaks," -b ",
                            HOCOMOCO ,TF,"_pwmscan.bed -wa -wb > ",TFBSinPeaks_raw)

  if (file.exists(TFBSinPeaks_raw)) {
    print(paste0(TFBSinPeaks_raw," is already created!"))
  } else {
    system(GetTFBSFromPeaks,intern = TRUE)
  }
  assertFileExists(TFBSinPeaks_raw)

  # # check the output don't have the 5th column
  TFBSinPeaks_corr = paste0(dir,TF,".",expirement,".peaks.corr.ext.bed")
  CorrectTFBSinPeaks = paste0("less ",dir,TF,".",expirement,
                              ".peaks.bed | cut -f4,5,6,7,8,11 | uniq  | awk \'{OFS=\"\t\"};{ print $3,$4-",
                              EXT,",$5+", EXT,",$1,$2,$6}\' > ", TFBSinPeaks_corr)

  if (file.exists(TFBSinPeaks_corr)) {
    print(paste0(TFBSinPeaks_corr," is already created!"))
  } else {
    system(CorrectTFBSinPeaks,intern = TRUE)
  }
  assertFileExists(TFBSinPeaks_corr)

  # # coverage for the TF motifs from the Bam file
  for (rep in names) {
    CovTFBSfile = paste0(dir,rep,".",TF,".",expirement,".PEAKS.bed")
    CovTFBSinPEAKS = paste0("/g/software/bin/multiBamCov -bams ",bamwd,rep,bam_ext," -bed ",
                            dir,TF,".",expirement,".peaks.corr.ext.bed > ", CovTFBSfile )

    if (file.exists(CovTFBSfile)) {
      print(paste0(CovTFBSfile," is already created!"))
    } else {
      system(CovTFBSinPEAKS,intern = TRUE)
    }
    assertFileExists(CovTFBSfile)
  }

  ## trhis is for only one site for TF
  setwd(dir)
  TF_df = list()
  TF.table = NULL
  TF_df1 = NULL
  TF_df2 = NULL
  for (rep in names) {
    #rep = "KO1"
    CovTFBSfile = paste0(dir,rep,".",TF,".",expirement,".PEAKS.bed")

    TF_df[[rep]] = read.table(CovTFBSfile) # ? changed for ENH
    colnames(TF_df[[rep]]) = c("chr","MSS","MES","ORI","ID","strand","COV")
    TF_df[[rep]]$identifier = factor(paste0(TF_df[[rep]]$chr,":",TF_df[[rep]]$MSS,"-",TF_df[[rep]]$MES),
                                     levels = unique(paste0(TF_df[[rep]]$chr,":",TF_df[[rep]]$MSS,"-",TF_df[[rep]]$MES)))
    # get only unique identifiers

    TF_df[[rep]] = TF_df[[rep]] %>% distinct(identifier, .keep_all = TRUE)

    TF_df2 = TF_df2 %>%
      dplyr::bind_cols(TF_df[[rep]][,"COV",drop=FALSE])

  }

  TF_df2 = TF_df2 %>%
    dplyr::bind_cols(TF_df[[Cond1Rep1]][,c("chr","MSS","MES","ORI","ID","strand","identifier")])

  TF_df2$mean = apply(TF_df2[,1:length(names)],1,mean)
  colnames(TF_df2) = c(names,"chr","MSS","MES","ORI","ID","strand","identifier","mean")

  TF_df1 = TF_df2 %>%
    dplyr::group_by(ID) %>%
    dplyr::slice(which.max(mean))

  TF.table = as.matrix(TF_df1[,names])


  colnames(TF.table) = names
  rownames(TF.table) = TF_df1$identifier

  # for each TF create coldata (parallelization)
  coldata = data.frame(row.names = names,
                       condition = c(rep(nameOfCond1,length(cond1)),
                                     rep(nameOfCond2,length(cond2)))
  )
  coldata$condition = factor(x = coldata$condition,levels = c(nameOfCond1,nameOfCond2))
  # for batch
  coldata$batch = batches_seq
  coldata$batch = factor(x = coldata$batch,levels = c(batch1,batch2))
  # create Deseq object from the TF specific data
  TF.cds <- DESeqDataSetFromMatrix(countData = TF.table,
                                   colData = coldata,
                                   design = ~ batch + condition)

  # normalize with normFacs
  normalizationFactors(TF.cds) <- normFacs[TF_df1$ID,]
  # low RC
  # check by rowMean
  TF.cds.filt = TF.cds[rowMeans(counts(TF.cds)) >= 0.1, ]

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
  # df.TF = final.df.TF[complete.cases(final.df.TF[,7]),]
  # df.TF_row = df.TF_row[complete.cases(final.df.TF[,7]),]
  # 23.02 talk with Armando and Judith decided to bring NAs back
  df.TF = final.df.TF
  df.TF_row = df.TF_row

  ################### Grapical output #######################################
    #ouput material

  TF_output = tbl_df(merge(res.TF, df.TF_row,  by = 0, all = F))
  save_col = c("identifier","baseMean","log2FoldChange","lfcSE","stat",
               "pvalue","padj","VST_diff")
  rm_col = c("COV.x","chr.y","ORI.y","COV.y","identifier.y" )
  order = c("TF","chr","MSS","MES","strand","PSS","PES","ORI","ID",
            "identifier","baseMean", "log2FoldChange","lfcSE","stat",
            "pvalue","padj","VST_diff")

  TF_output = TF_output %>%
    dplyr::rename(VST_diff = diff, identifier = Row.names) %>%
    dplyr::select(one_of(save_col)) %>%
    dplyr::arrange(identifier) %>%
    dplyr::full_join( TF_df[[1]],by = c("identifier")) %>%
    dplyr::full_join(peaks_df[[1]], by = "ID") %>%
    dplyr::filter(!is.na(baseMean)) %>%
    dplyr::rename(ORI = ORI.x, chr = chr.x, identifier = identifier.x) %>%
    dplyr::select(-one_of(rm_col)) %>%
    dplyr::mutate(TF = TF) %>%
    dplyr::select(one_of(order)) %>%
    dplyr::arrange(chr)

  write.table(TF_output, file = paste0(TF,".output.",expirement,".",EXT,".",extenstion,".csv"),quote = FALSE, sep = "\t",
              dec = ".", row.names = FALSE, col.names = TRUE)


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
  saveRDS(df.output,file = paste0(TF,".df.T.",extenstion,".output"))

  # GRAPHICS
  # MA plot for TF
  pdf(paste0("MAplot.",TF,".",extenstion,".realcounts.pdf") )
  DESeq2::plotMA(TF.cds.De)
  dev.off()

  # transformation for the TF
  notAllZeroTF <- (rowSums(counts( TF.cds.De)) > 0)
  pdf(paste0("VSTvsReal.",TF,".",extenstion,".peaks.pdf"))
  meanSdPlot(assay( TF.cds.De[notAllZeroTF,]))
  meanSdPlot(assay( vsd.TF[notAllZeroTF,]))
  dev.off()


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

  pdf(paste0("DENS.",TF,".",extenstion,".log2.dens.VST.pdf"), width = 4, height = 4)
  print(TF_dens_VST)
  dev.off()
  #

  # create ecdf plots for each TF
  ECDF_TF = ggplot() + stat_ecdf(aes(x = df.TF$D2_l2FC,colour= paste0("VST ", TF))) +
   stat_ecdf(aes(x = df.peaks$D2_l2FC,colour="VST peaks" )) +
    xlab("Log2FC WT/KO") + guides(colour=guide_legend(title="ORIGIN"))

  pdf(paste0(TF,".",expirement,".",EXT,".",extenstion,".ECDF.pdf"), width = 6, height = 4)
  print(ECDF_TF)
  dev.off()

}

# #
df.tf.T = data.frame()

for (TF in ret) {

  #
  dir = paste0(wd.HOC,expirement,"/",TF,"/")
  setwd(dir)
  StatFile = paste0(TF,".df.T.",extenstion,".output")

  if (!is.na(file.info(StatFile)$size) ) {
    stat = readRDS(StatFile)
    df.tf.T[TF,1] = stat$TF
    df.tf.T[TF,2] = stat$Pos_l2FC_VST
    df.tf.T[TF,3] = stat$Mean_l2FC_VST
    df.tf.T[TF,4] = stat$Median_l2FC_VST
    df.tf.T[TF,5] = stat$Mode_l2FC_VST
    df.tf.T[TF,6] = stat$Ttest_pval_VST
    df.tf.T[TF,7] = stat$Modeskewness
    df.tf.T[TF,8] = stat$T_statistic
    df.tf.T[TF,9] = stat$TFBS_num
  } else {
    print(paste0(StatFile, " is empty!"))
  }

}
# #
# #
# # #
colnames(df.tf.T) = c("TF_name","Pos_l2FC_VST","Mean_l2FC_VST","Median_l2FC_VST","Mode_l2FC_VST","Ttest_pval_VST",
                      "Mode_skewness","T_statistic","TFBS_num")
#
#
# #
# CGcontent = read.table(file = "/g/scb/zaugg/berest/TET2/data/HOCOMOCO.mouse/doc/CG.log", sep = "\t")
# colnames(CGcontent) = c("TF_name", "CG", "noCG")
# CGcontent = rbind(CGcontent, data.frame(TF_name = "CTCF.CG", CG = 83547, noCG = 0))
# CGcontent = rbind(CGcontent, data.frame(TF_name = "CTCF.noCG", CG = 0, noCG = 188677))
#
# CGcontent$CG_ratio = CGcontent$CG/(CGcontent$CG + CGcontent$noCG)
# CG_1 = CGcontent[order(CGcontent$TF_name),]
#
# #
# d = merge(df.tf.T, CGcontent, by = "TF_name", all.y = T, sort = F)
d = df.tf.T
#
#
#
d$Ttest_pval_VST[d$Ttest_pval_VST == 0] = 9.881313e-324
d$adj_pvalue = p.adjust(d$Ttest_pval_VST, method = "BH")

mode_peaks =   modeNum_VST = mlv(round(df.peaks$D2_l2FC, 2), method = "mfv", na.rm = T)

d$Diff_mean = d$Mean_l2FC_VST - mean(df.peaks$D2_l2FC, na.rm = T)
d$Diff_median = d$Median_l2FC_VST - median(df.peaks$D2_l2FC, na.rm = T)
d$Diff_mode = d$Mode_l2FC_VST - mode_peaks[[1]]
d$Diff_skew = d$Mode_skewness - mode_peaks[[2]]

d = na.omit(d)
#
# # add d fro mthe custom file
#d = read.table("/scratch/berest/PAH/result/PAH.NEW/PAH.map.MP.PAH.NEW.100.csv",
#               sep = "\t", header = T)
#
# d = d[which(d$Diff_mean < 0.2),]
#
# #
#
plot_thr_df = d[d$Diff_mean > 0.02 | d$Diff_mean < -0.02, ]
plot_thr_df = plot_thr_df[abs(plot_thr_df$T_statistic) > 1.9603, ]
#
TF_volcano_VST = ggplot() +
  geom_point(aes(x = d$Diff_mean,
                 y = abs(d$T_statistic),
                 label = d$TF_name),size = 1)   +
  geom_vline(xintercept = 0, size = 0.7,
             linetype = "longdash", color = "blue") +
  geom_hline(yintercept = 1.9603, size = 0.7,
             linetype = "longdash", color = "red") +
  geom_text_repel(aes(x = plot_thr_df$Diff_mean,
                      y = abs(plot_thr_df$T_statistic),
                      label = plot_thr_df$TF_name),size = 2.5,
                  segment.size = 0.5,box.padding = unit(0.05,"lines")) +
  ylab("Absolute T-statistic") +
  xlab(paste0("Mean(TF distr) - mean(peaks)")) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(face = "bold", color = "black", size = 20),
        axis.title.x = element_text(face = "bold", colour = "black", size = 24,margin = margin(25,0,0,0)),
        axis.title.y = element_text(face = "bold", colour = "black", size = 24,margin = margin(0,25,0,0)),
        axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.1,0.9),
        legend.justification = "center",
        legend.title = element_blank())


setwd(wd)
pdf(paste0(result,expirement,"/PAH.map.",expirement,".",EXT,".",extenstion,".VST.1TF1P.pdf"), width = 12, height = 8)
print(TF_volcano_VST)
dev.off()
# #
# # #
write.table(d, file = paste0(result,expirement,"/PAH.map.",expirement,".",EXT,".",extenstion,".csv"),quote = FALSE, sep = "\t",
            dec = ".", row.names = FALSE, col.names = TRUE)


end.time  <-  Sys.time()
message(" Finished execution using ",nCores," cores. TOTAL RUNNING TIME: ",
        round(end.time - start.time, 1), " ", units(end.time - start.time),"\n")



## addition to disect TFs of interest
# merge data frame from Deseq and VST transformation

###########
### END ###
###########
