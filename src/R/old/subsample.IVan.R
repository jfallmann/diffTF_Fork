start.time  <-  Sys.time()

#nTF=640
#nTF=100
#nsubSamples = 80
#nsubSamples * 7 + nsubSamples * 3 * nTF  + nsubSamples * 2  + nsubSamples * 1  + nsubSamples * 2  + nsubSamples * 2  + nsubSamples * 2 * nTF  + nsubSamples * 2
# 5776 vs 577600
# reduce the number of files for the analyzeTF part and for the permutations part!


#########################
# LIBRARY AND FUNCTIONS #
#########################

source("/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/functions.R")

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
packages = c("tidyverse", "futile.logger", "DESeq2", "vsn", "modeest", "csaw", "checkmate", "limma", "tools", "methods", "lsr", "ggrepel", "jsonlite", "foreach", "doParallel", "purrr", "reshape2")
checkAndLoadPackages(packages, verbose = TRUE)


###################
#### PARAMETERS ###
###################

par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"
par.l$signThresholdMAPlot = 0.1
par.l$FDR_threshold = 0.05
par.l$volcanoPlot_height = 16
par.l$volcanoPlot_width  = 24
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.2, 0.5, 0.8)
par.l$FDR_threshold = c(0.01, 0.05,0.1,0.2, 0.5, 0.7)
par.l$volcanoPlot_dpi  = 600
par.l$excludeSmallCohensPrinting = FALSE

# Generate subsamples of various sizes
par.l$stepsize = 5
par.l$nRepetitionsPerStepMax = 10


##############
# PARAMETERS #
##############

par.l$rootFolder                 = "/scratch/carnold/CLL/TF_act"

par.l$file_input_config          = paste0(par.l$rootFolder, "/input/CLL.config.json")

par.l$dir_peaks                  = paste0(par.l$rootFolder, "/output/PEAKS")
par.l$pattern_peaks              = ".overlapPeaks.bed"

par.l$file_output_summaryStats   = "TODO"

par.l$nCores  = 32

par.l$file_TF_list = paste0(par.l$rootFolder, "/input/names.TF.int.tsv")
par.l$file_TF_list = paste0(par.l$rootFolder, "/input/names.TF.all.tsv")

par.l$file_output_log            = paste0(par.l$rootFolder, "/output/Logs_and_Benchmarks/subsample.R.log")



#####################
# VERIFY PARAMETERS #
#####################

assertDirectoryExists(par.l$dir_peaks, access = "r")

assertCharacter(par.l$pattern_peaks, len = 1)

assertFileExists(par.l$file_TF_list, access = "r")

assertDirectoryExists(dirname(par.l$file_output_log), access = "w")

stopifnot(length(par.l$classes_CohensD) == length(par.l$thresholds_CohensD) + 1)


assertFileExists(par.l$file_input_config)
#config.l = read_json(par.l$file_input_config)

config.l = fromJSON(file(par.l$file_input_config), simplifyVector = FALSE)

assertCharacter(config.l$samples$summaryFile, len = 1, min.chars = 1)
file_sampleData = config.l$samples$summaryFile
assertFileExists(file_sampleData)

assertIntegerish(config.l$par_general$regionExtension, lower = 1, len = 1)
extensionSize = config.l$par_general$regionExtension

assertCharacter(config.l$par_general$analysisName, len = 1, min.chars = 1)
analysisName = config.l$par_general$analysisName

assertCharacter(config.l$par_general$conditionComparison, len = 1, min.chars = 1)
conditionComparison = config.l$par_general$conditionComparison
assertCharacter(conditionComparison, len = 1, min.chars = 3)

assertCharacter(config.l$additionalInputFiles$refGenome_fasta, len = 1, min.chars = 1)
refGenome = config.l$additionalInputFiles$refGenome_fasta
assertFileExists(refGenome, access = "r")


# Check if contrasts have been specified correctly
conditionsContrast = strsplit(conditionComparison, split = ",", fixed = TRUE)[[1]]
assertVector(conditionsContrast, len = 2)

assertCharacter(config.l$par_general$designContrast, len = 1, min.chars = 1)
designFormula  = config.l$par_general$designContrast


par.l$file_output_permResultsAll        = paste0(par.l$rootFolder, "/output/FINAL_OUTPUT/extension", extensionSize, "/", analysisName, ".permResultsAll.rds")
par.l$file_output_permResultsAllSummary = paste0(par.l$rootFolder, "/output/FINAL_OUTPUT/extension", extensionSize, "/", analysisName, ".permResultsAllSummary.rds")



######################
# FINAL PREPARATIONS #
startLogger(par.l$file_output_log, par.l$log_minlevel,  removeOldLog = TRUE)
printParametersLog(par.l)

allTF.df = read_tsv(par.l$file_TF_list, col_names = FALSE)
allTF = allTF.df$X1 # [1:30]

#allTF = allTF
nTF = length(allTF)


# #setup parallel backend to use n processors
cl <- makeCluster(par.l$nCores)
registerDoParallel(cl)



#####################
# CREATE SUBSAMPLES #
#####################

subsamples.l = permuteSampleTable(file_sampleData, conditionComparison, "Treatment", par.l$stepsize, par.l$nRepetitionsPerStepMax)

# Sanity check uniqueness of sample permutations
sampleSubsets = sapply(subsamples.l, function(x) { paste0(sort(unique(x$SampleID)), collapse = ",")})
stopifnot(length(unique(sampleSubsets)) == length(subsamples.l))

# Test if the factorial variable for the model contains at least two unique values for each case
batchesPerm = sapply(subsamples.l, function(x) { sort(unique(x$Treatment))})
stopifnot(length(which(sapply(batchesPerm, length) == 1)) == 0)


perm.l = list()

#subsamples.l = subsamples.l[1]

foreach(subsampleCounter = 1:length(subsamples.l),.packages = packages) %dopar% {
#for (nameSubsampleCur in names(subsamples.l)) {

  nameSubsampleCur = names(subsamples.l)[subsampleCounter]

  startLogger(paste0(par.l$file_output_log, "_", subsampleCounter), par.l$log_minlevel,  removeOldLog = TRUE)


  flog.info(paste0("Running for subsample ", nameSubsampleCur))


  sampleData.df = subsamples.l[[nameSubsampleCur]]

  assertSubset(c("bamReads", "conditionSummary"), colnames(sampleData.df))
  assertSubset(conditionsContrast, sampleData.df$conditionSummary)

  bamFileBasenames = file_path_sans_ext(basename(sampleData.df$bamReads))



  # Read and modify samples metadata
  sampleData.df = sampleData.df %>%
    mutate(name = file_path_sans_ext(basename(sampleData.df$bamReads))) %>%
    filter(conditionSummary %in% conditionsContrast) %>%
    mutate(conditionSummary = as.factor(conditionSummary)) %>%
    mutate(Condition = as.factor(Condition)) %>%
    mutate(Treatment = as.factor(Treatment))

  ##############################
  # ITERATE THROUGH PEAK FILES #
  ##############################

  #On this step it depends what waas the initial peaks file(how many columns)
  peaks.l = list()
  coverageAll.m = NULL

  flog.info(paste0(" Iterating over ", length(sampleData.df$name), " peak files in ", par.l$dir_peaks))

  for (nameCur in sampleData.df$name) {

    file_output_PEAKS = paste0(par.l$dir_peaks, "/", nameCur, par.l$pattern_peaks)
    flog.info(paste0("  Parsing file ", file_output_PEAKS))
    assertFileExists(file_output_PEAKS)

    peaks.l[[nameCur]] =  read_tsv(file_output_PEAKS, col_names = c("chr", "PSS", "PES", "annotation", "ID", "coverage"), col_types = cols())
    peaks.l[[nameCur]]$identifier = paste0(peaks.l[[nameCur]]$chr,":", peaks.l[[nameCur]]$PSS,"-", peaks.l[[nameCur]]$PES)

    # Filter and retain only unique identifiers
    peaks.filtered.df = distinct(peaks.l[[nameCur]], identifier, .keep_all = TRUE)

    nRowsFiltered = nrow(peaks.l[[nameCur]]) - nrow(peaks.filtered.df)
    if (par.l$verbose & nRowsFiltered  > 0) flog.info(paste0("   Filtered ", nRowsFiltered, " non-unique positions out of ", nrow(peaks.l[[nameCur]]), " from peaks table."))

    peaks.l[[nameCur]] = peaks.filtered.df

    # concatenate results from COV from each iteration
    coverageAll.m = cbind(coverageAll.m, peaks.l[[nameCur]]$coverage)
  }

  # TODO: Should contain all of peaks.l[[nameCur]] and not only the last one


  ## transform as matrix data frame with counts
  coverageAll.m = as.matrix(coverageAll.m)
  colnames(coverageAll.m) = sampleData.df$name
  rownames(coverageAll.m) = peaks.l[[1]]$identifier # Take the first element as nameCurresentative, they are all identical anyway

  #############
  # RUN DESEQ #
  #############

  designFormula = convertToFormula(designFormula, colnames(sampleData.df))

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


  cds.peaks.filt = cds.peaks[rowMeans(counts(cds.peaks)) > 0, ]
  #cds.peaks.filt$conditionSummary = factor(cds.peaks.filt$conditionSummary, levels = unique(sampleData.df$conditionSummary))

  cds.peaks <- DESeq(cds.peaks.filt, fitType = 'local', quiet = TRUE)

  # Save the comparison that DeSeq made for later scripts
  comparisonDESeq = getComparisonFromDeSeqObject(cds.peaks, designFormula)


  cds.peaks.df <- as.data.frame(DESeq2::results(cds.peaks))

  final.peaks.df = data_frame("position"    = rownames(cds.peaks.df),
                              "D2_baseMean" = cds.peaks.df$baseMean,
                              "D2_l2FC"     = cds.peaks.df$log2FoldChange,
                              "D2_ldcSE"    = cds.peaks.df$lfcSE,
                              "D2_stat"     = cds.peaks.df$stat,
                              "D2_pval"     =  cds.peaks.df$pvalue,
                              "D2_padj"     =  cds.peaks.df$padj #,
                              #"vst_diff"    =  peaks.df_row$diff
  )

  #########
  # PLOTS #
  #########

  filename = paste0(par.l$dir_peaks, "/", analysisName, ".subsample", nameSubsampleCur , ".allPlots.pdf")

  pdf(filename)

  DESeq2::plotMA(cds.peaks)

  # see how change the distr after transformation
  notAllZeroPeaks <- (rowSums(counts(cds.peaks)) > 0)

  meanSdPlot(assay(cds.peaks[notAllZeroPeaks,]))
  # meanSdPlot(assay(vsd.peaks.raw[notAllZeroPeaks,]))
  dev.off()


                                                                  ####################
                                                                  ####################
                                                                  #### P A R T  2 ####
                                                                  ####################
                                                                  ####################

  flog.info(paste0(" Iterate through all ", length(allTF), " TFs"))

  TF_motifsAll.df = c()

  TFCounter = 0
  for (TFCur in allTF) {

    TFCounter = TFCounter  + 1
    flog.info(paste0("  Run for TF ", TFCur, " (", TFCounter, " out of ", nTF, ")"))

    coverageAll.df = NULL

    fileList = paste0(par.l$rootFolder, "/output/TF-SPECIFIC/", TFCur, "/extension",
                      extensionSize, "/", TFCur, ".", analysisName, ".", bamFileBasenames, par.l$pattern_peaks)


    #flog.info(paste0("  Iterating over ", length(fileList), " TF-peak overlap files"))

    for (fileCur in fileList) {

      assertFileExists(fileCur, access = "r")

      peaksCur.df = read_tsv(fileCur, col_names = FALSE, col_types = cols())
      assertDataFrame(peaksCur.df, ncols = 7)

      colnames(peaksCur.df) = c("chr","MSS","MES","annotation","ID","strand","coverage")

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
    coverageAll_grouped.df = coverageAll.df %>%
      dplyr::group_by(ID) %>%
      dplyr::slice(which.max(mean))

    TF.table.m = as.matrix(coverageAll_grouped.df[,sampleData.df$name])
    colnames(TF.table.m) = sampleData.df$name
    rownames(TF.table.m) = coverageAll_grouped.df$identifier


    # Create formula based on user-defined design
    designFormula = convertToFormula(designFormula, colnames(sampleData.df))


    # create Deseq object from the TF specific data
    TF.cds <- DESeqDataSetFromMatrix(countData = TF.table.m,
                                     colData = sampleData.df,
                                     design = designFormula)

    # normalize with normFacs
    normalizationFactors(TF.cds) <- normFacs[as.numeric(coverageAll_grouped.df$ID),]
    # low RC, check by rowMean
    TF.cds.filt = TF.cds[rowMeans(counts(TF.cds)) > 0, ]



    # Deseq 2 functions
    # off/on shrinkage of the log2FC - betaPrior
    # with the simulations of negative binomial distribution increase the sample size
    res_DESeq <- DESeq(TF.cds.filt,fitType = 'local', quiet = TRUE)
    res_DESeq.df <- as.data.frame(DESeq2::results(res_DESeq))

    # addition  02.06
    final.TF.df = data_frame("position"    = rownames(res_DESeq.df),
                             "D2_baseMean" = res_DESeq.df$baseMean,
                             "D2_l2FC"     = res_DESeq.df$log2FoldChange,
                             "D2_ldcSE"    = res_DESeq.df$lfcSE,
                             "D2_stat"     = res_DESeq.df$stat,
                             "D2_pval"     = res_DESeq.df$pvalue,
                             "D2_padj"     = res_DESeq.df$padj
    )

    # assign final.peaks.df to the peaks.df and filter away NAs at the p.adjust

    assertSubset(rownames(res_DESeq.df), peaksCur.df$identifier)

    rm_col = c("coverage.x","chr.y","annotation.y","coverage.y","identifier.y" )
    order = c("TF","chr","MSS","MES","strand", "PSS","PES","annotation","ID", "identifier","baseMean", "log2FoldChange","lfcSE","stat", "pvalue","padj")



    TF_output.df = res_DESeq.df %>%
      rownames_to_column(var = "identifier") %>%
      dplyr::full_join(peaksCur.df,by = c("identifier")) %>%
      dplyr::full_join(peaks.filtered.df, by = "ID") %>%
      dplyr::filter(!is.na(baseMean)) %>%
      dplyr::rename(annotation = annotation.x, chr = chr.x, identifier = identifier.x) %>%
      dplyr::select(-one_of(rm_col)) %>%
      dplyr::mutate(TF = TFCur) %>%
      dplyr::select(one_of(order)) %>%
      dplyr::arrange(chr)


    # Concetanate all motifs from all TF
    if (is.null(nrow(TF_motifsAll.df))) {
      TF_motifsAll.df = TF_output.df
    } else {
      TF_motifsAll.df = rbind(TF_motifsAll.df, TF_output.df)
    }



    # d) Comparisons between peaks and binding sites


    modeNum     = mlv(round(final.TF.df$D2_l2FC, 2), method = "mfv", na.rm = TRUE)
    Ttest       = t.test(final.TF.df$D2_l2FC, final.peaks.df$D2_l2FC)

    output.df = data_frame(TF              = TFCur,
                           Pos_l2FC    = nrow(final.TF.df[final.TF.df$D2_l2FC > 0,]) / nrow(final.TF.df),
                           Mean_l2FC   = mean(final.TF.df$D2_l2FC, na.rm = TRUE),
                           Median_l2FC = median(final.TF.df$D2_l2FC, na.rm = TRUE),
                           Mode_l2FC   = modeNum[[1]],
                           Ttest_pval  = Ttest$p.value,
                           Modeskewness    = modeNum[[2]],
                           T_statistic     = Ttest$statistic[[1]],
                           TFBS_num        = nrow(final.TF.df)
    )

    # needed in the volcano plot part, which is skipepd for now
    # saveRDS(output.df,file = par.l$file_output_summaryStats)

    ############
    ############
    ## GRAPHS ##
    ############
    ############


    filename = paste0(par.l$rootFolder, "/output/TF-SPECIFIC/", TFCur, "/extension",
                      extensionSize, "/", TFCur, ".", analysisName, ".subsample", nameSubsampleCur , ".allPlots.pdf")
    pdf(filename)

    # MA plot for TF
    DESeq2::plotMA(res_DESeq, main = "MA plot", alpha = par.l$signThresholdMAPlot)

    # transformation for the TF
    notAllZeroTF <- (rowSums(counts(res_DESeq)) > 0)

    title = "not supported for meanSdPlot"
    meanSdPlot(assay(res_DESeq[notAllZeroTF,]))


    xlabLabel = paste0(" log2 FC ", comparisonDESeq)

    # density plot ## nice addition 27.04
    TF_dens = ggplot() + geom_density(aes(x = final.peaks.df$D2_l2FC,fill = "A" ),
                                      alpha = .5, color = "black") +
      geom_density(aes(x = final.TF.df$D2_l2FC,fill = "B"),size = 1, alpha = .7) +
      xlab(xlabLabel) +
      ggtitle("TF density") +
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
      scale_fill_manual(values = c("A" = "grey50" , "B" = "blue"), labels = c("PEAKS", TFCur))

    plot(TF_dens)



    # create ecdf plots for each TF
    ECDF_TF = ggplot() +
      stat_ecdf(aes(x = final.TF.df$D2_l2FC,colour = paste0("", TFCur))) +
      stat_ecdf(aes(x = final.peaks.df$D2_l2FC,colour = "Peaks" )) +
      ggtitle("ECDF") +
      xlab("Log2FC WT/KO") +
      guides(colour = guide_legend(title = "ORIGIN"))

    plot(ECDF_TF)

    dev.off()


  } # end for all TF

  flog.info(paste0("Finished for all TF. Now prepare and do permutations. "))


                                                                          ####################
                                                                          ####################
                                                                          #### P A R T  3 ####
                                                                          ####################
                                                                          ####################

  file_allMotifs = paste0(par.l$rootFolder, "/output/FINAL_OUTPUT/extension", extensionSize, "/", analysisName, ".subsample", nameSubsampleCur , ".allMotifs.tsv")
  write_tsv(TF_motifsAll.df, path = file_allMotifs, col_names = TRUE)

  # only used for bedtools nuc
  file_allMotifsTemp = paste0(par.l$rootFolder, "/output/TEMP/extension", extensionSize, "/", analysisName, ".subsample", nameSubsampleCur , ".allMotifs_noHeader.tsv")


  file_nucContent = paste0(par.l$rootFolder, "/output/TEMP/extension", extensionSize, "/", analysisName, ".subsample", nameSubsampleCur, ".motifs.coord.nucContent.bed")



  commandCall = paste("cat", file_allMotifs, '| awk \'{OFS=\"\\t\"};NR>1{print $2,$3,$4,$5,$1}\' >', file_allMotifsTemp)
  system(commandCall)

  commandCall = paste("bedtools nuc -fi", refGenome, "-bed", file_allMotifsTemp, " >", file_nucContent)
  system(commandCall)


                                                                          ####################
                                                                          ####################
                                                                          #### P A R T  4 ####
                                                                          ####################
                                                                          ####################

   flog.info(paste0(" Permute and summarize"))

    # import original dataframe with motifs
    TF.motifs.ori = read_tsv(file_allMotifs  , col_names = TRUE, col_types = cols())
    TF.motifs.CG  = read_tsv(file_nucContent , col_names = TRUE, col_types = cols())

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
      dplyr::rename(TF = TF.x, chr = chr.x, MSS = MSS.x, MES = MES.x)

    # Not needed anymore, delete
    rm(TF.motifs.CG)
    rm(TF.motifs.ori)

    # remove duplicated TFBS from different TFs to use in the permuations
    TF.motifs.all.unique = TF.motifs.all[!duplicated(TF.motifs.all[,"identifier"]),]

    perm.l[[nameSubsampleCur]] = list()

    permSummaryAll.df = c()

    for (TFCur in allTF) {



      output.global.TFs = data_frame(TF              = TFCur,
                                     weighted_mean   = numeric(length(TFCur)),
                                     weighted_Tstat  = numeric(length(TFCur)),
                                     weighted_CD     = numeric(length(TFCur)),
                                     weighted_median = numeric(length(TFCur)),
                                     TFBS            = numeric(length(TFCur))
      )

      uniqueBins = unique(TF.motifs.all$CG.bins)
      TFCounter = 1


      perm.l[[nameSubsampleCur]][[TFCur]] = list()
      nameCur = paste0(TFCur,"_summary.df") # delete
      perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]] = data_frame(bin        = character(),
                                              mean       = numeric(),
                                              pval       = numeric(),
                                              ratio_TFBS = numeric(),
                                              cohensD    = numeric(),
                                              median     = numeric()
      )

      nCol = ncol(perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]])

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
       # if (TFBS.quota.bin  <= 1 |  nRowsCur == 0) {
        if (TFBS.quota.bin  < 2 |  nRowsCur < 2) { # at least two observations needed each

          #flog.warn(paste0("0 rows for TF ", TFCur, " and bin ", bin))

          perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]][bin, "bin"] = bin
          perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]][bin, 2:nCol] = NA

        } else {

          l2FC_data = as.vector(binned_TF_all$l2FC)
          statistical.test = t.test(l2FC_data, TF.subset.df$l2FC)


          perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]][bin,1] = bin

          perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]][bin,"mean"] = mean(TF.subset.df$l2FC) - mean(l2FC_data)

          perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]][bin, "pval"] = statistical.test$statistic[[1]]

          perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]][bin, "ratio_TFBS"] = TFBS.quota.bin/nRowsCurTF

          perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]][bin, "cohensD"] = cohensD(l2FC_data, TF.subset.df$l2FC)

          perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]][bin, "median"] = median(TF.subset.df$l2FC, na.rm = TRUE) - median(l2FC_data, na.rm = TRUE)

        }

      }


      output.global.TFs[TFCounter, "TF"]   = TFCur

      output.global.TFs[TFCounter, "TFBS"] = nRowsCurTF

      output.global.TFs[TFCounter, "weighted_mean"] = weighted.mean(perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]]$mean ,
                                                                    perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]]$ratio_TFBS, na.rm = TRUE)
      output.global.TFs[TFCounter, "weighted_Tstat"] = weighted.mean(perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]]$pval,
                                                                     perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]]$ratio_TFBS, na.rm = TRUE)
      output.global.TFs[TFCounter, "weighted_CD"] = weighted.mean(perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]]$cohensD,
                                                                  perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]]$ratio_TFBS, na.rm = TRUE)
      output.global.TFs[TFCounter, "weighted_median"] = weighted.mean(perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]]$median,
                                                                      perm.l[[nameSubsampleCur]][[TFCur]][[nameCur]]$ratio_TFBS, na.rm = TRUE)


      output.global.TFs$Cohend_factor = ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[1], par.l$classes_CohensD[1],
                                               ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[2] , par.l$classes_CohensD[2],
                                                      ifelse(output.global.TFs$weighted_CD < par.l$thresholds_CohensD[3], par.l$classes_CohensD[3], par.l$classes_CohensD[4])))

      output.global.TFs$Cohend_factor = factor(output.global.TFs$Cohend_factor, levels = par.l$classes_CohensD, labels = seq_len(length(par.l$classes_CohensD)))

      # Concetanate all motifs from all TF
      if (is.null(nrow(permSummaryAll.df))) {
        permSummaryAll.df = output.global.TFs[1,]
      } else {
        permSummaryAll.df = rbind(permSummaryAll.df, output.global.TFs[1,])
      }


    } # end for (TFCur in allTF)


    #############
    # SUMMARIZE #
    #############

    permSummaryAll.df$Cohend_factor = ifelse(permSummaryAll.df$weighted_CD < par.l$thresholds_CohensD[1], par.l$classes_CohensD[1],
                                             ifelse(permSummaryAll.df$weighted_CD < par.l$thresholds_CohensD[2] , par.l$classes_CohensD[2],
                                                    ifelse(permSummaryAll.df$weighted_CD < par.l$thresholds_CohensD[3], par.l$classes_CohensD[3], par.l$classes_CohensD[4])))

    permSummaryAll.df$Cohend_factor = factor(permSummaryAll.df$Cohend_factor, levels = par.l$classes_CohensD, labels = seq_len(length(par.l$classes_CohensD)))


    # add this to save output to the rule
    file_StatsAllTFs = paste0(par.l$rootFolder, "/output/FINAL_OUTPUT/extension", extensionSize, "/", analysisName, ".subsample", nameSubsampleCur , ".StatsAllTF.tsv")
    write_tsv(permSummaryAll.df, path = file_StatsAllTFs, col_names = TRUE)

    #permSummaryAll.l[[nameSubsampleCur]] = permSummaryAll.df



    ########
    # PLOT #
    ########

    # Automatically calculate the significance thresholds
    # Reverse Ivans heuristic approach earlier

    T_stat = c()
    for (FDRCur in par.l$FDR_threshold) {
      threshold1 = FDRCur /  nTF
      threshold1 = FDRCur
      T_stat = c(T_stat, qt(threshold1/2, median(permSummaryAll.df$TFBS,  na.rm = TRUE), lower.tail = FALSE))

      permSummaryAll.df[, paste0("FDR_", FDRCur)] = qt(threshold1/2, permSummaryAll.df$TFBS, lower.tail = FALSE)
    }



    min_T_stat = min(T_stat)
    min_T_stat = 1.96
    if (par.l$excludeSmallCohensPrinting) {
      plot_thr_df = permSummaryAll.df[permSummaryAll.df$Cohend_factor != "1", ]
    } else {
      plot_thr_df = permSummaryAll.df
    }

    plot_thr_df = plot_thr_df[abs(plot_thr_df$weighted_Tstat) > min_T_stat, ]


    ## plot a distribtuion
    axis_colour  = "black"
    panel_colour = "white"


    labelsCohensD = c()
    for (i in seq_len(length(par.l$classes_CohensD))) {
      if (i != length(par.l$classes_CohensD)) {
        threshold = par.l$thresholds_CohensD[i]
        signCur = "<"
      } else {
        threshold = par.l$thresholds_CohensD[i - 1]
        signCur = ">="
      }

      labelCur = paste0(par.l$classes_CohensD[i], " (", signCur, " ", threshold, ")")
      labelsCohensD = c(labelsCohensD, labelCur)
    }

    classes_CohenD = sort(as.numeric(unique(permSummaryAll.df$Cohend_factor)))


    TF_volcano1 = ggplot() +
      geom_point(aes(x = permSummaryAll.df$weighted_mean,
                     y = abs(permSummaryAll.df$weighted_Tstat),
                     size = permSummaryAll.df$Cohend_factor))  +
      geom_hline(yintercept = T_stat, size = 0.5,
                 linetype = "longdash", color = "red") +
      annotate(geom="text", label= par.l$FDR_threshold, x=Inf, y = T_stat, vjust = 0, color = "red", size = 3, hjust = 1) +

      geom_hline(yintercept = min_T_stat, size = 0.5, linetype = "longdash", color = "darkgreen") +

      geom_label_repel(aes(x = plot_thr_df$weighted_mean,
                           y = abs(plot_thr_df$weighted_Tstat),
                           size = permSummaryAll.df$Cohend_factor,
                           label = plot_thr_df$TF),
                       size = 1.75, segment.size = 0.25, box.padding = unit(0.05,"lines")) +
      ylab("Absolute T-statistic\n(with various FDR thresholds)") +
      xlab(paste0("TF 'activity' ", conditionComparison)) +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
            axis.text.y = element_text(face = "bold", color = "black", size = 10),
            axis.title.x = element_text(face = "bold", colour = "black", size = 15, margin = margin(25,0,0,0)),
            axis.title.y = element_text(face = "bold", colour = "black", size = 15, margin = margin(0,25,0,0)),
            axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
            panel.grid.major = element_line(color = panel_colour),
            panel.grid.minor = element_line(color = panel_colour),
            panel.background =  element_rect(fill = panel_colour, colour = panel_colour),
            legend.justification = "center",
            legend.background = element_rect(fill = panel_colour, colour = panel_colour),
            legend.key = element_rect(fill = panel_colour, colour = panel_colour))  +
      scale_size_manual(values = classes_CohenD, labels = labelsCohensD[classes_CohenD]) +
      #  scale_colour_gradient2(low = "blue", mid = "green" , high = "red", midpoint = 0.5) +
      labs(size = "Cohen's D")


    filename = paste0(par.l$rootFolder, "/output/FINAL_OUTPUT/extension",
                      extensionSize, "/", analysisName, ".all.volcano.l2FC.removedCGbias.AR.subsample", nameSubsampleCur , ".pdf")

    ggsave(plot = TF_volcano1, height = par.l$volcanoPlot_height, width = par.l$volcanoPlot_width, dpi = par.l$volcanoPlot_dpi, filename = filename, useDingbats = FALSE)


    message("  Finished with subsample \n")

    #permSummaryAll.df

} # end loop through all subsamples


# PLOT
# add one for loop to concatenate results from all subsamples
permSummaryAll.l = list()

for (nameSubsampleCur in names(subsamples.l)) {
  stats_filename = paste0(par.l$rootFolder, "/output/FINAL_OUTPUT/extension", extensionSize, "/", analysisName, ".subsample", nameSubsampleCur , ".StatsAllTF.tsv")
  stats.df.subsample = read_tsv(stats_filename, col_names = T)
  permSummaryAll.l[[nameSubsampleCur]] = stats.df.subsample
}
names(permSummaryAll.l) = names(subsamples.l)

saveRDS(permSummaryAll.l, file = par.l$file_output_permResultsAllSummary)
saveRDS(perm.l, file = par.l$file_output_permResultsAll)



#sel = list.select(permSummaryAll.l, "weighted_mean")

sampleSizeStr = strsplit(names(permSummaryAll.l), split = "_")

sampleSizes = as.numeric(unlist(map(sampleSizeStr, 1)))

repetitions = as.numeric(unlist(map(sampleSizeStr, 2)))


pdf("subSamplingResults.pdf")

for (TFCur in allTF) {

  rowTF = which(permSummaryAll.l[[1]]$TF == TFCur)
  stopifnot(length(rowTF) == 1)

  meanValues = map(permSummaryAll.l,"weighted_mean")

  result.df = tibble(sampleSize = sampleSizes,
                     weightedMean = unlist(map(meanValues, rowTF)))


  p <- ggplot(result.df , aes(x = sampleSize, y = weightedMean, group = sampleSize))
  p <- p + geom_boxplot() + geom_jitter(alpha = 0.1, height = 0, color = "black")
  p <- p + .getThemeForGGPlot()
  p <- p + ggtitle(TFCur)
  #p <- p + stat_summary(mapping = aes(group = sampleSize), fun.y = mean, col = "red", geom = "line")
  p <- p + ylab("TF activity: Weighted mean")
  p <- p + xlab("Sample size")
  plot(p)
}

dev.off()



end.time  <-  Sys.time()
message(" Finished execution of script. TOTAL RUNNING TIME: ", round(end.time - start.time, 1), " ", units(end.time - start.time),"\n")
