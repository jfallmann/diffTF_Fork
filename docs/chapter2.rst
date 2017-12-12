.. _workflow:

Workflow
************************************************************

The workflow is illustrated by the following two Figures. First, we show a schematic of the diffTF workflow to illustrate concepts:


   .. figure:: Figures/Workflow_transparent.png
      :scale: 50 %
      :alt: Schematic of the diffTF workflow
      :align: center

      Schematic of the diffTF workflow, with input and output of the pipeline highlighted.


We now show which rules are executed by Snakemake for a specific example (see the caption of the image):


   .. figure:: Figures/DAG_transparent.png
         :scale: 70 %
         :alt: Directed acyclic graph of an example workflow
         :align: center

         Exact workflow (a so-called directed acyclic graph, or DAG) that is executed when calling Snakemake for an easy of example with two TFs (CEBPB and CTCF) for the two samples GMP.WT1 and MPP.WT1. Each node represents a rule name as defined in the Snakefile, and each arrow a dependency.

diffTF is implemented as a Snakemake pipeline. For a gentle introduction about Snakemake, see Section :ref:`workingWithPipeline`. As you can see, the workflow consists of the following steps or *rules*:

- ``checkParameterValidity``: R script that checks whether the specified peak file has the correct format, whether the provided *fasta* file and the *BAM* files are compatible, and other checks
- ``produceConsensusPeaks``:  R script that generates the consensus peaks if none are provided
- ``filterSexChromosomesAndSortPeaks``: Filters various chromosomes 8sex, unassembled ones, contigs, etc) from the peak file.
- ``sortTFBS``: Sort the TFBS lists by position
- ``resortBAM``: Sort the *BAM* file for optimized processing
- ``intersectPeaksAndBAM``: Count all reads for peak regions across all input files
- ``intersectPeaksAndTFBS``: Intersect all TFBS with peak regions to retain only TFBS in peak regions
- ``intersectTFBSAndBAM``: Count all reads from all TFBS across all input files in a TF-specific manner
- ``DESeqPeaks``: R script that performs a differential accessibility analysis for the peak regions as well as sample permutations
- ``analyzeTF``: R script that performs a TF-specific differential accessibility analysis
- ``summary1``:  R script that summarizes the previous script for all TFs
- ``concatenateMotifs``: Concatenates previous results (TFBS motives)
- ``calcNucleotideContent``: Calculates the GC content for all TFBS
- ``prepareBinning``:  R script that prepares the binning procedure
- ``binningTF``:  R script that performs the binning approach in a TF-specific manner
- ``summaryFinal``:  R script that summarizes the analysis and calculates final statistics
- ``cleanUpLogFiles``: Cleans up the ``LOGS_AND_BENCHMARKS`` directory (mostly relevant if run in cluster mode)


Input
************************************************************


Summary
==============================

As input for diffTF for your own analysis, the following data are needed:

- *BAM* file with aligned reads for each sample (see :ref:`parameter_summaryFile`)
- genome reference *fasta* that has been used to produce the *BAM* files (see :ref:`parameter_refGenome_fasta`)
- Optionally: corresponding RNA-Seq data (see :ref:`parameter_RNASeqCounts`)

In addition, the following files are need, all of which we provide already for human hg19, hg38 and mouse mm10:

- TF-specific list of TFBS (see :ref:`parameter_dir_TFBS`)
- mapping table (see :ref:`parameter_HOCOMOCO_mapping`)

.. _configurationFile:

General configuration file
==============================

To run the pipeline, a configuration file that defines various parameters of the pipeline is required.

.. note:: Please note the following important points:

  - the name of this file is irrelevant, but it must be in the right format (JSON) and it must be referenced correctly when calling Snakemake (via the ``--configfile`` parameter). We recommend naming it ``config.json``
  - neither section nor parameter names must be changed.
  - For parameters that specify a path, both absolute and relative paths are possible.  We recommend specifying an absolute path. Relative paths must be specified relative to the Snakemake working directory.
  - For parameters that specify a directory, there should be no trailing slash.

In the following, we explain all parameters in detail, organized by section names.


SECTION ``par_general``
--------------------------------------------

.. _parameter_outdir:


PARAMETER ``outdir``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default "output". Root output directory.

Details
  The root output directory where all output is stored.

.. _parameter_regionExtension:


PARAMETER ``regionExtension``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer >= 0. Default 100. Target region extension in base pairs.

Details
  Specifies the number of base pairs each target region (from the peaks file) should be extended in both 5’ and 3’ direction.

.. _parameter_comparisonType:


PARAMETER ``comparisonType``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default "".

Details
  This parameter helps to organize complex analysis for which multiple different types of comparisons should be done. Set it to a short but descriptive name that summarizes the type of comparison you are making or the types of cells you compare. The value of this parameter appears as prefix in most output files created by the pipeline. It may also be empty.

.. _parameter_designContrast:


PARAMETER ``designContrast``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default  *conditionSummary*. Design formula for the differential accessibility analysis in *DESeq2*.

Details
  This important parameter defines the actual contrast that is done in the differential analysis. That is, which groups of samples are being compared? Examples include mutant vs wild type, mutated vs. unmutated, etc. The last element in the formula must always be *conditionSummary*, which defines the two groups that are being compared. This name is currently hard-coded and required by the pipeline. Our pipeline allows including additional variables to model potential confounding variables, like gender, batches etc. For each additional variable that is part of the formula, a corresponding and identically named column in the sample summary file must be specified. For example, for an analysis that also includes the batch number of the samples, you may specify this as "*~ Treatment + conditionSummary*".

.. _parameter_designVariableTypes:


PARAMETER ``designVariableTypes``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default  *conditionSummary:factor*.   The data types of all elements listed in ``designContrast`` (:ref:`parameter_designContrast`).

Details
Names must be separated by commas, spaces are allowed and will be eliminated automatically. The data type must be specified with a “:”, followed by either “numeric”, “integer”, “logical”, or “factor”. For example, if ``designContrast`` (:ref:`parameter_designContrast`) is specified as "*~ Treatment + conditionSummary*", the corresponding types might be "Treatment:factor, conditionSummary:factor". If a data type is specified as either "logical" or "factor", the variable will be treated as a discrete variable with a finite number of distinct possibilities (something like batch, for example). *conditionSummary* is usually specified as factor because you want to make a pairwise comparison of exactly two conditions. If *conditionSummary* is specified as "integer" or "numeric", however, the variable is treated as continuously-scaled, which changes the interpretation of the results, see the note below.

  .. note:: Importantly, if the variable of interest is continuous-valued (i.e., marked as being integer or numeric), then the reported log2 fold change is per unit of change of that variable. That is, in the final circular plot, TFs displayed in the left side have a negative slope  per unit of change of that variable, while TFs at the right side have a positive one.

.. _parameter_nPermutations:


PARAMETER ``nPermutations``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer >=  0. Default 0. The number of random sample permutations.

Details
  If set to a value > 0, in addition to the real and non-permuted data, the sample conditions as specified in the sample table will be randomly permuted nPermutations times. Specifically, for the condition with fewer samples, 50% will be randomly chosen and switched with the same number of samples from the other condition. This procedure maximizes randomized conditions.

  .. note:: The running time of the pipeline will be increased when using permutations because parts of the pipeline are run for each permutation.

.. _parameter_TFs:


PARAMETER ``TFs``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default "all". Either "all" or a comma-separated list of TF names of TFs to include. If set to "all", all TFs that are found in the directory as specified in ``dir_TFBS`` (:ref:`parameter_dir_TFBS`) will be used.

Details
  If the analysis should be restricted to a subset of TFs, list the names of the TF to include in a comma-separated manner here.

  .. note:: For each TF ``{TF}``, a corresponding file ``{TF}_TFBS.bed`` needs to be present in the directory that is specified by ``dir_TFBS`` (:ref:`parameter_dir_TFBS`).

  .. warning:: We strongly recommending running diffTF with as many TF as possible due to our statistical model that we use that compares against a background model.

.. _parameter_dir_scripts:


PARAMETER ``dir_scripts``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. The path to the directory where the R scripts for running the pipeline are stored.

Details
  .. warning:: The folder name must be ``R``, and it has to be located in the same folder as the ``Snakefile``.

.. _parameter_RNASeqIntegration:


PARAMETER ``RNASeqIntegration``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Logical. true or false. Default false. Should RNA-Seq data be integrated into the pipeline?

Details
  If set to true, RNA-Seq counts as specified in ``RNASeqCounts`` (:ref:`parameter_RNASeqCounts`) will be used to classify each TF into either “activator”, “repressor”, “unknown”, or “not-expressed” for the final circular visualization and the summary table.

  .. note::RNA-Seq integration is only included in the very last step of the pipeline, so it can also be easily integrated later.


SECTION ``samples``
--------------------------------------------

.. _parameter_summaryFile:


PARAMETER ``summaryFile``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  String. Default "samples.tsv". Path to the sample metadata file.

Details
  Path to a tab-separated file that summarizes the input data. See the section :ref:`section_metadata` and the example file for how this file should look like.


SECTION ``peaks``
--------------------------------------------

.. _parameter_consensusPeaks:


PARAMETER ``consensusPeaks``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  String. Default "" (empty). Path to the consensus peak file.

Details
  If set to the empty string "", the pipeline will generate a consensus peaks out of the peak files from each individual sample. For this, you need to provide the following two things:

  - a peak file for each sample in the metadata file in the column *peaks*, see the section :ref:`section_metadata` for details.
  - The format of the peak files, as specified in ``peakType`` (:ref:`parameter_peakType`)

  If a file is provided, it must be a valid *BED* file with at least 3 columns:

  - tab-separated columns
  - no column names in the first row
  - Columns 1 to 3:

    1. Chromosome
    2. Start position
    3. End position

  - Optional:

    4. Identifier (will be made unique for each if this is not the case already)
    5.  Score
    6. Strand

.. _parameter_peakType:


PARAMETER ``peakType``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  String. Default ``narrow``. Format of the peaks. Only relevant if no consensus peak file has been provided (i.e., the :ref:`parameter_consensusPeaks` is empty).

Details
  Only needed if no consensus peak set has been provided. All individual peak files must be in the same format. See the help for *DiffBind dba* for a full list of supported formats, the most common ones include:

  - ``raw``: text file file; peak score is in fourth column
  - ``bed``: .bed file; peak score is in fifth column
  - ``narrow``: default peak.format: narrowPeaks file (from *MACS2*)

.. _parameter_minOverlap:


PARAMETER ``minOverlap``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer >= 0 or Float between 0 and 1. Default 2. Minimum overlap for peak files for a peak to be considered into the consensus peak set. Corresponds to the ``minOverlap`` argument in the *dba* function of *DiffBind*. Only relevant if no consensus peak file has been provided (i.e., ``consensusPeaks``, :ref:`parameter_consensusPeaks`, is empty).

Details
  Only include peaks in at least this many peak sets in the main binding matrix. If set to a value between zero and one, peak will be included from at least this proportion of peak sets. For more information, see the ``minOverlap`` argument in the *dba* function of *DiffBind*  `(see here) <http://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf>`_.


SECTION ``additionalInputFiles``
--------------------------------------------

.. _parameter_refGenome_fasta:


PARAMETER ``refGenome_fasta``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default ‘hg19.fasta’. Path to the reference genome *fasta* file.
Details

  .. warning:: You need write access to the directory in which the *fasta* file is stored, make sure this is the case or copy the *fasta* file to a different directory. The reason is that the pipeline produces a *fasta* index file, which is put in the same directory as the corresponding *fasta* file. This is a limitation of *samtools faidx* and not our pipeline.

  .. note:: This file has to be in concordance with the input data; that is, the exact same genome assembly version must be used. In the first step of the pipeline, this is checked explicitly, and any mismatches will result in an error.

.. _parameter_dir_TFBS:


PARAMETER ``dir_TFBS``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Path to the directory where the TF-specific files for TFBS results are stored.

Details
  Each TF *{TF}* has to have one *BED* file, in the format *{TF}.bed*.  Each file must be a valid *BED6* file with 6 columns, as follows:

  1. chromosome
  2. start
  3. end
  4. ID (or sequence)
  5. score or any other numeric column
  6. strand

  For user convenience, we provide these files as described in the publication as a separate download:

  - hg19: For a pre-compiled list of 620 human TF with in-silico predicted TFBS based on the *HOCOMOCO 10* database and *PWMScan* for hg19, `download this file: <https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_hg19_PWMScan_HOCOMOCOv10.tar.gz>`_
  - hg38: For a pre-compiled list of 771 human TF with in-silico predicted TFBS based on the *HOCOMOCO 11* database and *FIMO* from the MEME suite1 for hg38, `download this file: <https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_hg38_FIMO_HOCOMOCOv11.tar.gz>`_
  - mm10: For a pre-compiled list of 423 mouse TF with in-silico predicted TFBS based on the *HOCOMOCO 10* database and *PWMScan* for mm10, `download this file: <https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_mm10_PWMScan_HOCOMOCOv10.tar.gz>`_

  However, you may also manually create these files to include additional TF of your choice or to be more or less stringent with the predicted TFBS. For this, you only need PWMs for the TF of interest and then a motif prediction tool such as *FIMO* or *MOODS*.

.. _parameter_RNASeqCounts:


PARAMETER ``RNASeqCounts``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default “”. Path to the file with RNA-Seq counts.

Details
  If no RNA-Seq data is included, set to the empty string “”. Otherwise, if ``RNASeqIntegration`` (:ref:`parameter_RNASeqIntegration`) is set to true,  specify the path to a tab-separated file with normalized RNA-Seq counts. It does not matter whether the values have been variance-stabilized or not, as long as values across samples are comparable. Also, consider filtering lowly expressed genes. For guidance, you may want to read `Question 4 here <https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html>`_.

  The first line must be used for labeling the samples, with column names being identical to the sample names as specific in the sample summary table (``summaryFile``, :ref:`parameter_summaryFile`). If you have RNA-Seq data for only a subset of the input samples, this is no problem - the classification will then naturally only be based on the subset. The first column must be named ENSEMBL and it must contain ENSEMBL IDs (e.g., *ENSG00000028277*) without dots. The IDs are then matched to the IDs from the file as specified in ``HOCOMOCO_mapping`` (:ref:`parameter_HOCOMOCO_mapping`).

.. _parameter_HOCOMOCO_mapping:


PARAMETER ``HOCOMOCO_mapping``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Path to the TF-Gene translation table.

Details
  If RNA-Seq integration shall be used, a translation table to associate TFs and ENSEMBL genes is needed. For convenience, we provide such a translation table compatible with the pre-provided TFBS lists. Specifically, for each of the currently three TFBS lists, we provide corresponding translation tables for:

  1. hg19 with HOCOMOCO 10
  2. hg38 with HOCOMOCO 11
  3. mm10 with HOCOMOCO 10

  If you want to create your own version, check the example translation tables and construct one with an identical structure.

.. _section_metadata:


Input metadata
=============================================

This file summarizes the data and corresponding available metadata  that should be used for the analysis. The format is flexible and may contain additional columns that are ignored by the pipeline, so it can be used to capture all available information in a single place. Importantly, the file must be saved as tab-separated, the exact name does not matter as long as it is correctly specified in the configuration file.
It must contain at least contain the following columns (the exact names do matter):

- ``sampleID``: The ID of the sample
- ``bamReads``:  path to the *BAM* file corresponding to the sample.

  .. warning:: All *BAM* files must meet *SAM* format specifications. You may use the program *ValidateSamFile* from the *Picard tools* to check and identify problems with your file. Chromosome names must have a "*chr*" as prefix, otherwise diffTF may crash.

- ``peaks``: absolute path to the sample-specific peak file, in the format as given by ``peakType`` (:ref:`parameter_peakType`). Only needed if no consensus peak file is provided.
- ``conditionSummary``: String with an arbitrary condition name that defines which condition the sample belongs to. There must be only exactly two different conditions across all samples (e.g., *mutated and unmutated*, *day0 and day10*, ...)
- if applicable, all additional variables from the design formula except ``conditionSummary`` must also be present as a separate column.


.. warning:: Do not change the samples data after you started an analysis. You may introduce inconsistencies that will result in error messages. If you need to alter the sample data, we strongly advise to recalculate all steps in the pipeline.

Output
************************************************************

The pipeline produces quite a large number of output files, only some of which are however relevant for the regular user.

.. note:: In the following, the directory structure and the files are briefly outlined. As some directory or file names depend on specific parameters in the configuration file, curly brackets will be used to denote that the filename depends on a particular parameter or name. For example, ``{comparisonType}`` and ``{regionExtension}`` refer to ``comparisonType`` (:ref:`parameter_comparisonType`) and ``regionExtension`` ( :ref:`parameter_regionExtension`) as specified in the configuration file.

Most files have one of the following file formats:

- .bed.gz (gzipped bed file)
- .tsv (tab-separated value, text file with tab as column separators)
- .rds (binary R format, read into with the function ``readRDS``)
- .pdf (PDF format)
- .log (text format)

FOLDER ``FINAL_OUTPUT``
=============================================

In this folder, the final output files are stored. Most users want to examine the files in here for further analysis.


Sub-folder ``extension{regionExtension}``
----------------------------------------------

Stores results related to the user-specified extension size (``regionExtension``, :ref:`parameter_regionExtension`)

.. note:: In all output files, in the column ``permutation``, 0 always refers to the non-permuted, real data, while permutations > 0 reflect real permutations.

FILE ``{comparisonType}.allMotifs.tsv.gz``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Summary table for each TFBS

Details
  Columns are as follows:

  - *permutation*: The number of the permutation.
  - *TF*: name of the TF
  - *chr*, *MSS*, *MES*, *strand*, *TFBSID*: Genomic location and identifier of the (extended) TFBS
  - *PSS*, *PES*, *peakID*:  Genomic location and annotation of the overlapping peak region
  - *baseMean*, *log2FoldChange*, *lfcSE*, *stat*, *pvalue*, *padj*: Results from the *DESeq2* analysis. See the `DESeq2 documentation <https://www.bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html>`_ for details.

FILE ``{comparisonType}.TF_vs_peak_distribution.tsv``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  This summary table contains various results regarding TFs, their log2 fold change distribution across all TFBS and differences between all TFBS and the peaks

Details
  Columns are as follows:

  - *TF*: name of the TF
  - *permutation*: The number of the permutation.
  - *Pos_l2FC*, *Mean_l2FC*, *Median_l2FC*, *sd_l2FC*, *Mode_l2FC*, *skewness_l2FC*: fraction of positive values, mean, median, standard deviation, mode value and Bickel's measure of skewness of the log2 fold change distribution across all TFBS
  - *pvalue_raw* and *pvalue_adj*: raw and adjusted (fdr) p-value of the  t-test
  - *T_statistic*: the value of the T statistic from the t-test
  - *TFBS_num*: number of TFBS
  - *Diff_mean*, *Diff_median*, *Diff_mode*, *Diff_skew*: Difference of the mean, median, mode, and skewness between the log2 fold-change distribution across all TFBS and the peaks, respectively

FILE ``{comparisonType}.summary.tsv``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  The final summary table that is also used for the final circular visualization

Details
  The columns are as follows:

  - *TF*: name of the TF
  - *permutation*: The number of the permutation.
  - *weighted_meanDifference*: the weighted mean difference of the real and background distribution across all CG bins. This value is the basis for the final calculation of the x-axis position for the circular plot. Without permutations, this value is shown, while for permutations, the weighted_meanDifference_enrichment is used.
  - *weighted_Tstat*: the weighted value of the T statistic across all CG bins.
  - *TFBS*: The number of TF binding sites that overlap with the peaks
  - *variance*: Estimated variance for the weighted_Tstat estimate that is incorporated to calculate the p-value
  - *weighted_meanDifference_enrichment*: If permutations are used, the enrichment over the background is calculated based on the weighted_meanDifference values for the real and permuted data, respectively. Specifically, by default, the minimum and maximum of the permuted weighted_meanDifference values is taken as permutation thresholds, and the real weighted_meanDifference values are then divided by these thresholds to calculate an enrichment.
  - *weighted_Tstat_centralized*: a centralized version of the weighted_Tstat values that is used to calculate raw p-values while including the variance also. If possible, the MLE estimate of the distribution median is used for centralization.
  - *pvalue*: raw p-values
  - *pvalueAdj*: adjusted p-values (FDR, Benjamini & Hochberg (1995) correction)
  - *yValue*: the y-value for the plotting
  - *median.cor.tfs*: the median correlation, related to the RNA-Seq classification
  - *classification*: RNA-Seq classification (either activator, undetermined, repressor or not-expressed)
  - *Cohend_factor*	, *weighted_CD*, *weighted_median*, *weighted_sd*: Not used.

FILE ``{comparisonType}.summary.circular.pdf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  The final visualization of the diffTF results

Details
  The PDF contains multiple pages, the structure of which varies depending on the parameters:

  - Number of permutations > 0

    - RNA-Seq integration: Within each of the two following sections, the structure is identical (varying combinations of FDR threshold and RNA-Seq categories, see below)

      - Pages 1-15: Visualization of the results including the permutations
      - Pages 16-30: Visualization of the results excluding the permutations

    - No RNA-Seq integration: Within each of the two following section, the structure is identical (varying combinations of FDR threshold and RNA-Seq categories, see below)

      - Pages 1-5: Visualization of the results including the permutations for different FDR thresholds
      - Pages 6-10: Visualization of the results excluding the permutations for different FDR thresholds

  - Number of permutations = 0

    - RNA-Seq integration: Pages 1-15 show the results in a format as described below for varying combinations of FDR threshold and RNA-Seq categories
    - No RNA-Seq integration: Pages 1-5 show the results for different FDR thresholds

  - If RNA-Seq data is integrated, different combinations of categories are shown on each page (1: activator-undetermined-repressor-not-expressed, 2: activator-undetermined-repressor, 3: activator-repressor), which along with different FDR thresholds gives rise to multiple individual plots. The values on the x-axis denote the effect size (if permutations are incorporated enrichment over background, otherwise the mean difference between the two conditions), with higher values indicating a higher differential TF activity between the two conditions. TFs that are similarly active between the two conditions are close to 0, while TFs more active in either condition are located in the left and right part of the plot. The y-axis (radial position) denotes the statistical significance (adjusted p-values). The significance threshold is indicated as red circular line. TFs that pass the significance threshold are labeled and, if RNA-Seq data is integrated, colored according to their predicted role (see above).


FILE ``{comparisonType}.diagnosticPlots.pdf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Various diagnostic plots for the final TF activity values

Details
  If the number of permutations is larger than 0, the first three pages show various versions of the permuted weighted_meanDifference values and how they relate to the real ones. Permutation 0, as used everywhere throughout the pipeline, contains the real values, while any permutation > 0 refers to an actual permutation. Page 1 shows real and permuted values, page 2 only permuted ones, page 3 a density plot of the real values with the permutation thresholds as dashed lines, inside of which TFs are not labeled as they fall within the permutation and therefore noise area. The next page shows various diagnostic plots from the *locfdr* package to estimate the distribution median, while the remaining plots show histograms of all relevant columns in the final output table for different sets of TFs depending on a specific FDR threshold.

FOLDER ``PEAKS``
=============================================

Stores peak-associated files.


FILES ``{comparisonType}.consensusPeaks.bed`` and ``consensusPeaks_lengthDistribution.pdf``
----------------------------------------------------------------------------------------------

Summary
  Only present if no consensus peak file was provided (``consensusPeaks``, :ref:`parameter_consensusPeaks`). Produced in rule ``filterSexChromosomesAndSortPeaks``. Generated consensus peaks, before filtering (see below) as well as a diagnostic plot showing the length distribution of the peaks.

Details
  Filtered consensus peaks (removal of peaks from one of the following chromosomes: chrX, chrY, chrM, chrUn\*, \*random*, \*hap|_gl\*


FILE ``{comparisonType}.allBams.peaks.overlaps.bed``
----------------------------------------------------

Summary
  Produced in rule ``intersectPeaksAndBAM``. Counts for each consensus peak with each of the input *BAM* files.

Details
  No details provided yet.

FILE ``{comparisonType}.sampleMetadata.rds``
--------------------------------------------

Summary
  Produced in rule ``DESeqPeaks``. Stores data for the input data (similar to the input sample table), for both the real data and the permutations.

Details
  No details provided yet.


FILE ``{comparisonType}.peaks.rds``
--------------------------------------------

Summary
  Produced in rule ``DESeqPeaks``. Stores all peaks that will be used in the analysis.

Details
  No details provided yet.

FILE ``{comparisonType}.peaks.tsv``
--------------------------------------------

Summary
  Produced in rule ``DESeqPeaks``. Stores the *DESeq2* results of the differential accessibility analysis for the peaks.

Details
  No details provided yet.

FILE ``{comparisonType}.normFacs.rds``
--------------------------------------------

Summary
  Produced in rule ``DESeqPeaks``. Gene-specific normalization factors for each sample and peak.

Details
  This file is produces after the differential accessibility analysis for the peaks. The normalization factors will be used for the TF-specific differential accessibility analysis.


FILES ``{comparisonType}.diagnosticPlots.peaks.pdf`` and ``{comparisonType}.diagnosticPlots.peaks_permutation{perm}.pdf`` for each permutation ``{perm}``
---------------------------------------------------------------------------------------------------------------------------------------------------------------------
Summary
  Produced in rule ``DESeqPeaks``. Various diagnostic plots for the differential accessibility peak analysis for the real and permuted data, respectively

Details
  The pages are as follows:

  (1) MA plots
  (2) density plots of normalized and non-normalized counts
  (3) mean-average plots (average of the log-transformed counts vs the fold-change per peak) for each of the sample pairs
  (4) mean SD plots (row standard deviations versus row means)


FILE ``{comparisonType}.DESeq.object.rds``
--------------------------------------------

Summary
  Produced in rule ``DESeqPeaks``. The *DESeq2* object from the differential accessibility peak analysis.

Details
  No details provided yet.

FOLDER ``TF-SPECIFIC``
=============================================

Stores TF-specific files. For each TF ``{TF}``, a separate sub-folder ``{TF}`` is created by the pipeline. Within this folder, the following structure is created:

Sub-folder ``extension{regionExtension}``
----------------------------------------------

FILES ``{TF}.{comparisonType}.allBAMs.overlaps.bed.gz`` and ``{TF}.{comparisonType}.allBAMs.overlaps.bed.summary``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Overlap and *featureCounts* summary file of read counts across all TFBS for all input *BAM* files.

Details
  No details provided yet.


FILE ``{TF}.{comparisonType}.output.tsv``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Produced in rule ``analyzeTF``. A summary table for the *DESeq2* analysis.

Details
  See the file ``{comparisonType}.allMotifs.tsv.gz`` in the ``FINAL_OUTPUT`` folder for a column description.


FILE ``{TF}.{comparisonType}.summary.rds``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
   Produced in rule ``analyzeTF``. A summary table for the log2 fold-changes across all TFBS *DESeq2* results.

Details
  No details provided yet.


FILES ``{TF}.{comparisonType}.diagnosticPlots.pdf`` and ``{TF}.{comparisonType}.diagnosticPlots_permutation{perm}.pdf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Produced in rule ``analyzeTF``. Various diagnostic plots for the differential accessibility TFBS analysis for the real and permuted data, respectively.

Details
  See the description of the file ``{comparisonType}.diagnosticPlots.peaks.pdf`` in the ``PEAKS`` folder, which has an identical structure.


FILES ``{TF}.{comparisonType}.summaryPlots.pdf`` and ``{TF}.{comparisonType}.summaryPlots_permutation{perm}.pdf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Produced in rule ``analyzeTF``. A PDF with a summary of the *DESeq2* analysis for the real and permuted data, respectively.

Details
  Page 1 shows a density plot of the log2 fold-changes for the specific pairwise condition that the user selected, separately for the peaks only and across all TFBS from the specific TF. Page 2 shows the same but in a cumulative representation.


FILE ``{TF}.{comparisonType}.DESeq.object.rds``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Produced in rule ``analyzeTF``. Original *DESeq2* object.

Details
  No details provided yet.

FILE ``{TF}.{comparisonType}.permutationResults.rds``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Produced in rule ``binningTF``. Contains a data frame that stores the results of bin-specific results.

Details
  No details provided yet.

FILE ``{TF}.{comparisonType}.permutationSummary.tsv``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Produced in rule ``binningTF``. A final summary table that summarizes the results across bins by calculating weighted means.

Details
  The data of this table are used for the final visualization.


FILE ``{TF}.{comparisonType}.covarianceResults.rds``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Produced in rule ``binningTF``. Contains a data frame that stores the results of the pairwise bin covariances and the bin-specific weights.

Details
  .. note:: Covariances are only computed for the real data but not the permuted ones.


FOLDER ``LOGS_AND_BENCHMARKS``
=============================================

Stores various log and error files.

- ``*.log`` files from R scripts: Each log file is produced by the corresponding R script and contains debugging information as well as warnings and errors:

  - ``1.produceConsensusPeaks.R.log``
  - ``2.DESeqPeaks.R.log``
  - ``3.analyzeTF.{TF}.R.log`` for each TF ``{TF}``
  - ``4.summary1.R.log``
  - ``5.prepareBinning.log``
  - ``6.binningTF.{TF}.log``  for each TF ``{TF}``
  - ``7.summaryFinal.R.log``

- ``*.log`` summary files: Summary logs for user convenience, produced at very end of the pipeline only. They should contain all errors and warnings from the pipeline run.

  - ``all.errors.log``
  - ``all.warnings.log``

FOLDER ``TEMP``
=============================================

Stores temporary and intermediate files. Since they are usually not relevant for the user, they are explained only very briefly here.

Sub-folder ``SortedBAM``
------------------------------

Stores sorted versions of the original *BAMs* that are optimized for fast count retrieval using *featureCounts*.

- ``{basenameBAM}.bam`` for each input *BAM* file: Produced in rule ``resortBAM``. Resorted *BAM* file

Sub-folder ``extension{regionExtension}``
----------------------------------------------

Stores results related to the user-specified extension size (``regionExtension``, :ref:`parameter_regionExtension`)

- ``{comparisonType}.allTFBS.peaks.bed.gz``: Produced in rule ``intersectPeaksAndTFBS``. *BED* file containing all TFBS from all TF that overlap with the peaks after motif extension
- ``{comparisonType}.{TF}.allTFBS.peaks.extension.saf`` for each TF {TF}: Produced in rule ``intersectTFBSAndBAM``. TF-specific file in SAF format needed for *featureCounts*
- ``conditionComparison.rds``: Produced in rule ``DESeqPeaks``. Stores the condition comparison as a string. Some steps in diffTF need this file as input.
- ``{comparisonType}.motifs.coord.permutation{perm}.bed.gz`` and ``{comparisonType}.motifs.coord.nucContent.permutation{perm}.bed.gz`` for each permutation ``{perm}``: Produced in rule ``calcNucleotideContent``, and needed subsequently for the binning. Temporary and result file of *bedtools nuc*, respectively. The latter contains the GC content for all TFBS.
- ``{comparisonType}.allTFData_processedForPermutations.rd`` and ``{comparisonType}.allTFUniqueData_processedForPermutations.rds``. Produced in rule ``prepareBinning``, and needed subsequently for each ``6.binningTF.R`` step.
- ``{comparisonType}.consensusPeaks.*``: Produced in rule ``filterSexChromosomesAndSortPeaks``. Various temporary files related to the consensus peaks
- ``{comparisonType}.checkParameterValidity.done``: temporary flag file
- ``{TF}_TFBS.sorted.bed`` for each TF ``{TF}``: Produced in rule ``sortPWM``. Coordinate-sorted version of the input TFBS.
- ``{comparisonType}.allTFBS.peaks.bed.gz``: Produced in rule ``intersectPeaksAndTFBS``. *BED* file containing all TFBS from all TF that overlap with the peaks before motif extension

.. _workingWithPipeline:

Working with *diffTF* and FAQs
************************************************************

General remarks
==============================

`diffTF` is programmed as a Snakemake pipeline. Snakemake is a bioinformatics workflow manager that uses workflows that are described via a human readable, Python based language. It offers many advantages to the user because each step can easily be modified, parts of the pipeline can be rerun, and workflows can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition or only minimal modifications. However, with great flexibility comes a price: the learning curve to work with the pipeline might be a bit higher, especially if you have no Snakemake experience. For a deeper understanding and troubleshooting errors, some knowledge of Snakemake is invaluable.

Simply put, Snakemake executes various *rules*. Each *rule* can be thought of as a single *recipe* or task such as sorting a file, running an R script, etc. Each rule has, among other features, a name, an input, an output, and the command that is executed. You can see in the ``Snakefile`` what these rules are and what they do. During the execution, the rule name is displayed, so you know exactly at which step the pipeline is at the given moment. Different rules are connected through their input and output files, so that the output of one rule becomes the input for a subsequent rule, thereby creating *dependencies*, which ultimately leads to the directed acyclic graph (*DAG*) that describes the whole workflow. You have seen such a graph in Section :ref:`workflow`.

In diffTF, a rule is typically executed separately for each TF. One example for a particular rule is sorting the TFBS list for the TF CTCF.

In diffTF, the total number of *jobs* or rules to execute can roughly be calculated as 4 * ``nTF``, where ``nTF`` stands for the number of TFs that are included in the analysis. For each TF, four rules are executed:

1. Sorting the TFBS list (rule ``sortTFBS``)
2. Calculating read counts for each TFBS within the peak regions (rule ``intersectTFBSAndBAM``)
3. Differential accessibility analysis  (rule ``analyzeTF``)
4. Binning step (rule ``binningTF``)

In addition, a few other rules are executed that however do not add up much more to the overall rule count.


.. _timeMemoryRequirements:

Executing diffTF - Running times and memory requirements
===============================================================

*diffTF* is computationally demanding, and depending on the sample size and the number of peaks, a rather high amount of resources may be required to run it. In the following, we discuss various issues related to time and memory requirements and we provide some general guidelines that worked well for us.

.. warning:: We generally advise to run diffTF in a cluster environment. For small analysis, a local analysis on your machine might work just fine (see the example analysis in the Git repository), but running time increases substantially due to limited amount of available cores.

Analysis size
---------------

We now provide a *very rough* classification into small, medium and large with respect to the sample size and the number of peaks:

- Small: Fewer than 10-15 samples, number of peaks not exceeding 50,000-80,000, normal read depth per sample
- Large: Number of samples larger than say 20 or number of peaks clearly exceeds 100,000, or very high read depth per sample
- Medium: Anything between small and large

Memory
---------------

Some notes regarding memory:

- Disk space: Make sure you have enough space left on your device. As a guideline, analysis with 8 samples need around 12 GB of disk space, while a large analysis with 84 samples needs around 45 GB.
- Machine memory: Although most steps of the pipeline have a modest memory footprint, depending on the analysis size, some others (e.g., rule ``prepareBinning``) may need a lot of RAM during execution. We recommend having at least 16 GB available for large analysis (see above)

Number of cores
-----------------

Some notes regarding the number of available cores:

- diffTF can be invoked in a highly parallelized manner, so the more CPUs are available, the better.
- you can use the ``--cores`` option when invoking Snakemake to specify the number of cores that are available for the analysis. If you specify 4 cores, for example, up to 4 rules can be run in parallel (if each of them occupies only 1 core), or 1 rule can use up to 4 cores.
- we strongly recommend running diffTF in a cluster environment due to the massive parallelization. With Snakemake, it is easy to run diffTF in a cluster setting. Simply do the following:

  - write a cluster configuration file that specifies which resources each rule needs. For guidance and user convenience, we provide different cluster configuration files for a small and large analysis. See the folder ``src/clusterConfigurationTemplates`` for examples. Note that these are rough estimates only. See the `Snakemake documentation <http://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html#cluster-configuration>`_ for details for how to use cluster configuration files.
  - invoke Snakemake with one of the available cluster modes, which will depend on your cluster system. We used ``--cluster`` and tested the pipeline extensively with *LSF/BSUB* and *SLURM*. For more details, see the `Snakemake documentation <http://snakemake.readthedocs.io/en/latest/executable.html#cluster-execution>`_

Total running time
--------------------

Some notes regarding the total running time:

- the total running time is very difficult to estimate beforehand and depends on many parameters, most importantly the number of samples, their read depth, the number of peaks, and the number of TF included in the analysis.
- for small analysis such as the example analysis in the Git repository, running times are roughly 30 minutes with 2 cores for 50 TF and a few hours with all 640 TF.
- for large analysis, running time will be up to a day or so when executed on a cluster machine


Frequently asked questions
==============================

Here a few typical use cases, which we will extend regularly in the future if the need arises:

1. I received an error, and the pipeline did not finish.

  As explained in Section :ref:`docs-errors`, you first have to identify and fix the error. Rerunning then becomes trivially easy: just restart Snakemake, it will start off where it left off: at the step that produced that error.

2. I want to rerun a specific part of the pipeline only.

  This common scenario is also easy to solve: Just invoke Snakemake with ``--forcerun {rulename}``, where ``{rulename}`` is the name of the rule as defined in the Snakefile. Snakemake will then rerun the  specified run and all parts downstream of the rule. If you want to avoid rerunning downstream parts (think carefully about it, as there might be changes from the rerunning that might have consequences for downstream parts also), you can combine ``--forcerun`` with ``--until`` and specify the same rule name for both.

3. I want to modify the workflow.

  Simply add or modify rules to the Snakefile, it is as easy as that.


**If you feel that a particular use case is missing, let us know and we will add it here!**



.. _docs-errors:


Handling errors
************************************************************

Error types
==============================

Errors occur during the Snakemake run can principally be divided into:

- Temporary errors (often when running in a cluster setting)

  * might occur due to temporary problems such as bad nodes, file system issues or latencies
  * rerunning usually fixes the problem already. Consider using the option ``--restart-times`` in Snakemake.

- Permanent errors

  * indicates a real error related to the specific command that is executed
  * rerunning does not fix the problem as they are systematic (such as a missing tool)

Identify the cause
==============================

To troubleshoot errors, you have to first locate the exact error. Depending on how you run Snakemake (i.e., in a cluster setting or not), check the following places:

- in locale mode: the Snakemake output appears on the console. Check the output before the line "Error in rule", and try to identify what went wrong.  Errors from R script should in addition be written to the corresponding R log files in the in the ``LOGS_AND_BENCHMARKS`` directory.
- in cluster mode: either error, output or log file of the corresponding rule that threw the error in the ``LOGS_AND_BENCHMARKS`` directory. If you are unsure in which file to look, identify the rule name that caused the error and search for files that contain the rule name in it


Fixing the error
==============================

After locating the error, fix it accordingly. We here provide some guidelines of different error types that may help you fixing the errors you receive:

- Errors related to erroneous input: These errors are easy to fix, and the error message should be indicative. If not, please let us know, and we improve the error message in the pipeline.
- Errors of technical nature: Errors related to memory, missing programs, R libraries etc can be fixed easily by making sure the necessary tools are installed and by executing the pipeline in an environment that provides the required technical requirements. For example, if you receive a memory-related error, try to increase the available memory. In a cluster setting, adjust the cluster configuration file accordingly by either increasing the default memory or (preferably) or by overriding the default values for the specific rule.
- Errors related to the input data: Error messages that indicate the problem might be located in the data are more difficult to fix, and we cannot provide guidelines here. Feel free to contact us.

After fixing the error, rerun Snakemake. Snakemake will continue at the point at which the error message occurred, without rerunning already successfully computed previous steps (unless specified otherwise).
