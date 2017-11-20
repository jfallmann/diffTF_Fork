![logo|small](/doc/logo.png "diffTF logo")




Genome-wide quantification of differential transcription factor activity: diffTF         
=====================================================



Summary
-------
Transcription factor (TF) activity constitutes an important readout of cellular signalling pathways and thus for assessing regulatory differences across conditions. However, current technologies lack the ability to simultaneously assessing activity changes for multiple TFs and surprisingly little is known about whether a TF acts as repressor or activator. To this end, we introduce the widely applicable genome-wide method diffTF to assess differential TF binding activity and classifying TFs as activator or repressor by integrating any type of genome-wide chromatin with RNA-Seq data and in-silico predicted TF binding sites

Documentation
-------

[A detailed Documentation is available here](https://git.embl.de/grp-zaugg/diffTF/blob/master/doc/Documentation.pdf)

Installation and Quick Start
-------

The following quick start briefly summarizes the necessary steps to use our pipelines.

1. Install the necessary tools (Snakemake, samtools, and bedtools)
2. Clone the Git repository: 
    ``git clone https://git.embl.de/grp-zaugg/diffTF``
3. To run the example analysis, simply perform the following steps:
  * Change into the *example/input* directory within the Git repository
  
    ``cd diffTF/example/input``
  * Download the data via the download script
  
    ``sh downloadAllData.sh``
  * To test if the setup is correct, start a dryrun via the first helper script
  
    ``sh startAnalysisDryRun.sh``
  * Once the dryrun is successful, start the analysis via the second helper script
  
    ``sh startAnalysis.sh``
4. Modify the example analysis to run your own analysis. See Section 3 and 4 in the [Documentation](https://git.embl.de/grp-zaugg/diffTF/blob/master/doc/Documentation.pdf) for further details if you run into troubles or questions.


Citation
--------
*Please cite the following article if you use diffTF in your research*:

Ivan Berest*, Christian Arnold*, Armando Reyes-Palomares, Kasper Rassmussen & Judith B. Zaugg. Genome-wide quantification of differential transcription factor activity: diffTF. 2017. submitted.
