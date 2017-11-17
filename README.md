![logo|small](/doc/logo.png "diffTF logo")

Genome-wide quantification of differential transcription factor activity: diffTF         
=====================================================




Summary
-------
Transcription factor (TF) activity constitutes an important readout of cellular signalling pathways and thus for assessing regulatory differences across conditions. However, current technologies lack the ability to simultaneously assessing activity changes for multiple TFs and surprisingly little is known about whether a TF acts as repressor or activator. To this end, we introduce the widely applicable genome-wide method diffTF to assess differential TF binding activity and classifying TFs as activator or repressor by integrating any type of genome-wide chromatin with RNA-Seq data and in-silico predicted TF binding sites

Documentation
-------

A detailed Documentation is available here:

[Documentation](https://git.embl.de/carnold/diffTF/blob/master/doc/Documentation.pdf)

Instructions to run the example analysis
-----------
  * Change into the example/input directory within the Git repository
    * ``cd diffTF/example/input``
  * Download the data via the download script
    * ``sh downloadAllData.sh``
  * Start a dryrun via the helper script
    * ``sh startAnalysis.sh``

  * Once the dryrun is successful, change the startAnalysis.sh script and remove the dryrun directive to start the analysis and restart the script.
    * ``{EDIT THE FILE}``
    * ``sh startAnalysis.sh``



Citation
--------
*Please cite the following article if you use diffTF in your research*:
  * Ivan Berest*, Christian Arnold*, Armando Reyes-Palomares, Kasper Rassmussen & Judith B. Zaugg. Genome-wide quantification of differential transcription factor activity: diffTF. 2017. submitted.
