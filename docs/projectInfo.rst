.. _docs-project:

Biological motivation
============================
Transcription factor (TF) activity constitutes an important readout of cellular signalling pathways and thus for assessing regulatory differences across conditions. However, current technologies lack the ability to simultaneously assessing activity changes for multiple TFs and surprisingly little is known about whether a TF acts as repressor or activator. To this end, we introduce the widely applicable genome-wide method diffTF to assess differential TF binding activity and classifying TFs as activator or repressor by integrating any type of genome-wide chromatin with RNA-Seq data and in-silico predicted TF binding sites.

For a graphical summary of the idea, see the section :ref:`workflow`

We also put the paper on *bioRxiv*, please read all methodological details here:
`Quantification of differential transcription factor activity and multiomic-based classification into activators and repressors: diffTF <https://www.biorxiv.org/content/early/2018/07/13/368498>`_.


Help, contribute and contact
============================

If you have questions or comments, feel free to contact us. We will be happy to answer any questions related to this project as well as questions related to the software implementation. For method-related questions, contact Judith B. Zaugg (judith.zaugg@embl.de) or Ivan Berest (berest@embl.de). For technical questions, contact Christian Arnold (christian.arnold@embl.de).

If you have questions, doubts, ideas or problems, please use the `Bitbucket Issue Tracker <https://bitbucket.org/chrarnold/diffTF>`_. We will respond in a timely manner.

Citation
============================

If you use this software, please cite the following reference:

Ivan Berest*, Christian Arnold*, Armando Reyes-Palomares, Giovanni Palla, Kasper Dindler Rassmussen, Kristian Helin & Judith B. Zaugg. *Quantification of differential transcription factor activity and multiomic-based classification into activators and repressors: diffTF*. 2018. in review

We also put the paper on *bioRxiv*, please read all methodological details here:
`Quantification of differential transcription factor activity and multiomic-based classification into activators and repressors: diffTF <https://www.biorxiv.org/content/early/2018/07/13/368498>`_.


Change log
============================

Version 1.1 (2018-07-27)
    - added a new parameter *dir_TFBS_sorted* in the config file to specify that the TFBS input files are already sorted, which saves some computation time by not resorting them
    - updated the TFBS files that are available via download (some files were not presorted correctly)
    - added support for single-end BAM files. There is a new parameter *pairedEnd* in the config file that specifies whether reads are paired-end or not.
    - restructured some of the permutation-related output files to save space and computation time. The rule *concatenateMotifsPerm* should now be much faster, and the TF-specific *...outputPerm.tsv.gz* files are now much smaller due to an improved column structure

Version 1.0.1 (2018-07-25)
    - fixed a bug in *2.DiffPeaks.R* that sometimes caused the step to fail, thanks to Jonas Ungerbeck for letting us know
    - fixed a bug in *3.analyzeTF* for rare corner cases when *DESeq* fails

Version 1.0 (2018-07-01)
    - released stable version



License
============================


diffTF is licensed under the MIT License:

.. literalinclude:: ../LICENSE.md
    :language: text
