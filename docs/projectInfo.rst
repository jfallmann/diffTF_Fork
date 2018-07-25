.. _docs-project:

Biological motivation
============================
Transcription factor (TF) activity constitutes an important readout of cellular signalling pathways and thus for assessing regulatory differences across conditions. However, current technologies lack the ability to simultaneously assessing activity changes for multiple TFs and surprisingly little is known about whether a TF acts as repressor or activator. To this end, we introduce the widely applicable genome-wide method diffTF to assess differential TF binding activity and classifying TFs as activator or repressor by integrating any type of genome-wide chromatin with RNA-Seq data and in-silico predicted TF binding sites.

For a graphical summary of the idea, see the section :ref:`workflow`


Help, contribute and contact
============================

If you have questions or comments, feel free to contact us. We will be happy to answer any questions related to this project as well as questions related to the software implementation. For method-related questions, contact Judith B. Zaugg (judith.zaugg@embl.de) or Ivan Berest (berest@embl.de). For technical questions, contact Christian Arnold (christian.arnold@embl.de).

If you have questions, doubts, ideas or problems, please use the `Bitbucket Issue Tracker <https://bitbucket.org/chrarnold/diffTF>`_. We will respond in a timely manner.

Citation
============================

If you use this software, please cite the following reference:

Ivan Berest*, Christian Arnold*, Armando Reyes-Palomares, Giovanni Palla, Kasper Dindler Rassmussen, Kristian Helin & Judith B. Zaugg. *Quantification of differential transcription factor activity and multiomic-based classification into activators and repressors: diffTF*. 2018. submitted.


Change log
============================

Version 1.0.1 (2018-07-25)
    - fixed a bug in 2.DiffPeaks.R that sometimes caused the step to fail, thanks to Jonas Ungerbeck for letting us know
    - fixed a bug in 3.analyzeTF for rare corner cases when DeSeq fails

Version 1.0 (2018-07-01)
    - released stable version



License
============================


diffTF is licensed under the MIT License:

.. literalinclude:: ../LICENSE.md
    :language: text
