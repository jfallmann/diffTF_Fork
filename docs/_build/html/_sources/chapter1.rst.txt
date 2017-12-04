.. _docs-quickstart:

============================================================
Try it out now!
============================================================

The following quick start briefly summarizes the necessary steps to use our pipeline:

1. Install the necessary tools (Snakemake, samtools, bedtools, and Subread). We recommend installing them via conda, in which case the installation is as easy as

    ``conda install -c bioconda snakemake bedtools samtools subread``

    If conda is not yet installed, follow the `installation instructions <https://conda.io/docs/user-guide/install/index.html>`_. If you want to install the tools manually and outside of the conda framework, see the following instructions for each of the tools: `snakemake  <http://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_, `samtools <http://www.htslib.org/download>`_, `bedtools <http://bedtools.readthedocs.io/en/latest/content/installation.html>`_, `Subread <http://subread.sourceforge.net>`_.
2. Clone the Git repository:

    ``git clone https://git.embl.de/grp-zaugg/diffTF``
3. To run the example analysis for 50 TF, simply perform the following steps:
  * Change into the ``example/input`` directory within the Git repository

        ``cd diffTF/example/input``
  * Download the data via the download script

        ``sh downloadAllData.sh``
  * To test if the setup is correct, start a dryrun via the first helper script

        ``sh startAnalysisDryRun.sh``
  * Once the dryrun is successful, start the analysis via the second helper script

        ``sh startAnalysis.sh``
4. To run your own analysis, modify the files ``config.json`` and ``sampleData.ts``. See the instructions in the section `Run your own analysis`_ for more details.
5. If your analysis finished successfully, take a look into the ``FINAL_OUTPUT`` folder within your specified output directory, which contains the summary tables and visualization of your analysis. If you received an error, take a look in Section :ref:`docs-errors` to troubleshoot.

.. _docs-prerequisites:

============================================================
Prerequisites
============================================================

This section lists the required software and how to install them. As outlined in Section :ref:`docs-quickstart`, the easiest way is to install all of them via ``conda``. However, it is of course also possible to install the tools separately.

--------------------------
Snakemake
--------------------------

Please ensure that you have at least version 4.3 installed. Principally, there are `multiple ways to install Snakemake <http://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_. We recommend installing it, along with all the other required software, via conda.

----------------------------
samtools, bedtools, Subread
----------------------------

In addition, `samtools <http://www.htslib.org/download>`_, `bedtools <http://bedtools.readthedocs.io>`_ and `Subread <http://subread.sourceforge.net>`_ are needed to run *diffTF*. We recommend installing them, along with all the other required software, via conda.


--------------------------
R and R packages
--------------------------

A working ``R`` installation is needed and a number of packages from either CRAN or Bioconductor have to be installed.  Type the following in ``R`` to install them:

.. code-block:: R

  install.packages(c("checkmate", "futile.logger", "tidyverse", "reshape2", "gridExtra", "scales", "jsonlite", "RcolorBrewer", "rlist", "ggrepel", "lsr", "modeest", "locfdr", "boot"))
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("limma", "vsn", "csaw", "DESeq2", "DiffBind", "geneplotter", "Rsamtools"))



.. _docs-runOwnAnalysis:

============================================================
Run your own analysis
============================================================
