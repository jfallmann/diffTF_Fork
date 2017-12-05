.. _docs-quickstart:

============================================================
Try it out now!
============================================================

The following quick start briefly summarizes the necessary steps to use our pipeline:

1. Install the necessary tools (Snakemake, samtools, bedtools, and Subread).
  .. note:: Note that all tools require Python 3.

  We recommend installing them via conda, in which case the installation is as easy as

  .. code-block:: Bash

    conda install -c bioconda snakemake bedtools samtools subread

  If conda is not yet installed, follow the `installation instructions <https://conda.io/docs/user-guide/install/index.html>`_. Installation is quick and easy.

  .. note:: You do not need to uninstall other Python installations or packages in order to use conda. Even if you already have a system Python, another Python installation from a source such as the macOS Homebrew package manager and globally installed packages from pip such as pandas and NumPy, you do not need to uninstall, remove, or change any of them before using conda.

  If you want to install the tools manually and outside of the conda framework, see the following instructions for each of the tools: `snakemake  <http://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_, `samtools <http://www.htslib.org/download>`_, `bedtools <http://bedtools.readthedocs.io/en/latest/content/installation.html>`_, `Subread <http://subread.sourceforge.net>`_.
2. Clone the Git repository:

    .. code-block:: Bash

      git clone https://git.embl.de/grp-zaugg/diffTF

3. To run the example analysis for 50 TF, simply perform the following steps:

  * Change into the ``example/input`` directory within the Git repository

      .. code-block:: Bash

        cd diffTF/example/input

  * Download the data via the download script

        .. code-block:: Bash

          sh downloadAllData.sh

  * To test if the setup is correct, start a dryrun via the first helper script

        .. code-block:: Bash

          sh startAnalysisDryRun.sh

  * Once the dryrun is successful, start the analysis via the second helper script

        .. code-block:: Bash

          sh startAnalysis.sh

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

Running your own analysis is almost as easy as running the example analysis. Carefully read and follow the following steps and notes:

1. Copy the files ``config.json`` and ``startAnalysis.sh`` to a directory of your choice.
2. Modify the file ``config.json`` accordingly. For example, we strongly recommend running the analysis for all TF instead of just 50 as for the example analysis. For this, simply change the parameter “TFs” to “all”. See Section 4 (TODO) for details about the meaning of the parameters. Do not delete or rename any parameters or sections.
3. Create a tab-separated file that defines the input data, in analogy to the file ``sampleData.tsv`` from the example analysis, and refer to that in the file ``config.json`` (parameter ``summaryFile``)
4. Adapt the file ``startAnalysis.sh`` if necessary (the exact command line call to Snakemake and the various Snakemake-related parameters)

.. note::
  - For Snakemake to run properly, the R folder with the scripts has to be in the same folder as the Snakefile.
  - Since running the pipeline might be computationally demanding, make sure you have enough space left on your device. As a guideline, analysis with 8 samples need around 12 GB of disk space, while a large analysis with 84 samples needs around 45 GB. Also, adjust the number of available cores accordingly. The pipeline can be invoked in a highly parallelized manner, so the more cores are available, the better!
  - The pipeline is written in Snakemake, so for a deeper understanding and troubleshooting errors, some knowledge of Snakemake is invaluable. The same holds true for running the pipeline in a cluster setting. We recommend using a proper cluster configuration file in addition. For guidance and user convenience, we provide different cluster configuration files for a small (up to 10-15 samples) and large (>15 samples) analysis. See the folder ``src/clusterConfigurationTemplates`` for examples. Note that the sample number guidelines above are very rough estimates only. See the `Snakemake documentation <http://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html#cluster-configuration>`_ for details for how to use cluster configuration files.
