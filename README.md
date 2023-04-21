# TMP

Welcome to the TRE, massively parallel reporter assay, pipeline (TMP) github. This code was designed for users to be able to analyze their own fastq after performing MPRA studies using the TRE library from the english lab.


## Downloading instructions

TMP has a few dependencies required to run. If you do not want to download these software, we are currenltly working on a docker container that will include all of the needed dependencies.

### Required dependencies

* [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
  * Miniconda will be used to download Starcode and Trim-Galore
  * (See set up section for more details on setting up a conda environment)
* [Starcode](https://anaconda.org/bioconda/starcode)
* [Trim-galore](https://anaconda.org/bioconda/trim-galore)
* [ea-utils](https://expressionanalysis.github.io/ea-utils/)
  * Once it has been downloaded and built (see additional [instructions](https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/Compiling.md))you must add it to your path.
  * [Here](https://gist.github.com/nex3/c395b2f8fd4b02068be37c961301caa7) are some instructions on how to add something to your path if you are unfamiliar.
* [R](https://www.r-project.org/)
  * R will require the following packages
    * tidyverse
    * dplyr
    * testit
    * stringr
    * BiocParallel
    * MPRAnalyze
* Python
  * Python will require the following packages
    * numpy
    * pandas
    * biopython




