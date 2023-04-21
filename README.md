# TMP

Welcome to the TRE, massively parallel reporter assay, pipeline (TMP) github. This code was designed for users to be able to analyze their own fastq after performing MPRA studies using the TRE library from the english lab.


## Downloading Instructions

TMP has a few dependencies required to run. If you do not want to download these software, we are currenltly working on a docker container that will include all of the needed dependencies.

### Required Dependencies

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
* Python 3.9.7
  * Python will require the following packages
    * numpy
    * pandas
    * biopython
    
### Set Up

* General 
  * Clone the repository into the desired location
  * In the cloned repo go to the /barcode_map_data directory and unzip the finalBarcodeMap.csv.zip file
* Setting up a conda environment
  * Set up a new environment `conda create -y -n tmp_env`
  * Activate your new environment `conda activate tmp_env`
  * Run the following commands
    * `conda install python=3.9.7`
    * `conda install -c bioconda trim-galore`
    * `conda install -c bioconda starcode`
    
## Running TMP

### Workflow 

The general work flow of TMP is as follows. For each of your runs, you will run TMP_empirical.py. This will generate the alpha results for all of the treatment typs of that particular run. After running TMP_empirical.py for each of your runs, run TMP_comparative.py to create and run comparisons. TMP_comparative.py allows you to run pairwise comparisons (a baseline vs. a stimulated treatement) and multivareate comparisons (more than 2 treatmens). The results from the comparisons will be found in the multi_comparitive_results/ and pairwise_comparitive_results/ directories.

### Flags TMP_empirical.py 



