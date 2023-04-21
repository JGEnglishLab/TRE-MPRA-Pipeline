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

#### -r --run-name
* Used to specify the name of the run
* The same string you enter here will be seen in the /runs directory
* The name may not incluc "_" or "/"
* EXAMPLE "-r 19919"

#### -p --path_to_fastq
* Used to specify the path to the dirctory containing the fastq files for your run
* The path may be absolute or direct
* The fastq files inside the directory may be ".fq", ".fastq", ".fq.gz", or ".fastq.gz" files
* EXAMPLE "-r ../sequence_data/19919_fq_files"

#### -t --path_to_treatment_tsv
- Specify the path to the treatment TSV
- Treatment TSV is used to pair sample numbers with treatment names
- The first column should correspond to sample numbers
- The second column should correspond to treatment names
- At least one DNA sample must be present for the pipeline to work. All DNA samples should simply be labeled as "DNA"
- If you have multiple replicates of the same treatment they must be the same name
- Do not include a header for your tsv file
- See example of TSV below

```
1   Serum Free
2   Serum Free
3   Serum Free
4   ATP
5   Forskolin
6   DNA
7   DNA
8   DNA
```

#### -dt --path_to_dna_tsv
- Specify the path to the DNA TSV
- DNA TSV is only needed if more than one DNA sample is present in treatment TSV
- The DNA TSV is used to specify which DNA samples correspond to which RNA samples
- Do not include a header for your tsv file
- The first column corresponds to DNA sample numbers, the second column corresponds to RNA sample numbers
- Each RNA sample must correspond to 1 and only 1 DNA sample
- If you had a treatment TSV that looked like the one above (see -t) you would creat a DNA TSV like the one below
```
6   1
6   2
6   3
7   4
8   5
```
- This would mean the DNA from sample number 6 correspond to the Serum Free treatments, the DNA from sample 7 corresponds to the ATP treatment and the DNA from sample 8 corresponds to the Forskolin treatment.

#### -sr --sample_number_regex
- Used to specify a regex to detect sample numbers from fastq file names
- By default, and number directly following "S" will be counted as the sample number.
- If your file names looked like this 19919X42_220722_A00421_0459_AHH3JFDRX2_S42_L001_R2_001.fastq.gz
- The default would call this file name to be sample number 42
- If your file names looked like this 19919_R1_L001_sample_number_42.fastq.gz
- You would put " -sr sample_number_" to automatically detect the 42 as the sample number 

#### -d --path_to_DNA_fastq
- Sometimes you may not keep the DNA fastq files in the same location that your RNA fastq files
- Use this flag to specify the path to the directory containing the fastq files for your DNA samples

### Flags TMP_comparative.py 

#### -p --pairwise_tsv
- Used to specify a path to a tsv file containg pairwise comparisons
- If this flag is left blank TMP_comparative.py will ask for user input to create the TSV
- If you make your own tsv you must use this header 
- `id	base_treatment	stim_treatment	base_run	stim_run`
- The id column should just contain a unique integer for each row
- See example below

```
id	base_treatment	stim_treatment	base_run	stim_run
1	Serum-Free	Marin-1	19664	20250
2	ADRB2-2	ADRB2-5	20250	20250
3	Marin-1	ADRB2-5	20250	20250
4	FBS	Marin-1	19664	20250
```

#### -m --multi_tsv
- Used to specify a path to a tsv file containg multivariate comparisons
- If this flag is left blank TMP_comparative.py will ask for user input to create the TSV
- If you make your own tsv you must use this header 
- `id	treatment	run`
- The id column should be used to specify the multivariate comparisons, for example if you wanted to compare treatments A, B and C against each other, you would give them all the same id
- See example below

```
id	treatment	run
1	Marin-1	20250
1	Marin-2	20250
1	Marin-3	20250
2	Serum Free	19664
2	ADRB2-2	20250
2	ADRB2-5	20250
3	Serum Free	19664
3	ADRB2-5	20250
3	Marin-3	20250
3	FBS	19664
```

#### -n --n_workers
- Used to specify the number of threads
- Default is 6


