
def r():
    return """
    Takes a string to specify the name of the current run
    The name may not include "_" or "/" characters
    
    e.g. "-r 19919"
    """

def f():
    return """
    Takes a path to the directory that holds the fastq files
    The path may be absolute or direct
    The fastq files inside the directory may be ".fq", ".fastq", ".fq.gz", or ".fastq.gz" files
    
    e.g. "-f ../sequence_data/19919_fq_files"
    """

def t():
    return """
    Takes a path to the treatment TSV
    The path may be absolute or direct
    Treatment TSV is used to pair sample numbers with treatment names
    The first column should correspond to sample numbers
    The second column should correspond to treatment names
    At least one DNA sample must be present for the pipeline to work. All DNA samples should simply be labeled as "DNA"
    If you have multiple replicates of the same treatment they must be the same name
    The sample names may not include "_" "/" or "|" characters 
    The sample names may not include more than one consecutive space
    Do not include a header for your tsv file
    See example of TSV below
    
    1   Serum Free
    2   Serum Free
    3   Serum Free
    4   ATP
    5   Forskolin
    6   DNA
    7   DNA
    8   DNA
    
    e.g. "-t ../19919_treatments.tsv

    """


def dt():
    return """
    Takes a path to the dna TSV
    The path may be absolute or direct
    DNA TSV is only required if more than one DNA sample is present in treatment TSV
    The DNA TSV is used to specify which DNA samples correspond to which RNA samples
    Do not include a header for your tsv file
    The first column corresponds to DNA sample numbers, the second column corresponds to RNA sample numbers
    Each RNA sample must correspond to 1 and only 1 DNA sample
    If you had a treatment TSV that looked like the one above (see -t) you would create a DNA TSV like the one below

    See example of TSV below

    6   1
    6   2
    6   3
    7   4
    8   5
    
    This would mean the DNA from sample number 6 correspond to the Serum Free treatments, 
    the DNA from sample 7 corresponds to the ATP treatment 
    and the DNA from sample 8 corresponds to the Forskolin treatment.
    
    e.g. "-dt ../19919_dna_map.tsv"
    """

def sr():
    return """
    Takes a string that will be used as a regular expression to extract sample numbers
    Default = "S"
    If your file names looked like this 19919X42_S42_L001_R2_001.fastq.gz
    The default would call this file name to be sample number 42 because 42 directly follows "S"
    If your file names looked like this 19919_sample_number_42_R1_L001.fastq.gz
    You would put " -sr sample_number_" to automatically detect the 42 as the sample number 
    Any numbers directly following the string that you provide will be read as the sample number
    In the case where not every fastq file follows the same pattern, TMP_empirical.py will ask for manual user input for files that do not follow the pattern
    
    e.g. "-sr sample_number_"
    """

def s():
    return """
    Takes a path to a Spike-in file
    The path may be absolute or direct
    Spike-in file should contain one spike in per-line
    Spike-in should be 24 nucleotides and only contain ATGC characters
    See example of a spike-in file
    
    TAAATATGCCTCAGCACCCTGCTG
    AAGACGCGTCACAGACTTATAGAC
    CGGAGACACTTAATAGCCTCTAAC
    ATGTTAGTGAGTGTGCGAAGTAGG
    
    e.g. -s ../spike-in.txt
    """

def d():
    return """
    Takes a path to DNA fastq files
    The path may be absolute or direct    
    Sometimes you may not keep the DNA fastq files in the same location as your RNA fastq files
    This flag is only necessary if the directory specified by -f doesn't contain all the fastq files
    Use this flag to specify the path to the directory containing the fastq files for your DNA fastq files
    
    e.g. "-d ../sequence_data/19919_dna_fq_files
    """
def n():
    return """
    Takes an integer to specify the number of threads to be used by MPRAnalyze
    Default = 1
    
    e.g. "-n 40"
    """


def i():
    return """
    Takes a path to an ignore file
    The path may be absolute or direct
    An ignore file may be used to specify which fastq files in the directories specified by -f and -d should be ignored
    Each file name should appear on a separate line
    See example of ignore file
    
    19919X5_S5_L001_R2_001.fastq.gz
    19919X5_S5_L001_R2_001.fastq.gz
    undetermined.fq
    undetermined.fq

    e.g. "-i ../ignore.txt"
    """


def p():
    return """
    Takes a path to a tsv file containing pairwise comparisons
    The path may be absolute or direct
    If this flag is left blank TMP_comparative.py will ask for user input to create the TSV
    If you make your own tsv you must use this header 
    
    id  base_treatment	stim_treatment	base_run	stim_run`
    
    The id column should just contain a unique integer for each row
    
    See example below
    
    id	base_treatment	stim_treatment	base_run	stim_run
    1	Serum-Free	Marin-1	19664	20250
    2	ADRB2-2	ADRB2-5	20250	20250
    3	Marin-1	ADRB2-5	20250	20250
    4	FBS	Marin-1	19664	20250
    
    e.g. "-p ../pairwis_comparisons.tsv"
    """