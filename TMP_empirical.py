#!/usr/bin/env python3

import os
import sys
import re
import numpy as np
import pandas as pd
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
import argparse
from argparse import RawTextHelpFormatter
import time as t
from datetime import datetime
import subprocess
import help_txt as ht

print("Running TMP_empirical")
DEFAULT_THREADS = 1
BAD_CHARS = ["#","%","&","{","}","<",">","*","?","/","\\","$","'","!","\""," ",":","@","+","`","|","=",","]


def check_dir(
    path,
):  # Checks a directory to make sure that fastq files exist in directory
    for f in os.listdir(path):
        if (
            f.endswith(".fastq")
            or f.endswith(".fq")
            or f.endswith(".fq.gz")
            or f.endswith(".fastq.gz")
        ):
            return True
    return False


def check_y_n_inp(inp, correction_message="Try again"):
    inp = inp.lower()
    while inp != "y" and inp != "n":
        print(f"{inp} is an invalid response")
        print(correction_message)
        inp = input('Enter "y" or "n": ')
    return inp


def check_regex(reg_res, pattern):
    """
    Checks to make sure that candidate sample numbers returned by this regex "S\d+"
    all map to the same sample number
    So if the regex return ["S11", "S11"] or ["S11"] we can assume that the sample number is 11
    if the regex returns [] or ["S11", "S13"] we don't know what the sample number is
    """

    if len(reg_res) == 0:
        return False

    all_candidates = []
    for i in reg_res:
        candidate_sample_num = i.split(pattern)[1]
        all_candidates.append(candidate_sample_num)
    np_ar = np.array(all_candidates)
    return len(np.unique(np_ar)) == 1


def check_sample_pattern(entered_samples):
    """
    checks to make sure that entered samples follow this pattern
    1,4,5 ... etc.
    """
    pattern = re.compile(r"^\d+\s*(\s*,\s*\d+\s*)*$")
    while (
        not pattern.match(entered_samples)
        or entered_samples != pattern.match(entered_samples).group()
    ):
        print(f'"{entered_samples}" doesn\'t match the necessary pattern')
        entered_samples = input("Enter the sample numbers separated by commas: ")
    return entered_samples


def check_sample_numbers(entered_samples, remaining_samples):
    """
    Checks if the user accidentally enters a wrong sample number
    Returns the samples in a list
    """
    entered_samples = check_sample_pattern(entered_samples)

    sep_samples = []

    legal_samples = False
    while not legal_samples:
        for i in entered_samples.split(","):
            i = i.strip()
            if i not in remaining_samples:
                legal_samples = False
                print(f"Sample {i} is not in the remaining samples")
                entered_samples = input("Check samples and enter again: ")
                entered_samples = check_sample_pattern(entered_samples)
                continue
            else:
                sep_samples.append(i)

        legal_samples = True

    return sep_samples


def diff_letters(a, b):
    m = min([len(a), len(b)])
    s = sum(a[i] != b[i] for i in range(m))
    return s


def diff_reads(a, b):
    """
    makes sure the only difference between two file names is R1 and R2
    """
    pattern = re.compile(r"[r|R][1|2]")
    return pattern.search(a).group() != pattern.search(b).group()


def check_treatment_tsv(tsv):
    """
    :param tsv: a path to the tsv input by the user containing sample numbers and treatment types
    :return: a tuple
        a bool if correct or not
        a treatment dictionary
        a list of all sample numbers
        and an error or success message
    """

    df = pd.read_csv(tsv, sep='\t', dtype=str)

    #Make sure it includes the required columns
    if not "sample_number" in df.columns or not "treatment" in df.columns:
        return (
            False,None,None,
            "Treatment TSV must contain \"sample_number\" and \"treatment\" columns",
            None, None, None, None, None, None
        )

    #Make sure that the sample_number column is only numeric
    # if not isdigit(df['sample_number']):
    if not df['sample_number'].str.isdigit().all():

        return (
            False,None,None,
            "Error in Treatment TSV, \"sample_number\" column must be numeric",
            None, None, None, None, None, None
        )

    #Make sure that the sample numbers are all unique
    if not all(df['sample_number'].unique()) == all(df['sample_number']):
        return (
            False,None,None,
            "Error in Treatment TSV, \"sample_number\" column must have no repeated values",
            None, None, None, None, None, None
        )

    #Make sure that the short names don't contain any illegal chars
    for name in df["treatment"]:
        if not all(x not in name for x in BAD_CHARS):
            return (
                False,None,None,
                "Error in Treatment TSV, \"treatment\" column may not contain any of the following characters\n#_%&{}<>*?/\ $'!\":@+`|=,",
                None, None, None, None, None, None
            )
    #Make sure there are no commas
    if any(any(df[col].astype(str).str.contains(",")) for col in df.columns):
        return (
            False,None,None,
            "Error in Treatment TSV, no commas (,) are allowed in the Treatment TSV",
            None, None, None, None, None, None
        )
    # If cell type exists, make sure that they are the same for all replicates
    try:
        s = df.groupby(by="treatment")["cell_type"].nunique()
        if not all(s.loc[~(s.index == "DNA")] <=1):
            return (
                False,None,None,
                "Error in Treatment TSV, cell types must be the same for replicates",
                None, None, None, None, None, None
            )
    except KeyError:
        pass

    # If concentration exists, make sure that they are the same for all replicates
    try:
        s = df.groupby(by="treatment")["concentration"].nunique()
        if not all(s.loc[~(s.index == "DNA")] <=1):
            return (
                False,None,None,
                "Error in Treatment TSV, concentration must be the same for replicates of the same treatment (short name)",
                None, None, None, None, None, None
            )
    except KeyError:
        pass

    # If time exists, make sure that they are the same for all replicates
    try:
        s = df.groupby(by="treatment")["time"].nunique()
        if not all(s.loc[~(s.index == "DNA")] <= 1):
            return (
                False, None, None,
                "Error in Treatment TSV, time must be the same for replicates of the same treatment (short name)",
                None, None, None, None, None, None
            )
    except KeyError:
        pass

    # If time exists, make sure that they are the same for all replicates
    try:
        s = df.groupby(by="treatment")["time"].nunique()
        if not all(s.loc[~(s.index == "DNA")] <= 1):
            return (
                False, None, None,
                "Error in Treatment TSV, time must be the same for replicates of the same treatment (short name)",
                None, None, None, None, None, None
            )
    except KeyError:
        pass

    # If long_name exists, make sure that they are the same for all replicates
    try:
        s = df.groupby(by="treatment")["long_name"].nunique()
        if not all(s.loc[~(s.index == "DNA")] <= 1):
            return (
                False, None, None,
                "Error in Treatment TSV, long name must be the same for replicates of the same treatment (short name)",
                None, None, None, None, None, None
            )
    except KeyError:
        pass

    # If anonymous exists, make sure that they are the same for all replicates
    try:
        s = df.groupby(by="treatment")["anonymous_name"].nunique()
        if not all(s.loc[~(s.index == "DNA")] <= 1):
            return (
                False, None, None,
                "Error in Treatment TSV, anonymous name must be the same for replicates of the same treatment (short name)",
                None, None, None, None, None, None
            )
    except KeyError:
        pass

    # If tag exists, make sure that they are the same for all replicates
    try:
        s = df.groupby(by="treatment")["tag"].nunique()
        if not all(s.loc[~(s.index == "DNA")] <= 1):
            return (
                False, None, None,
                "Error in Treatment TSV, tag must be the same for replicates of the same treatment (short name)",
                None, None, None, None, None, None
            )
    except KeyError:
        pass

    #Make dictionary with sample number and short names
    treatments = dict(zip(df.sample_number, df.treatment))
    sample_list = list(df["sample_number"])

    #Make dictionaries out of optional columns
    try:
        long_names = dict(zip(df.sample_number, df.long_name))
    except:
        long_names = dict.fromkeys(df.sample_number)

    try:
        concentrations = dict(zip(df.sample_number, df.concentration))
    except:
        concentrations = dict.fromkeys(df.sample_number)

    try:
        times = dict(zip(df.sample_number, df.time))
    except:
        times = dict.fromkeys(df.sample_number)

    try:
        cell_types = dict(zip(df.sample_number, df.cell_type))
    except:
        cell_types = dict.fromkeys(df.sample_number)

    try:
        tags = dict(zip(df.sample_number, df.tag))
    except:
        tags = dict.fromkeys(df.sample_number)

    try:
        anonymous_names = dict(zip(df.sample_number, df.anonymous_name))
    except:
        anonymous_names = dict.fromkeys(df.sample_number)

    #All these lines do is replace NaN with "None" for each dict
    anonymous_names = {key: "None" if isinstance(value, float) and np.isnan(value) else value for key, value in anonymous_names.items()}
    tags = {key: "None" if isinstance(value, float) and np.isnan(value) else value for key, value in tags.items()}
    cell_types = {key: "None" if isinstance(value, float) and np.isnan(value) else value for key, value in cell_types.items()}
    times = {key: "None" if isinstance(value, float) and np.isnan(value) else value for key, value in times.items()}
    concentrations = {key: "None" if isinstance(value, float) and np.isnan(value) else value for key, value in concentrations.items()}
    long_names = {key: "None" if isinstance(value, float) and np.isnan(value) else value for key, value in long_names.items()}

    return True, treatments, sample_list, "Treatment TSV looks good", long_names, concentrations, times, cell_types, tags, anonymous_names


def check_DNA_tsv(tsv, dna_samples, rna_samples):
    """
    :param tsv: the tsv input by the user containing DNA sample numbers and corresponding RNA sample Numbers
    :param dna_samples: All of the DNA sample numbers
    :param rna_samples: All of the RNA sample numbers
    :return: A tuple
        A boolean if the DNA TSV is correct
        An error message
    """

    df = pd.read_csv(tsv, sep='\t')

    if not "DNA_sample_number" in df.columns and not "RNA_sample_number" in df.columns:
        error_msg = "Error in DNA TSV, the header must be \nDNA_sample_number\tRNA_sample_number"
        return False, error_msg

    if not is_numeric_dtype(df['DNA_sample_number']) and not is_numeric_dtype(df["RNA_sample_number"]):
        error_msg = "Columns must be numeric"
        return False, error_msg

    if not sorted(dict.fromkeys(list(df["DNA_sample_number"]))) == sorted([int(i) for i in dna_samples]):
        error_msg = f"Error in DNA TSV, the DNA sample numbers in DNA TSV don't match DNA sample numbers in treatment TSV"
        return False, error_msg

    if not sorted(dict.fromkeys(list(df["RNA_sample_number"]))) ==sorted([int(i) for i in rna_samples]):
        error_msg = f"Error in DNA TSV, the RNA sample numbers in DNA TSV don't match RNA sample numbers in treatment TSV"
        return False, error_msg


    return True, "DNA TSV looks good!"


###########################################################################
###########################################################################
#                   		GET ARGUMENTS
###########################################################################
###########################################################################
parser = argparse.ArgumentParser(
    formatter_class=RawTextHelpFormatter,
    prog="Setup snakemake",
    description="Sets up and writes a snakemake workflow to run MPRAnalyze",
)
parser.add_argument("-r", required=False, type=str, help=ht.r_empirical())
parser.add_argument("-f", required=True, type=str, help=ht.f())
parser.add_argument("-t", required=True, type=str, help=ht.t())
parser.add_argument("-dt", required=False, type=str, help=ht.dt())
parser.add_argument("-sr", required=False, type=str, default="S", help=ht.sr())
parser.add_argument("-d", required=False, type=str, help=ht.d())
parser.add_argument("-s", required=False, type=str, help=ht.s())
parser.add_argument("-n", required=False, type=int, default=DEFAULT_THREADS, help=ht.n())
parser.add_argument("-i", required=False, type=str, help=ht.i())

args = parser.parse_args()

dir_name = vars(args)["r"]
fastq_path = vars(args)["f"]
treatment_tsv_path = vars(args)["t"]
dna_tsv_path = vars(args)["dt"]
pattern = vars(args)["sr"]
dna_path = vars(args)["d"]
spike_path = vars(args)["s"]
ignore_path = vars(args)["i"]
threads = vars(args)["n"]


###########################################################################
###########################################################################
#                CHECK SPIKE-IN FILE -s
###########################################################################
###########################################################################

# Will just be an empty list if no spike path is provided, or if all spikes are invalid
valid_spikes = []
if spike_path:
    legal_path = os.path.exists(spike_path)

    while not legal_path:
        print(f"No spike-in file found at {spike_path}")
        spike_path = input("Check path and enter it again: ")
        legal_path = os.path.exists(spike_path)
    print("Spike-in file found")
    abs_spike_path = os.popen(f"readlink -f {spike_path}").read().strip()

    spike_file = open(abs_spike_path, "r+")
    nucleotides = "ACTG"
    for spike in spike_file:
        spike = spike.strip()

        # Errors will be dealt with in the R Scripts make_sample_files_SM.R and make_mpra_input_SM.R
        # We just draw attention to the errors here
        if len(spike) != 24:
            print(f"Error with spike-in {spike}")
            print("Length of spike-in must be 24")
            print("This spike-in will be ignored. If this should not happen start over")

        if not all(i in nucleotides for i in spike):
            print(f"Error with spike-in {spike}")
            print("Spike-ins must only contain A,T,G, and C characters")
            print("This spike-in will be ignored. If this should not happen start over")

    spike_file.close()
else:
    abs_spike_path = "None"

###########################################################################
###########################################################################
#                   		GET NAME OF RUN
###########################################################################
###########################################################################
wd = os.getcwd()


if not dir_name:
    dir_name = input("Enter the name of the run: ")

while not all(x not in dir_name for x in BAD_CHARS):
    print('The name of the run cannot include spaces or any of the following characters \n#_%&{}<>*?/\ $\'!":@+`|=,')
    dir_name = input(
        'Change the name of the run so it does not include spaces or any of the following characters \n#_%&{}<>*?/\ $\'!":@+`|=,\n:'
    )

while os.path.exists(f"./runs/{dir_name}"):
    print(f"The run {dir_name} already exists!")
    print("Adding date time to run name")
    now = datetime.now()
    dt_string = now.strftime("---%d-%m-%Y--%H-%M-%S")
    dir_name = dir_name + dt_string

###########################################################################
###########################################################################
#                   	GET TREATMENT AND DNA TSV
###########################################################################
###########################################################################

legal_path = os.path.exists(treatment_tsv_path) and bool(
    re.search(r".tsv(/)?$", treatment_tsv_path)
)
while not legal_path:
    print(f"Treatment TSV not found in {treatment_tsv_path}")
    treatment_tsv_path = input("Check path and enter it again: ")
    legal_path = os.path.exists(treatment_tsv_path) and bool(
        re.search(r".tsv(/)?$", treatment_tsv_path)
    )
print("Treatment TSV found!")

# treatment_tsv = open(treatment_tsv_path, "r")
correct_tsv, treatments, tsv_sample_numbers, msg, long_names, concentrations, times, cell_types,tags, anonymous_names= check_treatment_tsv(treatment_tsv_path)
while not correct_tsv:
    print("Error with treatment TSV provided")
    print(msg)
    _ = input('Check format and and try again (Enter "y" to continue)')
    correct_tsv, treatments, tsv_sample_numbers, msg, long_names, concentrations, times, cell_types,tags, anonymous_names= check_treatment_tsv(treatment_tsv_path)

noDNAProvided = True
DNA_sample_num = []
RNA_sample_num = []
for t in treatments:
    if treatments[t].lower() == "dna":
        noDNAProvided = False
        treatments[t] = "DNA"
        if t not in DNA_sample_num:
            DNA_sample_num.append(t)
    else:
        if t not in RNA_sample_num:
            RNA_sample_num.append(t)

if noDNAProvided:  # No DNA sample was auto detected in the TSV!
    print("No DNA sample was auto detected in the TSV!")
    print('The treatment TSV must contain at least 1 DNA sample labeled as "DNA"')
    print('If more than one DNA sample is present for run, label all of them as "DNA"')
    exit()

if len(DNA_sample_num) > 1:  # Multiple DNA sample detected.
    print("Multiple DNA samples detected!")

    if not dna_tsv_path:
        print(
            "If multiple DNA samples are present, you must provide a DNA tsv using the -dt flag."
        )
        print("see --help flag for more info")
        exit()

    legal_path = os.path.exists(dna_tsv_path) and bool(re.search(r".tsv(/)?$", dna_tsv_path))
    while not legal_path:
        print(f"DNA tsv not found in {dna_tsv_path}")
        dna_tsv_path = input("Check path and enter it again: ")
        legal_path = os.path.exists(dna_tsv_path) and bool(re.search(r".tsv(/)?$", dna_tsv_path))
    print("DNA tsv found!")

    # DNA_tsv = open(dna_tsv_path, "r")
    correct_tsv, error_msg = check_DNA_tsv(dna_tsv_path, DNA_sample_num, RNA_sample_num)
    while not correct_tsv:
        # DNA_tsv.close()
        print(error_msg)
        _ = input('Check format and and try again (Enter "y" to continue)')
        # DNA_tsv = open(dna_tsv_path, "r")
        correct_tsv, error_msg = check_DNA_tsv(dna_tsv_path, DNA_sample_num, RNA_sample_num)
    print(error_msg)
    abs_dna_tsv_path = os.popen(f"readlink -f {dna_tsv_path}").read().strip()

else:
    abs_dna_tsv_path = None

###########################################################################
###########################################################################
#                   	GET PATH TO FASTQS
###########################################################################
###########################################################################
legal_path = os.path.exists(fastq_path) and check_dir(fastq_path)

while not legal_path:
    if not os.path.exists(fastq_path):
        print(f'Path does not exist: "{fastq_path}"')
        fastq_path = input("Check path and enter it again: ")
    if os.path.exists(fastq_path) and not check_dir(fastq_path):
        print(f'No Fastq files found in path: "{fastq_path}"')
        fastq_path = input("Check path and enter it again: ")
    legal_path = os.path.exists(fastq_path) and check_dir(fastq_path)

# Print out the files
print(f"Files found in {fastq_path}")
for f in os.listdir(fastq_path):
    if bool(re.search(r".f(q|astq)(.gz)?$", f)):
        print(f)

# Get the absolute path to the directory
os.chdir(fastq_path)
fastq_path = os.getcwd()
if not fastq_path.endswith("/"):
    fastq_path = fastq_path + "/"
# Then change pack to the correct dir
os.chdir(wd)

if not dna_path:
    print(
        '\nNOTE!!!\nIf the DNA Fastq files are not contained in the directory above, you must provide a path to them using the "-d" flag'
    )
else:
    legal_path = os.path.exists(dna_path) and check_dir(dna_path)

    while not legal_path:
        if not os.path.exists(dna_path):
            print(f'Path to DNA samples, "{dna_path}" does not exist ')
            dna_path = input("Check path and enter it again: ")
        if os.path.exists(dna_path) and not check_dir(dna_path):
            print(f'No Fastq files found in path to DNA samples "{dna_path}"')
            dna_path = input("Check path and enter it again: ")
        legal_path = os.path.exists(dna_path) and check_dir(dna_path)

    # Print out the files
    print(f"Files found in {dna_path}")
    for f in os.listdir(dna_path):
        if bool(re.search(r".f(q|astq)(.gz)?$", f)):
            print(f)

    # Get the absolute path to the directory
    os.chdir(dna_path)
    dna_path = os.getcwd()
    if not dna_path.endswith("/"):
        dna_path = dna_path + "/"
    # Then change pack to the correct dir
    os.chdir(wd)

###########################################################################
###########################################################################
#                   	CREATE DIRECTORIES
###########################################################################
###########################################################################
os.system(f"mkdir -p ./runs/{dir_name}/raw_counts")
os.system(f"mkdir -p ./runs/{dir_name}/run_stats")
os.system(f"mkdir -p ./runs/{dir_name}/star_code")

# Used for files that need to be merged
os.system(f"mkdir -p ./runs/{dir_name}/trimmed_files")
os.system(f"mkdir -p ./runs/{dir_name}/joined_files")
os.system(f"mkdir -p ./runs/{dir_name}/joined_files/un_files")


###########################################################################
###########################################################################
#                   GET LIST OF FQ FILES TO BE IGNORED
###########################################################################
###########################################################################
ignore_files = []
if ignore_path:
    legal_path = os.path.exists(ignore_path)

    while not legal_path:
        print(f"No ignore file found at {ignore_path}")
        ignore_path = input("Check path and enter it again: ")
        legal_path = os.path.exists(ignore_path)
    print("Ignore file found")
    abs_ignore_path = os.popen(f"readlink -f {ignore_path}").read().strip()

    ignore_file = open(abs_ignore_path, "r+")

    for f in ignore_file:
        ignore_files.append(f.strip())


###########################################################################
###########################################################################
#                   GATHER FILES TO MAKE META DATA
###########################################################################
###########################################################################
file_info = open(f"./runs/{dir_name}/file_info.csv", "w+")
file_info.write("path_to_file,needs_to_join\n")

raw_file_paths = {}
files = []

for filename in os.listdir(fastq_path):
    if bool(re.search(r".f(q|astq)(.gz)?$", filename)) and filename not in ignore_files:
        files.append(filename)
        raw_file_paths[filename] = fastq_path

if dna_path and os.path.exists(dna_path):
    for filename in os.listdir(dna_path):
        if (
            bool(re.search(r".f(q|astq)(.gz)?$", filename))
            and filename not in ignore_files
        ):
            files.append(filename)
            raw_file_paths[filename] = dna_path

# If the files need to be merged, the names will look slightly different

"""
Loops through all the files in the cwd
If the file names end with .fq.gz it will get merge with the file that 
shares its name.

These are example file names that would get merged. 
20250X39_221115_A00421_0496_AHVN35DRX2_S39_L002_R2_001.fq.gz
20250X39_221115_A00421_0496_AHVN35DRX2_S39_L002_R1_001.fq.gz
Note in order to pair files, the only difference can be between R1 and R2
"""
f = open(f"./runs/{dir_name}/joined_files/matching_names.csv", "w+")
# new_files = []
numFiles = len(files)
numPaired = 0

print("Here are the auto detected pairs\n")
print("read1,\tread2,\tjoined name")
# These two will be used for meta_data.csv
# They will contain the altered final names
un_paired_files = []
paired_files = []
# These two will contain un-altered filenames
un_paired_files_raw = []
paired_files_raw = []

file1s = []
file2s = []
non_matching_files = []

for file in files:
    if bool(re.search(r"[r|R]1", file)):
        file1s.append(file)
    elif bool(re.search(r"[r|R]2", file)):
        file2s.append(file)
    else:
        non_matching_files.append(file)

file1s_copy = file1s.copy()
file2s_copy = file2s.copy()

for file1 in file1s:
    for file2 in file2s:
        if diff_letters(file1, file2) == 1 and diff_reads(file1, file2):
            file1s_copy.remove(file1)
            file2s_copy.remove(file2)

            paired_files_raw.append(file1)
            paired_files_raw.append(file2)
            # Get the paths to the files
            file1_path = raw_file_paths[file1]
            file2_path = raw_file_paths[file2]
            file_info.write(f"{file1_path + file1},True\n")
            file_info.write(f"{file2_path + file2},True\n")

            # The file names will have the .24bp_5prim added from trimGalore
            file1 = file1.split(".", 1)[0] + ".24bp_5prime." + file1.split(".", 1)[1]
            file2 = file2.split(".", 1)[0] + ".24bp_5prime." + file2.split(".", 1)[1]
            # trim_galore subs fq for fastq. So we need to change the file names
            file1 = re.sub(".fastq", ".fq", file1)
            file2 = re.sub(".fastq", ".fq", file2)

            finalName = re.sub(r"[r|R]1", "R1-2", file1)
            finalName = finalName.split(".")[0]
            fastq_join_name = f"{finalName}_%.fq"
            f.write(f"{file1},{file2},{fastq_join_name}\n")
            print(file1, "\t", file2, "\t", fastq_join_name)

            finalName = f"{finalName}_join.csv"
            paired_files.append(finalName)

for file1 in file1s_copy:
    file1_path = raw_file_paths[file1]
    file_info.write(f"{file1_path + file1},False\n")

    un_paired_files_raw.append(file1)
    file1 = re.sub(".f(q|astq)(.gz)?", ".csv", file1)
    un_paired_files.append(file1)

if len(file2s_copy) != 0:
    print("Read2 files detected without a matching Read1 file")
    print("These files will be ignored:")
    for file2 in file2s_copy:
        print(f"\t{file2}")

for non_match_file in non_matching_files:
    file1_path = raw_file_paths[non_match_file]
    file_info.write(f"{file1_path + non_match_file},False\n")

    un_paired_files_raw.append(non_match_file)
    non_match_file = re.sub(".f(q|astq)(.gz)?", ".csv", non_match_file)
    un_paired_files.append(non_match_file)

f.close()
file_info.close()

if len(un_paired_files) == 0:
    print("\nAuto detected all pairs!")
    merge = True

elif len(file1) == len(files):
    print("\nNo pairs were auto detected!")
    print(
        'In order to pair files together the only difference in the file names should be between "R1" and "R2"'
    )
    print('The patterns "R1" and "R2" may only be present once in each file name.')
    print("Run will continue without pairing files")
    print("If files must be paired, fix file names and start over.")
    merge = False

else:
    print("\nSome files were not paired")
    print(
        'In order to pair files together the only difference in the file names should be between "R1" and "R2"'
    )
    print('The patterns "R1" and "R2" may only be present once in each file name.')
    print("__________________________")
    for f in un_paired_files_raw:
        print(f)
    print("These files will not be paired")
    print("If they should be paired, fix them to match correct pattern and start over")
    print("__________________________")
    for f in paired_files_raw:
        print(f)
    print("These files will be paired")
    merge = True

files = res_list = [y for x in [paired_files, un_paired_files] for y in x]

fileSamplePair = {}
all_sample_numbers = []

for f in files:
    reg_res = re.findall(pattern + "\d+", f)
    if check_regex(reg_res, pattern):
        sample_number = reg_res[0].replace(pattern, "")  # Get the actual sample number
        fileSamplePair[f] = sample_number
        if sample_number not in all_sample_numbers:
            all_sample_numbers.append(sample_number)
    else:
        print("\n****************************************")
        print('Different file format detected for "' + f + '"')
        sample_number = input('Enter sample number e.g. "14" for sample 14: ')
        while not sample_number.isnumeric():
            print("Sample number must be a number.")
            input(f"Enter sample number for {f}")
        fileSamplePair[f] = sample_number
        if sample_number not in all_sample_numbers:
            all_sample_numbers.append(sample_number)

# Check to see if sample number detected from files and sample numbers from treatment TSV match
if tsv_sample_numbers.sort() != list(map(int, all_sample_numbers)).sort():
    print(
        "The sample numbers provided in treatment TSV and the sample numbers detected from files do not match!\n"
    )
    print(
        "These are the sample numbers from treatment TSV that are NOT detected in the files provided"
    )
    print("-")
    for i in tsv_sample_numbers:
        if i not in all_sample_numbers:
            print("\t", i)
    print("-")
    print(
        "These are the sample numbers from the files provided that are NOT in the treatment TSV"
    )
    print("-")
    for i in all_sample_numbers:
        if i not in tsv_sample_numbers:
            print("\t", i)
    print("-")
    print("Only these samples numbers will be used")
    joined_sample_numbers = list(set(all_sample_numbers) & set(tsv_sample_numbers))
    joined_sample_numbers.sort()
    print(joined_sample_numbers)

    con = input(
        "If you only want to use the above files enter y, otherwise enter n to exit and fixe your files (y/n): "
    )
    con = check_y_n_inp(con)

    if con == "n":
        exit()

    for file in fileSamplePair:
        if fileSamplePair[file] not in joined_sample_numbers:
            del fileSamplePair[file]

f = open(f"./runs/{dir_name}/metaData.tsv", "w+")
f.write("fileName\tsampleNumber\ttreatment\trun_name\tlong_name\tconcentration\ttime\tcell_type\ttag\tanonymous_name\n")

# Write TSV
for file in fileSamplePair:
    #File, sample#, treatment name, run_name
    line = (f"{file}\t{fileSamplePair[file]}\t{treatments[fileSamplePair[file]]}\t{dir_name}\t{long_names[fileSamplePair[file]]}\t{concentrations[fileSamplePair[file]]}\t{times[fileSamplePair[file]]}\t{cell_types[fileSamplePair[file]]}\t{tags[fileSamplePair[file]]}\t{anonymous_names[fileSamplePair[file]]}\n")
    print(line)
    f.write(line)

f.close()

###########################################################################
###########################################################################
#                   	WRITE SNAKE MAKE FILE START
###########################################################################
###########################################################################


sf = open(f"./runs/{dir_name}/Snakefile", "w+")

if merge:
    sf.write('MERGE_DIR = "merged_files/"')

sf.write(
    f"""
TRIM_DIR = "trimmed_files/"
JOIN_DIR = "joined_files/"
STARCODE_DIR = "star_code/"
SCRIPTS_DIR = "../../scripts/"
BARCODE_MAP_DIR = "../../barcode_map_data/"
DESCRIPTIVE_STATS_DIR = "run_stats/"
RAW_COUNTS = "raw_counts/"


rule run_quantitative_analysis:
	message: "Creating stats about run and running MPRAnalyze"
	input: 
		sc = dynamic(STARCODE_DIR + "analyzed_out_sample\u007bn3\u007d_mapped_sc_out.tsv"),
		md = "metaData.tsv",
		bcm = BARCODE_MAP_DIR + "finalBarcodeMap.csv"
	output: 
		DESCRIPTIVE_STATS_DIR + "filtering_ratios.png",
		DESCRIPTIVE_STATS_DIR + "type_ratios.png",
		DESCRIPTIVE_STATS_DIR + "dna_per_barcode.png",
		DESCRIPTIVE_STATS_DIR + "pre_filter_data.csv",
		DESCRIPTIVE_STATS_DIR + "run_summary.csv",
		"{dir_name}__empirical_results.csv",
		"MPRA_data.csv"

	shell: "Rscript " + SCRIPTS_DIR + "run_quantitative_analysis_SM.R \u007binput.md\u007d \u007binput.bcm\u007d {abs_spike_path} {abs_dna_tsv_path} {threads} \u007binput.sc\u007d"

rule analyze_starcode:
	message: "Analyzing Starcode"
	input: dynamic(STARCODE_DIR + "sample\u007bn2\u007d_mapped_sc_out.tsv")
	output: dynamic(STARCODE_DIR + "analyzed_out_sample\u007bn3\u007d_mapped_sc_out.tsv")
	shell: "python3 " + SCRIPTS_DIR + "analyze_star_code_SM.py \u007binput\u007d"

rule run_starcode:
	message: "Running Star Code"
	input: dynamic(STARCODE_DIR + "sample\u007bn1\u007d_mapped.tsv")
	output: dynamic(STARCODE_DIR + "sample\u007bn2\u007d_mapped_sc_out.tsv")
	shell: "bash " + SCRIPTS_DIR + "starcode_SM.sh"

rule make_starcode_input:
	message: "Processing input csv's before starcode"
	input: 
		md = "metaData.tsv",
		bcm = BARCODE_MAP_DIR + "finalBarcodeMap.csv",
		rc = dynamic(RAW_COUNTS + "\u007bn0\u007d.csv")
	output: dynamic(STARCODE_DIR + "sample\u007bn1\u007d_mapped.tsv")
	shell: "Rscript " + SCRIPTS_DIR + "pre_process_SM.R \u007binput.md\u007d \u007binput.bcm\u007d"
"""
)

if not merge:
    sf.write(
        f"""
rule count_barcodes:
	message: "Counting Barcodes"
	output: dynamic(RAW_COUNTS + "\u007bn0\u007d.csv")
	shell: "python3 " + SCRIPTS_DIR + "count_barcodes_SM.py {fastq_path}"
	"""
    )

if merge:
    sf.write(
        f"""

rule count_barcodes:
	message: "Counting Barcodes"
	input: "fastqs_joined.txt"
	output: dynamic(RAW_COUNTS + "\u007bn0\u007d.csv")
	run: 
	    shell("python3 " + SCRIPTS_DIR + "count_barcodes_SM.py file_info.csv")
	    shell("touch barcodes_counted.txt") 
	
rule join_fastqs:
	message: "Joining Fastq's"
	output: "fastqs_joined.txt"
	input: "fastqs_trimmed.txt"
	run: 
		shell("python3 " + SCRIPTS_DIR + "fastq_join_SM.py")
		shell("touch fastqs_joined.txt")

rule trim_fastqs:
	message: "Trimming Fastq's"
	input: "file_info.csv"
	output: "fastqs_trimmed.txt"
	run: 
		shell("python3 " + SCRIPTS_DIR + "fastq_trim_SM.py file_info.csv")
		shell("touch fastqs_trimmed.txt")
	"""
    )

sf.close()

print("Starting snakemake")
sys.stdout.flush()


command = f"snakemake -s runs/{dir_name}/Snakefile -d runs/{dir_name}/ -j"
os.system(command)

print("Empirical analysis finished!")
