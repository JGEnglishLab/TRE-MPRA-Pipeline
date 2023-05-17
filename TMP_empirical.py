import os
import sys
import re
import numpy as np
import argparse
import time as t
from datetime import datetime
import subprocess

print("Running TMP_empirical")


def check_dir(path): #Checks a directory to make sure that fastq files exist in directory
	for f in os.listdir(path):
		if f.endswith(".fastq") or f.endswith(".fq") or f.endswith(".fq.gz") or f.endswith(".fastq.gz"):
			return True
	return False

def print_new_section():
	print('\n')
	print("***********************************************************************")
	print("-----------------------------------------------------------------------")
	print("-----------------------------------------------------------------------")
	print("***********************************************************************")
	print("\n")

def check_y_n_inp(inp, correction_message = "Try again"):
	inp = inp.lower()
	while inp != "y" and inp != "n":
		print(f"{inp} is an invalid response")
		print(correction_message)
		inp = input("Enter \"y\" or \"n\": ")
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
	while not pattern.match(entered_samples) or entered_samples != pattern.match(entered_samples).group():
		print(f"\"{entered_samples}\" doesn't match the necessary pattern")
		entered_samples = input("Enter the sample numbers seperated by commas: ")
	return entered_samples


def check_sample_numbers(entered_samples, remaining_samples):
	"""
	Checks if the user accidentally enters a wronge sample number
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


def diff_letters(a,b):
	m = min([len(a), len(b)])
	s = sum ( a[i] != b[i] for i in range(m) )
	return s

def diff_reads(a,b):
	"""
	makes sure the only difference between two file names is R1 and R2
	"""
	pattern = re.compile(r"[r|R][1|2]")
	return pattern.search(a).group() != pattern.search(b).group()

def check_treatment_tsv(tsv):
	"""
	:param tsv: the tsv input by the user containing sample numbers and treatment types
	:return: a bool if correct or not, a treatment dictionary, a list of all sample numbers, and an error or success message
	"""
	treatments = {}
	sample_list = []
	for line in tsv:
		if line.isspace():
			continue
		try:
			sample_number  = line.split("\t")[0]
			treatment_type = line.split("\t")[1]
		except:
			return False, treatments, sample_list, "Treatment TSV must have two columns separated by a tab"

		if not sample_number.isnumeric():
			return False, treatments, sample_list,  "Treatment TSV must contain numeric sample numbers in first column"

		if sample_number in sample_list:
			return False, treatments, sample_list, "Treatment TSV may not contain the same sample number more than once."

		sample_list.append(sample_number)
		treatments[sample_number.strip()] = treatment_type.strip()
	return True, treatments, sample_list, "Treatment TSV looks good"

def check_DNA_tsv(tsv, dna_samples, rna_samples):
	"""
	:param tsv: the tsv input by the user containing DNA sample numbers and corresponding RNA sample Numbers
	:param dna_samples: All of the DNA sample numbers
	:param rna_samples: All of the RNA sample numbers
	:return:
	"""

	error_msg = ""
	for line in tsv:
		if line.isspace():
			continue
		try:
			DNA_num = line.split("\t")[0].strip()
			RNA_num = line.split("\t")[1].strip()
		except:
			error_msg = "Error in DNA TSV, must be tab separated."
			return False, error_msg

		if not DNA_num.isnumeric():
			error_msg = "Error in DNA TSV, DNA column (col 1) must only contain numbers."
			return False, error_msg
		if not RNA_num.isnumeric():
			error_msg = "Error in DNA TSV, RNA column (col 2) must only contain numbers."
			return False, error_msg
		if DNA_num not in dna_samples:
			error_msg = f"Error in DNA TSV, sample {DNA_num} found in DNA column (col 1). This isn't a DNA sample according to sample TSV"
			return False, error_msg
		if RNA_num not in rna_samples:
			error_msg = f"Error in DNA TSV, sample {RNA_num} found in RNA column (col 2). This isn't an RNA sample according to sample TSV"
			return False, error_msg

	return True, "DNA TSV looks good!"




###########################################################################
###########################################################################
#                   		GET ARGUMENTS
###########################################################################
###########################################################################
parser = argparse.ArgumentParser(
                    prog = 'Setup snakemake',
                    description = 'Sets up and writes a snakemake workflow to run MPRAnalyze')

parser.add_argument('-r', '--run_name', dest="dir_name", required=False)
parser.add_argument('-p', '--path_to_fq', dest="fastq_path", required=False)
# parser.add_argument('-j', '--join_fastq', dest="merge", required=False)
parser.add_argument('-t', '--path_to_treatment_tsv', dest="treatment_tsv_path", required=False)
parser.add_argument('-dt', '--path_to_dna_tsv', dest="dna_tsv_path", required=False)
parser.add_argument('-sr', '--sample_number_regex', dest="pattern", required=False)
parser.add_argument('-d', '--path_to_DNA_fastq', dest="dna_path", required=False)
parser.add_argument('-s', '--path_to_spikein_file', dest="spike_path", required=False)

#TODO finish these flags
parser.add_argument('-pr', '--read_pairing_regex', dest="read_pattern", required=False)

args = parser.parse_args()
dir_name=args.dir_name
fastq_path=args.fastq_path
# merge=args.merge
treatment_tsv_path=args.treatment_tsv_path
dna_tsv_path=args.dna_tsv_path
pattern=args.pattern
dna_path=args.dna_path
spike_path=args.spike_path

###########################################################################
###########################################################################
#                CHECK SPIKE-IN FILE IF SPECIFIED
###########################################################################
###########################################################################

#Will just be an empty list if no spike path is provided, or if all spikes are invalid
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
	nucleotides = 'ACTG'
	for spike in spike_file:
		spike = spike.strip()

		#Errors will be dealt with in the R Scripts make_sample_files_SM.R and make_mpra_input_SM.R
		#We just draw attention to the errors here
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
print_new_section()
wd = os.getcwd()
print(f"WORKING DIR {wd}")

if not dir_name:
	dir_name = input("Enter the name of the run: ")

while "_" in dir_name:
	print("The name of the run cannot include \"_\"")
	dir_name = input("Change the name of the run so it does not include any underscores: ")

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
if not treatment_tsv_path:
	print("\n")
	print("You are required to supply a .tsv containing the sample numbers and treatment names")
	print("The tsv should follow the following format")
	print("\n1\tSerum Free")
	print("2\tSerum Free")
	print("3\tSerum Free")
	print("4\tATP")
	print("5\tForskolin")
	print("6\tDNA\n")
	print("NOTES:")
	print("\t1: A DNA sample must be included for this pipeline to work. And it must be labeled as \"DNA\"")
	print("\t\tIf your run contains more than one DNA sample label all of them as \"DNA\"")
	print("\t2: Do not include a header.")
	print("\t3: If you have multiple replicates of the same treatment they must be the same name")
	print("\tfor example if you put")
	print("\t\t1\tSerum Free 1")
	print("\t\t2\tSerum Free 2")
	print("\t\t3\tSerum Free 3")
	print("\t Serum Free 1 2 and 3 would all be treated as separate treatment types. ")
	print("\n")

	treatment_tsv_path = input("\nEnter path to tsv containing sample numbers and corresponding treatments: ")

legal_path = os.path.exists(treatment_tsv_path) and bool(re.search(r".tsv(/)?$", treatment_tsv_path))
while not legal_path:
	print(f"tsv not found in {treatment_tsv_path}")
	treatment_tsv_path = input("Check path and enter it again: ")
	legal_path = os.path.exists(treatment_tsv_path) and bool(re.search(r".tsv(/)?$", treatment_tsv_path))
print("tsv found!")

treatment_tsv = open(treatment_tsv_path, "r")
correct_tsv, treatments, tsv_sample_numbers, msg = check_treatment_tsv(treatment_tsv)
while not correct_tsv:
	treatment_tsv.close()
	print("Error with treatment TSV provided")
	print(msg)
	_ = input("Check format and and try again (Enter \"y\" to continue)")
	treatment_tsv = open(treatment_tsv_path, "r")
	correct_tsv, treatments, tsv_sample_numbers, msg = check_treatment_tsv(treatment_tsv)

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

if noDNAProvided: #No DNA sample was auto detected in the TSV!
	print("No DNA sample was auto detected in the TSV!")
	print("The treatment TSV must contain at least 1 DNA sample labeled as \"DNA\"")
	print("If more than one DNA sample is present for run, label all of them as \"DNA\"")
	exit()
	#TODO maybe Make sure that the user inputs DNA

if len(DNA_sample_num) > 1: #Multiple DNA sample detected.
	print_new_section()
	print("Multiple DNA samples detected!")
	print("You must provide a TSV containing the DNA sample numbers in the first column")
	print("and the corresponding RNA samples in the second column.")
	print("For example if you had\n")
	print("\t1\tSerum Free")
	print("\t2\tSerum Free")
	print("\t3\tATP")
	print("\t4\tForskolin")
	print("\t5\tDNA")
	print("\t6\tDNA\n")
	print("And DNA 5 corresponded with the Serum Free samples and DNA 6 corresponded with the others")
	print("You would input a TSV of the following format\n")
	print("\t5\t1")
	print("\t5\t2")
	print("\t6\t3")
	print("\t6\t4\n")
	print("NOTES:")
	print("\t1: Do not include a header")
	print("\t2: Each RNA sample must correspond to 1 and only 1 DNA sample")
	if not dna_tsv_path:
		dna_tsv_path = input("Provide the path to the a TSV containing the DNA sample numbers and corresponding RNA sample numbers: ")

	legal_path = os.path.exists(dna_tsv_path) and bool(re.search(r".tsv(/)?$", dna_tsv_path))
	while not legal_path:
		print(f"tsv not found in {dna_tsv_path}")
		dna_tsv_path = input("Check path and enter it again: ")
		legal_path = os.path.exists(dna_tsv_path) and bool(re.search(r".tsv(/)?$", dna_tsv_path))
	print("tsv found!")

	DNA_tsv = open(dna_tsv_path, "r")
	correct_tsv, error_msg = check_DNA_tsv(DNA_tsv, DNA_sample_num, RNA_sample_num)
	while not correct_tsv:
		DNA_tsv.close()
		print(error_msg)
		_ = input("Check format and and try again (Enter \"y\" to continue)")
		DNA_tsv = open(dna_tsv_path, "r")
		correct_tsv, error_msg = check_DNA_tsv(DNA_tsv, DNA_sample_num, RNA_sample_num)
	print(error_msg)
	# abs_dna_tsv_path = os.system(f"readlink -f {dna_tsv_path}")
	# abs_dna_tsv_path = subprocess.Popen(f"readlink -f {dna_tsv_path}")
	abs_dna_tsv_path = os.popen(f"readlink -f {dna_tsv_path}").read().strip()

else:
	abs_dna_tsv_path = None




###########################################################################
###########################################################################
#                   	GET PATH TO FASTQS
###########################################################################
###########################################################################
print_new_section()
if not fastq_path:
	fastq_path = input("Enter the path to your fastq files: ")

legal_path = os.path.exists(fastq_path) and check_dir(fastq_path)

while not legal_path:
	if not os.path.exists(fastq_path):
		print(f"Path does not exist: \"{fastq_path}\"")
		fastq_path = input("Check path and enter it again: ")
	if os.path.exists(fastq_path) and not check_dir(fastq_path):
		print(f"No Fastq files found in path: \"{fastq_path}\"")
		fastq_path = input("Check path and enter it again: ")
	legal_path = os.path.exists(fastq_path) and check_dir(fastq_path)

#Print out the files
print(f"Files found in {fastq_path}")
for f in os.listdir(fastq_path):
	if bool(re.search(r".f(q|astq)(.gz)?$", f)):
		print(f)

#Get the absolute path to the directory
os.chdir(fastq_path)
fastq_path = os.getcwd()
if not fastq_path.endswith("/"):
	fastq_path = fastq_path + "/"
#Then change pack to the correct dir
os.chdir(wd)

if not dna_path:
	print("\nNOTE!!!\nIf the DNA Fastq files are not contained in the directory above, you must provide a path to them using the \"-d\" flag")
else:
	legal_path = os.path.exists(dna_path) and check_dir(dna_path)

	while not legal_path:
		if not os.path.exists(dna_path):
			print(f"Path to DNA samples, \"{dna_path}\" does not exist ")
			dna_path = input("Check path and enter it again: ")
		if os.path.exists(dna_path) and not check_dir(dna_path):
			print(f"No Fastq files found in path to DNA samples \"{dna_path}\"")
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
os.system(f"mkdir -p ./runs/{dir_name}/mpra_input")
os.system(f"mkdir -p ./runs/{dir_name}/raw_counts")
os.system(f"mkdir -p ./runs/{dir_name}/rna_dna_samples")
os.system(f"mkdir -p ./runs/{dir_name}/run_descriptive_stats")
os.system(f"mkdir -p ./runs/{dir_name}/star_code")
#Used for files that need to be merged
os.system(f"mkdir -p ./runs/{dir_name}/trimmed_files")
os.system(f"mkdir -p ./runs/{dir_name}/joined_files")
os.system(f"mkdir -p ./runs/{dir_name}/joined_files/un_files")

print_new_section()

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
	if bool(re.search(r".f(q|astq)(.gz)?$", filename)):
		files.append(filename)
		raw_file_paths[filename] = fastq_path

if dna_path and os.path.exists(dna_path):
	for filename in os.listdir(dna_path):
		if bool(re.search(r".f(q|astq)(.gz)?$", filename)):
			files.append(filename)
			raw_file_paths[filename] = dna_path

#If the files need to be merged, the names will look slightly different

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
count = 0
max_iter = len(files)
#These two will be used for meta_data.csv
#They will contain the altered final names
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
	if file1 == "run1X42_220722_A00421_0459_AHH3JFDRX2_S42_L001_R1_001_subset.fastq.gz":
		print("STOP!")
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
			file1 = file1.split('.', 1)[0] + ".24bp_5prime." + file1.split('.', 1)[1]
			file2 = file2.split('.', 1)[0] + ".24bp_5prime." + file2.split('.', 1)[1]
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
	files.remove(file1)
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
	files.remove(non_match_file)
	non_match_file = re.sub(".f(q|astq)(.gz)?", ".csv", non_match_file)
	un_paired_files.append(non_match_file)


f.close()
file_info.close()

if len(un_paired_files) == 0:
	print("\nAuto detected all pairs!")
	merge = True

elif len(file1) == len(files):
	print("\nNo pairs were auto detected!")
	print("In order to pair files together the only difference in the file names should be between \"R1\" and \"R2\"")
	print("The patterns \"R1\" and \"R2\" may only be present once in each file name.")
	print("Run will continue without pairing files")
	print("If files must be paired, fix file names and start over.")
	merge = False

else:
	print("\nSome files were not paired")
	print("In order to pair files together the only difference in the file names should be between \"R1\" and \"R2\"")
	print("The patterns \"R1\" and \"R2\" may only be present once in each file name.")
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

print_new_section()

files = res_list = [y for x in [paired_files, un_paired_files] for y in x]

fileSamplePair = {}
all_sample_numbers = []

if not pattern:
	pattern = "S"
	print("Each file name will be searched for a sample number automatically search for sample numbers.")
	print("The default is to search for a number following \"S\" IE S{number}")
	print("You can change this using the -sr flag")
	print("If you want to change it use -sr {pattern}")
	print("Enter the exact pattern that precedes the sample number in each of the file names.")
	print("For example. If your file names looked like this:")
	print("\ttest_lane1_sn1.fastq\n\ttest_lane2_sn1.fastq\n\ttest_lane1_sn2.fastq\n\ttest_lane2_sn2.fastq\n\ttest_lane1_sn3.fastq\n\ttest_lane2_sn3.fastq\n")
	print("and the sample number is indicated by the number following \"_sn\"")
	print("You would enter \"_sn\" as the pattern")

print_new_section()

for f in files:
	reg_res = re.findall(pattern + '\d+', f)
	if check_regex(reg_res, pattern):
		sample_number = reg_res[0].replace(pattern, '')  #Get the actual sample number
		fileSamplePair[f] = sample_number
		if sample_number not in all_sample_numbers:
			all_sample_numbers.append(sample_number)
	else:
		print("\n****************************************")
		print("Different file format detected for \"" + f + "\"")
		sample_number = input("Enter sample number e.g. \"14\" for sample 14: ")
		while not sample_number.isnumeric():
			print("Sample number must be a number.")
			input(f"Enter sample number for {f}")
		fileSamplePair[f] = sample_number
		if sample_number not in all_sample_numbers:
			all_sample_numbers.append(sample_number)

#Check to see if sample number detected from files and sample numbers from treatment TSV match
tsv_sample_numbers.sort()
all_sample_numbers.sort()
if tsv_sample_numbers != all_sample_numbers:
	print("The sample numbers provided in treatment TSV and the sample numbers detected from files do not match!\n")
	print("These are the sample numbers from treatment TSV that are NOT detected in the files provided")
	print("-")
	for i in tsv_sample_numbers:
		if i not in all_sample_numbers:
			print("\t", i)
	print("-")
	print("These are the sample numbers from the files provided that are NOT in the treatment TSV")
	print("-")
	for i in all_sample_numbers:
		if i not in tsv_sample_numbers:
			print("\t", i)
	print("-")
	print("Only these samples numbers will be used")
	joined_sample_numbers = list(set(all_sample_numbers) & set(tsv_sample_numbers))
	joined_sample_numbers.sort()
	print(joined_sample_numbers)

	con = input("If you only want to use the above files enter y, otherwise enter n to exit and fixe your files (y/n): ")
	con = check_y_n_inp(con)

	if con == "n":
		exit()

	for file in fileSamplePair:
		if fileSamplePair[file] not in joined_sample_numbers:
			del fileSamplePair[file]

print_new_section()


f = open(f"./runs/{dir_name}/metaData.csv", "w+")
f.write("fileName,sampleNumber,treatment,run_name\n")

#Write CSV
for file in fileSamplePair:
	line = file + "," + fileSamplePair[file] + "," + treatments[fileSamplePair[file]] + "," + dir_name + "\n"
	f.write(line)

f.close()


###########################################################################
###########################################################################
#                   	WRITE SNAKE MAKE FILE START
###########################################################################
###########################################################################


sf = open(f"./runs/{dir_name}/Snakefile", "w+")

if merge:
	sf.write("MERGE_DIR = \"merged_files/\"")

sf.write(f"""
STARCODE_DIR = "star_code/"
SCRIPTS_DIR = "../../scripts/"
BARCODE_MAP_DIR = "../../barcode_map_data/"
DESCRIPTIVE_STATS_DIR = "run_descriptive_stats/"
RNA_DNA_DIR = "rna_dna_samples/"
RAW_COUNTS = "raw_counts/"
MPRA_INPUT_DIR = "mpra_input/"

rule create_alphas:
	message: "Running MPRAnalyze to create alpha files"
	input: 
		rc = MPRA_INPUT_DIR + "rna_counts.csv",
		dc = MPRA_INPUT_DIR + "dna_counts.csv",
		rdv = MPRA_INPUT_DIR + "rna_depth_vals.csv",
		ddv = MPRA_INPUT_DIR + "dna_depth_vals.csv",
		ti = MPRA_INPUT_DIR + "treatment_id.csv",
		ca = MPRA_INPUT_DIR + "col_annotations.csv",
	output: "empirical_results.csv"
	shell: "Rscript " + SCRIPTS_DIR + "run_empirical_SM.R \u007binput.rc\u007d \u007binput.dc\u007d \u007binput.ti\u007d \u007binput.ca\u007d \u007binput.rdv\u007d \u007binput.ddv\u007d"

rule make_mpra_input:
	message: "Creating input for MPRAnalyze"
	input: 
		ad = RNA_DNA_DIR + "all_data_filtered.csv", 
		ds = RNA_DNA_DIR + "dna_samples.csv", 
		rs = RNA_DNA_DIR + "rna_samples.csv", 
		md = "metaData.csv",
		bcm = BARCODE_MAP_DIR + "finalBarcodeMap.csv"
	output:
		MPRA_INPUT_DIR + "dna_depth.csv",
		MPRA_INPUT_DIR + "rna_counts.csv",
		MPRA_INPUT_DIR + "rna_depth_vals.csv",
		MPRA_INPUT_DIR + "dna_depth_vals.csv",
		MPRA_INPUT_DIR + "dna_counts.csv",
		MPRA_INPUT_DIR + "treatment_id.csv",
		MPRA_INPUT_DIR + "rna_depth.csv",
		MPRA_INPUT_DIR + "col_annotations.csv"
	shell: "Rscript " + SCRIPTS_DIR + "make_mpra_input_SM.R \u007binput.bcm\u007d \u007binput.ad\u007d \u007binput.rs\u007d \u007binput.ds\u007d \u007binput.md\u007d {abs_dna_tsv_path} {abs_spike_path}"

rule make_dna_and_rna_samples:
	message: "Creating images about run and generating dna and rna samples"
	input: 
		sc = dynamic(STARCODE_DIR + "analyzed_out_sample\u007bn3\u007d_mapped_sc_out.tsv"),
		md = "metaData.csv",
		bcm = BARCODE_MAP_DIR + "finalBarcodeMap.csv"
	output: 
		DESCRIPTIVE_STATS_DIR + "filtering_ratios.png",
		DESCRIPTIVE_STATS_DIR + "type_ratios.png",
		DESCRIPTIVE_STATS_DIR + "spike_in_stats.png",
		DESCRIPTIVE_STATS_DIR + "dna_per_barcode.png",
		DESCRIPTIVE_STATS_DIR + "pre_filtering_barcodes.csv",
		DESCRIPTIVE_STATS_DIR + "spike_in_stats.csv",
		RNA_DNA_DIR + "all_data_filtered.csv", 
		RNA_DNA_DIR + "dna_samples.csv", 
		RNA_DNA_DIR + "rna_samples.csv", 
	shell: "Rscript " + SCRIPTS_DIR + "make_sample_files_SM.R \u007binput.md\u007d \u007binput.bcm\u007d {abs_spike_path} \u007binput.sc\u007d"

rule analyze_starcode:
	message: "Analyzing Starcode"
	input: dynamic(STARCODE_DIR + "sample\u007bn2\u007d_mapped_sc_out.tsv")
	output: dynamic(STARCODE_DIR + "analyzed_out_sample\u007bn3\u007d_mapped_sc_out.tsv")
	shell: "python " + SCRIPTS_DIR + "analyze_star_code_SM.py \u007binput\u007d"

rule run_starcode:
	message: "Running Star Code"
	input: dynamic(STARCODE_DIR + "sample\u007bn1\u007d_mapped.tsv")
	output: dynamic(STARCODE_DIR + "sample\u007bn2\u007d_mapped_sc_out.tsv")
	shell: "bash " + SCRIPTS_DIR + "starcode_SM.sh"

rule make_starcode_input:
	message: "Processing input csv's before starcode"
	input: 
		md = "metaData.csv",
		bcm = BARCODE_MAP_DIR + "finalBarcodeMap.csv",
		rc = dynamic(RAW_COUNTS + "\u007bn0\u007d.csv")
	output: dynamic(STARCODE_DIR + "sample\u007bn1\u007d_mapped.tsv")
	shell: "Rscript " + SCRIPTS_DIR + "pre_process_SM.R \u007binput.md\u007d \u007binput.bcm\u007d"
""")

if not merge:
	sf.write(f"""
rule count_barcodes:
	message: "Counting Barcodes"
	output: dynamic(RAW_COUNTS + "\u007bn0\u007d.csv")
	shell: "python " + SCRIPTS_DIR + "count_barcodes_SM.py {fastq_path}"
	""")

if merge:
	sf.write(f"""
rule count_barcodes:
	message: "Counting Barcodes"
	input: "fastqs_joined.txt"
	output: dynamic(RAW_COUNTS + "\u007bn0\u007d.csv")
	shell: "python " + SCRIPTS_DIR + "count_barcodes_SM.py file_info.csv" 
	
rule join_fastqs:
	message: "Joining Fastq's"
	output: "fastqs_joined.txt"
	input: "fastqs_trimmed.txt"
	run: 
		shell("python " + SCRIPTS_DIR + "fastq_join_SM.py")
		shell("touch fastqs_joined.txt")

rule trim_fastqs:
	message: "Trimming Fastq's"
	input: "file_info.csv"
	output: "fastqs_trimmed.txt"
	run: 
		shell("python " + SCRIPTS_DIR + "fastq_trim_SM.py file_info.csv")
		shell("touch fastqs_trimmed.txt")
	""")

sf.close()


print_new_section()

start = input("Would you like to start the run now? (y/n): ")
start = check_y_n_inp(start)

if start == "y":
	sbatch = input("Would you like to run this with SBATCH? (y/n): ")
	sbatch = check_y_n_inp(sbatch)
	if sbatch == "y":
		print("STILL NEED TO ADD CODE FOR RUNNING WITH SBATCH")
	elif sbatch == "n":
		cores = input("Enter the number of Cores to be used: ")
		while not cores.isnumeric():
			print("Must be a number")
			cores = input("Enter the number of Cores to be used: ")
		command = f"snakemake -s runs/{dir_name}/Snakefile -d runs/{dir_name}/ -j8"
		os.system(command)




#TODO
# Add regex for r[2|1] (line 459)
# Check multi dna
# Fix R scripts to use multi DNA (you might need to write it to a tsv or just use user tsv as input)
# Spike in stuff (Change R script too)
# Check to see if the entered sample numbers (TSV) match the found sample numbers (looping through files/regex)




