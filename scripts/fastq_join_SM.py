import os
import pandas as pd

"""
Reads names file "matching_names.csv" (Written in TMP_empirical.py)
Names file contains file1_name, file2_name, and joine_name
File1 and 2 look like this.
run1X8_220722_A00421_0459_AHH3JFDRX2_S8_L002_R1_001.24bp_5prime.fq
The .24bp_fprime.fq extension comes from trim-galore

The joined name looks like this
run1X8_220722_A00421_0459_AHH3JFDRX2_S8_L002_R1-2_001_%.fq
The % is replaced automatically by fastq-join with either 'join' 'un1' or 'un2'
The 'join' file contains the joined fastq reads
The 'un1' and 'un2' contain the unpaired reads

Runs fastq-join on read1 and read2 with default settings

Reads the number of reads in each file (joined, un1, and un2) so the user can see how many reads are being joined. 
"""

MERGE_DIR = "joined_files/"
TRIM_DIR = "trimmed_files/"

names_file = pd.read_csv(f"{MERGE_DIR}matching_names.csv")
merge_stats = open(f"merge_stats.csv", "w+")

merge_stats.write("n_joined,n_un1,n_un2,name\n")

for index, row in names_file.iterrows():
	file1 = row["file1_name"]
	file2 = row["file2_name"]
	joinName = row["joined_name"]

	print(f"joining {file1} \nand {file2} \nto make {joinName}")
	os.system(f"fastq-join {TRIM_DIR + file1} {TRIM_DIR+ file2} -o {MERGE_DIR + joinName}")

	finalJoinName = joinName.replace("%", "join")
	un1Name = joinName.replace("%", "un1")
	un2Name = joinName.replace("%", "un2")

	# Count the number of reads for each of the files
	join_count = os.popen(f'echo $(cat {MERGE_DIR}{finalJoinName} |wc -l)/4|bc').read().strip()
	un1_count = os.popen(f'echo $(cat {MERGE_DIR}{un1Name} |wc -l)/4|bc').read().strip()
	un2_count = os.popen(f'echo $(cat {MERGE_DIR}{un2Name} |wc -l)/4|bc').read().strip()
	merge_stats.write(f"{join_count},{un1_count},{un2_count},{finalJoinName}\n")

	os.system(f"mv {MERGE_DIR}{un1Name} {MERGE_DIR}un_files")
	os.system(f"mv {MERGE_DIR}{un2Name} {MERGE_DIR}un_files")

merge_stats.close()
# Remove all the trimmed fastq files
os.system(f"rm -r {TRIM_DIR}")

	

