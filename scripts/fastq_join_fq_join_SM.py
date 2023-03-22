import os

MERGE_DIR = "merged_files/"
TRIM_DIR = "trimmed_files/"

f = open(f"{MERGE_DIR}matching_names.csv", "r+")
# l = open(f"{MERGE_DIR}merge_stats.csv", "w+")

# l.write("n_joined,n_un1,n_un2,n_discarded,name\n")

for line in f:
	file1, file2, joinName = line.split(",")
	joinName = joinName.strip()
	print(f"joining {file1} \nand {file2} \nto make {joinName}")
	print(f"fastq-join {TRIM_DIR + file1} {TRIM_DIR+ file2} -o {MERGE_DIR + joinName}")
	os.system(f"fastq-join {TRIM_DIR + file1} {TRIM_DIR+ file2} -o {MERGE_DIR + joinName}")

	#Count the number of reads for each of the files
	# paired_count = os.popen(f'echo $(cat {MERGE_DIR}{joinName}.assembled.fastq |wc -l)/4|bc').read().strip()
	# un_paired_count_1 = os.popen(f'echo $(cat {MERGE_DIR}*unassembled.forward.fastq |wc -l)/4|bc').read().strip()
	# un_paired_count_2 = os.popen(f'echo $(cat {MERGE_DIR}*unassembled.reverse.fastq |wc -l)/4|bc').read().strip()
	# discarded_count = os.popen(f'echo $(cat {MERGE_DIR}*discarded.fastq |wc -l)/4|bc').read().strip()
	# l.write(f"{paired_count},{un_paired_count_1},{un_paired_count_2},{discarded_count},{joinName}\n")

	# os.system(f"mv {MERGE_DIR}*unassembled* {MERGE_DIR}un_files")
	# os.system(f"mv {MERGE_DIR}*discarded* {MERGE_DIR}un_files")

# l.close()
# f.close()
	

