import os

MERGE_DIR = "joined_files/"
TRIM_DIR = "trimmed_files/"

f = open(f"{MERGE_DIR}matching_names.csv", "r+")
l = open(f"merge_stats.csv", "w+")


l.write("n_joined,n_un1,n_un2,name\n")

for line in f:
	file1, file2, joinName = line.split(",")
	joinName = joinName.strip()
	print(f"joining {file1} \nand {file2} \nto make {joinName}")
	# print(f"pear -f {TRIM_DIR + file1} -r {TRIM_DIR+ file2} -o {MERGE_DIR + joinName} -m 24 -n 24")
	os.system(f"fastq-join {TRIM_DIR + file1} {TRIM_DIR+ file2} -o {MERGE_DIR + joinName}")

	#Join_name = 19919X42_220722_A00421_0459_AHH3JFDRX2_S42_L001_R1-2_001_%.fq
	finalJoinName = joinName.replace("%", "join")
	un1Name = joinName.replace("%", "un1")
	un2Name = joinName.replace("%", "un2")

	# #Count the number of reads for each of the files
	join_count = os.popen(f'echo $(cat {MERGE_DIR}{finalJoinName} |wc -l)/4|bc').read().strip()
	un1_count = os.popen(f'echo $(cat {MERGE_DIR}{un1Name} |wc -l)/4|bc').read().strip()
	un2_count = os.popen(f'echo $(cat {MERGE_DIR}{un2Name} |wc -l)/4|bc').read().strip()
	l.write(f"{join_count},{un1_count},{un2_count},{finalJoinName}\n")

	os.system(f"mv {MERGE_DIR}{un1Name} {MERGE_DIR}un_files")
	os.system(f"mv {MERGE_DIR}{un2Name} {MERGE_DIR}un_files")

l.close()
f.close()

os.system(f"rm -r {TRIM_DIR}")

	

