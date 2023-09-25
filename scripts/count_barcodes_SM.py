from Bio import SeqIO
import sys
import os
import subprocess
import re
MERGED_DIR = "joined_files/"

file_info = open(sys.argv[1], "r")
files_to_count = []

#Add files that didn't need to be trimmed and joined
for line in file_info:
	need_to_trim = line.split(",")[1].strip()
	file = line.split(",")[0].strip()
	if need_to_trim == "False":
		files_to_count.append(file)

#Add files that did need to be trimmed and joined
for file in os.listdir(MERGED_DIR):
	if bool(re.search(r".f(q|astq)(.gz)?$", file)):
		files_to_count.append(MERGED_DIR + file)

#Count all fastq files
for file in files_to_count:
	file_name = file.split("/")[-1]


	if file.endswith(".gz"):
		zipped = True
	else:
		zipped = False

	curDict = {}

	if zipped:
		process = subprocess.Popen(("gunzip","-c",file), stdout=subprocess.PIPE,text=True)
		curRecord = SeqIO.parse(process.stdout, "fastq")
	else:
		curRecord = SeqIO.parse(file, "fastq")

	print(f"counting barcodes of {file_name}")
	for i in curRecord:		
		curSeq = i.seq[0:24]
		try:
			curDict[curSeq] += 1
		except KeyError:
			curDict[curSeq] = 1

	file_name = file.split("/")[-1]
	csvName = re.split(".fastq|.fq", file_name)[0] + ".csv"
	csvName = "./raw_counts/" + csvName
	f = open(csvName, "w+")
	for key in curDict:
		line = str(key) + "," +str(curDict[key])+ "\n"
		f.write(line)
	f.close()

os.system(f"rm -r {MERGED_DIR}")







	
	

	



