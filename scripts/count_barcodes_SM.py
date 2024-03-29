from Bio import SeqIO
import sys
import os
import subprocess
import re
import pandas as pd

"""
Loops through the reads of every fastq file.
Counts the barcodes for each
Every fastq will hav a csv written that will look like this. 

AGC....,251
GCT....,52
"""

MERGED_DIR = "joined_files/"

file_info = pd.read_csv(sys.argv[1])
files_to_count = [] #Add all the file paths the will be counted

#Add files that didn't need to be trimmed and joined
for index, row in file_info.iterrows():
	need_to_trim = row["needs_to_join"]
	file = row["path_to_file"]
	if not need_to_trim:
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

	#Loop through every sequence.
	#If sequence already in dictionary add 1 to its count
	#Otherwise add it to the dictionary at a count of 1
	print(f"counting barcodes of {file_name}")
	for i in curRecord:		
		curSeq = i.seq[0:24]
		try:
			curDict[curSeq] += 1
		except KeyError:
			curDict[curSeq] = 1

	csvName = re.split(".fastq|.fq", file_name)[0] + ".csv" #Replace that .fq/.fastq with ".csv"
	csvName = "./raw_counts/" + csvName
	f = open(csvName, "w+")
	for key in curDict:
		line = str(key) + "," +str(curDict[key])+ "\n"
		f.write(line)
	f.close()

#Removed the merged fastq files
os.system(f"rm -r {MERGED_DIR}")







	
	

	



