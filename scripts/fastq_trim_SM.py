import os
import sys
import pandas as pd

"""
Loops through the rows of the file_info.csv (Written in TMP_empirical.py)
If the files need to be joined (I.E. They are paired end reads) they will be trimmed
Trimmed files are written to the trimmed_files/ directory
"""

file_info = pd.read_csv(sys.argv[1])
TRIMMED_DIR = "trimmed_files/"

for index, row in file_info.iterrows():
	need_to_trim = row["needs_to_join"]
	file = row["path_to_file"]

	if need_to_trim:
		os.system(f"trim_galore --hardtrim5 24 {file} -o {TRIMMED_DIR} -j 4")