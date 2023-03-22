import os
import sys

file_info = open(sys.argv[1], "r")

TRIMMED_DIR = "trimmed_files/"

for line in file_info:
	need_to_trim = line.split(",")[1].strip()
	file = line.split(",")[0].strip()

	if need_to_trim == "True":
		os.system(f"trim_galore --hardtrim5 24 {file} -o {TRIMMED_DIR} -j 4")