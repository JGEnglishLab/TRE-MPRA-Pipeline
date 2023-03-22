import os
import statistics
import sys

#Read in barcode map
#TODO consider passing the barcode map as snake make input
bc_map = "finalBarcodeMap.csv"
bc_map = "../../barcode_map_data/" + bc_map
bc_map_file = open(bc_map, "r+")
bc_dict = {}
for line in bc_map_file:
	if not line.startswith("id"):
		line = line.split(",")
		bc = line[6]
		architecture = line[0] + line[1] + line[2] + line[3] + line[4] + line[5]
		bc_dict[bc] = architecture

for starCodeOutputName in sys.argv:

	#*************************************
	#	Get all file names configured
	#*************************************


	#Starcode name will be in this format
	#star_code/sample4_mapped_sc_out.tsv
	if starCodeOutputName.endswith(".py"):
		continue
	print(starCodeOutputName)
	sc_out_file = open(starCodeOutputName, "r+")

	#Need output from pre_process_SM.R (pre starcode)
	#Will be in this format
	#star_code/sample4_mapped.tsv
	preStarcodeName = starCodeOutputName.replace("_sc_out.tsv", ".tsv")
	print(preStarcodeName)

	#Get the output name
	outputName = starCodeOutputName.replace("sample", "analyzed_out_sample")
	print(outputName)

	#*********************************************
	#	Make dictionaries on pre-starcode file
	#*********************************************

	pre_sc_file = open(preStarcodeName, "r+")
	
	curCountDict = {}
	curMappedDict = {}

	for line in pre_sc_file:
		line = line.split()
		bc = line[0]
		count = line[1]
		mapped = line[2]
		curCountDict[bc] = count
		curMappedDict[bc] = mapped
	pre_sc_file.close()

	#************************************
	#	Loop through starcode output
	#************************************

	mean_v = []
	median_v = []
	centroid_v = []
	totalCollapsed_v = []
	totalCount_v = []
	centroidCount_v = []
	centroidInMapped_v = []
	collapseBarCodes_v = []
	collapsedInMapped_v = []
	collapseArchMatchesCentroid_v = []

	for line in sc_out_file:
		line = line.split()
		allBarCodes = line[2]
		allBarCodes = allBarCodes.split(",")

		centroid = line[0]
		totalCount = line[1]
		totalCollapsed = len(allBarCodes)
		centroidCount = curCountDict[centroid]
		centroidInMapped = curMappedDict[centroid]

		collapsedCounts = []
		collapsedInMapped = []
		collapsedBarCodes = []
		for bc in allBarCodes:
			collapsedCounts.append(int(curCountDict[bc]))
			if bc != centroid:
				collapsedBarCodes.append(bc)
				if curMappedDict[bc] == "TRUE":
					collapsedInMapped.append(bc)

		
		if len(collapsedInMapped) == 0:
			collapsedInMapped.append("NULL")

		if len(collapsedBarCodes) == 0:
			collapsedBarCodes.append("NULL")
			collapseArchMatchesCentroid_v.append("Ignore")
		elif centroidInMapped == "TRUE":
			centroid_arch = bc_dict[centroid]
			collapse_same_arch = True

			if collapsedInMapped[0] != "NULL":
				for bc in collapsedInMapped:
					if bc_dict[bc] != centroid_arch:
						collapse_same_arch = False
				if collapse_same_arch:
					collapseArchMatchesCentroid_v.append("TRUE")
				else:
					collapseArchMatchesCentroid_v.append("FALSE")
			else:
				collapseArchMatchesCentroid_v.append("Ignore")
		else:
			collapseArchMatchesCentroid_v.append("Ignore")


		mean = statistics.mean(collapsedCounts)
		median = statistics.median(collapsedCounts)

		mean_v.append(mean)
		median_v.append(median)
		centroid_v.append(centroid)
		totalCollapsed_v.append(totalCollapsed)
		totalCount_v.append(totalCount)
		centroidCount_v.append(centroidCount)
		centroidInMapped_v.append(centroidInMapped)
		collapseBarCodes_v.append(",".join(collapsedBarCodes))
		collapsedInMapped_v.append(",".join(collapsedInMapped))

	sc_out_file.close()

	outPut = open(outputName, "w+")
	outPut.write("means\tmedians\ttotalCollapsed\ttotalCounts\tcentroid\tcentroidCounts\tcentroidInMapped\tcollapsedBarcodes\tarchitecturesMatch\tcollapsedInMapped\n")
	for i in range(len(mean_v)):
		newLine = (str(mean_v[i]) + "\t" +
					str(median_v[i]) + "\t" + 
					str(totalCollapsed_v[i]) + "\t" +
					str(totalCount_v[i]) + "\t" + 
					centroid_v[i] + "\t" + 
					centroidCount_v[i] + "\t" + 
					centroidInMapped_v[i] + "\t" + 
					collapseBarCodes_v[i] + "\t"+ 
					collapseArchMatchesCentroid_v[i] + "\t"+ 
					collapsedInMapped_v[i] + "\n")

		outPut.write(newLine)

	outPut.close()
					
