import os
import statistics
import sys

#Read in barcode map
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
	#star_code/sample4_map_sc_out.tsv
	if starCodeOutputName.endswith(".py"):
		continue
	print(starCodeOutputName)
	sc_out_file = open(starCodeOutputName, "r+")

	#Need output from pre_process_SM.R (pre starcode)
	#Will be in this format
	#star_code/sample4_map.tsv
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
	curMapDict = {}

	for line in pre_sc_file:
		line = line.split()
		bc = line[0]
		count = line[1]
		map = line[2]
		curCountDict[bc] = count
		curMapDict[bc] = map
	pre_sc_file.close()

	#************************************
	#	Loop through starcode output
	#************************************

	mean_v = []
	median_v = []
	centroid_v = []
	totalClustered_v = []
	totalCount_v = []
	centroidCount_v = []
	centroidInMap_v = []
	clusteredBarCodes_v = []
	clusteredInMap_v = []
	clusteredArchMatchesCentroid_v = []

	for line in sc_out_file:
		line = line.split()
		allBarCodes = line[2]
		allBarCodes = allBarCodes.split(",")

		centroid = line[0]
		totalCount = line[1]
		totalClustered = len(allBarCodes)
		centroidCount = curCountDict[centroid]
		centroidInMap = curMapDict[centroid]

		clusteredCounts = [] #Keep track of the total counts of all clustered barcodes (not the
		clusteredInMap = [] #Keep track of all the clustered barcodes (not the centroid) that are also in the map
		clusteredBarCodes = [] #Keep track of all clustered barcodes (not the centroid)
		for bc in allBarCodes:
			clusteredCounts.append(int(curCountDict[bc]))
			if bc != centroid:
				clusteredBarCodes.append(bc)
				if curMapDict[bc] == "TRUE":
					clusteredInMap.append(bc)

		
		if len(clusteredInMap) == 0:
			clusteredInMap.append("NULL")

		if len(clusteredBarCodes) == 0:
			clusteredBarCodes.append("NULL")
			clusteredArchMatchesCentroid_v.append("NULL")
		elif centroidInMap == "TRUE":
			centroid_arch = bc_dict[centroid]
			clustered_same_arch = True

			if clusteredInMap[0] != "NULL":
				for bc in clusteredInMap:
					if bc_dict[bc] != centroid_arch:
						clustered_same_arch = False
				if clustered_same_arch:
					clusteredArchMatchesCentroid_v.append("TRUE")
				else:
					clusteredArchMatchesCentroid_v.append("FALSE")
			else:
				clusteredArchMatchesCentroid_v.append("NULL")
		else:
			clusteredArchMatchesCentroid_v.append("NULL")


		mean = statistics.mean(clusteredCounts)
		median = statistics.median(clusteredCounts)

		mean_v.append(mean)
		median_v.append(median)
		centroid_v.append(centroid)
		totalClustered_v.append(totalClustered)
		totalCount_v.append(totalCount)
		centroidCount_v.append(centroidCount)
		centroidInMap_v.append(centroidInMap)
		clusteredBarCodes_v.append(",".join(clusteredBarCodes))
		clusteredInMap_v.append(",".join(clusteredInMap))

	sc_out_file.close()

	outPut = open(outputName, "w+")
	outPut.write("mean_cluster_count\tmedian_cluster_count\ttotal_clustered\ttotal_counts\tcentroid\tcentroid_counts\tcentroid_in_map\tclustered_barcodes\tarchitectures_match\tclustered_in_map\n")
	for i in range(len(mean_v)):
		newLine = (str(mean_v[i]) + "\t" +
					str(median_v[i]) + "\t" + 
					str(totalClustered_v[i]) + "\t" +
					str(totalCount_v[i]) + "\t" + 
					centroid_v[i] + "\t" + 
					centroidCount_v[i] + "\t" + 
					centroidInMap_v[i] + "\t" +
					clusteredBarCodes_v[i] + "\t"+
					clusteredArchMatchesCentroid_v[i] + "\t"+
					clusteredInMap_v[i] + "\n")

		outPut.write(newLine)

	outPut.close()
					
