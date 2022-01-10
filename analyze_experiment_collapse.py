import os
import statistics

for fileName in os.listdir("."):
	if fileName.endswith("_mapped_bcm2.tsv") and not fileName.startswith("sc"):
		scName = "sc_out_" + fileName
		finalName = "analyzed_out_" + fileName

		curFile = open(fileName, "r+")
		
		curCountDict = {}
		curMappedDict = {}

		#create dictionaries
		for line in curFile:
			line = line.split()
			bc = line[0]
			count = line[1]
			mapped = line[2]
			curCountDict[bc] = count
			curMappedDict[bc] = mapped
		curFile.close()

		starCodeOutput = open(scName, "r+")

		mean_v = []
		median_v = []
		centroid_v = []
		totalCollapsed_v = []
		totalCount_v = []
		centroidCount_v = []
		centroidInMapped_v = []
		collapseBarCodes_v = []
		collapsedInMapped_v = []

		for line in starCodeOutput:
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

			if len(collapsedBarCodes) == 0:
				collapsedBarCodes.append("NULL")
			if len(collapsedInMapped) == 0:
				collapsedInMapped.append("NULL")

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

		starCodeOutput.close()

		outPut = open(finalName, "w+")
		outPut.write("means\tmedians\ttotalCollapsed\ttotalCounts\tcentroid\tcentroidCounts\tcentroidInMapped\tcollapsedBarCodes\tcollapsedInMapped\n")
		for i in range(len(mean_v)):
			newLine = (str(mean_v[i]) + "\t" +
						str(median_v[i]) + "\t" + 
						str(totalCollapsed_v[i]) + "\t" +
						str(totalCount_v[i]) + "\t" + 
						centroid_v[i] + "\t" + 
						centroidCount_v[i] + "\t" + 
						centroidInMapped_v[i] + "\t" + 
						collapseBarCodes_v[i] + "\t"+ 
						collapsedInMapped_v[i] + "\n")

			outPut.write(newLine)

		outPut.close()
					