import os
import statistics
import sys
import pandas as pd
#Read in barcode map
bc_map_location = "../../barcode_map_data/finalBarcodeMap.csv"
bc_map = pd.read_csv(bc_map_location)

#Creat a dictionary that has the barcode as a key, and the architecture as the value
#This will be used to see if barcodes map to the same architecture
bc_dict =  dict(zip(bc_map['barcode'], bc_map['id'].astype(str)+bc_map['motif']+bc_map['spacer'].astype(str)+bc_map['period'].astype(str)+bc_map['promoter']))

class ClusterMetricsCalculator:
	def __init__(self):
		self.cluster_mean = []
		self.cluster_median = []
		self.cluster_centroid = []
		self.cluster_size = []
		self.cluster_count = []
		self.cluster_centroid_count = []
		self.cluster_centroid_in_map = []
		self.cluster_non_centroid_barcodes = []  # All barcodes that aren't the centroid
		self.cluster_non_centroid_barcodes_in_map = []  # All barcodes that aren't the centroid and are in the map
		self.cluster_non_centroid_match_centroid_architecture = []  # List of boolean values indicating whether the centroid barcode
																	# and its clustered non-centroid barcodes share the same mapped architecture.
																	# True indicates that all non-centroid barcodes that map to an architecture
																	# have the same architecture as the centroid, False otherwise.

	def calculate_metrics(self, starCodeOutputName, curCountDict, curMapDict, bc_dict):
		sc_out_file = pd.read_csv(starCodeOutputName, delimiter="\t", header=None)

		def row_metrics(row):
			cur_centroid = row[0]
			cur_cluster_count = row[1]
			cur_barcodes = row[2].split(",")

			cur_cluster_size = len(cur_barcodes)
			cur_centroid_count = curCountDict[cur_centroid]
			cur_centroid_in_map = curMapDict[cur_centroid]

			cluster_counts = []
			non_centroid_barcodes = []
			non_centroid_barcodes_in_map = []

			for bc in cur_barcodes:
				cluster_counts.append(int(curCountDict[bc]))
				if bc != cur_centroid:
					non_centroid_barcodes.append(bc)
					if curMapDict[bc]:
						non_centroid_barcodes_in_map.append(bc)

			if len(non_centroid_barcodes_in_map) == 0:
				self.cluster_non_centroid_match_centroid_architecture.append("FALSE")
				non_centroid_barcodes_in_map.append("NULL")

			elif cur_centroid_in_map:
				centroid_architecture = bc_dict[cur_centroid]
				all_same_architecture = True

				for bc in non_centroid_barcodes_in_map:
					if bc_dict[bc] != centroid_architecture:
						all_same_architecture = False
				if all_same_architecture:
					self.cluster_non_centroid_match_centroid_architecture.append("TRUE")
				else:
					self.cluster_non_centroid_match_centroid_architecture.append("FALSE")

			else:
				self.cluster_non_centroid_match_centroid_architecture.append("FALSE")

			cur_mean = statistics.mean(cluster_counts)
			cur_median = statistics.median(cluster_counts)

			self.cluster_mean.append(cur_mean)
			self.cluster_median.append(cur_median)
			self.cluster_centroid.append(cur_centroid)
			self.cluster_size.append(cur_cluster_size)
			self.cluster_count.append(cur_cluster_count)
			self.cluster_centroid_count.append(cur_centroid_count)
			self.cluster_centroid_in_map.append(cur_centroid_in_map)
			self.cluster_non_centroid_barcodes.append(",".join(non_centroid_barcodes))
			self.cluster_non_centroid_barcodes_in_map.append(",".join(non_centroid_barcodes_in_map))

		sc_out_file.apply(row_metrics, axis=1)
	def printTSV(self,outputName):


		# dictionary of lists
		metrics_dict = {'mean_cluster_count': self.cluster_mean,
				'median_cluster_count': self.cluster_median,
				'total_clustered': self.cluster_size,
				'total_counts': self.cluster_count,
				'centroid': self.cluster_centroid,
				'centroid_counts': self.cluster_centroid_count,
				'centroid_in_map': self.cluster_centroid_in_map,
				'clustered_barcodes': self.cluster_non_centroid_barcodes,
				'architectures_match': self.cluster_non_centroid_match_centroid_architecture,
				'clustered_in_map': self.cluster_non_centroid_barcodes_in_map}

		df = pd.DataFrame(metrics_dict)
		df.to_csv(outputName, sep="\t")

for starCodeOutputName in sys.argv:
	if starCodeOutputName.endswith(".py"): #First arg will be the name of the script
		continue

	#*************************************
	#	Get all file names configured
	#*************************************

	#Starcode name will be in this format
	#star_code/sample4_map_sc_out.tsv
	print(starCodeOutputName)

	#Need output from pre_process_SM.R (pre starcode)
	#Will be in this format
	#star_code/sample4_map.tsv
	preStarcodeName = starCodeOutputName.replace("_sc_out.tsv", ".tsv")
	print(preStarcodeName)

	#Get the output name
	# star_code/analyzed_out_sample4_map_sc_out.tsv
	outputName = starCodeOutputName.replace("sample", "analyzed_out_sample")
	print(outputName)

	#*********************************************
	#	Make dictionaries on pre-starcode file
	#*********************************************
	colnames = ["barcode","count","map"]
	pre_sc_file = pd.read_csv(preStarcodeName, sep = "\t", names=colnames)

	curCountDict = dict(zip(pre_sc_file['barcode'], pre_sc_file['count']))
	#'map' = a boolean stating if the barcode is in the barcode map or not
	curMapDict = dict(zip(pre_sc_file['barcode'], pre_sc_file['map']))


	cur_metrics = ClusterMetricsCalculator()
	cur_metrics.calculate_metrics(starCodeOutputName, curCountDict, curMapDict, bc_dict)
	cur_metrics.printTSV(outputName)