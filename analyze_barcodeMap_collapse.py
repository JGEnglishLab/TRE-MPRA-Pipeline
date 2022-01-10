
scOutput = open("sc_out_filteredBarcodeMap2.tsv", "r+")
analyzedOutput = open("analyzed_filteredBarcodeMap2.tsv", "w+")

analyzedOutput.write("barcode\tnumCollapsed\tisCentroid\n")

for line in scOutput:
	line = line.split()
	centroid = line[0]
	numCollapsed = line[1]
	collapsedBarcodes = line[2].split(",")

	for bc in collapsedBarcodes:
		bc = bc.strip()
		isCentroid = (centroid.strip() == bc)
		analyzedOutput.write(bc +"\t" + numCollapsed + "\t" + str(isCentroid) + "\n")

analyzedOutput.close()
scOutput.close()