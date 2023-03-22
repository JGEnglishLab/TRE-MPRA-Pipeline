import os

gatherFiles = True
getFileName = True

files = []
fileSamplePair = {}

for filename in os.listdir("./raw_counts/"):
	files.append(filename)

#Gets the name of the 
run_name = os.getcwd().split("/")[-1]

while gatherFiles:
	print("Files gathered")
	for f in files:
		print(f)
	print("Total Files = " + str(len(files)))
	complete = input("Are these all the files that you want to work with? (y/n): ")

	if complete == "y":
		gatherFiles = False
		break
	else:
		manualNames = True
		print("Enter file names exactly as they are written, type \"done\" when finished.")
		while manualNames:
			manualNames = input("Enter File Name: ")
			if manualNames.lower() == "done":
				manualNames = False
				break
			else:
				files.append(manualNames)

for f in files:

	splitFile = f.split("_")
	if len(splitFile) == 9:
		sample = splitFile[5].replace('S', '') 
		fileSamplePair[f] = sample

	else:
		print("\n****************************************")
		print("Different file format detected for \"" + f + "\"")
		curSampleName = input("Enter sample name e.g. \"14\" for sample 14: ")
		#FIXME make sure that "S**"
		fileSamplePair[f] = curSampleName


print("\n")
print("Now, enter sample names, and the corresponding treatment name.")
print("The names you enter will be exactly the same in all the images generated")
print("Enter the treatment type first e.g. \"Serum Free\"")
print("Then enter the corresponding sample numbers sperated by commas e.g. \"1,2,3,4\"")
print("Enter \"done\" when finished")
print("\n")


getTreatments = True
treatments = {}

while getTreatments:
	treatment = input("Enter treatment name e.g. \"GPR41 No Drug\" (Type \"done\" when finished): ")
	if treatment.lower() == "done":
		break
	samples = input("Enter corresponding sample numbers seperated by commas e.g. \"19,20,21\": ")
	for i in samples.split(","):
		treatments[i] = treatment

print(treatments)

#FIXME Create a way for them to fix it if they need


f = open("metaData.csv", "w+")
f.write("fileName,sampleNumber,treatment,run_name\n")

#Make big CSV
for file in fileSamplePair:
	line = file + "," + fileSamplePair[file] + "," + treatments[fileSamplePair[file]] + "," + run_name + "\n"
	f.write(line)








