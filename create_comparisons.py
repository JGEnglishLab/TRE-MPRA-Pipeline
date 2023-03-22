import os
import pandas as pd
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(
                    prog = 'set up comparisons',
                    description = 'Takes optional input of comparison csvs')

parser.add_argument('-p', '--pairwise_csv', dest="pairwise_comps", required=False)
parser.add_argument('-m', '--multi_csv', dest="multi_comps", required=False)
args = parser.parse_args()

pairwise_comps=args.pairwise_comps
multi_comps=args.multi_comps


treatment_list = [] #Used for user manually creating comparisons
run_treatment_dict = {} #Used for checking an input CSV


for d in os.listdir("./runs"): #Loop through directories (ie run names)
	if (d != ".DS_Store" and d != "19664"):
		if os.path.exists(f"./runs/{d}/metaData.csv"):
			md = pd.read_csv(f"./runs/{d}/metaData.csv")
			cur_treatments = md.treatment.unique()
			for t in cur_treatments:
				if t != "DNA":
					treatment_list.append(t + " || " + d)
					run_treatment_dict[d] = t
				
treatment_dict = {} #Used for user manually creating comparisons
counter = 0
for i in treatment_list:
	counter+=1
	treatment_dict[counter] = i


# *****************************************************************************
# *****************************************************************************
#                                 Pairwise comparisons
# *****************************************************************************
# *****************************************************************************


if not pairwise_comps: #They will manually enter the pairwise comparisons
	print("id\ttreatment\trun")
	def print_menu():
		for key in treatment_dict:
			t = treatment_dict[key].split(" || ")[0]
			r = treatment_dict[key].split(" || ")[1]
			print(f"{key}\t{t}\t{r}")
	print_menu()

	print("***********************************")
	print("***********************************")
	print("***********************************")
	print("Now, enter the pairwise comparisons that you want to make")
	print("Enter the numbers corresponding to the treatments you want to compare")
	print("The first treatment entered will be the base, the second with be the stimulated condition")
	print("Seperate the two numbers with a comma e.g. \"1,2\"")
	print("Enter \"done\" when finished")
	print("Enter \"menu\" to see list again")
	print("\n")

	getTreatments = True
	now = datetime.now()
	dt_string = now.strftime("%d-%m-%Y_%H-%M-%S")
	pairwise_comparison_name = f"pairwise_comparisons_{dt_string}.csv"
	pairwise_f = open(pairwise_comparison_name, "w+")
	head="id,base_treatment,stim_treatment,base_run,stim_run\n"
	pairwise_f.write(head)

	iter_id = 0
	while getTreatments:
		comp = input("Enter pairwise comparison starting with basal treatment e.g. \"1,3\" (Type \"done\" after entering all comparisons): ")
		if comp.lower() == "done":
			break
		if comp.lower() == "menu":
			print_menu()
			continue
			
		iter_id+=1

		comp = comp.split(",")
		base = int(comp[0].strip())
		stim = int(comp[1].strip())

		base_treatment = treatment_dict[base].split(" || ")[0]
		base_run = treatment_dict[base].split(" || ")[1]
		stim_treatment = treatment_dict[stim].split(" || ")[0]
		stim_run = treatment_dict[stim].split(" || ")[1]

		line = f"{iter_id},{base_treatment},{stim_treatment},{base_run},{stim_run}\n"
		pairwise_f.write(line)

	os.system(f"Rscript run_pairwise_comparisons.R ./{pairwise_comparison_name}")

else: #Check their input file
	if not os.path.exists(pairwise_comps):
		print(f"file {pairwise_comps} not found!")
		exit()
	f = open(pairwise_comps, "r+")
	id_list = []
	line_number = 0
	for line in f:
		line = line.strip()
		line_number+=1
		if line_number == 1:
			if line != "id,base_treatment,stim_treatment,base_run,stim_run":
				print(f"pairwise input file doesn't have correct header")
		if len(line.split(",")) == 5:
			print(f"{line} from pairwise input doesn't follow the correct format")
		cur_id, cur_base_treatment, cur_stim_treatment, cur_base_run, cur_stim_run = line.split(",")

		if not cur_id.isnumeric():
			print("error in pairwise tsv on line {line} \nid column must be numeric")
			exit()
		if cur_id in id_list:
			print("error in pairwise tsv on line {line} \nall ids must be unique")
			exit()
		if run_treatment_dict[cur_base_run] != cur_base_treatment:
			print("error in pairwise tsv on line {line} \nbasal run name and treatment name do not match")
			exit()
		if run_treatment_dict[cur_stim_run] != cur_stim_treatment:
			print("error in pairwise tsv on line {line} \nbasal run name and treatment name do not match")
			exit()
	os.system(f"Rscript run_pairwise_comparisons.R {pairwise_comps}")



# *****************************************************************************
# *****************************************************************************
#                                 Multi comparisons
# *****************************************************************************
# *****************************************************************************













