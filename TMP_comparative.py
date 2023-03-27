import os
import pandas as pd
import argparse
from datetime import datetime


def check_y_n_inp(inp, correction_message = "Try again"):
	inp = inp.lower()
	while inp != "y" and inp != "n":
		print(f"{inp} is an invalid response")
		print(correction_message)
		inp = input("Enter \"y\" or \"n\": ")
	return inp

DEFAULT_THREADS = 6

parser = argparse.ArgumentParser(
                    prog = 'set up comparisons',
                    description = 'Takes optional input of comparison Tsvs')

parser.add_argument('-p', '--pairwise_tsv', dest="pairwise_comps", required=False)
parser.add_argument('-m', '--multi_tsv', dest="multi_comps", required=False)
parser.add_argument('-n', '--n_threads', dest="threads", required=False)

args = parser.parse_args()

now = datetime.now()
dt_string = now.strftime("%d-%m-%Y_%H-%M")

pairwise_comps=args.pairwise_comps
multi_comps=args.multi_comps
threads = args.threads

if not threads:
	threads = DEFAULT_THREADS
elif not threads.isnumeric():
	print(f"Error in -n flag. Must be numeric. Using default of {DEFAULT_THREADS}")
	threads = DEFAULT_THREADS


treatment_list = [] #Used for user manually creating comparisons
run_treatment_dict = {} #Used for checking an input TSV


for d in os.listdir("./runs"): #Loop through directories (ie run names)
	if (d != ".DS_Store"):
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


def print_menu():
	for key in treatment_dict:
		t = treatment_dict[key].split(" || ")[0]
		r = treatment_dict[key].split(" || ")[1]
		print(f"{key}\t{t}\t{r}")


# *****************************************************************************
# *****************************************************************************
#                                 Pairwise comparisons
# *****************************************************************************
# *****************************************************************************
run_pairwise = False


if not pairwise_comps: #They will manually enter the pairwise comparisons
	manual_pairwise = input("Would you like to manually input pairwise comparisons? (y/n): ")
	check_y_n_inp(manual_pairwise)
	if manual_pairwise == "y":
		print("id\ttreatment\trun")
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
		pairwise_comparison_name = f"entered_pairwise_comparisons/pairwise_comparisons_{dt_string}.tsv"
		pairwise_f = open(pairwise_comparison_name, "w+")
		head="id\tbase_treatment\tstim_treatment\tbase_run\tstim_run\n"
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
			if len(comp) != 2:
				print("For pairwise comparisons please enter exactly 2 ids")
				continue

			base = int(comp[0].strip())
			stim = int(comp[1].strip())

			base_treatment = treatment_dict[base].split(" || ")[0]
			base_run = treatment_dict[base].split(" || ")[1]
			stim_treatment = treatment_dict[stim].split(" || ")[0]
			stim_run = treatment_dict[stim].split(" || ")[1]

			line = f"{iter_id}\t{base_treatment}\t{stim_treatment}\t{base_run}\t{stim_run}\n"
			pairwise_f.write(line)
		pairwise_f.close()
		run_pairwise = True


else: #Check their pairwise input file
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
			if line != "id\tbase_treatment\tstim_treatment\tbase_run\tstim_run":
				print(f"pairwise input file, {pairwise_comps} doesn't have correct header")
			continue
		if len(line.split("\t")) != 5:
			print(f"\"{line}\" from pairwise input doesn't follow the correct format")
			exit()
		cur_id, cur_base_treatment, cur_stim_treatment, cur_base_run, cur_stim_run = line.split("\t")

		if not cur_id.isnumeric():
			print(f"error in pairwise tsv on line \"{line}\" \nid column must be numeric")
			exit()
		if cur_id in id_list:
			print(f"error in pairwise tsv on line \"{line}\" \nall ids must be unique")
			exit()

		id_list.append(cur_id)

		try:
			if run_treatment_dict[cur_base_run] != cur_base_treatment:
				print(f"error in pairwise tsv on line \"{line}\" \n{cur_base_treatment} doesn't belong to the {cur_base_run} run.")
				exit()
		except KeyError:
			print(f"error in pairwise tsv on line \"{line}\" \n{cur_base_run} doesn't exist")
			exit()

		try:
			if run_treatment_dict[cur_stim_run] != cur_stim_treatment:
				print(f"error in pairwise tsv on line \"{line}\" \n{cur_stim_treatment} doesn't belong to the {cur_stim_run} run.")
				exit()
		except KeyError:
			print(f"error in pairwise tsv on line \"{line}\" \n{cur_stim_run} doesn't exist")
			exit()

	run_pairwise = True



# *****************************************************************************
# *****************************************************************************
#                                 Multi comparisons
# *****************************************************************************
# *****************************************************************************

#If they didn't enter a multi TSV
run_multi = False
if not multi_comps:
	manual_pairwise = input("\nWould you like to manually input multi comparisons? (y/n): ")
	check_y_n_inp(manual_pairwise)
	if manual_pairwise == "y":
		print("id\ttreatment\trun")
		print_menu()

		print("***********************************")
		print("***********************************")
		print("***********************************")
		print("Now, enter the multi comparisons that you want to make")
		print("Enter the numbers corresponding to the treatments you want to compare")
		print("Separate the two numbers with a comma e.g. \"1,2,4,5\"")
		print("Enter \"done\" when finished")
		print("Enter \"menu\" to see list again")
		print("\n")

		getTreatments = True
		multi_comparison_name = f"entered_multi_comparisons/multi_comparisons_{dt_string}.tsv"
		multi_f = open(multi_comparison_name, "w+")
		head = "id\ttreatment\trun\n"
		multi_f.write(head)

		iter_id = 0
		while getTreatments:
			comp = input(
				"Enter multi comparisons separated by commas \"1,2,4,5\" (Type \"done\" after entering all comparisons): ")
			if comp.lower() == "done":
				break
			if comp.lower() == "menu":
				print_menu()
				continue

			iter_id += 1

			comp = comp.split(",")
			if not len(comp) >= 3:
				print("For multi comparisons please enter 3 or more")
				continue


			for i in comp:
				i = int(i)
				cur_treatment = treatment_dict[i].split(" || ")[0]
				cur_run = treatment_dict[i].split(" || ")[1]
				line = f"{iter_id}\t{cur_treatment}\t{cur_run}\n"
				multi_f.write(line)
		run_multi = True
		multi_f.close()


else: #If their multi-input file
	if not os.path.exists(multi_comps):
		print(f"file {multi_comps} not found!")
		exit()
	f = open(multi_comps, "r+")
	id_list = []
	line_number = 0
	for line in f:
		line = line.strip()
		line_number+=1
		if line_number == 1:
			if line != "id\ttreatment\trun":
				print(f"multi input file, {multi_comps} doesn't have correct header")
			continue
		if len(line.split("\t")) != 3:
			print(f"\"{line}\" from multi input file doesn't follow the correct format")
			exit()
		cur_id, cur_treament, cur_run = line.split("\t")
		if not cur_id.isnumeric():
			print(f"error in multi tsv on line \"{line}\" \nid column must be numeric")
			exit()
		try:
			if run_treatment_dict[cur_run] != cur_treament:
				print(f"error in multi tsv on line \"{line}\" \n{cur_treament} doesn't belong to the {cur_run} run.")
				exit()
		except KeyError:
			print(f"error in multi tsv on line \"{line}\" \n{cur_run} doesn't exist")
	run_multi = True


# *****************************************************************************
# *****************************************************************************
#                                 Run All
# *****************************************************************************
# *****************************************************************************

# if run_pairwise:
# 	if not pairwise_comps: #IE they didn't enter a TSV
# 		os.system(f"Rscript scripts/run_pairwise_comparisons.R ./{pairwise_comparison_name} {threads}")
# 	else: #They did enter a TSV
# 		os.system(f"Rscript scripts/run_pairwise_comparisons.R {pairwise_comps} {threads}")
#
# if run_multi:
# 	if not pairwise_comps: #IE they didn't enter a TSV
# 		os.system(f"Rscript scripts/run_multi_comparisons.R ./{multi_comparison_name} {threads}")
# 	else: #They did enter a TSV
# 		os.system(f"Rscript scripts/run_multi_comparisons.R {multi_comps} {threads}")
#




#TODO
# Fix run treatment dict
# Test flags
# In R scripts
# Spit out things that missing architectures



