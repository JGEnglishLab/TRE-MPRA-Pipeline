library(tidyverse)
library(dplyr)
PATH_TO_INPUT_CSVS = "./raw_counts/"
PATH_TO_STARCODE = "./star_code/"

wd = getwd()

# Get the command-line arguments.
args = commandArgs(trailingOnly=TRUE)

#Read in meta data (From snakemake)
mData=read_tsv(args[1])

#Read in barcode map (From snakemake)
barcodeMap=read_csv(args[2])

#Set up location where files will be written for starcode
STARCODE_DIR = paste0(wd,"/star_code")

# Loop through all of the sample numbers
# For every sample number read every file of that sample number
# Group by barcode, and summarise by the sum of the counts 
# We do this to add all the counts across lanes. 
# Add a boolean that states if each barcode is in the dictionary or not
# Write file to starcode directory

for (sampleNumber in unique(mData$sampleNumber)){
  
  rows = mData %>%
    filter(sampleNumber == !!sampleNumber)
  files = rows$fileName
  cur_sample_number_file = data.frame(matrix(ncol = 2, nrow = 0))
  
  for (file in files){ #Read every file of the same sample number
    cur_file = read_csv(paste0(PATH_TO_INPUT_CSVS,file), col_names = F)
    cur_sample_number_file = rbind(cur_file, cur_sample_number_file)
  }
  
  output_name = paste0(PATH_TO_STARCODE, "sample", sampleNumber,"_mapped.tsv")
  
  cur_sample_number_file %>%
    group_by(X1)%>% #Group by barcode
    summarise(X2 = sum(X2)) %>% #Summarise the counts as the sum of the counts
    mutate(mapped = X1 %in% barcodeMap$barcode) %>% #Get a boolean that says if the barcode is in map or not
    write_tsv(output_name, col_names = F)
}

