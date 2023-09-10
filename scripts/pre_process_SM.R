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
print(STARCODE_DIR)

#Adds a sample lane column to meta data
#Samples of the same number but different file name (different lanes) will be added
mData %>% 
  select(sampleNumber, fileName) %>% 
  unique() %>%
  group_by(sampleNumber) %>% 
  reframe(lane = row_number(), fileName = fileName) %>%
  right_join(mData) -> mData

 
#Created new merged column
mData$merged = paste0(mData$sampleNumber,"_",mData$lane)

#Create empty file to hold all the raw data
longFile <- data.frame(matrix(ncol = 3, nrow = 0))

#Read in all the raw data
for (file in mData$fileName){
  row = mData[mData$fileName == file,]
  n = paste0(row$sampleNumber,"_",row$lane)
  file_loc = paste0(PATH_TO_INPUT_CSVS, file)
  curFile = read_csv(file_loc, col_names = F) %>% mutate(name=n)
  longFile = rbind(longFile,curFile)
}

#Pivot the data 
longFile %>% pivot_wider(names_from = name, values_from = X2) %>%
  rename(barcode = X1) -> wideFile
wideFile[is.na(wideFile)] <-0

#Writes a tsv for each sample
#Now, each sample will be written in the PATH_TO_STARCODE dir with info a boolean column specifying if the barcode is in the map or not
tables_written = 0
for (i in unique(mData$sampleNumber)){
  tables_written = tables_written + 1
  sampleName = paste0("sample",i)
  curFile = wideFile %>% select("barcode", ((mData %>% filter(sampleNumber == i))$merged))
  curFile[, sampleName]= rowSums(curFile[ ,-1]) #Sum up lanes

  curFile %>%
    select(barcode,sampleName) %>%
    filter(!!sym(sampleName) !=0) %>%
    mutate(mapped = barcode %in% barcodeMap$barcode) -> curFile

  #Write all the files in prep for running starcode
  loc = paste0(PATH_TO_STARCODE, sampleName,"_mapped.tsv")
  write.table(curFile, loc, sep = "\t", col.names = FALSE, quote = F, row.names = F)
}


write_tsv(mData, "metaData.tsv")
