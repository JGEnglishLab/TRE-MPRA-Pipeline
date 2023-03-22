library(tidyverse)
library(dplyr)
PATH_TO_INPUT_CSVS = "./raw_counts/"
PATH_TO_STARCODE = "./star_code/"

wd = getwd()

# Get the command-line arguments.
args = commandArgs(trailingOnly=TRUE)

#Read in meta data
mData=read_csv(args[1])
print("INPUT CSV")
print(mData)

#Read in barcode map
barcodeMap=read_csv(args[2])

#Set up location where files will be written for starcode
STARCODE_DIR = paste0(wd,"/star_code")
print(STARCODE_DIR)

#Adds a sample replicate column to meta data
#The first time a treatment type is in mData it will be labled as replicate 1
#Each subsequent time that treatment is seen it will be labled as the next replecate
mData$sampleReplicate = 999
sampleNumbers = as.array(mData$sampleNumber)
for (file in mData$fileName){
  curSampleNumber = (mData %>% filter(fileName == file))$sampleNumber
  curReplicateSum = sum(sampleNumbers == curSampleNumber)

  for (i in 1:length(sampleNumbers)){
    if (sampleNumbers[i] == curSampleNumber){
      sampleNumbers = sampleNumbers[-i]
      break
    }
  }
  mData[mData$fileName == file,]$sampleReplicate = curReplicateSum
}
 
#Created new merged column
mData$merged = paste0(mData$sampleNumber,"_",mData$sampleReplicate)

#Read in all the raw data
longFile <- data.frame(matrix(ncol = 3, nrow = 0))

for (file in mData$fileName){
  row = mData[mData$fileName == file,]
  n = paste0(row$sampleNumber,"_",row$sampleReplicate)
  file_loc = paste0(PATH_TO_INPUT_CSVS, file)
  curFile = read_csv(file_loc, col_names = F) %>% mutate(name=n)
  longFile = rbind(longFile,curFile)
}

longFile %>% pivot_wider(names_from = name, values_from = X2) %>%
  rename(barcode = X1) -> wideFile
wideFile[is.na(wideFile)] <-0

#Writes a tsv for each 
tables_written = 0
for (i in unique(mData$sampleNumber)){
  tables_written = tables_written + 1
  sampleName = paste0("sample",i)
  curFile = wideFile %>% select("barcode", ((mData %>% filter(sampleNumber == i))$merged))
  curFile[, sampleName]= rowSums(curFile[ ,-1])

  curFile %>%
    select(barcode,sampleName) %>%
    filter(!!sym(sampleName) !=0) %>%
    mutate(mapped = barcode %in% barcodeMap$barcode) -> curFile

  #Write all the files in prep for running starcode
  loc = paste0(PATH_TO_STARCODE, sampleName,"_mapped.tsv")
  write.table(curFile, loc, sep = "\t", col.names = FALSE, quote = F, row.names = F)
}

print(tables_written)
print(mData)
write_csv(mData, "metaData.csv")
