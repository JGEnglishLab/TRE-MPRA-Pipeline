library(tidyverse)
library(dplyr)
library(stringr)
PATH_TO_BARCODE_MAP = "./barcode_map_data/finalBarcodeMap.csv"
PATH_TO_INPUT_CSVS = "./raw_counts/"
PATH_TO_STARCODE = "./star_code/"
PATH_TO_STATS = "./run_descriptive_stats/"
PATH_TO_RNA_DNA = "./rna_dna_samples/"

 #l = c("star_code/analyzed_out_sample43_mapped_sc_out.tsv", "star_code/analyzed_out_sample2_mapped_sc_out.tsv", "star_code/analyzed_out_sample42_mapped_sc_out.tsv", "star_code/analyzed_out_sample8_mapped_sc_out.tsv", "star_code/analyzed_out_sample5_mapped_sc_out.tsv", "star_code/analyzed_out_sample1_mapped_sc_out.tsv", "star_code/analyzed_out_sample6_mapped_sc_out.tsv", "star_code/analyzed_out_sample3_mapped_sc_out.tsv")


check_valid_nucleotide <- function(string){
  valid_nuc = c("A", "a", "T", "t", "G", "g", "C", "c")
  for (i in strsplit(string, "")[[1]]){
    if (!(i %in% valid_nuc)){
      return(FALSE)
    }
  }
  return(TRUE)
}

wd = getwd()

# Get the command-line arguments.
args = commandArgs(trailingOnly=TRUE)

#Read in meta data
mData=read_csv(args[1])

#read in barcode map
barcodeMap=read_csv(args[2])

#read in spike-in file
spike_ins = c()
if (file.exists(args[3])){
  spikeInFile=read_table(args[3], col_names = F)
  for (i in spikeInFile$X1){
    if (nchar(i) == 24 && check_valid_nucleotide(i)) {
      
      spike_ins = c(spike_ins, i)
    }
  }
}

addType <- function(collapsedAnalyzed, inpName){
  collapsedAnalyzed$name = inpName
  collapsedAnalyzed$type = "NONE"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == TRUE & collapsedAnalyzed$collapsedInMapped == "NULL"] = "type1"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == FALSE & collapsedAnalyzed$collapsedInMapped == "NULL"] = "type2"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == TRUE & collapsedAnalyzed$collapsedInMapped != "NULL"] = "type3" 
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == FALSE & collapsedAnalyzed$collapsedInMapped != "NULL"] = "type4"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == TRUE & collapsedAnalyzed$collapsedInMapped != "NULL" & collapsedAnalyzed$architecturesMatch == T] = "type1.2"
  
  #If it is the spike in, keep as type1
  collapsedAnalyzed$type[collapsedAnalyzed$centroid %in% spike_ins] = "type1"
  
  # collapsedAnalyzed$type[collapsedAnalyzed$centroid == "TAAATATGCCTCAGCACCCTGCTG"] = "type1"
  # collapsedAnalyzed$type[collapsedAnalyzed$centroid == "AAGACGCGTCACAGACTTATAGAC"] = "type1"
  # collapsedAnalyzed$type[collapsedAnalyzed$centroid == "CGGAGACACTTAATAGCCTCTAAC"] = "type1"
  # collapsedAnalyzed$type[collapsedAnalyzed$centroid == "ATGTTAGTGAGTGTGCGAAGTAGG"] = "type1"
  collapsedAnalyzed %>%
    mutate(
      collapse_stat = 0,
      collapse_stat = if_else(totalCollapsed!=1, centroidCounts/totalCollapsed, collapse_stat)) -> collapsedAnalyzed
  return(collapsedAnalyzed)
}

#Make a new long file for all the filtered/analyzed samples
longFile <- data.frame(matrix(ncol = 14, nrow = 0))

#Read in the analyzed starcode data
#Give each cluster a type from addType function
#Type 1 = centroid is in map & no clustered barcodes in map
#Type 1.2 = centroid is in map & all clustered barcodes that are in the map belong to same architecture as the centroid
#Type 2 = centroid not in map & none of the clustered barcodes are in the map
#Type 3 = centroid is in map & one or more clustered barcode is in map
#Type 4 = centroid not in map & one or more clustered barcode is in map
#We only keep type 1 and 1.2

for (i in 1:length(args)){
  if (i != 1 && i != 2 && i != 3){ #First three arguments are metadata, barcode map, spike in file
    print(args[i])
    curFile = read_tsv(args[i])
    sampleName = str_split(args[i], "_", n = Inf, simplify = FALSE)[[1]][4] #Gets the sample name ie "sample2"
    sampleNumber = substring(sampleName, 7) #Removes "sample" and just gets "2"
    print("SAMPLE NUMBER")
    print(sampleNumber)
    curFile = addType(curFile, sampleName)
    curFile$treatment = unique(mData[mData$sampleNumber == sampleNumber,]$treatment)
    longFile = rbind(longFile,curFile)
  }
}
unique(mData[mData$sampleNumber == i, ]$treatment)
#Calculate the RPM.
longFile %>%
  group_by(name) %>%
  summarise(norm = sum(totalCounts) / 1000000) %>%
  right_join(longFile) %>%
  mutate(rpm = totalCounts/norm) -> longFile


#Get the ratios of the filtered reads, to reads that made it through
longFile %>%
  mutate(filtered = ifelse(type == "type1" | type == "type1.2", "Not Filtered", "Filtered")) %>%
  group_by(name, filtered) %>%
  summarise(sum = sum(totalCounts)) %>%
  ggplot(aes(x = name, y = sum, fill = filtered)) +
  geom_bar(stat = "identity") + coord_flip() + labs( y = "Total Counts") -> filtered_ratios

png(paste0(PATH_TO_STATS,"filtering_ratios.png"))
print(filtered_ratios)
dev.off()

write_csv(longFile, paste0(PATH_TO_STATS,"pre_filtering_barcodes.csv"))

# #DELETE LATER
# setwd("/Volumes/external_disk/english_lab/mpra_snake_make/runs/full-test")
# longFile = read_csv(paste0(PATH_TO_STATS,"pre_filtering_barcodes.csv"))
# mData = read_csv("metaData.csv")

# ####

#Get all types
longFile %>%
  group_by(name, type) %>%
  summarise(sum = sum(totalCounts)) %>%
  ggplot(aes(x = name, y = sum, fill = type)) +
  geom_bar(stat = "identity") + coord_flip() + labs(y = "Total Counts") +
  labs(title = "The number of types of collapses.\n(2,3 and 4 are filtered)") -> type_ratios

png(paste0(PATH_TO_STATS,"type_ratios.png"))
print(type_ratios)
dev.off()

longFile %>% filter(type == "type1" | type == "type1.2") -> longFileFiltered

longFileFiltered %>%
  group_by(name, treatment)%>%
  summarise(totalReads = sum(totalCounts),
            numberBarcodes = n()) -> summaryStats

longFileFiltered %>%
  filter(centroid %in% spike_ins)%>%
  mutate(spike = centroid) %>%
  select(totalCounts, name, spike, treatment) %>%
  rename(spikeInCount = totalCounts) -> spikeTotals

#Join all data for a summary of the run.
left_join(summaryStats, spikeTotals) %>%
  mutate(spike_percent = spikeInCount/totalReads) -> completeSummaryStats

if (length(spike_ins) > 0){
  completeSummaryStats %>%
    mutate(i = paste0(name,", ", treatment)) %>%
    ggplot(aes(x = i, y = spike_percent, fill = spike)) +
    geom_bar(stat = "Identity", position="dodge") +
    labs(x="Sample", y = "Spike Reads / Total Reads") +
    coord_flip() -> spike_in_stats
  
  png(paste0(PATH_TO_STATS,"spike_in_stats.png"))
  print(spike_in_stats)
  dev.off()
} else{ #Create an empty png if there are no spikes
  png::writePNG(array(0, dim = c(1,1,4)), paste0(PATH_TO_STATS,"spike_in_stats.png"))
  
}

write_csv(completeSummaryStats, paste0(PATH_TO_STATS,"spike_in_stats.csv"))


# ******************************************************************************
#                             MAKE DNA AND RNA SAMPLES
# ******************************************************************************

longFileFiltered %>% filter(treatment == "DNA") %>%
  rename(barcode = centroid) -> dnaSamples

longFileFiltered %>% filter(treatment != "DNA") %>%
  rename(barcode = centroid) -> rnaSamples

longFileFiltered %>%
  rename(barcode = centroid) -> allData

#DNA Counts per barcode
dnaSamples %>%
  ggplot(aes(x = totalCounts)) +
  geom_histogram(binwidth = 2)+
  labs(title = "DNA counts per barcode", x = "DNA count", y = "number of barcodes") -> dna_per_barcode

png(paste0(PATH_TO_STATS,"dna_per_barcode.png"))
print(dna_per_barcode)
dev.off()

fileConn<-file(paste0(PATH_TO_STATS,"dna_per_barcode.txt"))
writeLines(c("The mean number of DNA counts per barcode is ",as.character(dnaSamples$totalCounts %>% mean())), fileConn)
writeLines(c("The mean number of DNA counts per barcode is ",
             as.character(dnaSamples$totalCounts %>% mean()),
             "The median number of DNA counts per barcode is ",
             as.character(dnaSamples$totalCounts %>% median())), fileConn)
close(fileConn)



dnaSamples %>%
  select(barcode, totalCounts, rpm) %>%
  rename(DNA_count = totalCounts, DNA_rpm = rpm) -> dna_data

left_join(rnaSamples, dna_data) %>%
  mutate(rna_dna_rpm_ratio = rpm/DNA_rpm) %>%
  select(-type,
         -collapsedInMapped,
         -collapsedBarcodes,
         -medians,
         -means,
         -norm,
         -centroidCounts,
         -collapse_stat,
         -totalCollapsed,
         -centroidInMapped) %>%
  filter(!(barcode %in% spike_ins)) -> allDataFiltered #Get rid of all the spike ins

allDataFiltered %>%
  left_join(barcodeMap, by = "barcode") -> allDataFilteredJoined


allDataFilteredJoined %>%
  mutate(architecture = paste0(motif,":", id,", ", period,", ", spacer,", ", promoter)) -> allDataFilteredJoined

write_csv(allDataFilteredJoined, paste0(PATH_TO_RNA_DNA,"all_data_filtered.csv"))
write_csv(dnaSamples, paste0(PATH_TO_RNA_DNA,"dna_samples.csv"))
write_csv(rnaSamples, paste0(PATH_TO_RNA_DNA,"rna_samples.csv"))

