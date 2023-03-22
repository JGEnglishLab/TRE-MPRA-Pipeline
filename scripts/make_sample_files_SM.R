library(tidyverse)
library(dplyr)
library(stringr)
PATH_TO_BARCODE_MAP = "./barcode_map_data/finalBarcodeMap.csv"
PATH_TO_INPUT_CSVS = "./raw_counts/"
PATH_TO_STARCODE = "./star_code/"
PATH_TO_STATS = "./run_descriptive_stats/"
PATH_TO_RNA_DNA = "./rna_dna_samples/"

wd = getwd()

# Get the command-line arguments.
args = commandArgs(trailingOnly=TRUE)

#Read in meta data
mData=read_csv(args[1])
print("INPUT CSV")
print(mData)

#read in barcode map
barcodeMap=read_csv(args[2])

print("NUM ARGS")
print(length(args))

addType <- function(collapsedAnalyzed, inpName){
  collapsedAnalyzed$name = inpName
  collapsedAnalyzed$type = "NONE"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == TRUE & collapsedAnalyzed$collapsedInMapped == "NULL"] = "type1"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == FALSE & collapsedAnalyzed$collapsedInMapped == "NULL"] = "type2"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == TRUE & collapsedAnalyzed$collapsedInMapped != "NULL"] = "type3" 
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == FALSE & collapsedAnalyzed$collapsedInMapped != "NULL"] = "type4"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == TRUE & collapsedAnalyzed$collapsedInMapped != "NULL" & collapsedAnalyzed$architecturesMatch == T] = "type1.2"
  
  #If it is the spike in, keep as type1
  collapsedAnalyzed$type[collapsedAnalyzed$centroid == "TAAATATGCCTCAGCACCCTGCTG"] = "type1"
  collapsedAnalyzed$type[collapsedAnalyzed$centroid == "AAGACGCGTCACAGACTTATAGAC"] = "type1"
  collapsedAnalyzed$type[collapsedAnalyzed$centroid == "CGGAGACACTTAATAGCCTCTAAC"] = "type1"
  collapsedAnalyzed$type[collapsedAnalyzed$centroid == "ATGTTAGTGAGTGTGCGAAGTAGG"] = "type1"
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
  if (i != 1 && i != 2){ #First two arguments are metadata and barcode map
    print(args[i])
    curFile = read_tsv(args[i])
    sampleName = str_split(args[i], "_", n = Inf, simplify = FALSE)[[1]][4] #Gets the sample name ie "sample2"
    sampleNumber = substring(sampleName, 7) #Removes "sample" and just gets "2"
    print("SAMPLE NUMBER")
    print(sampleNumber)
    curFile = addType(curFile, sampleName)
    # curFile$treatment = unique(mData[mData$sampleNumber == i - 2,]$treatment)
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
# 
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
  summarise(totalCount = sum(totalCounts),
            nrow = n()) -> summaryStats

longFileFiltered %>%
  filter(centroid == "TAAATATGCCTCAGCACCCTGCTG" |
           centroid == "AAGACGCGTCACAGACTTATAGAC" |
           centroid == "CGGAGACACTTAATAGCCTCTAAC" |
           centroid == "ATGTTAGTGAGTGTGCGAAGTAGG")%>%
  mutate(spike = case_when(centroid == "TAAATATGCCTCAGCACCCTGCTG" ~ "Spike1",
                           centroid == "AAGACGCGTCACAGACTTATAGAC" ~ "Spike2",
                           centroid == "CGGAGACACTTAATAGCCTCTAAC" ~ "Spike3",
                           centroid == "ATGTTAGTGAGTGTGCGAAGTAGG" ~ "Spike4")) %>%
  select(totalCounts, name, spike, treatment) -> spikeTotals

#Join all data for a summary of the run.
left_join(summaryStats, spikeTotals) %>%
  pivot_wider(names_from = spike, values_from = totalCounts) %>%
  mutate(Spike1_percent = Spike1 / totalCount,
         Spike2_percent = Spike2 / totalCount,
         Spike3_percent = Spike3 / totalCount,
         Spike4_percent = Spike4 / totalCount) -> completeSummaryStats


completeSummaryStats[is.na(completeSummaryStats)] = 0
completeSummaryStats$treatment <- factor(completeSummaryStats$treatment, levels=unique(mData$treatment))

completeSummaryStats %>%
  select(name,Spike1_percent, Spike2_percent, Spike3_percent, Spike4_percent) %>%
  pivot_longer(!name, names_to = "spikeIn") %>%
  left_join(completeSummaryStats) %>%
  ggplot(aes(x = name, y = value, fill = spikeIn)) +
  geom_bar(stat = "Identity", position="dodge") +
  labs(x="Sample Name", y = "Spike Count / Total Count") +
  facet_grid(cols = vars(treatment), scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.x = element_text(size = 5)) -> spike_in_stats

png(paste0(PATH_TO_STATS,"spike_in_stats.png"))
print(spike_in_stats)
dev.off()

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
  filter(barcode != "TAAATATGCCTCAGCACCCTGCTG", #Get rid of all the spike ins
         barcode != "AAGACGCGTCACAGACTTATAGAC",
         barcode != "CGGAGACACTTAATAGCCTCTAAC",
         barcode != "ATGTTAGTGAGTGTGCGAAGTAGG")-> allDataFiltered



allDataFiltered %>%
  left_join(barcodeMap, by = "barcode") -> allDataFilteredJoined


allDataFilteredJoined %>%
  mutate(architecture = paste0(motif,":", id,", ", period,", ", spacer,", ", promoter)) -> allDataFilteredJoined

write_csv(allDataFilteredJoined, paste0(PATH_TO_RNA_DNA,"all_data_filtered.csv"))
write_csv(dnaSamples, paste0(PATH_TO_RNA_DNA,"dna_samples.csv"))
write_csv(rnaSamples, paste0(PATH_TO_RNA_DNA,"rna_samples.csv"))
#
# allDataFilteredJoined %>%
#   pivot_wider(id_cols=c(barcode,architecture),names_from = treatment, values_from = c(rpm,totalCounts),  names_glue = "_{treatment}_{.value}") -> wide_all
#
# wide_all[is.na(wide_all)] = 0
#
# write_csv(wide_all,  paste0(RUN_NAME, "/all_data_filtered_joined_wide_no_sc_", RUN_NAME, ".csv"))
#
# dna_data %>%
#   left_join(barcodeMap, by = "barcode") %>%
#   mutate(architecture = paste0(motif,":", id,", ", period,", ", spacer,", ", promoter)) %>%
#   select(barcode, architecture, DNA_count, DNA_rpm)-> dna_csv
#
# write_csv(dna_csv,  paste0(RUN_NAME, "/dna_data_no_sc_", RUN_NAME, ".csv"))


