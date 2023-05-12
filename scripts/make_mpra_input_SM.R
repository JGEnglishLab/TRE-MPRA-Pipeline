library(tidyverse)
library(dplyr)
library(testit)
PATH_TO_MPRA_INPUT = "./mpra_input/"

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

#read in input
barcodeMap=read_csv(args[1])
allDataFilteredJoined=read_csv(args[2])
rnaSamples=read_csv(args[3])
dnaSamples=read_csv(args[4])
mData=read_csv(args[5])

#read in spike-in file
spike_ins = c()
if (file.exists(args[7])){
  spikeInFile=read_table(args[7], col_names = F)
  for (i in spikeInFile$X1){
    if (nchar(i) == 24 && check_valid_nucleotide(i)) {

      spike_ins = c(spike_ins, i)
    }
  }
}

RUN_NAME = unique(mData$run_name)

# #####
# #DELETE LATER!
#setwd("/Volumes/external_disk/english_lab/TRE-MPRA-Pipeline/runs/run1t")
# 
#barcodeMap = read_csv("../../barcode_map_data/finalBarcodeMap.csv")
#allDataFilteredJoined = read_csv("rna_dna_samples/all_data_filtered.csv")
#rnaSamples = read_csv("rna_dna_samples/rna_samples.csv")
#dnaSamples = read_csv("rna_dna_samples/dna_samples.csv")
#mData = read_csv("./metaData.csv")
#####

# #Get DNA depth factor
# dnaSamples %>% 
#   group_by(name) %>%
#   summarise(uq = quantile(totalCounts, .75)) -> df
# df$run_name = RUN_NAME
# write_csv(df, paste0(PATH_TO_MPRA_INPUT, "dna_depth.csv"))

#Get a unique id for each treatment
run_name = unique(mData$run_name)
allDataFilteredJoined %>%
  group_by(treatment) %>%
  summarise(n = n()) %>%
  select(-n) %>%
  mutate(run = RUN_NAME) %>%
  mutate(t_id = paste0(run_name,'_t_',row_number())) -> treatment_ids

#Filter 
dnaSamples %>%
  select(totalCounts, barcode, name)  -> current_dna_samples

current_dna_samples %>%
  left_join(barcodeMap, by = "barcode", "name") %>%
  mutate(architecture = paste0(motif,":", id,", ", period,", ", spacer,", ", promoter))-> dna_joined 

#Get the top 100 (by DNA count) barcodes for each architecture
dna_joined  %>% 
  group_by(architecture, name) %>% 
  slice_max(order_by = totalCounts, n = 100,  with_ties = FALSE) -> dna_top

dna_top %>% group_by(architecture, name) %>%
  summarise(n = n()) -> expected_n

cbind(dna_top %>%
        group_by(architecture, name) %>%
        summarise(bc = row_number()), dna_top) %>%
  mutate(architecture = architecture...1, name = name...6)%>%
  select(-architecture...14,-architecture...1,-occurrence, -name...6,-name...2)-> dna_top

dna_top %>%
  select(bc, barcode, totalCounts, architecture, name) %>%
  rename(dna_count = totalCounts) -> dna_top


#Join with RNA
rnaSamples %>% 
  filter(!(barcode %in% spike_ins)) %>%
  select(barcode, totalCounts, name, treatment) %>%
  left_join(barcodeMap, by = "barcode") %>%
  mutate(architecture = paste0(motif,":", id,", ", period,", ", spacer,", ", promoter)) %>%
  select(architecture, totalCounts, barcode, name, class, treatment) -> rna_joined 

#This is where we check to see if > 1 DNA sample exists. If it does, we must join the DNA samples to the corresponding RNA samples
if (args[6] != "None"){ #If a sixth argument is passed in that is not "None" we have > 1 dna sample
  dna_sample_map = read_tsv(args[6], col_names = F)
  ###
  #DELETE LATER
  #dna_sample_map = read_tsv("/Volumes/external_disk/english_lab/TRE-MPRA-Pipeline/test_data/input/run1_dna_map.tsv", col_names = F)
  ####
  dna_sample_map$X1 = paste0("sample",dna_sample_map$X1) # X1 is DNA name
  dna_sample_map$X2 = paste0("sample",dna_sample_map$X2) # X2 is RNA name
  
  rna_joined = left_join(rna_joined, dna_sample_map,by = c("name" = "X2")) %>% rename("DNA_name" = X1)
  dna_top %>% rename("DNA_name" = name) -> dna_top
  mpra_base = left_join(dna_top, rna_joined)
  
} else{
  only_dna_name = (current_dna_samples)$name %>% unique()
  assert(length(only_dna_name)==1)
  X2 = unique(rnaSamples$name)
  X1 = rep(only_dna_name, length(X2))
  dna_sample_map <- data.frame (X1,X2)
  left_join(dna_top %>% rename("DNA_name" = name), rna_joined, by = c("barcode", "architecture")) -> mpra_base
}

allDataFilteredJoined %>% 
  group_by(treatment, name) %>%
  summarise(batch = n()) %>%
  mutate(j = paste(treatment, name)) -> original_samples

original_samples %>%
  group_by(treatment) %>%
  summarise(batch = row_number()) -> batchNums

cbind(batchNums, original_samples %>% ungroup() %>% select(-batch, -treatment)) %>%
  select(name, treatment, batch) -> batch_info

left_join(mpra_base, batch_info) ->mpra_base
mpra_base %>% filter(!is.na(totalCounts)) -> mpra_base

mpra_base %>% 
  left_join(treatment_ids) %>%
  select(-treatment) %>%
  rename(treatment = t_id) -> mpra_base


#Get RNA depth factors
rna_joined %>% 
  group_by(name, treatment) %>%
  summarise(uq = quantile(totalCounts, .75)) %>% 
  left_join(batch_info) %>%
  mutate(norm = uq / first(.$uq)) %>%
  select(batch, norm, treatment, uq) %>%
  left_join(treatment_ids) %>%
  select(-treatment) %>% 
  rename(treatment = t_id) -> depth_factors

depth_factors %>% 
  left_join(dna_sample_map, by = c("name" = "X2")) %>%
  rename(dna_name = X1) ->  depth_factors
write_csv(depth_factors, paste0(PATH_TO_MPRA_INPUT, "rna_depth.csv"))


#Get DNA depth factors
dnaSamples %>% 
  group_by(name) %>%
  summarise(uq = quantile(totalCounts, .75)) %>%
  mutate(norm = uq / first(.$uq))-> dna_depth_factors
dna_depth_factors$run_name = RUN_NAME
write_csv(dna_depth_factors, paste0(PATH_TO_MPRA_INPUT, "dna_depth.csv"))


mpra_base %>%
  pivot_wider(values_from = totalCounts, id_cols = architecture, names_from = c(treatment, batch, bc), names_sep = ":") %>% 
  select(-architecture) %>%
  colnames()-> depth_cols


# Get depth columns
rna_depth = c()
dna_depth = c()
for (i in depth_cols){

  cur_condition = str_split(i, ":")[[1]][1]
  cur_batch = str_split(i, ":")[[1]][2]
  
  cur_sample_name = (depth_factors %>% filter(treatment == cur_condition, batch == cur_batch))$name
  cur_dna_sample = (dna_sample_map %>% filter(X2 == cur_sample_name))$X1
  cur_dna_depth = as.double((dna_depth_factors %>% filter(name == cur_dna_sample))$norm)
  dna_depth = append(dna_depth, cur_dna_depth)
  cur_depth_factor = as.double((depth_factors %>%
                                  filter(treatment == cur_condition, batch == as.character(cur_batch)))$norm)
  rna_depth = append(rna_depth, cur_depth_factor)
}
  
as.data.frame(depth_cols) -> col_ano
col_ano[c('conditions', 'batch', "barcode")] <- str_split_fixed(col_ano$depth_cols, ':', 3)
col_ano <- data.frame(col_ano%>%select(-depth_cols), row.names=col_ano$depth_cols)
col_ano$condition <- as.factor(col_ano$condition)
col_ano =col_ano %>%select(-conditions)
  
write_csv(col_ano, paste0(PATH_TO_MPRA_INPUT, "col_annotations.csv"))

#RNA
mpra_base %>%
  pivot_wider(values_from = totalCounts, id_cols = architecture, names_from = c(treatment, batch, bc), names_sep = ".") -> rna_mat

#DNA
mpra_base %>%
  pivot_wider(values_from = dna_count, id_cols = architecture, names_from = c(treatment, batch, bc), names_sep = ".") -> dna_mat

rna_mat[is.na(rna_mat)] <- 0
dna_mat[is.na(dna_mat)] <- 0

write_csv(as.data.frame(rna_depth), paste0(PATH_TO_MPRA_INPUT, "rna_depth_vals.csv"))
write_csv(as.data.frame(dna_depth), paste0(PATH_TO_MPRA_INPUT, "dna_depth_vals.csv"))

write_csv(rna_mat, paste0(PATH_TO_MPRA_INPUT, "rna_counts.csv"))
write_csv(dna_mat, paste0(PATH_TO_MPRA_INPUT, "dna_counts.csv"))
write_csv(treatment_ids, paste0(PATH_TO_MPRA_INPUT, "treatment_id.csv"))

