setwd("/Volumes/external_disk/english_lab/mpra_pipeline")

library(tidyverse)
library(dplyr)
library(MPRAnalyze)
library(BiocParallel)
library(ggiraph)
library(jsonlite)

RESULTS_DIR = './comparative_results/'


id <- c(1, 2) #should just be a sequential list of numbers for the # of comparisons you have
names <- c("Name 1", "Name 2") # The string values of the comparison set
#Put all the comparisons in here 
#Seperate the name of each treatment with " # "
comp_list <- c(
  "20250_t_1 # 20250_t_2 # 20250_t_3 # 20250_t_5", 
  "19919_t_1 # 20250_t_2 # 20250_t_3") 
comparisons <- data.frame(id,names, comp_list)

all_treatments = c()
regex = "architecture|"
for (id_iter in comparisons$id){
  
  comparisons %>% 
    filter(id == id_iter) -> cur_row
  
  comps = cur_row$comp_list
  all_treatments = c(str_split(comp_list, " # ")[[1]], all_treatments)
  
  
  comps = paste0(str_replace_all(comps, " # ", "\\\\.|"), "\\.")
  regex = paste0(regex, comps)
}

rna_counts = read_csv("combined_rna_counts.csv")
dna_counts = read_csv("combined_dna_counts.csv")
depth_factor_data = read_csv("combined_depth_factor_data.csv")
mpra_base_data = read_csv("combined_mpra_base_data.csv")
treatment_ids = read_csv("combined_treatment_ids.csv")
dna_depth_df = read_csv("combined_dna_depth_df.csv")


grep(regex, colnames(rna_counts)) -> split_indices_rna
grep(regex, colnames(dna_counts)) -> split_indices_dna

rna_counts = rna_counts[split_indices_rna]
dna_counts = dna_counts[split_indices_dna]
mpra_base_data = mpra_base_data %>% filter(treatment %in% all_treatments)

#Make the normalization for depth factors RNA
#Get the normalizer
norm = depth_factor_data[1,]$uq
depth_factor_data$norm = depth_factor_data$uq / norm

#Make depth factor stuff for DNA 
norm = dna_depth_df[1,]$uq
dna_depth_df$norm = dna_depth_df$uq / norm

#Get bool vector of if each row is a control
controls = grepl('^Spacer|^Scramble', rna_counts$architecture)

keep_archi_order <- rna_counts
#Configure the row names

rna_counts <- column_to_rownames(rna_counts, var = "architecture")
dna_counts <- column_to_rownames(dna_counts, var = "architecture")


mpra_base_data %>%
  pivot_wider(values_from = totalCounts, id_cols = architecture, names_from = c(treatment, batch, bc), names_sep = ".") %>% 
  select(-architecture) %>%
  colnames()-> depth_cols


iter = 0
for (id_iter in comparisons$id){

  param <- BatchtoolsParam(workers = 24)
  


  iter = iter + 1
  print("__________")
  print("iter")
  print(id_iter)
  comparisons %>% 
    filter(id == id_iter) -> cur_row
  
  name = cur_row$name
  treatments = cur_row$comp_list
  
  regex = c()
  for (t in str_split(treatments, pattern = " # ")){
    for (i in t){
      regex = paste0(regex, i, "\\.|")
    }
  }



  name_csv = paste0(RESULTS_DIR,RUN_NAME,"__", name)
  name_csv = gsub(' ','_',name_csv)
  
  
  grep(regex, colnames(rna_counts)) -> split_indices
  cur_depth_cols <- depth_cols[split_indices]
  
  # Get depth columns
  rna_depth = c()
  dna_depth = c()
  for (i in cur_depth_cols){
    cur_run = str_split(i, "_",  n = 2)[[1]][1]
    sample_batch_info = str_split(i, "_",  n = 2)[[1]][2]
    cur_condition = str_split(i, "\\.", n = 2)[[1]][1]
    cur_batch = str_split(sample_batch_info, "\\.")[[1]][2]

    cur_depth_factor = as.double((depth_factor_data %>%
                                    filter(treatment == cur_condition, batch == as.character(cur_batch)))$norm)
    rna_depth = c(rna_depth, cur_depth_factor)
    
    cur_dna_depth =  (dna_depth_df %>% filter(run_name == cur_run))$norm
    dna_depth = c(dna_depth, cur_dna_depth)
    
  }

  as.data.frame(cur_depth_cols) -> new_col_ano
  new_col_ano[c('conditions', 'batch', "barcode")] <- str_split_fixed(new_col_ano$cur_depth_cols, '\\.', 3)
  new_col_ano <- data.frame(new_col_ano%>%select(-cur_depth_cols), row.names=new_col_ano$cur_depth_cols)
  new_col_ano$condition <- as.factor(new_col_ano$condition)
  new_col_ano=new_col_ano %>%select(-conditions)
  complete_obj <- MpraObject(dnaCounts = as.matrix(dna_counts[split_indices]),
                             rnaCounts = as.matrix(rna_counts[split_indices]),
                             dnaAnnot = as.data.frame(new_col_ano),
                             rnaAnnot = as.data.frame(new_col_ano),
                             controls = controls)


 
  #We will have to fix this to get the dna depth for two comparisons..
  complete_obj = setDepthFactors(complete_obj, dnaDepth = dna_depth, rnaDepth = rna_depth)
  

  if (length(unique(new_col_ano$batch)) > 1){ #If more than one batch is present us it in comparative
    comp_obj <- analyzeComparative(obj = complete_obj,
                              dnaDesign = ~ barcode + batch + condition,
                              rnaDesign = ~ condition,
                              reducedDesign = ~ 1, 
                              BPPARAM = param)
  } else{
    comp_obj <- analyzeComparative(obj = complete_obj, 
                                   dnaDesign = ~ barcode + condition, 
                                   rnaDesign = ~ condition, 
                                   reducedDesign = ~ 1,
                                   BPPARAM = param)
  }
  
  res <- testLrt(comp_obj)
  res$architecture = row.names(res)

  write_csv(res, paste0(name_csv,'.csv'))
  
}

