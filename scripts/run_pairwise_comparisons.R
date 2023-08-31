library(tidyverse)
library(dplyr)
library(MPRAnalyze)
library(BiocParallel)


RESULTS_DIR = "pairwise_results/"

wd = getwd()
args = commandArgs(trailingOnly=TRUE)

comps = read_tsv(args[1])
param <- BatchtoolsParam(workers = args[2])
runs_dir = args[3]

#For testing
# comps = read_tsv("../entered_pairwise_comparisons/pairwise_comparisons_30-08-2023_14-29.tsv")
# param <- BatchtoolsParam(workers = 6)
# runs_dir = "../runs/"
##

all_runs = unique(c(comps$base_run, comps$stim_run))

all_run_data = data.frame()

count = 0
#Read in all the data from all the runs
for (i in all_runs){
  count = count + 1
  if (count == 1){ #If its the first iteration, just initialize the dataframe
    all_run_data = read_csv(paste0(runs_dir,i,"/MPRA_data.csv")) %>% mutate(run = i)
    
  }
  else{ #If we are beyond the first iteration, join the data frames
    all_run_data = rbind(read_csv(paste0(runs_dir,i,"/MPRA_data.csv")) %>% mutate(run=i), all_run_data)
  }
}



for (comp_id in comps$id){

  comps %>% 
    filter(id == comp_id) -> cur_row
  
  base_treatment = cur_row$base_treatment
  base_run = cur_row$base_run
  
  stim_treatment = cur_row$stim_treatment
  stim_run = cur_row$stim_run
  
  #Filter the data for just the two treatments in the comparison
  all_run_data %>% 
    filter((treatment == base_treatment | treatment == stim_treatment) &(run == stim_run | run == base_run)) -> cur_comp_data
  
  #Get rna replicate numbers
  cur_comp_data %>% ungroup() %>% 
    select(RNA_sample_number, treatment, run) %>% 
    unique() %>%
    group_by(treatment, run) %>% 
    reframe(replicate_n = row_number(), RNA_sample_number = RNA_sample_number) %>%
    right_join(cur_comp_data) -> cur_comp_data
  
  #Get treatment numbers
  treatment = c(base_treatment, stim_treatment)
  run = c(base_run, stim_run)
  treatment_n = c(1,2) #There will only be two treatments, and we want the base treatment to be 1
  treatment_numbers = tibble(treatment, run, treatment_n)
  
  treatment_numbers %>%
    right_join(cur_comp_data)  -> cur_comp_data 
  
  #Check to see if DNA from the two treatments is the same
  dna1 = cur_comp_data %>%
    filter(treatment == base_treatment, run == base_run) %>%
    select(DNA_count, architecture, class, barcode) %>%
    arrange(DNA_count, architecture, class, barcode)
  
  dna2 = cur_comp_data %>%
    filter(treatment == stim_treatment, run == stim_run) %>%
    select(DNA_count, DNA_rpm, architecture, class, barcode) %>%
    arrange(DNA_count, DNA_rpm, architecture, class, barcode)
  
  same_DNA = identical(dna1,dna2)
  
  if (same_DNA){
    #Pivot all data to get DNA_mat for MPRAnalyze
    cur_comp_data %>% 
      select(DNA_count, architecture, barcode_n, class) %>% 
      unique() %>%
      pivot_wider(id_cols = c("architecture", "class"), 
                  values_from = c("DNA_count"), 
                  names_from = c("barcode_n"),
                  names_sep = ":")-> dna_mat 
    
   
  } else{
    
    cur_comp_data %>% 
      select(DNA_count, architecture, barcode_n, class, treatment_n) %>% 
      unique() %>%
      pivot_wider(id_cols = c("architecture", "class"), 
                  values_from = c("DNA_count"), 
                  names_from = c("barcode_n", "treatment_n"),
                  names_sep = ":")-> dna_mat 
    
  }
  
  #Pivot all data to get RNA_mat for MPRAnalyze
  cur_comp_data %>% 
    select(RNA_count, architecture, barcode_n, replicate_n, treatment_n, class) %>%
    pivot_wider(id_cols = c("architecture", "class"), 
                values_from = c("RNA_count"), 
                names_from = c("barcode_n", "replicate_n", "treatment_n"), 
                names_sep = ":")-> rna_mat  
  
  
  #Extract control booleans
  controls = rna_mat$class == "scramble" | rna_mat$class == "spacer"
  
  #Make the rna col annotations
  cur_comp_data %>% 
    select(barcode_n, replicate_n, treatment_n) %>% 
    unique() -> rna_col_ano
  
  rna_col_ano$barcode = as.factor(rna_col_ano$barcode_n)
  rna_col_ano$batch = as.factor(rna_col_ano$replicate_n)
  rna_col_ano$condition = as.factor(rna_col_ano$treatment_n)

  rna_col_ano %>% select(barcode,batch,condition) -> rna_col_ano
  
  #Make dna col annotations
  if(same_DNA){
    dna_col_ano = as.data.frame(cur_comp_data %>% select(barcode_n) %>% unique())
    dna_col_ano$barcode = as.factor(dna_col_ano$barcode_n)
    dna_col_ano %>% select(barcode) -> dna_col_ano
  } else{
    cur_comp_data %>% 
      select(barcode_n, treatment_n) %>% 
      unique() -> dna_col_ano
    
    dna_col_ano$barcode = as.factor(dna_col_ano$barcode_n)
    dna_col_ano$condition = as.factor(dna_col_ano$treatment_n)
    
    dna_col_ano %>% select(barcode,condition) -> dna_col_ano
  }
  
  rna_mat %>% select(-class) %>% column_to_rownames("architecture") -> ready_rna_mat
  dna_mat %>% select(-class) %>% column_to_rownames("architecture") -> ready_dna_mat
  
  #Replace all NAs with 0
  ready_rna_mat[is.na(ready_rna_mat)] <- 0
  ready_dna_mat[is.na(ready_dna_mat)] <- 0
  
  
  ################################################################################
  #                              RUN MPRAnalyze
  ################################################################################
  obj <- MpraObject(dnaCounts = as.matrix(head(ready_dna_mat)), rnaCounts = as.matrix(head(ready_rna_mat)), 
                    dnaAnnot = as.data.frame(dna_col_ano), rnaAnnot = as.data.frame(rna_col_ano), 
                    controls = controls, 
                    BPPARAM = param)
  
  obj <- estimateDepthFactors(obj, lib.factor = c("condition", "batch"),
                              which.lib = "rna", 
                              depth.estimator = "uq")
  
  if (!same_DNA){
    obj <- estimateDepthFactors(obj, lib.factor = c("condition"),
                                which.lib = "dna", 
                                depth.estimator = "uq")
  }
  if (same_DNA){
    obj <- analyzeComparative(obj = obj,
                                  dnaDesign = ~ barcode,
                                  rnaDesign = ~ condition,
                                  BPPARAM = param,
                                  reducedDesign = ~ 1)
  } else{
    obj <- analyzeComparative(obj = obj,
                              dnaDesign = ~ barcode + condition,
                              rnaDesign = ~ condition,
                              BPPARAM = param,
                              reducedDesign = ~ 1)
  }
  
  lrt <- testLrt(obj)
  lrt$architecture = row.names(lrt)
  
  cur_comp_data %>% 
    group_by(architecture, treatment) %>% 
    summarise(DNA_barcodes = n(), RNA_barcodes = sum(!is.na(RNA_count))) %>% 
    pivot_wider(values_from = c(RNA_barcodes, DNA_barcodes), names_from = treatment, names_prefix = "n_") -> architecture_summary
  
    left_join(architecture_summary, lrt) -> lrt
    
    write_csv(lrt, paste0(RESULTS_DIR,base_treatment,"__", base_run, "_vs_",stim_treatment, "__", stim_run,'.csv'))
    
}
