library(tidyverse)
library(dplyr)
library(MPRAnalyze)
library(BiocParallel)

RESULTS_DIR = "pairwise_results/"

print(getwd())
wd = getwd()
args = commandArgs(trailingOnly=TRUE)

comps = read_tsv(args[1])
param <- BatchtoolsParam(workers = args[2])
runs_dir = args[3]


all_runs = unique(c(comps$base_run, comps$stim_run))

rna_counts = data.frame()
dna_counts = data.frame()
treatment_ids = data.frame()
rna_depth_factor = data.frame()
dna_depth_factor = data.frame()
col_ano = data.frame()
# depth_cols = data.frame()

count = 0
for (i in all_runs){
  count = count + 1
  if (count == 1){ #If its the first iteration, just initialize the dataframe
    rna_counts = read_csv(paste0(runs_dir,i,"/mpra_input/rna_counts.csv"))
    dna_counts = read_csv(paste0(runs_dir,i,"/mpra_input/dna_counts.csv"))
    treatment_ids = read_csv(paste0(runs_dir,i,"/mpra_input/treatment_id.csv"))
    rna_depth_factor = read_csv(paste0(runs_dir,i,"/mpra_input/rna_depth.csv"))
    dna_depth_factor = read_csv(paste0(runs_dir,i,"/mpra_input/dna_depth.csv"))
    # depth_cols = read_csv(paste0("runs/",i,"/mpra_input/depth_cols.csv"))
    col_ano = read_csv(paste0(runs_dir,i,"/mpra_input/col_annotations.csv"))
    
  }
  else{ #If we are beyond the first iteration, join the data frames
    rna_counts = full_join(read_csv(paste0(runs_dir,i,"/mpra_input/rna_counts.csv")), rna_counts)
    dna_counts = full_join(read_csv(paste0(runs_dir,i,"/mpra_input/dna_counts.csv")), dna_counts)
    treatment_ids = rbind(read_csv(paste0(runs_dir,i,"/mpra_input/treatment_id.csv")), treatment_ids)
    rna_depth_factor = rbind(read_csv(paste0(runs_dir,i,"/mpra_input/rna_depth.csv")), rna_depth_factor)
    dna_depth_factor = rbind(read_csv(paste0(runs_dir,i,"/mpra_input/dna_depth.csv")), dna_depth_factor)
    # depth_cols = rbind(read_csv(paste0("runs/",i,"/mpra_input/depth_cols.csv")), depth_cols)
    col_ano = rbind(read_csv(paste0(runs_dir,i,"/mpra_input/col_annotations.csv")), col_ano)
    
  }
}

#Re-calculate the RNA and DNA depth factors now that we've combined 
rna_depth_factor %>%
  mutate(norm = uq / first(.$uq)) -> rna_depth_factor

dna_depth_factor %>%
  mutate(norm = uq / first(.$uq)) -> dna_depth_factor


#Get bool vector of if each row is a control
controls = grepl('^Spacer|^Scramble', rna_counts$architecture)

stopifnot(nrow(dna_counts) == nrow(rna_counts))
architecture_order = dna_counts$architecture

#Configure the row names
dna_counts %>% remove_rownames %>% column_to_rownames(var="architecture") -> dna_counts
rna_counts %>% remove_rownames %>% column_to_rownames(var="architecture") -> rna_counts

for (id_iter in comps$id){

  comps %>% 
    filter(id == id_iter) -> cur_row
  
  base_treatment = cur_row$base_treatment
  base_run = cur_row$base_run
  
  stim_treatment = cur_row$stim_treatment
  stim_run = cur_row$stim_run
  
  base_t_id = (treatment_ids %>% filter(run == base_run, treatment == base_treatment))$t_id
  stim_t_id = (treatment_ids %>% filter(run == stim_run, treatment == stim_treatment))$t_id
  regex = paste0(base_t_id,"\\.|",stim_t_id,"\\.") 
  
  name_csv = paste0(RESULTS_DIR, base_treatment, "__",base_run ,"_vs_",stim_treatment, "__",stim_run)
  name_csv = gsub(' ','_',name_csv)
  
  grep(regex, colnames(rna_counts)) -> split_indices
  cur_depth_cols <- colnames(rna_counts[split_indices])
  
  # Get depth columns
  rna_depth = c()
  dna_depth = c()
  for (i in cur_depth_cols){
    cur_run = str_split(i, "_",  n = 2)[[1]][1]
    sample_batch_info = str_split(i, "_",  n = 2)[[1]][2]
    cur_condition = str_split(i, "\\.", n = 2)[[1]][1]
    cur_batch = str_split(sample_batch_info, "\\.")[[1]][2]
    
    cur_rna_depth = as.double((rna_depth_factor %>%
                                    filter(treatment == cur_condition, batch == as.character(cur_batch)))$norm)
    cur_dna_name = (rna_depth_factor %>%
                                 filter(treatment == cur_condition, batch == as.character(cur_batch)))$dna_name
    rna_depth = c(rna_depth, cur_rna_depth)
    
    cur_dna_depth =  (dna_depth_factor %>% filter(run_name == cur_run, name == cur_dna_name))$norm
    dna_depth = c(dna_depth, cur_dna_depth)
  }
  
  as.data.frame(cur_depth_cols) -> new_col_ano
  new_col_ano[c('conditions', 'batch', "barcode")] <- str_split_fixed(new_col_ano$cur_depth_cols, '\\.', 3)
  new_col_ano <- data.frame(new_col_ano%>%select(-cur_depth_cols), row.names=new_col_ano$cur_depth_cols)
  new_col_ano$condition <- as.factor(new_col_ano$condition)
  new_col_ano$condition <- relevel(new_col_ano$condition, base_t_id) #The condition must be a factor with base first
  new_col_ano=new_col_ano %>%select(-conditions)
  complete_obj <- MpraObject(dnaCounts = as.matrix(dna_counts[split_indices]),
                             rnaCounts = as.matrix(rna_counts[split_indices]),
                             dnaAnnot = as.data.frame(new_col_ano),
                             rnaAnnot = as.data.frame(new_col_ano),
                             controls = controls)
  
  
  
  #We will have to fix this to get the dna depth for two comparisons..
  complete_obj = setDepthFactors(complete_obj, dnaDepth = dna_depth, rnaDepth = rna_depth)
  
  
  if (length(unique(new_col_ano$batch)) > 1){ #If more than one batch is present us it in comparison
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
  
  #Get the architectures that are completely missing from each treatment type
  rna_counts[split_indices] %>% 
    rownames_to_column(var="architecture")  %>% 
    pivot_longer(cols = -architecture) %>% 
    mutate(treatment = str_extract(name, "^[^.]+"))   %>%
    select(-name) %>%
    group_by(architecture, treatment) %>% 
    summarise(sum = sum(value))  %>%
    filter(sum == 0) %>% 
    left_join(treatment_ids, by=c("treatment" = "t_id")) %>%
    select(-treatment, -sum) %>%
    rename(treatment = treatment.y) %>%
    arrange(treatment, run) -> missing_architectures
  
  write_csv(missing_architectures, paste0(name_csv,'__missing_architectures.csv'))
  
  
}
