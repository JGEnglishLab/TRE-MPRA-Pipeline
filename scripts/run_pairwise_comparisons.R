library(tidyverse)
library(dplyr)
library(MPRAnalyze)
library(BiocParallel)


RESULTS_DIR = "pairwise_results/"

wd = getwd()
args = commandArgs(trailingOnly=TRUE)

param <- MulticoreParam(workers = args[1])
runs_dir = args[2]
base_treatment = args[3]
base_run = args[4]
stim_treatment = args[5]
stim_run = args[6]


###For testing
# RESULTS_DIR = "../../../mpra_final_data/pairwise_res_controls_correct_false"
# # RESULTS_DIR = "../../../mpra_final_data/pairwise_results/"
# comps = read_tsv("pairwise_comps.tsv")
# param <- MulticoreParam(workers = 6)
# runs_dir = "../../../mpra_final_data/"
# base_treatment = "NPFFR1-NP"
# base_run = "22482"
# stim_treatment = "NPFFR1mut-NP"
# stim_run = "22482"
###

if (base_run != stim_run){
  rbind(
    read_csv(paste0(runs_dir,base_run,"/MPRA_data.csv")) %>% mutate(run = base_run),
    read_csv(paste0(runs_dir,stim_run,"/MPRA_data.csv")) %>% mutate(run=stim_run)) %>% 
      filter((treatment == base_treatment & run == base_run ) | (run == stim_run & treatment == stim_treatment)) -> cur_comp_data
} else{
  read_csv(paste0(runs_dir,stim_run,"/MPRA_data.csv")) %>% mutate(run=stim_run) %>% 
    filter((treatment == base_treatment & run == base_run ) | (run == stim_run & treatment == stim_treatment)) -> cur_comp_data
}



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
  unique() %>%
  arrange(barcode, class, architecture, DNA_count)

dna2 = cur_comp_data %>%
  filter(treatment == stim_treatment, run == stim_run) %>%
  select(DNA_count, architecture, class, barcode) %>% 
  unique() %>%
  arrange(barcode, class, architecture, DNA_count)

same_DNA = identical(dna1,dna2)
dna1<-NULL
dna2<-NULL
gc() #return memory to the system


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


#Get the 
cur_comp_data %>% 
  group_by(architecture, treatment) %>% 
  summarise(DNA_barcodes = n(), RNA_barcodes = sum(!is.na(RNA_count))) %>% 
  pivot_wider(values_from = c(RNA_barcodes, DNA_barcodes), names_from = treatment) -> architecture_summary


#Make sure that all barcode numbers will be counted as barcodes per replicate
n_dna_base_name = paste0("DNA_barcodes_", base_treatment)
n_rna_base_name = paste0("RNA_barcodes_", base_treatment)

n_dna_stim_name = paste0("DNA_barcodes_", stim_treatment)
n_rna_stim_name = paste0("RNA_barcodes_", stim_treatment)

cur_comp_data %>% 
  group_by(treatment) %>%
  summarise(n_replicates = max(replicate_n)) -> n_replicates

#Divide the number of barcodes by the number of replicates for each condition
architecture_summary[n_dna_base_name] = architecture_summary[n_dna_base_name] / (n_replicates%>%filter(treatment == base_treatment))$n_replicates
architecture_summary[n_rna_base_name] = architecture_summary[n_rna_base_name] / (n_replicates%>%filter(treatment == base_treatment))$n_replicates
architecture_summary[n_dna_stim_name] = architecture_summary[n_dna_stim_name] / (n_replicates%>%filter(treatment == stim_treatment))$n_replicates
architecture_summary[n_rna_stim_name] = architecture_summary[n_rna_stim_name] / (n_replicates%>%filter(treatment == stim_treatment))$n_replicates


rna_mat <- NULL
dna_mat <- NULL
cur_comp_data <- NULL
gc() #Return memory to system

################################################################################
#                              RUN MPRAnalyze
################################################################################
obj <- MpraObject(dnaCounts = as.matrix(ready_dna_mat), rnaCounts = as.matrix(ready_rna_mat),
                  dnaAnnot = as.data.frame(dna_col_ano), rnaAnnot = as.data.frame(rna_col_ano), 
                  controls = controls, 
                  BPPARAM = param
                  )

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
                                correctControls=F,
                                reducedDesign = ~ 1)
} else{
  obj <- analyzeComparative(obj = obj,
                            dnaDesign = ~ barcode + condition,
                            rnaDesign = ~ condition,
                            BPPARAM = param,
                            correctControls=F,
                            reducedDesign = ~ 1)
}

lrt <- testLrt(obj)
lrt$architecture = row.names(lrt)
left_join(architecture_summary, lrt) -> lrt
lrt$log2FC = log2(exp(lrt$logFC))
lrt %>% select(-logFC) -> lrt
write_csv(lrt, paste0(RESULTS_DIR,base_treatment,"__", base_run, "_vs_",stim_treatment, "__", stim_run,'.csv'))



obj <- NULL
lrt <- NULL
ready_dna_mat <- NULL
ready_rna_mat <- NULL
gc() #Return memory to system


