library("tidyverse")
library("MPRAnalyze")
library("dplyr")
library("BiocParallel")
library("testit")

wd = getwd()

# Get the command-line arguments.
args = commandArgs(trailingOnly=TRUE)

####
#DELETE AND CHANGE LATER!
# setwd("/Volumes/external_disk/english_lab/TRE-MPRA/runs/run1")
# rna_counts = read_csv("./mpra_input/rna_counts.csv")
# dna_counts = read_csv("./mpra_input/dna_counts.csv")
# treatment_ids = read_csv("./mpra_input/treatment_id.csv")
# col_ano = read_csv("./mpra_input/col_annotations.csv")
# rna_depth = read_csv("./mpra_input/rna_depth_vals.csv")$rna_depth
# dna_depth = read_csv("./mpra_input/dna_depth_vals.csv")$dna_depth

rna_counts = read_csv(args[1])
dna_counts = read_csv(args[2])
treatment_ids = read_csv(args[3])
col_ano = read_csv(args[4])
rna_depth = read_csv(args[5])$rna_depth
dna_depth = read_csv(args[6])$dna_depth

controls = grepl('^Spacer|^Scramble', rna_counts$architecture)

param <- BatchtoolsParam(workers = 5)

assert("RNA and DNA counts same length", nrow(dna_counts) == nrow(rna_counts))
architecture_order = dna_counts$architecture

#Remove architecture and set it as row name
dna_counts %>% remove_rownames %>% column_to_rownames(var="architecture") -> dna_counts
rna_counts %>% remove_rownames %>% column_to_rownames(var="architecture") -> rna_counts

col_ano %>% mutate(rowName = paste(condition, batch, barcode, sep = ":")) -> col_ano
col_ano <- data.frame(col_ano%>%select(-rowName), row.names=col_ano$rowName)
col_ano$condition <- as.factor(col_ano$condition)

complete_obj <- MpraObject(dnaCounts = as.matrix(dna_counts),
                           rnaCounts = as.matrix(rna_counts),
                           dnaAnnot = as.data.frame(col_ano),
                           rnaAnnot = as.data.frame(col_ano),
                           controls = controls)

complete_obj = setDepthFactors(complete_obj, dnaDepth = dna_depth, rnaDepth = rna_depth)

quant_obj <- analyzeQuantification(obj = complete_obj,
                                   dnaDesign = ~ batch + condition,
                                   rnaDesign = ~ condition,
                                   BPPARAM = param)


alpha <- getAlpha(quant_obj, by.factor = "condition")

#Get empircal results
all_emp_results = alpha %>% rownames_to_column(var = "architecture") %>% select(architecture)
for (col in names(alpha)){
  if (sum(is.na(alpha[[col]])) != length(alpha[[col]])){
    emp_res <- testEmpirical(obj = quant_obj,statistic = alpha[[col]])
    colnames(emp_res) <- paste(col, colnames(emp_res), sep = "_")
    all_emp_results = cbind(all_emp_results, emp_res)
  }
}

#All_emp_result has column names of the following format t_1_(result)
#This function we replace t_1 (the treatment id) with the actual treatment name
modify_col_names <- function(cols_to_modify) {
  mod_cols = c()
  for (c in cols_to_modify){
    if (c != "architecture"){
      run_name = str_split(c, "_")[[1]][1]
      treatment_id = paste0(run_name,"_t_", str_split(c, "_")[[1]][3])
      column_type = str_split(c, "_")[[1]][4]
      treatment = (treatment_ids %>% filter(t_id == treatment_id))$treatment
      
      new_col = paste(column_type, treatment, sep = "_")
      mod_cols = c(mod_cols, new_col)
    }
    else{
      mod_cols = c(mod_cols, c)
    }
  }
  return(mod_cols)
}	

colnames(all_emp_results) <- modify_col_names(colnames(all_emp_results))
new_order = sort(colnames(all_emp_results))
all_emp_results_correct_names <- all_emp_results[, new_order]

#The all_emp_results comes with multiple control columns. These are redundent.
#This creates one controls column
all_emp_results_correct_names %>%
  rename("controls_cols" = matches("^control_*")) %>%
  rename("controls" = "controls_cols1")%>%
  select(-matches("^controls_cols*")) -> all_emp_results_correct_names


write_csv(all_emp_results_correct_names, "empirical_results.csv")
