# sam-tidy-bc-norm.R
# William Owens and Joshua Bloom and Sam Himes
library(tidyverse)
library(magrittr)
library(furrr)
library(stringdist)
library("ggpubr")


# filtering values
minimum_RNA_count <- 10
minimum_DNA_count <- 10

#In this version I will only work with barcodmap2
#barcode_map_fname <- "MPRA_barcode_map.tsv.gz"
barcode_map_fname_2 <- "MPRA_barcode_map_distance2.tsv.gz"
experiment_fname <- "MPRA_exp1_counts.tsv.gz"

barcode_map2 <- read_tsv(barcode_map_fname_2) 
#barcode_map <- read_tsv(barcode_map_fname)
experiment_data_o <- read_tsv(experiment_fname)

#filter to make sure all of barcodes are between 22 and 24 bps
#also, barcode_map2_filtered will only have barcodes with an occurrence >= 20
barcode_map2$barcodeLength = nchar(barcode_map2$barcode)
barcode_map2 %>% filter(barcodeLength >= 22 & barcodeLength <=24) -> barcode_map2

barcode_map2 %>%
  filter(occurrence <50 ) %>%
  ggplot(aes(y=occurrence)) +
  geom_histogram(binwidth = 1) + 
  coord_flip() +
  labs(title = "Barcode occurrance in barcode map (0,50)")

barcode_map2 %>% filter(occurrence >= 20) -> barcode_map2_filtered

# We tested 2 variations of the "library", full and bottlenecked. We also
# ran each biological replicate ("sample") on two different indexes, so we
# need to combine those technical replicates. You shouldn't have to do this.
experiment_data <-
  experiment_data_o %>% #Select for Bottlenecked library, and only ones treated with DNA, DMSO, and Forsk
  filter(library == "BN", treatment %in% c("DNA", "DMSO", "Forsk")) %>% #Collapse the replicates and sum the counts
  group_by(sample, treatment, barcode) %>%
  summarize(count = sum(count))

#***********************************************************************************************************
#                                 Run starcode on the barcode map
#***********************************************************************************************************

#write a file with just the barcodes to be run in starcode
barcode_map2_filtered$barcode -> barcode_map2_filtered_justBarCodes
write.table(barcode_map2_filtered_justBarCodes, "filteredBarcodeMap2.tsv", sep = "\t", col.names = FALSE, quote = F, row.names = F)

#Ran this command "starcode -i filteredBarcodeMap2.tsv --print-clusters -s -d 1 -o sc_out_filteredBarcodeMap2.tsv"
#Then I ran analyzeMappedSC.py to get the # of barcodes that are clusterd with a particular cluster and if each barcode is a centroid or not

#This is to check and see how many architechtures lose barccode representation by filtering for only centroids and > 20 occurence
#This is in comparison to the original barcode map 2
barcode_map2_filtered_clustered = read_tsv("analyzed_filteredBarcodeMap2.tsv")
left_join(barcode_map2, barcode_map2_filtered_clustered, by="barcode") -> barcode_map2_numCollapsed

#This gets the number of barcodes representing each architechtrue before and after filtering
left_join(barcode_map2_numCollapsed %>% 
            group_by(id, motif, period, promoter, spacer, class) %>%
            summarise(includingAll = n()),
          
          barcode_map2_numCollapsed %>% 
            filter(isCentroid == TRUE) %>%  #Only filter for barcodes that are the the centroid for now
            group_by(id, motif, period, promoter, spacer, class) %>%
            summarise(filteredBarcodeRep = n())
) -> architectureComparison

ggplot(data = pivot_longer(architectureComparison, cols = c(includingAll, filteredBarcodeRep)), 
       aes(x = value, fill = name, color = name)) + 
  geom_histogram(color="#e9ecef", alpha=0.5, binwidth = 50, position = "identity") + 
  labs(x = "Number of barcodes for each architecture", title = "Original barcodeMap2 vs.barcodeMap2 with\noccurrence>20 AND only centroids")

#Gives median barcode representation
median(filter(architectureComparison, !is.na(filteredBarcodeRep))$filteredBarcodeRep)

#This gives the number of architechtures that lose all representation
sum(is.na(architectureComparison$filteredBarcodeRep))

#Get the final filtered and cleaned barcode map
left_join(barcode_map2, barcode_map2_filtered_clustered, by="barcode") %>%
  filter(isCentroid == TRUE) %>%
  select(-barcodeLength, -numCollapsed, -isCentroid) -> finalBarcodeMap



#***********************************************************************************************************
#                                 Run starcode on the experimental data
#***********************************************************************************************************

#Seperate experiment data into each individual sample so they can all be collapsed in starcode
experiment_data %>% ungroup() %>% 
  filter(treatment == "Forsk", sample == "Sample7") %>% 
  group_by(barcode) %>%
  summarise(sum(count)) -> f7_collapse

experiment_data %>% ungroup() %>% 
  filter(treatment == "Forsk", sample == "Sample8") %>% 
  group_by(barcode) %>%
  summarise(sum(count)) -> f8_collapse

experiment_data %>% ungroup() %>% 
  filter(treatment == "DMSO", sample == "Sample4") %>% 
  group_by(barcode) %>%
  summarise(sum(count)) -> d4_collapse

experiment_data %>% ungroup() %>% 
  filter(treatment == "DMSO", sample == "Sample3") %>% 
  group_by(barcode) %>%
  summarise(sum(count)) -> d3_collapse

#FIXME!!!! ADDED
experiment_data %>% ungroup() %>% 
  filter(treatment == "DNA", sample == "Sample12") %>% 
  group_by(barcode) %>%
  summarise(sum(count)) -> dna12_collapse

#add a column to each to see if it is in the barcode map
f7_collapse$mapped = f7_collapse$barcode %in% finalBarcodeMap$barcode
f8_collapse$mapped = f8_collapse$barcode %in% finalBarcodeMap$barcode
d4_collapse$mapped = d4_collapse$barcode %in% finalBarcodeMap$barcode
d3_collapse$mapped = d3_collapse$barcode %in% finalBarcodeMap$barcode
#FIXME ADDED
dna12_collapse$mapped = dna12_collapse$barcode %in% finalBarcodeMap$barcode


write.table(f7_collapse, "f7_mapped_bcm2.tsv", sep = "\t", col.names = FALSE, quote = F, row.names = F)
write.table(f8_collapse, "f8_mapped_bcm2.tsv", sep = "\t", col.names = FALSE, quote = F, row.names = F)
write.table(d4_collapse, "d4_mapped_bcm2.tsv", sep = "\t", col.names = FALSE, quote = F, row.names = F)
write.table(d3_collapse, "d3_mapped_bcm2.tsv", sep = "\t", col.names = FALSE, quote = F, row.names = F)
#FIXME!!! ADDED
write.table(dna12_collapse, "dna12_mapped_bcm2.tsv", sep = "\t", col.names = FALSE, quote = F, row.names = F)

#Run this in the terminal
"for file in *_bcm2.tsv
do echo starcode -i $file --print-clusters -s -d 1 -o sc_out_${file}
starcode -i $file --print-clusters -s -d 1 -o sc_out_${file}
done"

#Then run python3 analyzeStarcode_revised.py 

#***********************************************************************************************************
#                                 open all of the collapsed experimental data
#                                 and visualize the outcome of the collapse
#***********************************************************************************************************

read_tsv("analyzed_out_f7_mapped_bcm2.tsv") -> f7_analyzed
read_tsv("analyzed_out_f8_mapped_bcm2.tsv") -> f8_analyzed
read_tsv("analyzed_out_d3_mapped_bcm2.tsv") -> d3_analyzed
read_tsv("analyzed_out_d4_mapped_bcm2.tsv") -> d4_analyzed
#FIXME!!! ADDED
read_tsv("analyzed_out_dna12_mapped_bcm2.tsv") -> dna12_analyzed


#add a column to each dataframe stating the type of collapse and the of the sampel
addType <- function(collapsedAnalyzed, inpName){
  collapsedAnalyzed$name = inpName
  collapsedAnalyzed$type = "NONE"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == TRUE & collapsedAnalyzed$collapsedInMapped == "NULL"] = "type1"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == FALSE & collapsedAnalyzed$collapsedInMapped == "NULL"] = "type2"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == TRUE & collapsedAnalyzed$collapsedInMapped != "NULL"] = "type3"
  collapsedAnalyzed$type[collapsedAnalyzed$centroidInMapped == FALSE & collapsedAnalyzed$collapsedInMapped != "NULL"] = "type4"
  
  return(collapsedAnalyzed)
}

addType(f7_analyzed, "Forskalin Sample 7") -> f7_analyzed
addType(f8_analyzed, "Forskalin Sample 8") -> f8_analyzed
addType(d3_analyzed, "DMSO Sample 3") -> d3_analyzed
addType(d4_analyzed, "DMSO Sample 4") -> d4_analyzed
#FIXME ADDED
addType(dna12_analyzed, "DNA Sample 12") -> dna12_analyzed


joined_experiment_data = rbind(f7_analyzed, f8_analyzed, d3_analyzed, d4_analyzed)
#FIXME!!! ADDED
joined_experiment_data = rbind(f7_analyzed, f8_analyzed, d3_analyzed, d4_analyzed, dna12_analyzed)


#This visualizes the types of collapses for each sample
ggplot(data = joined_experiment_data, aes(x = name, fill = type)) +
  geom_bar(position = position_stack(reverse = TRUE))+
  labs(x="Barcode Map", fill = "Experiment Collapse Type", title="Experimental Collapses")

#***********************************************************************************************************
#                               Filter for type 1 collapses
#                               Take the total count as the RNA count for each centroid
#                               Make new experiment data
#***********************************************************************************************************

#Get the stats on filtering for type 1
joined_experiment_data %>% group_by(name) %>%
  summarise(originalCount = sum(totalCollapsed), remaining = sum(type=="type1"), percent= remaining/originalCount)

#select only type 1 collapsed
joined_experiment_data %>% filter(type== "type1") -> joined_experiment_data_type1

#Adds Adam's stat
joined_experiment_data_type1 %>%
  mutate(
    collapse_stat = Inf,
    collapse_stat = if_else(totalCollapsed!=1, centroidCounts/totalCollapsed, collapse_stat)) -> joined_experiment_data_type1

#Plots adam's stat
joined_experiment_data_type1 %>%
  ggplot(aes(x = collapse_stat))+
  geom_histogram(bins = 60)+
  labs(x = "Centroid count / Number of barcodes collapsed") +
  facet_wrap(~ name, scales = "free")

joined_experiment_data_type1 %>%
  ggplot(aes(x = collapse_stat, y = totalCollapsed, color = medians))+
  geom_point()+
  labs(x = "Centroid count / Number of barcodes collapsed") +
  facet_wrap(~ name, scales = "free")+
  scale_color_gradientn(colours = rainbow(10)) 

joined_experiment_data_type1 %>%
  arrange(collapse_stat) %>%
  ggplot(aes(x = medians, y = totalCollapsed, color = collapse_stat))+
  geom_point()+
  labs(x = "medians") +
  facet_wrap(~ name, scales = "free")+
  scale_color_gradientn(colours = rainbow(50)) 

#Summarise Adam's stat
joined_experiment_data_type1 %>%
  filter(collapse_stat!= Inf)%>%
  group_by(name) %>%
  summarise(median_collapse_stat = median(collapse_stat),
            mean_collapse_stat = mean(collapse_stat))
  
#Get the increase in barcode count before/after collapse
ggplot(data = pivot_longer(joined_experiment_data_type1 %>% rename(dummyName = name), cols = c(centroidCounts, totalCounts)), 
       aes(y = value, x = name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=20, linetype="dashed", color = "red")+
  scale_y_continuous(trans='log10') +
  labs(y = "Barcode Counts (log10)",x="before/after collapse", title = "LD = 2, Barcode Counts Before (centroidCounts) and after (totalCounts) collapse")

#mean difference
mean(joined_experiment_data_type1$totalCounts - joined_experiment_data_type1$centroidCounts)

#median difference
median(joined_experiment_data_type1$totalCounts- joined_experiment_data_type1$centroidCounts)

#Get number rescued
sum(joined_experiment_data_type1$totalCounts >= 10 & joined_experiment_data_type1$centroidCounts < 10)

#Make the joined data_experiment data look like the original experiment data.
joined_experiment_data_type1 %>%
  select(centroid, totalCounts, name, collapse_stat) %>%
  rename("barcode" = centroid, "count" = totalCounts) %>%
  mutate(treatment = case_when(
    startsWith(name, "F") ~ "Forskalin",
    startsWith(name, "DMSO") ~ "DMSO", 
    startsWith(name, "DNA") ~ "DNA"
  )) %>%
  mutate(sample = case_when(
    endsWith(name, "7") ~ "Sample7",
    endsWith(name, "8") ~ "Sample8",
    endsWith(name, "3") ~ "Sample3",
    endsWith(name, "4") ~ "Sample4",
    endsWith(name, "12") ~ "Sample12"
  )) %>%
  select(-name) -> newExperimentData



#***********************************************************************************************************
#                                 Run Will's analysis with the new Files
#***********************************************************************************************************

# Map the barcodes observed in the count data to TREs based on the associations
# observed in the previous barcode mapping experiment.
mapped_data <-
  newExperimentData %>%
  left_join(finalBarcodeMap, by = "barcode")  

###############################################################################
# Quality control: how many times does a given barcode appear in each sample?
# Intuition holds that for low (< 10) values, Poisson effects will massively
# increase your noise.
# The dashed lines on these plots show the cutoffs for RNA/DNA count filters.

# only look at mapped reads for the QC
mapped_only <-
  mapped_data %>%
  filter(!is.na(id)) %>%
  mutate(
    cutoff = if_else(treatment == "DNA", minimum_DNA_count, minimum_RNA_count)
  )

mapped_only %>%  #Gets the raw counts of DNA/RNA_with_treatment by sample and treatment
  group_by(sample, treatment) %>%
  # cut off outliers that mess up the plot
  filter(count <= quantile(count, probs = 0.95)) %>%
  ggplot(aes(x = count)) +
  geom_histogram(binwidth = 5) +
  geom_vline(aes(xintercept = cutoff), linetype = "dashed") +
  labs(
    x = "barcode count", y = "frequency",
    title = "Per sample barcode count distributions with out Collapse"
  ) +
  facet_wrap(~ paste(sample, treatment))

#Get the number of barcodes above the cutoff
mapped_only %>% group_by(sample, treatment) %>% summarise(numAboveCutoff = sum(count > cutoff), avgCount = mean(count))


mapped_only %>% #Gets the number of counts per barcode in each sample
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), y = count)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = cutoff), linetype = "dashed") +
  scale_y_log10() +
  labs(x = "Sample (Treatment)", y = "Counts per barcode (log10)") +
  ggtitle("Total counts per mapped barcode (LD = 1)")

mapped_data %>% #Gets the counts per TRE 
  count(sample, treatment, id, motif, spacer, period, promoter) %>%
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), y = n)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(
    x = "Sample (Treatment)", y = "Number of Barcodes Per Architecture (log10)",
    title = "Barcodes per architecture (LD = 1)"
  )
#FIXME! I think that the above graph gets the # of barcodes per architecture. I think the name is slightly misleading 
#this is what labs was before
#
#labs(
#  x = "Sample (Treatment)", y = "Counts per TRE (log10)",
#  total = "Total counts per mapped TRE"
#I changed it to what it is now. 


mapped_data %>% count(sample, treatment, id, motif, spacer, period, promoter) -> bcPerA

#Gets the number of represented architechtures (across every set)
length(bcPerA$n)
#Gets avg barcode representation
mean(bcPerA$n)

# Computing RNA/DNA ratios
ratios <- 
  left_join(
    mapped_data %>%
      filter(treatment != "DNA") %>%
      rename(RNA_count = count),
    mapped_data %>%
      ungroup() %>%
      # only look at mapped reads
      filter(treatment == "DNA", !is.na(id)) %>%
      transmute(id, motif, spacer, period, promoter, barcode, DNA_count = count)
  ) %>%
  mutate(log_ratio = log2((RNA_count + 1) / (DNA_count + 1)), expression = RNA_count / DNA_count)

# Quality control: how many barcodes have DNA representation?
# At least for TREMPRArun1, this shows that the majority of barcodes without a
# DNA representation are just low-representation (count < 3), likely
# sequencing errors there are a few high-appearing orphans, but this can
# probably be solved by barcode mapping with levenshtein = 2 instead of 0. 
ratios %>% #Shows the counts per barcode for each sample with and without DNA representation
  mutate(has_dna = if_else(is.na(DNA_count), "no DNA representation", "DNA representation")) %>%
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), y = RNA_count, color = has_dna)) +
  geom_boxplot() +
  geom_hline(yintercept = minimum_RNA_count, linetype = "dashed") +
  scale_y_log10() +
  labs(
    x = "Sample (Treatment)", y = "RNA counts per barcode (log10)",
    title = "Total counts per barcode (with / without corresponding DNA rep.) (LD = 1)"
  )

ratios %>% 
  group_by(sample, treatment) %>% 
  summarise(numberMissingDNA = sum(is.na(DNA_count)), numberWithDNA = sum(!is.na(DNA_count)))


  

# Formally introduce filtering reasons
ratios_to_filter <-
  mutate(
    ratios,
    filter_reason = NA_character_,
    filter_reason = if_else(is.na(DNA_count) | DNA_count < minimum_DNA_count, "DNA low", filter_reason),
    filter_reason = if_else(RNA_count < minimum_RNA_count, "RNA low", filter_reason),
    filter_reason = if_else(RNA_count < minimum_RNA_count &(is.na(DNA_count) | DNA_count < minimum_DNA_count), "RNA & DNA low", filter_reason),
    filter_reason = if_else(is.na(id), "Unmapped barcode", filter_reason)
    # add other filters here (in order of ascending precedence)
  )

# Visualize why each each barcode is being filtered (Grouped by sample)
ratios_to_filter %>%
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), fill = filter_reason)) +
  geom_bar() +
  labs(
    x = "Sample (Treatment)", y = "Total barcodes",
    title = "Barcodes filtered (by sample / filter reason) (LD = 1)"
  )

ratios_to_filter %>%
  group_by(sample, treatment, filter_reason) %>%
  summarize(sum = sum(RNA_count)) %>%
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), y = sum, fill = filter_reason)) +
  geom_col() +
  labs(
    x = "Sample (Treatment)", y = "Total reads",
    title = "Reads filtered (by sample / filter reason) (LD = 1)"
  )

ratios_to_filter %>%
  group_by(treatment, sample, filter_reason) %>%
  summarise(sum = n())


# apply the filter
ratios_filtered <- filter(ratios_to_filter, is.na(filter_reason))

#Get's the number of barcodes remains
length(ratios_filtered$barcode)

###############################################################################
# Normalization

compute_raw_ratio <- function(df) {
  dmso_mean <- filter(df, treated == 0) %$% mean(log_ratio)
  treated_mean <- filter(df, treated == 1) %$% mean(log_ratio)
  treated_mean / dmso_mean
}

# Find the median RNA / DNA ratio of the spacer set, which should be
# representative of a non-responsive condition
spacer_stats <-
  filter(ratios_filtered, class == "spacer") %>%
  group_by(sample) %>%
  summarize(median_spacer = median(log_ratio), count_spacer = sum(RNA_count))


#I added this

#JUST NEGATIVE CONTROLS
negativeControls <-
  ratios_filtered %>%
  inner_join(spacer_stats) %>%
  mutate(
    treated = treatment != "DMSO"
  )%>%
  filter(class == "scramble" | class =="spacer") 

ggplot(data = negativeControls, aes(x = log_ratio, color = treated, fill = treated))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = .2, position = "identity")+
  facet_wrap(~ class)+
  labs(x = "\nExpression of each negative control\n\nexpression = log2((RNA_count + 1) / (DNA_count + 1))")

#Gets the difference of the means expression between treatment and not treated for scramble
mean((negativeControls %>% 
        filter(treated == T) %>%
        filter(class == "scramble"))$log_ratio) - 
mean((negativeControls %>% 
        filter(treated == F)%>%
        filter(class == "scramble"))$log_ratio)

#Gets the difference of the means expression between treatment and not treated for spacer
mean((negativeControls %>% 
        filter(treated == T) %>%
        filter(class == "spacer"))$log_ratio) - 
  mean((negativeControls %>% 
          filter(treated == F)%>%
          filter(class == "spacer"))$log_ratio)

#Gets the difference in expression before and after for matching barcodes
negativeControls %>%
  group_by(barcode) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  group_by(barcode, treated)%>%
  summarise(expression = mean(expression)) %>%
  pivot_wider(names_from = treated, values_from = expression, id_cols = c(treated, barcode))%>%
  mutate(diff = `TRUE`  - `FALSE`, fc =log2((`TRUE` + 1 )/ (`FALSE` + 1))) %>%
  filter(!is.na(diff)) -> negativeControlsMatchingBarcodes

negativeControlsMatchingBarcodes %>%
  ggplot(aes(x = fc))+
  geom_histogram(binwidth = .08)+
  geom_vline(xintercept = mean(negativeControlsMatchingBarcodes$fc), linetype="dashed", color = "red")+
  labs(x = "Difference in expression between treated and non-treated barcodes", title = "Negative Controls")

#Gets the difference in expression before and after for architectures
negativeControls %>%
  group_by(id, spacer, period, promoter) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  group_by(treated, motif, spacer, period, promoter, id)%>%
  summarise(expression = mean(expression)) %>%
  pivot_wider(names_from = treated, values_from = expression, id_cols = c(treated, spacer, period, promoter,id, motif))%>%
  mutate(diff = `TRUE`  - `FALSE`, fc =log2((`TRUE` + 1 )/ (`FALSE` + 1))) %>%
  filter(!is.na(diff)) -> negativeControlsArchitectures

negativeControlsArchitectures %>%
  ggplot(aes(x = fc))+
  geom_histogram(binwidth = .08)+
  geom_vline(xintercept = mean(negativeControlsArchitectures$fc), linetype="dashed", color = "red")+
  labs(x = "Difference in expression between treated and non-treated barcodes", title = "Negative Controls")


#CANDIDATE AND PROMEGA
test <-
  ratios_filtered %>%
  inner_join(spacer_stats) %>%
  mutate(
    treated = treatment != "DMSO")%>%
  filter(class == "candidate" | class == "promega") 

ggplot(data = test, aes(x = log_ratio, color = treated, fill = treated))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = .2, position = "identity")+
  labs(x = "\nExpression of each architecture\n\nexpression = log2((RNA_count + 1) / (DNA_count + 1))")

#Gets the difference of the means expression between treatment and not treated for scramble
mean((test %>% 
        filter(treated == T))$log_ratio) - 
  mean((test %>% 
          filter(treated == F))$log_ratio)


#Gets the difference in expression before and after treatment for individual barcodes
test %>%
  group_by(barcode) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  group_by(barcode, treated, motif)%>%
  summarise(expression = mean(expression)) %>%
  pivot_wider(names_from = treated, values_from = expression, id_cols = c(treated, barcode, motif))%>%
  mutate(diff = `TRUE` - `FALSE`, fc =log2((`TRUE` + 1 )/ (`FALSE` + 1))) %>%
  filter(!is.na(diff)) -> test2

#Plot the test for the individual barcodes
ggplot(data = test2, aes(x = fc, y = `FALSE`)) +
  geom_point() +
  geom_text_repel(data = test2[test2$fc>2,], 
                  aes(x=fc,y=`FALSE`,label=motif),
                  force =7,
                  max.overlaps = 390)+
  labs(y = "Expression (RNA_count / DNA_count) of the non-treated", 
       x = "Log 2 Fold change from non-treated expression to treated expression",
       title = "Expression of candidate/promega barcodes")


#Gets the difference in expression before and after treatment for architechtures
test %>%
  group_by(id, spacer, period, promoter) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  group_by(treated, motif, spacer, period, promoter, id)%>%
  summarise(expression = mean(expression), nBarcodes = n()) %>%
  pivot_wider(names_from = treated, values_from = c(expression,nBarcodes),  id_cols = c(treated, spacer, period, promoter,id, motif))%>%
  mutate(diff_expression = expression_TRUE - expression_FALSE, 
         fc =log2((expression_TRUE + 1 )/ (expression_FALSE + 1)), 
         diff_barcodes = nBarcodes_TRUE - nBarcodes_FALSE) %>%
  filter(!is.na(fc)) -> test3

#Gets the difference in expression before and after treatment for architechtures UPDATED 
test %>%
  group_by(id, spacer, period, promoter) %>%
  group_by(treated, motif, spacer, period, promoter, id)%>%
  summarise(expression = mean(expression), nBarcodes = n()) %>%
  pivot_wider(names_from = treated, values_from = c(expression,nBarcodes),  id_cols = c(treated, spacer, period, promoter,id, motif))%>%
  mutate(diff_expression = expression_TRUE - expression_FALSE, 
         fc =log2((expression_TRUE + 1 )/ (expression_FALSE + 1)), 
         diff_barcodes = nBarcodes_TRUE - nBarcodes_FALSE)  -> test3new


test3 %>%
  ggplot(aes(x = diff_barcodes)) + 
  geom_histogram(binwidth = 5) + 
  labs(x = "Difference in # of barcodes between treated and non-treated for the same architecture",
       title = "Candidate/promega architectures")



#Plot the test for the architectures
ggplot(data = test3, aes(x = fc, y = expression_FALSE, color = diff_barcodes)) +
  geom_point() +
  geom_text_repel(data = test3[test3$fc>2,], 
                  aes(x=fc,y=expression_FALSE,label=motif),
                  force =7,
                  max.overlaps = 60)+
  labs(y = "Expression (RNA_count / DNA_count) of the non-treated", 
       x = "Log 2 Fold change from non-treated expression to treated expression",
       title = "Expression of candidate/promega architectures")+
  scale_color_gradientn(colours = rainbow(100)) 



ggplot(data = test3, aes(x = diff_barcodes, y = fc))+
  geom_point()+
  labs(x = "Difference in of number barcodes between treated and non-treated\n (N barcodes treated) - (N barcodes not treated)", y = "Fold Change")

ggscatter(test3, x = "diff_barcodes", y = "fc", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          cor.coef.size = 8,
          xlab = "Difference in of number barcodes between treated and non-treated\n (N barcodes treated) - (N barcodes not treated)", 
          ylab = "Fold Change",
          title = "Spearman Correlation")

ggplot(data = test3, aes(x=expression_FALSE, y = expression_TRUE))+
  geom_point()+
  geom_text_repel(data = test3[test3$expression_TRUE>50 & test3$expression_FALSE < 25,], 
                               aes(x=expression_FALSE,y=expression_TRUE,label=motif),
                               force =7,
                               max.overlaps = 40)+
  labs(x = "Expression of non-treated", y = "Expression on treated", title = "Expression of candidate/promega architectures")





library(RColorBrewer)
library(wordcloud)

test2 %>%
  filter(fc > 5) %>%
  group_by(motif) %>%
  summarise(freq = n()) -> test3

test3 %>%
  ggplot(aes(y = freq, x = motif))+
  geom_boxplot()

#OR

wordcloud(words = test2$motif, freq = test3$freq, min.freq = 1,
          max.words=200, 
          random.order=FALSE, 
          rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))

test2 %>%
  ggplot(aes(x = fc))+
  geom_histogram(binwidth = .25)+
  geom_vline(xintercept = mean(test2$fc), linetype="dashed", color = "red")+
  labs(x = "Difference in expression between treated and non-treated barcodes", title = "Histogram ")



#ALL TOGETHER
allTest <-
  ratios_filtered %>%
  inner_join(spacer_stats) %>%
  mutate(
    treated = treatment != "DMSO") 

allTest %>%
  group_by(id, spacer, period, promoter) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  group_by(treated, motif, spacer, period, promoter, id)%>%
  summarise(expression = median(expression), nBarcodes = n()) %>% #FIXME! THIS IS WHERE I CHANGED TO MEDIAN
  pivot_wider(names_from = treated, values_from = c(expression,nBarcodes),  id_cols = c(treated, spacer, period, promoter,id, motif))%>%
  mutate(diff_expression = expression_TRUE - expression_FALSE, 
         fc =log2((expression_TRUE + 1 )/ (expression_FALSE + 1)), 
         diff_barcodes = nBarcodes_TRUE - nBarcodes_FALSE) %>%
  filter(!is.na(fc)) -> allTestArchitectures

ggplot(data = allTestArchitectures, aes(x = fc))+
  geom_histogram(binwidth = .05)+
  labs(x = "Log2 Fold change from treated to non-treated expression for all architectures \n(Spacer, Scramble, Candidate, and Promega)")


allTestArchitectures %>% 
  mutate(Zscore = (fc - mean(allTestArchitectures$fc)) /sd(allTestArchitectures$fc),
         PVal = 2*pnorm(abs(Zscore), lower.tail=FALSE, mean  = mean(allTestArchitectures$fc) , sd = sd(allTestArchitectures$fc)),
         PValCorrected = p.adjust(PVal,method="fdr",n=length(allTestArchitectures)),
         significant = (PValCorrected < .05)) -> t

ggplot(data = t, aes(x = PValCorrected))+
  geom_histogram(binwidth = .025)+
  geom_vline(xintercept = .05)



ggplot(data = t,(aes(y = expression_FALSE, x = fc, color = significant)))+
  geom_point()+
  labs(x = "Fold Change", 
       y = "Expression (RNA_count / DNA_count) of the non-treated", 
       color = "Corrected P.Value < .05?",
       title = "All architectures")

t %>% 
  filter(PValCorrected < .05)%>%
  group_by(motif) %>%
  summarise(freq = n()) %>%
  arrange(freq)%>%
  ggplot(aes(y = reorder(motif, -freq), x = freq))+
  geom_bar(stat = "identity")+
  labs(y = "Motif", 
       x = "Number of architectures that achieved significance (with correction)",
       title = "Using median expression to summarize treated/not-treated")
  

t %>% filter(grepl( "Spacer", motif, fixed = TRUE) | grepl( "Scramble", motif, fixed = TRUE)) %>% 
  ggplot(aes(x = fc)) +      
  geom_histogram(binwidth = .08)

wordcloud(words = s$motif, freq = (s %>% group_by(motif) %>% summarise(n()))$`n()`, min.freq = 1,
          max.words=200, 
          random.order=FALSE, 
          rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))

t %>% filter(significant == T) -> sigs
View(sigs)



#NEW AND IMPROVED!


#mafb testing
allTest %>% filter(spacer == "set1", period == 4, promoter == "minProm", id == 628, motif == "RARA") -> mafBtest
mafBtest %>% 
  group_by(barcode, treated, motif) %>%
  summarise(sd = sd(expression), expression = median(expression)) %>%
  pivot_wider(names_from = treated, values_from = c(expression, sd), id_cols = c(treated, barcode, motif)) -> mafBtest2
mafBtest2[is.na(mafBtest2$expression_FALSE),]$expression_FALSE = 0
mafBtest2 %>%
  mutate(expression_FALSE = expression_FALSE + .5,expression_TRUE = expression_TRUE + .5 , fc = log2(expression_TRUE/expression_FALSE)) -> mafBtest2

#Is it alwayw the untreated that we be na
allTest %>% 
  group_by(barcode, treated, motif) %>%
  summarise(sd = sd(expression), expression = median(expression)) %>%
  pivot_wider(names_from = treated, values_from = c(expression, sd), id_cols = c(treated, barcode, motif)) -> expTest













#Take a closer look at the spcer/scramble that achieved significance
allTest %>% filter(period == 4, spacer == "set2", promoter == "minTK", motif == "Spacer2") -> spacer2
spacer2 %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity")+
  labs(title = "period == 4, spacer == \"set2\", promoter == \"minTK\", motif == \"Spacer2\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 4, spacer == "set2", promoter == "minCMV", motif == "Spacer2") -> spacer2B
spacer2B %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 4, spacer == \"set2\", promoter == \"minCMV\", motif == \"Spacer2B\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 4, spacer == "set1", promoter == "minTK", motif == "Spacer3") -> spacer3
spacer3 %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 4, spacer == \"set1\", promoter == \"minTK\", motif == \"Spacer3\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 8, spacer == "set2", promoter == "minTK", motif == "Scramble8") -> scramble8
scramble8 %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 8, spacer == \"set2\", promoter == \"minTK\", motif == \"Scramble8\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 4, spacer == "set2", promoter == "minCMV", motif == "VENTX") -> ventx
ventx %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 4, spacer == \"set2\", promoter == \"minCMV\", motif == \"VENTX\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 4, spacer == "set1", promoter == "minCMV", motif == "CREB3") -> creb3_lowestP
creb3_lowestP %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 4, spacer == \"set1\", promoter == \"minCMV\", motif == \"CREB3\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 8, spacer == "set2", promoter == "minCMV", motif == "XBP1") -> XBP1_lowest
XBP1_lowest %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 8, spacer == \"set2\", promoter == \"minCMV\", motif == \"XBP1\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 1, spacer == "set1", promoter == "minCMV", motif == "Mafb") -> mafb_lowest
mafb_lowest %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 1, spacer == \"set1\", promoter == \"minCMV\", motif == \"Mafb\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 4, spacer == "set2", promoter == "minCMV", motif == "ATF7") -> atf7_lowest
atf7_lowest %>%
ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 4, spacer == \"set2\", promoter == \"minCMV\", motif == \"ATF7\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 4, spacer == "set1", promoter == "minCMV", motif == "GMEB2") -> gmeb2_lowest
gmeb2_lowest %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 4, spacer == \"set1\", promoter == \"minCMV\", motif == \"GMEB2\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 1, spacer == "set1", promoter == "minProm", motif == "GMEB2") -> gmeb2_highest
gmeb2_highest %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 4, spacer == \"set1\", promoter == \"minCMV\", motif == \"GMEB2 (Highest p.val)\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 1, spacer == "set2", promoter == "minTK", motif == "MYF6") -> myf6
myf6 %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 1, spacer == \"set2\", promoter == \"minTK\", motif == \"MYF6\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))


allTest %>% filter(period == 4, spacer == "set1", promoter == "minProm", motif == "NR2F1") -> NR2F1
NR2F1 %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 4, spacer == \"set1\", promoter == \"minProm\", motif == \"NR2F1\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 4, is.na(spacer) , promoter == "minCMV", motif == "Promega_CRE") -> promega_cre
promega_cre %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 5, position = "identity") +
  labs(title = "period == 4, is.na(spacer) , promoter == \"minCMV\", motif == \"Promega_CRE\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 8, is.na(spacer), promoter == "minTK", motif == "Promega_SIE") -> promega_sie
promega_sie %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 5, position = "identity") +
  labs(title = "period == 8, is.na(spacer), promoter == \"minTK\", motif == \"Promega_SIE\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))


allTest %>% filter(period == 8, spacer == "set1" , promoter == "minProm", motif == "TFCP2") -> tfpc2
tfpc2 %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 8, spacer == \"set1\" , promoter == \"minProm\", motif == \"TFCP2\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 8, spacer == "set1" , promoter == "minCMV", motif == "TBX4") -> tbx4
tbx4 %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 8, spacer == \"set1\" , promoter == \"minCMV\", motif == \"TBX4\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))

allTest %>% filter(period == 1, spacer == "set2" , promoter == "minCMV", motif == "YY2") -> yy2
yy2 %>%
  ggplot(aes(x = expression, fill = treatment, color = treatment))+
  geom_histogram(color="#e9ecef", alpha=0.45, binwidth = 1, position = "identity") +
  labs(title = "period == 1, spacer == \"set2\" , promoter == \"minCMV\", motif == \"YY2\"",
       y = "number of barcodes")+
  theme(text = element_text(size = 16))


scramble_stats <-
  filter(ratios_filtered, class == "scramble") %>%
  group_by(sample) %>%
  summarize(median_scramble = median(log_ratio), count_scramble = sum(RNA_count))

##
# My understanding of this. This is where they group my architecture. 

##
# Fit each TRE under each condition to a linear model:
# log_ratio ~ treated + median_spacer
linear_fits <-
  ratios_filtered %>%
  inner_join(spacer_stats) %>%
  mutate(
    # We just have one condition here, but you'll need to do some filtering /
    # joining to model different conditions to the same negative control.
    # Change this if you change the name of your control treatment.
    treated = as.numeric(treatment != "DMSO"),
  ) %>%
  group_by(id, motif, spacer, period, promoter, class) %>%
  group_nest() %>%
  mutate(
    fits = future_map(data, function(d) {
      # we're just pulling out the values for the "treated" term,
      # but feel free to use something different
      # lm(dependentVa ~ independentVar1 +independentVar2 etc.)
      # Should I add median scramble?
      broom::tidy((lm(log_ratio ~ treated + median_spacer, data = d)))
    }),
    # For comparison purposes, get the raw ratio of ratios
    raw_ratio = future_map_dbl(data, compute_raw_ratio)
  )

# Extract the estimate, std.error, statistic, and p-value of the "treated" term.
# The estimate
normalized <-
  linear_fits %>%
  unnest(fits) %>%
  filter(term == "treated") %>%
  transmute(
    id, motif, spacer, period, class, raw_ratio, 
    norm_ratio = estimate, std_error = std.error, z_score = statistic, p_value = p.value
  )

#CONCERNS
# 1 what they are calling the z_score is actually just the t value
# 2 Why is the median spacer being used as a predictor of log_ratio?

# Compute the "empirical range" based on the spacer and scramble controls.
# Any candidate TREs with scores in this range are probably unresponsive.
empirical_range <-
  normalized %>%
  filter(class == "spacer" | class == "scramble") %$%
  range(z_score)

normalized %>%
  pivot_longer(
    c(raw_ratio, norm_ratio),
    names_to = "normalization",
    values_to = "ratio"
  ) %>%
  ggplot(aes(x = class, y = ratio, color = normalization)) +
  geom_boxplot() +
  labs(y = "Estimated Forskolin log2(RNA/DNA) / DMSO log2(RNA/DNA)", title="LD = 1") +
  scale_y_continuous(limits = c(-5, 8))

# Export a TSV for further analysis
normalized %>%
  write_tsv(file = "MPRA_exp1_norm_collapsed.tsv.gz")



