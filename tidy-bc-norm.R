# tidy-bc-norm.R
# William Owens and Joshua Bloom
library(tidyverse)
library(magrittr)
# for parallel linear modeling
library(furrr)

# filtering values
minimum_RNA_count <- 10
minimum_DNA_count <- 10

barcode_map_fname <- "MPRA_barcode_map.tsv.gz"
experiment_fname <- "MPRA_exp1_counts.tsv.gz"

barcode_map <- read_tsv(barcode_map_fname)
experiment_data <- read_tsv(experiment_fname)

# We tested 2 variations of the "library", full and bottlenecked. We also
# ran each biological replicate ("sample") on two different indexes, so we
# need to combine those technical replicates. You shouldn't have to do this.
experiment_data <-
  experiment_data %>%
  filter(library == "BN", treatment %in% c("DNA", "DMSO", "Forsk")) %>%
  group_by(sample, treatment, barcode) %>%
  summarize(count = sum(count))

# Map the barcodes observed in the count data to TREs based on the associations
# observed in the previous barcode mapping experiment.
mapped_data <-
  experiment_data %>%
  left_join(barcode_map, by = "barcode")

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

mapped_only %>%
  group_by(sample, treatment) %>%
  # cut off outliers that mess up the plot
  filter(count <= quantile(count, probs = 0.95)) %>%
  ggplot(aes(x = count)) +
  geom_histogram(binwidth = 5) +
  geom_vline(aes(xintercept = cutoff), linetype = "dashed") +
  labs(
    x = "barcode count", y = "frequency",
    title = "Per sample barcode count distributions"
  ) +
  facet_wrap(~ paste(sample, treatment))

mapped_only %>%
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), y = count)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = cutoff), linetype = "dashed") +
  scale_y_log10() +
  labs(x = "Sample (Treatment)", y = "Counts per barcode") +
  ggtitle("Total counts per mapped barcode")

mapped_data %>%
  count(sample, treatment, id, motif, spacer, period, promoter) %>%
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), y = n)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(
    x = "Sample (Treatment)", y = "Counts per TRE",
    total = "Total counts per mapped TRE"
  )

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
  mutate(log_ratio = log2((RNA_count + 1) / (DNA_count + 1)))

# Quality control: how many barcodes have DNA representation?
# At least for TREMPRArun1, this shows that the majority of barcodes without a
# DNA representation are just low-representation (count < 3), likely
# sequencing errors there are a few high-appearing orphans, but this can
# probably be solved by barcode mapping with levenshtein = 2 instead of 0.
ratios %>%
  mutate(has_dna = if_else(is.na(DNA_count), "no DNA representation", "DNA representation")) %>%
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), y = RNA_count, color = has_dna)) +
  geom_boxplot() +
  geom_hline(yintercept = minimum_RNA_count, linetype = "dashed") +
  scale_y_log10() +
  labs(
    x = "Sample (Treatment)", y = "Counts per barcode",
    title = "Total counts per barcode (with / without corresponding DNA rep.)"
  )

# Formally introduce filtering reasons
ratios_to_filter <-
  mutate(
    ratios,
    filter_reason = NA_character_,
    filter_reason = if_else(is.na(DNA_count) | DNA_count < 10, "DNA low", filter_reason),
    filter_reason = if_else(RNA_count < 10, "RNA low", filter_reason),
    filter_reason = if_else(is.na(id), "Unmapped barcode", filter_reason)
    # add other filters here (in order of ascending precedence)
  )

# Visualize why each each barcode / read is being filtered
ratios_to_filter %>%
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), fill = filter_reason)) +
  geom_bar() +
  labs(
    x = "Sample (Treatment)", y = "Total barcodes",
    title = "Barcodes filtered (by sample / filter reason)"
  )

ratios_to_filter %>%
  group_by(sample, treatment, filter_reason) %>%
  summarize(sum = sum(RNA_count)) %>%
  ggplot(aes(x = sprintf("%s\n(%s)", sample, treatment), y = sum, fill = filter_reason)) +
  geom_col() +
  labs(
    x = "Sample (Treatment)", y = "Total reads",
    title = "Reads filtered (by sample / filter reason)"
  )

# apply the filter
ratios_filtered <- filter(ratios_to_filter, is.na(filter_reason))

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
  labs(y = "Estimated Forskolin log2(RNA/DNA) / DMSO log2(RNA/DNA)") +
  scale_y_continuous(limits = c(-5, 8))

# Export a TSV for further analysis
normalized %>%
  write_tsv(file = "MPRA_exp1_norm.tsv.gz")
