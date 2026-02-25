# Compare library build methods: SCR vs BEST

# Set working directory ----
#setwd("")

# Load packages ----
library(tidyverse)
library(jsonlite)

# 1. Load data ----
## 1.2. Stats (paleomix, noMicro, damageProfiler read length averages)
libs_stats <- read_csv("paleomix_libCompare_noMicro.csv", col_names = T)
# 1.3. Read length distributions (damageProfiler)
# 1.3.1. Get list of all json files with length distribution info to process (do once, then load at step 1.4)
json_files <- list.files(path = "damageProfiler",
                        pattern = "dmgprof\\.json$",
                        recursive = TRUE,
                        full.names = TRUE)
# 1.3.2. Define function to get read length distribution from jsons
process_dp_json <- function(json_files) {
  
  dp <- fromJSON(json_files)
  
  folder <- basename(dirname(normalizePath(json_files)))
  SeqID  <- str_remove(folder, "\\.damageProfiler$")
  SampleID <- str_extract(SeqID, "^[^_]+")
  LibraryID <- str_remove(SeqID, "^[^_]+_")
  
  # Convert forward and reverse length distributions
  lendist_fw <- enframe(dp$lendist_fw, name = "Length", value = "Occurrences")
  lendist_rv <- enframe(dp$lendist_rv, name = "Length", value = "Occurrences")
  
  lendist <- bind_rows(lendist_fw, lendist_rv) %>%
    mutate(Length = as.numeric(Length),
           Occurrences = as.numeric(unlist(Occurrences))) %>%
    group_by(Length) %>%
    summarise(Occurrences = sum(Occurrences), .groups = "drop") %>%
    mutate(SeqID = SeqID) %>%
    mutate(SampleID = SampleID) %>%
    mutate(LibraryID = LibraryID)
}

# 1.3.3. Process jsons
dp_lg_results <- map(json_files, process_dp_json)
# 1.3.4. Store all length distributions in tibble
lg_dist <- bind_rows(dp_lg_results)
# 1.3.5. Save as csv
write_csv(lg_dist, "lg_dist.csv", col_names = T)

# 1.4. Load csv
lg_dist <- read_csv("lg_dist.csv", col_names = T)

# 2. Compare read length distributions 
## 2.1. Weighted Kolmogorov-Smirnov test ----
#Nonparametric
#Tests the entire distribution shape
#Appropriate when you do not want to assume normality

### 2.1.1. Expand counts into individual read lengths
expanded <- lg_dist %>%
  uncount(Occurrences)
### 2.1.2. Run KS test per SampleID
ks_results <- expanded %>%
  group_by(SampleID) %>%
  group_modify(~ {
    libs <- unique(.x$LibraryID)
    if (length(libs) != 2) return(tibble())
    
    x <- .x %>% filter(LibraryID == libs[1]) %>% pull(Length)
    y <- .x %>% filter(LibraryID == libs[2]) %>% pull(Length)
    
    ks <- ks.test(x, y)
    
    tibble(
      Library1 = libs[1],
      Library2 = libs[2],
      statistic = unname(ks$statistic),
      p_value   = ks$p.value,
      ks_raw    = list(ks)   # <- full htest object
    )
  }) %>%
  ungroup()
write_csv(ks_results, "ks_test_libs.csv", col_names = T)
view(ks_results) # Read length distributions are significantly different
ks_results$ks_raw[[1]]
ks_results$ks_raw[[2]]

### 2.1.2. Downsample reads to equal size (conservative) and repeat KS test
set.seed(1)

balanced <- expanded %>%
  group_by(SampleID) %>%
  group_modify(~ {
    min_reads <- .x %>%
      group_by(LibraryID) %>%
      summarise(n_reads = n(), .groups = "drop") %>%
      summarise(min_n = min(n_reads)) %>%
      pull(min_n)
    
    # Downsample each library to min_reads
    .x %>%
      group_by(LibraryID) %>%
      slice_sample(n = min_reads) %>%
      ungroup()
  }) %>%
  ungroup()

ks_results_bal <- balanced %>%
  group_by(SampleID) %>%
  group_modify(~ {
    libs <- unique(.x$LibraryID)
    if (length(libs) != 2) return(tibble())
    
    x <- .x %>% filter(LibraryID == libs[1]) %>% pull(Length)
    y <- .x %>% filter(LibraryID == libs[2]) %>% pull(Length)
    
    ks <- ks.test(x, y)
    
    tibble(
      Library1 = libs[1],
      Library2 = libs[2],
      statistic = unname(ks$statistic),
      p_value   = ks$p.value,
      ks_raw    = list(ks)   # <- full htest object
    )
  }) %>%
  ungroup()
write_csv(ks_results_bal, "ks_test_libs_equalN.csv", col_names = T)
view(ks_results_bal) # Read length distributions are still significantly different
ks_results_bal$ks_raw[[1]]
ks_results_bal$ks_raw[[2]]

## 2.2. Weighted Wilcoxon (median shift) ----
#Tests difference in central tendency (not full distribution)
#wilcox.test(Length ~ LibraryID,
#            data = lg_dist %>% filter(SampleID == "RxA007"),
#            weights = "Occurrences") # Medians not significantly different
#wilcox.test(Length ~ LibraryID,
#            data = lg_dist %>% filter(SampleID == "ScA005"),
#            weights = "Occurrences") # Medians not significantly different
## But note that this is a less powerful test than the KS test


# 3. Plot read length distributions ----
# Prepare annotation dataframe for KS results
ks_annot <- ks_results %>%
  mutate(label = paste0("D = ", round(statistic, 3), 
                   ", p = ", signif(p_value, 3))) %>%
  select(SampleID, label)
# Prepare plot labels
## 3.1. Barplot - absoulte counts
lg_dist %>%
  filter(Length <= 100) %>%
  ggplot() +
  geom_col(aes(Length, Occurrences, fill = LibraryID), position = "dodge") +
  facet_wrap(~SampleID, scales = "free", ncol = 1, nrow = 2)
## 3.2. Line plot - absolute counts
lg_dist %>%
  filter(Length <= 100) %>%
  ggplot(aes(Length, Occurrences, colour = LibraryID, group = LibraryID)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(breaks = c("SCR_002", "BEST_012", "SCR_007", "BEST_017"),
                      values = c("#1E5EA8", "#ADD8E6", "#2E8B28", "#C7E89A")) +
  facet_wrap(~SampleID, scales = "free_y", ncol = 1) +
  geom_text(data = ks_annot, aes(x = 80, y = Inf, label = label), 
            inherit.aes = FALSE, vjust = 1.2, hjust = 1, size = 4.5) +
  labs(y = "Read count", x = "Length (bp)") +
  theme_classic(base_size = 16) +
  theme(strip.background = element_rect(fill = "white", colour = NA)) 
ggsave("plot_readLength_count_lines.png", dpi = 300)
ggsave("plot_readLength_count_lines.pdf", width = 6.11, height = 4.19, units = "in")

## 3.3. Line plot - normalized frequencies - Figure 3
lg_dist %>%
  filter(Length <= 100) %>%
  group_by(SampleID, LibraryID) %>%
  mutate(freq = Occurrences / sum(Occurrences)) %>%
  mutate(letter = case_when(SampleID == "RxA007" ~ "A", SampleID == "ScA005" ~ "B")) %>%
  ggplot(aes(Length, freq, colour = LibraryID, group = LibraryID)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(breaks = c("SCR_002", "BEST_012", "SCR_007", "BEST_017"),
                      values = c("#1E5EA8", "#ADD8E6", "#2E8B28", "#C7E89A")) +
  facet_wrap(~SampleID, scales = "free_y", ncol = 1) +
  geom_text(data = ks_annot, aes(x = 80, y = Inf, label = label), 
            inherit.aes = FALSE, vjust = 1.2, hjust = 1, size = 4.5) +
  labs(y = "Proportion of reads", x = "Length (bp)") +
  theme_classic(base_size = 16) +
  theme(strip.background = element_rect(fill = "white", colour = NA))
ggsave("plot_readLength_freq_lines.png", dpi = 300)
ggsave("plot_readLength_freq_lines.pdf", width = 6.11, height = 4.19, units = "in") # Add letters in inkscape


## 3.4. Empirical CDF plots (directly representing the KS test)
#(empirical cumulative distribution function)
lg_dist %>%
  filter(Length <= 150) %>%
  group_by(SampleID, LibraryID) %>%
  summarise(
    Length = Length,
    Occurrences = Occurrences,
    .groups = "drop"
  ) %>%
  uncount(Occurrences) %>%
  ggplot(aes(Length, colour = LibraryID)) +
  stat_ecdf(size = 1) +
  facet_wrap(~SampleID, scales = "free_y", ncol = 1) +
  labs(y = "Empirical CDF")


