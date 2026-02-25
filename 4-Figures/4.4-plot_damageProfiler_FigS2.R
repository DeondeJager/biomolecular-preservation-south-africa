# Get and plot damage patterns and read lengths from DamageProfiler JSON output

# Set working directory ----
#setwd("")

# Load packages ----
library(tidyverse)
library(jsonlite)
library(ggpubr)
library(cowplot)

# Get misincorporation data ----
# Get damage profiles - Only need to do steps 1-7 once, after which just need 
# to load dp_misincorporation_summary_R.csv (step 7).
# 1. Get list of all misincorporation files to process
mis_files <- list.files(path = "damageProfiler",
                        pattern = "misincorporation\\.txt$",
                        recursive = TRUE,
                        full.names = TRUE)
## 1.1. Build a tibble and extract the folder name
mis_tbl <- tibble(path = mis_files) %>%
  mutate(folder = basename(dirname(path)))
# 1.2. Extract SeqID (everything before the period)
mis_tbl <- mis_tbl %>%
  mutate(SeqID = str_remove(folder, "\\.damageProfiler$"))
## 1.3. Split into SampleID and LibraryID
mis_tbl <- mis_tbl %>%
  mutate(SampleID  = str_extract(SeqID, "^[^_]+"),
         LibraryID = str_remove(SeqID, "^[^_]+_")) %>%
  mutate(LibraryID = if_else(LibraryID == SeqID, NA_character_, LibraryID))
## 1.4. Save as csv and manually add missing LibraryID, then load back in
write_csv(mis_tbl, "mis_tbl.csv")
mis_tbl <- read_csv("mis_tbl.csv", col_names = T)

# 2. Add "species_identifiable" column to filter on - only plotting libraries with true endo aDNA
gt_data <- read_csv("Table S3 - paleomix nuclear & collagen.csv", 
                    col_select = c(SampleID, species_identifiable))
## 2.1. Join with misincorporation table
mis_tbl <- mis_tbl %>%
  left_join(gt_data, by = "SampleID")
# 2.2. Filter to include only libraries where the species was genetically identifiable (true endo DNA)
mis_tbl_yes <- mis_tbl %>%
  filter(species_identifiable == "Yes")
## 2.3. Update mis_files to include only those with true endo DNA
mis_files_yes <- as.vector(mis_tbl_yes$path)

# 3. Initialize an empty list to store each file's processed data
results_list <- list()

# 4. Loop through each file
for (file in mis_files_yes) {
  
  folder <- basename(dirname(normalizePath(file)))
  SeqID  <- str_remove(folder, "\\.damageProfiler$")
  
  # Find the line number where "Chr" first appears in *misincorporation.txt file
  first_chr_line <- grep("^Chr", readLines(file), ignore.case = FALSE)[1]
  # Read the TSV starting from that line
  # subtract 1 because skip counts lines *before* reading
  data <- read_tsv(file, skip = first_chr_line - 1, col_names = TRUE)
  
  # Perform calculations on each row
  processed_data <- data %>%
    #rowwise() %>%
    mutate(GtoA = `G>A`/G) %>% 
    mutate(CtoT = `C>T`/C) %>%
    #ungroup() %>%
    select(End, Std, Pos, GtoA, CtoT) %>%
    mutate(SeqID = SeqID) # Add a columns with the SeqID
  
  # Add the processed data to the results list
  results_list[[length(results_list) + 1]] <- processed_data
}

# 5. Combine all the processed data into a single tibble by row-binding
dp_mis <- bind_rows(results_list)
# View the final combined tibble
print(dp_mis)

# 6. Average the frequencies per position across the two strands (+ and -)
tmp1 <- dp_mis %>% 
  group_by(End, Pos, SeqID) %>% 
  summarise(GtoA_mean = mean(GtoA)) %>%
  arrange(SeqID)
tmp2 <- dp_mis %>% 
  group_by(End, Pos, SeqID) %>% 
  summarise(CtoT_mean = mean(CtoT)) %>%
  arrange(SeqID)
# Overwrite previous dataframe with new one with means across positions
dp_mis <- inner_join(tmp1, tmp2)
print(dp_mis)
# Combine with mis_tbl_yes on SeqID column, to get SampleID and LibraryID
dp_mis <- dp_mis %>%
  left_join(mis_tbl_yes, by = "SeqID")

# 7. Save as file and load again
write_csv(dp_mis, "dp_misincorporation_summary_R.csv")
mis_yes <- read_csv("dp_misincorporation_summary_R.csv", col_names = T)

# Get read length data ----
# Only need to do steps 8-13 once, after that just load dp_length_dist_summary_r.csv (step 12)
# and dp_length_summary_r.csv (step 13)
# 8. Get list of all json files with length distribution info to process
lg_files <- list.files(path = "damageProfiler",
                        pattern = "dmgprof\\.json$",
                        recursive = TRUE,
                        full.names = TRUE)
## 8.1. Build a tibble and extract the folder name
lg_tbl <- tibble(path = lg_files) %>%
  mutate(folder = basename(dirname(path)))
## 8.2 Combine with mis_tbl (load from step 1.4. if necessary) to get SeqID, SampleID and LibraryID
lg_tbl <- lg_tbl %>%
  left_join(mis_tbl %>% select(-path), by = "folder")
### Save as csv and load again
write_csv(lg_tbl, "lg_tbl.csv")
lg_tbl <- read_csv("lg_tbl.csv", col_names = T)

# 9. Same as step 2 for misincorporation data: Add "species_identifiable" column to filter on; 
# only plotting libraries with true endo aDNA
gt_data <- read_csv("Table S3 - paleomix nuclear & collagen.csv", 
                    col_select = c(SampleID, species_identifiable))
## 9.1. Join with length table
lg_tbl <- lg_tbl %>%
  left_join(gt_data, by = "SampleID")
## 9.2. Filter to include only libraries where the species was genetically identifiable (true endo DNA)
lg_tbl_yes <- lg_tbl %>%
  filter(species_identifiable == "Yes")
## 9.3. Make a vector of length files to include only those with true endo DNA
lg_files_yes <- as.vector(lg_tbl_yes$path)

# 10. Process json function
process_dp_json <- function(json_files) {
  
  dp <- fromJSON(json_files)
  
  folder <- basename(dirname(normalizePath(json_files)))
  SeqID  <- str_remove(folder, "\\.damageProfiler$")
  
  # Convert forward and reverse length distributions
  lendist_fw <- enframe(dp$lendist_fw, name = "Length", value = "Occurrences")
  lendist_rv <- enframe(dp$lendist_rv, name = "Length", value = "Occurrences")
  
  lendist <- bind_rows(lendist_fw, lendist_rv) %>%
    mutate(
      Length = as.numeric(Length),
      Occurrences = as.numeric(unlist(Occurrences))
    ) %>%
    group_by(Length) %>%
    summarise(Occurrences = sum(Occurrences), .groups = "drop") %>%
    mutate(SeqID = SeqID)
  
  median_length <- dp$summary_stats$median
  mean_length <- dp$summary_stats$mean
  
  list(
    lendist = lendist,
    median = median_length,
    mean = mean_length
  )
}


# 11. Process json files for read length dist (only those wit true endo DNA)
dp_lg_results <- map(lg_files_yes, process_dp_json)

# 12. Combine all length distributions
lg_dist_yes <- map_dfr(dp_lg_results, "lendist")
# 13. Combine with lg_tbl_yes to get SampleID and LibraryID
lg_dist_yes <- lg_dist_yes %>%
  left_join(lg_tbl_yes, by = "SeqID")
## Save as csv and load again
write_csv(lg_dist_yes, "dp_length_dist_summary_r.csv")
lg_dist_yes <- read_csv("dp_length_dist_summary_r.csv", col_names = T)

# 13. Get library-level (SeqID) summary stats of read lengths
sample_lg_stats <- map_dfr(dp_lg_results, 
                        ~ tibble(SeqID = unique(.x$lendist$SeqID),
                                                median_lg = .x$median,
                                                mean_lg = .x$mean))
## Combine with lg_tbl_yes to get SampleID, LibraryID, etc.
sample_lg_stats <- sample_lg_stats %>%
  left_join(lg_tbl_yes, by = "SeqID")
## Save as csv and load again
write_csv(sample_lg_stats, "dp_length_summary_r.csv")
lg_stats_yes <- read_csv("dp_length_summary_r.csv", col_names = T)

# Plot damage patterns and read length dist together ----
## Make individual plots and then add them together with cowplot
# 14. Convert to factors so we can explicitly maintain the order that LibraryID is currently in
mis_yes$LibraryID <- factor(mis_yes$LibraryID, levels = unique(mis_yes$LibraryID))
## Get max value to set y-axis limits
filter(mis_yes, Pos <= 25) %>% ungroup() %>% summarise(max(CtoT_mean, na.rm = T)) # 0.405

# 15. Make an empty list to store each library's plots in
results_list <- list()
# 16. Loop through each SeqID in the dataframe and make plots - note we use "levels" instead of "unique" in the for loop, as unique() does not respect the order of levels
for (lib_name in levels(mis_yes$LibraryID)) {
  # Make plots
  ## 5p end (CtoT)
  p1 <- ggplot(filter(mis_yes, Pos <= 25) %>% filter(LibraryID == lib_name) %>% filter(End == "5p")) +
    theme_pubr() +
    geom_line(aes(Pos, CtoT_mean), colour = "red", linewidth = 1) +
    geom_line(aes(Pos, GtoA_mean), colour = "blue", linewidth = 1) +
    scale_y_continuous(limits = c(0,0.41), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4), name = NULL) +
    scale_x_continuous(breaks = seq(from = 0, to = 25, by = 5), name = NULL)
  ## 3p end (GtoA)
  p2 <- ggplot(filter(mis_yes, Pos <= 25) %>% filter(LibraryID == lib_name) %>% filter(End == "3p")) +
    theme_pubr() +
    geom_line(aes(Pos, CtoT_mean), colour = "red", linewidth = 1) +
    geom_line(aes(Pos, GtoA_mean), colour = "blue", linewidth = 1) +
    scale_x_reverse(breaks = seq(from = 0, to = 25, by = 5), labels = c("0", "-5", "-10", "-15", "-20", "-25"), name = NULL) +
    scale_y_continuous(limits = c(0,0.41), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4), name = NULL, position = "right", labels = NULL)
  ## Join together in a side-by-side plot with joint title
  p3 <- cowplot::plot_grid(p1, p2, ncol = 2, align = "hv")
  ## Get SampleID of the current lib_name
  sample_id <- mis_yes %>%
    filter(LibraryID == lib_name) %>%
    pull(SampleID) %>%
    unique()
  ## Safety check
  stopifnot(length(sample_id) == 1)
  ## Add label
  p4 <- ggdraw(p3) + 
    draw_label(paste0(sample_id,"\n",lib_name), fontface = "bold", size = 14, x = 0.5, y = 0.9)
  ## Plot length
  ### Get length stats first
  lg_stat <- lg_stats_yes %>%
    filter(LibraryID == lib_name)
  stopifnot(nrow(lg_stat) == 1)
  ### Then plot
  p_lg <- ggplot(filter(lg_dist_yes, LibraryID == lib_name)) +
    theme_pubr() +
    geom_col(aes(Length, Occurrences)) +
    ## Mean read length (vertical dashed line)
    geom_vline(data = lg_stat, aes(xintercept = mean_lg),
               linetype = "dashed") +
    ## Mean label
    geom_text(data = lg_stat,
              aes(x = mean_lg, y = Inf, label = paste0("Mean=", round(mean_lg, 1))),
              vjust = 1.2, hjust = -0.1, size = 3.5) +
    scale_x_continuous(breaks = seq(30, 150, by = 10),
                       labels = c("30", "", "", "60", "", "", "90", "", "", "120", "", "", "150"),
                       name = "Read length (bp)") +
    coord_cartesian(xlim = c(30, 150)) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                       n.breaks = 3,
                       name = "Count")
  # Join misincorporation and read length plots
  p4_lg <- ggdraw() +
    draw_plot(p4) +
    draw_plot(p_lg, width = 0.5, height = 0.5, hjust = -0.5, vjust = -0.5)
  
  # Add all the plots to the results list
  results_list[[length(results_list) + 1]] <- p4_lg
}
results_list[[1]]
results_list

# 17. Combine plots into groups of 10 in 2x5 grids
plots_per_page <- 10
results_chunks <- split(results_list,
                        ceiling(seq_along(results_list) / plots_per_page))
p5_list <- map(results_chunks, 
               ~ cowplot::plot_grid(plotlist = .x,
                                    ncol = 2,
                                    nrow = 5,
                                    align = "hv"))
# p5_list[[1]] = plots 1-10
# p5_list[[2]] = plots 11-20, etc.

# 18. Get legend to use as common legend - extract it from another plot
## Transform dataframe into long format for End (5p, 3p)
mis_yes_long <- pivot_longer(mis_yes, cols = c("GtoA_mean", "CtoT_mean"), names_to = "Damage", values_to = "Frequency")
mis_yes_long$End <- factor(mis_yes_long$End, levels = c("5p", "3p"))
p_legend <- ggplot() +
  theme_pubr() +
  geom_line(data = filter(mis_yes_long, Pos <= 10) %>% filter(SeqID == "PaA014"), 
            mapping = aes(Pos, Frequency, colour = Damage), linewidth = 1) +
  scale_colour_manual(values = c("red", "blue"), labels = c("T", "A"), name = "Nucleotide:") +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        legend.margin = margin(c(0,0,0,0))) +
  facet_wrap(~End)
p_legend
leg <- ggpubr::get_legend(p_legend)
leg

# 19. Arrange all plots in a grid with legend
## Test
#p5_list[[1]] <- ggarrange(p5_list[[1]], ncol = 1, nrow = 1, labels = NULL, align = "hv",
#                          common.legend = TRUE, legend.grob = leg, legend = "top")
#p5_list[[1]]

## Apply (re-run step 17 before doing this, if you ran the test with p5_list[[1]] above) 
p5_list <- map(p5_list,
               ~ ggarrange(.x, ncol = 1, nrow = 1, labels = NULL, align = "hv",
                           common.legend = TRUE, legend.grob = leg, legend = "top"))
p5_list[[1]] # Check it worked

# 20. Add x and y titles (common for all plots)
p5_list <- map(p5_list, 
               ~ annotate_figure(.x, left = text_grob("Frequency", size = 16, rot = 90), 
                                 bottom = text_grob("Read position", size = 16)))
p5_list[[7]] # Check it worked

# 21. Save final plots

## As pdf, 10 plots per page
pdf("damage_plots_10perpage.pdf", width = 8.07, height = 11.496)
walk(p5_list, print)
dev.off()

## As sets of PNGs
walk2(p5_list, seq_along(p5_list),
      ~ ggsave(filename = sprintf("damage_plots_pages_png/damage_plots_page_%02d.png", .y),
               plot = .x, width = 205, height = 292, units = "mm", dpi = 300))
