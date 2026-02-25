# Plot screening data mapped to mito_panel for species ID - Creates Figure S1

# Set working directory
#setwd("")

## Load packages
library(tidyverse)
library(reshape2)
library(writexl)
library(ggforce)
library(ggh4x) # to duplicate discrete axis in ggplot (currently only natively available for continuous scales)
library(grid)
library(cowplot)
library(ggpubr)

# Read data
screening <- read_csv("Table S2 - paleomix mito panel.csv", col_names = T)

### 3.2. Plot by sample
#### 3.2.1. Convert SampleID to factors so we can explicitly maintain the order
screening$SampleID <- factor(screening$SampleID, levels = unique(screening$SampleID))
#### 3.2.2. Make an empty list to store each sample's plots in
results_list <- list()
#### 3.2.2. Plot
for (sample_id in levels(screening$SampleID)){
  # Get LibraryID of the current sample_id
  library_id <- screening %>%
    filter(SampleID == sample_id) %>%
    pull(LibraryID) %>%
    unique()
  ## Safety check
  stopifnot(length(library_id) == 1)
  # Plot
  p1 <- ggplot(filter(screening, SampleID == sample_id) %>% filter(variable == "hits_unique")) +
    geom_col(aes(x = reference, y = value, fill = species_identifiable)) +
    scale_fill_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "darkgrey"),
                      name = "Species identifiable?") +
    geom_text(aes(x = reference, y = value, label = value), 
              size = 3.5, vjust = -0.2) +
    #labs(x = "Reference mitogenome", y = "No. of unique reads mapped") +
    #theme(axis.title = element_text(size = 16),
    #      axis.text = element_text(size = 14),
    #      strip.text = element_text(size = 16)) +
    #facet_wrap_paginate(~SampleID + LibraryID, scales = "free_y", ncol = 1, nrow = 1, page = i) +
    scale_x_discrete(guide = guide_axis(angle = 45), name = NULL) +
    scale_y_continuous(name = NULL) +
    theme_classic(base_size = 14) +
    labs(subtitle = paste0(sample_id, ":", " ", library_id)) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 10),
          plot.subtitle = element_text(size = 12, face = "bold"))
  ## Get max y-value
  #max_y <- screening %>%
  #  filter(SampleID == sample_id) %>%
  #  filter(variable == "hits_unique") %>%
  #  pull("value") %>%
  #  max()
  ## Safety check
  #stopifnot(length(max_y) == 1)
  ## Add label
  #p2 <- ggdraw(p1) + 
  #  draw_label(paste0(sample_id,"\n",library_id), 
  #             fontface = "bold", 
  #             size = 10, 
  #             x = 10, y = max_y*0.9)
  # Add all the plots to the results list
  results_list[[length(results_list) + 1]] <- p1
}
results_list[[1]]
#results_list # plots all subplots in plotting area, so can take a while

### 3.3 Combine plots into groups of 10 in 2x5 grids
plots_per_page <- 10
results_chunks <- split(results_list,
                        ceiling(seq_along(results_list) / plots_per_page))
p1_list <- map(results_chunks, 
               ~ cowplot::plot_grid(plotlist = .x,
                                    ncol = 2,
                                    nrow = 5,
                                    align = "hv"))
p1_list[[1]]

### 3.4. Get legend to use as common legend - extract it from another plot
p_legend <- screening  %>%
  filter(SampleID %in% c("PaA001", "PaA014")) %>%
  filter(variable == "hits_unique") %>%
  ggplot() +
  geom_col(aes(x = reference, y = value, fill = species_identifiable)) +
  scale_fill_manual(breaks = c("Yes", "No"), 
                    values = c("#2A9E72", "darkgrey"),
                    name = "Species identifiable?") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")
p_legend
leg <- ggpubr::get_legend(p_legend)
leg

### 3.5. Arrange all plots in a grid with legend
p1_list <- map(p1_list,
               ~ ggarrange(.x, ncol = 1, nrow = 1, labels = NULL, align = "hv",
                           common.legend = TRUE, legend.grob = leg, legend = "top"))
p1_list[[1]]

### 3.6. Add x and y titles (common for all plots)
p1_list <- map(p1_list, 
               ~ annotate_figure(.x, left = text_grob("Unique reads mapped", size = 16, rot = 90), 
                                 bottom = text_grob("Reference mitogenome", size = 16)))
p1_list[[7]] # Check it worked


### 3.7. Save plots
#### As one PDF with 10 plots per page
pdf("FigureS1.pdf", width = 8.2, height = 11.496)
#print(p1_list[[1]]) # Test with one set of 10 
walk(p1_list, print) # Plot all
dev.off()

#### As sets of PNGs
#walk2(p1_list, seq_along(p1_list),
#      ~ ggsave(filename = sprintf("mito_panel_plots_pages_png/mito_panel_plots_page_%02d.png", .y),
#               plot = .x, width = 205, height = 292, units = "mm", dpi = 300))
