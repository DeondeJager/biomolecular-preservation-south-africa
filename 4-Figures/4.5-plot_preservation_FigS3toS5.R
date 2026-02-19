# Plot preservation of DNA and collagen site, taxon, material - Makes Fig S3-S5 (additional edits done in Inkscape)

# Set working directory as necessary
#setwd("")

## Load packages
library(tidyverse)
library(jsonlite)
library(rstatix)
library(lme4)
library(MASS)
library(scales)
library(ggpubr)

# 1. Read data for plotting ----
Data <- read_csv("Table S3 - paleomix nuclear & collagen.csv")

## Remove negative controls (not mapped to nuclear genomes)
Data <- Data %>% filter(Taxon_plots_common != "Negative control")

## Remove individual libraries for samples with multiple libs, keeping only the "combined" sample
Data <- Data %>% filter(!LibraryID %in% c("BEST_013", "BEST_014", "BEST_015", "SCR_003", "SCR_004", "SCR_005",
                                          "BEST_011", "SCR_001", "BEST_012", "SCR_002", "BEST_017", "SCR_007",
                                          "BEST_018", "SCR_008", "BEST_019", "SCR_009"))

## Pivot data into long format for ggplot
Data_long <- Data %>%
  pivot_longer(cols = seq_reads_pairs:last_col(),
               names_to = "variable", values_to = "value") %>%
  separate(col = "variable", into = c("variable", "reference"), sep = "x") # Split variable column on "x" to get a column with reference species names (i.e. "Target" and "Human")


# SI Figures ----
## Fig. S3: Endo DNA and collagen by site
### Define sites to plot for collagen plots (Phase 1), so it also includes Byneskranskop 1 (n = 0)
all_sites <- c("Boomplaas Cave","Klasies River Mouth",
               "Nelson Bay Cave","Byneskranskop 1",
               "Die Kelders Cave 1","Elands Bay Cave")
### Endo DNA
p_S3A <- Data_long %>% 
  filter(reference == "Target") %>%
  filter(variable == "noMicro_endoProp") %>%
  #  filter(value != "NA") %>% count(species_identifiable) # Yes = 65/144 = 45%, No = 79/144 = 55%
  #  filter(value != "NA") %>% filter(latest_chrono_age < 11700) %>% count(species_identifiable) # Holocene: Yes = 61/97 = 63%, No = 36/97 = 37%
  #filter(value != "NA") %>% filter(latest_chrono_age >= 11700) %>% count(species_identifiable) # Pleistocene: Yes = 4/47 = 8.5%, No = 91.5%
  ggplot() +
  geom_point(aes(latest_chrono_age, value*100, colour = species_identifiable), shape = 17, size = 6, alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_log10(labels = label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 10),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_continuous("Endogenous DNA (%)") +
  coord_cartesian(ylim = c(0, 30.5), clip = "off") +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  #annotate(geom = "text", x = 11700, y = 34, label = "11.7", size = 4) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14)) +
  facet_wrap(~Site_long)
p_S3A

### Collagen
p_S3B <- Data_long %>% 
  filter(Site_long %in% all_sites) %>%
  filter(StudyPhase == "Phase 1") %>%
  filter(variable == "Coll_Yield_perc") %>%
  complete(Site_long = all_sites, Taxon_plots_common) %>%  # fill missing combinations
  ggplot() +
  geom_point(aes(latest_chrono_age, value, colour = species_identifiable), shape = 17, size = 6, alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_log10(labels = label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 10),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_continuous("Collagen yield (%)") +
  coord_cartesian(ylim = c(0, 8.5), clip = "off") +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  #annotate(geom = "text", x = 11700, y = 34, label = "11.7", size = 4) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14)) +
  facet_wrap(~Site_long)
p_S3B

### Combine
p_S3 <- cowplot::plot_grid(p_S3A, p_S3B, nrow = 2, ncol = 1, labels = "AUTO")
p_S3
ggsave("FigureS3.png", p_S3, width = 200, height = 290, units = "mm", dpi = 300)
ggsave("FigureS3.pdf", p_S3, width = 200, height = 290, units = "mm")

## Fig. S4: Endo DNA and collagen by taxon
### Endo DNA
p_S4A <- Data_long %>% 
  filter(reference == "Target") %>%
  filter(variable == "noMicro_endoProp") %>%
  #  filter(value != "NA") %>% count(species_identifiable) # Yes = 65/144 = 45%, No = 79/144 = 55%
  #  filter(value != "NA") %>% filter(latest_chrono_age < 11700) %>% count(species_identifiable) # Holocene: Yes = 61/97 = 63%, No = 36/97 = 37%
  #filter(value != "NA") %>% filter(latest_chrono_age >= 11700) %>% count(species_identifiable) # Pleistocene: Yes = 4/47 = 8.5%, No = 91.5%
  ggplot() +
  geom_point(aes(latest_chrono_age, value*100, colour = species_identifiable), shape = 17, size = 6, alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_log10(labels = label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 10),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_continuous("Endogenous DNA (%)") +
  coord_cartesian(ylim = c(0, 30.5), clip = "off") +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  #annotate(geom = "text", x = 11700, y = 34, label = "11.7", size = 4) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14)) +
  facet_wrap(~Taxon_plots_common)
p_S4A

### Collagen
p_S4B <- Data_long %>% 
  filter(StudyPhase == "Phase 1") %>%
  filter(variable == "Coll_Yield_perc") %>%
  ggplot() +
  geom_point(aes(latest_chrono_age, value, colour = species_identifiable), shape = 17, size = 6, alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_log10(labels = label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 10),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_continuous("Collagen yield (%)") +
  coord_cartesian(ylim = c(0, 8.5), clip = "off") +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  #annotate(geom = "text", x = 11700, y = 34, label = "11.7", size = 4) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14)) +
  facet_wrap(~Taxon_plots_common)
p_S4B

### Combine
p_S4 <- cowplot::plot_grid(p_S4A, p_S4B, nrow = 2, ncol = 1, labels = "AUTO")
p_S4
ggsave("FigureS4.png", p_S4, width = 200, height = 290, units = "mm", dpi = 300)
ggsave("FigureS4.pdf", p_S4, width = 200, height = 290, units = "mm")


## Fig. S5: Endo DNA and collagen by material group
### Endo DNA
p_S5A <- Data_long %>% 
  filter(reference == "Target") %>%
  filter(variable == "noMicro_endoProp") %>%
  #  filter(value != "NA") %>% count(species_identifiable) # Yes = 65/144 = 45%, No = 79/144 = 55%
  #  filter(value != "NA") %>% filter(latest_chrono_age < 11700) %>% count(species_identifiable) # Holocene: Yes = 61/97 = 63%, No = 36/97 = 37%
  #filter(value != "NA") %>% filter(latest_chrono_age >= 11700) %>% count(species_identifiable) # Pleistocene: Yes = 4/47 = 8.5%, No = 91.5%
  ggplot() +
  geom_point(aes(latest_chrono_age, value*100, colour = species_identifiable), shape = 17, size = 6, alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_log10(labels = label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 10),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_continuous("Endogenous DNA (%)") +
  coord_cartesian(ylim = c(0, 30.5), clip = "off") +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  #annotate(geom = "text", x = 11700, y = 34, label = "11.7", size = 4) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14)) +
  facet_wrap(~material_group)
p_S5A

### Collagen
p_S5B <- Data_long %>% 
  filter(StudyPhase == "Phase 1") %>%
  filter(variable == "Coll_Yield_perc") %>%
  #  filter(value != "NA") %>% count(species_identifiable) # Yes = 65/144 = 45%, No = 79/144 = 55%
  #  filter(value != "NA") %>% filter(latest_chrono_age < 11700) %>% count(species_identifiable) # Holocene: Yes = 61/97 = 63%, No = 36/97 = 37%
  #filter(value != "NA") %>% filter(latest_chrono_age >= 11700) %>% count(species_identifiable) # Pleistocene: Yes = 4/47 = 8.5%, No = 91.5%
  ggplot() +
  geom_point(aes(latest_chrono_age, value, colour = species_identifiable), shape = 17, size = 6, alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_log10(labels = label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 10),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_continuous("Collagen yield (%)") +
  coord_cartesian(ylim = c(0, 8.5), clip = "off") +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  #annotate(geom = "text", x = 11700, y = 34, label = "11.7", size = 4) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14)) +
  facet_wrap(~material_group)
p_S5B

### Combine
p_S5 <- cowplot::plot_grid(p_S5A, p_S5B, nrow = 2, ncol = 1, labels = "AUTO")
p_S5
ggsave("FigureS5.png", p_S5, width = 200, height = 290, units = "mm", dpi = 300)
ggsave("FigureS5.pdf", p_S5, width = 200, height = 290, units = "mm")
