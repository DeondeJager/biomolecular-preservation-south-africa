# Plot preservation of DNA and collagen - Makes Fig2A-D (additional edits done in Inkscape)

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

# 2. Figure 2 ----
## 2.1. Fig2A - Endo vs age, species identifiable
p_2A <- Data_long %>% 
  filter(reference == "Target") %>%
  filter(variable == "noMicro_endoProp") %>%
  #  filter(value != "NA") %>% count(species_identifiable) # Yes = 65/144 = 45%, No = 79/144 = 55%
  #  filter(value != "NA") %>% filter(latest_chrono_age < 11700) %>% count(species_identifiable) # Holocene: Yes = 61/97 = 63%, No = 36/97 = 37%
  #filter(value != "NA") %>% filter(latest_chrono_age >= 11700) %>% count(species_identifiable) # Pleistocene: Yes = 4/47 = 8.5%, No = 91.5%
  ggplot() +
  geom_point(aes(latest_chrono_age, value*100, 
                 colour = species_identifiable), 
             shape = 17, size = 6,  alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_log10(labels = label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 10),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_continuous("Endogenous DNA (%)") +
  coord_cartesian(ylim = c(0, 30.5), clip = "off") +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  annotate(geom = "text", x = 40000, y = 28, label = "Yes = 4/47 (8.5%)", size = 4) + # Pleistocene
  annotate(geom = "text", x = 4000, y = 28, label = "Yes = 61/97 (63%)", size = 4) + # Holocene
  annotate(geom = "text", x = 4000, y = 34, label = "Holocene", size = 4) +
  annotate(geom = "text", x = 40000, y = 34, label = "Pleistocene", size = 4) +
  annotate(geom = "text", x = 11700, y = 34, label = "11.7k", size = 4) +
  #annotate("segment", x = 9500, xend = 7000, y = 0.37, yend = 0.37, size = 1, arrow=arrow(length = unit(0.05, "inches"))) +
  #annotate("segment", x = 14250, xend = 18750, y = 0.37, yend = 0.37, size = 1, arrow=arrow(length = unit(0.05, "inches"))) +
  #annotate(geom = "text", x = 40000, y = 0.05, label = "R^2 = 0.16, p-value = 1.74e-05\nDeviance explained = 17%") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14))
p_2A
### Get legend and then re-run plot without legend
#leg <- ggpubr::get_legend(p_2A)

## 2.2. Fig 2B - Collagen vs age, species identifiable
p_2B <- Data_long %>% 
  filter(StudyPhase == "Phase 1") %>%
  filter(variable == "Coll_Yield_perc") %>%
  #  filter(value != "NA") %>% count(species_identifiable) # Yes = 21/58 = 36%, No = 37/58 = 64%
  ggplot() +
  geom_point(aes(latest_chrono_age, value, colour = species_identifiable), shape = 17, size = 6, alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_log10(labels = label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 10),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_continuous("Collagen yield (%)") +
  coord_cartesian(ylim = c(0, 8.1), clip = "off") +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 14))
p_2B

## 2.3. Fig 2C - Endo DNA vs collagen yield
p_2C <- Data %>% 
  filter(StudyPhase == "Phase 1") %>%
  ggplot() +
  geom_point(aes(Coll_Yield_perc, noMicro_endoPropxTarget*100, colour = species_identifiable), shape = 17, size = 6, alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_continuous("Collagen yield (%)") +
  scale_y_continuous("Endogenous DNA (%)") +
  coord_cartesian(ylim = c(0, 30.5), clip = "off") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 14))
p_2C

## 2.4. Fig 2D - Mean read length vs age
p_2D <- Data_long %>% 
  filter(species_identifiable == "Yes") %>%
  filter(variable == "noMicro_meanLength") %>%
  ggplot() +
  geom_point(aes(latest_chrono_age, value, colour = species_identifiable), shape = 17, size = 6, alpha = 0.65) +
  scale_colour_manual(breaks = c("Yes", "No"), 
                      values = c("#2A9E72", "#A9A9A9"),
                      name = "Species identifiable?") +
  scale_x_log10(labels = label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 10),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_continuous("Mean read length (bp)") +
  coord_cartesian(ylim = c(0, 80), clip = "off") +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 14))
p_2D

## 2.5. Fig2 - Combine plots and legend (legend obtained in step 2.1)
p_2 <- ggpubr::ggarrange(p_2A, p_2C, p_2B, p_2D,
                         ncol = 2, nrow = 2,
                         common.legend = TRUE,
                         legend = "top",
                         legend.grob = leg,
                         align = "hv",
                         labels = "AUTO")
p_2

## 2.6. Save combined figures
ggsave("Figure2.png", p_2, width = 205, height = 292/2, units = "mm", dpi = 300)
ggsave("Figure2.pdf", p_2, width = 205, height = 292/2, units = "mm")

