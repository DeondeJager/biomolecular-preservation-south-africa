# Figure 1 - Sampling overview

# Set working directory
#setwd("")

# Load packages
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)

# Load data
samples <- read_csv("Table S1 - Specimen metadata.csv", col_names = T, 
                    name_repair = "universal", # Remove spaces from column names and replaces with "."
                    n_max = 321) # Prevents reading of additional info after table

## Rename name.of.site.or.location.where.sample.originated column to Site_long
samples <- rename(samples, Site_long = name.of.site.or.location.where.sample.originated)

# Plot samples ----
all_sites <- c("Boomplaas Cave","Klasies River Mouth","Nelson Bay Cave","Byneskranskop 1","Die Kelders Cave 1","Elands Bay Cave")
site_labels <- c("Boomplaas\nCave","Klasies River\nMouth","Nelson Bay\nCave","Byneskranskop\n1","Die Kelders\nCave 1","Elands Bay\nCave")

## Check age ranges to set x-axis limits (do for Phase 1 and 2 separately)
#samples %>%
#  filter(StudyPhase == "Phase 2") %>%
#  pull(latest.chronometric.age) %>%
#  range()

## Plot Phase 1
p_phase1 <- samples %>%
  filter(Site_long %in% all_sites) %>%
  filter(StudyPhase == "Phase 1") %>%
  complete(Site_long = all_sites, Taxon_plots_common) %>%  # fill missing combinations
  ggplot() +
  geom_point(aes(latest.chronometric.age, y = Taxon_plots_common, colour = Taxon_plots_common), 
             shape = 2, position = "jitter", size = 2.5) +
  facet_wrap(~Site_long, nrow = 6, ncol = 1, strip.position = "right", axes = "margins",
             labeller = labeller(Site_long = site_labels)) +
  scale_x_log10(limits = c(1400, 112000),
                labels = scales::label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 5),
                name = "Age (thousand years ago) [log10 scale]") +
  scale_y_discrete(name = "Taxon",
                   guide = guide_axis(angle = 45)) +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  coord_flip() +
  labs(title = "Phase 1") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
p_phase1  
ggsave("plot_phase1.pdf", p_phase1, width = 73, height = 290/2, units = "mm")

## Plot Phase 2
p_phase2 <- samples %>%
  filter(Site_long %in% all_sites) %>%
  filter(StudyPhase == "Phase 2") %>%
  complete(Site_long = all_sites, Taxon_plots_common) %>%  # fill missing combinations
  ggplot() +
  geom_point(aes(latest.chronometric.age, y = Taxon_plots_common, colour = Taxon_plots_common), 
             shape = 2, position = "jitter", size = 2.5) +
  facet_wrap(~Site_long, nrow = 6, ncol = 1, strip.position = "right", axes = "margins") +
  scale_x_log10(limits = c(1400, 112000),
                labels = scales::label_number(scale = 1e-3),
                breaks = scales::breaks_log(n = 5),
                name = "") +
  scale_y_discrete(name = "Taxon",
                   guide = guide_axis(angle = 45)) +
  geom_vline(xintercept = 11700, linetype = "dashed") +
  coord_flip() +
  labs(title = "Phase 2") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
p_phase2  
ggsave("plot_phase2.pdf", p_phase2, width = 73, height = 290/2, units = "mm")

# Combine plots
p_comb <- cowplot::plot_grid(p_phase1, p_phase2, nrow = 1, ncol = 2, 
            align = "v", 
            axis = "lr", 
            labels = c("A", "B"))
p_comb
ggsave("plot_phases.pdf", p_comb, width = 146, height = 292/2, units = "mm")



# Plot sites map ----
sites <- read_csv("Sites_GPS.csv", col_names = T)
map_sites <- c("BPA", "KRM", "NBC", "BNK1", "DK1", "EBC")
## Prep site data
sites_df <- samples %>%
  left_join(sites, by = "Site_short") %>%
  filter(Site_short %in% map_sites) %>%
  select(Site_short, Lat, Lon) %>%
  distinct()   # one row per site
## Convert to sf object
sites_sf <- st_as_sf(sites_df, coords = c("Lon", "Lat"), crs = 4326)
## Get the world map as an sf object
world <- ne_countries(scale = "medium", returnclass = "sf")
# Compute centroids for country labels
world_centroids <- st_centroid(world)

## Plot map with labels for use in later editing in Inkscape
#ggplot() +
#  geom_sf(data = world, fill = "gray95", color = "gray70") +
#  geom_sf(data = sites_sf, aes(), shape = 24, size = 4, fill = "#23232380", colour = "gray10") +
#  geom_sf_text(data = world_centroids, aes(label = name),              # country labels
#               size = 4, color = "gray30") +
#  coord_sf(xlim = c(min(sites_df$Lon) - 2, max(sites_df$Lon) + 2),
#           ylim = c(min(sites_df$Lat) - 2, max(sites_df$Lat) + 2)) + # zoom to region
#  theme_minimal(base_size = 14) +
#  labs(colour = "Site", x = NULL, y = NULL)
#ggsave("plot_Sites_map_labels.png")
# Determine zoom region
xlim_region <- c(min(sites_df$Lon) - 1, max(sites_df$Lon) + 2)
ylim_region <- c(min(sites_df$Lat) - 1, max(sites_df$Lat) + 2)

# Main map
main_map <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70") +
  geom_sf(data = sites_sf, aes(), shape = 24, size = 4, fill = "#23232380", colour = "gray10") +
  geom_sf_text(data = world_centroids, aes(label = name),
               size = 4, color = "gray30") +
  coord_sf(xlim = xlim_region, ylim = ylim_region) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14)
main_map
#ggsave("plot_Sites_map_labels.png", main_map)

# Africa map
africa <- world %>% filter(continent == "Africa")

# Bounding box for the zoomed region
bbox <- st_as_sfc(st_bbox(c(xmin = xlim_region[1], xmax = xlim_region[2],
                            ymin = ylim_region[1], ymax = ylim_region[2]),
                          crs = st_crs(africa)))

inset_map <- ggplot() +
  geom_sf(data = africa, fill = "gray90", color = "gray50") +
  geom_sf(data = bbox, fill = NA, color = "red", linewidth = 0.5) +   # highlight main map region
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.spacing = unit(0.0, "cm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  theme_void()
inset_map

# Place inset in the top-right corner of main map
ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.65, y = 0.6, width = 0.35, height = 0.35)
ggsave("plot_map.pdf", width = 6.11, height = 4.19, units = "in")
ggsave("plot_map.svg", width = 6.11, height = 4.19, units = "in")
