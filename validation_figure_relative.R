# Alternative version of validation figure showing relative differences

library(tidyverse)
library(USAboundaries)
library(sf)
library(cowplot)
library(grid)

load(file = '~/temp/grandtotals_finalfig.RData')

label_names <- c('This study', 'U.S. EPA', 'Microsoft', 'Facebook', 'Block group area weighting')
estimate_labels <- data.frame(estimate = c('our_dasy', 'epa_dasy', 'huang', 'fb', 'blockgroup'),
                              estimate_label = factor(label_names, levels = label_names))

wildfire_risks_long <- wildfire_grandtotals %>%
  select(county_label, fips, env_class, our_dasy:blockgroup) %>%
  pivot_longer(our_dasy:blockgroup, names_to = 'estimate', values_to = 'population') %>%
  group_by(county_label, fips, estimate) %>%
  mutate(proportion = population/sum(population, na.rm = TRUE)) %>%
  left_join(estimate_labels)

flood_risks_long <- flood_grandtotals %>%
  select(county_label, fips, env_class, our_dasy:blockgroup) %>%
  pivot_longer(our_dasy:blockgroup, names_to = 'estimate', values_to = 'population') %>%
  group_by(county_label, fips, estimate) %>%
  mutate(proportion = population/sum(population, na.rm = TRUE)) %>%
  left_join(estimate_labels)

counties_use_wf <- c('49009', '35057', '49013', '53007', '16001')
wildfire_risks_reduced <- wildfire_risks_long %>%
  filter(env_class %in% as.character(3:5), fips %in% counties_use_wf) %>%
  group_by(county_label, estimate, estimate_label) %>%
  summarize(proportion = sum(proportion), population = sum(population))

counties_use_fl <- c('24019', '13029', '12035', '01003', '24003')
flood_risks_reduced <- flood_risks_long %>%
  filter(env_class %in% '1', fips %in% counties_use_fl) 


# Convert to relative difference vs naive
flood_risks_reduced <- flood_risks_reduced %>%
  group_by(county_label) %>%
  mutate(diff = population/population[estimate == 'blockgroup'] - 1) %>%
  filter(!estimate %in% 'blockgroup')

wildfire_risks_reduced <- wildfire_risks_reduced %>%
  group_by(county_label) %>%
  mutate(diff = population/population[estimate == 'blockgroup'] - 1) %>%
  filter(!estimate %in% 'blockgroup')

# Make plots --------------------------------------------------------------

okabe_colors <- palette.colors(7, palette = 'Okabe-Ito')[c(7, 3, 4, 2)]

zero_line <- geom_hline(yintercept = 0, linetype = 'dashed')

p_flood <- flood_risks_reduced %>%
  ungroup %>%
  mutate(county_label = factor(county_label, labels = c('Dorchester, Maryland\n(32,500)', 'Bryan, Georgia\n(34,100)', 'Flagler, Florida\n(103,000)', 'Baldwin, Alabama\n(200,000)', 'Anne Arundel,Maryland\n(560,000)')) %>% fct_rev) %>%
  ggplot(aes(x = county_label, y = diff, color = estimate_label)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_y_continuous(name = 'Change in flood risk estimate relative to naïve method',
                     limits = c(-1.04, 1.04), labels = scales::percent, position = 'right') + 
  zero_line +
  coord_flip() +
  scale_color_manual(name = 'Dasymetric method', values = unname(okabe_colors)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none',
        plot.margin = unit(c(1, 0.2, 0.1, 2), 'cm'))

p_wildfire <- wildfire_risks_reduced %>%
  ungroup %>%
  mutate(county_label = factor(county_label, labels = c('Daggett, Utah\n(750)', '  Torrance, New Mexico\n(15,600)', 'Duchesne, Utah\n(20,100)', 'Chelan, Washington\n(74,800)', 'Ada, Idaho\n(426,000)')) %>% fct_rev) %>%
  ggplot(aes(x = county_label, y = diff, color = estimate_label)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_y_continuous(name = 'Change in wildfire risk estimate relative to naïve method',
                     limits = c(-1.04, 1.04), labels = scales::percent) + 
  zero_line +
  labs(y = 'Change in wildfire risk estimate relative to naive method') +
  coord_flip() +
  scale_color_manual(name = 'Dasymetric method', values = unname(okabe_colors)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.8, 0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.5)),
        legend.key.size = unit(0.34, 'cm'),
        legend.margin = margin(0, 0.1, 0.1, 0.1, unit = 'cm'),
        legend.background = element_rect(fill = 'white', colour = 'gray50', size = 0.5), 
        plot.margin = unit(c(0.1, 0.2, 1, 2), 'cm'))

### add fancy map thingies
# Also note state_name is in there twice
countymaps <- us_counties(resolution = 'high') %>%
  st_transform(3857) %>%
  set_tidy_names %>%
  mutate(fips = paste0(statefp, countyfp))

make_inset <- function(cofips) {
  statefips <- substr(cofips, 1, 2)
  state <- countymaps %>% filter(statefp == statefips)
  county <- state %>% filter(fips == cofips)
  p <- ggplot() + 
    geom_sf(data = st_geometry(state), fill = NA, size = 0.3) + 
    geom_sf(data = st_geometry(county), fill = 'red', color = NA) +
    theme_void()
  p
}

# Loop through and create the maps
insets_flood <- map(counties_use_fl, make_inset)
insets_wildfire <- map(counties_use_wf, make_inset)


inset_size <- 0.12
inset_x <- 0.03
p_top <- ggdraw(p_flood) +
  draw_plot(insets_flood[[1]], x = inset_x, y = 0.6, width = 0.1, height = 0.1) +
  draw_plot(insets_flood[[2]], x = inset_x, y = 0.45, width = inset_size, height = inset_size) +
  draw_plot(insets_flood[[3]], x = inset_x, y = 0.3, width = inset_size, height = inset_size) +
  draw_plot(insets_flood[[4]], x = inset_x, y = 0.16, width = inset_size, height = inset_size) +
  draw_plot(insets_flood[[5]], x = inset_x, y = 0.03, width = 0.1, height = 0.1)
p_bottom <- ggdraw(p_wildfire) +
  draw_plot(insets_wildfire[[1]], x = inset_x, y = 0.85, width = inset_size, height = inset_size) +
  draw_plot(insets_wildfire[[2]], x = inset_x, y = 0.7, width = inset_size, height = inset_size) +
  draw_plot(insets_wildfire[[3]], x = inset_x, y = 0.55, width = inset_size, height = inset_size) +
  draw_plot(insets_wildfire[[4]], x = inset_x, y = 0.4, width = inset_size, height = inset_size) +
  draw_plot(insets_wildfire[[5]], x = inset_x, y = 0.27, width = inset_size, height = inset_size)

p_all <- plot_grid(p_top, p_bottom, nrow = 2)

ggsave('~/temp/logscale_pop_fig_with_insets.png', p_all, height = 6, width = 6, dpi = 400)
