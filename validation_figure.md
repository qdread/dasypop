---
title: "Validation figure code"
author: "Quentin D. Read"
date: "10/25/2021"
output: html_document
---

For the case study, we explored how the method used to dasymetrically estimate population distribution might affect policy-relevant inference. We investigated two key environmental hazards: wildfire and coastal flooding. Raster data products are available for both wildfire and coastal flooding. The U.S. Forest Service produced a map of wildfire hazard potential (WHP) for the contiguous United States at 270-meter pixel resolution, with five risk categories (Dillon et al. 2015). The U.S. Federal Emergency Management Agency (FEMA) has released flood risk data products for many U.S. counties, including the water surface elevation (WSE) for a 1% flood event (expected to occur once every 100 years). The WSE product is provided at 10-meter pixel resolution (FEMA 2021).

For each of the two hazard categories, we compared estimation methods from four different sources: the present study, the U.S. Environmental Protection Agency (EPA; Pickard et al. 2015), Huang et al. (2021), and Facebook/CIESIN (Tiecke et al. 2017). In addition, we compared all the methods to a fifth baseline method: assuming that individuals are evenly distributed across the entire geographical area of a Census block group (Zandbergen et al. 2011).

The procedure is:

- Take a random sample of the counties of interest, stratified by population.
- Clip population rasters and environmental rasters to the boundaries of each selected county. 
- Convert the environmental rasters to polygons, merging all adjacent pixels with the same class
- For each polygon, estimate population totals for each method by summing the dasymetric population pixels in that polygon.
- Sum all polygons in each risk class in each county.
- Qualitatively compare risk estimates across methods for each county.

## Selection of counties for case study

For the wildfire case study, we took a random sample, stratified by population, of 15 counties in the eleven western states of the contiguous U.S. (WA, OR, CA, ID, NV, MT, WY, CO, UT, AZ, NM), sampling three counties in each population quintile. We chose five counties to display in the final visualization that have sufficient spatial variation in wildfire risk to differentiate between the population methods. For the flood case study, we took a population-stratified random sample from all counties bordering a coastline in the contiguous U.S. However, FEMA does not provide flood risk data products for all of these counties, so we continued sampling until we had ten counties for which flood risk data products were available. Not all states bordering a coastline have flood risk data products available for any counties. Again, we chose five to display that maximize differentiation between the population methods.

## Data sources

We obtained the WHP raster product for the entire contiguous United States, and the 1% flood event WSE product for each of the counties in the case study. We obtained the gridded population estimates for the three comparison methods as raster layers covering the contiguous United States. Finally, we used the previously obtained population estimates for 2016 from the American Community Survey for each of the counties chosen for the case study, as well as the boundaries of each block group as a polygon layer. Locations where we downloaded each data product are in the supplemental table below.

## Procedure

The following sections include a description of the methods as well as the R code. The random sampling of counties and downloading of the external data sources are not included in this notebook. Only the final section of the code (final data processing and figure generation) is run at the time the notebook is rendered. See the table at the bottom of this document for URLs where you can download the other external data. The only external data object provided is the USA county boundaries file geopackage, projected to Albers equal-area.

### Note on downloading flood risk rasters

Note in particular that the FEMA 1% flood water surface elevation raster files for each county need to be downloaded from the [FEMA flood map service center](https://msc.fema.gov/portal/advanceSearch). For each county, select the appropriate state name and county name, then select "All Jurisdictions" for the community name and press `Search`. Expand the `Flood Risk Products` folder, then the `Flood Risk Database` submenu. Download the compressed folder called `FRD_xxxxx_GeoTIFFs` where `xxxxx` is the county code, then extract it to `input_data/FRD` in the current working directory.

This code block defines the county names and FIPS codes used for the case studies.

```
counties_wildfire <- read.delim(textConnection('
county_name	state_name	state	state_code	county_code
Daggett County	Utah	UT	49	009
Torrance County	New Mexico	NM	35	057
Duchesne County	Utah	UT	49	013
Chelan County	Washington	WA	53	007
Ada County	Idaho	ID	16	001'), colClasses = 'character')
counties_flood <- read.delim(textConnection('
state_name	county_name	fips	state_code	county_code	wse_filename
Maryland	Anne Arundel	24003	24	003	input_data/FRD/FRD_24003C_Coastal_GeoTIFFs_20150909/WSE_01pct.tif
Georgia	Bryan	13029	13	029	input_data/FRD/FRD_13029C_GeoTiffs_20131216/WSE_01pctTif.tif
Florida	Flagler	12035	12	035	input_data/FRD/FRD_12035C_GeoTIFF_20180117/FRD_12035C_GeoTIFF_20180117/WSE_01pct.tif
Maryland	Dorchester	24019	24	019	input_data/FRD/FRD_24019C_Coastal_GeoTIFFs_20160120/WSE_01pct.tif
Alabama	Baldwin	01003	01	003	input_data/FRD/FRD_01003C_GeoTIFFS_20190321/FRD_01097C_GeoTIFFS_20190321/wse_01pct.tif'), colClasses = 'character')
```

This code block creates the directory structure to hold the external data objects. Download and extract them to the appropriate directory.

```
dir.create('input_data/epadasy')
dir.create('input_data/huang_grid')
dir.create('input_data/facebookpop')
dir.create('input_data/WHP')
dir.create('input_data/FRD')
dir.create('input_data/countybounds')
dir.create('input_data/countyrasters')
```


### Initial raster processing

We clipped the WHP raster layer and the population raster layers (U.S. EPA, Huang et al., and Facebook) to the extent of each county in the case study. The water surface elevation rasters were already provided at the single county level. For simplicity, we converted both environmental raster layers to binary form (i.e., at risk and not at risk). For the wildfire layer, we treated all pixels in the medium, high, and very high risk categories as being at risk, and the remainder as not at risk. For the flooding layer, we treated all pixels with water surface elevation > 0 as being at risk. Next, we converted the wildfire and flooding rasters to polygons by merging all adjacent pixels with the same value into a polygon. Finally, we transformed these polygon layers into the coordinate reference system of each of the population rasters.

This code block loads the R packages needed.

```
library(stars)
library(sf)
library(purrr)
library(dplyr)
library(glue)
library(tidycensus)
library(readr)
library(ggplot2)
library(gdalUtils)
library(rslurm)
library(tidyr)
library(grid)
library(gridExtra)
library(USAboundaries)
library(cowplot)
library(forcats)

counties_wildfire <- counties_wildfire %>% mutate(fips = paste0(state_code, county_code))
```

This code block creates individual GeoPackage shapefiles for the boundaries of each county (Albers equal-area projection). In addition a shapefile in unprojected latitude-longitude coordinates is created for each county for use with the unprojected population raster provided by Facebook.

```
gpkgpath <- "input_data/countybounds"
raster_proj <- gdalsrsinfo('input_data/facebookpop/population_usa_2019-07-01.vrt', as.CRS = TRUE)

all_fips <- bind_rows(counties_flood, counties_wildfire) %>%
  select(state_code, county_code)

pwalk(all_fips, function(state_code, county_code) {
  fips <- paste0(state_code, county_code)
  # Extract only one county from the all county geopackage
  ogr2ogr(src_datasource_name = 'input_data/USA_county_2014_aea.gpkg',
          dst_datasource_name = glue('{gpkgpath}/county{fips}.gpkg'),
          overwrite = TRUE, f = 'GPKG', where = glue('fips=\'{fips}\''))
  # Project the single county GPKG to the geographic CRS used by Facebook
  ogr2ogr(src_datasource_name = glue('{gpkgpath}/county{fips}.gpkg'),
          dst_datasource_name = glue('{gpkgpath}/county{fips}longlat.gpkg'),
          overwrite = TRUE, f  = 'GPKG', t_srs = raster_proj)
})
```

Our dasymetric rasters and the flood risk rasters are already in separate files for each county. The following code block creates separate raster files for the wildfire raster, the EPA dasymetric raster, the Facebook raster, and the Huang et al. population grid. We also ensure that our dasymetric raster matches the long-lat geographic CRS by creating a new version.

```
tifpath <- "input_data/countyrasters"
gpkgpath <- "input_data/countybounds"

pwalk(all_fips, function(state_code, county_code) {
  fips <- paste0(state_code, county_code)
  # Wildfire
  gdalwarp(srcfile = 'input_data/WHP/Data/whp2020_GeoTIF/whp2020_cls_conus.tif',
           dstfile = glue('{tifpath}/county{fips}_wildfire.tif'),
           overwrite = TRUE, crop_to_cutline = TRUE,
           cutline = glue('{gpkgpath}/county{fips}.gpkg'))
  # EPA Dasymetric
  gdalwarp(srcfile = 'input_data/epadasy/dasymetric_us_20160208/dasymetric_us_20160208.tif',
           dstfile = glue('{tifpath}/county{fips}_epadasy.tif'),
           overwrite = TRUE, crop_to_cutline = TRUE,
           cutline = glue('{gpkgpath}/county{fips}.gpkg'))
  # Facebook Dasymetric
  gdalwarp(srcfile = 'input_data/facebookpop/population_usa_2019-07-01.vrt',
           dstfile = glue('{tifpath}/county{fips}_fblonglat.tif'),
           overwrite = TRUE, crop_to_cutline = TRUE,
           cutline = glue('{gpkgpath}/county{fips}longlat.gpkg'))
  # Huang Dasymetric
  gdalwarp(srcfile = 'input_data/huang_grid/PopGrid.tif',
           dstfile = glue('{tifpath}/county{fips}_huangdasy.tif'),
           overwrite = TRUE, crop_to_cutline = TRUE,
           cutline = glue('{gpkgpath}/county{fips}.gpkg'))
  # Our dasymetric, crop to cutline with the longlat object
  gdalwarp(srcfile = glue('output_tifs/neon-dasy-{state_code}-{county_code}.tif'),
           dstfile = glue('{tifpath}/county{fips}_dasy.tif'),
           overwrite = TRUE, crop_to_cutline = TRUE, tr = c(30, 30),
           cutline = glue('{gpkgpath}/county{fips}longlat.gpkg'))
})
```

In this code block, the environmental and population rasters are read into R as `stars` objects.

```
tifpath <- 'input_data/countyrasters'
dasypath <- 'output_tifs'

wildfire_rasters <- map(counties_wildfire$fips, ~ read_stars(glue('{tifpath}/county{.}_wildfire.tif')))
ourdasy_rasters_wf <- map(counties_wildfire$fips, ~ read_stars(glue('{dasypath}/neon-dasy-{substr(., 1, 2)}-{substr(., 3, 5)}.tif')))
epadasy_rasters_wf <- map(counties_wildfire$fips, ~ read_stars(glue('{tifpath}/county{.}_epadasy.tif')))
huang_rasters_wf <- map(counties_wildfire$fips, ~ read_stars(glue('{tifpath}/county{.}_huangdasy.tif')))
fbpop_rasters_wf <- map(counties_wildfire$fips, ~ read_stars(glue('{tifpath}/county{.}_fblonglat.tif')))

flood_rasters <- map(counties_flood$wse_filename, read_stars)
ourdasy_rasters_fl <- map(counties_flood$fips, ~ read_stars(glue('{dasypath}/neon-dasy-{substr(., 1, 2)}-{substr(., 3, 5)}.tif')))
epadasy_rasters_fl <- map(counties_flood$fips, ~ read_stars(glue('{tifpath}/county{.}_epadasy.tif')))
huang_rasters_fl <- map(counties_flood$fips, ~ read_stars(glue('{tifpath}/county{.}_huangdasy.tif')))
fbpop_rasters_fl <- map(counties_flood$fips, ~ read_stars(glue('{tifpath}/county{.}_fblonglat.tif')))
```

In this code block, we obtain the census block group population estimates and geographies from the U.S. Census Bureau API. Note that this requires a valid Census API key. These will be used for the block group area-weighting population estimate. The flood risk rasters have a different coordinate reference system for each county, so the block group population polygons are transformed each to a different one.

```
get_blockgroup_pop <- function(FIPS, crs) {
 
  pop <- get_acs(geography = "block group", variables = "B01003_001", 
                 year = 2016, state= substr(FIPS, 1, 2), county = substr(FIPS, 3, 5), 
                 geometry = TRUE)   
  
  # Data QC: remove empty geometries from pop
  pop <- pop[!is.na(st_dimension(pop)), ]
  
  st_transform(pop, crs = crs)
}

census_api_key(readLines('censusapikey.txt')) 
crs_wildfire <- st_crs(wildfire_rasters[[1]])
blockgroup_pops_wf <- map(counties_wildfire$fips, get_blockgroup_pop, crs = crs_wildfire)

crs_flood <- map(flood_rasters, st_crs)
blockgroup_pops_fl <- map2(counties_flood$fips, crs_flood, get_blockgroup_pop)
```

### Conversion of rasters to polygons

In the following code block, we use `st_as_sf()` to convert the wildfire raster and flood raster (converted to binary) to polygons. The job is run in parallel on a Slurm cluster due to high memory requirements.

```
sjob_wild <- slurm_map(wildfire_rasters, st_as_sf, jobname = 'wf_to_poly', 
                       nodes = 2, cpus_per_node = 4,
                       as_points = FALSE, merge = TRUE)
wildfire_polygons <- get_slurm_out(sjob_wild)
cleanup_files(sjob_wild)

make_poly_and_raster <- function(wse_filename, fips) {
  flood_raster <- raster(wse_filename)
  file_name <- glue('input_data/countyrasters/county{fips}_floodbinary.tif')
  if (!file.exists(file_name)) {
    flood_notna_raster <- st_as_stars(!is.na(flood_raster))
    write_stars(flood_notna_raster, file_name)
  } else {
    flood_notna_raster <- read_stars(file_name)
  }
  flood_notna_poly <- st_as_sf(flood_notna_raster, as_points = FALSE, merge = TRUE)
  return(flood_notna_poly)
}

sjob_flood <- slurm_apply(make_poly_and_raster, counties_flood[,c('wse_filename', 'fips')],
                          pkgs = c('raster', 'stars', 'sf', 'glue'), jobname = 'fl_to_poly',
                          nodes = 4, 
                          slurm_options = list(partition = 'sesync', mem = '100gb'))

flood_polygons <- get_slurm_out(sjob_flood)
cleanup_files(sjob_flood)
```

In this code block, we convert the polygons to each of the rasters' coordinate reference systems. We correct invalid geometries in the lat-long polygons using `st_make_valid()`.

```
wildfire_polygons_ourdasycrs <- map(wildfire_polygons, st_transform, crs = st_crs(ourdasy_rasters_wf[[1]]))
wildfire_polygons_epadasycrs <- map(wildfire_polygons, st_transform, crs = st_crs(epadasy_rasters_wf[[1]]))
wildfire_polygons_huangcrs <- map(wildfire_polygons, st_transform, crs = st_crs(huang_rasters_wf[[1]]))
wildfire_polygons_fbpopcrs <- map(wildfire_polygons, st_transform, crs = st_crs(fbpop_rasters_wf[[1]])) %>%
  map(st_make_valid)

flood_polygons_ourdasycrs <- map(flood_polygons, st_transform, crs = st_crs(ourdasy_rasters_fl[[1]]))
flood_polygons_epadasycrs <- map(flood_polygons, st_transform, crs = st_crs(epadasy_rasters_fl[[1]]))
flood_polygons_huangcrs <- map(flood_polygons, st_transform, crs = st_crs(huang_rasters_fl[[1]]))
flood_polygons_fbpopcrs <- map(flood_polygons, st_transform, crs = st_crs(fbpop_rasters_fl[[1]])) %>%
  map(st_make_valid)

```


### Finding population totals in each risk category

We overlaid the polygonized wildfire and flood layers onto the population raster layers for each dasymetric estimation method. For each wildfire and flood polygon, we summed the counts of individuals in all pixels of the population raster contained within that polygon, then calculated the grand totals for each risk category in each county.

In this code block, we define a function to aggregate the raster pixels by polygon for each population estimation method, then apply it in parallel on a Slurm cluster to each county.

```
aggregate_by_poly <- function(i, type, suffix) {
  sums_ourdasy <- aggregate(get(paste0('ourdasy_rasters_', suffix))[[i]], get(paste0(type,'_polygons_ourdasycrs'))[[i]], FUN = sum, na.rm = TRUE) %>%
    st_as_sf() %>% st_drop_geometry()
  sums_epadasy <- aggregate(get(paste0('epadasy_rasters_', suffix))[[i]], get(paste0(type,'_polygons_epadasycrs'))[[i]], FUN = sum, na.rm = TRUE) %>%
    st_as_sf() %>% st_drop_geometry()
  sums_huang <- aggregate(get(paste0('huang_rasters_', suffix))[[i]], get(paste0(type,'_polygons_huangcrs'))[[i]], FUN = sum, na.rm = TRUE) %>%
    st_as_sf() %>% st_drop_geometry()
  sums_fbpop <- aggregate(get(paste0('fbpop_rasters_', suffix))[[i]], get(paste0(type,'_polygons_fbpopcrs'))[[i]], FUN = sum, na.rm = TRUE) %>%
    st_as_sf() %>% st_drop_geometry()
  cbind(sums_ourdasy, sums_epadasy, sums_huang, sums_fbpop) %>%
    setNames(c('our_dasy', 'epa_dasy', 'huang', 'fb'))
}

sjob_wild_sums <- slurm_apply(aggregate_by_poly, data.frame(i = 1:nrow(counties_wildfire)),
                              type = 'wildfire', suffix = 'wf',
                              global_objects = c('wildfire_polygons_ourdasycrs', 'wildfire_polygons_epadasycrs',
                                                 'wildfire_polygons_huangcrs', 'wildfire_polygons_fbpopcrs',
                                                 'ourdasy_rasters_wf', 'epadasy_rasters_wf',
                                                 'huang_rasters_wf', 'fbpop_rasters_wf'), 
                              jobname = 'agg_wf', nodes = 4, cpus_per_node = 2,
                              pkgs = c('stars', 'sf'))

wildfire_sums <- get_slurm_out(sjob_wild_sums)
cleanup_files(sjob_wild_sums)

sjob_flood_sums <- slurm_apply(aggregate_by_poly, data.frame(i = 1:nrow(counties_flood)),
                              type = 'flood', suffix = 'fl',
                              global_objects = c('flood_polygons_ourdasycrs', 'flood_polygons_epadasycrs',
                                                 'flood_polygons_huangcrs', 'flood_polygons_fbpopcrs',
                                                 'ourdasy_rasters_fl', 'epadasy_rasters_fl',
                                                 'huang_rasters_fl', 'fbpop_rasters_fl'), 
                              jobname = 'agg_fl', nodes = 4, cpus_per_node = 2,
                              pkgs = c('stars', 'sf'))

flood_sums <- get_slurm_out(sjob_flood_sums)
cleanup_files(sjob_flood_sums)
```

For the block group population polygons, we calculated the areas of intersection of each block group polygon with each wildfire or flood polygon. We multiplied the population total of the block group polygon by the proportional area of overlap of each environmental risk polygon to yield the population total at risk within each block group, then calculated the grand totals for each county.

The method was slightly different for wildfire risk and flood risk. For wildfire, we used the rasters to generate tables of pixel counts in each class for each polygon. For flood risk, because there were only two categories, we intersected the polygonized flood risk raster and the block group polygons.

In the following code block, we define the function to generate the tables for each block group polygon and apply it to the wildfire raster for each county in the wildfire case study.

```
table_by_polygon <- function(env_raster, bg_pop_poly) {
  map(st_geometry(bg_pop_poly), ~ as.data.frame(table(st_crop(env_raster, .)[[1]], useNA = 'always')))
}

wildfire_tables_blockgroups <- map2(wildfire_rasters, blockgroup_pops_wf, table_by_polygon)

sum_classes_blockgroup <- function(pop_table, bg_poly) {
  map_dfr(1:length(pop_table), ~ data.frame(st_drop_geometry(bg_poly[., c('GEOID', 'estimate')]), pop_table[[.]])) %>%
    group_by(GEOID) %>%
    mutate(pop = estimate * Freq / sum(Freq))
}

wildfire_sums_blockgroups <- map2(wildfire_tables_blockgroups, blockgroup_pops_wf, sum_classes_blockgroup)
```

In this code block, we define the function to calculate area of intersection of the flood risk polygons and block group polygons, and apply it in parallel on a Slurm cluster to each of the counties in the flood risk case study.

```
get_blockgroup_sums <- function(bg_pop, env_poly, env_crs) {
  env_poly <- st_make_valid(env_poly)
  bg_pop %>%
    st_transform(env_crs) %>%
    mutate(area = st_area(.)) %>%
    st_intersection(env_poly) %>%
    mutate(area_int = st_area(.))
}

sjob_flood_block <- slurm_apply(function(i) get_blockgroup_sums(blockgroup_pops_fl[[i]],
                                                                flood_polygons_ourdasycrs[[i]],
                                                                st_crs(flood_polygons_ourdasycrs[[1]])),
                                params = data.frame(i = 1:length(blockgroup_pops_fl)),
                                global_objects = c('get_blockgroup_sums', 'blockgroup_pops_fl', 'flood_polygons_ourdasycrs'),
                                jobname = 'floodblock', nodes = 4, cpus_per_node = 2,
                                pkgs = c('stars', 'sf', 'dplyr'))

flood_sums_blockgroups <- get_slurm_out(sjob_flood_block)
cleanup_files(sjob_flood_block)
```

### Grand totals by risk category

In this code block, we find the total Census population for each county and join it to the county lookup table, and generate informative labels for the final figure.

```
wildfire_pop_totals <- map_dbl(blockgroup_pops_wf, ~ sum(.$estimate))
flood_pop_totals <- map_dbl(blockgroup_pops_fl, ~ sum(.$estimate))

counties_wildfire <- mutate(counties_wildfire, total_pop = wildfire_pop_totals)
counties_flood <- mutate(counties_flood, total_pop = flood_pop_totals)

# Ordered factor facet labels for counties
counties_wildfire <- counties_wildfire %>%
  mutate(total_pop = wildfire_pop_totals,
         county_label = glue('{county_name}, {state_name} ({as.character(as.integer(signif(total_pop, 3)))})'),
         county_label = factor(county_label, levels = county_label[order(total_pop)]))

counties_flood <- counties_flood %>%
  mutate(total_pop = flood_pop_totals,
         county_label = glue('{county_name}, {state_name} ({as.character(as.integer(signif(total_pop, 3)))})'),
         county_label = factor(county_label, levels = county_label[order(total_pop)]))
```

In this code block, we define a function to calculate grand totals of the risk classes for each county, then apply it to each county in the wildfire case study. Grand totals for the block group area-weighted method are calculated with a different function.

```
grandtotal_classes <- function(env_poly, pop_sums) {
  data.frame(env_class = env_poly[[1]], pop_sums) %>%
    mutate(env_class = as.character(env_class)) %>%
    group_by(env_class) %>%
    summarize_all(sum, na.rm = TRUE)
}

wildfire_grandtotals <- map2(wildfire_polygons, wildfire_sums, grandtotal_classes)

wildfire_grandtotals_blockgroups <- map(wildfire_sums_blockgroups, function(dat) {
  dat %>%
    rename(env_class = Var1) %>%
    mutate(env_class = as.character(env_class)) %>%
    group_by(env_class) %>%
    summarize(blockgroup = sum(pop))
})

wildfire_grandtotals <- map2(wildfire_grandtotals, wildfire_grandtotals_blockgroups, full_join)

wildfire_grandtotals <- counties_wildfire %>% 
  mutate(listcol = wildfire_grandtotals) %>%
  unnest(cols = listcol)
```

In this code block, we do the same for the flood risk case study, again using a different function for the block group area-weighted method.

```
flood_grandtotals <- map2(flood_polygons, flood_sums, grandtotal_classes)

get_blockgroup_grandtotals <- function(env_sums) {
  env_sums %>%
    st_drop_geometry() %>%
    rename_with(~ 'env_class', ends_with('tif') | contains('layer')) %>%
    mutate(env_class = as.character(env_class)) %>%
    group_by(env_class) %>%
    summarize(blockgroup = as.numeric(sum(estimate * area_int/area)))
}

flood_grandtotals_blockgroups <- map(flood_sums_blockgroups, get_blockgroup_grandtotals)

flood_grandtotals <- map2(flood_grandtotals, flood_grandtotals_blockgroups, full_join)

flood_grandtotals <- counties_flood %>% 
  mutate(listcol = flood_grandtotals) %>%
  unnest(cols = listcol)

save(list = c('wildfire_grandtotals', 'flood_grandtotals'), file = 'temp_files/grandtotals_finalfig.RData')
```

## Figure

In this code block, we prepare the data for figure generation. First, we reshape the data to long form for plotting. We also simplify the wildfire risk classes into binary by summing all individuals in risk classes 3, 4, and 5 (medium to very high risk), considering them all to be "at risk."

```
load(file = 'temp_files/grandtotals_finalfig.RData')

label_names <- c('This study', 'U.S. EPA', 'Huang et al.', 'Facebook', 'Block group area weighting')
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
```

The following code block creates the two dot plots, one for flood risk and one for wildfire risk.

```
okabe_colors <- palette.colors(7, palette = 'Okabe-Ito')[c(7, 3, 4, 2, 1)]

p_flood <- flood_risks_reduced %>%
  ungroup %>%
  mutate(county_label = factor(county_label, labels = c('Dorchester, Maryland\n(32,500)', 'Bryan, Georgia\n(34,100)', 'Flagler, Florida\n(103,000)', 'Baldwin, Alabama\n(200,000)', 'Anne Arundel,Maryland\n(560,000)')) %>% fct_rev) %>%
  ggplot(aes(x = county_label, y = population, color = estimate_label)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_y_log10(name = 'Population at risk from flooding',
                limits = c(110, 44000),
                expand = expansion(mult = 0.02), position = 'right') +
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
ggplot(aes(x = county_label, y = population, color = estimate_label)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_y_log10(name = 'Population at risk from wildfire',
                limits = c(110, 44000),
                expand = expansion(mult = 0.02)) +
  coord_flip() +
  scale_color_manual(name = 'Dasymetric method', values = unname(okabe_colors)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.2, 0.25),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.5)),
        legend.key.size = unit(0.34, 'cm'),
        legend.margin = margin(0, 0.1, 0.1, 0.1, unit = 'cm'),
        legend.background = element_rect(fill = 'white', colour = 'gray50', size = 0.5), 
        plot.margin = unit(c(0.1, 0.2, 1, 2), 'cm'))
```

The following code generates the inset maps and draws them along the margins of the plot panels.

```
countymaps <- us_counties(resolution = 'high') %>%
  st_transform(3857) %>%
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

grid.draw(p_all)
```

## Supplemental table: data sources for case studies

| Data product | Provider | Data year | Resolution | Coverage | Native CRS | Source URL | Date accessed |
| ------------ | -------- | --------- | ---------- | -------- | ---------- | ---------- | ------------- |
| wildfire hazard potential (WHP) | U.S. Forest Service | 2020 | 270 m | contiguous United States | CONUS Albers equal-area | https://www.fs.usda.gov/rds/archive/catalog/RDS-2015-0047-3 | 4 June 2021 |
| water surface elevation (WSE) for 1% flood event | U.S. Federal Emergency Management Agency | varies by county | 10 m | individual counties | varies by county | https://www.fema.gov/flood-maps/tools-resources/risk-map/products | 3 August 2021 |
| dasymetric population raster | U.S. Environmental Protection Agency | 2016 | 30 m | contiguous United States | CONUS Albers equal-area | https://www.epa.gov/enviroatlas/data-download-step-2#Individual-datasets | 4 June 2021 |
| dasymetric population raster | Huang et al. 2021 | 2017 | 100 m | contiguous United States | CONUS Albers equal-area | https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/DLGP7Y | 12 July 2021 |
| dasymetric population raster | Facebook / CIESIN | 2019 | 30 m | global (by country) | WGS84 lat-long | https://data.humdata.org/dataset/united-states-high-resolution-population-density-maps-demographic-estimates | 7 May 2021 |

## Works cited

Dillon, G., Menakis, J. & Fay, F. (2015). Wildland fire potential: A tool for assessing wildfire risk and fuels management needs. In: Keane, Robert E.; Jolly, Matt; Parsons, Russell; Riley, Karin. Proceedings of the large wildland fires conference; May 19-23, 2014; Missoula, MT. Proc. RMRS-P-73. Fort Collins, CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. p. 60-76., 73, 60–76.

Huang, X., Wang, C., Li, Z. & Ning, H. (2021). A 100 m population grid in the CONUS by disaggregating census data with open-source Microsoft building footprints. Big Earth Data, 5, 112–133.

Pickard, B.R., Daniel, J., Mehaffey, M., Jackson, L.E. & Neale, A. (2015). EnviroAtlas: A new geospatial tool to foster ecosystem services science and resource management. Ecosystem Services, 14, 45–55.

Tiecke, T.G., Liu, X., Zhang, A., Gros, A., Li, N., Yetman, G., et al. (2017). Mapping the world population one building at a time. arXiv:1712.05839 [cs].

U.S. Federal Emergency Management Agency (FEMA). (2021). Risk MAP Products | FEMA.gov. Available at: https://www.fema.gov/flood-maps/tools-resources/risk-map/products. Last accessed 10 August 2021.

Zandbergen, P.A. (2011). Dasymetric Mapping Using High Resolution Address Point Datasets. Transactions in GIS, 15, 5–27.

