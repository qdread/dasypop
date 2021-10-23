# Methods figure with intermediate steps

# 1. Block group populations
# 2. Superimpose impervious surface on this
# 3. Remove roads from impervious surface
# 4. Final dasymetric map

# Modified 12 Oct 2021 to meet journal formatting guidelines

library(tidycensus)
library(tidyverse)
library(raster)
library(sf)
library(glue)
library(gdalUtils)
library(stars)

### State and county ID used for the example (can be changed)
stid <- '24'
ctyid <- '003'


# Produce data for map (adapted from get_dasy_data.R) ---------------------

# Albers equal-area projection
aea <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

census_api_key(readLines('/research-home/qread/censusapikey.txt')) 
# This is done to not have the API key in your environment or scripts (good practice)

pop <- get_acs(geography = "block group", variables = "B01003_001", 
               year = 2016, state= stid, county = ctyid, 
               geometry = TRUE)   

# Data QC: remove empty geometries from pop
pop <- pop[!is.na(st_dimension(pop)), ]

# Project population to Albers equal-area
pop.projected <- st_transform(pop, crs = aea)

# Update: use gdalwarp to extract the county area, from the NLCD impervious raster, already in Albers projection
# Use a temporary directory I created for the purpose to write the county polygon for extraction. 
nlcd_imp_vrt <- '/nfs/public-data/NLCD/VRTs/NLCD_2016_Impervious_L48_20190405.vrt'
temp_polygon_filename <- as.character(glue("/nfs/qread-data/temp/county-{stid}-{ctyid}.gpkg"))
temp_nlcdraster_filename <- as.character(glue("/nfs/qread-data/temp/countynlcd-{stid}-{ctyid}.tif"))
if (!file.exists(temp_nlcdraster_filename)) {
  st_write(st_union(pop.projected), dsn = temp_polygon_filename, driver = 'GPKG')
  gdalwarp(srcfile = nlcd_imp_vrt, dstfile = temp_nlcdraster_filename, cutline = temp_polygon_filename, crop_to_cutline = TRUE, tr = c(30, 30), dstnodata = "None")
}
lu <- raster(temp_nlcdraster_filename)

#download 2010 block-level data, filter for only the blocks with 0 pop
zero.pop <- get_decennial(geography = "block", variables = "P001001", 
                          year = 2010, state = stid, county = ctyid, 
                          geometry = TRUE) %>% filter(value == 0) %>% st_transform(., aea)

#Remove NLCD data <=1% (masking no longer necessary as it's already masked)
lu[lu <= 1] <- NA

#create lu ratio
lu.ratio <- lu/100

#mask out zero pop blocks
lu.ratio.zp <- mask(lu.ratio, as(zero.pop, "Spatial"), inverse=TRUE)

#get the impervious surface descriptor dataset from: https://www.mrlc.gov/data?f%5B0%5D=category%3ALand%20Cover&f%5B1%5D=category%3AUrban%20Imperviousness&f%5B2%5D=year%3A2016
# Now VRT is used.
imp.surf.desc <- raster("/nfs/public-data/NLCD/NLCD_2016_Impervious_descriptor_L48_20190405/NLCD_2016_Impervious_descriptor_L48_20190405.img")
#mask out primary, secondary, and urban tertiary roads
imp.surf.crop <- raster::crop(imp.surf.desc, as(pop.projected, "Spatial")) #crop imp surface to county
#plot(imp.surf.crop)
imp.surf.mask <- raster::mask(imp.surf.crop, as(pop.projected, "Spatial")) #mask all non-county values to NA

reclass.table <- matrix(c(1,6,1,7,14,NA), ncol = 3, byrow = TRUE) #reclassify values 1-6 into 1 for keep drop the rest

imp.roads <- reclassify(imp.surf.mask, reclass.table, right = NA)
imp.roads.p <- projectRaster(imp.roads, lu.ratio.zp, method = 'ngb')#have to reproject the descriptor file
#Mask out roads (i.e, all NonNA values in imp.roads.p)
RISA <- overlay(lu.ratio.zp, imp.roads.p, fun = function(x, y) {
  x[!is.na(y[])] <- NA
  return(x)
})

#get the block-group level sum of the remaining impervious surface pixels
RISA.sum <- raster::extract(RISA, as(pop.projected,"Spatial"), fun=sum, na.rm=TRUE,df=TRUE)

pop.df <- cbind(pop.projected, RISA.sum$layer)
bg.sum.pop <- fasterize::fasterize(pop.projected, RISA, field = "estimate")
bg.sum.RISA <- fasterize::fasterize(pop.df, RISA, field = "RISA.sum.layer")

#generate density (people/30 m pixel)
dasy.pop <- (bg.sum.pop/bg.sum.RISA) * RISA



# Data objects for maps -------------------------------------------------------------

# Reclassify road versus not road
reclass.table2 <- matrix(c(1,6,1,7,11,2), ncol = 3, byrow = TRUE) #reclassify values 1-6 into 1 for keep drop the rest

imp.roads2 <- reclassify(imp.surf.mask, reclass.table2, right = NA)
imp.roads.p2 <- projectRaster(imp.roads2, lu.ratio.zp, method = 'ngb')#have to reproject the descriptor file

lu_stars <- st_as_stars(lu)
imp_stars <- st_as_stars(imp.roads.p2)
imp_stars[[1]] <- as.character(imp_stars[[1]])
imp_stars[imp_stars == 0] <- NA
dasy_stars <- st_as_stars(dasy.pop)

# Create maps with NEON imagery -------------------------------------------

library(ggspatial)
library(gridExtra)

box_x <- c(-76.60, -76.505)
box_y <- c(38.90, 38.93)

sercimage_raster <- stack('/nfs/qread-data/DASY/neon_imagery/2017_SERC_3_all_5m_UTM_geo.tif')

plot_box <- st_as_sf(data.frame(x = box_x, y = box_y), coords = c('x', 'y'), crs = 4326) %>%
  st_transform(crs = aea) %>%
  st_bbox

maptheme <- theme(legend.position = 'bottom', legend.text = element_text(size = rel(0.5)), legend.title = element_text(size = rel(0.6)),
                  axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
mapcoord <- coord_sf(crs = aea, xlim = plot_box[c(1,3)], ylim = plot_box[c(2,4)], expand = FALSE)
image_alpha <- 1

# Get block group population limits inside the plot_box so that we can set appropriate limits on the fill scale.
range_pop <- st_crop(pop.projected, plot_box) %>% pull(estimate) %>% range # 690 to 3008

# Get dasymetric population limits inside the plot_box so that we can set appropriate limits on the fill scale.
plot_box_dasy <- st_as_sf(data.frame(x = box_x, y = box_y), coords = c('x', 'y'), crs = 4326) %>%
  st_transform(crs = st_crs(dasy_stars)) %>%
  st_bbox
max_dasy <- ceiling(max(dasy_stars[plot_box_dasy][[1]], na.rm = TRUE)) # 23

# Create scale bar manually, draw segment 1 km long
scale_bar_line <- annotate(geom = 'errorbarh', xmin = plot_box['xmax'] - 1500, xmax = plot_box['xmax'] - 500,
                           y = plot_box['ymin'] + 500, color = 'white', size = 0.75, height = 150)
scale_bar_text <- annotate(geom = 'text', label = '1 km', x = plot_box['xmax'] - 2000, y = plot_box['ymin'] + 500, color = 'white', size = rel(2))

letter_label <- function(letter) annotate(geom = 'text', label = letter, hjust = -0.1, vjust = 1.1, x = -Inf, y = Inf, color = 'white')

p1map <- ggplot() +
  annotation_spatial(data = sercimage_raster, alpha = image_alpha) +
  geom_sf(data = pop.projected, color = 'white', alpha = 0.3, aes(fill = estimate)) +
  letter_label('a') +
  scale_fill_viridis_c(name = 'block group  \npopulation   ', option = 'B', limits = range_pop) +
  mapcoord + maptheme

p2map <- ggplot() +
  annotation_spatial(data = sercimage_raster, alpha = image_alpha) +
  geom_stars(data = lu_stars/100) +
  geom_sf(data = pop.projected, color = 'white', alpha = 0.3, fill = NA, size = 0.3) +
  letter_label('b') +
  scale_fill_viridis_c(name = 'impervious surface     \npercentage', labels = scales::percent, na.value = 'transparent', option = 'D') +
  mapcoord + maptheme

p3map <- ggplot() +
  annotation_spatial(data = sercimage_raster, alpha = image_alpha) +
  geom_stars(data = imp_stars) +
  geom_sf(data = pop.projected, color = 'white', alpha = 0.3, fill = NA, size = 0.3) +
  letter_label('c') +
  scale_fill_manual(name = 'surface type    ', labels = c('road (discarded)    ', 'non-road (kept)    '), na.value = 'transparent', na.translate = FALSE, values = rev(colorspace::diverging_hcl(palette='Berlin',n=2))) +
  mapcoord + maptheme 

p4map <- ggplot() +
  annotation_spatial(data = sercimage_raster, alpha = image_alpha) +
  geom_stars(data = dasy_stars) +
  geom_sf(data = pop.projected, color = 'white', alpha = 0.3, fill = NA, size = 0.3) +
  scale_bar_line + scale_bar_text + 
  letter_label('d') +
  scale_fill_viridis_c(name = 'dasymetric\npopulation     ', na.value = 'transparent', option = 'B', limits = c(0.5, max_dasy), trans = 'log10') +
  mapcoord + maptheme

# Bind together so that main panels line up (2x2)
# Note: extra spaces added at the end of all the legend names and values as a workaround to the odd overlapping of text and boxes
allmaps <- gtable_rbind(
  gtable_cbind(ggplotGrob(p1map), ggplotGrob(p2map)), 
  gtable_cbind(ggplotGrob(p3map), ggplotGrob(p4map))
)
ggsave('/nfs/qread-data/DASY/figs/methods_figure_formatted.png', allmaps, height = 5.5*.9, width = 6*.9, dpi = 300)

