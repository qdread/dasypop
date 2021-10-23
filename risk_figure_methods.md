# Methods for case study: how does the dasymetric method affect our inference about environmental hazards?

For the case study, we explored how the method used to dasymetrically estimate population distribution might affect policy-relevant inference. We investigated two key environmental hazards: wildfire and coastal flooding. Raster data products are available for both wildfire and coastal flooding. The U.S. Forest Service produced a map of wildfire hazard potential (WHP) for the contiguous United States at 270-meter pixel resolution, with five risk categories (Dillon et al. 2015). The U.S. Federal Emergency Management Agency (FEMA) has released flood risk data products for many U.S. counties, including the water surface elevation (WSE) for a 1% flood event (expected to occur once every 100 years). The WSE product is provided at 10-meter pixel resolution (FEMA 2021).

For each of the two hazard categories, we compared estimation methods from four different sources: the present study, the U.S. Environmental Protection Agency (EPA; Pickard et al. 2015), Huang et al. (2021), and Facebook/CIESIN (Tiecke et al. 2017). In addition, we compared all the methods to a fifth baseline method: assuming that individuals are evenly distributed across the entire geographical area of a Census block group (Zandbergen et al. 2011).

### Selection of counties for case study

For the wildfire case study, we took a random sample, stratified by population, of 15 counties in the eleven western states of the contiguous U.S. (WA, OR, CA, ID, NV, MT, WY, CO, UT, AZ, NM), sampling three counties in each population quintile. We chose five counties to display in the final visualization that have sufficient spatial variation in wildfire risk to differentiate between the population methods. For the flood case study, we took a population-stratified random sample from all counties bordering a coastline in the contiguous U.S. However, FEMA does not provide flood risk data products for all of these counties, so we continued sampling until we had ten counties for which flood risk data products were available. Not all states bordering a coastline have flood risk data products available for any counties. Again, we chose five to display that maximize differentiation between the population methods.

### Data sources

We obtained the WHP raster product for the entire contiguous United States, and the 1% flood event WSE product for each of the counties in the case study. We obtained the gridded population estimates for the three comparison methods as raster layers covering the contiguous United States. Finally, we used the previously obtained population estimates for 2016 from the American Community Survey for each of the counties chosen for the case study, as well as the boundaries of each block group as a polygon layer (see method section in main text for references). Locations where we downloaded each data product are in the supplemental table below.

### Initial raster processing

We clipped the WHP raster layer and the population raster layers (U.S. EPA, Huang et al., and Facebook) to the extent of each county in the case study. The water surface elevation rasters were already provided at the single county level. For simplicity, we converted both environmental raster layers to binary form (i.e., at risk and not at risk). For the wildfire layer, we treated all pixels in the medium, high, and very high risk categories as being at risk, and the remainder as not at risk. For the flooding layer, we treated all pixels with water surface elevation > 0 as being at risk. Next, we converted the wildfire and flooding rasters to polygons by merging all adjacent pixels with the same value into a polygon. Finally, we transformed these polygon layers into the coordinate reference system of each of the population rasters.

### Finding population totals in each risk category

We overlaid the polygonized wildfire and flood layers onto the population raster layers for each dasymetric estimation method. For each wildfire and flood polygon, we summed the counts of individuals in all pixels of the population raster contained within that polygon, then calculated the grand totals for each risk category in each county.

For the block group population polygons, we calculated the areas of intersection of each block group polygon with each wildfire or flood polygon. We multiplied the population total of the block group polygon by the proportional area of overlap of each environmental risk polygon to yield the population total at risk within each block group, then calculated the grand totals for each county.

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
