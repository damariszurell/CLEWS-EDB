# Session 4: Environmental data in R

## ------------------------------------------------------------------------------------------------
library(raster)


clim <- getData("worldclim", var="bio", res=10, download=T, path="data")
 
# Now, let's look at the data:
clim
# Can you explain, what a raster stack is?
plot(clim[[1:2]])

# Please note that you have to set download=T if you haven't downloaded the data before:
clim_fut <- getData('CMIP5', var='bio', download=F, res=10,
                  rcp=45, model='NO', year=50, path=my_filepath)

# Rename to hold the same bioclim names as the current climate layers
names(clim_fut) <- names(clim)

# Plot climate layers (only the first 4 layers for simplicity)
plot(clim_fut[[1:4]])


# Read in the CLC data for Germany:
corine2018 <- raster('data/CLC2018_DE.tif')

# Read in the legend:
legend <- read.table('data/CLC2018_CLC2018_V2018_20_QGIS.txt',sep=',')  # columns 2-4 indicate the RGB colours for reproducing the standard CLC colour maps. Column 6 provides the land use/cover classes.

legend[,6]

# Plot the CLC map
plot(corine2018)


## ------------------------------------------------------------------------------------------------
# Number of cells in x/y direction that will be aggregated
aggregation_factor <- 10

# Calculate proportional cover for class 23 (broad-leaved forest) at 1km resolution
corine2018_1km <- aggregate(corine2018,aggregation_factor, 
                            fun=function(x,...){ sum(x==23)/(aggregation_factor^2) })

# Calculate proportional cover for class 24 (Coniferous forest) and add to raster stack
corine2018_1km <- addLayer(corine2018_1km, 
                           aggregate(corine2018,aggregation_factor, 
                                     fun=function(x,...){ sum(x==24)/(aggregation_factor^2) }))

# Calculate proportional cover for class 25 (Coniferous forest) and add to raster stack
corine2018_1km <- addLayer(corine2018_1km, 
                           aggregate(corine2018,aggregation_factor, 
                                     fun=function(x,...){ sum(x==25)/(aggregation_factor^2) }))

# Add names to layers
names(corine2018_1km) <- c('broadleaved_forest','coniferous_forest','mixed_forest')

# Plot the proportional covers
plot(corine2018_1km)

# And don't forget to save the resulting raster
writeRaster(corine2018_1km, "data/corine2018_forestcover_1km.grd")


# Read in raster stack
ndvi <- stack('data/NDVI_1km_DE.grd')

# plot NDVI layers
plot(ndvi)

