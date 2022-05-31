# Session 1: Spatial data in R

## install.packages(c('sp','raster','dismo'), dep=T)


## ------------------------------------------------------------------------------------------------
# We set a seed for the random number generator, so that we all obtain the same results
set.seed(12345)

coords <- cbind(
  x <- rnorm(50, sd=2),
  y <- rnorm(50, sd=2)
)

# sort coords data along x coordinates
coords <- coords[order(coords[,1]),]

# look at data structure
str(coords)

# plot the coordinates
plot(coords)


## ------------------------------------------------------------------------------------------------
library(sp)

# Convert into SpatialPoints
sp <- SpatialPoints(coords)

# Check out the object
class(sp)

# Inspect the (spatial) information in the object:
showDefault(sp)

# The raster package also provides a nicer print summary, so we load it as well:
library(raster)
sp


## ------------------------------------------------------------------------------------------------
sp <- SpatialPoints(coords, proj4string=CRS('+proj=longlat +datum=WGS84'))
sp


## ------------------------------------------------------------------------------------------------
# Create attribute table
data <- data.frame(ID=1:50,species= NA)
# Define the 15 western trees as birch
data[1:15, 'species'] <- 'birch'
# Randomly assign the species beach or oak to the remaining 35 trees
data[-c(1:15), 'species'] <- sample(c('beech','oak'),35,replace=T)
                  
# Create SpatialPointsDataFrame
(spdf <- SpatialPointsDataFrame(sp, data))

# To see what is inside
str(spdf)


## ----warning=F-----------------------------------------------------------------------------------
# Create SpatialLines through all oak trees in the data
lns <- spLines(subset(spdf,species=='oak'), crs=CRS('+proj=longlat +datum=WGS84'))

# Create SpatialPolygons around all birch trees in the data
# Hull around birch trees:
birch_hull <- chull(subset(spdf,species=='birch')@coords)
# SpatialPolygons:
pols <- spPolygons(subset(spdf,species=='birch')[birch_hull,], crs=CRS('+proj=longlat +datum=WGS84'))


## ------------------------------------------------------------------------------------------------
plot(pols, border='blue',col='yellow', axes=T)
plot(lns, add=T, col='red', lwd=2)
points(sp, lwd=2)



## ------------------------------------------------------------------------------------------------
(shrew <- shapefile('data/IUCN_Sorex_alpinus.shp'))

# Plot Central Europe
library(maps)
map('world',xlim=c(5,30), ylim=c(40,55))

# Overlay the range of the Alpine Shrew
plot(shrew, col='red', add=T)



## ------------------------------------------------------------------------------------------------
(r1 <- raster(ncol=10, nrow=10, xmx=-80, xmn=-150, ymn=20, ymx=60))


## ------------------------------------------------------------------------------------------------
summary(values(r1))
values(r1) <- rnorm(ncell(r1))

# plot the raster
plot(r1)


## ------------------------------------------------------------------------------------------------
# Create another RasterLayer and assign values
r2 <- r1
values(r2) <- 1:ncell(r2)

# Stack the raster layers
(r <- stack(r1,r2))
plot(r)


## ------------------------------------------------------------------------------------------------
(temp <- raster('data/UK_temp.tif'))

plot(temp)


## ------------------------------------------------------------------------------------------------
(bioclim <- stack('data/UK_bioclim.grd'))

plot(bioclim)


# # Download global bioclimatic data from worldclim:
# (clim <- getData("worldclim", var="bio", res=10, download=T, path="data"))

## ------------------------------------------------------------------------------------------------
plot(clim)


## ------------------------------------------------------------------------------------------------
# Crop the temperature layer (bio1) to roughly European extent
temp_eur <- crop(clim[[1]], c(-15,45,35,72))

# Aggregate to one-degree and two-degree resolution
temp_eur_onedeg <- aggregate(temp_eur, 6)
temp_eur_twodeg <- aggregate(temp_eur, 12)

par(mfrow=c(1,3))
plot(temp_eur)
plot(temp_eur_onedeg)
plot(temp_eur_twodeg)


## ------------------------------------------------------------------------------------------------
# Generate random locations
library(dismo)
lonlat <- randomPoints(temp_eur,10)

# Map temperature and locations
plot(temp_eur)
points(lonlat, cex=2, pch=19)

# Extract temperature values at these locations
raster::extract(temp_eur, lonlat)

