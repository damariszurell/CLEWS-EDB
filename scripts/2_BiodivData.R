# Session 2: Biodiversity data

library(raster)

# Load the shapefile
(shrew <- shapefile('data/IUCN_Sorex_alpinus.shp'))

# Plot the Central Europe
library(maps)
map('world',xlim=c(5,30), ylim=c(40,55))

# Overlay the range of the Alpine Shrew
plot(shrew, col='red', add=T)


# Read shapefile for all mammals using the raster package:
mammals <- shapefile('data/MAMMTERR.shp')


## ------------------------------------------------------------------------------------------------
mammals

# Inspect attribute table
head(mammals@data)

 
# Show all entries for the species 'Lynx lynx'
subset(mammals@data, BINOMIAL=='Lynx lynx')
 
# Show all entries for the species with the word 'Lynx' in their name
grep('Lynx',mammals@data$BINOMIAL)
mammals@data[grep('Lynx',mammals@data$BINOMIAL),]


## ------------------------------------------------------------------------------------------------
# Range map of the Eurasian lynx
lynx_lynx <- mammals[mammals@data$BINOMIAL=='Lynx lynx',]

# Map the range
library(maps)
map('world')
plot(lynx_lynx, col='red', add=T)


## ----message=F-----------------------------------------------------------------------------------
library(BIEN)
library(sp)

# Load the range map for the monkey puzzle
(monkey_puzzle <- BIEN_ranges_load_species('Araucaria_araucana'))

# Map
map('world',xlim = c(-180,-20),ylim = c(-60,80))
plot(monkey_puzzle,col='red',add=T)


## ------------------------------------------------------------------------------------------------
# Range area of alpine shrew in square meters:
area(shrew)

# Range area of alpine shrew in square kilometers:
area(shrew)/1000000

# Range area of monkey puzzle in square meters:
area(monkey_puzzle)

# Range area of monkey puzzle n square kilometers:
area(monkey_puzzle)/1000000


## ------------------------------------------------------------------------------------------------
library(rgeos)

# Range centroid:
gCentroid(shrew)

# Map the species range and add the centroid to the map
map('world',xlim=c(5,30), ylim=c(40,55))
plot(shrew, col='red', add=T)
plot(gCentroid(shrew), col='blue',add=T,lwd=3)


# ------------------------------------------------------------------------
# Let's crop the range polygons to South America
monkey_puzzle_SAm <- gIntersection(monkey_puzzle, as(extent(-85, -30, -55, 5),"SpatialPolygons"))

# Map the range and range centroid
map('world',xlim = c(-100,-10),ylim = c(-60,15))
plot(monkey_puzzle_SAm,col='red',add=T)
plot(gCentroid(monkey_puzzle_SAm), col='blue',add=T,lwd=3)


## ------------------------------------------------------------------------------------------------
# By default, raster() will create a 1° resolution map in the *WGS 84* coordinate system (lon/lat).
(r_1deg <- raster())

(shrew_1deg <- rasterize(shrew, r_1deg))

map('world',xlim=c(5,30), ylim=c(40,55))
plot(shrew, col='red', add=T)
plot(shrew_1deg, add=T, alpha=0.6, legend=F)


## ------------------------------------------------------------------------------------------------
# Define an empty raster of the world at 2° spatial resolution
(r_2deg <- raster(res=2))

# Rasterize the eurasian lynx data
lynx_lynx_2deg <- rasterize(lynx_lynx, r_2deg)

# Map the occupied grid cells
plot(lynx_lynx_2deg)


## ------------------------------------------------------------------------------------------------
values(lynx_lynx_2deg)[!is.na(values(lynx_lynx_2deg))] <- 1

# Map the occupied cells
map('world')
plot(lynx_lynx_2deg, add=T, legend=F)



library(letsR)

# The lets.presab() function expects specific column names in the Polygons data frame
colnames(monkey_puzzle@data) <- "binomial"

# We set the resolution to 1 degree (the default) and restrict the spatial extent to South America
r_monkey_puzzle <- lets.presab(monkey_puzzle, resol=1, xmn = -100, xmx = -10, ymn = -57, ymx = 15)

# Map the range and range centroid
map('world',xlim = c(-100,-10),ylim = c(-60,15))
plot(monkey_puzzle_SAm,col='blue',add=T)
plot(r_monkey_puzzle, add=T, alpha=0.6, legend=F)



# Extract the available Pinus species names
(pinus_names <- BIEN_ranges_genus("Pinus",match_names_only = T)[,1])

# Download the range maps for all Pinus species
pinus <- BIEN_ranges_load_species(pinus_names)

# Format the column names and rasterise
colnames(pinus@data) <- "binomial"
r_pinus <- lets.presab(pinus, resol=1, xmn = -170, xmx = -10, ymn = -57, ymx = 60)

# Plot species richness
plot(r_pinus)


## ------------------------------------------------------------------------------------------------
# Subset the SpatialPolygonsDataFrame
artibeus_spp <- mammals[grep('Artibeus',mammals@data$BINOMIAL),]

# Rasterize the ranges using the letsR package
r_artibeus_spp <- lets.presab(artibeus_spp, resol=2, 
                              presence = 1, origin = 1, seasonal = 1)

# Map the species richness
plot(r_artibeus_spp)

# Map single species - here, just the first two
par(mfrow=c(1,2))
plot(r_artibeus_spp, name = "Artibeus amplus")
plot(r_artibeus_spp, name = "Artibeus anderseni")

# Look at structure of the object and at the presence-absence matrix
str(r_artibeus_spp, 1)
head(r_artibeus_spp$Presence_and_Absence_Matrix)


library(rgbif)

# Check out the number of occurrences found in GBIF:
occ_count()

# number of observations:
occ_count(basisOfRecord='OBSERVATION')

# number of occurrences reported for Germany:
occ_count(country=isocodes[grep("Germany", isocodes$name), "code"])

# number of observations reported for Germany:
occ_count(country=isocodes[grep("Germany", isocodes$name), "code"],basisOfRecord='OBSERVATION')


# Check for synonyms
name_suggest(q='Sorex alpinus', rank='species')

# Check number of records - here filtered to those with coordinate information
occ_search(scientificName = "Sorex alpinus", hasCoordinate=T, limit = 10)


## ------------------------------------------------------------------------------------------------
occ_search(scientificName = "Sorex alpinus", hasCoordinate=T, basisOfRecord='HUMAN_OBSERVATION', limit = 10)


gbif_shrew <- occ_search(scientificName = "Sorex alpinus", hasCoordinate=T, basisOfRecord='HUMAN_OBSERVATION', limit = 600) 

# We are just interested in the data frame containing the records
gbif_shrew <- gbif_shrew$data

library(maptools)
data(wrld_simpl)
plot(wrld_simpl,xlim=c(5,50), ylim=c(40,55))
points(gbif_shrew$decimalLongitude, gbif_shrew$decimalLatitude, col='red',  pch=19)


library(CoordinateCleaner)

# We use only those data entries with coordinate information - Note that you don't need this if you have used the hasCoordinate=T in the occ_search() function:
gbif_shrew <- subset(gbif_shrew, !is.na(decimalLatitude))

# We now clean the coordinates and check for outliers - see ?clean_coordinates for more options
cl_gbif_shrew <- clean_coordinates(gbif_shrew, lon="decimalLongitude", lat="decimalLatitude", countries="countryCode", tests=c("centroids", "outliers", "duplicates", "institutions"), inst_rad = 1000)

plot(wrld_simpl,xlim=c(5,50), ylim=c(40,55))
points(gbif_shrew$decimalLongitude, gbif_shrew$decimalLatitude, col='red',  pch=19)
points(gbif_shrew$decimalLongitude[cl_gbif_shrew$.summary], gbif_shrew$decimalLatitude[cl_gbif_shrew$.summary], col='blue',  pch=18)

gbif_shrew <- gbif_shrew[cl_gbif_shrew$.summary,]


save(gbif_shrew,file='data/gbif_shrew.RData')

