# Session 3: Species threat data 

# Read in mammal shapefile
mammals <- shapefile('data/MAMMTERR.shp')

# rough European extent
extent_eur <- c(-15,45,35,72)

# Crop SpatialPolygonsDataFrame to European extent
mammals_eur <- crop(mammals, extent_eur)

# Rasterise mammal range maps for European extent
r_mammals_eur <- lets.presab(mammals_eur, resol=1, presence = 1, origin = 1, seasonal = 1)

# Save the PresenceAbsence object
save(extent_eur, r_mammals_eur,file='data/r_mammals_eur.RData')


## ------------------------------------------------------------------------------------------------
load('data/r_mammals_eur.RData')

# Remember how the PresenceAbsence objects look like:
str(r_mammals_eur,2)

# Plot the richness raster:
plot(crop(r_mammals_eur$Richness_Raster,extent_eur))


## ------------------------------------------------------------------------------------------------
# Extract the species names from the PresenceAbsence object
names_mammals_eur <- colnames(r_mammals_eur$Presence_and_Absence_Matrix)[-c(1:2)]


# library(rredlist)
# 
# # Generate your personal API token
# rl_use_iucn()


----
# Download red list category for single species using your personal API token "MY_IUCN_REDLIST_KEY"
# (rl_search('Lynx lynx', key= MY_IUCN_REDLIST_KEY))


# # Download red list categories for all species
# redlist_status <- do.call(rbind,lapply(names_mammals_eur,FUN=function(sp){
#   rl_search(sp, key= MY_IUCN_REDLIST_KEY)$result
#   }
#   ))


## ------------------------------------------------------------------------------------------------
redlist_status <- read.table('data/mammals_eur_redlist_status.csv', header=T, sep=';')


## ------------------------------------------------------------------------------------------------
redlist_status[1:10,10:20]


## ------------------------------------------------------------------------------------------------
# Conservation status
table(redlist_status$category)

# Population trends
table(redlist_status$population_trend)


## ------------------------------------------------------------------------------------------------
# Download red list threats for single species
# rl_threats('Lynx lynx', key= MY_IUCN_REDLIST_KEY)


## # Download red list threats for all species
## redlist_threats <- do.call(rbind,
##                            lapply(seq_len(length(names_mammals_eur)),FUN=function(i){
##                              xi <- rl_threats(names_mammals_eur[i], key= MY_IUCN_REDLIST_KEY);
##                              if(length(xi$result)) {
##                                data.frame(species=names_mammals_eur[i],xi$result)
##                                }
##                              }
##                              ))


## ------------------------------------------------------------------------------------------------
redlist_threats <- read.table('data/mammals_eur_redlist_threats.csv', header=T, sep=';')


## ------------------------------------------------------------------------------------------------
redlist_threats[sample(nrow(redlist_threats),10),-c(1:2)]


## ------------------------------------------------------------------------------------------------
table(redlist_threats$timing)


## ------------------------------------------------------------------------------------------------
(subset(redlist_status,category=='VU')$scientific_name)


## ----echo=T, warning=F, message=F----------------------------------------------------------------
library(raster)

# Identify all vulnerable species 
vu_spp <- subset(redlist_status,category=='VU')$scientific_name

# Identify all least concern species
lc_spp <- subset(redlist_status,category=='LC')$scientific_name

# Now, we extract the distribution data for the VU and LC species groups, make rasters, stack these and plot
# Calculate species richness
richness_RL_status <- data.frame(r_mammals_eur$Presence_and_Absence_Matrix[,1:2],
  richness_vu = rowSums(r_mammals_eur$Presence_and_Absence_Matrix[,vu_spp]),
  richness_lc = rowSums(r_mammals_eur$Presence_and_Absence_Matrix[,lc_spp]))

# Make raster and plot
spplot( rasterFromXYZ(richness_RL_status))


# Which ongoing threats are the most common ?
sort(table(subset(redlist_threats, species %in% names_mammals_eur & timing=='Ongoing')$title), decreasing=T)[1:10]

# Identify the species experiencing threats from hunting
mammals_threat1 <- subset(redlist_threats,title=="Hunting & trapping terrestrial animals" & species %in% names_mammals_eur)$species
# Identify the species experiencing threats from agriculture, specifically non-timber crops 
mammals_threat2 <- subset(redlist_threats,title=="Annual & perennial non-timber crops" & species %in% names_mammals_eur)$species

# Map species experiencing threats from hunting
plot(rasterFromXYZ(
  data.frame(
    r_mammals_eur$Presence_and_Absence_Matrix[,1:2],
    rowSums(r_mammals_eur$Presence_and_Absence_Matrix[,mammals_threat1]))), 
  main="Hunting & trapping terrestrial animals")

# species experiencing threats from agriculture, specifically non-timber crops 
plot(rasterFromXYZ(
  data.frame(
    r_mammals_eur$Presence_and_Absence_Matrix[,1:2],
    rowSums(r_mammals_eur$Presence_and_Absence_Matrix[,mammals_threat2]))), 
  main="Annual & perennial non-timber crops")

