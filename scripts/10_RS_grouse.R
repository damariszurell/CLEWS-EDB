# Session 10: Black grouse range dynamics in RangeShiftR

## ------------------------------------------------------------------------------------------
# Set path to model directory
dirpath = "models/RS_BlackGrouse/"


## ------------------------------------------------------------------------------
# load required packages
library(RangeShiftR)
library(raster)
library(rasterVis)
library(viridis)
library(tidyverse)
library(ggplot2) 

# Create sub-folders (if not already existing)
if(!file.exists(paste0(dirpath,"Inputs"))) {
  dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE) }
if(!file.exists(paste0(dirpath,"Outputs"))) {
  dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE) }
if(!file.exists(paste0(dirpath,"Output_Maps"))) {
  dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE) }


## ------------------------------------------------------------------------------
# Read in landscape files
r_lc2000 <- raster(paste0(dirpath,'Inputs/BlackGrouse_SDMpredictions2000.asc'))
r_lc2010 <- raster(paste0(dirpath,'Inputs/BlackGrouse_SDMpredictions2010.asc'))
r_lc2020 <- raster(paste0(dirpath,'Inputs/BlackGrouse_SDMpredictions2020.asc'))
r_lc2030 <- raster(paste0(dirpath,'Inputs/BlackGrouse_SDMpredictions2030.asc'))
r_lc2040 <- raster(paste0(dirpath,'Inputs/BlackGrouse_SDMpredictions2040.asc'))
r_lc2050 <- raster(paste0(dirpath,'Inputs/BlackGrouse_SDMpredictions2050.asc'))
r_lc2060 <- raster(paste0(dirpath,'Inputs/BlackGrouse_SDMpredictions2060.asc'))
r_lc2070 <- raster(paste0(dirpath,'Inputs/BlackGrouse_SDMpredictions2070.asc'))


## ----------------------------------------------------------------------------
# Read in border mask of Switzerland
CH_mask <- raster(paste0(dirpath,'Inputs/CH_mask.grd'))

# For displaying the rasters, we stack them and then mask them using the Swiss border mask
plot(mask(stack(mget(grep('r_lc',ls(),value=T))), CH_mask))


## ------------------------------------------------------------------------------------------------
# Set static landscape parameters
land <- ImportedLandscape(LandscapeFile = "BlackGrouse_SDMpredictions2000.asc", 
        Resolution = 1000, # resolution of one grid cell in m
        HabPercent = TRUE, # habitat map includes percentage of suitability
        K_or_DensDep = 0.5) # strength of density dependence 1/b


## ------------------------------------------------------------------------------------------------
# Numbers of spinup years in dynamic landscape
spinup <- 25

# Set dynamic landscape parameters
land_cc <- ImportedLandscape(LandscapeFile = c("BlackGrouse_SDMpredictions2000.asc",
                                               "BlackGrouse_SDMpredictions2010.asc",
                                               "BlackGrouse_SDMpredictions2020.asc",
                                               "BlackGrouse_SDMpredictions2030.asc",
                                               "BlackGrouse_SDMpredictions2040.asc",
                                               "BlackGrouse_SDMpredictions2050.asc",
                                               "BlackGrouse_SDMpredictions2060.asc",
                                               "BlackGrouse_SDMpredictions2070.asc"), 
          DynamicLandYears = c(0, spinup, spinup+10, spinup+20, spinup+30, spinup+40, spinup+50, spinup+60), # years to switch landscapes
          Resolution = 1000, # resolution of one grid cell
          HabPercent = TRUE, # habitat map includes percentage of suitability
          K_or_DensDep = 0.5)  # strength of density dependence 1/b



## ------------------------------------------------------------------------------------------------
# Define transition matrix
(trans_mat <- matrix(c(0, 0.6,1.646,   0.5), nrow = 2, byrow = F)) 
# stage 0 : dispersing juveniles; stage 1: adult, reproducing hens


## ------------------------------------------------------------------------------------------------
# Define stage structured module
stg <- StageStructure(Stages = 2, # 2 stages
                      TransMatrix = trans_mat,
                      MaxAge = 10,
                      FecDensDep = F , 
                      SurvDensDep = T , # Only mortality increases with density
                      SurvDensCoeff = 1 # Density dependence relative to 1/b
)

# Define demography module
demo <- Demography(StageStruct = stg, ReproductionType = 0) # asexual model


## ------------------------------------------------------------------------------------------------
plotProbs(stg)


## ------------------------------------------------------------------------------------------------
# Define dispersal module
disp <-  Dispersal(
  # Emigration phase: stage 0 has constant emigration probability of 0.81
  Emigration = Emigration(StageDep=T, EmigProb = cbind(0:1,c(0.81,0)) ),
  # Transfer phase: negative exponential dispersal kernel with mean dispersal distance of 8km
  Transfer = DispersalKernel(Distances = 8000),
  # Settlement: if individual arrives in unsuitable cells, it can randomly chose a suitable neighbouring cell or will die
  Settlement = Settlement(Settle = 2) 
  )


## ------------------------------------------------------------------------------------------------
# Define initial conditions for simulations
init <- Initialise( InitType = 0, FreeType = 0, NrCells = 8000, # free initialization in 8000 random cells 
                    InitDens = 2, # start population at given density (in individuals per cell; next line)
                    IndsHaCell = 1, # set initial density to 2 inds/cell
                    PropStages = c(0, 1) # we initialise only adult hens
)


## ------------------------------------------------------------------------------------------------
# Set the number of replicate simulations to run - for quick illustration, we use few replicates here, typically should be >20 
RepNb <- 10
sim_years <- 100


## ------------------------------------------------------------------------------------------------
sim <- Simulation(Simulation = 0, 
                  Replicates = RepNb, # number of replicates to be run
                  Years = spinup + sim_years, # number of years to be simulated
                  EnvStoch = 2, #include local environmental stochasticity
                  EnvStochType = 0,  # Environmental Stochasticity included in fecundity (in Zurell et al. 2012, in pleadYoung)
                  std=0.15,minR=1.28,maxR=2.53, # standard deviation defined to correspond roughly to pleadYoung variability in Zurell et al. (2012)
                  ac=0, # no temporal autocorrelation
                  OutIntPop = 1,
                  OutIntOcc = 1,
                  OutIntRange = 1)


## ------------------------------------------------------------------------------------------------
# RangeShifter parameter master object for static landscape (we set a seed for replicable results)
s <- RSsim(batchnum = 1, land = land, demog = demo, dispersal = disp, simul = sim, init = init, seed = 324135)


## ------------------------------------------------------------------------------------------------
# RangeShifter parameter master object for dynamic landscape
s_cc <- RSsim(batchnum = 2, land = land_cc, demog = demo, dispersal = disp, simul = sim, init = init, seed = 324135)


## ---------------------------------------------------------------------------------------
# Run simulations with 25 years spin up
RunRS(s, dirpath)


## ------------------------------------------------------------------------------------------------
# plot the resulting abundances and occupancy
par(mfrow=c(1,2))

# Abundances
plotAbundance(s, dirpath)

# Occupancy
plotOccupancy(s, dirpath)


## --------------------------------------------------------------------------------
# Map mean abundances

# Read population output file into a data frame
pop_static <- readPop(s, dirpath)

# Make data frame with mean abundances per cell and selected year:
pop_static_wide <- pop_static %>%
  # Select years you want to map
  filter(Year %in% c(spinup, spinup+sim_years) ) %>%
  # Define grouping variables
  group_by(x,y,Year) %>%
  # Average individuals over replicate runs
  summarise(meanInd = mean(NInd)) %>%
  # Make separate columns for different years
  pivot_wider(names_from=Year, values_from=meanInd)

pop_static_wide


## ------------------------------------------------------------------------------------------------
# Make raster stacks

# Transfer into correct coordinate system
pop_static_wide$x <- pop_static_wide$x + xmin(CH_mask)
pop_static_wide$y <- pop_static_wide$y + ymin(CH_mask)

# Make raster stack
r_pop_static <- rasterFromXYZ(pop_static_wide, crs=crs(CH_mask))
# Extend to original spatial extent
r_pop_static <- extend(r_pop_static, CH_mask)
# Change names
names(r_pop_static) <- sub('X','Year.',names(r_pop_static))
# Set empty cells within Switzerland to zero (instead of NA)
values(r_pop_static)[!is.na(values(CH_mask))&is.na(values(r_pop_static))] <- 0

# Map (while masking to Swiss borders)
plot(mask(r_pop_static, CH_mask))


## --------------------------------------------------------------------------------
# Run simulations with 25 years spin up
RunRS(s_cc, dirpath)


## ------------------------------------------------------------------------------------------------
# plot the resulting abundances and occupancy
par(mfrow=c(1,2))

# Abundances
plotAbundance(s_cc, dirpath)

# Occupancy
plotOccupancy(s_cc, dirpath)


## -----------------------------------------------------------------------------------
# Read population output file into a data frame
pop_dynamic <- readPop(s_cc, dirpath)

# Make data frame with mean abundances per cell and selected year:
pop_dynamic_wide <- pop_dynamic %>%
  # Select years you want to map: every ten years after spinup
  filter(Year %in% c(spinup, spinup+seq(10,sim_years,by=10)) ) %>%
  # Define grouping variables
  group_by(x,y,Year) %>%
  # Average individuals over replicate runs
  summarise(meanInd = mean(NInd), .groups='keep') %>%
  # Make separate columns for different years
  pivot_wider(names_from=Year, values_from=meanInd)


## ------------------------------------------------------------------------------------------------
# Make raster stacks

# Transfer into correct coordinate system
pop_dynamic_wide$x <- pop_dynamic_wide$x + xmin(CH_mask)
pop_dynamic_wide$y <- pop_dynamic_wide$y + ymin(CH_mask)

# Make raster stack
r_pop_dynamic <- rasterFromXYZ(pop_dynamic_wide, crs=crs(CH_mask))
# Extend to original spatial extent
r_pop_dynamic <- extend(r_pop_dynamic, CH_mask)
# Change names
names(r_pop_dynamic) <- sub('X','Year.',names(r_pop_dynamic))
# Set empty cells within Switzerland to zero (instead of NA)
values(r_pop_dynamic)[!is.na(values(CH_mask))&is.na(values(r_pop_dynamic))] <- 0

# Map
plot(mask(r_pop_dynamic, CH_mask))


## ------------------------------------------------------------------------------------------------
# We reduce juvenile survival
(trans_mat2 <- matrix(c(0, 0.6,1.646,   0.45), nrow = 2, byrow = F)) 

# Define stage structured module
stg2 <- StageStructure(Stages = 2, # 2 stages plus stage 0
                      TransMatrix = trans_mat2,
                      MaxAge = 10,
                      FecDensDep = F , 
                      SurvDensDep = T , # Only mortality increases with density
                      SurvDensCoeff = 1 # Density dependence relative to 1/b
)

# Define demography module
demo2 <- Demography(StageStruct = stg2, ReproductionType = 0) # asexual model

# Set new parameter master object for static landscape
s2 <- RSsim(batchnum = 12, land = land, demog = demo2, dispersal = disp, simul = sim, init = init, seed = 324135)


## --------------------------------------------------------------------------------------
# Run simulations
RunRS(s2, dirpath)


## ---------------------------------------------------------------------------------
# plot the resulting abundances and occupancy
par(mfrow=c(1,2))
plotAbundance(s2, dirpath)
plotOccupancy(s2, dirpath)


## ------------------------------------------------------------------------------------------------
# We reduce juvenile survival
(trans_mat3 <- matrix(c(0, 0.6,1.646,   0.55), nrow = 2, byrow = F)) 

# Define stage structured module
stg3 <- StageStructure(Stages = 2, # 2 stages plus stage 0
                      TransMatrix = trans_mat3,
                      MaxAge = 10,
                      FecDensDep = F , 
                      SurvDensDep = T , # Only mortality increases with density
                      SurvDensCoeff = 1 # Density dependence relative to 1/b
)

# Define demography module
demo3 <- Demography(StageStruct = stg3, ReproductionType = 0) # asexual model

# Set new parameter master object for static landscape
s3 <- RSsim(batchnum = 13, land = land, demog = demo3, dispersal = disp, simul = sim, init = init, seed = 324135)


## -------------------------------------------------------------------------------------
# Run simulations
RunRS(s3, dirpath)


## -------------------------------------------------------------------------------------
# plot the resulting abundances and occupancy
par(mfrow=c(1,2))
plotAbundance(s3, dirpath)
plotOccupancy(s3, dirpath)


## ------------------------------------------------------------------------------------------------
# Set static landscape parameters
land2 <- ImportedLandscape(LandscapeFile = "BlackGrouse_SDMpredictions2000.asc", 
        Resolution = 1000, # resolution of one grid cell in m
        HabPercent = TRUE, # habitat map includes percentage of suitability
        K_or_DensDep = 0.1) # strength of density dependence 1/b

# Set new parameter master object for static landscape
s4 <- RSsim(batchnum = 14, land = land2, demog = demo, dispersal = disp, simul = sim, init = init, seed = 324135)


## ----------------------------------------------------------------------------------
# Run simulations
RunRS(s4, dirpath)


## ---------------------------------------------------------------------------------
# plot the resulting abundances and occupancy
par(mfrow=c(1,2))
plotAbundance(s4, dirpath)
plotOccupancy(s4, dirpath)


## --------------------------------------------------------------------------

# Join mean (and sd) of abundances over all static scenarios in a single data frame
abund_sens <- bind_rows(
  # Scenario static 1:
  # The function readRange() runs in the background of plotAbundance(). Here, we extract abundances per scenario by hand.
  readRange(s,dirpath) %>%
    group_by(Year) %>%
    summarise(Abundance = mean(NInds), sd = sd(NInds)) %>% add_column(Scenario = "1 - Default"), 
  # Scenario static 2
  readRange(s2,dirpath) %>%
    group_by(Year) %>%
    summarise(Abundance = mean(NInds), sd = sd(NInds)) %>% add_column(Scenario = "2 - AdultSurv -5%"), 
  # Scenario static 3
  readRange(s3,dirpath) %>%
    group_by(Year) %>%
    summarise(Abundance = mean(NInds), sd = sd(NInds)) %>% add_column(Scenario = "3 - AdultSurv +5%"),
  # Scenario static 4
  readRange(s4,dirpath) %>%
    group_by(Year) %>%
    summarise(Abundance = mean(NInds), sd = sd(NInds)) %>% add_column(Scenario = "4 - TerritorySize * 5")) 

# Plot abundance
ggplot(data = abund_sens, mapping = aes(x = Year, y = log(Abundance),  color=Scenario)) + 
  geom_line(size=2) +
  geom_ribbon(aes(ymin=log(Abundance-sd), ymax=log(Abundance+sd)), linetype=2, alpha=0.1)

