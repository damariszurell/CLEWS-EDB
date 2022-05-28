# Session 9: Eurasian lynx reintroduction in RangeShiftR

## ----------------------------------------------------------------------------------------
# Set path to model directory
dirpath = "models/RS_Lynx/"


## ---------------------------------------------------------------------
# load required packages
library(RangeShiftR)
library(raster)
library(rasterVis)
library(viridis)
library(RColorBrewer)
library(latticeExtra)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggplot2) 

# Create sub-folders (if not already existing)
if(!file.exists(paste0(dirpath,"Inputs"))) {
  dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE) }
if(!file.exists(paste0(dirpath,"Outputs"))) {
  dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE) }
if(!file.exists(paste0(dirpath,"Output_Maps"))) {
  dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE) }


## --------------------------------------------------------------------------
# Land cover map of scotland
landsc <- raster(paste0(dirpath, "Inputs/LCM_Scotland_2015_1000.asc"))
plot(landsc)

# Plot land cover map and highlight cells with initial species distribution - option 2 with categorical legend:
landsc.f <- as.factor(landsc)

# add the land cover classes to the raster attribute table
(rat <- levels(landsc.f)[[1]])
rat[["land-use"]] <- c("salt water", "arable + horticulture", "freshwater", "built-up areas + gardens", "inland rock", "grasslands", "woodland", "supra-/littoral sediment", "marsh, swamp", "heath")
levels(landsc.f) <- rat

levelplot(landsc.f, margin=F, scales=list(draw=FALSE), col.regions=brewer.pal(n = 10, name = "Spectral"))


## --------------------------------------------------------------------------
# Read in woodland patch files
patches <- raster(paste0(dirpath,'Inputs/woodland_patchIDs_1000.asc'))

# Plot the patches in different colours:
plot(patches,axes=F, legend=F, col = c('grey',rep(brewer.pal(n = 10, name = "Paired"),4)))

# In total, we have 39 patches of different size
values(patches) %>% table


## -------------------------------------------------------------------------------

# Plot patches and add labels
plot(patches,axes=F, legend=F, col = c('grey',rep(brewer.pal(n = 10, name = "Paired"),4)))
patches_pol <- rasterToPolygons(patches, dissolve=T) # Makes spatial polygons
text(patches_pol,labels=patches_pol@data$woodland_patchIDs_1000)


## ------------------------------------------------------------------------------
land <- ImportedLandscape(LandscapeFile = "LCM_Scotland_2015_1000.asc",
                          PatchFile = "woodland_patchIDs_1000.asc", 
                          Resolution = 1000,
                          Nhabitats = 10,
                          K_or_DensDep = c(0, 0, 0, 0, 0, 0, 0.000285, 0, 0, 0)
                          )


## ------------------------------------------------------------------------------------------------
# Transition matrix
(trans_mat <- matrix(c(0,1,0,0,0,0, 0.53, 0,0, 0, 0, 0.63, 5,0, 0, 0.8), nrow = 4, byrow = F))  # stage 0: dispersing newborns; stage 1 : first-year juveniles; stage 2: non-breeding sub-adults; stage 3: breeding adults

# Define the stage structure
stg <- StageStructure(Stages = 4, # four stages including dispersing juveniles
                      TransMatrix = trans_mat,  # transition matrix
                      MaxAge = 17,  # maximum age
                      FecDensDep = T   # density dependence in fecundity
)

# Plot the vital rates against different density levels
plotProbs(stg)

# Define the demography module
demo <- Demography(StageStruct = stg, 
                   ReproductionType = 1, # simple sexual model
                   PropMales = 0.5) 


## ------------------------------------------------------------------------------------------------
emig <- Emigration(DensDep=T, StageDep=T, SexDep = T,
                   EmigProb = cbind(c(0,0,1,1,2,2,3,3),
                                    c(0,1,0,1,0,1,0,1),
                                    c(0.4, 0.9, 0, 0, 0, 0, 0, 0),
                                    c(10, 10, 0, 0, 0, 0, 0, 0), 
                                    c(1, 1, 0, 0, 0, 0, 0, 0)) ) # only emigration of juveniles, females higher than males


## ------------------------------------------------------------------------------------------------
transfer <- SMS(PR = 1, # Perceptual range in number of cells
                PRMethod = 2, # Harmonic mean used to quantify movement cost within perceptual range
                MemSize = 5, # number of steps remembered when applying directional persistence.
                DP = 5,  # directional persistence.
                Costs = c(100000, 30, 100, 1000, 1000, 10, 1, 10, 10, 7), # movement cost per habitat class
                StepMort = c(0.9999, 0.0002, 0.0005, 0.007, 0.00001, 0.00001, 0, 0.00001, 0.00001, 0.00001)) # per step mortality per habitat class


## ------------------------------------------------------------------------------------------------
settle <- Settlement(StageDep = F,
                     SexDep = T,
                     Settle = cbind(c(0, 1), c(1.0, 1.0), c(-10, -10), c(1, 1)), # here no difference between sexes
                     FindMate = c(F, T), # only males need to find a female
                     DensDep = T,
                     MaxSteps = 500
                     )


## ------------------------------------------------------------------------------------------------
# Dispersal
disp <-  Dispersal(Emigration = emig,
                   Transfer = transfer,
                   Settlement = settle) 

# plot parameters
plotProbs(disp@Emigration)




## ------------------------------------------------------------------------------------------------
# prepare dataframe for InitIndsFile
(init_df_29 <- data.frame(Year=0,Species=0,PatchID=29,Ninds=c(5,5),Sex=c(0,1),Age=3,Stage=3))


## ------------------------------------------------------------------------------------------------
# write InitIndsFile to file
write.table(init_df_29, file=paste0(dirpath,'Inputs/InitInds_29.txt'), sep='\t', row.names=F, quote=F)

# Set initialisation
init_29 <- Initialise(InitType = 2,       # Initialise from initial individuals list file
                      InitIndsFile = 'InitInds_29.txt')



## ------------------------------------------------------------------------------------------------
# prepare dataframe for InitIndsFile
(init_df_29_15 <- data.frame(Year=0,Species=0,PatchID=rep(c(29,15),each=2),Ninds=rep(c(9,7),each=2),Sex=c(0,1),Age=3,Stage=3))

# write InitIndsFile to fil
write.table(init_df_29_15, file=paste0(dirpath,'Inputs/InitInds_29_15.txt'), sep='\t', row.names=F, quote=F)

# Set initialisation
init_29_15 <- Initialise(InitType = 2,       # from loaded species distribution map
                      InitIndsFile = 'InitInds_29_15.txt')



## ------------------------------------------------------------------------------------------------
# Simulations
RepNb <- 100

sim <- Simulation(Simulation = 0, # ID
                  Replicates = RepNb, # number of replicate runs
                  Years = 100, # number of simulated years
                  OutIntPop = 1, # output interval
                  OutIntOcc = 1,
                  OutIntRange = 1)



## ------------------------------------------------------------------------------------------------
# RangeShifter parameter master object for single-site reintroduction
s_29 <- RSsim(batchnum = 1, land = land, demog = demo, dispersal = disp, simul = sim, init = init_29, seed = 324135)

# RangeShifter parameter master object for multi-site reintroduction
s_29_15 <- RSsim(batchnum = 3, land = land, demog = demo, dispersal = disp, simul = sim, init = init_29_15, seed = 324135)


## --------------------------------------------------------------------------------------
# Run simulations
RunRS(s_29, dirpath)


## ------------------------------------------------------------------------------------------------
# general output of population size + occupancy
par(mfrow=c(1,2))
plotAbundance(s_29,dirpath,sd=T, rep=F)
plotOccupancy(s_29, dirpath, sd=T, rep=F)



## ------------------------------------------------------------------------------------------------
# Colonisation metrics
col_stats_29 <- ColonisationStats(s_29, dirpath, years = 100, maps = T)


## ------------------------------------------------------------------------------------------------
# mean occupancy probability in year 100
head(col_stats_29$occ_prob)

# map occupancy probability
mycol_occprob <- colorRampPalette(c('blue','orangered','gold'))
levelplot(col_stats_29$map_occ_prob, margin=F, scales=list(draw=FALSE), at=seq(0,1,length=11), col.regions=mycol_occprob(11))


## ------------------------------------------------------------------------------------------------
# Some useful functions

# Store underlying landscape map display for later:
# Note: Not sure if this is really suitable, as we have 10 different landscape classes
bg <- function(main=NULL){
  levelplot(landsc, margin=F, scales=list(draw=FALSE), colorkey=F,
            col.regions = rev(colorRampPalette(brewer.pal(3, "Greys"))(10)), main=main)
}


# map occupancy probability on landscape background. For this, we first define a colorkey function
col.key <- function(mycol, at, space='top',pos=0.95, height=0.6, width=1) {
  key <- draw.colorkey(
    list(space=space, at=at, height=height, width=width,
         col=mycol)
  )
  key$framevp$y <- unit(pos, "npc")
  return(key)
}


## ------------------------------------------------------------------------------------------------
# map occupancy probability after 100 years
bg() + levelplot(col_stats_29$map_occ_prob, margin=F, scales=list(draw=FALSE), at=seq(0,1,length=11), col.regions=mycol_occprob(11))
grid.draw(col.key(mycol_occprob(11),at=seq(0,1,length=11)))

# map colonisation time
mycol_coltime <- colorRampPalette(c('orangered','gold','yellow','PowderBlue','LightSeaGreen'))
levelplot(col_stats_29$map_col_time, margin=F, scales=list(draw=FALSE), at=c(-9,seq(-.001,100,length=11)), col.regions=c('blue',mycol_coltime(11)))

# map colonisation time on landscape background
bg() + levelplot(col_stats_29$map_col_time, margin=F, scales=list(draw=FALSE), at=c(-9,seq(-.001,100,length=11)), col.regions=c('blue',mycol_coltime(11)))
grid.draw(col.key(c('blue',mycol_coltime(11)), c(-9,seq(-.001,100,length=11))))


## -----------------------------------------------------------------------------------
# read population output file into a data frame
pop_29 <- readPop(s_29, dirpath)


# Calculate survival probability as number of replicate with surviving individuals per year
# Extinction probability is 1- survival probability:

# Calculate  extinction probability by hand:
pop_29 %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate
  summarise(sumPop = sum(NInd), .groups='keep') %>%
  group_by(Year) %>%
  # Average extinction probability (1 minus the proportion of replicates with surviving populations)
  summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb)



## ------------------------------------------------------------------------------------------------
# Mean time to extinction:
pop_29 %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate    
  summarise(sumPop = sum(NInd), .groups='keep') %>% 
  # Identify in which year they go extinct
  filter(sumPop==0) %>% 
  pull(Year) %>% mean



## ------------------------------------------------------------------------------------------------
# Define a function for calculating extinction probability
Calc_ExtProb <- function(pop_df,s) {
  require(dplyr)
  
  pop_df %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate
  summarise(sumPop = sum(NInd), .groups='keep') %>%
  group_by(Year) %>%
  # Average extinction probability (1 minus the proportion of replicates with surviving populations)
  summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb) %>%
  # Make sure that data frame is filled until last year of simulation
  right_join(tibble(Year = seq_len(s@simul@Years)), by='Year') %>% replace_na(list(extProb=1))
}

# Define a function for calculating mean time to extinction
Calc_ExtTime <- function(pop_df) {
  require(dplyr)
  
  pop_df %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate    
  summarise(sumPop = sum(NInd), .groups='keep') %>% 
  # Identify in which year they go extinct
  filter(sumPop==0) %>% 
  pull(Year) %>% mean
}


## ------------------------------------------------------------------------------------------------
# extinction probability
extProb_29 <- Calc_ExtProb(pop_29,s_29)

# Plot extinction probabilities
ggplot(data = extProb_29, mapping = aes(x = Year, y = extProb)) + 
  geom_line() +
  ylim(0,1)
  
# mean time to extinction
Calc_ExtTime(pop_29)


## ------------------------------------------------------------------------------------
# Run simulations
RunRS(s_29_15, dirpath)


## ------------------------------------------------------------------------------------------------
# general output of population size + occupancy
par(mfrow=c(1,2))
plotAbundance(s_29_15,dirpath,sd=T, rep=F)
plotOccupancy(s_29_15, dirpath, sd=T, rep=F)



## ------------------------------------------------------------------------------------------------
# Colonisation metrics
col_stats_29_15 <- ColonisationStats(s_29_15, dirpath, years = 100, maps = T)



## ------------------------------------------------------------------------------------------------
# map occupancy probability after 100 years
bg() + levelplot(col_stats_29_15$map_occ_prob, margin=F, scales=list(draw=FALSE), at=seq(0,1,length=11), col.regions=mycol_occprob(11))
grid.draw(col.key(mycol_occprob(11),at=seq(0,1,length=11)))

# map colonisation time 
bg() + levelplot(col_stats_29_15$map_col_time, margin=F, scales=list(draw=FALSE), at=c(-9,seq(-.001,100,length=11)), col.regions=c('blue',mycol_coltime(11)))
grid.draw(col.key(c('blue',mycol_coltime(11)), c(-9,seq(-.001,100,length=11))))


## ------------------------------------------------------------------------------------------------
# read population output file into a data frame
pop_29_15 <- readPop(s_29_15, dirpath)

# extinction probability
extProb_29_15 <- Calc_ExtProb(pop_29_15,s_29_15)

# mean time to extinction
Calc_ExtTime(pop_29_15)

# Compare extinction probabilities for single- and multi-site reintroduction

# Join extinction probabilities in a single data frame
extProb_scens <- bind_rows(extProb_29 %>% add_column(Scenario = "10 ind. Kintyre "),
                    extProb_29_15 %>% add_column(Scenario = "32 ind. Kintyre + Aberdeenshire"))

ggplot(data = extProb_scens, mapping = aes(x = Year, y = extProb, color=Scenario)) + 
  geom_line(size=2) +
  ylim(0,1)

