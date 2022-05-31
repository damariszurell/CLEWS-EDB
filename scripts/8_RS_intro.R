# Session 8: Getting starter with RangeShiftR

## ------------------------------------------------------------------------------------------------
library(RangeShiftR)


## ------------------------------------------------------------------------------------------------
s <- RSsim()


dirpath = "models/RS_GettingStarted/"


# Create sub-folders (if not already existing)
if(!file.exists(paste0(dirpath,"Inputs"))) {
  dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE) }
if(!file.exists(paste0(dirpath,"Outputs"))) {
  dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE) }
if(!file.exists(paste0(dirpath,"Output_Maps"))) {
  dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE) }


RunRS(s,dirpath)


## ------------------------------------------------------------------------------------------------
s

## ------------------------------------------------------------------------------------------------
sim <- Simulation(Simulation = 2,
                  Years = 50,
                  Replicates = 2,
                  OutIntPop = 50)


?Simulation


# land <- ImportedLandscape()
# land <- ArtificialLandscape()


## ------------------------------------------------------------------------------------------------
land <- ArtificialLandscape(Resolution = 10,  # in meters
                            K_or_DensDep = 1500,  # ~ 15 inds/cell
                            propSuit = 0.2,
                            dimX = 129, dimY = 257, 
                            fractal = T, hurst = 0.3,
                            continuous = F)


## ------------------------------------------------------------------------------------------------
demo <- Demography(Rmax = 2.2, ReproductionType = 1, PropMales = 0.45)


## ------------------------------------------------------------------------------------------------
stg <- StageStructure(Stages = 3,
                      TransMatrix = matrix(c(0,1,0,5.7,.5,.4,3.4,0,.9),nrow = 3),
                      FecDensDep = T,
                      SurvDensDep = T)
demo <- demo + stg


## ------------------------------------------------------------------------------------------------
demo <- Demography(StageStruct = stg, ReproductionType = 1, PropMales = 0.45)


## ------------------------------------------------------------------------------------------------
plotProbs(stg)


## ------------------------------------------------------------------------------------------------
disp <-  Dispersal(Emigration = Emigration(EmigProb = 0.2), 
                   Transfer   = DispersalKernel(Distances = 50),
                   Settlement = Settlement() )


## ------------------------------------------------------------------------------------------------
plotProbs(DispersalKernel(Distances = matrix(c(0,1,2,70,50,30),nrow = 3), StageDep = T))


## ------------------------------------------------------------------------------------------------
disp <-  disp + Settlement(SexDep = T, 
                           Settle = matrix(c(0,1,1,0), nrow = 2))


## ------------------------------------------------------------------------------------------------
init <- Initialise(FreeType = 0, 
                   NrCells = 2250,
                   InitDens = 2, 
                   IndsHaCell = 3, 
                   PropStages = c(0,0.7,0.3))
init


## ------------------------------------------------------------------------------------------------
s <- RSsim(simul = sim, land = land, demog = demo, dispersal = disp, init = init)

# Alternative notation:
# s <- RSsim() + land + demo + disp + sim + init + gene


## ------------------------------------------------------------------------------------------------
validateRSparams(s)


RunRS(s, dirpath)


## ---------------------------------------------------------------------------------------
range_df <- readRange(s, dirpath)

# ...with replicates:
par(mfrow=c(1,2))
plotAbundance(range_df)
plotOccupancy(range_df)

# ...with standard deviation:
par(mfrow=c(1,2))
plotAbundance(range_df, sd=T, replicates = F)
plotOccupancy(range_df, sd=T, replicates = F)


## ---E-------------------------------------------------------------------
# read population output file into a dataframe
pop_df <- readPop(s, dirpath)

# Make data frame with number of individuals per output year - for only one repetition (Rep==0):
pop_wide_rep0 <- reshape(subset(pop_df,Rep==0)[,c('Year','x','y','NInd')], timevar='Year', v.names=c('NInd'), idvar=c('x','y'), direction='wide')

head(pop_wide_rep0)

# use raster package to make a raster from the data frame
library(raster)
stack_years_rep0 <- rasterFromXYZ(pop_wide_rep0)
names(stack_years_rep0) <- c('Year.0', 'Year.50')
spplot(stack_years_rep0, zlim = c(0,10))

