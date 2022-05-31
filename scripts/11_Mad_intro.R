# Session 11: Running MadingleyR


# library(devtools)
# 
# # Install the MadingleyR package
# install_github('MadingleyR/MadingleyR', subdir='Package', build_vignettes = TRUE)


## ------------------------------------------------------------------------------------------------
# Load MadingleyR package 
library('MadingleyR')

# Get version MadingleyR and C++ source code
madingley_version( )

## -------------------------------------------------------------------------------------
# View the MadingleyR tutorial vignette
vignette('MadingleyR')


## ------------------------------------------------------------------------------------------------
dirpath = paste0(getwd(),"/models/Mad_intro")


## -------------------------------------------------------------------------------
# Initialise MadingleyR data with default values
mdata = madingley_init()

# Structure of initialissed data
str(mdata,1)

# Initialised spatial window
plot_spatialwindow(mdata$spatial_window)


## ------------------------------------------------------------------------------------------------
# Check potential inputs
madingley_inputs()

# Load default inputs for spatial layers
sptl_inp = madingley_inputs("spatial inputs")

# Look at structure of spatial input layers
str(sptl_inp,1)


## ------------------------------------------------------------------------------------------------
# Inspect raster layer for endothermic omnivores
sptl_inp$Endo_O_max

# Restrict size to 150 kg (150'000 g)
sptl_inp$Endo_O_max[sptl_inp$Endo_O_max>150000] = 150000

# Plot initial body mass distribution of endothermic omnivores
plot(sptl_inp$Endo_O_max)

# Reinitialise the model
mdata = madingley_init(spatial_inputs = sptl_inp)


## -----------------------------------------------------------------------------
# Run spin-up of 50 years
mres = madingley_run(madingley_data = mdata,
                    years = 20,
                    max_cohort = 200,
                    out_dir=dirpath)

# save model object
save(mres,file=paste0(dirpath,'/mres.RData'))


## ------------------------------------------------------------------------------------------------
# Look at structure of output data sets
str(mres,1)


## -----------------------------------------------------------------------------------
# Plot biomass time series
plot_timelines(mres)

# Plot body mass distributions of heterotroph functional groups 
# (note that semelparous and iteroparous ectotherms are pooled together)
plot_densities(mres)

# Plot trophic pyramid
plot_trophicpyramid(mres)

# Plot foodweb plot
plot_foodweb(mres, max_flows = 5)

# Plot total abundance per grid cell and abundance distribution per functional group
plot_spatialabundances(mres)
plot_spatialabundances(mres, functional_filter = T)

# Plot total biomass per grid cell and biomass distribution per functional group
plot_spatialbiomass(mres)
plot_spatialbiomass(mres, functional_filter = T)


