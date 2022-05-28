# Session 12: MadingleyR, role of large carnivores

## ------------------------------------------------------------------------------------------------
dirpath = paste0(getwd(),"/models/Mad_carnivores")


## --------------------------------------------------------------------------------
library(MadingleyR)

# Spatial window Bialowieza:
sptl_bialowieza = c(20, 25, 50, 54) 

# Spatial window Serengeti:
sptl_serengeti = c(32, 36, -4, 0)

# Initialise models for the two locations
mdat_B = madingley_init(spatial_window = sptl_bialowieza)
mdat_S = madingley_init(spatial_window = sptl_serengeti)



## ----------------------------------------------------------------------------
# Run spin-up of 50 years for Bialowieza and for Serengeti
mres_spinup_B = madingley_run(madingley_data = mdat_B,
                      years = 50,
                    out_dir=dirpath)

mres_spinup_S = madingley_run(madingley_data = mdat_S,
                      years = 50,
                    out_dir=dirpath)


## ------------------------------------------------------------------------------------
# save model objects
save(mres_spinup_B, mres_spinup_S, file=paste0(dirpath,'/mres_carnivores.RData'))


## --------------------------------------------------------------------------
# Run control scenario for 50 years.
# Note that simulations start from spinup results:
mres_control_B = madingley_run(madingley_data = mres_spinup_B,
                      years = 50,
                      out_dir=dirpath)

mres_control_S = madingley_run(madingley_data = mres_spinup_S,
                      years = 50,
                      out_dir=dirpath)


## -------------------------------------------------------------------------------------
# save model objects
save(mres_spinup_B, mres_spinup_S, mres_control_B, mres_control_S, file=paste0(dirpath,'/mres_carnivores.RData'))


## -----------------------------------------------------------------------------------
# Inspect adult body mass of endothermic carnivores in Bialowieza
summary(subset(mres_spinup_B$cohorts, FunctionalGroupIndex==1)$AdultMass)

# Inspect adult body mass of endothermic carnivores in Serengeti
summary(subset(mres_spinup_S$cohorts, FunctionalGroupIndex==1)$AdultMass)


## ----------------------------------------------------------------------------
# Define new output
mres_removal_B = mres_spinup_B
mres_removal_S = mres_spinup_S

# Remove large carnivores in Bialowieza
remove_idx_B = which(mres_removal_B$cohorts$AdultMass > 21000 &
                    mres_removal_B$cohorts$FunctionalGroupIndex == 1)
mres_removal_B$cohorts = mres_removal_B$cohorts[-remove_idx_B, ]

# Remove large carnivores in Serengeti
remove_idx_S = which(mres_removal_S$cohorts$AdultMass > 21000 &
                    mres_removal_S$cohorts$FunctionalGroupIndex == 1)
mres_removal_S$cohorts = mres_removal_S$cohorts[-remove_idx_S, ]

# Define maximum body size of endotherm carnivores in spatial input layers:
sptl_inp = madingley_inputs('spatial inputs')
sptl_inp$Endo_C_max[ ] = 21000

# Run simulation:
# Bialowieza
mres_removal_B = madingley_run(madingley_data = mres_removal_B,
                      spatial_inputs = sptl_inp,
                      years = 50,
                      out_dir=dirpath)

# Serengeti
mres_removal_S = madingley_run(madingley_data = mres_removal_S,
                      spatial_inputs = sptl_inp,
                      years = 50,
                      out_dir=dirpath)


## --------------------------------------------------------------------------------------
## # save model objects
## save(mres_spinup_B, mres_spinup_S, mres_control_B, mres_control_S, mres_removal_B, mres_removal_S, file=paste0(dirpath,'/mres_carnivores.RData'))


## ----warning=F, message=F------------------------------------------------------------------------
par(mfrow=c(2,2))

# Plot time series for Carnivores
dens_change_carn_B <- mres_removal_B$time_line_cohorts$Biomass_FG_1 / mres_control_B$time_line_cohorts$Biomass_FG_1
dens_change_carn_S <- mres_removal_S$time_line_cohorts$Biomass_FG_1 / mres_control_S$time_line_cohorts$Biomass_FG_1

plot(dens_change_carn_B,type='l', col="cyan4",lwd=2, xlab='Months', ylab='Change in density',ylim=c(0,2), main ='Carnivores')
abline(1,0, lty='dotted', col='grey')
lines(dens_change_carn_S,col='khaki',lwd=2)
legend('bottom',ncol=2,col=c('cyan4','khaki'),lwd=2,legend=c('Bialowieza','Serengeti'), bty='n')

# Plot time series for Herbivores
dens_change_herb_B <- mres_removal_B$time_line_cohorts$Biomass_FG_0 / mres_control_B$time_line_cohorts$Biomass_FG_0
dens_change_herb_S <- mres_removal_S$time_line_cohorts$Biomass_FG_0 / mres_control_S$time_line_cohorts$Biomass_FG_0

plot(dens_change_herb_B,type='l', col="cyan4",lwd=2, xlab='Months', ylab='Change in density',ylim=c(0.5,5), main='Herbivores')
abline(1,0, lty='dotted', col='grey')
lines(dens_change_herb_S,col='khaki',lwd=2)

# Plot time series for Omnivores
dens_change_omn_B <- mres_removal_B$time_line_cohorts$Biomass_FG_2 / mres_control_B$time_line_cohorts$Biomass_FG_2
dens_change_omn_S <- mres_removal_S$time_line_cohorts$Biomass_FG_2 / mres_control_S$time_line_cohorts$Biomass_FG_2

plot(dens_change_omn_B,type='l', col="cyan4",lwd=2, xlab='Months', ylab='Change in density',ylim=c(0.5,3), main='Omnivores')
abline(1,0, lty='dotted', col='grey')
lines(dens_change_omn_S,col='khaki',lwd=2)

# Plot time series for Autotrophs
dens_change_aut_B <- mres_removal_B$time_line_stocks$TotalStockBiomass / mres_control_B$time_line_stocks$TotalStockBiomass
dens_change_aut_S <- mres_removal_S$time_line_stocks$TotalStockBiomass / mres_control_S$time_line_stocks$TotalStockBiomass

plot(dens_change_aut_B,type='l', col="cyan4",lwd=2, xlab='Months', ylab='Change in density',ylim=c(0,1.5), main='Autotrophs')
abline(1,0, lty='dotted', col='grey')
lines(dens_change_aut_S,col='khaki',lwd=2)



## ------------------------------------------------------------------------------------------------
# Plot foodweb plot
plot_foodweb(mres_control_B, max_flows = 5)
title('Bialowieza - control scenario')
plot_foodweb(mres_removal_B, max_flows = 5)
title('Bialowieza - removal scenario')
plot_foodweb(mres_control_S, max_flows = 5)
title('Serengeti - control scenario')
plot_foodweb(mres_removal_S, max_flows = 5)
title('Serengeti - removal scenario')

# Plot trophic pyramid
plot_trophicpyramid(mres_control_B)
title('Bialowieza - control scenario')
plot_trophicpyramid(mres_removal_B)
title('Bialowieza - removal scenario')
plot_trophicpyramid(mres_control_S)
title('Serengeti - control scenario')
plot_trophicpyramid(mres_removal_S)
title('Serengeti - removal scenario')

