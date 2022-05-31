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
# Define maximum body size of endotherms in the cohort definitions:
c_defs_B =c_defs_S = madingley_inputs('cohort definition')

c_defs_B$PROPERTY_Maximum.mass[1] = 500000 # Largest herbivore Bialowieza: Bison bonasus - bison (500 kg)
c_defs_B$PROPERTY_Maximum.mass[2] = 32000 # Largest carnivore Bialowieza: Canis lupus -wolve (32 kg)
c_defs_B$PROPERTY_Maximum.mass[3] = 181000 # Largest omnivore Bialowieza: Ursus arctos - brown bear (181 kg)

c_defs_S$PROPERTY_Maximum.mass[1] = 3940000 # Largest herbivore Serengeti: Loxodonta africana - elefant (3940 kg)
c_defs_S$PROPERTY_Maximum.mass[2] = 162000 # Largest carnivore Serengeti: Panthera leo - lion (162 kg)
c_defs_S$PROPERTY_Maximum.mass[3] = 83000 # Largest omnivore Serengeti: Phacochoerus africanus - warthog (83 kg)


# Run spin-up of 50 years for Bialowieza and for Serengeti
mres_spinup_B = madingley_run(madingley_data = mdat_B,
                    years = 50,
                    cohort_def = c_defs_B,
                    out_dir=dirpath)

mres_spinup_S = madingley_run(madingley_data = mdat_S,
                    years = 50,
                    cohort_def = c_defs_S,
                    out_dir=dirpath)


## ------------------------------------------------------------------------------------
# save model objects
save(mres_spinup_B, mres_spinup_S, file=paste0(dirpath,'/mres_carnivores.RData'))


## --------------------------------------------------------------------------
# Run control scenario for 50 years.
# Note that simulations start from spinup results:
mres_control_B = madingley_run(madingley_data = mres_spinup_B,
                      years = 50,
                      cohort_def = c_defs_B,
                      out_dir=dirpath)

mres_control_S = madingley_run(madingley_data = mres_spinup_S,
                      years = 50,
                      cohort_def = c_defs_B,
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

# Remove large carnivores >21 kg in Bialowieza
remove_idx_B = which(mres_removal_B$cohorts$AdultMass > 21000 &
                    mres_removal_B$cohorts$FunctionalGroupIndex == 1)
mres_removal_B$cohorts = mres_removal_B$cohorts[-remove_idx_B, ]

# Remove large carnivores >21 kg in Serengeti
remove_idx_S = which(mres_removal_S$cohorts$AdultMass > 21000 &
                    mres_removal_S$cohorts$FunctionalGroupIndex == 1)
mres_removal_S$cohorts = mres_removal_S$cohorts[-remove_idx_S, ]


# Restrict maximum body size of endotherm carnivores in the cohort definitions to <21kg:
c_defs_B_removal = c_defs_B
c_defs_B_removal$PROPERTY_Maximum.mass[2] = 21000
c_defs_S_removal = c_defs_S
c_defs_S_removal$PROPERTY_Maximum.mass[2] = 21000


# Run simulation:
# Bialowieza
mres_removal_B = madingley_run(madingley_data = mres_removal_B,
                      cohort_def = c_defs_B_removal,
                      years = 50,
                      out_dir=dirpath)

# Serengeti
mres_removal_S = madingley_run(madingley_data = mres_removal_S,
                      cohort_def = c_defs_S_removal,
                      years = 50,
                      out_dir=dirpath)


## --------------------------------------------------------------------------------------
# save model objects
save(mres_spinup_B, mres_spinup_S, mres_control_B, mres_control_S, mres_removal_B, mres_removal_S, file=paste0(dirpath,'/mres_carnivores.RData'))


## ---------------------------------------------------------------------------
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
par(mfrow=c(1,1))
plot_foodweb(mres_control_B, max_flows = 5)
title('Bialowieza - control scenario')
plot_foodweb(mres_removal_B, max_flows = 5)
title('Bialowieza - removal scenario')
plot_foodweb(mres_control_S, max_flows = 5)
title('Serengeti - control scenario')
plot_foodweb(mres_removal_S, max_flows = 5)
title('Serengeti - removal scenario')

# Plot body mass densities
par(mfcol=c(3,2))
hist(subset(mres_control_B$cohorts, FunctionalGroupIndex==0)$AdultMass/1000, main='Bialowieza control - herbivores', xlab='Body mass kg')
hist(subset(mres_control_B$cohorts, FunctionalGroupIndex==1)$AdultMass/1000, main='Bialowieza control - carnivores', xlab='Body mass kg')
hist(subset(mres_control_B$cohorts, FunctionalGroupIndex==2)$AdultMass/1000, main='Bialowieza control - omnivores', xlab='Body mass kg')
hist(subset(mres_removal_B$cohorts, FunctionalGroupIndex==0)$AdultMass/1000, main='Bialowieza removal - herbivores', xlab='Body mass kg')
hist(subset(mres_removal_B$cohorts, FunctionalGroupIndex==1)$AdultMass/1000, main='Bialowieza removal - carnivores', xlab='Body mass kg')
hist(subset(mres_removal_B$cohorts, FunctionalGroupIndex==2)$AdultMass/1000, main='Bialowieza removal - omnivores', xlab='Body mass kg')

par(mfcol=c(3,2))
hist(subset(mres_control_S$cohorts, FunctionalGroupIndex==0)$AdultMass/1000, main='Serengeti control - herbivores', xlab='Body mass kg')
hist(subset(mres_control_S$cohorts, FunctionalGroupIndex==1)$AdultMass/1000, main='Serengeti control - carnivores', xlab='Body mass kg')
hist(subset(mres_control_S$cohorts, FunctionalGroupIndex==2)$AdultMass/1000, main='Serengeti control - omnivores', xlab='Body mass kg')
hist(subset(mres_removal_S$cohorts, FunctionalGroupIndex==0)$AdultMass/1000, main='Serengeti removal - herbivores', xlab='Body mass kg')
hist(subset(mres_removal_S$cohorts, FunctionalGroupIndex==1)$AdultMass/1000, main='Serengeti removal - carnivores', xlab='Body mass kg')
hist(subset(mres_removal_S$cohorts, FunctionalGroupIndex==2)$AdultMass/1000, main='Serengeti removal - omnivores', xlab='Body mass kg')

