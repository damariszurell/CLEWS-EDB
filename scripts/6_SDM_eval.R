# Session 6: SDM assessment and prediction

library(raster)

load('data/lynx_thinned.RData')


m_step <- step(
  glm( occ ~ bio11 + I(bio11^2) + bio10 + I(bio10^2), 
               family='binomial', data=lynx_thinned)
  )

summary(m_step)

# If we do not provide "newdata", then predict() should simply return the fitted values: 
head(predict(m_step, type='response'))
head(m_step$fitted)


## ------------------------------------------------------------------------------------------------
# We want to make predictions for all combinations of the two predictor variables
# and along their entire environmental gradients:
xyz <- expand.grid(
  # We produce a sequence of environmental values within the predictor ranges:
	bio11 = seq(min(lynx_thinned$bio11),max(lynx_thinned$bio11),length=50),
	bio10 = seq(min(lynx_thinned$bio10),max(lynx_thinned$bio10),length=50)
	)

# Now we can make predictions to this new data frame
xyz$z <- predict(m_step, newdata=xyz, type='response')
summary(xyz)

# As result, we have a 3D data structure and want to visualise this.
# Here, I first set a color palette
library(RColorBrewer)
cls <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)

# Finally, we plot the response surface using the wireframe function from the lattice package
library(lattice)
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1))

# We can also rotate the axes to better see the surface
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1), 
          screen=list(z = -120, x = -70, y = 3))


## ------------------------------------------------------------------------------------------------
library(mecofun)

# Names of our variables:
my_preds <- c('bio11', 'bio10')

# We want two panels next to each other:
par(mfrow=c(1,2))

# Plot the partial responses
partial_response(m_step, predictors = lynx_thinned[,my_preds])


expl_deviance(obs = lynx_thinned$occ,
              pred = m_step$fitted)


## ------------------------------------------------------------------------------------------------
preds_cv <- crossvalSDM(m_step, traindat = lynx_thinned, colname_species = 'occ', colname_pred = my_preds)


## ------------------------------------------------------------------------------------------------
plot(m_step$fitted.values, preds_cv, xlab='Fitted values', ylab='Predicted values from CV')
abline(0,1,col='red',lwd=2)


library(PresenceAbsence)

# We first prepare our data:
# Prepare cross-validated predictions:
thresh_dat <- data.frame(
  ID = seq_len(nrow(lynx_thinned)), 
	obs = lynx_thinned$occ,
	pred = preds_cv)
		
# Then, we find the optimal thresholds:		
(thresh_cv <- PresenceAbsence::optimal.thresholds(DATA= thresh_dat))


## ------------------------------------------------------------------------------------------------
(cmx_maxSSS <- PresenceAbsence::cmx(DATA= thresh_dat, threshold=thresh_cv[3,2]))


## ------------------------------------------------------------------------------------------------
# Proportion of correctly classified observations
PresenceAbsence::pcc(cmx_maxSSS, st.dev=F)

# Sensitivity = true positive rate
PresenceAbsence::sensitivity(cmx_maxSSS, st.dev=F)

# Specificity = true negative rate
PresenceAbsence::specificity(cmx_maxSSS, st.dev=F)


# Kappa
PresenceAbsence::Kappa(cmx_maxSSS, st.dev=F)

# True skill statistic
TSS(cmx_maxSSS)	


library(AUC)

# Let's have a look a the ROC curve:
roc_cv <- roc(preds_cv, as.factor(lynx_thinned$occ))
plot(roc_cv, col = "grey70", lwd = 2)

# Compute the AUC:
AUC::auc(roc_cv)


evalSDM(lynx_thinned$occ, preds_cv)


# First set your file path that you want to download climate data
# to or where you have stored climate data already.
my_filepath <- 'data'

library(raster)

# rough European extent
extent_eur <- c(-15,45,35,72)

# Please note that you have to set download=T if you haven't downloaded the data before:
bio_curr <- getData("worldclim", var="bio", res=10, download=F, path=my_filepath)

# Crop data to European extent
bio_curr <- crop(bio_curr, extent_eur)

# Aggregate to 1° resolution
bio_curr <- aggregate(bio_curr, 6)


# Please note that you have to set download=T if you haven't downloaded the data before:
bio_fut <- getData('CMIP5', var='bio', download=F, res=10, 
                   rcp=45, model='NO', year=50, path=my_filepath)

# Rename to hold the same bioclim names as the current climate layers
names(bio_fut) <- names(bio_curr)

# Crop data to European extent
bio_fut <- crop(bio_fut, extent_eur)

# Aggregate to 1° resolution
bio_fut <- aggregate(bio_fut, 6)


## ------------------------------------------------------------------------------------------------
# Prepare data frames
bio_curr_df <- data.frame(rasterToPoints(bio_curr))
bio_fut_df <- data.frame(rasterToPoints(bio_fut))

# Make continuous predictions:
bio_curr_df$pred_glm <- predict(m_step, newdata= bio_curr_df, type="response")
bio_fut_df$pred_glm <- predict(m_step, newdata= bio_fut_df, type="response")


## ------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))

# Make raster of predictions to current environment:
r_pred_curr <- rasterFromXYZ(bio_curr_df[,c('x','y','pred_glm')])
plot(r_pred_curr, axes=F, main='Occ. prob. - today')

# Make raster stack of predictions to future environment:
r_pred_fut <- rasterFromXYZ(bio_fut_df[,c('x','y','pred_glm')])
plot(r_pred_fut, axes=F, main='Occ. prob. - 2050')


## ------------------------------------------------------------------------------------------------
# Make binary predictions:
bio_curr_df$bin_glm <- ifelse(bio_curr_df$pred_glm >= thresh_cv[3,2], 1, 0)
bio_fut_df$bin_glm <- ifelse(bio_fut_df$pred_glm >= thresh_cv[3,2], 1, 0)

# Make raster stack of predictions to current environment:
r_pred_curr <- rasterFromXYZ(bio_curr_df[,c('x','y','pred_glm','bin_glm')])
plot(r_pred_curr, axes=F)

# Make raster stack of predictions to future environment:
r_pred_fut <- rasterFromXYZ(bio_fut_df[,c('x','y','pred_glm','bin_glm')])
plot(r_pred_fut, axes=F)

library(dismo)

# MESS maps from the dismo package:
r.mess <- mess(bio_fut[[my_preds]], lynx_thinned[,my_preds])
plot(r.mess, axes=F)

# Negative values indicate dissimilar=novel environments:
r.mess.mask <- r.mess<0
r.mess.mask <- mask(r.mess.mask,bg)
plot(r.mess.mask, axes=F)

# Predictions to analogous climates:
r_analog_fut <- r_pred_fut
values(r_analog_fut)[values(r.mess)<0] <- NA
plot(r_analog_fut, axes=F)

# Predictions to novel climates:
r_novel_fut <- r_pred_fut
values(r_novel_fut)[values(r.mess)>=0] <- NA
plot(r_novel_fut, axes=F)

