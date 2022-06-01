# Session 7: SDM algorithms and ensembles

library(raster)
load('data/lynx_thinned.RData')

my_preds <- c('bio11', 'bio10')


## ----eval=T--------------------------------------------------------------------------------------
set.seed(54321)

# First, we randomly select 70% of the rows that will be used as training data
train_i <- sample(seq_len(nrow(lynx_thinned)), size=round(0.7*nrow(lynx_thinned)))

# Then, we can subset the training and testing data
lynx_train <- lynx_thinned[train_i,]
lynx_test <- lynx_thinned[-train_i,]

# We store the split information for later:
write(train_i, file='data/indices_traindata.txt')


bio_curr <- getData("worldclim", var="bio", res=10, download=F, path='data')


## ------------------------------------------------------------------------------------------------
# rough European extent
extent_eur <- c(-15,45,35,72)

# Crop data to European extent
bio_curr <- crop(bio_curr, extent_eur)

# Aggregate to 1° resolution
bio_curr <- aggregate(bio_curr, 6)


## ----message=F, warning=F------------------------------------------------------------------------
library(dismo)

# Fit BIOCLIM model
m_bc <- bioclim(bio_curr[[my_preds]], # provide environmental raster
                lynx_train[lynx_train$occ==1,1:2]) # provide coordinates of presence locations
plot(m_bc)

# by hand:
lynx_env <- extract(bio_curr[[my_preds]],lynx_train[lynx_train$occ==1,1:2])
plot(lynx_env)
polygon(quantile(lynx_env[,1],c(.05,.95,.95,.05)), quantile(lynx_env[,2],c(.05,.05,.95,.95)))
polygon(quantile(lynx_env[,1],c(.025,.975,.975,.025)), quantile(lynx_env[,2],c(.025,.025,.975,.975)),lty='dotted')

# For the response surface, we first prepare the 3D-grid with environmental gradient and predictions
xyz <- expand.grid(
	seq(min(lynx_train[,my_preds[1]]),max(lynx_train[,my_preds[1]]),length=50),
	seq(min(lynx_train[,my_preds[2]]),max(lynx_train[,my_preds[2]]),length=50))
names(xyz) <- my_preds
# Make predictions to gradients:
xyz$z <- predict(m_bc, xyz)

# Define colour palette:
library(RColorBrewer)
cls <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)

# Plot response surface:
library(lattice)
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='BIOCLIM', xlab='bio11', ylab='bio10', 
          screen=list(z = -120, x = -70, y = 3))


library(mecofun)

# Plot partial response curves:
par(mfrow=c(1,2)) 
partial_response(m_bc, predictors = lynx_train[,my_preds], main='BIOCLIM')


# We use the default MaxSens+Spec threshold:
(perf_bc <- evalSDM(lynx_test$occ, predict(m_bc, lynx_test[,my_preds])))


## ------------------------------------------------------------------------------------------------
# Map predictions:
r_bc_pred <- r_bc_bin <- predict(m_bc,bio_curr)

# Threshold predictions using the maxTSS threshold (max sens+spec)
values(r_bc_bin) <- ifelse(values(r_bc_pred)>=perf_bc$thresh, 1, 0)

# plot the maps
plot(stack(r_bc_pred, r_bc_bin),main=c('BIOCLIM prob.','BIOCLIM bin.'), axes=F)	


# Fit GLM
m_glm <- step(glm( occ ~ bio11 + I(bio11^2) + bio10 + I(bio10^2),
	family='binomial', data=lynx_train))

# Now, we plot the response surface:
xyz$z <- predict(m_glm, xyz, type='response')

wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='GLM', xlab='bio11', ylab='bio10', 
          screen=list(z = -120, x = -70, y = 3))
	
# Plot partial response curves:
par(mfrow=c(1,2)) 
partial_response(m_glm, predictors = lynx_train[,my_preds], main='GLM')

# Performance measures
(perf_glm <- evalSDM(lynx_test$occ, predict(m_glm, lynx_test[,my_preds], type='response') ))

# Map predictions:
bio_curr_df <- data.frame(rasterToPoints(bio_curr[[my_preds]]))
r_glm_bin <- r_glm_pred <- rasterFromXYZ(cbind(bio_curr_df[,1:2],
                                  predict(m_glm, bio_curr_df, type='response')))
values(r_glm_bin) <- ifelse(values(r_glm_pred)>=perf_glm$thresh, 1, 0)
plot(stack(r_glm_pred, r_glm_bin),main=c('GLM prob.','GLM bin.'), axes=F)	


library(randomForest)

# Fit RF
m_rf <- randomForest( x=lynx_train[,my_preds], y=lynx_train$occ, 
	ntree=1000, importance =T)
	
# Variable importance:
importance(m_rf,type=1)
varImpPlot(m_rf)
	
# Look at single trees:
head(getTree(m_rf,1,T))

# Now, we plot the response surface:
xyz$z <- predict(m_rf, xyz, type='response')
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='Random Forest', xlab='bio11', ylab='bio10', 
          screen=list(z = -120, x = -70, y = 3))

	
# Plot partial response curves:
par(mfrow=c(1,2)) 
partial_response(m_rf, predictors = lynx_train[,my_preds], main='Random Forest')

# Performance measures of RF
(perf_rf <- evalSDM(lynx_test$occ, predict(m_rf, lynx_test[,my_preds],  type='response') ))

# Map predictions:
r_rf_bin <- r_rf_pred <- rasterFromXYZ(cbind(bio_curr_df[,1:2],
                                 predict(m_rf, bio_curr_df,type='response')))
values(r_rf_bin) <- ifelse(values(r_rf_pred)>=perf_rf$thresh, 1, 0)
plot(stack(r_rf_pred, r_rf_bin),main=c('RF prob.','RF bin.'), axes=F)	


library(gbm)

# Fit BRT
m_brt <- gbm.step(data = lynx_train, 
	gbm.x = my_preds,
	gbm.y = 'occ', 
	family = 'bernoulli',
	tree.complexity = 2,
	bag.fraction = 0.75,
	learning.rate = 0.001,
	verbose=F)
	
# Variable importance:
m_brt$contributions

# Interactions (not very meaningful here with only 2 predictors):
gbm.interactions(m_brt)$interactions
gbm.interactions(m_brt)$rank.list

# Now, we plot the response surface:
xyz$z <- predict.gbm(m_brt, xyz, n.trees=m_brt$gbm.call$best.trees, type="response")
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='Boosted regression trees', xlab='bio11', 
          ylab='bio10', screen=list(z = -120, x = -70, y = 3))
	
# Plot partial response curves:
par(mfrow=c(1,2)) 
partial_response(m_brt, predictors = lynx_train[,my_preds], main='BRT')

# Performance measures of BRT
(perf_brt <- evalSDM(lynx_test$occ, predict.gbm(m_brt, lynx_test[,my_preds], n.trees=m_brt$gbm.call$best.trees, type='response') ))

# Map predictions:
r_brt_bin <- r_brt_pred <- rasterFromXYZ(cbind(bio_curr_df[,1:2],
                                  predict.gbm(m_brt, bio_curr_df,
                                              n.trees=m_brt$gbm.call$best.trees, 
                                              type="response")))
values(r_brt_bin) <- ifelse(values(r_brt_pred)>=perf_brt$thresh, 1, 0)
plot(stack(r_brt_pred, r_brt_bin),main=c('BRT prob.','BRT bin.'), axes=F)	


(comp_perf <- rbind(bc = perf_bc, glm = perf_glm, rf = perf_rf, brt = perf_brt))

# We add a column containing the names of the algorithm
comp_perf <- data.frame(alg=row.names(comp_perf),comp_perf)

# Adapt the file path to your folder structure
write.table(comp_perf, file='data/SDM_alg_performances.txt', row.names=F)



pred_testdata <- data.frame(
  bc = predict(m_bc, lynx_test[,my_preds]),
  glm = predict(m_glm, lynx_test[,my_preds], type='response'),
  rf = predict(m_rf, lynx_test[,my_preds], type='response'),
  brt = predict.gbm(m_brt, lynx_test[,my_preds], 
		                 n.trees=m_brt$gbm.call$best.trees, type="response")
)

summary(pred_testdata)

# We use the sapply function to apply the thresholding to all columns 
# in the prediction data frame. You could also use a loop or construct
# the data frame by hand.
binpred_testdata <- sapply(names(pred_testdata), FUN=function(alg){ 
  ifelse(pred_testdata[,alg]>=comp_perf[comp_perf$alg==alg,'thresh'],1,0)
})

summary(binpred_testdata)


## ------------------------------------------------------------------------------------------------
# Mean of probabilities
mean_prob <- rowMeans(pred_testdata)

# Median of probabilities
median_prob <- apply(pred_testdata, 1, median)

# Weighted mean of probabilities, weighted by TSS 
# (Make sure that order of models is the same in df for predictions and performance!!)
wmean_prob <- apply(pred_testdata, 1, weighted.mean, w=comp_perf[,'TSS'])

# Committee averaging of binary predictions: calculates the proportion of
# models that predict the species to be present.
committee_av <- rowSums(binpred_testdata)/ncol(binpred_testdata)

# We can also calculate uncertainty measures, 
# e.g. the standard deviation when making ensembles of mean probabilities.
sd_prob <- apply(pred_testdata, 1, sd)


## ------------------------------------------------------------------------------------------------
# performance measures for "mean of probabilities"
(perf_mean_prob <- evalSDM(lynx_test$occ, mean_prob))

# performance measures for "median of probabilities":
(perf_median_prob <- evalSDM(lynx_test$occ, median_prob))

# performance measures for "weighted mean of probabilities":
(perf_wmean_prob <- evalSDM(lynx_test$occ, wmean_prob))

# Compare:
(ens_perf <- rbind(mean_prob = perf_mean_prob, median_prob = perf_median_prob, 
                    wmean_prob = perf_mean_prob))

# Note that we do not assess the performance of the committee average as the interpretation of this ensemble is quite different: it gives the proportion of models that agree on predicted presence


## ------------------------------------------------------------------------------------------------
# Response surfaces:
# Make predictions of all models to hypothetical grid:
xyz_preds <- data.frame(
  bc = predict(m_bc, xyz),
  glm = predict(m_glm, xyz, type='response'),
  rf = predict(m_rf, xyz, type='response'),
  brt = predict.gbm(m_brt, xyz, 
		                 n.trees=m_brt$gbm.call$best.trees, type="response")
)

# Make binary predictions
xyz_bin <- sapply(names(xyz_preds), FUN=function(alg){ 
  ifelse(xyz_preds[,alg]>=comp_perf[comp_perf$alg==alg,'thresh'],1,0)
})

# Make ensembles:
xyz_ensemble <- data.frame(
  mean_prob = rowMeans(xyz_preds),
  median_prob = apply(xyz_preds, 1, median),
  wmean_prob = apply(xyz_preds, 1, weighted.mean, w=comp_perf[,'TSS']),
  committee_av = rowSums(xyz_bin)/ncol(xyz_bin),
  sd_prob = apply(xyz_preds, 1, sd)
)	

# Plot ensemble of mean probabilities:
xyz$z <- xyz_ensemble[,'mean_prob']
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='Ensemble: mean prob', xlab='bio11', 
          ylab='bio10', screen=list(z = -120, x = -70, y = 3))

# Plot ensemble of median probabilities:
xyz$z <- xyz_ensemble[,'median_prob']
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='Ensemble: median prob', xlab='bio11', 
          ylab='bio10', screen=list(z = -120, x = -70, y = 3))
	
# Plot ensemble of weighted mean probabilities:
xyz$z <- xyz_ensemble[,'wmean_prob']
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='Ensemble: weighted mean prob', xlab='bio11', 
          ylab='bio10', screen=list(z = -120, x = -70, y = 3))

# Plot ensemble of committee average. This provides the proportion of models that agree on predicted presence:
xyz$z <- xyz_ensemble[,'committee_av']
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='Ensemble: committee average', xlab='bio11', 
          ylab='bio10', screen=list(z = -120, x = -70, y = 3))

# Plot standard deviation of mean probabilities. This gives us an indication where in environmental space we have highest uncertainty:
xyz$z <- xyz_ensemble[,'sd_prob']
wireframe(z ~ bio11 + bio10, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='Ensemble: sd', xlab='bio11', ylab='bio10', 
          screen=list(z = -120, x = -70, y = 3))


## ------------------------------------------------------------------------------------------------
# Prepare data frame with environmental data
bio_curr_df <- data.frame(rasterToPoints(bio_curr[[my_preds]]))

# We make predictions of all models:
env_preds <- data.frame(bio_curr_df[,1:2], 
  bc = predict(m_bc, bio_curr_df),
  glm = predict(m_glm, bio_curr_df, type='response'),
  rf = predict(m_rf, bio_curr_df, type='response'),
  brt = predict.gbm(m_brt, bio_curr_df, 
		                 n.trees=m_brt$gbm.call$best.trees, type="response"))
                        

# Binarise predictions of all algorithms
env_preds_bin <- data.frame(bio_curr_df[,1:2],
  sapply(names(env_preds[,-c(1:2)]), FUN=function(alg){ 
    ifelse(env_preds[,alg]>=comp_perf[comp_perf$alg==alg,'thresh'],1,0)
  }))

# Make rasters from predictions:
r_preds <- rasterFromXYZ(env_preds, crs=crs(bg))
r_preds_bin <- rasterFromXYZ(env_preds_bin, crs=crs(bg))

# Map predicted occurrence probabilities:
spplot(r_preds)

# Map predicted presences:
spplot(r_preds_bin)


## ------------------------------------------------------------------------------------------------
# We make ensembles:	
env_ensemble <- data.frame(bio_curr_df[,1:2], 
  mean_prob = rowMeans(env_preds[,-c(1:2)]),
  median_prob = apply(env_preds[,-c(1:2)], 1, median),
  wmean_prob = apply(env_preds[,-c(1:2)], 1, weighted.mean, w=comp_perf[,'TSS']),
  committee_av = rowSums(env_preds_bin[,-c(1:2)])/ncol(env_preds_bin[,-c(1:2)]),
  sd_prob = apply(env_preds[,-c(1:2)], 1, sd))

# Make rasters from ensemble predictions:
r_ens <- rasterFromXYZ(env_ensemble, crs=crs(bg))

# Map continuous ensemble predictions:
spplot(r_ens[[1:4]])


## ------------------------------------------------------------------------------------------------
# Map standard deviation across model algorithms:
plot(r_ens[['sd_prob']])


## ------------------------------------------------------------------------------------------------
# Binarise ensemble predictions
env_ensemble_bin <- data.frame(bio_curr_df[,1:2], 
    sapply(c('mean_prob', 'median_prob', 'wmean_prob'), 
           FUN=function(x){ifelse(env_ensemble[,x]>= ens_perf[x,'thresh'],1,0)}))

# Make rasters:
r_ens_bin <- rasterFromXYZ(env_ensemble_bin, crs=crs(bg))

# Map predicted presence from ensembles:
spplot(r_ens_bin)	

