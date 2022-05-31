# Session 5: SDM introduction

library(raster)
library(letsR)

load('data/r_mammals_eur.RData')
lynx_dist <- data.frame(r_mammals_eur$Presence_and_Absence_Matrix[,1:2],
                        occ = r_mammals_eur$Presence_and_Absence_Matrix[,'Lynx lynx'])

# plot the distribution
plot(rasterFromXYZ(lynx_dist))

# Set download=T if you haven't downloaded the data in the previous practical
clim <- getData("worldclim", var="bio", res=10, download=F, path="data")


## ------------------------------------------------------------------------------------------------
# Aggregate to 1° resolution
clim_1deg <- aggregate(clim, 6)

# Merge with Lynx data
lynx_dist <- data.frame(lynx_dist,
                        raster::extract(x = clim_1deg, y = lynx_dist[,1:2]))

# Exclude NAs (cells without climate data)
lynx_dist <- na.exclude(lynx_dist)


## ------------------------------------------------------------------------------------------------
library(dismo)

# We create a background map at 1° spatial resolution with European land masses - this is only for plotting purposes
bg <- rasterFromXYZ(lynx_dist[,1:3])
values(bg)[!is.na(values(bg))] <- 0

# We only want to have one record per 2°x2° and use a mask to achieve this
mask <- aggregate(bg,2)

# Thinned coordinates
set.seed(54321)
xy <- gridSample(lynx_dist[,1:2], mask, n=1)

# Checking number of data points 
nrow(xy)
nrow(lynx_dist)

# Reducing the Lynx data
lynx_thinned <- merge(xy,lynx_dist)

# Plot European land mass
plot(bg,col='grey',axes=F,legend=F)

# Plot presences in red and absences in black
points(lynx_thinned[,1:2],pch=19,col=c('black','red')[as.factor(lynx_thinned$occ)])


## ------------------------------------------------------------------------------------------------
# We first fit a GLM for the bio11 variable assuming a linear relationship:
m1 <- glm(occ ~ bio11, family="binomial", data= lynx_thinned)
	
# We can get a summary of the model:
summary(m1)	


# Fit a quadratic relationshop with bio1:
m1_q <- glm(occ ~ bio11 + I(bio11^2), family="binomial", data= lynx_thinned)

# Fit a quadratic relationship with bio11:
m1_q <- glm(occ ~ bio11 + I(bio11^2), family="binomial", data= lynx_thinned)
summary(m1_q)
	
# Or use the poly() function:
summary( glm(occ ~ poly(bio11,2) , family="binomial", data= lynx_thinned) )

# Fit two linear variables:
summary( glm(occ ~ bio11 + bio8, family="binomial", data= lynx_thinned) )

# Fit three linear variables:
summary( glm(occ ~ bio11 + bio8 + bio17, family="binomial", data= lynx_thinned) )

# Fit three linear variables with up to three-way interactions
summary( glm(occ ~ bio11 * bio8 * bio17, family="binomial", data= lynx_thinned) )

# Fit three linear variables with up to two-way interactions
summary( glm(occ ~ bio11 + bio8 + bio17 +
              bio11:bio8 + bio11:bio17 + bio8:bio17,
            family="binomial", data= lynx_thinned) )


library(corrplot)

# We first estimate a correlation matrix from the predictors. 
# We use Spearman rank correlation coefficient, as we do not know 
# whether all variables are normally distributed.
cor_mat <- cor(lynx_thinned[,-c(1:3)], method='spearman')

# We can visualise this correlation matrix. For better visibility, 
# we plot the correlation coefficients as percentages.
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)


# # Installing mecofun
# library(devtools)
# devtools::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")


library(mecofun)

# Run select07()
var_sel <- select07(X=lynx_thinned[,-c(1:3)], 
                    y=lynx_thinned$occ, 
                    threshold=0.7)

# Check out the structure of the resulting object:
str(var_sel)

# We extract the names of the weakly correlated predictors ordered by the univariate variable importance in terms of AIC:
pred_sel <- var_sel$pred_sel


## ------------------------------------------------------------------------------------------------
# How many presence points do we have?
sum(lynx_thinned$occ)


## ------------------------------------------------------------------------------------------------
# Fit the full model:
m1 <- glm( occ ~ bio11 + I(bio11^2) + bio10 + I(bio10^2), 
               family='binomial', data=lynx_thinned)

# Inspect the model:
summary(m1)


# Explained deviance:
expl_deviance(obs = lynx_thinned$occ,
              pred = m1$fitted)


m_step <- step(m1)


# Inspect the model:
summary(m_step)

# Explained deviance:
expl_deviance(obs = lynx_thinned$occ,
              pred = m_step$fitted)

save(bg, clim_1deg, lynx_thinned, pred_sel, file='data/lynx_thinned.RData')

