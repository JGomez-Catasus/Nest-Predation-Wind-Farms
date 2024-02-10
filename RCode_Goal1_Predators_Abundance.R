
########################################################################
# Housekeeping
# Load packages from R and support function
library(readxl)
library(lattice)  
library(ggplot2)
library(MASS)
library(mgcv)
library(plyr)
library(INLA)
library(spdep)
library(ncf)
source("HighstatLibV10.R")
########################################################################

#### *A* TERRESTRIAL PREDATORS
######## *A.2* TERRESTRIAL PREDATORS - Spatial Autocorrelation
#### *B* AVIAN PREDATORS
######## *B.2* AVIAN PREDATORS - Spatial Autocorrelation

#########################################################################################################
#########################################################################################################
#########################################################################################################
##### *A* TERRESTRIAL PREDATORS
Datos <- read_excel ("Excrementos2016_20190205.xlsx")
names(Datos)
str (Datos)

# Prepare data
Datos$Potentials <- Datos$Zorro + Datos$Perro + Datos$Jabali
# We use only the potentials, because we believe that this is not a good method to estimate lagarto abundance
# On the other hand, home range is bigger for foxes, dogs and wild boars, whereas in lagarto
# (although they use roads to get warmer) they dont use them to move. We can find them if they are in rocks next to a path
# but not if they are far away

Datos$fWF <- factor (Datos$WF, levels = c(0,1), labels = c("No", "Yes"))
str(Datos)

## Standardized variables
Datos$Suelo_no_edificado_estd <- scale (Datos$Suelo_no_edificado)
Datos$Crop_surface_estd <- scale (Datos$Crop_surface)
Datos$Tree_surface_estd <- scale (Datos$Tree_surface)

#########################################################
# Start INLA

# TP_i ~ Poisson(mu_i)
# E(TP_i) = mu_i 
# var(TP_i) = mu_i

#          Sex_i + Location_i + Length_i + Location_i x Length_i
# mu_i = e

M1 <- inla (Potentials ~ fWF,
           control.compute = list(dic = TRUE),
           family = "poisson",
           data = Datos)
summary(M1)

M1 <- inla (Potentials ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd,
            control.compute = list(dic = TRUE),
            family = "poisson",
            data = Datos)
summary(M1)

#Check for overdispersion frequentist style
mu1 <- M1$summary.fitted.values[,"mean"]
E1  <- (Datos$Potentials - mu1) / sqrt(mu1)
N   <- nrow(Datos) # + 1 In negative binomial due to the dispersion parameter K
p   <- nrow(M1$summary.fixed)
Dispersion <- sum(E1^2) / (N - p)
Dispersion
# In a frequentist analysis this would indicate 
# severe overdispersion.

# But it would be nice to come to the conclusion 
# that the model is overdispersed using Bayesian 
# tools.

############################################
# We will now implement the 7-step protocol
# for assessing whether the Poisson GLM
# is over- or underdispersed.


# Step 1: Apply the model in INLA
# This may sound obvious as we already 
# executed the Poisson GLM in INLA, but we 
# need to do it again, with a small modification 
# this time. The config = TRUE option allows 
# us to simulate regression parameters in the next 
# step.

M2 <- inla(Potentials ~ fWF,
           control.compute=list(config = TRUE, dic = TRUE),
           family = "poisson",
           data = Datos)

M2 <- inla (Potentials ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd,
            control.compute=list(config = TRUE, dic = TRUE),
            family = "poisson",
            data = Datos)


# Step 2: Simulate regression parameters
# We use the function inla.posterior.sample 
# to simulate from the model. The output is 
# stored in the Sim object.

set.seed(12345)
Sim <- inla.posterior.sample(n = 1, result = M2)

Sim[[1]]$latent
# This gives an object with 162 rows for this 
# specific data set and model. The first 18 
# rows are simulated values for eta = X * beta, 
# where X is the matrix with covariates, and 
# the last 2 rows are simulated regression parameters. 
# This is just one set of simulated values 
# (due to n = 1 in the function above). 

# Step 3: Calculate predicted values
# We have multiple options to do this step. 
# We can either take the first 18 rows and 
# exponent them to get mu = exp(eta), 
# or we access the last 2 (or 4 for land uses) rows and calculate 
# the fitted values via mu = exp(X * beta).
RowNum <- 19:20  # Last 2 rows are the betas
RowNum <- 19:22 # Land uses

X      <- model.matrix(~fWF, data = Datos)
X      <- model.matrix(~Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd, data = Datos)

Betas  <- Sim[[1]]$latent[RowNum]
mu     <- exp(X %*% Betas)

# Step 4: Simulate count data
# We use the rpois function for this.

Ysim <- rpois(n = nrow(Datos), lambda = mu)

# We now have 18 simulated values that 
# correspond to the observed covariate values. 
# We can calculate how many zeros this simulated 
# data set has, determine the maximum value, 
# the minimum value, the range, etc.

# Step 5: Calculate summary statistic
# Instead of the numbers of zeros, or the 
# maximum value, we can also calculate the 
# Pearson residuals for the simulated data, 
# and square and sum them.
Es <- (Ysim - mu) /sqrt(mu) # Predicted - Mean
sum(Es^2)

# By the way, the sum of squared Pearson 
# residuals for the original data and model is:
E1 <- (Datos$Potentials - mu1) /sqrt(mu1)
SS <- sum(E1^2)
SS    #bigger!

# Step 6: Repeat steps 2 to 5 thousand times
#Repeating the simulation step a large 
# number of times requires only minor modification 
# of the code presented in step 2

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = M2)

# Now we have thousand simulated mean values 
# from the model. Processing the information 
# requires more coding. We first determine the 
# number of rows in the data (N), and then 
# create a matrix Ysim that can hold thousand 
# simulated data sets (one in each column of Ysim). 
# Then we start a loop and in each iteration we 
# extract the simulated betas and predict the 
# response variable.

N    <- nrow(Datos)
Ysim   <- matrix(nrow = N, ncol = NSim)
mu.sim <- matrix(nrow = N, ncol = NSim)

for (i in 1:NSim){
  Betas <- SimData[[i]]$latent[RowNum]
  mu.sim[,i] <- exp(X %*% Betas)
  Ysim[,i]   <- rpois(n = N, lambda = mu.sim [,i])
}

# We now have thousand simulated data sets with 
# numbers of potential predators.

# Step 7: Compare simulation results and observed data
# In this step we compare the simulation results 
# with the observed data. There are many things that 
# we can do with the simulated data sets. 
# For example, for each simulated data set we can 
# count the number of zeros, and see whether the 
# percentages of zeros are comparable to the percentage 
# of zeros in the observed data. If the model does 
# not produce enough zeros then we can start thinking 
# about zero inflated models. However, zero inflation 
# is not an issue for this data set.

# The sum of squared Pearson residuals 
# for the observed data is calculated as follows.
E1 <- (Datos$Potentials - mu1) /sqrt(mu1)
SS <- sum(E1^2)
SS

# We calculate the same statistic for 
# each of the simulated data sets.

# POISSON
SS.sim <- NULL
for(i in 1:NSim){
  e2 <- (Ysim[,i] - mu.sim[,i]) /sqrt(mu.sim[,i])
  SS.sim[i] <- sum(e2^2) 
}
sum(SS > SS.sim) / NSim 

# We would like to have a value here saying 50%, 40%
# If the result is 10% then we would have UNDERDISPERSION
# If the resulta is 90%,100$, then we have OVERDISPERSION
# Admisible values: 50%, 60&, 70%, 40%, 30%
# Cuases? No correct distribution, missing variables, etc.

# If it is 100%: The sum of squared Pearson residuals for the 
# observed data is always larger than that of 
# the simulate data sets, and this indicates that 
# the variation in the observed data is larger 
# than is allowed for by a Poisson distribution. 
# And that is called overdispersion.

#########################################
# What else can we do with the 1,000 simulated data sets?
# We can also calculate the dispersion statistic for each
# simulated data set

#RUN FROM HERE........
p        <- length(Betas)
Disp.Sim <- vector(length = NSim) #Create space

# Calculate the dispersion for each simulated data set

# POISSON
for(i in 1:NSim){
  e2 <- (Ysim[,i] - mu.sim[,i]) /sqrt(mu.sim[,i])
  Disp.Sim[i] <- sum(e2^2) / (N - p)
}

# Plot this as a table
hist(Disp.Sim,
     xlab = "Dispersion",
     ylab = "Frequency",
     xlim = c(0.6, 5),
     breaks = 25)

# And visualize the dispersion for the original Poisson GLMM
points(x = Dispersion, 
       y = 0, 
       pch = 16, 
       cex = 3, 
       col = 2)

# TO HERE ............
# The red dot is the overdispersion in the original data set.
# We have serious overdispersion!
########################################



########################
# Do we need zero inflated models?
# Calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim[,i] == 0)
}
table(zeros)

#Let's plot this as a table
plot(table(zeros), 
     #axes = FALSE,
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Number of zeros in 1000 simulated data sets",
     xlim = c(0, 20),
     main = "Simulation results")
points(x = sum(Datos$Potentials == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#The red dot is the number of zeros in the original data set.
#The data simulated from the Poisson model
# Zero inflation is not an issue here.

#########################################

###############################
#NB GLM
M3 <- inla(Potentials ~ fWF,
           control.compute=list(config = TRUE, dic = TRUE),
           family = "nbinomial",
           quantiles = c(0.025, 0.975),
           data = Datos)

M3 <- inla(Potentials ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd,
           control.compute=list(config = TRUE, dic = TRUE),
           family = "nbinomial",
           quantiles = c(0.025, 0.975),
           data = Datos)

M3 <- inla(Potentials ~ Suelo_no_edificado_estd,
           control.compute=list(config = TRUE, dic = TRUE),
           family = "nbinomial",
           quantiles = c(0.025, 0.975),
           data = Datos)

summary(M3)
# Task: Write down the fitted model

# Totalparasites_i ~ NB(mu_i, k)
# E(Totalparasites_i) = mu_i
# var(Totalparasites_i) = mu_i + mu_i^2 / k
# log(mu_i) = Covariate stuff


# We go on with model M3 (the NB GLM in which k is
# estimated).            
summary(M3)            

# We can get some info on the hyperparameters:
# Posterior mean of k
k.pd <- M3$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm <- inla.emarginal(function(x) x, k.pd)
k.pm


# We can also plot the posterior 
# distribution of k. Because we get 
# the k in marginals.hyperpar there is
# no need for things like: sigma = 1/sqrt(tau)
# Therefore the info in $marginals.hyperpar
# is similar to that in tmarginal:

k.marg <- inla.tmarginal(function(x) x, k.pd)

par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
plot(k.pd, 
     type = "l",
     xlim = c(0,10),
     xlab = expression(paste(k)),
     ylab = expression(paste("Pr(", k," | data)")))
text(1.1, 2.2, "A", cex = 1.5)

plot(k.marg, 
     type = "l",
     xlim = c(0,10),
     xlab = expression(paste(k)),
     ylab = expression(paste("Pr(", k," | data)")))
text(1.1, 2.2, "B", cex = 1.5)
# Second one is just a little bit smoother.

########################################################################################################
### MODEL VALIDATION NBINOMIAL

# We should still check whether the NB GLM
# is under- or overdispersed.
# We can do it the classical way:
mu3 <- M3$summary.fitted.values[1:18,"mean"]
E3  <- (Datos$Potentials - mu3) /sqrt(mu3 + mu3^2 / k.pm)
N <- nrow(Datos)
p <- nrow(M3$summary.fixed) + 1 # 4 Betas (modelo suelos) 2 Betas (modelo WF), and 1 is the hyperparameter k
Dispersion <- sum(E3^2) / (N - p)
Dispersion  

# We can also do the simulation study again
# to check for over- or underdispersion.

# The code below is copy-paste from above. 
set.seed(12345)

# Model WF
RowNum <- 19:20
X <- model.matrix(~fWF, data = Datos)

# Model Usos
RowNum <- 19:22  #These are the rows we want
X <- model.matrix(~Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd, data = Datos)
X <- model.matrix(~Suelo_no_edificado_estd, data = Datos)

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = M3)
N      <- nrow(Datos)
Ysim   <- matrix(nrow = N, ncol = NSim)
mu.sim <- matrix(nrow = N, ncol = NSim)

for (i in 1: NSim){
  k <- exp(SimData[[i]]$logdens$hyperpar)
  Betas <- SimData[[i]]$latent[RowNum]
  mu.sim[,i] <- exp(X %*% Betas)
  Ysim[,i] <- rnegbin(n = N, # Negative Binomial distributed simulated data
                      mu = mu.sim[,i],
                      theta = k)
}

# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?

#RUN FROM HERE........
Disp.Sim <- vector(length = NSim) #Create space

# Calculate the dispersion for each simulated data set
for(i in 1:NSim){
  e3 <- (Ysim[,i] - mu.sim[,i]) /sqrt(mu.sim[,i] + mu.sim[,i]^2 / k.pm)
  Disp.Sim[i] <- sum(e3^2) #/ (N - p)
}

# This is the summary statistic for the observed data:
SS  <- sum(E3^2) 
SS

#How often is this value larger than a simulated value
mean(SS > Disp.Sim)

# This means that the NB distribution complies with the data.

# Do we need ZERO inflated models?
# Calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim[,i] == 0)
}

#Let's plot this as a table
par (mfrow = c(1,1))
plot(table(zeros), 
     #axes = FALSE,
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Number of zeros in 1000 simulated data sets",
     xlim = c(4, 20),
     main = "Simulation results")
points(x = sum(Datos$Potentials == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#The red dot is the number of zeros in the original data set.
#The data simulated from the Poisson model
# Zero inflation is not an issue here.

############################################
# Model validation of the NB GLM
# Get fitted values and Pearson residuals
mu3 <- M3$summary.fitted.values[1:18,"mean"]

# Negative Binomial
E3  <- (Datos$Potentials - mu3) / sqrt(mu3 + mu3^2 / k.pm)

# Outliers?
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = mu3, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     
# No clear outliers.
# These bands of lines are due to the 
# discrete nature of the data.

# Plot residuls vs each covariate
Datos$E3 <- E3 # From the Negative Binomial
MyVar <- c("fWF", "Suelo_no_edificado", "Crop_surface", "Tree_surface")
MyMultipanel.ggp2(Z=Datos,
                  varx = MyVar,
                  vary = "E3",
                  ylab = "Pearson Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

## MORAN'S I TEST
Loc <- cbind(Datos$XCOORD/1000, Datos$YCOORD/1000)
D <- dist(Loc)

w <- 1/as.matrix(dist(Loc))
diag(w) <- 0

moran.test (E3, mat2listw(w))

leadI <- spline.correlog(x=Loc[,1], y=Loc[,2],
                         z=E3, resamp=100, quiet=TRUE)

par (mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
par (mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(density(E3))
plot (leadI)

## VARIOGRAM
MyData <- data.frame(E3 = E3, 
                     X  = Datos$XCOORD/1000, 
                     Y  = Datos$YCOORD/1000)
coordinates(MyData)  <- c("X", "Y")
V1 <- variogram(E3 ~ 1, 
                MyData, 
                cressie = TRUE)

#Plot both variograms in one graph
p <- ggplot()
p <- p + geom_point(data = V1,
                    aes(x = dist, 
                        y = gamma)
)
p <- p + geom_line(data = V1,
                   aes(x = dist, 
                       y = gamma),
                   col = "red")

p <- p + xlab("Distance") + ylab("Semi-variogram")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p <- p + ylim(0,1)
p

# A different GRAPH
plot(V1, 
     main = "", 
     xlab = list(label = "Distance", cex = 1.5), 
     ylab = list(label = "Semi-variogram", cex = 1.5),
     pch = 16,
     col = 1,
     cex = 1.5
)

# CONCLUSION: MORAN'S I does not show Spatial Autocorrelation, but the Variograms are not trully nice

### PLOT RESIDUALS!
# Option 1: Use different point sizes and symbols
#           based on the values of the residuals.
E3
MyCex <- 3 * abs(E3) / max(E3) + 0.5
Sign  <- as.numeric(E3 >=0) + 1
MyPch <- c(1, 16)[Sign]

utmcoor <- SpatialPoints (coords = cbind(Datos$XCOORD, Datos$YCOORD), proj4string = CRS("+proj=utm +zone=30N"))
longlat <- spTransform(utmcoor, CRS ("+proj=longlat"))
Datos$Longitude <- coordinates (longlat)[,1]
Datos$Latitude <- coordinates (longlat)[,2]

xyplot (Longitude ~ Latitude,
        data = Datos, aspect = "iso",
        cex = MyCex, pch = MyPch, col = 1, 
        xlab = list(label = "Easting", cex = 1.5),
        ylab = list(label = "Northing", cex = 1.5))

############################################################################################################
############################################################################################################
# MODEL SELECTION with the NB

k.pm #Posterior mean value of k

M5 <- inla(Potentials ~ fWF,
           control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
           family = "nbinomial",
           data = Datos)

M5 <- inla(Potentials ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd,
           control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
           family = "nbinomial",
           data = Datos)
# TASK: Write down the model that we just fitted.

#Drop fWF
M5a <- inla(Potentials ~  1,
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)

M5b <- inla(Potentials ~  Crop_surface_estd + Tree_surface_estd,
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)
M5c <- inla(Potentials ~  Suelo_no_edificado_estd + Tree_surface_estd,
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)
M5d <- inla(Potentials ~  Suelo_no_edificado_estd + Crop_surface_estd,
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)

dic  <- c(M5$dic$dic, M5a$dic$dic)
dic  <- c(M5$dic$dic, M5b$dic$dic, M5c$dic$dic, M5d$dic$dic)
waic <- c(M5$waic$waic, M5a$waic$waic)
waic <- c(M5$waic$waic, M5b$waic$waic, M5c$waic$waic, M5d$waic$waic)
Z <- cbind(dic, waic)
rownames(Z) <- c("fWF", "Null")
rownames(Z) <- c("Full", "Drop Suelo no edificado", "Drop Crops", "Drop Trees")
Z

## According to DIC --> Best model FULL MODEL
## According to wAIC --> Best model the one that drops crop surface!

M5$summary.fixed
M5c$summary.fixed

M5e <- inla(Potentials ~  Tree_surface_estd,
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)
M5f <- inla(Potentials ~  Suelo_no_edificado_estd,
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)

dic  <- c(M5c$dic$dic, M5e$dic$dic, M5f$dic$dic)
waic <- c(M5c$waic$waic, M5e$waic$waic, M5f$waic$waic)

Z <- cbind(dic, waic)
rownames(Z) <- c("Both", "Drop Suelo no edificado", "Drop Trees")
Z

M5c$summary.fixed
M5f$summary.fixed

########################################################################################################
############## PRESENTING THE RESULTS!
M5$summary.fixed

# Sketch the model fit.
# The last thing we do in this analysis 
# is to make a visualization of the optimal 
# model. We discussed in detail how to do 
# this in Section 9.8. We will use the 
# inla.make.lincombs approach. Recall that 
# these were the steps that we have to execute.

# 1. Create a grid of covariate values.
# 2. Make a design matrix using the model.matrix function.
# 3. Run the model in INLA with the inla.make.lincombs option.
# 4. Visualize the results using ggplot2.

# We follow the four steps. We first make a grid of 
# covariate values.


# 1. Create a grid of covariate values.
MyData <- ddply(Datos, 
                .(fWF))

MyData <- data.frame (Suelo_no_edificado = seq(min(Datos$Suelo_no_edificado), 
                                             max(Datos$Suelo_no_edificado), 
                                             length = 25),
                      Crop_surface = seq(min(Datos$Crop_surface), 
                                         max(Datos$Crop_surface), 
                                         length = 25),
                      Tree_surface = seq(min(Datos$Tree_surface), 
                                         max(Datos$Tree_surface), 
                                         length = 25))

MyData <- data.frame (Suelo_no_edificado = seq(min(Datos$Suelo_no_edificado), 
                                               max(Datos$Suelo_no_edificado), 
                                               length = 25),
                      Crop_surface = rep(mean (Datos$Crop_surface), 
                                         length = 25),
                      Tree_surface = rep(mean (Datos$Tree_surface), 
                                         length = 25),
                      Crop_surface_estd = rep(mean (Datos$Crop_surface_estd), 
                                         length = 25),
                      Tree_surface_estd = rep(mean (Datos$Tree_surface_estd), 
                                         length = 25))

MyData <- data.frame (Suelo_no_edificado = seq(min(Datos$Suelo_no_edificado), 
                                               max(Datos$Suelo_no_edificado), 
                                               length = 25))
                   
MyData$Suelo_no_edificado_estd <- scale (MyData$Suelo_no_edificado)
Crop_surface_estd <- scale (MyData$Crop_surface)
MyData$Tree_surface_estd <- scale (MyData$Tree_surface)

head(MyData)

# 2. Make a design matrix
Xmat <- model.matrix(~ fWF,
                     data = MyData)

Xmat <- model.matrix(~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd,
                     data = MyData)

Xmat <- model.matrix(~ Suelo_no_edificado_estd,
                     data = MyData)
head(Xmat)
Xmat <- as.data.frame(Xmat)


# 3. Run the model in INLA with the 
#    inla.make.lincombs option
lcb <- inla.make.lincombs(Xmat)
Hyper.NB <- list(size = list(initial = 1, fixed = TRUE))

M7 <- inla(Potentials ~ fWF,
           lincomb = lcb,
           control.inla = list(lincomb.derived.only = FALSE),
           control.predictor = list(#link = 1,
             compute = TRUE),
           family = "nbinomial",
           control.family = list(hyper = Hyper.NB),
           data = Datos)

M7 <- inla(Potentials ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd,
           lincomb = lcb,
           control.inla = list(lincomb.derived.only = FALSE),
           control.predictor = list(#link = 1,
             compute = TRUE),
           family = "nbinomial",
           control.family = list(hyper = Hyper.NB),
           data = Datos)

M7 <- inla(Potentials ~ Suelo_no_edificado_estd,
           lincomb = lcb,
           control.inla = list(lincomb.derived.only = FALSE),
           control.predictor = list(#link = 1,
             compute = TRUE),
           family = "nbinomial",
           control.family = list(hyper = Hyper.NB),
           data = Datos)

# 4. Visualize the results using ggplot2.

# get the predicted values
Pred.pm <- M7$summary.lincomb[,c("mean", 
                                 "0.025quant", 
                                 "0.975quant")]

print(head(M7$summary.lincomb.derived), digits = 2)
Out1 <- M7$summary.lincomb.derived

# The output above is on the log-scale; 
# the exponent has not been taken yet. 
# We could do this ourselves via:

MyData$mu.wrong   <- exp(Out1$'mean')
MyData$selo.wrong <- exp(Out1$'0.025quant')
MyData$seup.wrong <- exp(Out1$'0.975quant')

# and now we could run the ggplot2 
# code from Section 9.8 again. But strictly 
# speaking this is not correct because 
# E(exp(x)) !=  exp(E(x)) if x is a stochastic 
# variable. And we are doing a Bayesian analysis, 
# so x is stochastic. The ‘E()’ stands for the 
# expected value and exp for the exponential function. 
# This was not a problem in Chapter 9 with the 
# identity link function, but with the log link 
# function it is.

# The solution is to use support functions 
# from INLA that convert the distribution of x 
# into the distribution of exp(x). Starting point 
# is the $marginals.lincomb.derived, which 
# contains the marginal posterior distribution for 
# each predicted value, and we have 75 of them.

M7$marginals.lincomb.derived
Pred.marg <- M7$marginals.lincomb.derived

# It is slightly confusing to access the 75 
# marginal posterior distributions. The marginal 
# posterior distribution for the first predicted 
# value is obtained via:
Pred.marg[[1]] 

# This contains an entire posterior distribution, 
# which can be plotted if you want to. What we 
# really want is its posterior mean (or median) 
# and a 95% credible interval. That is obtained with
inla.qmarginal(c(0.025, 0.5, 0.975), 
               inla.tmarginal(exp, Pred.marg[[1]]))

# For the posterior mean use:
inla.emarginal(exp, Pred.marg[[1]])

# Here is some code that does it for all 75 points
# at once:
MyData$mu <- unlist( 
  lapply(
    Pred.marg,
    function(x) inla.emarginal(exp,x)))

MyData$selo <- unlist( 
  lapply(
    Pred.marg,
    function(x) 
      inla.qmarginal(c(0.025), 
                     inla.tmarginal(exp, x))
  ))

MyData$seup <- unlist( 
  lapply(
    Pred.marg,
    function(x) 
      inla.qmarginal(c(0.975), 
                     inla.tmarginal(exp, x))
  ))

# If you don’t like fancy coding (like us), 
# just run the following loop. It repeats the 
# steps that we just explained for each predicted 
# value. Using a loop means that we have to think 
# less. The results are exactly the same.

for (i in 1:18){
  MyData$mu2[i]  <- inla.emarginal(exp, Pred.marg[[i]])
  lo.up <- inla.qmarginal(c(0.025, 0.975), 
                          inla.tmarginal(exp, Pred.marg[[i]]))
  MyData$selo2[i] <- lo.up[1]
  MyData$seup2[i] <- lo.up[2]    	
}               

MyData

# And the rest is a matter of plotting it all. 
p <- ggplot(Datos, aes(x=fWF, y=Potentials)) + 
  geom_boxplot()
# Box plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=2)
# Box plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0), size = 5) + xlab ("Wind farms presence") + ylab ("Predator abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(family = "Times", size = 20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 20, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text = element_text(family = "Times",size=18))


######################################################################################################
# And the rest is a matter of plotting it all. 
p <- ggplot()
p <- p + geom_point(data = Datos, 
                    aes(y = Potentials, x = Suelo_no_edificado),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Suelo no edificado") + ylab("Potential predators")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_line(data = MyData, 
                   aes(x = Suelo_no_edificado, 
                       y = mu), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Suelo_no_edificado, 
                         ymax = seup, 
                         ymin = selo),
                     alpha = 0.2)
p

p <- ggplot()
p <- p + geom_point(data = Datos, 
                    aes(y = Potentials, x = Tree_surface),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Tree surface") + ylab("Potential predators")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_line(data = MyData, 
                   aes(x = Tree_surface, 
                       y = mu), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Tree_surface, 
                         ymax = seup, 
                         ymin = selo),
                     alpha = 0.2)
p

range(MyData$Tree_surface)
range(Datos$Tree_surface)
range(MyData$Suelo_no_edificado)


#########################################################################################################
##### *A.2* OVERDISPERSION! INCORPORATE SPATIAL AUTOCORRELATION
library(lattice)
library(sp)
library(raster)
library(dismo)
library(splancs)
library(INLA)
library(reshape)
library(gstat)
library(ggplot2)
library(ggmap)
library(rworldmap)
library(rgdal)
library(rgeos)

##################################
# Transform coordinates from UTM to LONG/LAT
utmcoor <- SpatialPoints (coords = cbind(Datos$XCOORD, Datos$YCOORD), proj4string = CRS("+proj=utm +zone=30N"))
longlat <- spTransform(utmcoor, CRS ("+proj=longlat"))
Datos$Longitude <- coordinates (longlat)[,1]
Datos$Latitude <- coordinates (longlat)[,2]

####################################################
# We will implement the following 8 steps.
# 1. Make a mesh.
# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
# 3. Define the SPDE.
# 4. Define the spatial field.
# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.
# 6. Specify the model formula in terms of the 
#    response variable, covariates and the 
#    spatial correlated term.
# 7. Run the spatial model in INLA.
# 8. Inspect the results.



########################
#1. Make a mesh.
#   Step 1 of making a mesh:  Get a 
#   sense for the distribution of 
#   distances between sampling locations. 
Loc <- cbind(Datos$XCOORD/1000, Datos$YCOORD/1000)
D <- dist(Loc)
max(D)
mean(D)

par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
plot (x = sort(D),
      y = (1:length(D))/length(D),
      type = "l",
      xlab = "Distance between sites (km)",
      ylab = "Cumulative proportion")

Position <- which((1:length(D))/length(D) > 0.49 & (1:length(D))/length(D) < 0.5)
sort(D)[Position]
# El 50% de las zonas estan a menos de 8 km

# We need a grid on top of our sampling points
# Mesh basada en las sampling locations

Bound <- inla.nonconvex.hull(Loc, convex= -0.15)
mesh    <- inla.mesh.2d(boundary = Bound, 
                        max.edge=  c(1,2) # Maxima distancia permitida
)
mesh$n

# A different approach
# RangeGuess <- 8 # Greater frequency
# MaxEdge <- RangeGuess / 5
# mesh    <- inla.mesh.2d(boundary = ConvHull, 
#                         max.edge=  c(1,5) * MaxEdge,
#                         cutoff = MaxEdge / 5) # Without outer part

# Bellow different script without outer part
# Recommended settings 
# See Section 3.1 at: https://haakonbakka.bitbucket.io/btopic104.html
# mesh    <- inla.mesh.2d(boundary = ConvHull, 
# max.edge=  c(1.5))
# max.edge: maximum allowed triangle edge lengths in 
#           the inner domain and in the outer extension
# cutoff: minimum allowed distance between points. Points 
#         at a closer distance than the supplied value are 
#         replaced by a single vertex

par(mfrow=c(1,1), mar=c(0,0,2,0))
plot(mesh, asp= 1, main= "")
points(Loc, col = 1, pch = 16, cex = 1)
mesh$n #We can cope with a finer mesh.

#########################################
# Step 2. Define the weighting factors a_ik (also called 
#         the projector matrix).
A2 <- inla.spde.make.A(mesh, loc = Loc)
dim(A2)  #18 observations on a 2485 grid
#        A is a weight matrix
head(A2)
############################################



############################################
# Step 3. Define the SPDE.
#### Basic Style
spde <- inla.spde2.matern(mesh, alpha=2)

# Use PC priors
# P(Range < range0) = alpha  and P(sigma > sigma0) = alpha

# We are going to use:
# P(Range < 1) = 0.5  # DIFFUSE PRIOR: No idea...it can be anything..though we think the range is small-ish
# and P(sigma > 1) = 0.05

# The 1km comes from our desire to avoid overfitting.
# The sigma....let's assume that the covariates are not important
M1 <- glm(Potentials ~ 1 , data = Datos, family = "poisson")
summary(M1)

# This means
# E[NEggs] = exp(Intercept) = exp(0.2007) = 1.22
range(Datos$Potentials)
# How big do the spatial correlated residuals u_iu need to be 
# if we want to cover 0 to 6 and:

# E[Potentials] = exp(0.2007 + u)
#               = exp(0.2007) *  exp(u)
#               = 1.22* exp(u)

# The exp(u) should not be larger than 5.
# Hence, the u should not be not bigger than log(6) = 1.8-ish.
# Because we assume u_i ~ N(0, sigma^2), using sigma = 1 should do 
# the job.

# P(sigma > 1) = 0.05

# Below an uninformative prior: most likely the range is larger than 1 km
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(1, 0.05), 
                            prior.sigma = c(1, 0.05))

# Diffuse prior. No Idea, but we think that the range is small
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(1, 0.5), 
                            prior.sigma = c(1, 0.05))

# Diffuse prior. No Idea, but we think that the range is small <--- WE USE THIS ONE
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(7, 0.5), 
                            prior.sigma = c(3, 0.05))

##########################################
# Step 4. Define the spatial field.
# Next we set up a list for the spatial random 
# intercept u. As explained in Section 13.6, 
# u is rewritten internally as A * w. We need 
# to specify the w in INLA. This is done with 
# the inla.spde.make.index function

# The size of the mesh is: 
mesh$n

# This number is also in
spde$n.spde

# It is also the number of w_k values that we will get.
# For this mesh we use:
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde,
  n.group = 1,
  n.rep1 = 1)

str(w.index)
# Here is where spatial-temporal models can be defines.
# We will do that later. 
#####################################

#####################################
# Step 5.	Make a stack. 

# This is a confusing step. We discussed it
# in detail in Chapter 13.

# Make the X matrix

# Set up the model. 
# Create a data frame with an intercept and covariates.
N <- nrow(Datos)

# WIND FARMS - Factors
Xm <- model.matrix(~ WF, data= Datos)

X <- data.frame(Intercept = Xm[,1],
                fWF = Xm[,2])

# LAND USES - No Factors
X <- data.frame(Intercept = rep (1,N),
                Suelo_no_edificado_estd = Datos$Suelo_no_edificado_estd,
                Crop_surface_estd = Datos$Crop_surface_estd,
                Tree_surface_estd = Datos$Tree_surface_estd)

X <- as.data.frame(X) #Avoids problems

# And here is the stack.
StackFit <- inla.stack(
  tag  = "Fit",
  data = list(y = Datos$Potentials),  
  A    = list(1,A2),  # Intercept and covariates, spatial field                    
  effects = list( 
    X = X, #Covariates
    w = w.index          #Spatial field  
  )) 

#############################################
#6.	Specify the model formula in terms of the 
#   response variable, covariates and the 
#   spatial correlated term.

# These are the models that we will fit:
# Y_i ~ Poisson(mu_i)
# E(Y_i) = mu_i
# Model 1: log(mu_i) = Covariate stuff
# Model 2: log(mu_i) = Covariate stuff + u_i
f1 <- y ~ -1 + Intercept + fWF
f1 <- y ~ -1 + Intercept + Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd

f2 <- y ~ -1 + Intercept + fWF + 
  f(w, model = spde)
f2 <- y ~ -1 + Intercept + Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + 
  f(w, model = spde)

#############################################
# 7. Run the spatial model in INLA.
# First we run the model without spatial dependency.

I1.poisson <- inla(f1,
                   family = "poisson", 
                   data = inla.stack.data(StackFit),
                   control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
                   control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit)))

I1.nb <- inla(f1,
              family = "nbinomial", 
              data = inla.stack.data(StackFit),
              control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
              control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit)))


# And this is the model with the spatial field:
#  Y_i ~ Poisson(mu_i)
#  E(Y_i)   = mu_i
#  var(Y_i) =  mu_i
#  log(mu_i) = Intercept + Covariates + u_i
#  Where u_i is spatially correlated noise.

I2.poisson <- inla(f2,
                   family = "poisson", 
                   data = inla.stack.data(StackFit),
                   control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
                   control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit)))

I2.nb <- inla(f2,
              family = "nbinomial", 
              data = inla.stack.data(StackFit),
              control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
              control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit)))

# And compare the models with DICs and WAICs
dic  <- c(I1.poisson$dic$dic, I1.nb$dic$dic, I2.poisson$dic$dic, I2.nb$dic$dic)
waic <- c(I1.poisson$waic$waic, I1.nb$waic$waic, I2.poisson$waic$waic, I2.nb$waic$waic)
Z     <- cbind(dic, waic)
rownames(Z) <- c("Poisson GLM",  "NB GLM", 
                 "Poisson GLM + SPDE", "NB GLM + SPDE")
Z

# For the effect of wind farms
#                       dic      waic
# Poisson GLM        68.72457  74.23878
# NB GLM             63.90720  66.99557 ***
# Poisson GLM + SPDE 45.65335 501.69332
# NB GLM + SPDE      46.60376  45.44446 <---

# For the effect of land uses 
#                           dic         waic
# Poisson GLM           61.30920 6.993844e+01
# NB GLM                59.51933 6.483316e+01
# Poisson GLM + SPDE 35863.70553 2.770178e+73
# NB GLM + SPDE         46.67300 8.666628e+01

I1.poisson$summary.fixed
I2.nb$summary.fixed


## MODEL SELECTION
dic  <- c(I1.poisson$dic$dic, I1a.poisson$dic$dic, I1b.poisson$dic$dic, 
          I1.nb$dic$dic, I1a.nb$dic$dic, I1b.nb$dic$dic,
          I2.poisson$dic$dic, I2a.poisson$dic$dic, I2b.poisson$dic$dic,
          I2.nb$dic$dic, I2a.nb$dic$dic, I2b.nb$dic$dic)

waic <- c(I1.poisson$waic$waic, I1a.poisson$waic$waic, I1b.poisson$waic$waic,
          I1.nb$waic$waic, I1a.nb$waic$waic, I1b.nb$waic$waic,
          I2.poisson$waic$waic, I2a.poisson$waic$waic, I2b.poisson$waic$waic, 
          I2.nb$waic$waic, I2a.nb$waic$waic, I2b.nb$waic$waic)
Z     <- cbind(dic, waic)
rownames(Z) <- c("Poisson GLM Full",  "Poisson GLM Drop Crop", "Poisson GLM Only Suelo",
                 "NB GLM Full",  "NB GLM Drop Crop", "NB GLM Only Suelo", 
                 "Poisson GLM + SPDE Full", "Poisson GLM + SPDE Drop Crop", "Poisson GLM + SPDE Only Suelo",
                 "NB GLM + SPDE Full", "NB GLM + SPDE Drop Crop", "NB GLM + SPDE Only Suelo")
Z

#######################################
#8. Inspect the results of I2

# Fixed parameters
# Here are the results of the model without and with the
# spatial random effects

I1.poisson$summary.fixed[, c("mean", "0.025quant", "0.975quant")]
I2.poisson$summary.fixed[, c("mean", "0.025quant", "0.975quant")]

# Let's plot the results of the model, without and with
# the spatial correlation side by side.
# Better not try to understand all this R code.
# RUN FROM HERE......

Combined <- rbind(I1.poisson$summary.fixed[, c("mean", "0.025quant", "0.975quant")],
                  I2.poisson$summary.fixed[, c("mean", "0.025quant", "0.975quant")])

Combined$WhichModel <- rep(c("Non-spatial", "Spatial"), 
                           each = nrow(I2.poisson$summary.fixed))

rownames(I2.poisson$summary.fixed)

Combined$WhichVariable <- rep(c("Intercept", "Wind farms"), 2)
Combined$WhichVariable <- rep(c("Intercept", "S.Gravel roads", "S.Crops", "S.Trees"), 2)

colnames(Combined) <- c("Mean", "Lo", "Up", "Model", "Predictor")
positions <- unique(Combined$Predictor)

windows()
ggplot(Combined, aes(y=Mean, x= Predictor)) +
  geom_errorbar(
    aes(ymin = Lo, ymax = Up, color = Model),
    position = position_dodge(0.3), width = 0.2, size= 1) +
  geom_hline(yintercept = 0, size = 0.7, linetype = "dashed") +
  geom_point(aes(color = Model), position = position_dodge(0.3), size= 2) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme(legend.position="top",
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(family = "Times", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text = element_text(family = "Times",size=14)) + 
  scale_x_discrete(limits = positions)
# TO HERE ......

# For the NEGATIVE BINOMIAL MODELS
Combined <- rbind(I1.nb$summary.fixed[, c("mean", "0.025quant", "0.975quant")],
                  I2.nb$summary.fixed[, c("mean", "0.025quant", "0.975quant")])

Combined$WhichModel <- rep(c("Non-spatial", "Spatial"), 
                           each = nrow(I2.nb$summary.fixed))

rownames(I2.nb$summary.fixed)

Combined$WhichVariable <- rep(c("Intercept", "Wind farms"), 2)
Combined$WhichVariable <- rep(c("Intercept", "S.Gravel roads", "S.Crops", "S.Trees"), 2)

colnames(Combined) <- c("Mean", "Lo", "Up", "Model", "Predictor")
positions <- unique(Combined$Predictor)

windows()
ggplot(Combined, aes(y=Mean, x= Predictor)) +
  geom_errorbar(
    aes(ymin = Lo, ymax = Up, color = Model),
    position = position_dodge(0.3), width = 0.2, size= 1) +
  geom_hline(yintercept = 0, size = 0.7, linetype = "dashed") +
  geom_point(aes(color = Model), position = position_dodge(0.3), size= 2) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme(legend.position="top",
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(family = "Times", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text = element_text(family = "Times",size=14)) + 
  scale_x_discrete(limits = positions)

###############################################################################
## HYPERPARAMETERS:

# Let's focus on the hyper-parameters.
# The code is a 'little' bit on the ugly side.
# Just run it without asking why.

# The structure of the NB model is as follows
# TotalPredators_i ~ NB(mu_i, k)
# E(TotalPredators_i) = mu_i
# var(TotalPredators_i) = mu_i + mu_i^2 / k
# log(mu_i) = Intercept + Covariates + u_i
#  Where u_i is spatially correlated noise.

# We have the mean (mu_i) and shape (k) parameters
# As its name implies, the negative binomial shape parameter, k, 
# describes the shape of a negative binomial distribution. In other words, 
# k is only a reasonable measure to the extent that your data represent a 
# negative binomial distribution. The broader the k, The broader the distribution
I2 <- I2.poisson
I2 <- I2.nb

SpatField.w <- inla.spde2.result(inla = I2,
                                 name = "w",
                                 spde = spde,
                                 do.transfer = TRUE)

Kappa <- inla.emarginal(function(x) x, 
                        SpatField.w$marginals.kappa[[1]] )

Sigma_u <- inla.emarginal(function(x) sqrt(x), 
                          SpatField.w$marginals.variance.nominal[[1]] )

Range <- inla.emarginal(function(x) x, 
                        SpatField.w$marginals.range.nominal[[1]] )

Kappa
Sigma_u
Range       #Distance at which the correlation diminishes
I2.poisson$summary.hyperpar

# This is perhaps a nicer graph to make and present.
# Show correlation structure
# First we obtain the locations of each point of the mesh.
LocMesh <- mesh$loc[,1:2]

# And then we calculate the distance between each vertex.
D <- as.matrix(dist(LocMesh))

# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, max(D), length = 100)      
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1) 
Cor.M[1] <- 1

# Which we plot here:
par(mfrow=c(1,1), mar = c(5,5,2,2))
plot(x = d.vec, 
     y = Cor.M, 
     pch = 16, 
     type = "l", 
     cex.lab = 1.5,
     xlab = "Distance", 
     ylab = "Correlation",
     xlim = c(0, 30))
abline(h = 0.1, lty = 2)

# Negative Binomial: posterior mean of k
# As its name implies, the negative binomial shape parameter, k, 
# describes the shape of a negative binomial distribution. In other words, 
# k is only a reasonable measure to the extent that your data represent a 
# negative binomial distribution. The broader the k, The broader the distribution
NB <- I1.nb
NB <- I2.nb

k.pd <- NB$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm <- inla.emarginal(function(x) x, k.pd)
k.pm


############################################################################################################
##### MODEL VALIDATION
#### *1* Plot fitted values versus observed data
#### *2* Check the Pearson residuals for any remaining spatial dependency - Variogram and Moran's I
#### *3* PLOT RESIDUALS
#### *4* RESIDUALS VERSUS EACH COVARIATE
#### *5* Proportion of zeros and overdispersion

########################################################
Poisson <- I1.poisson # Poisson
NB <- I1.nb # Negative Binomial

Poisson <- I2.poisson # Poisson + SPDE
NB <- I2.nb # Negative Binomial + SPDE

########################################################
#### *1* Plot fitted values versus observed data
# Poisson
mu_poisson <- Poisson$summary.fitted.values[1:18, "mean"]
plot(x = mu_poisson,
     y = Datos$Potentials)

# Negative Binomial
mu_nb <- NB$summary.fitted.values[1:18, "mean"]
plot(x = mu_nb,
     y = Datos$Potentials)
# Hmm. Looks too good?
# Are we overfitting?

########################################################
#### *2* Check the Pearson residuals for any remaining spatial dependency.
# Variogram

# Poisson
E_poisson <- (Datos$Potentials - mu_poisson) / sqrt (mu_poisson) # Poisson Pearson Residuals

MyData <- data.frame(E1 = E_poisson, 
                     X  = Datos$XCOORD/1000, 
                     Y  = Datos$YCOORD/1000)
coordinates(MyData)  <- c("X", "Y")
V_poisson <- variogram(E_poisson ~ 1, 
                MyData, 
                cressie = TRUE)


# Negative Binomial
k.pd <- NB$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm <- inla.emarginal(function(x) x, k.pd)
k.pm

E_nb <- (Datos$Potentials - mu_nb) / sqrt (mu_nb + mu_nb^2/k.pm) # Negative BinomialPearson Residuals

MyData <- data.frame(E1 = E_nb, 
                     X  = Datos$XCOORD/1000, 
                     Y  = Datos$YCOORD/1000)
coordinates(MyData)  <- c("X", "Y")
V_nb <- variogram(E_nb ~ 1, 
                MyData, 
                cressie = TRUE)

# Plot both variograms in one graph
p <- ggplot()
p <- p + geom_point(data = V_poisson,
                    aes(x = dist, 
                        y = gamma)
)
p <- p + geom_line(data = V_poisson,
                   aes(x = dist, 
                       y = gamma),
                   col = "red")

p <- p + geom_point(data = V_nb,
                    aes(x = dist, 
                        y = gamma)
)
p <- p + geom_line(data = V_nb,
                   aes(x = dist, 
                       y = gamma),
                   col = "blue"
)

p <- p + xlab("Distance") + ylab("Semi-variogram")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p <- p + ylim(0,1)
p
# No clear patterns for the Poisson GLM with spatial correlation
# But strange patterns for the NB GLM with spatial correlation

### MORAN'S I
## MORAN'S I TEST
E <- E_poisson
E <- E_nb

Loc <- cbind(Datos$XCOORD/1000, Datos$YCOORD/1000)
D <- dist(Loc)
w <- 1/as.matrix(dist(Loc))
diag(w) <- 0

moran.test (E, mat2listw(w))
leadI <- spline.correlog(x=Loc[,1], y=Loc[,2],
                         z=E, resamp=100, quiet=TRUE)
par (mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
plot(density(E))
plot (leadI)

########################################################
#### 3* PLOT RESIDUALS
# Option 1: Use different point sizes and symbols
#           based on the values of the residuals.
E <- E_poisson
E <- E_nb

MyCex <- 3 * abs(E) / max(E) + 0.5
Sign  <- as.numeric(E >=0) + 1
MyPch <- c(1, 16)[Sign]

xyplot (Longitude ~ Latitude,
        data = Datos, aspect = "iso",
        cex = MyCex, pch = MyPch, col = 1, 
        xlab = list(label = "Easting", cex = 1.5),
        ylab = list(label = "Northing", cex = 1.5))

# Patterns in the residuals are present! 
# Black dots is for positive residuals (left side of the island) whereas white dots is for negative residuals (right side of the island)

# If you see patterns then we have to stop! Standard errors, t values and estimates are not true -> They are biased
# Standard errors are faluty, and then also de t-value and also the p-values (Domino!)
# If the p-value is really low (0.0000) the effects are going to be low. But if the p-values are near 0.01 then we have troubles!

########################################################
#### *4* RESIDUALS VERSUS EACH COVARIATE
# Plot residuls vs each covariate
Datos$E <- E_poisson # From the Poisson
Datos$E <- E_nb # From the Negative Binomial

MyVar <- c("fWF", "Suelo_no_edificado", "Crop_surface", "Tree_surface")
MyMultipanel.ggp2(Z=Datos,
                  varx = MyVar,
                  vary = "E",
                  ylab = "Pearson Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

########################################################
### *5* Proportion of zeros and overdispersion

# Simulating from the model

# In the Turco parasite example we implemented a
# 7-step protocol for simulating from the model.
# See also Powerpoint presentation: 
#    P15_SpatTempBook_Chapter10_V1.pdf


# We will now implement the 7-step protocol
# for assessing whether the Poisson GLM
# is over- or underdispersed.

# Step 1: Apply the model in INLA
# This may sound obvious as we already 
# executed the Poisson GLM in INLA, but we 
# need to do it again, with a small modification 
# this time. The config = TRUE option allows 
# us to simulate regression parameters in the next 
# step. ---> We did it before!

# WITHOUT SPDE
# Poisson
Poisson.sim <- I1.poisson
# Negative Binomial
NB.sim <- I1.nb

#################
# WITH SPDE
# Poisson
Poisson.sim <- I2.poisson
# Negative Binomial
NB.sim <- I2.nb

###########

k.pd <- NB.sim$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm <- inla.emarginal(function(x) x, k.pd)
k.pm
# Step 2: Simulate regression parameters
# We use the function inla.posterior.sample 
# to simulate from the model. The output is 
# stored in the Sim object.

# Here is some fancy code to grap the 1  betas (the same for both models)
Sim <- inla.posterior.sample(n = 1, result = Poisson.sim) # Simulate one dataset from the poisson

MyParams <- rownames(Poisson.sim$summary.fixed)
RowNum.Betas <- lapply (MyParams, 
                        function (x)
                          grep (x,
                                rownames (Sim[[1]]$latent),
                                fixed = TRUE))
RowNum.Betas <- as.numeric (RowNum.Betas)
RowNum.Betas # In this rows are our parameters!

# In case of spatial autocorrelation
# Here is some fancy code to grap the w 
Nw <- mesh$n 
# The w are labelled as w:1 to w:Nw   where Nw is the number of vertices
MyParams <- paste("w", seq(from = 1, to = Nw), sep = ":")
MyID <- function(x){ which(rownames(Sim[[1]]$latent) == x) }
RowNum.w <- lapply(MyParams, MyID)
RowNum.w <- as.numeric(RowNum.w)
RowNum.w

############# Simulate 1000 datasets
## Poisson without SPDE
NSim <- 1000
SimData_Poisson <- inla.posterior.sample(n = NSim, result = Poisson.sim)

N  <- nrow(Datos)
Ysim_poisson <- matrix(nrow = N, ncol = NSim)
mu.i_poisson <- matrix(nrow = N, ncol = NSim)

Xm <- as.matrix(X)

for (i in 1:NSim){
  Betas <- SimData_Poisson[[i]]$latent[RowNum.Betas]
  FixedPart   <- Xm %*% Betas
  mu.i_poisson[,i]    <- exp(FixedPart)
  Ysim_poisson[,i]    <- rpois(n = nrow(Datos), lambda = mu.i_poisson[,i])
}

table(Ysim_poisson)
par(mfrow= c(1,1), mar= c(5,5,5,5))
plot (table(Ysim_poisson),
      xlab= "Simulated predator abundance values",
      ylab= "Frequencies")

## Negative Binomial without SPDE
NSim <- 1000
SimData_NB <- inla.posterior.sample(n = NSim, result = NB.sim)

N  <- nrow(Datos)
Ysim_nb <- matrix(nrow = N, ncol = NSim)
mu.i_nb <- matrix(nrow = N, ncol = NSim)
k.i_nb <- vector (length = NSim)

Xm <- as.matrix(X)

for (i in 1:NSim){
  Betas <- SimData_NB[[i]]$latent[RowNum.Betas]
  FixedPart   <- Xm %*% Betas
  mu.i_nb[,i]    <- exp(FixedPart)
  k.i_nb [i]    <- exp(SimData_NB[[i]]$logdens$hyperpar)
  Ysim_nb[,i]    <- rnegbin(n = nrow(Datos), mu = mu.i_nb[,i], theta = k.i_nb[i])
}

table(Ysim_nb)
plot (table(Ysim_nb),
      xlab= "Simulated predator abundance values",
      ylab= "Frequencies")


#######################################
## Poisson with SPDE
NSim <- 1000
SimData_Poisson <- inla.posterior.sample(n = NSim, result = Poisson.sim)

N  <- nrow(Datos)
Ysim_poisson <- matrix(nrow = N, ncol = NSim)
mu.i_poisson <- matrix(nrow = N, ncol = NSim)

Xm <- as.matrix(X)
Am <- as.matrix(A2)

for (i in 1:NSim){
  Betas <- SimData_Poisson[[i]]$latent[RowNum.Betas]
  wk    <- SimData_Poisson[[i]]$latent[RowNum.w]
  FixedPart   <- Xm %*% Betas
  SpatialPart <- Am %*% wk
  mu.i_poisson[,i]    <- exp(FixedPart + SpatialPart)
  Ysim_poisson[,i]    <- rpois(n = N, lambda = mu.i_poisson[,i])
}

table(Ysim_poisson)
par(mfrow= c(1,1), mar= c(5,5,5,5))
plot (table(Ysim_poisson),
      xlab= "Simulated predator abundance values",
      ylab= "Frequencies")

## Negative Binomial with SPDE
NSim <- 1000
SimData_NB <- inla.posterior.sample(n = NSim, result = NB.sim)

N  <- nrow(Datos)
Ysim_nb <- matrix(nrow = N, ncol = NSim)
mu.i_nb <- matrix(nrow = N, ncol = NSim)
k.i_nb <- vector (length = NSim)

Xm <- as.matrix(X)
Am <- as.matrix(A2)

for (i in 1:NSim){
  Betas <- SimData_NB[[i]]$latent[RowNum.Betas]
  wk    <- SimData_NB[[i]]$latent[RowNum.w]
  FixedPart   <- Xm %*% Betas
  SpatialPart <- Am %*% wk
  mu.i_nb[,i]    <- exp(FixedPart + SpatialPart)
  k.i_nb [i]    <- exp(SimData_NB[[i]]$logdens$hyperpar)
  Ysim_nb[,i]    <- rnegbin(n = N, mu = mu.i_nb[,i], theta = k.i_nb[i])
}

table(Ysim_nb)
plot (table(Ysim_nb),
      xlab= "Simulated predator abundance values",
      ylab= "Frequencies")

########################################################

# OVERDISPERSION
######## Poisson
Dispersion_poisson <- vector(length = NSim) #Create space


# Calculate the dispersion for each simulated data set
for(i in 1:NSim){
  E_poisson <- (Ysim_poisson[,i] - mu.i_poisson[,i]) / sqrt(mu.i_poisson[,i])
  Dispersion_poisson[i] <- sum(E_poisson^2) # / (N - p)
}

Low95 <- max(which(sort(Dispersion_poisson) < sort(Dispersion_poisson)[0.025*length(Dispersion_poisson)]))
High95 <- min(which(sort(Dispersion_poisson) > sort(Dispersion_poisson)[0.975*length(Dispersion_poisson)]))

# Plot this as a table
hist(Dispersion_poisson,
     xlab = "Dispersion",
     ylab = "Frequency",
     breaks = 1000)

abline(v= sort(Dispersion_poisson)[Low95], col="red")
abline(v= sort(Dispersion_poisson)[High95], col="red")

# And visualize the dispersion for the original NB GLM SPDE
mu_pois_observed <- Poisson$summary.fitted.values[1:18, "mean"]
E_pois_observed <- (Datos$Potentials - mu_pois_observed) / sqrt (mu_pois_observed)
Dispersion_pois_obs <- sum(E_pois_observed^2) # / (N - p)

points(x = Dispersion_pois_obs, 
       y = 0, 
       pch = 16, 
       cex = 1, 
       col = 2)

# The red dot is the underdispersion in the original data set.
length(which(Dispersion_pois_obs > Dispersion_poisson))/length(Dispersion_poisson)

# GGPLOT2
C95 <- data.frame(Low95= sort(Dispersion_poisson)[Low95], High95= sort(Dispersion_poisson)[High95])
Dispersion_poisson <- data.frame(Dispersion_poisson= Dispersion_poisson)
Dispersion_pois_obs <- data.frame(Dispersion_pois_obs= Dispersion_pois_obs)

windows()
ggplot() + geom_histogram(data= Dispersion_poisson, aes(x= Dispersion_poisson, y= ..density..), color= "grey", alpha= 0.1) +
  geom_line(data= Dispersion_poisson, aes(x= Dispersion_poisson), size= 0.8, stat= "density") + 
  geom_vline(xintercept = C95$Low95, size = 0.7, colour = "#FF3721") +
  geom_vline(xintercept = C95$High95, size = 0.7, colour = "#FF3721") +
  geom_point(data= Dispersion_pois_obs, aes(x= Dispersion_pois_obs, y= 0), size= 3) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab ("Dispersion")

windows()
ggplot() + geom_histogram(data= Dispersion_poisson, aes(x= Dispersion_poisson, y= ..density..), color= "grey", alpha= 0.1) +
  geom_line(data= Dispersion_poisson, aes(x= Dispersion_poisson), size= 1, stat= "density") + 
   geom_point(data= Dispersion_pois_obs, aes(x= Dispersion_pois_obs, y= 0), size= 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab ("Dispersion")

Dispersion_pois_obs <- sum(E_pois_observed^2) # / (N - p)
windows()
ggplot() + geom_histogram(data= Dispersion_poisson, aes(x= Dispersion_poisson, y= ..density..), color= "grey", alpha= 0.1) +
  #geom_line(data= Dispersion_poisson, aes(x= Dispersion_poisson), size= 0.8, stat= "density") + 
  geom_vline(xintercept = Dispersion_pois_obs, size = 1, linetype= "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(family = "Times", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text = element_text(family = "Times",size=14)) + 
  ylab("Relative Frequency")+
  xlab ("Dispersion") + ggtitle("P(D|Data > D|model) = 0.969")
#+  xlim(c(0,100))

windows()
ggplot() + geom_density(data= Dispersion_poisson, aes(x= Dispersion_poisson), color="#E69F00", fill="#E69F00", position="stack", alpha= 0.8, size= 1) +
  #geom_line(data= Dispersion_poisson, aes(x= Dispersion_poisson), size= 0.8, stat= "density") + 
  geom_vline(xintercept = Dispersion_pois_obs, size = 1, linetype= "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, family = "Times", size= 18, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.x = element_text(family = "Times", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text = element_text(family = "Times",size=14),) + 
  ylab("Relative Frequency")+
  xlab ("Dispersion") +
  ggtitle("P(D|data > D|model) = 0.962") +  xlim(c(0,100))

# Save width 1000 and height 550

#####################################################################
#### Negative Binomial
Dispersion_nb <- vector(length = NSim) #Create space

# Calculate the dispersion for each simulated data set
for(i in 1:NSim){
  E_nb <- (Ysim_nb[,i] - mu.i_nb[,i]) / sqrt(mu.i_nb[,i] + mu.i_nb[,i]^2/k.pm)
  Dispersion_nb[i] <- sum(E_nb^2) # / (N - p)
}

Low95 <- max(which(sort(Dispersion_nb) < sort(Dispersion_nb)[0.025*length(Dispersion_nb)]))
High95 <- min(which(sort(Dispersion_nb) > sort(Dispersion_nb)[0.975*length(Dispersion_nb)]))

# Plot this as a table
range(Dispersion_nb)
hist(Dispersion_nb,
     xlab = "Dispersion",
     ylab = "Frequency",
     xlim = c(0,500),
     breaks = 100000)

abline(v= sort(Dispersion_nb)[Low95], col="red")
abline(v= sort(Dispersion_nb)[High95], col="red")

# And visualize the dispersion for the original NB GLM SPDE
k.pd_obs <- NB$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm_obs <- inla.emarginal(function(x) x, k.pd_obs)
k.pm_obs

mu_nb_observed <- NB$summary.fitted.values[1:18, "mean"]
E_nb_observed <- (Datos$Potentials - mu_nb_observed) / sqrt (mu_nb_observed + mu_nb_observed^2/k.pm_obs)
Dispersion_nb_obs <- sum(E_nb_observed^2) # / (N - p)
Dispersion_nb_obs

points(x = Dispersion_nb_obs, 
       y = 0, 
       pch = 16, 
       cex = 3, 
       col = 2)

# The red dot is the overdispersion in the original data set.
length(which(Dispersion_nb_obs > Dispersion_nb))/length(Dispersion_nb)
#How often is this value larger than a simulated value
mean(Dispersion_nb_obs > Dispersion_nb)

## GGPLOT 2
C95 <- data.frame(Low95= sort(Dispersion_nb)[Low95], High95= sort(Dispersion_nb)[High95])
Dispersion_nb <- data.frame(Dispersion_nb= Dispersion_nb)
Dispersion_nb_obs <- data.frame(Dispersion_nb_obs= Dispersion_nb_obs)

windows()
ggplot() + geom_histogram(data= Dispersion_nb, aes(x= Dispersion_nb, y= ..density..), color= "grey", alpha= 0.1) +
  geom_line(data= Dispersion_nb, aes(x= Dispersion_nb), size= 1, stat= "density") + 
  geom_vline(xintercept = C95$Low95, size = 0.5, colour = "#FF3721") +
  geom_vline(xintercept = C95$High95, size = 0.5, colour = "#FF3721") +
  geom_point(data= Dispersion_pois_obs, aes(x= Dispersion_pois_obs, y= 0), size= 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab ("Dispersion") + xlim(0, 500)

Dispersion_nb_obs <- sum(E_nb_observed^2) # / (N - p)
windows()
ggplot() + geom_histogram(data= Dispersion_nb, aes(x= Dispersion_nb, y= ..density..), color= "grey", alpha= 0.1) +
  #geom_line(data= Dispersion_poisson, aes(x= Dispersion_poisson), size= 0.8, stat= "density") + 
  geom_vline(xintercept = Dispersion_nb_obs, size = 1, linetype= "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(family = "Times", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text = element_text(family = "Times",size=14)) + 
  ylab("Relative Frequency")+
  xlab ("Dispersion") + 
  xlim(0,500)

windows()
ggplot() + geom_density(data= Dispersion_nb, aes(x= Dispersion_nb), color="#E69F00", fill="#E69F00", position="stack", alpha= 0.8, size= 1) +
  #geom_line(data= Dispersion_poisson, aes(x= Dispersion_poisson), size= 0.8, stat= "density") + 
  geom_vline(xintercept = Dispersion_nb_obs, size = 1, linetype= "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, family = "Times", size= 18, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.x = element_text(family = "Times", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text = element_text(family = "Times",size=14),) + 
  ylab("Relative Frequency")+
  xlab ("Dispersion") + ggtitle("P(D|data > D|model) = 0.51") +
  xlim (c(0,500))

########################################################
### How well is predicting the model??

Ysim <- Ysim_poisson
Ysim <- Ysim_nb

# RUN FROM HERE.....
windows()
Z <- matrix(nrow = max(Ysim)+1, ncol = NSim)
for (i in 1: NSim){
  zi <- table(Ysim[,i])
  I <- as.numeric(names(zi)) + 1
  Z[I,i] <- zi
}

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
Z[is.na(Z)] <- 0
Xi <- 0: max(Ysim)
AverageTable <- rowSums(Z) / NSim

AverageTable <- rowSums(Z) / NSim
apply(Z, 1, mean)
SD <- apply(Z, 1, sd)

plot(x = Xi, 
     y = AverageTable,
     type = "h",
     lwd = 5,
     xlab = "Simulated number of potential predators",
     ylab = "Frequencies",
     ylim = c(0, 20),
     xlim= c(0,20))

# Add SDs
Zs <- table(Datos$Potentials)
nx <- length(Zs)
NamesZ <- as.numeric(names(Zs))
nx <- length(NamesZ)

for (i in 1:nx){
  segments(x0 = NamesZ[i],
           x1 = NamesZ[i],
           y0 = AverageTable[i],
           y1 = AverageTable[i] + SD[i],
           lwd = 15,
           col = 1)
}


#And add the table for the observed data
for (i in 1:nx){
  segments(x0 = NamesZ[i] + 0.2,
           x1 = NamesZ[i] + 0.2,
           y0 = 0,
           y1 = Zs[i],
           lwd = 2,
           col = 2)
}
# TO HERE....		
# The red bars represent the frequency table
# for the observed data. The black lines is
# the average frequency table for the 1000
# simulated data sets. They should hopefully 
# match. If they do then the model simulates
# data that is comparable to the observed data.		

# It is not perfect. We are not predicting enough
# zeros, and too many ones.

#############################################
########################################################
Ysim <- Ysim_poisson
Ysim <- Ysim_nb
# Now that we have 1000 simulated data sets from the 
# model, what else shall we do with these simulated 
# data sets?

# We could calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim[,i] == 0)
}

table(zeros)

# From the 1,000 simulated data sets, in 2 simulated
# data sets we had 1141 zeros. In 2 simulated data sets
# we had 1148 zeros,etc......
# Your results will be different as mine.

#Let's plot this as a table
plot(table(zeros), 
     #axes = FALSE,
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Number of zeros in 1000 simulated data sets",
     xlim = c(9, 19),
     main = "Simulation results")
points(x = sum(Datos$Potentials == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)

# The red dot is the number of zeros in the original data set.
# The data simulated from the Poisson model
# contains to many zeros.

#########################################
####################################################
# We finally present the spatial component, the wks. 
# Their posterior mean values can be obtained via
w.pm.NB   <- I2.poisson$summary.random$w$mean  

# This is a vector of length 745 by 1. Each value 
# in w.pm belongs to a specific vertex on mesh 5. 
# We can either obtain the coordinates of the 745 
# vertices (LocMesh), match these with w.pm 
# (in the correct order) and use existing graphical 
# functions to plot the spatial random field, or we 
# can use INLA functions to do this for us. 
# We will go for the second approach.  



# This function is modified code from material on Haakon Bakka's website
# We will not explain what is inside this function. Just run it.
PlotField <- function(field, mesh, ContourMap, xlim, ylim, Add=FALSE, ...){
  stopifnot(length(field) == mesh$n)
  # Plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim <- ContourMap@bbox[1, ] 
  if (missing(ylim)) ylim <- ContourMap@bbox[2, ]
  
  # inla.mesh.projector: it creates a lattice using the mesh and specified ranges. 
  proj <- inla.mesh.projector(mesh, 
                              xlim = xlim, 
                              ylim = ylim, 
                              dims = c(300, 300))
  # The function inla.mesh.project can then 
  # be used to project the w's on this grid.
  field.proj <- inla.mesh.project(proj, field)
  
  # And plot the whole thing
  image.plot(list(x = proj$x, 
                  y = proj$y,
                  z = field.proj), 
             xlim = xlim, 
             ylim = ylim,
             asp = 1,
             add = Add,
             ...)  
}




# Plot the spatial random field 
library(fields)
par (mfrow = c(1,1), mar= c(3,3,3,3))

PlotField(field = w.pm.NB, mesh = mesh, xlim = range(mesh$loc[,1]), ylim = range(mesh$loc[,2]))

# Add the sampling locations (in UTM)
points(x = Loc[,1],
       y = Loc[,2], 
       cex = 0.5, 
       col = "black", 
       pch = 16)

# Posterior mean of the spatial random field
w.pm <- I2.nb$summary.random$w$mean
exp(I2.nb$summary.random$w$mean)

w.proj <- inla.mesh.projector(mesh)

w.pm100_100 <- inla.mesh.project(w.proj, w.pm)

grid <- expand.grid(x= w.proj$x/1000,
                    y= w.proj$y/1000)
grid$z <- as.vector(w.pm100_100)
levelplot(z ~ x * y,
          data= grid,
          scales= list(draw= TRUE),
          xlab= list("Easting(km)", cex= 1.5),
          ylab= list("Northing(km)", cex= 1.5),
          main= list("posterior mean spatial random field", cex= 1.5),
          panel= function(...){
            panel.levelplot(...)
            grid.points (x= Loc[,1]/1000,
                         y= Loc[,2]/1000,
                         pch= 1,
                         size= unit(0.5, "char"))
          })

#########################################################################################################
#########################################################################################################
#########################################################################################################
##### *B* AVIAN PREDATORS

Datos <- read_excel("Data_DEP_AEREOSmodf_20190212.xlsx")
str(Datos)

Datos$WF <- ifelse (Datos$fWF == "NO", 0, 1)
Datos$fWF <- factor (Datos$WF, levels = c(0,1), labels = c("No", "Yes"))
str(Datos)

## Standardized variables
Datos$Suelo_no_edificado_estd <- scale (Datos$Suelo_no_edificado)
Datos$Crop_surface_estd <- scale (Datos$Crop_surface)
Datos$Tree_surface_estd <- scale (Datos$Tree_surface)

#########################################################
# Start INLA

# TP_i ~ Poisson(mu_i)
# E(TP_i) = mu_i 
# var(TP_i) = mu_i

#          Sex_i + Location_i + Length_i + Location_i x Length_i
# mu_i = e
M1 <- inla (AvianDEP_Tot ~ fWF + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE),
            family = "poisson",
            data = Datos)
summary(M1)

# Open landscape species
M1 <- inla (AvianDEP_Tot_Open ~ fWF + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE),
            family = "poisson",
            data = Datos)
summary(M1)

# Forest species
M1 <- inla (AvianDEP_Tot_Forest ~ fWF + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE),
            family = "poisson",
            data = Datos)
summary(M1)

# LANDSCAPE
M1 <- inla (AvianDEP_Tot ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE),
            family = "poisson",
            data = Datos)
summary(M1)

# Open landscape species
M1 <- inla (AvianDEP_Tot_Open ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE),
            family = "poisson",
            data = Datos)
summary(M1)

# Forest species
M1 <- inla (AvianDEP_Tot_Forest ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE),
            family = "poisson",
            data = Datos)
summary(M1)

#Check for overdispersion frequentist style
mu1 <- M1$summary.fitted.values[,"mean"]
E1  <- (Datos$AvianDEP_Tot - mu1) / sqrt(mu1)
N   <- nrow(Datos) # + 1 In negative binomial due to the dispersion parameter K
p   <- nrow(M1$summary.fixed)
Dispersion <- sum(E1^2) / (N - p)
Dispersion
# In a frequentist analysis this would indicate 
# severe overdispersion.

# But it would be nice to come to the conclusion 
# that the model is overdispersed using Bayesian 
# tools.

############################################
# We will now implement the 7-step protocol
# for assessing whether the Poisson GLM
# is over- or underdispersed.


# Step 1: Apply the model in INLA
# This may sound obvious as we already 
# executed the Poisson GLM in INLA, but we 
# need to do it again, with a small modification 
# this time. The config = TRUE option allows 
# us to simulate regression parameters in the next 
# step.

M2 <- inla(AvianDEP_Tot ~ fWF + offset(log(Esfuerzo_horas)),
           control.compute=list(config = TRUE, dic = TRUE),
           family = "poisson",
           data = Datos)

# Open landscape species
M2 <- inla(AvianDEP_Tot_Open ~ fWF + offset(log(Esfuerzo_horas)),
           control.compute=list(config = TRUE, dic = TRUE),
           family = "poisson",
           data = Datos)

# Forest species
M2 <- inla(AvianDEP_Tot_Forest ~ fWF + offset(log(Esfuerzo_horas)),
           control.compute=list(config = TRUE, dic = TRUE),
           family = "poisson",
           data = Datos)

## LANDSCAPE
M2 <- inla (AvianDEP_Tot ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)),
            control.compute=list(config = TRUE, dic = TRUE),
            family = "poisson",
            data = Datos)


# Step 2: Simulate regression parameters
# We use the function inla.posterior.sample 
# to simulate from the model. The output is 
# stored in the Sim object.

set.seed(12345)
Sim <- inla.posterior.sample(n = 1, result = M2)

Sim[[1]]$latent
# This gives an object with 162 rows for this 
# specific data set and model. The first 18 
# rows are simulated values for eta = X * beta, 
# where X is the matrix with covariates, and 
# the last 2 rows are simulated regression parameters. 
# This is just one set of simulated values 
# (due to n = 1 in the function above). 

# Step 3: Calculate predicted values
# We have multiple options to do this step. 
# We can either take the first 155 rows and 
# exponent them to get mu = exp(eta), 
# or we access the last 7 rows and calculate 
# the fitted values via mu = exp(X * beta).
RowNum <- 19:20  # Last 2 rows are the betas
RowNum <- 19:22

X      <- model.matrix(~fWF + offset(log(Esfuerzo_horas)), data = Datos)
X      <- model.matrix(~Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)), data = Datos)

Betas  <- Sim[[1]]$latent[RowNum]
mu     <- exp(X %*% Betas)

# Step 4: Simulate count data
# We use the rpois function for this.

Ysim <- rpois(n = nrow(Datos), lambda = mu)

# We now have 18 simulated values that 
# correspond to the observed covariate values. 
# We can calculate how many zeros this simulated 
# data set has, determine the maximum value, 
# the minimum value, the range, etc.

# Step 5: Calculate summary statistic
# Instead of the numbers of zeros, or the 
# maximum value, we can also calculate the 
# Pearson residuals for the simulated data, 
# and square and sum them.
Es <- (Ysim - mu) /sqrt(mu) # Predicted - Mean
sum(Es^2)

# By the way, the sum of squared Pearson 
# residuals for the original data and model is:
E1 <- (Datos$AvianDEP_Tot - mu1) /sqrt(mu1)
SS <- sum(E1^2)
SS    #bigger!

# Step 6: Repeat steps 2 to 5 thousand times
#Repeating the simulation step a large 
# number of times requires only minor modification 
# of the code presented in step 2

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = M2)

# Now we have thousand simulated mean values 
# from the model. Processing the information 
# requires more coding. We first determine the 
# number of rows in the data (N), and then 
# create a matrix Ysim that can hold thousand 
# simulated data sets (one in each column of Ysim). 
# Then we start a loop and in each iteration we 
# extract the simulated betas and predict the 
# response variable.

N    <- nrow(Datos)
Ysim   <- matrix(nrow = N, ncol = NSim)
mu.sim <- matrix(nrow = N, ncol = NSim)

for (i in 1:NSim){
  Betas <- SimData[[i]]$latent[RowNum]
  mu.sim[,i] <- exp(X %*% Betas)
  Ysim[,i]   <- rpois(n = N, lambda = mu.sim [,i])
}

# We now have thousand simulated data sets with 
# numbers of potential predators.

# Step 7: Compare simulation results and observed data
# In this step we compare the simulation results 
# with the observed data. There are many things that 
# we can do with the simulated data sets. 
# For example, for each simulated data set we can 
# count the number of zeros, and see whether the 
# percentages of zeros are comparable to the percentage 
# of zeros in the observed data. If the model does 
# not produce enough zeros then we can start thinking 
# about zero inflated models. However, zero inflation 
# is not an issue for this data set.

# The sum of squared Pearson residuals 
# for the observed data is calculated as follows.
E1 <- (Datos$AvianDEP_Tot - mu1) /sqrt(mu1)
SS <- sum(E1^2)
SS

# We calculate the same statistic for 
# each of the simulated data sets.

# POISSON
SS.sim <- NULL
for(i in 1:NSim){
  e2 <- (Ysim[,i] - mu.sim[,i]) /sqrt(mu.sim[,i])
  SS.sim[i] <- sum(e2^2) 
}
sum(SS > SS.sim) / NSim 

# We would like to have a value here saying 50%, 40%
# If the result is 10% then we would have UNDERDISPERSION
# If the resulta is 90%,100$, then we have OVERDISPERSION
# Admisible values: 50%, 60&, 70%, 40%, 30%
# Cuases? No correct distribution, missing variables, etc.

# If it is 100%: The sum of squared Pearson residuals for the 
# observed data is always larger than that of 
# the simulate data sets, and this indicates that 
# the variation in the observed data is larger 
# than is allowed for by a Poisson distribution. 
# And that is called overdispersion.

#########################################
# What else can we do with the 1,000 simulated data sets?
# We can also calculate the dispersion statistic for each
# simulated data set

#RUN FROM HERE........
p        <- length(Betas)
Disp.Sim <- vector(length = NSim) #Create space

# Calculate the dispersion for each simulated data set

# POISSON
for(i in 1:NSim){
  e2 <- (Ysim[,i] - mu.sim[,i]) /sqrt(mu.sim[,i])
  Disp.Sim[i] <- sum(e2^2) / (N - p)
}

# Plot this as a table
hist(Disp.Sim,
     xlab = "Dispersion",
     ylab = "Frequency",
     xlim = c(0.6, 5),
     breaks = 25)

# And visualize the dispersion for the original Poisson GLMM
points(x = Dispersion, 
       y = 0, 
       pch = 16, 
       cex = 3, 
       col = 2)

# TO HERE ............
# The red dot is the overdispersion in the original data set.
# We have serious overdispersion!
########################################



########################
# Do we need zero inflated models?
# Calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim[,i] == 0)
}
table(zeros)

#Let's plot this as a table
plot(table(zeros), 
     #axes = FALSE,
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Number of zeros in 1000 simulated data sets",
     xlim = c(5, 19),
     main = "Simulation results")
points(x = sum(Datos$AvianDEP_Tot == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#The red dot is the number of zeros in the original data set.
#The data simulated from the Poisson model
# Zero inflation is not an issue here.

#########################################

###############################
#NB GLM
M3 <- inla(AvianDEP_Tot ~ fWF + offset(log(Esfuerzo_horas)),
           control.compute=list(config = TRUE, dic = TRUE),
           family = "nbinomial",
           control.predictor = list(#link = 1,
              compute = TRUE),
           quantiles = c(0.025, 0.975),
           data = Datos)

M3 <- inla(AvianDEP_Tot ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)),
           control.compute=list(config = TRUE, dic = TRUE),
           family = "nbinomial",
           quantiles = c(0.025, 0.975),
           data = Datos)

summary(M3)
# Task: Write down the fitted model

# Totalparasites_i ~ NB(mu_i, k)
# E(Totalparasites_i) = mu_i
# var(Totalparasites_i) = mu_i + mu_i^2 / k
# log(mu_i) = Covariate stuff


# We go on with model M3 (the NB GLM in which k is
# estimated).            
summary(M3)            

# We can get some info on the hyperparameters:
# Posterior mean of k
k.pd <- M3$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm <- inla.emarginal(function(x) x, k.pd)
k.pm


# We can also plot the posterior 
# distribution of k. Because we get 
# the k in marginals.hyperpar there is
# no need for things like: sigma = 1/sqrt(tau)
# Therefore the info in $marginals.hyperpar
# is similar to that in tmarginal:

k.marg <- inla.tmarginal(function(x) x, k.pd)

par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
plot(k.pd, 
     type = "l",
     xlim = c(0,10),
     xlab = expression(paste(k)),
     ylab = expression(paste("Pr(", k," | data)")))
text(1.1, 2.2, "A", cex = 1.5)

plot(k.marg, 
     type = "l",
     xlim = c(0,100),
     xlab = expression(paste(k)),
     ylab = expression(paste("Pr(", k," | data)")))
text(1.1, 2.2, "B", cex = 1.5)
# Second one is just a little bit smoother.

########################################################################################################
### MODEL VALIDATION NBINOMIAL

# We should still check whether the NB GLM
# is under- or overdispersed.
# We can do it the classical way:
mu3 <- M3$summary.fitted.values[1:18,"mean"]
E3  <- (Datos$AvianDEP_Tot - mu3) /sqrt(mu3 + mu3^2 / k.pm)
N <- nrow(Datos)
p <- nrow(M3$summary.fixed) + 1 # 4 Betas (modelo suelos) 2 Betas (modelo WF), and 1 is the hyperparameter k
Dispersion <- sum(E3^2) / (N - p)
Dispersion  

# We can also do the simulation study again
# to check for over- or underdispersion.

# The code below is copy-paste from above. 
set.seed(12345)

# Model WF
RowNum <- 19:20
X <- model.matrix(~fWF + offset(log(Esfuerzo_horas)), data = Datos)

# Model Usos
RowNum <- 19:22  #These are the rows we want
X <- model.matrix(~Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)), data = Datos)

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = M3)
N      <- nrow(Datos)
Ysim   <- matrix(nrow = N, ncol = NSim)
mu.sim <- matrix(nrow = N, ncol = NSim)

for (i in 1: NSim){
  k <- exp(SimData[[i]]$logdens$hyperpar)
  Betas <- SimData[[i]]$latent[RowNum]
  mu.sim[,i] <- exp(X %*% Betas)
  Ysim[,i] <- rnegbin(n = N, # Negative Binomial distributed simulated data
                      mu = mu.sim[,i],
                      theta = k)
}

# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?

#RUN FROM HERE........
Disp.Sim <- vector(length = NSim) #Create space

# Calculate the dispersion for each simulated data set
for(i in 1:NSim){
  e3 <- (Ysim[,i] - mu.sim[,i]) /sqrt(mu.sim[,i] + mu.sim[,i]^2 / k.pm)
  Disp.Sim[i] <- sum(e3^2)# / (N - p)
}

# This is the summary statistic for the observed data:
SS  <- sum(E3^2) 
SS

#How often is this value larger than a simulated value
mean(SS > Disp.Sim)

# Plot this as a table
par(mfrow = c(1,1))
hist(Disp.Sim,
     xlab = "Dispersion",
     ylab = "Frequency",
     xlim = c(0, 20),
     breaks = 100000)

# And visualize the dispersion for the original Poisson GLMM
points(x = Dispersion, 
       y = 0, 
       pch = 16, 
       cex = 3, 
       col = 2)
# This means that the NB distribution complies with the data.

# Do we need ZERO inflated models?
# Calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim[,i] == 0)
}

#Let's plot this as a table
par (mfrow = c(1,1))
plot(table(zeros), 
     #axes = FALSE,
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Number of zeros in 1000 simulated data sets",
     xlim = c(7, 19),
     main = "Simulation results")
points(x = sum(Datos$AvianDEP_Tot == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#The red dot is the number of zeros in the original data set.
#The data simulated from the Poisson model
# Zero inflation is not an issue here.

############################################
# Model validation of the NB GLM
# Get fitted values and Pearson residuals
mu3 <- M3$summary.fitted.values[1:18,"mean"]

# And the rest is a matter of plotting it all. 
p <- ggplot(Datos, aes(x=fWF, y=AvianDEP_Tot)) + 
  geom_boxplot()
# Box plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=2)
# Box plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0), size = 5) + xlab ("Wind farms presence") + ylab ("Predator abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(family = "Times", size = 20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 20, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text = element_text(family = "Times",size=18))


# Negative Binomial
E3  <- (Datos$AvianDEP_Tot - mu3) / sqrt(mu3 + mu3^2 / k.pm)

# Outliers?
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = mu3, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     
# No clear outliers.
# These bands of lines are due to the 
# discrete nature of the data.

# Plot residuls vs each covariate
Datos$E3 <- E3 # From the Negative Binomial
MyVar <- c("fWF", "Suelo_no_edificado", "Crop_surface", "Tree_surface")
MyMultipanel.ggp2(Z=Datos,
                  varx = MyVar,
                  vary = "E3",
                  ylab = "Pearson Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

## MORAN'S I TEST
library (spdep)
library (ncf)
Loc <- cbind(Datos$XCOORD/1000, Datos$YCOORD/1000)
D <- dist(Loc)

w <- 1/as.matrix(dist(Loc))
diag(w) <- 0

moran.test (E3, mat2listw(w))

leadI <- spline.correlog(x=Loc[,1], y=Loc[,2],
                         z=E3, resamp=100, quiet=TRUE)

par (mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
par (mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(density(E3))
plot (leadI)

## VARIOGRAM
MyData <- data.frame(E3 = E3, 
                     X  = Datos$XCOORD/1000, 
                     Y  = Datos$YCOORD/1000)
coordinates(MyData)  <- c("X", "Y")
V1 <- variogram(E3 ~ 1, 
                MyData, 
                cressie = TRUE)

#Plot both variograms in one graph
p <- ggplot()
p <- p + geom_point(data = V1,
                    aes(x = dist, 
                        y = gamma)
)
p <- p + geom_line(data = V1,
                   aes(x = dist, 
                       y = gamma),
                   col = "red")

p <- p + xlab("Distance") + ylab("Semi-variogram")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p <- p + ylim(0,1)
p

# A different GRAPH
plot(V1, 
     main = "", 
     xlab = list(label = "Distance", cex = 1.5), 
     ylab = list(label = "Semi-variogram", cex = 1.5),
     pch = 16,
     col = 1,
     cex = 1.5
)

# CONCLUSION: MORAN'S I does not show Spatial Autocorrelation, but the Variograms are not trully nice

### PLOT RESIDUALS!
# Option 1: Use different point sizes and symbols
#           based on the values of the residuals.
E3
MyCex <- 3 * abs(E3) / max(E3) + 0.5
Sign  <- as.numeric(E3 >=0) + 1
MyPch <- c(1, 16)[Sign]

utmcoor <- SpatialPoints (coords = cbind(Datos$XCOORD, Datos$YCOORD), proj4string = CRS("+proj=utm +zone=30N"))
longlat <- spTransform(utmcoor, CRS ("+proj=longlat"))
Datos$Longitude <- coordinates (longlat)[,1]
Datos$Latitude <- coordinates (longlat)[,2]

xyplot (Longitude ~ Latitude,
        data = Datos, aspect = "iso",
        cex = MyCex, pch = MyPch, col = 1, 
        xlab = list(label = "Easting", cex = 1.5),
        ylab = list(label = "Northing", cex = 1.5))

############################################################################################################
############################################################################################################
# MODEL SELECTION with the NB

k.pm #Posterior mean value of k

M5 <- inla(AvianDEP_Tot ~ fWF + offset(log(Esfuerzo_horas)),
           control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
           family = "nbinomial",
           data = Datos)

M5 <- inla(AvianDEP_Tot ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)),
           control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
           family = "nbinomial",
           data = Datos)
# TASK: Write down the model that we just fitted.

#Drop fWF
M5a <- inla(AvianDEP_Tot ~  1 + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)

M5b <- inla(AvianDEP_Tot ~  Crop_surface_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)
M5c <- inla(Potentials ~  Suelo_no_edificado_estd + Tree_surface_estd + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)
M5d <- inla(AvianDEP_Tot ~  Suelo_no_edificado_estd + Crop_surface_estd + offset(log(Esfuerzo_horas)),
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)

dic  <- c(M5$dic$dic, M5a$dic$dic)
dic  <- c(M5$dic$dic, M5b$dic$dic, M5c$dic$dic, M5d$dic$dic)
waic <- c(M5$waic$waic, M5a$waic$waic)
waic <- c(M5$waic$waic, M5b$waic$waic, M5c$waic$waic, M5d$waic$waic)
Z <- cbind(dic, waic)
rownames(Z) <- c("fWF", "Null")
rownames(Z) <- c("Full", "Drop Suelo no edificado", "Drop Crops", "Drop Trees")
Z

## According to DIC --> Best model FULL MODEL
## According to wAIC --> Best model the one that drops crop surface!

M5$summary.fixed
M5c$summary.fixed

M5e <- inla(Potentials ~  Tree_surface_estd,
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)
M5f <- inla(Potentials ~  Suelo_no_edificado_estd,
            control.compute = list(dic = TRUE, waic = TRUE),
            family = "nbinomial",
            data = Datos)

dic  <- c(M5c$dic$dic, M5e$dic$dic, M5f$dic$dic)
waic <- c(M5c$waic$waic, M5e$waic$waic, M5f$waic$waic)

Z <- cbind(dic, waic)
rownames(Z) <- c("Both", "Drop Suelo no edificado", "Drop Trees")
Z

M5c$summary.fixed
M5f$summary.fixed

########################################################################################################
############## PRESENTING THE RESULTS!
M5$summary.fixed

# Sketch the model fit.
# The last thing we do in this analysis 
# is to make a visualization of the optimal 
# model. We discussed in detail how to do 
# this in Section 9.8. We will use the 
# inla.make.lincombs approach. Recall that 
# these were the steps that we have to execute.

# 1. Create a grid of covariate values.
# 2. Make a design matrix using the model.matrix function.
# 3. Run the model in INLA with the inla.make.lincombs option.
# 4. Visualize the results using ggplot2.

# We follow the four steps. We first make a grid of 
# covariate values.


# 1. Create a grid of covariate values.
MyData <- ddply(Datos, 
                .(fWF))

MyData <- data.frame (Suelo_no_edificado = seq(min(Datos$Suelo_no_edificado), 
                                               max(Datos$Suelo_no_edificado), 
                                               length = 25),
                      Crop_surface = seq(min(Datos$Crop_surface), 
                                         max(Datos$Crop_surface), 
                                         length = 25),
                      Tree_surface = seq(min(Datos$Tree_surface), 
                                         max(Datos$Tree_surface), 
                                         length = 25))

MyData <- data.frame (Suelo_no_edificado = seq(min(Datos$Suelo_no_edificado), 
                                               max(Datos$Suelo_no_edificado), 
                                               length = 25),
                      Crop_surface = rep(mean (Datos$Crop_surface), 
                                         length = 25),
                      Tree_surface = rep(mean (Datos$Tree_surface), 
                                         length = 25),
                      Crop_surface_estd = rep(mean (Datos$Crop_surface_estd), 
                                              length = 25),
                      Tree_surface_estd = rep(mean (Datos$Tree_surface_estd), 
                                              length = 25))

MyData <- data.frame (Suelo_no_edificado = seq(min(Datos$Suelo_no_edificado), 
                                               max(Datos$Suelo_no_edificado), 
                                               length = 25))

MyData$Suelo_no_edificado_estd <- scale (MyData$Suelo_no_edificado)
Crop_surface_estd <- scale (MyData$Crop_surface)
MyData$Tree_surface_estd <- scale (MyData$Tree_surface)

head(MyData)

# 2. Make a design matrix
Xmat <- model.matrix(~ fWF,
                     data = MyData)

Xmat <- model.matrix(~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd,
                     data = MyData)

Xmat <- model.matrix(~ Suelo_no_edificado_estd,
                     data = MyData)
head(Xmat)
Xmat <- as.data.frame(Xmat)


# 3. Run the model in INLA with the 
#    inla.make.lincombs option
lcb <- inla.make.lincombs(Xmat)
Hyper.NB <- list(size = list(initial = 1, fixed = TRUE))

M7 <- inla(Potentials ~ fWF,
           lincomb = lcb,
           control.inla = list(lincomb.derived.only = FALSE),
           control.predictor = list(#link = 1,
             compute = TRUE),
           family = "nbinomial",
           control.family = list(hyper = Hyper.NB),
           data = Datos)

M7 <- inla(Potentials ~ Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd,
           lincomb = lcb,
           control.inla = list(lincomb.derived.only = FALSE),
           control.predictor = list(#link = 1,
             compute = TRUE),
           family = "nbinomial",
           control.family = list(hyper = Hyper.NB),
           data = Datos)

M7 <- inla(Potentials ~ Suelo_no_edificado_estd,
           lincomb = lcb,
           control.inla = list(lincomb.derived.only = FALSE),
           control.predictor = list(#link = 1,
             compute = TRUE),
           family = "nbinomial",
           control.family = list(hyper = Hyper.NB),
           data = Datos)

# 4. Visualize the results using ggplot2.

# get the predicted values
Pred.pm <- M7$summary.lincomb[,c("mean", 
                                 "0.025quant", 
                                 "0.975quant")]

print(head(M7$summary.lincomb.derived), digits = 2)
Out1 <- M7$summary.lincomb.derived

# The output above is on the log-scale; 
# the exponent has not been taken yet. 
# We could do this ourselves via:

MyData$mu.wrong   <- exp(Out1$'mean')
MyData$selo.wrong <- exp(Out1$'0.025quant')
MyData$seup.wrong <- exp(Out1$'0.975quant')

# and now we could run the ggplot2 
# code from Section 9.8 again. But strictly 
# speaking this is not correct because 
# E(exp(x)) !=  exp(E(x)) if x is a stochastic 
# variable. And we are doing a Bayesian analysis, 
# so x is stochastic. The ‘E()’ stands for the 
# expected value and exp for the exponential function. 
# This was not a problem in Chapter 9 with the 
# identity link function, but with the log link 
# function it is.

# The solution is to use support functions 
# from INLA that convert the distribution of x 
# into the distribution of exp(x). Starting point 
# is the $marginals.lincomb.derived, which 
# contains the marginal posterior distribution for 
# each predicted value, and we have 75 of them.

M7$marginals.lincomb.derived
Pred.marg <- M7$marginals.lincomb.derived

# It is slightly confusing to access the 75 
# marginal posterior distributions. The marginal 
# posterior distribution for the first predicted 
# value is obtained via:
Pred.marg[[1]] 

# This contains an entire posterior distribution, 
# which can be plotted if you want to. What we 
# really want is its posterior mean (or median) 
# and a 95% credible interval. That is obtained with
inla.qmarginal(c(0.025, 0.5, 0.975), 
               inla.tmarginal(exp, Pred.marg[[1]]))

# For the posterior mean use:
inla.emarginal(exp, Pred.marg[[1]])

# Here is some code that does it for all 75 points
# at once:
MyData$mu <- unlist( 
  lapply(
    Pred.marg,
    function(x) inla.emarginal(exp,x)))

MyData$selo <- unlist( 
  lapply(
    Pred.marg,
    function(x) 
      inla.qmarginal(c(0.025), 
                     inla.tmarginal(exp, x))
  ))

MyData$seup <- unlist( 
  lapply(
    Pred.marg,
    function(x) 
      inla.qmarginal(c(0.975), 
                     inla.tmarginal(exp, x))
  ))

# If you don’t like fancy coding (like us), 
# just run the following loop. It repeats the 
# steps that we just explained for each predicted 
# value. Using a loop means that we have to think 
# less. The results are exactly the same.

for (i in 1:18){
  MyData$mu2[i]  <- inla.emarginal(exp, Pred.marg[[i]])
  lo.up <- inla.qmarginal(c(0.025, 0.975), 
                          inla.tmarginal(exp, Pred.marg[[i]]))
  MyData$selo2[i] <- lo.up[1]
  MyData$seup2[i] <- lo.up[2]    	
}               

MyData

# And the rest is a matter of plotting it all. 
p <- ggplot(Datos, aes(x=fWF, y=Potentials)) + 
  geom_boxplot()
# Box plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=2)
# Box plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0), size = 5) + xlab ("Wind farms presence") + ylab ("Predator abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(family = "Times", size = 20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 20, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text = element_text(family = "Times",size=18))


######################################################################################################
# And the rest is a matter of plotting it all. 
p <- ggplot()
p <- p + geom_point(data = Datos, 
                    aes(y = Potentials, x = Suelo_no_edificado),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Suelo no edificado") + ylab("Potential predators")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_line(data = MyData, 
                   aes(x = Suelo_no_edificado, 
                       y = mu), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Suelo_no_edificado, 
                         ymax = seup, 
                         ymin = selo),
                     alpha = 0.2)
p

p <- ggplot()
p <- p + geom_point(data = Datos, 
                    aes(y = Potentials, x = Tree_surface),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Tree surface") + ylab("Potential predators")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_line(data = MyData, 
                   aes(x = Tree_surface, 
                       y = mu), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Tree_surface, 
                         ymax = seup, 
                         ymin = selo),
                     alpha = 0.2)
p

range(MyData$Tree_surface)
range(Datos$Tree_surface)
range(MyData$Suelo_no_edificado)


#########################################################################################################
##### *B.2* OVERDISPERSION! INCORPORATE SPATIAL AUTOCORRELATION
library(lattice)
library(sp)
library(raster)
library(dismo)
library(splancs)
library(INLA)
library(reshape)
library(gstat)
library(ggplot2)
library(ggmap)
library(rworldmap)
library(rgdal)
library(rgeos)

##################################
# Transform coordinates from UTM to LONG/LAT
utmcoor <- SpatialPoints (coords = cbind(Datos$XCOORD, Datos$YCOORD), proj4string = CRS("+proj=utm +zone=30N"))
longlat <- spTransform(utmcoor, CRS ("+proj=longlat"))
Datos$Longitude <- coordinates (longlat)[,1]
Datos$Latitude <- coordinates (longlat)[,2]


####################################################
# We will implement the following 8 steps.
# 1. Make a mesh.
# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
# 3. Define the SPDE.
# 4. Define the spatial field.
# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.
# 6. Specify the model formula in terms of the 
#    response variable, covariates and the 
#    spatial correlated term.
# 7. Run the spatial model in INLA.
# 8. Inspect the results.



########################
#1. Make a mesh.
#   Step 1 of making a mesh:  Get a 
#   sense for the distribution of 
#   distances between sampling locations. 

# Done before

# We need a grid on top of our sampling points
Loc <- cbind(Datos$XCOORD/1000, Datos$YCOORD/1000)
Bound <- inla.nonconvex.hull(Loc, convex= -0.15)
mesh    <- inla.mesh.2d(boundary = Bound, 
                        max.edge=  c(1,2) # Maxima distancia permitida
)
mesh$n

# Bellow different script without outer part
# Recommended settings 
# See Section 3.1 at: https://haakonbakka.bitbucket.io/btopic104.html
# mesh    <- inla.mesh.2d(boundary = ConvHull, 
# max.edge=  c(1.5))
# max.edge: maximum allowed triangle edge lengths in 
#           the inner domain and in the outer extension
# cutoff: minimum allowed distance between points. Points 
#         at a closer distance than the supplied value are 
#         replaced by a single vertex

par(mfrow=c(1,1), mar=c(0,0,2,0))
plot(mesh)
points(Loc, col = 1, pch = 16, cex = 1)
mesh$n #We can cope with a finer mesh.

#########################################
# Step 2. Define the weighting factors a_ik (also called 
#         the projector matrix).
A2 <- inla.spde.make.A(mesh, loc = Loc)
dim(A2)  #18 observations on a 2485 grid
#        A is a weight matrix
head(A2)
############################################



############################################
# Step 3. Define the SPDE.
#### Basic Style
spde <- inla.spde2.matern(mesh, alpha=2)

# Use PC priors
# P(Range < range0) = alpha  and P(sigma > sigma0) = alpha

# We are going to use:
# P(Range < 1) = 0.5  # DIFFUSE PRIOR: No idea...it can be anything..though we think the range is small-ish
# and P(sigma > 1) = 0.05

# The 1km comes from our desire to avoid overfitting.
# The sigma....let's assume that the covariates are not important
M1 <- glm(AvianDEP_Tot ~ 1 , data = Datos, family = "poisson")
summary(M1)

# This means
# E[NEggs] = exp(Intercept) = exp(0.8232) = 2.27

range(Datos$AvianDEP_Tot)
# How big do the spatial correlated residuals u_iu need to be 
# if we want to cover 0 to 8 and:

# E[Potentials] = exp(0.2007 + u)
#               = exp(0.2007) *  exp(u)
#               = 2.27* exp(u)

# The exp(u) should not be larger than 4-ish.
# Hence, the u should not be not bigger than log(6) = 1.4-ish.
# Because we assume u_i ~ N(0, sigma^2), using sigma = 1 should do 
# the job.

# P(sigma > 1) = 0.05

# Below an uninformative prior: most likely the range is larger than 1 km
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(1, 0.05), 
                            prior.sigma = c(1, 0.05))

# Diffuse prior. No Idea, but we think that the range is small
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(1, 0.5), 
                            prior.sigma = c(1, 0.05))

# Diffuse prior. No Idea, but we think that the range is small <--- WE USE THIS ONE
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(7, 0.5), 
                            prior.sigma = c(3, 0.05))
##########################################

##########################################
# Step 4. Define the spatial field.
# Next we set up a list for the spatial random 
# intercept u. As explained in Section 13.6, 
# u is rewritten internally as A * w. We need 
# to specify the w in INLA. This is done with 
# the inla.spde.make.index function

# The size of the mesh is: 
mesh$n

# This number is also in
spde$n.spde

# It is also the number of w_k values that we will get.
# For this mesh we use:
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde,
  n.group = 1,
  n.rep1 = 1)

str(w.index)
# Here is where spatial-temporal models can be defines.
# We will do that later. 
#####################################

#####################################
# Step 5.	Make a stack. 

# This is a confusing step. We discussed it
# in detail in Chapter 13.

# Make the X matrix

# Set up the model. 
# Create a data frame with an intercept and covariates.
N <- nrow(Datos)

# WIND FARMS - Factors
Xm <- model.matrix(~ fWF + offset(log(Esfuerzo_horas)), data= Datos)
X <- data.frame(Intercept = Xm[,1],
                fWF = Xm[,2])

# LAND USES
X <- data.frame(Intercept = rep (1,N),
                Suelo_no_edificado_estd = Datos$Suelo_no_edificado_estd,
                Crop_surface_estd = Datos$Crop_surface_estd,
                Tree_surface_estd = Datos$Tree_surface_estd)

X <- as.data.frame(X) #Avoids problems

# And here is the stack.
StackFit <- inla.stack(
  tag  = "Fit",
  data = list(y = Datos$AvianDEP_Tot),  
  A    = list(1, 1, A2),  # Intercept and covariates, spatial field                    
  effects = list( 
    X = X, #Covariates
    Esfuerzo_horas = Datos$Esfuerzo_horas, # The offset
    w = w.index          #Spatial field  
  )) 

# Open landscape species
StackFit <- inla.stack(
  tag  = "Fit",
  data = list(y = Datos$AvianDEP_Tot_Open),  
  A    = list(1, 1, A2),  # Intercept and covariates, spatial field                    
  effects = list( 
    X = X, #Covariates
    Esfuerzo_horas = Datos$Esfuerzo_horas, # The offset
    w = w.index          #Spatial field  
  )) 

# Forestry species
StackFit <- inla.stack(
  tag  = "Fit",
  data = list(y = Datos$AvianDEP_Tot_Forest),  
  A    = list(1, 1, A2),  # Intercept and covariates, spatial field                    
  effects = list( 
    X = X, #Covariates
    Esfuerzo_horas = Datos$Esfuerzo_horas, # The offset
    w = w.index          #Spatial field  
  )) 
#############################################
#6.	Specify the model formula in terms of the 
#   response variable, covariates and the 
#   spatial correlated term.

# These are the models that we will fit:
# Y_i ~ Poisson(mu_i)
# E(Y_i) = mu_i
# Model 1: log(mu_i) = Covariate stuff
# Model 2: log(mu_i) = Covariate stuff + u_i

f1 <- y ~ -1 + Intercept + fWF + offset (log(Esfuerzo_horas))
f1 <- y ~ -1 + Intercept + Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset (log(Esfuerzo_horas))

f2 <- y ~ -1 + Intercept + fWF + offset (log(Esfuerzo_horas)) +
  f(w, model = spde)
f2 <- y ~ -1 + Intercept + Suelo_no_edificado_estd + Crop_surface_estd + Tree_surface_estd + offset (log(Esfuerzo_horas))+ 
  f(w, model = spde)


#############################################
# 7. Run the spatial model in INLA.
# First we run the model without spatial dependency.

I1.poisson <- inla(f1,
                   family = "poisson", 
                   data = inla.stack.data(StackFit),
                   control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
                   control.predictor = list(
                     A = inla.stack.A(StackFit)))

I1.nb <- inla(f1,
              family = "nbinomial", 
              data = inla.stack.data(StackFit),
              control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
              control.predictor = list(
                A = inla.stack.A(StackFit)))


# And this is the model with the spatial field:
#  Y_i ~ Poisson(mu_i)
#  E(Y_i)   = mu_i
#  var(Y_i) =  mu_i
#  log(mu_i) = Intercept + Covariates + u_i
#  Where u_i is spatially correlated noise.

I2.poisson <- inla(f2,
                   family = "poisson", 
                   data = inla.stack.data(StackFit),
                   control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
                   control.predictor = list(
                     A = inla.stack.A(StackFit)))

I2.nb <- inla(f2,
              family = "nbinomial", 
              data = inla.stack.data(StackFit),
              control.compute = list(config= TRUE, config = TRUE, dic = TRUE, waic = TRUE),
              control.predictor = list(
                A = inla.stack.A(StackFit)))

# Priors for the intercept and the coefficients
inla.set.control.fixed.default()
# Priors for the hyperparameters

# And compare the models with DICs and WAICs
dic  <- c(I1.poisson$dic$dic, I1.nb$dic$dic, I2.poisson$dic$dic, I2.nb$dic$dic)
waic <- c(I1.poisson$waic$waic, I1.nb$waic$waic, I2.poisson$waic$waic, I2.nb$waic$waic)
Z2     <- cbind(dic, waic)
rownames(Z2) <- c("Poisson GLM",  "NB GLM", 
                  "Poisson GLM + SPDE", "NB GLM + SPDE")
Z2

## MODEL SELECTION
dic  <- c(I1.poisson$dic$dic, I1a.poisson$dic$dic, I1b.poisson$dic$dic, 
          I1.nb$dic$dic, I1a.nb$dic$dic, I1b.nb$dic$dic,
          I2.poisson$dic$dic, I2a.poisson$dic$dic, I2b.poisson$dic$dic,
          I2.nb$dic$dic, I2a.nb$dic$dic, I2b.nb$dic$dic)

waic <- c(I1.poisson$waic$waic, I1a.poisson$waic$waic, I1b.poisson$waic$waic,
          I1.nb$waic$waic, I1a.nb$waic$waic, I1b.nb$waic$waic,
          I2.poisson$waic$waic, I2a.poisson$waic$waic, I2b.poisson$waic$waic, 
          I2.nb$waic$waic, I2a.nb$waic$waic, I2b.nb$waic$waic)
Z     <- cbind(dic, waic)
rownames(Z) <- c("Poisson GLM Full",  "Poisson GLM Drop Crop", "Poisson GLM Only Suelo",
                 "NB GLM Full",  "NB GLM Drop Crop", "NB GLM Only Suelo", 
                 "Poisson GLM + SPDE Full", "Poisson GLM + SPDE Drop Crop", "Poisson GLM + SPDE Only Suelo",
                 "NB GLM + SPDE Full", "NB GLM + SPDE Drop Crop", "NB GLM + SPDE Only Suelo")
Z


## dic     waic
## Poisson GLM        70.82523 74.42231
## NB GLM             70.74562 73.89878
## Poisson GLM + SPDE 71.00877 73.68832
## NB GLM + SPDE      70.61799 71.69746

#######################################
#8. Inspect the results of I2

# Fixed parameters
# Here are the results of the model without and with the
# spatial random effects

# Fixed parameters
# Here are the results of the model without and with the
# spatial random effects

I1.poisson$summary.fixed[, c("mean", "0.025quant", "0.975quant")]
I2.poisson$summary.fixed[, c("mean", "0.025quant", "0.975quant")]

I2.nb$summary.fixed[, c("mean", "0.025quant", "0.975quant")]

# Let's plot the results of the model, without and with
# the spatial correlation side by side.
# Better not try to understand all this R code.
# RUN FROM HERE......

Combined <- rbind(I1.poisson$summary.fixed[, c("mean", "0.025quant", "0.975quant")],
                  I2.poisson$summary.fixed[, c("mean", "0.025quant", "0.975quant")])

Combined$WhichModel <- rep(c("Non-spatial", "Spatial"), 
                           each = nrow(I2.poisson$summary.fixed))

rownames(I2.poisson$summary.fixed)

Combined$WhichVariable <- rep(c("Intercept", "Wind farms"), 2)
Combined$WhichVariable <- rep(c("Intercept", "S.Gravel roads", "S.Crops", "S.Trees"), 2)

colnames(Combined) <- c("Mean", "Lo", "Up", "Model", "Predictor")
positions <- unique(Combined$Predictor)

windows()
ggplot(Combined, aes(y=Mean, x= Predictor)) +
  geom_errorbar(
    aes(ymin = Lo, ymax = Up, color = Model),
    position = position_dodge(0.3), width = 0.2, size= 1) +
  geom_hline(yintercept = 0, size = 0.7, linetype = "dashed") +
  geom_point(aes(color = Model), position = position_dodge(0.3), size= 2) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme(legend.position="top",
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(family = "Times", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text = element_text(family = "Times",size=14)) + 
  scale_x_discrete(limits = positions)
# TO HERE ......


############################################################################
## HYPERPARAMETERS
I2 <- I2.poisson
I2 <- I2.nb

SpatField.w <- inla.spde2.result(inla = I2,
                                 name = "w",
                                 spde = spde,
                                 do.transfer = TRUE)

Kappa <- inla.emarginal(function(x) x, 
                        SpatField.w$marginals.kappa[[1]] )

Sigma_u <- inla.emarginal(function(x) sqrt(x), 
                          SpatField.w$marginals.variance.nominal[[1]] )

Range <- inla.emarginal(function(x) x, 
                        SpatField.w$marginals.range.nominal[[1]] )

Kappa
Sigma_u
Range       #Distance at which the correlation diminishes

I2$summary.hyperpar
# This is perhaps a nicer graph to make and present.
# Show correlation structure
# First we obtain the locations of each point of the mesh.
LocMesh <- mesh$loc[,1:2]

# And then we calculate the distance between each vertex.
D <- as.matrix(dist(LocMesh))

# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, max(D), length = 100)      
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1) 
Cor.M[1] <- 1

# Which we plot here:
par(mfrow=c(1,1), mar = c(5,5,2,2))
plot(x = d.vec, 
     y = Cor.M, 
     pch = 16, 
     type = "l", 
     cex.lab = 1.5,
     xlab = "Distance", 
     ylab = "Correlation",
     xlim = c(0, 30))
abline(h = 0.1, lty = 2)

# Negative Binomial: posterior mean of k
# As its name implies, the negative binomial shape parameter, k, 
# describes the shape of a negative binomial distribution. In other words, 
# k is only a reasonable measure to the extent that your data represent a 
# negative binomial distribution. The broader the k, The broader the distribution
NB <- I1.nb
NB <- I2.nb

k.pd <- NB$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm <- inla.emarginal(function(x) x, k.pd)
k.pm

############################################################################################################
##### MODEL VALIDATION
#### *1* Plot fitted values versus observed data
#### *2* Check the Pearson residuals for any remaining spatial dependency - Variogram and Moran's I
#### *3* PLOT RESIDUALS
#### *4* RESIDUALS VERSUS EACH COVARIATE
#### *5* Proportion of zeros and overdispersion

########################################################
Poisson <- I1.poisson # Poisson
NB <- I1.nb # Negative Binomial

Poisson <- I2.poisson # Poisson + SPDE
NB <- I2.nb # Negative Binomial + SPDE

########################################################
#### *1* Plot fitted values versus observed data
# Poisson
mu_poisson <- Poisson$summary.fitted.values[1:18, "mean"]
plot(x = mu_poisson,
     y = Datos$AvianDEP_Tot_Open)

# Negative Binomial
mu_nb <- NB$summary.fitted.values[1:18, "mean"]
plot(x = mu_nb,
     y = Datos$AvianDEP_Tot)
# Hmm. Looks too good?
# Are we overfitting?

########################################################
#### *2* Check the Pearson residuals for any remaining spatial dependency.
# Variogram

# Poisson
E_poisson <- (Datos$AvianDEP_Tot - mu_poisson) / sqrt (mu_poisson) # Poisson Pearson Residuals

MyData <- data.frame(E1 = E_poisson, 
                     X  = Datos$XCOORD/1000, 
                     Y  = Datos$YCOORD/1000)
coordinates(MyData)  <- c("X", "Y")
V_poisson <- variogram(E_poisson ~ 1, 
                       MyData, 
                       cressie = TRUE)


# Negative Binomial
k.pd <- NB$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm <- inla.emarginal(function(x) x, k.pd)
k.pm

E_nb <- (Datos$AvianDEP_Tot - mu_nb) / sqrt (mu_nb + mu_nb^2/k.pm) # Negative BinomialPearson Residuals

MyData <- data.frame(E1 = E_nb, 
                     X  = Datos$XCOORD/1000, 
                     Y  = Datos$YCOORD/1000)
coordinates(MyData)  <- c("X", "Y")
V_nb <- variogram(E_nb ~ 1, 
                  MyData, 
                  cressie = TRUE)

# Plot both variograms in one graph
p <- ggplot()
p <- p + geom_point(data = V_poisson,
                    aes(x = dist, 
                        y = gamma)
)
p <- p + geom_line(data = V_poisson,
                   aes(x = dist, 
                       y = gamma),
                   col = "red")

p <- p + geom_point(data = V_nb,
                    aes(x = dist, 
                        y = gamma)
)
p <- p + geom_line(data = V_nb,
                   aes(x = dist, 
                       y = gamma),
                   col = "blue"
)

p <- p + xlab("Distance") + ylab("Semi-variogram")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p <- p
p
# No clear patterns for the Poisson GLM with spatial correlation
# But strange patterns for the NB GLM with spatial correlation

### MORAN'S I
## MORAN'S I TEST
E <- E_poisson
E <- E_nb

Loc <- cbind(Datos$XCOORD/1000, Datos$YCOORD/1000)
D <- dist(Loc)
w <- 1/as.matrix(dist(Loc))
diag(w) <- 0

moran.test (E, mat2listw(w))
leadI <- spline.correlog(x=Loc[,1], y=Loc[,2],
                         z=E, resamp=100, quiet=TRUE)
par (mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
plot(density(E))
plot (leadI)

########################################################
#### 3* PLOT RESIDUALS
# Option 1: Use different point sizes and symbols
#           based on the values of the residuals.
E <- E_poisson
E <- E_nb

MyCex <- 3 * abs(E) / max(E) + 0.5
Sign  <- as.numeric(E >=0) + 1
MyPch <- c(1, 16)[Sign]

xyplot (Longitude ~ Latitude,
        data = Datos, aspect = "iso",
        cex = MyCex, pch = MyPch, col = 1, 
        xlab = list(label = "Easting", cex = 1.5),
        ylab = list(label = "Northing", cex = 1.5))

# Patterns in the residuals are present! 
# Black dots is for positive residuals (left side of the island) whereas white dots is for negative residuals (right side of the island)

# If you see patterns then we have to stop! Standard errors, t values and estimates are not true -> They are biased
# Standard errors are faluty, and then also de t-value and also the p-values (Domino!)
# If the p-value is really low (0.0000) the effects are going to be low. But if the p-values are near 0.01 then we have troubles!

########################################################
#### *4* RESIDUALS VERSUS EACH COVARIATE
# Plot residuls vs each covariate
Datos$E <- E_poisson # From the Poisson
Datos$E <- E_nb # From the Negative Binomial

MyVar <- c("fWF", "Suelo_no_edificado", "Crop_surface", "Tree_surface")
MyMultipanel.ggp2(Z=Datos,
                  varx = MyVar,
                  vary = "E",
                  ylab = "Pearson Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

########################################################
### *5* Proportion of zeros and overdispersion

# Simulating from the model

# In the Turco parasite example we implemented a
# 7-step protocol for simulating from the model.
# See also Powerpoint presentation: 
#    P15_SpatTempBook_Chapter10_V1.pdf


# We will now implement the 7-step protocol
# for assessing whether the Poisson GLM
# is over- or underdispersed.

# Step 1: Apply the model in INLA
# This may sound obvious as we already 
# executed the Poisson GLM in INLA, but we 
# need to do it again, with a small modification 
# this time. The config = TRUE option allows 
# us to simulate regression parameters in the next 
# step.

# WITHOUT SPDE
# Poisson
Poisson.sim <- I1.poisson
# Negative Binomial
NB.sim <- I1.nb

#################
# WITH SPDE
# Poisson
Poisson.sim <- I2.poisson
# Negative Binomial
NB.sim <- I2.nb


###########

k.pd <- NB.sim$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm <- inla.emarginal(function(x) x, k.pd)
k.pm
# Step 2: Simulate regression parameters
# We use the function inla.posterior.sample 
# to simulate from the model. The output is 
# stored in the Sim object.

# Here is some fancy code to grap the 1  betas (the same for both models)
Sim <- inla.posterior.sample(n = 1, result = Poisson.sim) # Simulate one dataset from the poisson

MyParams <- rownames(Poisson$summary.fixed)
RowNum.Betas <- lapply (MyParams, 
                        function (x)
                          grep (x,
                                rownames (Sim[[1]]$latent),
                                fixed = TRUE))
RowNum.Betas <- as.numeric (RowNum.Betas)
RowNum.Betas # In this rows are our parameters!

# In case of spatial autocorrelation
# Here is some fancy code to grap the w 
Nw <- mesh$n 
# The w are labelled as w:1 to w:Nw   where Nw is the number of vertices
MyParams <- paste("w", seq(from = 1, to = Nw), sep = ":")
MyID <- function(x){ which(rownames(Sim[[1]]$latent) == x) }
RowNum.w <- lapply(MyParams, MyID)
RowNum.w <- as.numeric(RowNum.w)
RowNum.w

############# Simulate 1000 datasets
## Poisson without SPDE
NSim <- 1000
SimData_Poisson <- inla.posterior.sample(n = NSim, result = Poisson.sim)

N  <- nrow(Datos)
Ysim_poisson <- matrix(nrow = N, ncol = NSim)
mu.i_poisson <- matrix(nrow = N, ncol = NSim)

Xm <- as.matrix(X)

for (i in 1:NSim){
  Betas <- SimData_Poisson[[i]]$latent[RowNum.Betas]
  FixedPart   <- Xm %*% Betas
  mu.i_poisson[,i]    <- exp(FixedPart)
  Ysim_poisson[,i]    <- rpois(n = nrow(Datos), lambda = mu.i_poisson[,i])
}

table(Ysim_poisson)
par(mfrow= c(1,1), mar= c(5,5,5,5))
plot (table(Ysim_poisson),
      xlab= "Simulated predator abundance values",
      ylab= "Frequencies")

## Negative Binomial without SPDE
NSim <- 1000
SimData_NB <- inla.posterior.sample(n = NSim, result = NB.sim)

N  <- nrow(Datos)
Ysim_nb <- matrix(nrow = N, ncol = NSim)
mu.i_nb <- matrix(nrow = N, ncol = NSim)
k.i_nb <- vector (length = NSim)

Xm <- as.matrix(X)

for (i in 1:NSim){
  Betas <- SimData_NB[[i]]$latent[RowNum.Betas]
  FixedPart   <- Xm %*% Betas
  mu.i_nb[,i]    <- exp(FixedPart)
  k.i_nb [i]    <- exp(SimData_NB[[i]]$logdens$hyperpar)
  Ysim_nb[,i]    <- rnegbin(n = nrow(Datos), mu = mu.i_nb[,i], theta = k.i_nb[i])
}

table(Ysim_nb)
plot (table(Ysim_nb),
      xlab= "Simulated predator abundance values",
      ylab= "Frequencies")


#######################################
## Poisson with SPDE
NSim <- 1000
SimData_Poisson <- inla.posterior.sample(n = NSim, result = Poisson.sim)

N  <- nrow(Datos)
Ysim_poisson <- matrix(nrow = N, ncol = NSim)
mu.i_poisson <- matrix(nrow = N, ncol = NSim)

Xm <- as.matrix(X)
Am <- as.matrix(A2)

for (i in 1:NSim){
  Betas <- SimData_Poisson[[i]]$latent[RowNum.Betas]
  wk    <- SimData_Poisson[[i]]$latent[RowNum.w]
  FixedPart   <- Xm %*% Betas
  SpatialPart <- Am %*% wk
  mu.i_poisson[,i]    <- exp(FixedPart + SpatialPart)
  Ysim_poisson[,i]    <- rpois(n = nrow(Datos), lambda = mu.i_poisson[,i])
}

table(Ysim_poisson)
par(mfrow= c(1,1), mar= c(5,5,5,5))
plot (table(Ysim_poisson),
      xlab= "Simulated predator abundance values",
      ylab= "Frequencies")

## Negative Binomial with SPDE
NSim <- 1000
SimData_NB <- inla.posterior.sample(n = NSim, result = NB.sim)

N  <- nrow(Datos)
Ysim_nb <- matrix(nrow = N, ncol = NSim)
mu.i_nb <- matrix(nrow = N, ncol = NSim)
k.i_nb <- vector (length = NSim)

Xm <- as.matrix(X)
Am <- as.matrix(A2)

for (i in 1:NSim){
  Betas <- SimData_NB[[i]]$latent[RowNum.Betas]
  wk    <- SimData_NB[[i]]$latent[RowNum.w]
  FixedPart   <- Xm %*% Betas
  SpatialPart <- Am %*% wk
  mu.i_nb[,i]    <- exp(FixedPart + SpatialPart)
  k.i_nb [i]    <- exp(SimData_NB[[i]]$logdens$hyperpar)
  Ysim_nb[,i]    <- rnegbin(n = nrow(Datos), mu = mu.i_nb[,i], theta = k.i_nb[i])
}

table(Ysim_nb)
par(mfrow= c(1,1))
plot (table(Ysim_nb),
      xlab= "Simulated predator abundance values",
      ylab= "Frequencies")

########################################################

# OVERDISPERSION
######## Poisson
Dispersion_poisson <- vector(length = NSim) #Create space

# Calculate the dispersion for each simulated data set
for(i in 1:NSim){
  E_poisson <- (Ysim_poisson[,i] - mu.i_poisson[,i]) / sqrt(mu.i_poisson[,i])
  Dispersion_poisson[i] <- sum(E_poisson^2) # / (N - p)
}

Low95 <- max(which(sort(Dispersion_poisson) < sort(Dispersion_poisson)[0.025*length(Dispersion_poisson)]))
High95 <- min(which(sort(Dispersion_poisson) > sort(Dispersion_poisson)[0.975*length(Dispersion_poisson)]))

# Plot this as a table
hist(Dispersion_poisson,
     xlab = "Dispersion",
     ylab = "Frequency",
     breaks = 1000)

abline(v= sort(Dispersion_poisson)[Low95], col="red")
abline(v= sort(Dispersion_poisson)[High95], col="red")

# And visualize the dispersion for the original NB GLM SPDE
mu_pois_observed <- Poisson$summary.fitted.values[1:18, "mean"]
E_pois_observed <- (Datos$AvianDEP_Tot - mu_pois_observed) / sqrt (mu_pois_observed)
Dispersion_pois_obs <- sum(E_pois_observed^2) # / (N - p)

points(x = Dispersion_pois_obs, 
       y = 0, 
       pch = 16, 
       cex = 1, 
       col = 2)

# The red dot is the underdispersion in the original data set.
length(which(Dispersion_pois_obs > Dispersion_poisson))/1000

# GGPLOT2
C95 <- data.frame(Low95= sort(Dispersion_poisson)[Low95], High95= sort(Dispersion_poisson)[High95])
Dispersion_poisson <- data.frame(Dispersion_poisson= Dispersion_poisson)

windows()
ggplot() + geom_density(data= Dispersion_poisson, aes(x= Dispersion_poisson), color="#E69F00", fill="#E69F00", position="stack", alpha= 0.8, size= 1) +
  #geom_line(data= Dispersion_poisson, aes(x= Dispersion_poisson), size= 0.8, stat= "density") + 
  geom_vline(xintercept = Dispersion_pois_obs, size = 1, linetype= "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, family = "Times", size= 18, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.x = element_text(family = "Times", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text = element_text(family = "Times",size=14),) + 
  ylab("Relative Frequency")+
  xlab ("Dispersion") +
  ggtitle("P(D|data > D|model) = 0.431") +  xlim(c(0,200))


#####################################################################
#### Negative Binomial
Dispersion_nb <- vector(length = NSim) #Create space

# Calculate the dispersion for each simulated data set
for(i in 1:NSim){
  E_nb <- (Ysim_nb[,i] - mu.i_nb[,i]) / sqrt(mu.i_nb[,i] + mu.i_nb[,i]^2/k.pm)
  Dispersion_nb[i] <- sum(E_nb^2) # / (N - p)
}

Low95 <- max(which(sort(Dispersion_nb) < sort(Dispersion_nb)[0.025*length(Dispersion_nb)]))
High95 <- min(which(sort(Dispersion_nb) > sort(Dispersion_nb)[0.975*length(Dispersion_nb)]))

# Plot this as a table
range(Dispersion_nb)
hist(Dispersion_nb,
     xlab = "Dispersion",
     ylab = "Frequency",
     breaks = 100000)

# And visualize the dispersion for the original NB GLM SPDE
k.pd_obs <- NB$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
k.pm_obs <- inla.emarginal(function(x) x, k.pd_obs)
k.pm_obs

mu_nb_observed <- NB$summary.fitted.values[1:18, "mean"]
E_nb_observed <- (Datos$AvianDEP_Tot - mu_nb_observed) / sqrt (mu_nb_observed + mu_nb_observed^2/k.pm_obs)
Dispersion_nb_obs <- sum(E_nb_observed^2) # / (N - p)
Dispersion_nb_obs

points(x = Dispersion_nb_obs, 
       y = 0, 
       pch = 16, 
       cex = 3, 
       col = 2)

# The red dot is the overdispersion in the original data set.
length(which(Dispersion_nb > Dispersion_nb_obs))/length(Dispersion_nb)
#How often is this value larger than a simulated value
mean(Dispersion_nb_obs > Dispersion_nb)

## GGPLOT 2
C95 <- data.frame(Low95= sort(Dispersion_nb)[Low95], High95= sort(Dispersion_nb)[High95])
Dispersion_nb <- data.frame(Dispersion_nb= Dispersion_nb)

ggplot() + geom_histogram(data= Dispersion_nb, aes(x= Dispersion_nb, y= ..density..), color= "grey", alpha= 0.1) +
  geom_line(data= Dispersion_nb, aes(x= Dispersion_nb), size= 1, stat= "density") + 
  geom_vline(xintercept = C95$Low95, size = 0.5, colour = "#FF3721") +
  geom_vline(xintercept = C95$High95, size = 0.5, colour = "#FF3721") +
  geom_point(data= Dispersion_pois_obs, aes(x= Dispersion_pois_obs, y= 0), size= 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab ("Dispersion") + xlim(0, 100)

# GGPLOT2
C95 <- data.frame(Low95= sort(Dispersion_nb)[Low95], High95= sort(Dispersion_nb)[High95])
Dispersion_nb <- data.frame(Dispersion_nb= Dispersion_nb)

windows()
ggplot() + geom_density(data= Dispersion_nb, aes(x= Dispersion_nb), color="#E69F00", fill="#E69F00", position="stack", alpha= 0.8, size= 1) +
  #geom_line(data= Dispersion_poisson, aes(x= Dispersion_poisson), size= 0.8, stat= "density") + 
  geom_vline(xintercept = Dispersion_nb_obs, size = 1, linetype= "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, family = "Times", size= 18, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.x = element_text(family = "Times", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Times", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text = element_text(family = "Times",size=14),) + 
  ylab("Relative Frequency")+
  xlab ("Dispersion") +
  ggtitle("P(D|data > D|model) = 0.494") +  xlim(c(0,100))

########################################################
# Do we need ZERO inflated models?
# Calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim_nb[,i] == 0)
}

#Let's plot this as a table
par (mfrow = c(1,1))
plot(table(zeros), 
     #axes = FALSE,
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Number of zeros in 1000 simulated data sets",
     xlim = c(0, 20),
     main = "Simulation results")
points(x = sum(Datos$AvianDEP_Tot_Forest == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#The red dot is the number of zeros in the original data set.
#The data simulated from the Poisson model
# Zero inflation is not an issue here.

#######################################################3
### How well is predicting the model??

Ysim <- Ysim_poisson
Ysim <- Ysim_nb

# RUN FROM HERE.....
Z <- matrix(nrow = max(Ysim)+1, ncol = NSim)
for (i in 1: NSim){
  zi <- table(Ysim[,i])
  I <- as.numeric(names(zi)) + 1
  Z[I,i] <- zi
}

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
Z[is.na(Z)] <- 0
Xi <- 0: max(Ysim)
AverageTable <- rowSums(Z) / NSim

AverageTable <- rowSums(Z) / NSim
apply(Z, 1, mean)
SD <- apply(Z, 1, sd)

plot(x = Xi, 
     y = AverageTable,
     type = "h",
     lwd = 5,
     xlab = "Simulated number of potential predators",
     ylab = "Frequencies",
     ylim = c(0, 20))

# Add SDs
Zs <- table(Datos$AvianDEP_Tot)
nx <- length(Zs)
NamesZ <- as.numeric(names(Zs))
nx <- length(NamesZ)

for (i in 1:nx){
  segments(x0 = NamesZ[i],
           x1 = NamesZ[i],
           y0 = AverageTable[i],
           y1 = AverageTable[i] + SD[i],
           lwd = 15,
           col = 1)
}


#And add the table for the observed data
for (i in 1:nx){
  segments(x0 = NamesZ[i] + 0.2,
           x1 = NamesZ[i] + 0.2,
           y0 = 0,
           y1 = Zs[i],
           lwd = 2,
           col = 2)
}
# TO HERE....		
# The red bars represent the frequency table
# for the observed data. The black lines is
# the average frequency table for the 1000
# simulated data sets. They should hopefully 
# match. If they do then the model simulates
# data that is comparable to the observed data.		

# It is not perfect. We are not predicting enough
# zeros, and too many ones.

#############################################
########################################################
Ysim <- Ysim_poisson
Ysim <- Ysim_nb
# Now that we have 1000 simulated data sets from the 
# model, what else shall we do with these simulated 
# data sets?

# We could calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim[,i] == 0)
}

table(zeros)

# From the 1,000 simulated data sets, in 2 simulated
# data sets we had 1141 zeros. In 2 simulated data sets
# we had 1148 zeros,etc......
# Your results will be different as mine.

#Let's plot this as a table
plot(table(zeros), 
     #axes = FALSE,
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Number of zeros in 1000 simulated data sets",
     xlim = c(5, 19),
     main = "Simulation results")
points(x = sum(Datos$AvianDEP_Tot == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)

# The red dot is the number of zeros in the original data set.
# The data simulated from the Poisson model
# contains to many zeros.

#########################################
####################################################
# We finally present the spatial component, the wks. 
# Their posterior mean values can be obtained via
w.pm.NB   <- I2.poisson$summary.random$w$mean  

# This is a vector of length 745 by 1. Each value 
# in w.pm belongs to a specific vertex on mesh 5. 
# We can either obtain the coordinates of the 745 
# vertices (LocMesh), match these with w.pm 
# (in the correct order) and use existing graphical 
# functions to plot the spatial random field, or we 
# can use INLA functions to do this for us. 
# We will go for the second approach.  



# This function is modified code from material on Haakon Bakka's website
# We will not explain what is inside this function. Just run it.
PlotField <- function(field, mesh, ContourMap, xlim, ylim, Add=FALSE, ...){
  stopifnot(length(field) == mesh$n)
  # Plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim <- ContourMap@bbox[1, ] 
  if (missing(ylim)) ylim <- ContourMap@bbox[2, ]
  
  # inla.mesh.projector: it creates a lattice using the mesh and specified ranges. 
  proj <- inla.mesh.projector(mesh, 
                              xlim = xlim, 
                              ylim = ylim, 
                              dims = c(300, 300))
  # The function inla.mesh.project can then 
  # be used to project the w's on this grid.
  field.proj <- inla.mesh.project(proj, field)
  
  # And plot the whole thing
  image.plot(list(x = proj$x, 
                  y = proj$y,
                  z = field.proj), 
             xlim = xlim, 
             ylim = ylim,
             asp = 1,
             add = Add,
             ...)  
}




# Plot the spatial random field 
library(fields)
par (mfrow = c(1,1))

PlotField(field = w.pm.NB, mesh = mesh, xlim = range(mesh$loc[,1]), ylim = range(mesh$loc[,2]))

# Add the sampling locations (in UTM)
points(x = Loc[,1],
       y = Loc[,2], 
       cex = 0.5, 
       col = "black", 
       pch = 16)

#########################################################################################################
#########################################################################################################