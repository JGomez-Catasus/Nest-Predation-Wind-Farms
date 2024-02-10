###################################################
###         PROBABILITY OF NEST PREDATION        ##
###################################################

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
source("HighstatLibV10.R")
########################################################################
Datos <- read_excel("Datos_Depredacion_2016_20200317.xlsx")
Datos <- read_excel("Datos_Depredacion_TOTAL_20190207.xlsx")

names(Datos)

#### *A* DATA PREPARATION
str (Datos)

Datos$fWF <- factor (Datos$WF, levels = c(0,1), labels = c("No", "Yes"))
Datos$fTanda <- factor (Datos$Tanda, levels = c(0,1), labels = c("Tanda1", "Tanda2"))

Datos$GLOBSuelo_no_edificado.std <- scale (Datos$GLOBSuelo_no_edificado)
Datos$GLOBCrop_surface.std <- scale (Datos$GLOBCrop_surface)
Datos$GLOBTree_surface.std <- scale (Datos$GLOBTree_surface)
Datos$GLOBOtras_construcciones.std <- scale (Datos$GLOBOtras_construcciones)
Datos$PC1ver.std <- scale (Datos$PC1ver)
Datos$PC2ver.std <- scale (Datos$PC2ver)
Datos$PC1hor.std <- scale (Datos$PC1hor)
Datos$PC2hor.std <- scale (Datos$PC2hor)
Datos$PC3hor.std <- scale (Datos$PC3hor)
Datos$PC4hor.std <- scale (Datos$PC4hor)
Datos$Excrementos_mean.std <- scale (Datos$Excrementos_mean)
Datos$Excrementos_potentials.std <- scale (Datos$Excrementos_potentials)
Datos$Avian_Predators.std <- scale (Datos$Avian_Predators)
Datos$DIST_TREE.std <- scale (Datos$DIST_TREE)
Datos$Dist_WFm.std <- scale (Datos$Dist_WFm)
Datos$DIST_cropm.std <- scale (Datos$DIST_cropm)
Datos$DIST_paths.std <- scale (Datos$DIST_paths)
Datos$Avian_Predators_residuals.std <- scale (Datos$Avian_Predators_residuals)

str(Datos)

#########################################################################################################
#########################################################################################################
#########################################################################################################

## Memory Effect

# Nests predated in the first and second period
Predated1_ID <- Datos$Nest[Datos$DEP == 1 & Datos$Tanda == 0]
length(Predated1_ID)
Predated11 <- sapply(Predated1_ID, function(x) {length(Datos$FID[Datos$Nest == x & Datos$DEP == 1 & Datos$Tanda == 1])})
Predated11
sum(Predated11)/length(Predated1_ID)

# Nests predated in the first period but not in the second period
Predated10 <- sapply(Predated1_ID, function(x) {length(Datos$FID[Datos$Nest == x & Datos$DEP == 0 & Datos$Tanda == 1])})
Predated10
sum(Predated10)/length(Predated1_ID)

# Nests predated in the second period but not in the first period
Predated0_ID <- Datos$Nest[Datos$DEP == 0 & Datos$Tanda == 0]
length(Predated0_ID)
Predated01 <- sapply(Predated0_ID, function(x) {length(Datos$FID[Datos$Nest == x & Datos$DEP == 1 & Datos$Tanda == 1])})
Predated01
sum(Predated01)/length(Predated0_ID)

TotalPredated1 <- length(Predated1_ID)
TotalPredated2 <- length(Datos$Nest[Datos$DEP == 1 & Datos$Tanda == 1])

TotalPredated1 + TotalPredated2


# Only 9 nests out of the 22 nests predated in the first period, were predated again the second period.

# Now with sations:

# Estaciones con al menos un nido depredado en el primer periodo
Predated1_ID <- Datos$Station[Datos$DEP == 1 & Datos$Tanda == 0]
length(Predated1_ID)
Predated1_ID
length(unique(Predated1_ID)) # 10 estaciones de 18 tuvieron al menos 1 nido depredado en el primer periodo

# Estaciones con al menos un nido depredado en el primero y segundo periodo
Predated11 <- unlist(sapply(unique(Predated1_ID), function(x) {unique(Datos$Station[Datos$Station == x & Datos$DEP == 1 & Datos$Tanda == 1])}))
Predated11
length(Predated11) # 7 estaciones de 10 fueron revisitadas en el segundo periodo

# 3 Estaciones tuvieron al menos un nido depredado en el primer periodo, pero ninguno en el segundo

# Estaciones con ningun nido depredado en el primer periodo y al menos uno en el segundo periodo
Predated0_ID <- sapply(unique(Datos$Station), function(x){if (length(Datos$Station[Datos$Station == x & Datos$DEP == 0 & Datos$Tanda == 0]) == length(Datos$Station[Datos$Station == x & Datos$Tanda == 0])) {x} else {1}})
Predated0_ID
length(Predated0_ID)

 X <- which(Predated0_ID != 1) 
 Predated0_ID <- as.vector(Predated0_ID[X])

Predated01 <- unlist(sapply(unique(Predated0_ID), function(x) {unique(Datos$Station[Datos$Station == x & Datos$DEP == 1 & Datos$Tanda == 1])}))
Predated01 # 2 estaciones que no presentaron ningun nido depredado en el primer periodo lo presentaron en el segundo


#########################################################
# Start INLA

# TP_i ~ Poisson(mu_i)
# E(TP_i) = mu_i 
# var(TP_i) = mu_i

#          Sex_i + Location_i + Length_i + Location_i x Length_i
# mu_i = e
## WIND FARMS
M1 <- inla (DEP ~ fWF + Tanda,
            control.compute = list(dic = TRUE),
            family = "binomial",
            Ntrials = 1,
            data = Datos)

M1a <- inla (DEP ~ fWF + f(Station, model = "iid"),
           control.compute = list(dic = TRUE),
           family = "binomial",
           Ntrials = 1,
           data = Datos)

M1b <- inla (DEP ~ fWF + f(Station, model = "iid") + f(Nest, model = "iid"),
             control.compute = list(dic = TRUE),
             family = "binomial",
             Ntrials = 1,
             data = Datos)

M1c <- inla (DEP ~ fWF*Tanda + f(Station, model = "iid") + f(Nest, model = "iid"),
             control.compute = list(dic = TRUE),
             family = "binomial",
             Ntrials = 1,
             data = Datos)

M1d <- inla (DEP ~ fWF + Tanda + f(Station, model = "iid") + f(Nest, model = "iid"),
             control.compute = list(dic = TRUE),
             family = "binomial",
             Ntrials = 1,
             data = Datos)

summary(M1)

dic  <- c(M1$dic$dic, M1a$dic$dic, M1b$dic$dic, M1c$dic$dic, M1d$dic$dic)
waic <- c(M1$waic$waic, M1a$waic$waic, M1b$waic$waic, M1c$waic$waic, M1d$waic$waic)
Z <- cbind(dic, waic)
rownames(Z) <- c("Normal Bernoulli", "Bernoulli + Station random", "Bernoulli + Station random + Nest random", "Bernoulli:Tanda + Station random + Nest random", "Bernoulli + Tanda + Station random + Nest random")
Z

M1c$summary.fixed

## MICROHABITAT
M1 <- inla (DEP ~ PC1ver.std + PC2ver.std + PC1hor.std + PC2hor.std + PC3hor.std + PC4hor.std + f(Station, model = "iid") + f(Nest, model = "iid"),
            control.compute = list(dic = TRUE),
            family = "poisson",
            data = Datos)
summary(M1)

##############################################################################################################
#### MODEL VALIDATION
M1 <- M1c

# Get Pearson residuals for Bernoulli models
# Model 1:
Pi <- M1$summary.fitted.values[1:N, "mean"]
ExpY <- Pi
VarY <- (1 - Pi) * Pi
E1 <- (Datos$DEP - ExpY) / sqrt(VarY)

## MORAN'S I TEST
Loc <- cbind(Datos$X/1000, Datos$Y/1000)
D <- dist(Loc)

w <- 1/as.matrix(dist(Loc))
diag(w) <- 0

moran.test (E1, mat2listw(w))

leadI <- spline.correlog(x=Loc[,1], y=Loc[,2],
                         z=E1, resamp=100, quiet=TRUE)

par (mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
par (mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(density(E1))
plot (leadI)

## VARIOGRAM
MyData <- data.frame(E1 = E1, 
                     X  = Datos$X/1000, 
                     Y  = Datos$Y/1000)
coordinates(MyData)  <- c("X", "Y")
V1 <- variogram(E1 ~ 1, 
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

# CONCLUSION: Variogram, Spatial Autocorrelation

### PLOT RESIDUALS!
# Option 1: Use different point sizes and symbols
#           based on the values of the residuals.
E1
MyCex <- 3 * abs(E1) / max(E1) + 0.5
Sign  <- as.numeric(E1 >=0) + 1
MyPch <- c(1, 16)[Sign]

utmcoor <- SpatialPoints (coords = cbind(Datos$X, Datos$Y), proj4string = CRS("+proj=utm +zone=30N"))
longlat <- spTransform(utmcoor, CRS ("+proj=longlat"))
Datos$Longitude <- coordinates (longlat)[,1]
Datos$Latitude <- coordinates (longlat)[,2]

xyplot (Longitude ~ Latitude,
        data = Datos, aspect = "iso",
        cex = MyCex, pch = MyPch, col = 1, 
        xlab = list(label = "Easting", cex = 1.5),
        ylab = list(label = "Northing", cex = 1.5))

### RESIDUALS VERSUS EACH COVARIATE
# Plot residuls vs each covariate
Datos$E1 <- E1 # From the Negative Binomial
MyVar <- c("fWF", "Tanda", "GLOBSuelo_no_edificado", "GLOBCrop_surface", "GLOBTree_surface")
MyMultipanel.ggp2(Z=Datos,
                  varx = MyVar,
                  vary = "E1",
                  ylab = "Pearson Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

########################################################################################################
############## PRESENTING THE RESULTS!
M1$summary.fixed

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
utmcoor <- SpatialPoints (coords = cbind(Datos$X, Datos$Y), proj4string = CRS("+proj=utm +zone=30N"))
longlat <- spTransform(utmcoor, CRS ("+proj=longlat"))
Datos$Longitude <- coordinates (longlat)[,1]
Datos$Latitude <- coordinates (longlat)[,2]

range (Datos$Longitude)
range (Datos$Latitude)

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

Loc <- cbind(Datos$X/1000, Datos$Y/1000)
D <- dist(Loc)
max(D)
mean(D)
0.5*diff(range(Datos$Y))

# Create the function.
Datos2 <- Datos[Datos$fTanda == "Tanda1",]
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(as.vector(D))

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
Bound <- inla.nonconvex.hull(Loc, convex= -0.15)

mesh  <- inla.mesh.2d(boundary = Bound, 
                      max.edge=  c(0.2,2) # Maxima distancia permitida
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

#########################################
# Step 2. Define the weighting factors a_ik (also called 
#         the projector matrix).
A2 <- inla.spde.make.A(mesh, loc = Loc)
dim(A2)  # 312 observations on a 2485 grid
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
M1 <- glm(DEP ~ 1 , data = Datos, family = "binomial")
summary(M1)

# This means
# E[Predation] = exp(Intercept) / (1 + exp(Intercept)) = exp(-1.5198) / (1 + exp(-1.5198)) = 0.1794

range(Datos$Potentials)
# How big do the spatial correlated residuals u_iu need to be 
# if we want to cover 0 to 6 and:

# E[Potentials] = exp(-1.5198 + u) / (1+ exp(-1.5198 + u))
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
#*#
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

# Option 2: When we have factors!
## CAREFUL!! Dont introduce here the random effects!!
Xm <- model.matrix(~ fWF + fTanda + PC1hor.std + PC2hor.std + PC3hor.std + PC4hor.std + PC1ver.std + 
                     PC2ver.std + GLOBSuelo_no_edificado.std + GLOBCrop_surface.std + GLOBTree_surface.std +
                     GLOBOtras_construcciones.std + Avian_Predators.std + Excrementos_potentials.std +
                     DIST_TREE.std + DIST_cropm.std + DIST_paths.std + Dist_WFm.std, data= Datos)[,-1]

X <- data.frame(Intercept = rep(1,N),
                WF = Xm[,1],
                Tanda = Xm[,2],
                PC1hor.std = Xm[,3],
                PC2hor.std = Xm[,4],
                PC3hor.std = Xm[,5],
                PC4hor.std = Xm[,6],
                PC1ver.std = Xm[,7],
                PC2ver.std = Xm[,8],
                GLOBSuelo_no_edificado.std = Xm[,9],
                GLOBCrop_surface.std = Xm[,10],
                GLOBTree_surface.std = Xm[,11],
                GLOBOtras_construcciones.std = Xm[,12],
                Avian_Predators.std = Xm[,13],
                Excrementos_potentials.std = Xm[,14],
                DIST_TREE.std = Xm[,15],
                DIST_cropm.std = Xm[,16],
                DIST_paths.std = Xm[,17],
                Dist_WFm.std = Xm[,18])

#######################################################################################################
# And here is the stack.
# Opcion 1
StackFit <- inla.stack(
  tag  = "Fit",
  data = list(y = Datos$DEP),  
  A    = list(1,A2, 1, 1),  # Intercept and covariates, spatial field                    
  effects = list(X = X, #Intercept and Covariates
                 w = w.index,          #Spatial field  
                 Nest = Datos$Nest,
                 Station = Datos$Station
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

# WIND FARMS, TANDA AND PREDATORS
f1 <- y ~ -1 + Intercept + WF + Tanda + Avian_Predators.std + Excrementos_potentials.std + f(Station, model = "iid") + f(Nest, model = "iid")
f2 <- y ~ -1 + Intercept + WF + Tanda + Avian_Predators.std + Excrementos_potentials.std + f(w, model = spde)

# DISTANCES + LAND USES
f1 <- y ~ -1 + Intercept + DIST_TREE.std + DIST_cropm.std + DIST_paths.std + GLOBSuelo_no_edificado.std + GLOBCrop_surface.std + GLOBTree_surface.std + f(Station, model = "iid") + f(Nest, model = "iid")
f2 <- y ~ -1 + Intercept + DIST_TREE.std + DIST_cropm.std + DIST_paths.std + GLOBSuelo_no_edificado.std + GLOBCrop_surface.std + GLOBTree_surface.std + f(w, model = spde)

# MICROHABITAT
f1 <- y ~ -1 + Intercept + PC1hor.std + PC2hor.std + PC3hor.std + PC4hor.std + PC1ver.std + PC2ver.std + f(Station, model = "iid") + f(Nest, model = "iid")
f2 <- y ~ -1 + Intercept + PC1hor.std + PC2hor.std + PC3hor.std + PC4hor.std + PC1ver.std + PC2ver.std + f(w, model = spde)

# FINAL MODEL
f1 <- y ~ -1 + Intercept + Tanda + PC3hor.std + GLOBSuelo_no_edificado.std + GLOBCrop_surface.std + f(Station, model = "iid") + f(Nest, model = "iid")
f2 <- y ~ -1 + Intercept + Tanda + PC3hor.std + GLOBSuelo_no_edificado.std + GLOBCrop_surface.std + f(w, model = spde)

# Repeated measures design
f3 <- y ~ -1 + Intercept + Tanda + PC3hor.std + GLOBSuelo_no_edificado.std + GLOBCrop_surface.std + GLOBTree_surface.std + f(Nest, model = "iid") + f(w)
# This model has higher wAIC values because it is less parsimonous, it incorporates de spatial component and the random intercept for the Nest

#############################################
# 7. Run the spatial model in INLA.
# First we run the model without spatial dependency.

I1.bernoulli <- inla(f1,
                   family = "binomial", 
                   data = inla.stack.data(StackFit),
                   control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
                   control.predictor = list(A = inla.stack.A(StackFit)))

I2.bernoulli <- inla(f2,
                     family = "binomial", 
                     data = inla.stack.data(StackFit),
                     control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
                     control.predictor = list(A = inla.stack.A(StackFit)))

# And compare the models with DICs and WAICs
dic  <- c(I1.bernoulli$dic$dic, I2.bernoulli$dic$dic)
waic <- c(I1.bernoulli$waic$waic, I2.bernoulli$waic$waic)
Z2     <- cbind(dic, waic)
rownames(Z2) <- c("Mixed Effects (Station/Nest)", "SPDE")
Z2

# WITH SPDE PC PRIORS P(r < 1) = 0.5 and P(sigma > 1) = 0.05
#                                  dic     waic
# Mixed Effects (Station/Nest)     Inf 234.2109
# SPDE                         236.387 234.3780
I1.bernoulli$summary.fixed
I2.bernoulli$summary.fixed

# We get Inf values for DIC. As it is dicussed in the internet, this could be because
# the prior distributions has to be tight up. 

# https://groups.google.com/forum/embed/?parenturl=http%3A%2F%2Fwww.r-inla.org%2Fcomments-1&service=jotspot&ul=1&theme=default&place=forum%2Fr-inla-discussion-group&showpopout=true&showsearch=true#!searchin/r-inla-discussion-group/DIC/r-inla-discussion-group/JI6hFEqxlAk/2Newnet9BQAJ
# https://groups.google.com/forum/embed/?parenturl=http%3A%2F%2Fwww.r-inla.org%2Fcomments-1&service=jotspot&ul=1&theme=default&place=forum%2Fr-inla-discussion-group&showpopout=true&showsearch=true#!searchin/r-inla-discussion-group/DIC/r-inla-discussion-group/2eGHHN6xPIU/uQbzKYHHAwAJ

# We can see that the Inf values come from the sampling stations in AV2
I2.bernoulli$dic$local.dic
which(I2.bernoulli$dic$local.dic == Inf)
Infinites <- which(I2.bernoulli$dic$local.dic == Inf)
Datos$Station[Infinites]

which(I1.bernoulli$dic$local.dic == Inf) # The same
# Thus, something is going on that the random effect do not work in those sampling stations

# Change the priors
inla.models()$latent$iid$hyper # Default priors
I2.bernoulli$all.hyper$random[[1]]$hyperid
I2.bernoulli$all.hyper$random[[1]]$hyper

# For the Mixed Effect Model:
prec.prior <- list(prec = list(prior = "loggamma", param = c(2, 0.1)))
# We see that once DIC is different to Inf the change on the mean of the loggama distribution does not affect the estimate of the DIC
f2 <- y ~ -1 + Intercept + WF + Tanda + Avian_Predators_residuals.std + Excrementos_potentials.std + f(Station, model = "iid", hyper= prec.prior) + f(Nest, model = "iid", hyper= prec.prior)

# For the SPDE model
I3.bernoulli$all.hyper$random[[1]]$hyperid
I3.bernoulli$all.hyper$random[[1]]$hyper$theta1
I3.bernoulli$all.hyper$random[[1]]$group.hyper$theta

# Below an informative prior: Most likely the range is larger than 100 km.
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(100, 0.05), 
                            prior.sigma = c(1, 0.05))

# Below an uninformative prior: most likely the range is larger than 1 km
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(1, 0.05), 
                            prior.sigma = c(1, 0.05))
# Another try
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(0.1, 0.05), 
                            prior.sigma = c(1, 0.05))
# Another try
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(0.5, 0.05), 
                            prior.sigma = c(1, 0.05))
# Go to #*# and fit again the SPDE

#######################################
#8. Inspect the results of I2

# Let's plot the results of the model, without and with
# the spatial correlation side by side.
# Better not try to understand all this R code.
# RUN FROM HERE......

Combined <- rbind(I1.bernoulli$summary.fixed[, c("mean", "0.025quant", "0.975quant")],
                  I2.bernoulli$summary.fixed[, c("mean", "0.025quant", "0.975quant")])

Combined$WhichModel <- rep(c("Non-spatial", "Spatial"), 
                           each = nrow(I1.bernoulli$summary.fixed))

rownames(I1.bernoulli$summary.fixed)

Combined$WhichVariable <- rep(c("Intercept", "Wind farms", "Period", "Avian predators", "Terrestrial predators"), 2)
Combined$WhichVariable <- rep(c("Intercept", "PC1Hor", "PC2Hor", "PC3Hor", "PC4Hor", "PC1Ver", "PC2Ver"), 2)
Combined$WhichVariable <- rep(c("Intercept", "D.Tree", "D.Crop", "D.Path", "S.Gravel roads", "S.Crops", "S.Trees"), 2)
Combined$WhichVariable <- rep(c("Intercept", "Period", "PC3Hor", "S.Gravel roads", "S.Crops"), 2)

colnames(Combined) <- c("Mean", "Lo", "Up", "Model", "Predictor")
positions <- unique(Combined$Predictor)
positions <- c("Intercept", "Period", "S.Gravel roads", "S.Crops", "PC3Hor")

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
  
#################################################################################################
###### MODEL SELECTION

# Let's focus on the hyper-parameters.
# The code is a 'little' bit on the ugly side.
# Just run it without asking why.
I2 <- I2.bernoulli # JUST SPDE

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

I2.bernoulli$summary.hyper 

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

# Plot fitted values versus observed data
mu2 <- I2$summary.fitted.values[1:nrow(Datos), "mean"]
plot(x = mu2,
     y = Datos$DEP)
# Hmm. Looks too good?

# Check the Pearson residuals of I2 for
# any remaining spatial dependency.
# Variogram model I2
inla.stack.index(StackFit, tag = "Fit")$data

############################################################################################################
##### MODEL VALIDATION
#Variogram Model I1 (without Spatial Autocorrelationa)
# PEARSON RESIDUALS
N <- nrow (Datos)

I2 <- I1.bernoulli # Mixed Effects Model
Pi <- I2$summary.fitted.values[1:N, "mean"]
ExpY <- Pi
VarY <- (1 - Pi) * Pi
E2 <- (Datos$DEP - ExpY) / sqrt(VarY)

I3 <- I2.bernoulli # SPDE
N <- nrow (Datos)
Pi <- I3$summary.fitted.values[1:N, "mean"]
ExpY <- Pi
VarY <- (1 - Pi) * Pi
E3 <- (Datos$DEP - ExpY) / sqrt(VarY)

MyData2 <- data.frame(E1 = E1,
                      E2 = E2,
                      E3 = E3,
                      X  = Datos$X/1000, 
                      Y  = Datos$Y/1000)
coordinates(MyData2) <- c("X", "Y")
V1 <- variogram(E1 ~ 1, 
                MyData2, 
                cressie = TRUE)

V2 <- variogram(E2 ~ 1, 
                MyData2, 
                cressie = TRUE)

V3 <- variogram(E3 ~ 1, 
                MyData2, 
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

p <- p + geom_point(data = V2,
                    aes(x = dist, 
                        y = gamma)
)
p <- p + geom_line(data = V2,
                   aes(x = dist, 
                       y = gamma),
                   col = "blue"
)
p <- p + geom_point(data = V3,
                    aes(x = dist, 
                        y = gamma)
)
p <- p + geom_line(data = V3,
                   aes(x = dist, 
                       y = gamma),
                   col = "black"
)
p <- p + xlab("Distance") + ylab("Semi-variogram")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p <- p + ylim(0,1)
p
# No clear patterns for the Poisson GLM with spatial correlation.

### MORAN'S I
## MORAN'S I TEST
Loc <- cbind(Datos$X/1000, Datos$Y/1000)
D <- dist(Loc)

w <- 1/as.matrix(dist(Loc))
diag(w) <- 0

library (spdep)
moran.test (E3, mat2listw(w))

library (ncf)
leadI <- spline.correlog(x=Loc[,1], y=Loc[,2],
                         z=E3, resamp=100, quiet=TRUE)

par (mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
par (mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(density(E3))
plot (leadI)

### PLOT RESIDUALS!
# Option 1: Use different point sizes and symbols
#           based on the values of the residuals.
E2 <- E3

MyCex <- 3 * abs(E2) / max(E2) + 0.5
Sign  <- as.numeric(E2 >=0) + 1
MyPch <- c(1, 16)[Sign]

xyplot (Longitude ~ Latitude,
        data = Datos, aspect = "iso",
        cex = MyCex, pch = MyPch, col = 1, 
        xlab = list(label = "Easting", cex = 1.5),
        ylab = list(label = "Northing", cex = 1.5))

# Patterns in the residuals are present! 
# Rain comes from the west. 
# Black dots is for positive residuals (left side of the island) whereas white dots is for negative residuals (right side of the island)

# If you see patterns then we have to stop! Standard errors, t values and estimates are not true -> They are biased
# Standard errors are faluty, and then also de t-value and also the p-values (Domino!)
# If the p-value is really low (0.0000) the effects are going to be low. But if the p-values are near 0.01 then we have troubles!


### RESIDUALS VERSUS EACH COVARIATE
# Plot residuls vs each covariate
Datos$E2 <- E3 # SPDE Only
MyVar <- c("fWF", "Tanda", "GLOBSuelo_no_edificado", "GLOBCrop_surface", "GLOBTree_surface")
MyVar <- c("fWF", "Tanda", "PC1hor.std", "PC2hor.std", "PC3hor.std","PC1ver.std","PC2ver.std","GLOBSuelo_no_edificado", "GLOBCrop_surface", "GLOBTree_surface")
MyMultipanel.ggp2(Z=Datos,
                  varx = MyVar,
                  vary = "E2",
                  ylab = "Pearson Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

####################################################
# We finally present the spatial component, the wks. 
# Their posterior mean values can be obtained via
w.pm.NB   <- I2.bernoulli$summary.random$w$mean  

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

#############################
# Simulating from the model

# In the Turco parasite example we implemented a
# 7-step protocol for simulating from the model.
# See also Powerpoint presentation: 
#    P15_SpatTempBook_Chapter10_V1.pdf


# We will now implement the 7-step protocol
# for assessing whether the Poisson GLM
# is over- or underdispersed.

# Step 2: Simulate regression parameters
# We use the function inla.posterior.sample 
# to simulate from the model. The output is 
# stored in the Sim object.
I2.sim <- I2.bernoulli

NSim <- 1000
Sim <- inla.posterior.sample(n = NSim, result = I2.sim)

# This gives an object with 2140 rows for this 
# specific data set and model. The first 155 
# rows are simulated values for eta = X * beta, 
# where X is the matrix with covariates, and 
# the last 8 rows are simulated regression parameters. 
# The middle part contains the simulated spatial field w.
# This is just one set of simulated values 
# (due to n = 1 in the function above). 


# Here is some fancy code to grap the 1  betas
MyParams <- rownames(I2.sim$summary.fixed)
RowNum.Betas <- lapply (MyParams, 
                        function (x)
                          grep (x,
                                rownames (Sim[[1]]$latent),
                                fixed = TRUE))
RowNum.Betas <- as.numeric (RowNum.Betas)
RowNum.Betas # In this rows are our parameters!
Sim[[1]]$latent[RowNum.Betas,]

# Here is some fancy code to grap the w 
Nw <- mesh$n 
MyParams <- paste("w", seq(from = 1, to = Nw), sep = ":")
MyID <- function(x){ which(rownames(Sim[[1]]$latent) == x) }
RowNum.w <- lapply(MyParams, MyID)
RowNum.w <- as.numeric(RowNum.w)
RowNum.w

# Step 3: Calculate predicted values
# It is identical to the Turco parasite example, except
# there is now also the A * w component.

Betas <- Sim[[1]]$latent[as.numeric (RowNum.Betas)]
wk    <- Sim[[1]]$latent[RowNum.w]

X <- model.matrix (~ fWF + fTanda + Avian_Predators.std + Excrementos_potentials.std
                # + f(Station, model = "iid")
                # + f(Nest, model = "iid")
                , data = Datos)

X <- model.matrix (~ PC1hor.std + PC2hor.std + PC3hor.std + PC4hor.std + 
                     PC1ver.std + PC2ver.std
                   # + f(Station, model = "iid")
                   # + f(Nest, model = "iid")
                   , data = Datos)

X <- model.matrix (~ Dist_WFm.std + DIST_TREE.std + DIST_cropm.std + 
                     DIST_paths.std + GLOBSuelo_no_edificado.std + GLOBCrop_surface.std + 
                     GLOBTree_surface.std
                   # + f(Station, model = "iid")
                   # + f(Nest, model = "iid")
                   , data = Datos)

X <- model.matrix (~ fTanda + PC3hor.std + GLOBSuelo_no_edificado.std + GLOBCrop_surface.std + 
                     GLOBTree_surface.std
                   # + f(Station, model = "iid")
                   # + f(Nest, model = "iid")
                   , data = Datos)

# This is just one simulated data set.
# Repeat the whole thing 1000 times.

Xm <- as.matrix(X)
Am <- as.matrix(A2)

N  <- nrow(Datos)
Ysim <- matrix(nrow = N, ncol = NSim)
mu.i <- matrix(nrow = N, ncol = NSim)

for (i in 1:NSim){
  Betas <- Sim[[i]]$latent[RowNum.Betas]
  wk    <- Sim[[i]]$latent[RowNum.w]
  FixedPart   <- Xm %*% Betas
  SpatialPart <- Am %*% wk
  mu.i[,i]    <- exp(FixedPart + SpatialPart) / (1 + exp(FixedPart + SpatialPart))
  Ysim[,i]    <- rbinom (n = nrow(Datos), size = 1, prob = mu.i[,i])
}

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
     xlab = "Simulated probability of Nest Predation",
     ylab = "Frequencies",
     ylim = c(0, 260))

# Add SDs
Zs <- table(Datos$DEP)
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
  segments(x0 = NamesZ[i] + 0.01,
           x1 = NamesZ[i] + 0.01,
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
# Now that we have 1000 simulated data sets from the 
# model, what else shall we do with these simulated 
# data sets?

# BERNOULLI MODELS CAN NOT BE OVERDISPERSED OR ZERO-INFLETED!!

##################################################################
## RESULTS!

I3.bernoulli$summary.fixed

