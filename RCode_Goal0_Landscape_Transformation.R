
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

##### We load the same data than for the Abundance of Terrestrial predators, because there
# we have the information for each sampling station
Datos <- read_excel ("Excrementos2016_20190205.xlsx")
names(Datos)
str (Datos)

Datos$fWF <- factor (Datos$WF, levels = c(0,1), labels = c("No", "Yes"))

## Standardized variables
Datos$Suelo_no_edificado_estd <- scale (Datos$Suelo_no_edificado)
Datos$Crop_surface_estd <- scale (Datos$Crop_surface)
Datos$Tree_surface_estd <- scale (Datos$Tree_surface)

#########################################################################################################
##### INCORPORATE SPATIAL AUTOCORRELATION
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

# And here is the stack.
StackFit_suelo <- inla.stack(
  tag  = "Fit",
  data = list(y = Datos$Suelo_no_edificado_estd),  
  A    = list(1,A2),  # Intercept and covariates, spatial field                    
  effects = list( 
    X = X, #Covariates
    w = w.index          #Spatial field  
  )) 

StackFit_crops <- inla.stack(
  tag  = "Fit",
  data = list(y = Datos$Crop_surface_estd),  
  A    = list(1,A2),  # Intercept and covariates, spatial field                    
  effects = list( 
    X = X, #Covariates
    w = w.index          #Spatial field  
  )) 

StackFit_trees <- inla.stack(
  tag  = "Fit",
  data = list(y = Datos$Tree_surface_estd),  
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

f2 <- y ~ -1 + Intercept + fWF + 
  f(w, model = spde)

#############################################
# 7. Run the spatial model in INLA.
# First we run the model without spatial dependency.

#### PARA EL SUELO
I1_suelo <- inla(f1,
           family = "gaussian", 
           data = inla.stack.data(StackFit_suelo),
           control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
           control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit_suelo)))

I2_suelo <- inla(f2,
          family = "gaussian", 
          data = inla.stack.data(StackFit_suelo),
          control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
          control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit_suelo)))


# And compare the models with DICs and WAICs
dic  <- c(I1_suelo$dic$dic, I2_suelo$dic$dic)
waic <- c(I1_suelo$waic$waic, I2_suelo$waic$waic)
Zsuelo     <- cbind(dic, waic)
rownames(Zsuelo) <- c("Gaussian GLM",
                 "Gaussian GLM + SPDE")
Zsuelo

I2_suelo$summary.fixed[, c("mean", "0.025quant", "0.975quant")]

#### PARA CROPS
I1_Crops <- inla(f1,
                 family = "gaussian", 
                 data = inla.stack.data(StackFit_crops),
                 control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
                 control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit_crops)))

I2_Crops <- inla(f2,
                 family = "gaussian", 
                 data = inla.stack.data(StackFit_crops),
                 control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
                 control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit_crops)))


# And compare the models with DICs and WAICs
dic  <- c(I1_Crops$dic$dic, I2_Crops$dic$dic)
waic <- c(I1_Crops$waic$waic, I2_Crops$waic$waic)
Zcrops     <- cbind(dic, waic)
rownames(Zcrops) <- c("Gaussian GLM",
                 "Gaussian GLM + SPDE")
Zcrops

I2_Crops$summary.fixed[, c("mean", "0.025quant", "0.975quant")]

##### PARA TREES

I1_Trees <- inla(f1,
                 family = "gaussian", 
                 data = inla.stack.data(StackFit_trees),
                 control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
                 control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit_trees)))

I2_Trees <- inla(f2,
                 family = "gaussian", 
                 data = inla.stack.data(StackFit_trees),
                 control.compute = list(config= TRUE, dic = TRUE, waic = TRUE),
                 control.predictor = list(compute= TRUE, A = inla.stack.A(StackFit_trees)))


# And compare the models with DICs and WAICs
dic  <- c(I1_Trees$dic$dic, I2_Trees$dic$dic)
waic <- c(I1_Trees$waic$waic, I2_Trees$waic$waic)
Ztrees    <- cbind(dic, waic)
rownames(Ztrees) <- c("Gaussian GLM",
                      "Gaussian GLM + SPDE")
Ztrees

I2_Trees$summary.fixed[, c("mean", "0.025quant", "0.975quant")]

#######################################

# Let's plot the results of the model, without and with
# the spatial correlation side by side.
# Better not try to understand all this R code.
# RUN FROM HERE......
I1 <- I1_suelo
I1 <- I1_Crops
I1 <- I1_Trees

I2 <- I2_suelo
I2 <- I2_Crops
I2 <- I2_Trees

Combined <- rbind(I1$summary.fixed[, c("mean", "0.025quant", "0.975quant")],
                  I2$summary.fixed[, c("mean", "0.025quant", "0.975quant")])

Combined$WhichModel <- rep(c("Non-spatial", "Spatial"), 
                           each = nrow(I2$summary.fixed))

rownames(I2$summary.fixed)

Combined$WhichVariable <- rep(c("Intercept", "Wind farms"), 2)

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


############################################################################
## HYPERPARAMETERS
I2 <- I2_suelo

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