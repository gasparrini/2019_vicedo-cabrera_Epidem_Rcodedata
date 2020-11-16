################################################################################
# Updated version of the R code for an example of the analysis in:
#
#  "A hands-on tutorial on a modelling framework for projections of climate 
#    change impacts on health"
#   Ana M. Vicedo-Cabrera, Francesco Sera, Antonio Gasparrini
#  http://www.ag-myresearch.com/2019_vicedo-cabrera_Epidem.html
#
# This code reproduces the analysis described in the tutorial.
# R Code and data released under GNU General Public License (version 3).
#
# Update: 29 May 2019
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata
#
# Requirements:
#  - R version 3.5.1
#  - 'dlnm' package version 2.3.5
#  - 'splines' package version  3.5.1
#  - 'MASS' package version 7.3-50
#
################################################################################

# LOAD THE PACKAGES
library(dlnm) ; library(splines) ; library(MASS) ; library(tsModel)

################################################################################
# 01 ESTIMATION OF THE EXPOSURE-RESPONSE ASSOCIATIONS
################################################################################

# Use the observed daily temperature-mortality series to estimate
#     the coefficients defining the exposure-response association.
# In this example, we estimate the temperature-related mortality association 
#     using data from London (UK) between 1990 and 2012.


# LOAD OBSERVED DATA - DAILY TEMPERATURE-SERIES BETWEEN 1990-2012 IN LONDON
#  NB: file "lndn_obs.csv" provided as downloadable supplemental material  
#      along this code. Description of the variables available in 
#     "VarDescr_lndn_obs.csv".
obs <- read.csv("lndn_obs.csv")
obs$date <- as.Date(obs$date, format="%d/%m/%Y")

# DEFINITION OF THE CROSS-BASIS FOR TEMPERATURE
# - SPECIFICATION PARAMETERS OF THE EXPOSURE-RESPONSE DIMENSION OF THE CROSS-BASIS

# argvar: main model, cubic natural spline with three internal knots in 
#   the 10th, 75th, 90th percentiles of the temperature distribution

argvar <- list(fun="ns", knots = quantile(obs$tmean,c(10,75,90)/100, na.rm=T),
  Bound=range(obs$tmean,na.rm=T))

# argavar1 & argvar2: additional tests, argavar1 fits a linear function, and
#   argvar2 fits a double threshold, with linear relationship for temperatures 
#   below and above the 10th and 90th percentiles,
#   with null association in between.

argvar1 <- list(fun="lin")
argvar2 <- list(fun="thr", thr.value=quantile(obs$tmean,c(10,90)/100, na.rm=T),
  side="d")

# - SPECIFICATION PARAMETERS OF THE LAG-ASSOCIATION DIMENSION OF THE CROSS-BASIS
# Definition of the maximum lag, that is, 21 days
maxlag <- 21
# arglag: main model, it fits a cubic natural spline with three internal knots 
#   equally-spaced in the log-scale.
arglag <- list(fun="ns",knots=logknots(maxlag,nk=3))

# - CREATE CROSSBASIS OBJECTS
cb <- crossbasis(obs$tmean,maxlag,argvar,arglag)  
cb1 <- crossbasis(obs$tmean,maxlag,argvar1,arglag)   
cb2 <- crossbasis(obs$tmean,maxlag,argvar2,arglag)

# FIT THE MODEL
# Include in the model the crossbasis term of temperature, along with the 
#    indicator for day of day of the week (dow) and natural cubic spline of 
#    time with 8 df per year.

m <- glm(all ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
  data=obs, family=quasipoisson)

m1 <- glm(all ~ cb1 + dow + ns(date,df=round(8*length(date)/365.25)), 
  data=obs, family=quasipoisson)
m2 <- glm(all ~ cb2 + dow + ns(date,df=round(8*length(date)/365.25)), 
  data=obs, family=quasipoisson)

# GET PREDICTIONS & PLOT
# - DEFINE PROVISIONAL CENTERING POINT TO HAVE THE INITIAL PREDICTION
#    NB: 'varcen' is initially set to a provisional temperature value within 
#        the observed range. 
varcen <- 19 

# - ESTIMATE MMT FROM THE PREDICTED EXPOSURE-RESPONSE ASSOCIATION 
# MMT corresponds to the temperature of minimum mortality, which will be used as
#    as reference to estimate relative risks and as temperature threshold 
#    to differentiate the contribution of heat and cold to the total mortality 
#    attributable to non-optimal temperatures.

cp <- crosspred(cb,m,cen=varcen,by=0.1)
cen <- cp$predvar[which.min(cp$allRRfit)] 

# - RE-CENTER & GET PREDICTIONS FOR EACH MODEL CENTERING ON THE MMT 
pred <- crosspred(cb, m, cen=cen, by=1)   

pred1 <- crosspred(cb1, m1, cen=cen, by=0.1)   
pred2 <- crosspred(cb2, m2, cen=cen, by=0.1)   


# PLOT - FIGURE 1

xlab <- expression(paste("Temperature (",degree,"C)"))

#pdf("Figure 1.pdf",height=4.5,width=10)
layout(matrix(c(1,2,3),ncol=3,byrow=T))

# PLOT - 3D
par(mar=c(2,3,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"3d",ltheta=150,xlab="Temperature (C)",ylab="Lag",zlab="RR", 
  col=gray(0.9), main="Exposure-lag-response")

# OVERALL
# The plots show the cumulative exposure-response association, in terms of 
#    relative risks (RR) and centered in the MMT, across the 21 days of lag.
par(mar=c(5,4,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"overall",col="red",ylim=c(0.5,2.5),axes=T,lab=c(6,5,7),xlab=xlab,
  ylab="RR",main="Overall")

# OVERALL DIFFERENT FUNCTIONS
# See the different shapes of the exposure-response association using the
#   three functions (non-linear, linear, double threshold).
par(mar=c(5,4,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"overall",col="red",ylim=c(0.5,2.5),axes=T,lab=c(6,5,7),xlab=xlab,
  ylab="RR",ci="n",main="Overall - different functions")
lines(pred1, col="blue")
lines(pred2, col="forestgreen")

legend("top",c("Natural spline (Main)","Linear","Double threshold")
  ,xpd = TRUE,col=c("red","blue","forestgreen"),lwd=2,bg="white"
  ,cex=0.8,ncol=1,inset=0.04)

layout(1)
#dev.off()

# NB: run pdf() and dev.off() for saving the plot in pdf format.