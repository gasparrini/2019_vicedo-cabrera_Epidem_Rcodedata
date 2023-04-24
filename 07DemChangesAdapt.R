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

################################################################################
# 07 ACCOUNTING FOR DEMOGRAPHIC CHANGES AND ADAPTATION
################################################################################

################################################################################
# - 07.1 DEMOGRAPHIC CHANGES
################################################################################

# The method requires an extension of the method described in the previous steps,
#   basically performing age-specific projections using estimates of 
#   exposure-response relationships and mortality series specific for each age
#   group.

# To do so, we first (1) replace Rscript 01 with the following code to run
#   age-specific models and obtain the corresponding exposure-response curves

################################################################################
# DEFINITION OF THE CROSS-BASIS FOR TEMPERATURE
# NB: the same specifications of the main cross-basis defined in sript 01

argvar <- list(fun="ns", knots = quantile(obs$tmean,c(10,75,90)/100, na.rm=T),
  Bound=range(obs$tmean,na.rm=T))
maxlag <- 21
arglag <- list(fun="ns",knots=logknots(maxlag,nk=3))

cb <- crossbasis(obs$tmean,maxlag,argvar,arglag)  


# FIT THE MODEL
# Apply the same model defined in the script 01 to the age-specific counts:

m0_64 <- glm(all_0_64 ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
  data=obs, family=quasipoisson)
m65_74 <- glm(all_65_74 ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
  data=obs, family=quasipoisson)
m75_84 <- glm(all_75_84 ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
  data=obs, family=quasipoisson)
m85plus <- glm(all_85plus ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
  data=obs, family=quasipoisson)

# GET PREDICTIONS & PLOT
# - DEFINE PROVVISIONAL CENTERING POINT TO HAVE THE INITIAL PREDICTION
varcen <- 19 

# - ESTIMATE MMT FROM THE PREDICTED EXPOSURE-RESPONSE ASSOCIATION
#     FOR EACH AGE CATEGORY

cp0_64 <- crosspred(cb,m0_64,cen=varcen,by=0.1)
cen0_64 <- cp0_64$predvar[which.min(cp0_64$allRRfit)] 

cp65_74 <- crosspred(cb,m65_74,cen=varcen,by=0.1)
cen65_74 <- cp65_74$predvar[which.min(cp65_74$allRRfit)] 

cp75_84 <- crosspred(cb,m75_84,cen=varcen,by=0.1)
cen75_84 <- cp75_84$predvar[which.min(cp75_84$allRRfit)] 

cp85plus <- crosspred(cb,m85plus,cen=varcen,by=0.1)
cen85plus <- cp85plus$predvar[which.min(cp85plus$allRRfit)] 

# - RE-CENTER & GET PREDICTIONS FOR EACH MODEL CENTERING ON THE MMT 
pred0_64 <- crosspred(cb, m0_64, cen=cen0_64, by=1)   
pred65_74 <- crosspred(cb, m65_74, cen=cen65_74, by=1)   
pred75_84 <- crosspred(cb, m75_84, cen=cen75_84, by=1)  
pred85plus <- crosspred(cb, m85plus, cen=cen85plus, by=1)  

################################################################################

# The next steps consist in (2) estimating the projected mortality for each age
#   group. This step, not shown here, can rely on age-specific baseline population
#   and mortality rates built under specific scenarios. The same projected
#   temperature series can be used for all the age groups.
# Steps (3)-(4) (recalibration and extrapolation) remains the same. We can then
#   project the corresponding (5) mortality and the corresponding (6) uncertainty
#   for each age group. Counts in each category are then aggregated to estimate
#   the total projected mortality. 


################################################################################
# - 07.2 ADAPTATION
################################################################################

# To account for adaptation, an estimate of the modified exposure-response curve
#   that represents the attenuation in risk is needed.
# As a simple example, we show how to implement the assumption of a decrease of
#   30% in the risk associated to heat, by applying a scaling factor to log-risk
#   in the related part of the curve. 
# Below, we show how to plot the modified exposure-response curve corresponding 
#   to teh assumed reduction. To implement this step in the framework with the
#   corresponding reduction in excess mortality, the line 128 of the code
#   "04_05_06ExtCurveProjUncert.R" can be modified accordingly.


# PLOT - FIGURE 6

#pdf("Figure 6.pdf",height=5,width=9)

# FIGURE 6A DEMOGRAPHIC CHANGES
layout(matrix(c(1,2),ncol=2,byrow=T))
xlab <- expression(paste("Temperature (",degree,"C)"))

par(mar=c(5,4,4,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred0_64,"overall",col="magenta",ylim=c(0.5,2.5),axes=T,lab=c(6,5,7),xlab=xlab,
  ylab="RR",ci="n",main="Demographic changes", lwd=2)
lines(pred65_74, col="blue", lwd=2, lty=2)
lines(pred75_84, col="forestgreen", lwd=2, lty=3)
lines(pred85plus, col="orange", lwd=2, lty=4)

legend(8,2.4,c("<65 years","65-74 years","75-84 years", ">84 years")
  ,xpd = TRUE,col=c("magenta","blue","forestgreen", "orange"),lwd=2,bg="white"
  ,cex=0.9,ncol=1,inset=0.04,bty="n", lty=1:4)
text(12,2.5, "Age-specific relationships", cex=1)


# FIGURE 6B ADAPTATION

# INDICATOR FOR HEAT (ABOVE MMT)
ind <- pred$predvar>=cen 
par(mar=c(5,4,4,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)

plot(pred,"overall",col="red",ylim=c(0.5,2.5),axes=T,lab=c(6,5,7),xlab=xlab,
  ylab="RR",ci="n",main="Adaptation")

lines(pred$predvar,pred$allRRfit,col=2,lwd=2)
# PLOT LINE OF THE MODIFIED CURVE (30% REDUCTION OF HEAT) 
lines(pred$predvar[ind],exp(pred$allfit[ind]*0.7),col=2,lwd=2,lty=2)
legend(8,2.4,c("No adaptation","With adaptation"),lty=1:2,inset=0.005,
  cex=0.9,bty="n",col=2)
text(12,2.5, "Reduction in heat risk", cex=1)

layout(1)
#dev.off()

# NB: run pdf() and dev.off() for saving the plot in pdf format.





