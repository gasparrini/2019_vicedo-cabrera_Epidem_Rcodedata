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
# 03 BIAS-CORRECTION
################################################################################

# The modelled temperature included in the files "lndn_rcp4p5.csv" and 
#    "lndn_rcp8p5.csv" have been previously downscaled through bi-linear 
#    interpolation at a 0.5?X0.5?spatial resolution and linear interpolated 
#    by day of the year. 


# LOAD BIAS CORRECTION FUNCTION
# This is a function created to apply the bias-correction method developed 
#    within ISI-MIP (Hempel et al. 2013). More details on the calibration 
#    procedure are described in the "fhempel.R" code.

source("fhempel.R")

# RE-CALIBRATE USING THE BIAS CORRECTION FUNCTION
rcp4p5cal <- fhempel(obs[c("date","tmean")],rcp4p5)
rcp8p5cal <- fhempel(obs[c("date","tmean")],rcp8p5)


# PLOT - FIGURE 4 CALIBRATION
# The two plots compared the distribution between the observed and the modelled 
#    temperature obtained using one specific GCM, before and after applying the 
#    calibration. 

#pdf("Figure 4.pdf",height=5,width=9)
layout(matrix(c(1,2),ncol=2,byrow=T))
par(mar=c(5,4,3,1),mgp=c(3,1,0),las=1,cex.axis=0.8,cex.lab=1)
rgbcol <- col2rgb(c(grey(0.8),"lightgreen"))
col <- apply(rgbcol,2,
  function(x) do.call(rgb,c(as.list(x),alpha=255/2,max=255)))
ind <- rcp4p5$date %in% obs$date
dobs <- density(obs$tmean,na.rm=T)
dmod <- density(rcp4p5[ind,3],na.rm=T)
plot(c(obs$tmean,rcp4p5[ind,3]),c(obs$tmean,rcp4p5[ind,3]),type="n",bty="l",
  main="",ylim=range(c(dobs$y,dmod$y))*1.05,ylab="Distribution",
  xlab=expression(paste("Temperature (",degree,"C)")))
polygon(dobs,col=col[1],lty=4)
polygon(dmod,col=col[2],lty=5)
legend("topleft",c("Observed","GCM"),pch=22,pt.bg=col, cex=0.9,bty="n",
  inset=0.05)
abline(h=0)


par(mar=c(5,4,3,1),mgp=c(3,1,0),las=1,cex.axis=0.8,cex.lab=1)
ind <- rcp4p5$date %in% obs$date
x <- pretty(c(obs$tmean,rcp4p5[,3]),n=100)
plot(x,ecdf(obs$tmean)(x),col=1,main="",type="l",lwd=2,
  ylab="Cumulative distribution",
  xlab=expression(paste("Temperature (",degree,"C)")))
lines(x,ecdf(rcp4p5[ind,3])(x),col=3,lwd=2)
lines(x,ecdf(rcp4p5cal[ind,3])(x),col=3,lwd=2,lty=4)
legend("topleft",c("Observed","GCM","GCM calibrated"),lwd=2,lty=c(1,1,4),
  col=c(1,3,3), cex=0.9,bty="n",inset=0.05)

layout(1)
#dev.off()

# NB: run pdf() and dev.off() for saving the plot in pdf format.


# REPLACE THE ORIGINAL SERIES WITH THE RECALIBRATED
rcp4p5 <- rcp4p5cal
rcp8p5 <- rcp8p5cal
