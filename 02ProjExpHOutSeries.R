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

################################################################################
# 02 PROJECTED EXPOSURE AND HEALTH OUTCOME SERIES
################################################################################

################################################################################
# PROJECTED TEMPERATURE SERIES:

# LOAD MODELLED TEMPERATURE DATA (1971-2099)
# Climate data was obtained, processed and made available by 
#    the Inter-Sectoral Impact Model Intercomparison Project 
#    (ISI-MIP, https://www.isimip.org/) .
# It includes daily temperature series for 5 different GCMs and
#    2 (RCP4.5 and RCP8.5) climate change scenarios, 
#    for the historical (1971-2005) and future periods (2006-2100).

# - RCP4.5 SCENARIO 
rcp4p5 <- read.csv("lndn_rcp4p5.csv",colClasses=c(date="Date"))

# - RCP8.5 SCENARIO 
rcp8p5 <- read.csv("lndn_rcp8p5.csv",colClasses=c(date="Date"))


# PLOT - FIGURE 2 TEMPERATURE TRENDS
#  The plot shows the temporal trends annual temperature averaged across 
#     the 5 GCMs for each scenario (solid line), expressed as increase in 
#     annual mean temperature from the historical average temperature. 
#  It also shows the variability in the projected temperature across GCMs 
#     (shade area).


# - COMPUTE ANNUAL AVERAGES (CENTERED)
tsyear <- as.numeric(format(rcp4p5$date,format="%Y"))
year <- unique(tsyear)
rcp4p5avg <- apply(rcp4p5[,-1],2,tapply,tsyear,mean)
rcp8p5avg <- apply(rcp8p5[,-1],2,tapply,tsyear,mean)

# - COMPUTE THE DIFFERENCE IN TEMPERATURE FROM THE HISTORICAL AVERAGE
rcp4p5avg <- rcp4p5avg-mean(rcp4p5avg[year%in%1986:2005,])
rcp8p5avg <- rcp8p5avg-mean(rcp8p5avg[year%in%1986:2005,])

ylab <- expression(paste("Temperature increase (",degree,"C)"))

col1 <- do.call(rgb,c(as.list(col2rgb("blue")),alpha=255/6,max=255))
col2 <- do.call(rgb,c(as.list(col2rgb("red")),alpha=255/6,max=255))

#pdf("Figure 2.pdf",height=5,width=5)
par(mar=c(5,4,3,1),mgp=c(2,0.8,0),las=1,cex.axis=0.8,cex.lab=1)
plot(year,rowMeans(rcp8p5avg),type="n",bty="l",ylab=ylab,xlab="Year",
  ylim=c(-1,6),xlim=c(1970,2110),main="")
abline(h=0)
lines(year[year<2006],rowMeans(cbind(rcp4p5avg,rcp8p5avg)[year<2006,]),lwd=1.5)
polygon(c(year[year>=2006],rev(year[year>=2006])),
  c(apply(rcp4p5avg,1,max)[year>=2006],rev(apply(rcp4p5avg,1,min)[year>=2006])),
  col=col1,border=NA)
polygon(c(year[year>=2006],rev(year[year>=2006])),
  c(apply(rcp8p5avg,1,max)[year>=2006],rev(apply(rcp8p5avg,1,min)[year>=2006])),
  col=col2,border=NA)
lines(year[year>=2006],rowMeans(rcp4p5avg[year>=2006,]),lwd=1.5,col=4)
lines(year[year>=2006],rowMeans(rcp8p5avg[year>=2006,]),lwd=1.5,col=2)
abline(v=2005,lty=3)
rect(2103-1,mean(apply(rcp4p5avg,1,min)[as.character(2090:2099)]),2103+1,
  mean(apply(rcp4p5avg,1,max)[as.character(2090:2099)]),col=4,border=NA)
rect(2106-1,mean(apply(rcp8p5avg,1,min)[as.character(2090:2099)]),2106+1,
  mean(apply(rcp8p5avg,1,max)[as.character(2090:2099)]),col=2,border=NA)
legend("top",c("Historical","RCP4.5","RCP8.5"),lwd=1.5,col=c(1,4,2),cex=0.8,bty="n",inset=0.05)
#dev.off()


################################################################################
# PROJECTED MORTALITY SERIES:

# DEFINE PROJECTED MORTALITY SERIES
# It is computed as the average mortality for each day of the year 
#   from daily observed deaths, then repeated along the same projection period
#   of the modelled temperature series.

deathdoy <- tapply(obs$all,as.numeric(format(obs$date,"%j")),
  mean,na.rm=T)[seq(365)]
while(any(isna <- is.na(deathdoy)))
  deathdoy[isna] <- rowMeans(Lag(deathdoy,c(-1,1)),na.rm=T)[isna]
deathproj <- rep(deathdoy,length=nrow(rcp8p5))


# PLOT - FIGURE 3 SEASONAL MORTALITY TRENDS
# The plot depicts the seasonal trends in daily mortality in London.
#    The grey dots are the actual number of deaths in each day of the year 
#    between 1992-2005 and the blue line is the average daily mortality 
#    per day of the year, which is projected along the projection period.

#pdf("Figure 3.pdf",height=5,width=6)
par(mar=c(5,4,3,1),mgp=c(2.5,1,0),las=1,cex.axis=0.8,cex.lab=1)
plot(all ~ doy,obs,pch=19,col=grey(0.8),cex=0.5,ylab="Daily deaths",xlab="Day of the year",
  main="")
lines(1:366,tapply(obs$all,obs$doy,mean),lwd=1.5,col=4)
#dev.off()
