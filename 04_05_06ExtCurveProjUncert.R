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
# 04 EXTRAPOLATION OF THE EXPOSURE-RESPONSE CURVE
# 05 PROJECTION & QUANTIFICATION OF THE IMPACT
# 06 ENSEMBLE ESTIMATES & QUANTIFICATION OF THE UNCERTAINTY
################################################################################

# The three last steps of the analysis (extrapolation of the curve, 
#   impact projections and quantification of the uncertainty) can be performed 
#   sequentially using the following code.

# In brief, once we extrapolate the curve, we estimate the daily number of attributable 
#   deaths (AN) in each scenario, GCM and temperature range. 
# Then, we compute the sum of ANs per each 10-years period for the ensemble 
#   per each RCP and temperature range. We also estimate the difference in AN 
#   relative to the ANs estimated for the current days (2010-19). 
# By dividing between the total mortality, we estimate the corresponding 
#   attributable fractions (AFs).
# Uncertainty of the estimated impacts is expressed in terms of empirical 
#   confidence intervals, defined as the 2.5th and 97.5th percentiles of the 
#   empirical distribution of the impacts across coefficients samples and GCMs. 
#   The distribution is obtained through Monte Carlo simulation of the coefficients.


# EXTRACT COEFFICIENTS AND VCOV FROM THE MODEL IN STAGE 1
# With the crossreduce function we can reduce the fit of the bidimensional DLNM
#   (of the stage 1 using the observed temperature-mortality series)
#   to summaries defined in the exposure-response dimension. We can then 
#   extract the coefficients and the covariance matrix defining the overall 
#   exposure-response association.

red <- crossreduce(cb,m,cen=cen)
coef <- coef(red)
vcov <- vcov(red)

# STORE THE MORTALITY BY PERIOD
deathperiod <- sum(deathdoy)*10

# CREATE TEMPORARY OBJECT TO STORE THE ESTIMATED AN 
# We create an array with 6 dimensions to store the AN estimated in each 
#   simulation for each period, range of temperature, gcm and rcp. An additional
#   dimension is also included to store the absolute AN (estimated directly from
#   the formula below), and the difference in AN relative to current period.

# - DEFINE THE NAMES AND NUMBER OF LEVELS PER DIMENSION (ALSO USED TO DEFINE 
#     THE INDEX TO RUN THE LOOPS BELOW)

# (1) DIMENSION - 10-YEAR PERIOD 
#  *LABEL THE HISTORICAL PERIOD
histperiod <- "1971-1989"

#  *LABELS THE PROJECTION PERIODS
projperiod <- paste(199:209,0,"-",substr(199:209,3,3),9,sep="")

#  *DEFINE SEQUENCE OF PERIODS FOR THE PREDICTIONS (HISTORICAL & PROJECTED)
histseqperiod <- factor(rep(histperiod,length.out=365*length(seq(1971,1989))))
projseqperiod <- factor(rep(projperiod,each=365*10))
seqperiod <- factor(c(as.numeric(histseqperiod)-1,as.numeric(projseqperiod)))
levels(seqperiod) <- c(histperiod,projperiod)

# (2) DIMENSION - RANGE OF TEMPERATURES
temprange <- c("tot","cold","heat")

# (3) DIMENSION - ABSOLUTE AN/DIFFERENCE IN AN
absrel <- c("abs","rel")

# (4) DIMENSION - GENERAL CIRCULATION MODELS
# *LIST OF GCMs 
gcm <- c("GFDL-ESM2M"="tmean_gfdl","HadGEM2-ES"="tmean_hadgem2",
  "IPSL-CM5A-LR"="tmean_ipsl","MIROC-ESM-CHEM"="tmean_miroc",
  "NorESM1-M"="tmean_noresm1")

# (5) DIMENSION - SCENARIO DIMENSION
#  *LIST OF REPRESENTATIVE CONCENTRATION PATHWAYS SCENARIOS 
rcp <- c(RCP4.5="rcp4p5",RCP8.5="rcp8p5")

# (6) DIMENSION - NUMBER OF ITERATION IN THE MONTE-CARLO SIMULATION 
nsim <- 1000


# DEFINE THE ARRAY
ansim <- array(NA,dim=c(length(levels(seqperiod)),length(temprange),
  length(absrel), length(gcm),length(rcp),nsim+1), 
  dimnames=list(levels(seqperiod),temprange,absrel,
    names(gcm),names(rcp), c("est",paste0("sim",seq(nsim)))))


# RUN LOOP PER RCP
for (i in seq(rcp)) {
  
  # PRINT
  cat("\n\n",names(rcp)[i],"\n")
  
  # SELECTION OF THE PROJECTED TEMPERATURE SERIES FOR A SPECIFIC RCP SCENARIO
  tmeanproj <- get(rcp[[i]])
  
  # RUN LOOP PER GCM
  for(j in seq(gcm)) {
    
    # PRINT
    cat(gcm[j],"")
    
    # (4) EXTRAPOLATION OF THE CURVE: 
    # - DERIVE THE CENTERED BASIS USING THE PROJECTED TEMPERATURE SERIES
    #   AND EXTRACT PARAMETERS
    bvar <- do.call(onebasis,c(list(x=tmeanproj[,j+1]),argvar))
    cenvec <- do.call(onebasis,c(list(x=cen),argvar))
    bvarcen <- scale(bvar,center=cenvec,scale=F)
    
    # INDICATOR FOR COLD/HEAT DAYS
    indheat <- tmeanproj[,j+1]>cen
    
    # (5) IMPACT PROJECTIONS:
    # - COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
    an <- (1-exp(-bvarcen%*%coef))*deathproj
    
    # - SUM AN (ABS) BY TEMPERATURE RANGE AND PERIOD, STORE IN ARRAY BEFORE THE ITERATIONS
    # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
    ansim[,"tot","abs",j,i,1] <- tapply(an,seqperiod,sum)
    ansim[,"cold","abs",j,i,1] <- tapply(an[!indheat],factor(seqperiod[!indheat]),sum)
    ansim[,"heat","abs",j,i,1] <- tapply(an[indheat],factor(seqperiod[indheat]),sum)
    
    # (6) ESTIMATE UNCERTAINTY OF THE PROJECTED AN:
    # - SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
    set.seed(13041975+j)
    coefsim <- mvrnorm(nsim,coef,vcov)
    
    # - LOOP ACROSS ITERATIONS
    for(s in seq(nsim)) {
      
      # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
      an <- (1-exp(-bvarcen%*%coefsim[s,]))*deathproj
      
      # STORE THE ATTRIBUTABLE MORTALITY
      ansim[,"tot","abs",j,i,s+1] <- tapply(an,seqperiod,sum)
      ansim[,"cold","abs",j,i,s+1] <- tapply(an[!indheat],factor(seqperiod[!indheat]),sum)
      ansim[,"heat","abs",j,i,s+1] <- tapply(an[indheat],factor(seqperiod[indheat]),sum)
      
    }
  }
}


################################################################################
# ESTIMATE AN IN EACH PERIOD RELATIVE TO CURRENT DAYS (2010-19)

ansim[,,"rel",,,] <- ansim[,,"abs",,,] - ansim[rep("2010-19",length(levels(seqperiod))),,"abs",,,]


################################################################################
# SUMMARIZE THE RESULTS 
# COMPUTE AN/AF (95%CI) IN THE ENSEMBLE, BY RANGE & PERIOD & RCP

# CREATE NEW OBJECTS TO STORE RESULTS
# We now create 4 new arrays (2 arrays to store the abs and rel AN, and another 2
#   to store the estimated rel and abs AF) with 4 dimensions each to store the 
#   ensemble estimates (average impacts across GCMs) with the empirical 
#   95% confidence intervals. 
# In this case, the 4 dimensions correspond to the 10-year period, the point estimate 
#   and the CI, the temperature range and the scenario. 

estci <- c("est","ci.l","ci.u")

anabs <- afabs <- anrel <- afrel <- array(NA,dim=c(length(levels(seqperiod)),
  length(estci),length(temprange),length(rcp)), 
  dimnames=list(levels(seqperiod),estci,temprange,names(rcp)))


# ATTRIBUTABLE NUMBERS 
# ABSOLUTE
anabs[,"est",,"RCP4.5"] <- apply(ansim[,,"abs",,"RCP4.5",1],1:2,mean)
anabs[,"ci.l",,"RCP4.5"] <- apply(ansim[,,"abs",,"RCP4.5",-1],1:2,quantile,0.025)
anabs[,"ci.u",,"RCP4.5"] <- apply(ansim[,,"abs",,"RCP4.5",-1],1:2,quantile,0.975)

anabs[,"est",,"RCP8.5"] <- apply(ansim[,,"abs",,"RCP8.5",1],1:2,mean)
anabs[,"ci.l",,"RCP8.5"] <- apply(ansim[,,"abs",,"RCP8.5",-1],1:2,quantile,0.025)
anabs[,"ci.u",,"RCP8.5"] <- apply(ansim[,,"abs",,"RCP8.5",-1],1:2,quantile,0.975)

# RELATIVE
anrel[,"est",,"RCP4.5"] <- apply(ansim[,,"rel",,"RCP4.5",1],1:2,mean)
anrel[,"ci.l",,"RCP4.5"] <- apply(ansim[,,"rel",,"RCP4.5",-1],1:2,quantile,0.025)
anrel[,"ci.u",,"RCP4.5"] <- apply(ansim[,,"rel",,"RCP4.5",-1],1:2,quantile,0.975)

anrel[,"est",,"RCP8.5"] <- apply(ansim[,,"rel",,"RCP8.5",1],1:2,mean)
anrel[,"ci.l",,"RCP8.5"] <- apply(ansim[,,"rel",,"RCP8.5",-1],1:2,quantile,0.025)
anrel[,"ci.u",,"RCP8.5"] <- apply(ansim[,,"rel",,"RCP8.5",-1],1:2,quantile,0.975)

# ATTRIBUTABLE FRACTION
afabs[,,,] <- anabs[,,,]/deathperiod*100
afrel[,,,] <- anrel[,,,]/deathperiod*100

# REMOVE ansim
rm(ansim)



################################################################################
# FIGURE 5

# To obtained Figure 5 of the manuscript, we restrict the projections to one 
#    RCP (RCP8.5), one GCM (NorESM1-M) and two periods ("2010-19","2090-99").

# RESTRICT MODELLED TEMP TO 2010 AND 2090, NorESM1-M
tmean2010 <- rcp8p5cal[substr(rcp8p5cal$date,1,3)%in%"201",6]
tmean2090 <- rcp8p5cal[substr(rcp8p5cal$date,1,3)%in%"209",6]

# PREDICTIONS
bvar <- do.call(onebasis,c(list(x=c(tmean2010,tmean2090)),argvar))
#red <- crossreduce(cb,m,cen=19)
pred_exp <- crosspred(bvar,coef=coef,vcov=vcov,model.link="log",
  cen=cen,by=0.1)

# RESTRICT MORTALITY DATA
death2010 <- deathproj[substr(rcp8p5$date,1,3)%in%"201"]
death2090 <- deathproj[substr(rcp8p5$date,1,3)%in%"209"]

# COMPUTE ATTRIBUTABLE DEATHS
breaks <- seq(min(c(tmean2010,tmean2090),na.rm=T),
  max(c(tmean2010,tmean2090),na.rm=T),length=30)
cenvec <- do.call(onebasis,c(list(x=cen),argvar))
bvar <- do.call(onebasis,c(list(x=tmean2010),argvar))
bvarcen <- scale(bvar,center=cenvec,scale=F)
anhist <- (1-exp(-bvarcen%*%coef(red)))*death2010
bvar <- do.call(onebasis,c(list(x=tmean2090),argvar))
bvarcen <- scale(bvar,center=cenvec,scale=F)
anproj <- (1-exp(-bvarcen%*%coef(red)))*death2090

# PLOT THE FIGURE 
#pdf("Figure 5.pdf",height=6,width=4)
layout(matrix(1:4,ncol=1),heights=c(0.6,rep(0.3,2),0.20))

par(mar=c(0.2,5,2,1),mgp=c(3,1,0),las=1,cex.axis=0.8,cex.lab=1)

plot(pred_exp,type="n",ylim=c(0.8,3),axes=F,lab=c(6,5,7),xlab="",ylab="RR",
  main=" ")
ind1 <- pred_exp$predvar<=min(tmean2010)
ind2 <- pred_exp$predvar>=min(tmean2010) & pred_exp$predvar<=cen
ind3 <- pred_exp$predvar>=cen & pred_exp$predvar<=max(tmean2010)
ind4 <- pred_exp$predvar>=max(tmean2010)
if(any(ind1)) lines(pred_exp$predvar[ind1],pred_exp$allRRfit[ind1],col=4,lwd=1.5,lty=2)
lines(pred_exp$predvar[ind2],pred_exp$allRRfit[ind2],col=4,lwd=1.5)
lines(pred_exp$predvar[ind3],pred_exp$allRRfit[ind3],col=2,lwd=1.5)
lines(pred_exp$predvar[ind4],pred_exp$allRRfit[ind4],col=2,lwd=1.5,lty=2)
axis(2,at=0:5*4*0.1+1)
abline(v=cen,lty=3)
abline(v=max(tmean2010),lty=2)

par(mar=c(0.1,5,0.3,1))

rgbcol <- col2rgb(c(grey(0.8),"lightgreen"))
col <- apply(rgbcol,2,
  function(x) do.call(rgb,c(as.list(x),alpha=255/2,max=255)))

d00 <- density(tmean2010)
d90 <- density(tmean2090)
plot(d00,type="n",bty="l",main="",axes=F,xlim=range(pred_exp$predvar),
  ylim=c(0,max(d00$y,d90$y)*1.05),ylab="Temperature\ndistribution",xlab="")
polygon(d00,col=col[1],lty=4)
polygon(d90,col=col[2],lty=5)
axis(2,at=pretty(c(0,max(d00$y,d90$y)*1.05),n=4))
legend("topleft",c("2010-19","2090-99"),pch=22,inset=0.005,cex=0.9,bty="n",
  pt.bg=col)
abline(v=cen,lty=3)
abline(v=max(tmean2010),lty=2)
abline(h=0)

w00 <- abs(anhist)/sum(death2010)
d00 <- density(tmean2010,weights=w00/sum(w00))
w90 <- abs(anproj)/sum(death2090)
d90 <- density(tmean2090,weights=w90/sum(w90))
plot(d00,type="n",bty="l",main="",axes=F,xlim=range(pred_exp$predvar),
  ylim=c(0,max(d00$y,d90$y)*1.05),ylab="Excess mortality (%)",xlab="")
axis(2,at=pretty(c(0,max(d00$y,d90$y)*1.05),n=4))
polygon(d00,col=col[1],lty=4)
polygon(d90,col=col[2],lty=5)
abline(v=cen,lty=3)
abline(v=max(tmean2010),lty=2)
abline(h=0)

par(mar=c(4,5,0,1),mgp=c(2.5,1,0))

plot(pred_exp$predvar,rep(0,length(pred_exp$predvar)),type="n",yaxt="n",ylab="",
  xlab=xlab,frame.plot=F)
axis(1,at=-50:50,tck=0,labels=F)

layout(1)
#dev.off()

# NB: uncomment pdf() and dev.off() for saving the plot in pdf format.

################################################################################
# CHANGE IN AF (%) BETWEEN 2010-19 & 2090-99 (RCP8.5)

#   NB: These figures corresponds to the ones reported in section 5 of the main 
#       manuscript.

# HEAT
afrel["2090-99",,"heat","RCP8.5"]

# COLD
afrel["2090-99",,"cold","RCP8.5"]

