#############################################################
### Building a PLSR model to predict traits from spectra ####
### Anna K. Schweiger, April 2017 

# based on: Couture, J.J., Singh, A., Rubert‐Nason, K.F., Serbin, S.P., Lindroth, R.L. 
# and Townsend, P.A. (2016), Spectroscopic determination of ecologically relevant 
# plant secondary metabolites. Methods Ecol Evol, 7: 1402-1412. 
# doi:10.1111/2041-210X.12596

library(pls)
library(reshape)
library(agricolae)

dat <- read.csv("./example_dat.csv") 
###  Spectra are vector normalized using R package spectrolab:
### Meireles, J. E., Schweiger, A. K. & Cavender-Bares, J. spectrolab: Class and Methods for Hyperspectral Data. R package version 0.0.2.
### spectrolab::normalize()

# select wavelength range and sampling interval for model development
# e.g., NSC: 1200-2400 nm, 20 nm wvl interval
wvl <- names(dat)[grep("X", names(dat))] # wavelenghts
inVars <- names(dat)[1] # variable(s) to be predicted, example: non-structural carbohydrates (NSC)
range <- wvl[which(wvl=="X1200"):which(wvl=="X2400")]
inBands <- range[seq(1,1201,20)]

VIPjh=function(object, j, h)
{
  if (object$method != "oscorespls") stop("Only implemented for orthogonal scores algorithm. Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1) stop("Only implemented for single-response models")
  b=c(object$Yloadings)[1:h]
  T=object$scores[,1:h, drop = FALSE]
  SS=b^2 * colSums(T^2)
  W=object$loading.weights[,1:h, drop = FALSE]
  Wnorm2=colSums(W^2)
  VIP=sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
  return(VIP)
}

###### Part 1: Find number of components 
set.seed(1840)
for (i in 1:length(inVars)){
  inVar <- inVars[i]
  iForm <- paste(inBands,collapse="+")
  iForm <- as.formula(paste(inVar,"~",iForm))
  dir.create(paste("~/Desktop/", inVar, sep=""),recursive = T) # create folder to save output
  
  ### find number of components
  nCut <- floor(0.8*nrow(dat)) # 80% data for calibration, can be adjusted
  nComps <- 15 # maximum no of components, can be adjusted
  nsims <- 50 # no of simulations, can be adjusted/increased
  outMat <- matrix(data=NA,nrow=nsims,ncol=nComps)
  outMatRM <- matrix(data=NA,nrow=nsims,ncol=nComps)
  
  for (nsim in seq(nsims)){
    print(nsim)
    flush.console()
    dat$RAND <- order(runif(nrow(dat)))
    subData <- dat[dat["RAND"]<nCut,]
    resNCOMP <- plsr(iForm,data=subData,ncomp=nComps,validation="LOO",method="oscorespls")
    resPRESS <- as.vector(resNCOMP$validation$PRESS)
    outMat[nsim,seq(resNCOMP$validation$ncomp)] <-resPRESS
    resRMSEP <- as.numeric(RMSEP(resNCOMP,estimate="CV",intercept=F)$val)
    outMatRM[nsim,seq(resNCOMP$validation$ncomp)] <-resRMSEP
  }
  
  ### PRESS stat: Tukey test and plot
  pressDF <- as.data.frame(outMat)
  names(pressDF) <- as.character(seq(nComps))
  pressDFres <- melt(pressDF)
  
  modi <- lm (value~variable, pressDFres) 
  tuk <- HSD.test (modi,"variable")
  tuk_dat <- as.data.frame(tuk$groups)
  tuk_dat$var <- as.numeric(row.names(tuk_dat))
  tuk_dat <- tuk_dat[order(tuk_dat$var,decreasing = F),]
  letters <- as.character(tuk_dat$groups)
  
  jpeg(paste("~/Desktop/", inVar, "/", inVar,"_PRESS.jpg",sep=""),
       width=6,height=5,units="in",res=200)
  par(bty="l")
  boxplot(pressDFres$value~pressDFres$variable,
          xlab="n Components",ylab="PRESS",main=inVar)
  text(x=1:max(as.numeric(pressDFres$variable)), 
       y=rep(max(pressDFres$value),15),letters)
  dev.off()
  
  ### RMSEP: Tukey test and plot
  RMDF <- as.data.frame(outMatRM)
  names(RMDF) <- as.character(seq(nComps))
  RMDFres <- melt(RMDF)
  
  modi <- lm (value~variable, RMDFres) 
  tuk <- HSD.test (modi,"variable")
  tuk_dat <- as.data.frame(tuk$groups)
  tuk_dat$var <- as.numeric(row.names(tuk_dat))
  tuk_dat <- tuk_dat[order(tuk_dat$var,decreasing = F),]
  letters <- as.character(tuk_dat$groups)
  
  jpeg(paste("~/Desktop/", inVar, "/", inVar,"_RMSEP.jpg",sep=""),
       width=6,height=5,units="in",res=200)
  par(bty="l")
  boxplot(RMDFres$value~RMDFres$variable,
          xlab="n Components",ylab="RMSEP",main=inVar)
  text(x=1:max(as.numeric(RMDFres$variable)), 
       y=rep(max(RMDFres$value),15),letters)
  dev.off()
}

saveData <- dat  

### check output: select no of comps using RMSEP/PRESS criterion and Tukey tests

#######################################
###### Part 2: Final models ###########

(inVar <- inVars[1]) ### one trait at a time ... 
dat <- saveData[complete.cases(saveData[,inVar]),]
nCompss <- c(8) ### number of comps, muliple options can be selected

iForm <- paste(inBands,collapse="+")
iForm <- as.formula(paste(inVar,"~",iForm))

set.seed(1840)
for (i in 1:length(nCompss)){
  nsims <- 50 # no of simulations, can be adjusted
  nComps <- nCompss[i] 
  nCut <- floor(0.8*nrow(dat)) # 80% data for calibration, can be adjusted
  
  coefMat <- matrix(data=NA,nrow=nsims,ncol=length(inBands)+1)
  coefStd <- matrix(data=NA,nrow=nsims,ncol=length(inBands)+1)
  vipMat <- matrix(data=NA,ncol=nsims,nrow=length(inBands))
  statMat <- matrix(data=NA,nrow=nsims,ncol=6)
  
  for (nsim in seq(nsims)){
    print(nsim)
    flush.console()
    dat$RAND <- order(runif(nrow(dat)))
    subData <- dat[dat["RAND"]<nCut,]
    tstData <- dat[dat["RAND"]>=nCut,]
    
    resX <- plsr(iForm,data=subData,ncomp=nComps,method="oscorespls") ###
    resS <- plsr(iForm,data=subData,ncomp=nComps,method="oscorespls",scale=T) # scaled by SD
    
    ### Coefficients (raw and standardized)
    coefs <- as.vector(coef(resX,ncomp=nComps,intercept=T))
    zcoef <- as.vector(coef(resS,ncomp=nComps,intercept=T))
    coefMat[nsim,] <- coefs
    coefStd[nsim,] <- zcoef # standardized coeffis for importance of wvls
    
    ### VIP
    vip <- c()
    for (j in seq(length(inBands))){
      vip <- c(vip,VIPjh(resS,j,nComps))
    }
    vipMat[,nsim] <- vip
    
    ### Model stats
    fitX <- as.vector(unlist(resX$fitted.values[,1,nComps]))
    preX <- as.vector(unlist(predict(resX,ncomp=nComps,newdata=tstData)))
    fitBias <- mean(subData[,inVar]-fitX)
    valBias <- mean(tstData[,inVar]-preX)
    fitR2 <- summary(lm(subData[,inVar]~fitX))$r.squared
    valR2 <- summary(lm(tstData[,inVar]~preX))$r.squared
    fitRMSE <- sqrt(mean((subData[,inVar]-fitX)^2))
    valRMSE <- sqrt(mean((tstData[,inVar]-preX)^2))
    outVec <- c(fitR2,fitRMSE,fitBias,valR2,valRMSE,valBias)
    statMat[nsim,] <- outVec
  }
  
  statMat <- as.data.frame(statMat)
  names(statMat) <- c("fitR2","fitRMSE","fitBias","valR2","valRMSE","valBias")
  write.csv(statMat,paste("~/Desktop/", inVar, "/", inVar,"_", nComps,"comps_stats.csv", sep=""),row.names=FALSE)
  
  ### Coefficients
  coeffis <- data.frame(matrix(nrow = length(inBands)+1, ncol = 3))
  names(coeffis) <- c("bands", "mean","stdv")
  coeffis$bands <- c("Intercept",inBands)
  coeffis$mean <- apply(coefMat,MAR=2,FUN=mean)
  coeffis$stdv <- apply(coefMat,MAR=2,FUN=sd)
  write.table(coeffis,paste("~/Desktop/", inVar, "/", inVar, "_",nComps,"comps_coeffMEAN.csv", sep=""), 
              sep=",",col.names=T,row.names=F)
  
  specMat <- dat[,inBands]
  specMat <- cbind(rep(1,nrow(specMat)),specMat)
  specMat <- as.matrix(specMat)
  
  predMat <- specMat%*%t(coefMat)
  predMean <- apply(predMat,FUN=mean,MAR=1)
  predStdv <- apply(predMat,FUN=sd,MAR=1)
  
  preds <- dat[, !(names(dat) %in% inBands|names(dat)=="RAND")]
  preds[,paste("predMean_",inVar,sep="")] <- predMean
  preds[,paste("predStdv_",inVar,sep="")] <- predStdv
  write.csv(preds,paste("~/Desktop/", inVar, "/", inVar, "_", nComps,"comps_preds.csv", sep=""), row.names=FALSE)
  
  ### Plot Predictions 
  modCI <- quantile(statMat$fitR2, probs=c(0.05,0.95))
  
  formi <- as.formula(paste(paste(inVar," ~ predMean_",inVar, sep="")))
  lmformi <-   as.formula(paste(paste("predMean_",inVar," ~ ", inVar, sep="")))
  
  jpeg(paste("~/Desktop/", inVar,  "/", inVar,"_", nComps, "comp_predplot.jpg",sep=""),
       width=6,height=5,units="in",res=300)
  plot(formi, data= preds, pch=16,cex=0.8,ylab="measured",xlab="predicted",
       main=inVar,xlim=c(min(predMean-predStdv),max(predMean+predStdv)))
  abline(lm(lmformi, data= preds))
  abline(a = 0,b = 1, lty=2)
  arrows(predMean,dat[,inVar],predMean+predStdv,dat[,inVar],angle=90,length=0.05, lwd=0.8)
  arrows(predMean,dat[,inVar],predMean-predStdv,dat[,inVar],angle=90,length=0.05,lwd=0.8)
  legend("topleft", bty="n", cex=0.8,
         c(paste("R² = ", sprintf("%.2f",signif (mean(statMat$fitR2),3))," [",signif(modCI[1],2),",",signif(modCI[2],2),"]", sep = ""), 
           paste("RMSEP =",sprintf("%.2f",signif (mean(statMat$fitRMSE),3)), sep=" "),
           paste("ncomps =", nComps, sep=" ")))
  dev.off()
  
  ### VIP for plotting
  vipAggr <- as.data.frame(t(apply(vipMat,MAR=1,FUN=quantile,probs=c(0.05,0.5,0.95))))
  vipAggr$mean_VIP <- apply(vipMat,MAR=1,FUN=mean)
  vipAggr$stdv <- apply(vipMat,MAR=1,FUN=sd)
  serr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))  ## standard error of mean
  vipAggr$se <- apply(vipMat,MAR=1,FUN=serr)
  vipAggr$band <- inBands
  
  ### Standardized coefficients for plotting
  coeff_std <- data.frame(matrix(nrow = length(inBands)+1, ncol = 3))
  names(coeff_std) <- c("bands", "mean","stdv")
  coeff_std$bands <- c("Intercept",inBands)
  coeff_std$mean <- apply(coefStd,MAR=2,FUN=mean)
  coeff_std$stdv <- apply(coefStd,MAR=2,FUN=sd)
  
  ### Plot VIP and standardized coefficients
  jpeg(paste("~/Desktop/", inVar,  "/", inVar, "_", nComps, "comp_varimp.jpg",sep=""),
       width=6,height=7,units="in",res=300)
  par(mfrow=c(2,1), mar=c(1.5,4,2.5,1.5), oma=c(3,0,0,0))
  plot(coeff_std$mean[-1]~as.numeric(substr(coeff_std$bands[-1],2,nchar(coeff_std$bands[-1]))),
       type="p",pch=19, xlab="",ylab="coeff_stdmean",main=paste(inVar,nComps,"comps",sep = "_"),
       ylim=c(-max(abs(coeff_std$mean[-1])),max(abs(coeff_std$mean[-1]))), bty="l")
  abline(h=0)
  points(abs(coeff_std$mean)[-1]~as.numeric(substr(coeff_std$bands[-1],2,nchar(coeff_std$bands[-1]))),
         xlab="wvl",ylab="coeff_stdmean", col=2, pch=16, cex=0.8)
  
  lines(abs(coeff_std$mean)[-1]~as.numeric(substr(coeff_std$band[-1],2,nchar(coeff_std$band[-1]))), col=2)
  
  plot(as.numeric(substr(vipAggr$band,2,nchar(vipAggr$band))),vipAggr$mean_VIP, type="l",
       xlab = "wvl",ylab = "VIP", bty="l")
  polygon(x=c(as.numeric(substr(vipAggr$band,2,nchar(vipAggr$band))),
              rev(as.numeric(substr(vipAggr$band,2,nchar(vipAggr$band))))),
          y=c(vipAggr$mean_VIP+vipAggr$stdv*1.96, rev(vipAggr$mean_VIP-vipAggr$stdv*1.96)),  
          col =  adjustcolor("red", alpha.f = 0.2), border = NA)
  mtext("wavelength(nm)",1,outer = T,line = 1)
  dev.off()
}

### END: check output