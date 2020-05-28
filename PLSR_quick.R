###################################################################
### Building a QUICK PLSR model to predict traits from spectra ####
### to see if data are useful, without iterations, uncertainty calcs
### use PLSR_iter for final models ################################
### Anna Schweiger March 2016 

library(pls)

dada <- read.csv("./example_data_for_model.csv")
###  Spectra are vector normalized using R package spectrolab:
### Meireles, J. E., Schweiger, A. K. & Cavender-Bares, J. spectrolab: Class and Methods for Hyperspectral Data. R package version 0.0.2.
### spectrolab::normalize()

# select wavelength range and sampling interval for model development
# e.g., NSC: 1200-2400 nm, 20 nm wvl interval
wvl <- names(dada)[grep("X", names(dada))] # wavelenghts
range <- wvl[which(wvl=="X1200"):which(wvl=="X2400")]
inBands <- range[seq(1,1201,20)]

### set up new data matrix
dat <- data.frame (SampN = c(1:nrow(dada)))
dat$NSC <- dada$nonstructural_perc ## select trait
dat$spec <- as.matrix(dada[,inBands]) ## add spectral data
dat <- na.exclude(dat) 

### PLSR, select maximum no of comps and internal validation method 
mod <- plsr (NSC ~ spec, ncomp=15, data = dat, validation = "LOO") 
 
validationplot(mod, val.type="RMSEP", estimate="CV")
validationplot(mod, val.type="R2", estimate="CV")

compi <- 12 ## select no of components
summary(mod)

### Modelstats ################################
a <- predplot(mod,ncomp= compi, which="validation") # validationplot
b <- predplot(mod, ncomp= compi, which="train") # calibrationplot

dat_predtrue <- data.frame(abbrev=character(nrow(dat)))
dat_predtrue$abbrev <- dat$abbrev
dat_predtrue$measured <- a[, colnames(a)=="measured"]
dat_predtrue$predicted_val <- a[, colnames(a)=="predicted"]
dat_predtrue$predicted_train <- b[, colnames(b)=="predicted"]

### R2 validation ####
mod_fit<- lm(measured ~ predicted_val, data=dat_predtrue)
summary(mod_fit) 

plot (measured ~ predicted_val, data=dat_predtrue)
abline(mod_fit) # regression line predicted vs measured
abline(0,1,lty=2) # 1:1 line

### Test if slope is sign diff from 1
mod1 <- nls (measured ~ k*predicted_val+d, data= dat_predtrue, start=list(k=1, d=0)) ## mod1 ... alternative mod, fixed slope at 1 and intercept at 0
mod0 <- nls (measured ~ predicted_val+d, data= dat_predtrue,start=list(d=0)) ## mod0 ... our model, but intercept fixed at 0 (only look at deviation of slope)
anova (mod1,mod0) 

### Test if intercept is sign diff from 0
mod0d<- nls (measured ~ k*predicted_val, data= dat_predtrue, start=list(k=1)) ## mod0 ... our model, but slope fixed at 1 (only look at dev of intercept)
anova (mod1,mod0d) ## ideally no sign. diff

### Coefficients for predictions #####
coeffi <- coef(mod, ncomp=compi, intercept = T)

### Model stats ##########################
### RMSEP
rmse <- function(obs, pred) {sqrt(mean((obs-pred)^2))} 
(rm<- rmse(dat_predtrue$measured, dat_predtrue$predicted_val))

### Theil's U 
a <- sqrt(mean(dat_predtrue$measured^2))
b <- sqrt(mean(dat_predtrue$predicted_val^2))
(theil_u <- rm/(a+b))

### % samples within x % mean prediction error 
mpe <- (abs(dat_predtrue$predicted_val-dat_predtrue$measured))/(mean(dat_predtrue$measured))
((sum(mpe < 0.2))/nrow(dat)) ### % samples within 20% MPE (% can be adjusted)

#### Loadings plot 
wvl <- as.numeric(substr(colnames(dat$spec), 2, nchar(colnames(dat$spec))))
loadingplot(mod,comps = compi,xaxt="n", xlab="wavelength")
axis(1, seq(0, length(wvl),by = 10), labels = seq(wvl[1],wvl[length(wvl)], length.out = 7))

### END