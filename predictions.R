######################################################
####### Predict traits from PLSR models ##############
### Anna Schweiger March 2016 
library(spectrolab)

### example spectra NOT vector normalized
spec <- read.csv("./example_data_for_preds.csv", as.is = T)

### make spectra object
colnames(spec) <- gsub("X","", colnames(spec)) # colnames need to represent wavelengths
spec_ex <- as.spectra(spec)
plot_interactive(spec_ex) # plot data

spec_n <- normalize(spec_ex) # vector normalize

### read coefficients 
coefN <- read.csv("./PLSR/N_percent/12comps_coeffMEAN.csv", check.names = F)

#### check wavelength range and spectral resolution used in model
spec_res <- resample(spec_n,seq(1100,2400,20)) # resample such that new data matches wvl range and resolution of model
wavelengths(spec_res) # check result
plot(spec_res[1:5,])

#### predict trait and save as metadata
meta(spec_res, "nitrogen_perc") <- reflectance(spec_res) %*% coefN$mean[-1] + coefN$mean[1] # matrix multiplication plus intercept

### check output and convert to dataframe
meta(spec_res, "nitrogen_perc")
dat <- as.data.frame(spec_res)

### END 
