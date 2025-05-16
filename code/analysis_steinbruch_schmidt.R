#### Load packages ####
library(admtools)

#### Set seed ####
set.seed(42)

#### Constants ####
# from da silva et al 2020
targetE=c(130.7, 123.8, 98.9, 94.9)
astrochron::bergerPeriods(372, genplot = FALSE)
targetP=c(19.98, 16.87)

# from da silva https://doi.org/10.1038/s41598-020-69097-6, fig 2
# stratigraphic positions in m
h_ashbed = 1.55
h_top_lkw = 1.15
h_bottom_ukw = 3.6
h_top_ukw = 4.05
h_f_f_bdry = 4.05

# age constraints for the ff from da silva https://doi.org/10.1038/s41598-020-69097-6
# in Myr
ff_mean_dasilva = 371.870
ff_2sigma_dasilva = 0.108

# age constraints for the ff from Gradstein https://doi.org/10.1016/B978-0-12-824360-2.00022-X
# in Myr
ff_mean_gradstein = 371.1
ff_2sigma_gradstein = 1.1

# timing of ashbed based on Percival 2018 in kyr https://doi.org/10.1038/s41598-018-27847-7
# in kyr
t_ashbed_mean = 372.360 * 1000
uncertainty_radiometric <- 0.053 # [Ma], 2 sigma around the age
t_ashbed_sd = uncertainty_radiometric * 1000 / 2 # percival uses 2 sigma

# rate parameter for estimation
rate = 3
no_of_rep = 1000

#### Load data ####
# from da silva 2020, suppl material
data=read.csv("data/raw/SbS_XRF_forfactor3.csv",header = T, sep=";")

#### Clean data ####
Height=data[,c(1)]
mydata=data[,c(2:23)]
mydata[is.na(mydata)] <- 0
MS<-data[,c(1,23)]
MS=na.omit(MS)
MSi=astrochron::linterp(dat = MS, dt = 0.02, verbose = FALSE, genplot = FALSE)


#### run eTimeOpt ####
# taken from da silva 2020, suppl code https://doi.org/10.5281/zenodo.12516430
# test for precession amplitude modulation
etimeOptMS_prec=astrochron::eTimeOpt(MSi,
                                     win=0.02*100,
                                     step=0.02*10,
                                     sedmin=0.1,
                                     sedmax=0.6, 
                                     numsed=100,
                                     linLog=1,
                                     limit=T,
                                     fit=1,
                                     fitModPwr=T,
                                     flow=NULL,
                                     fhigh=NULL,
                                     roll=NULL,
                                     targetE=targetE,
                                     targetP=targetP,
                                     detrend=T,
                                     ydir=1,
                                     output=1,
                                     genplot=F,
                                     check=T,
                                     verbose=1)
# test for short eccentricty amplitude modulation
etimeOptMS_secc=astrochron::eTimeOpt(MSi,
                                     win=0.02*100,
                                     step=0.02*10,
                                     sedmin=0.1,
                                     sedmax=0.6,
                                     numsed=100,
                                     linLog=1,
                                     limit=T,
                                     fit=2,
                                     fitModPwr=T,
                                     flow=NULL,
                                     fhigh=NULL,
                                     roll=NULL,
                                     targetE=targetE,
                                     targetP=targetP,
                                     detrend=T,
                                     ydir=1,
                                     output=1,
                                     genplot=F,
                                     check=T,
                                     verbose=1)

par(mfrow = c(1,1))


#### Estimation procedure ####
#strat tie point: bentonite layer
h_tp = function(){
  return(h_ashbed)
}

# absolute age tie point with age from Percival et al 2018
t_tp_absolute = function(){
  t = rnorm(1, mean =  - t_ashbed_mean, sd = t_ashbed_sd)
  return(t)
}

#  absolute age tie point: U-Pb age without uncertainty
t_tp_anchor_no_error = function(){
  return(- t_ashbed_mean)
}

# extract r^2_opt from eTimeOpt
fa_prec = admtools::get_data_from_eTimeOpt(etimeOptMS_prec, index = 3)
# generate sed rate generator
se_prec = admtools::sed_rate_from_matrix(height = fa_prec$heights,
                               sedrate = fa_prec$sed_rate / 100, # convert cm/kyr to m/kyr
                               matrix = fa_prec$results,
                               mode = "poisson",
                               rate = rate)
fa_secc = admtools::get_data_from_eTimeOpt(etimeOptMS_secc, index = 3)
se_secc = admtools::sed_rate_from_matrix(height = fa_secc$heights,
                               sedrate = fa_secc$sed_rate / 100, # convert cm/kyr to m/kyr
                               matrix = fa_secc$results,
                               mode = "poisson",
                               rate = rate)

# heights of interest
h = seq(1.05, 4.25, by = 0.05)

#### estimate ADMs ####
adm_prec_abs = admtools::sedrate_to_multiadm(h_tp , 
                                             t_tp_absolute,
                                             sed_rate_gen = se_prec,
                                             h, 
                                             no_of_rep = no_of_rep )
adm_prec_abs_no_error = admtools::sedrate_to_multiadm(h_tp ,
                                                      t_tp_anchor_no_error,
                                                      sed_rate_gen = se_prec,
                                                      h, 
                                                      no_of_rep = no_of_rep )


adm_secc_abs = admtools::sedrate_to_multiadm(h_tp ,
                                             t_tp_absolute, 
                                             sed_rate_gen = se_secc,
                                             h = h, 
                                             no_of_rep = no_of_rep )
adm_secc_abs_no_error = admtools::sedrate_to_multiadm(h_tp ,
                                                      t_tp_anchor_no_error,
                                                      sed_rate_gen = se_secc,
                                                      h,
                                                      no_of_rep = no_of_rep )

save.image(file = "data/res/sbs_data.RData")


