# load data from Murphy et al. https://doi.org/10.1016/j.gca.2010.03.039
data = read.csv("data/raw/murphy_et_al_2010_1-s2.0-S0016703710003108-mmc3.csv", header = TRUE, sep = "\t")

# heights
h = data$Depth..mcd.


base_clay = 306.78 #Murphy et al., Table 1

library(admtools)
h_eval = seq(303.5, base_clay, by = 0.01) # heights where the ages are determined - cm resolution
subdiv = 10000 # numeric options for integration


# stratigraphic tie point at the base of the clay layer
h_tp = function(){
  return(base_clay)
}
# time tie point: measure time relative to the deposition of the base of the clay layer
t_tp = function(){
  return(0)
}

## 3He fluxes in the time domain
# based on the number in the main text
# seems to be a typo, not used further
time_const_gen_text = function(){
  eps = 0.0001
  r = max(eps,rnorm(1, mean = 0.37, sd = 0.06/2)) # murphy et al 2010, pcc/cm^-2/kyr
  f = approxfun(x = c(-1000, 1000), y = rep(r, 2), rule = 2)
  return(f)
}
# based on the number in the suppl. materials
mean_3He_flux = 0.48 # mu based on supplementary materials
twosigma_3H3_flux = 0.08 # 2 sigma based on supplementary materials
# constant flux in time domain
time_const_gen_supp = function(){
  eps = 0.0001
  r = max(eps,rnorm(1, mean = mean_3He_flux, sd = twosigma_3H3_flux/2)) # murphy et al 2010, pcc/cm^-2/kyr
  f = approxfun(x = c(-1000, 1000), y = rep(r, 2), rule = 2)
  return(f)
}
# increasing flux in time domain
time_inc_gen_supp = function(){
  eps = 0.0001
  r = max(eps,rnorm(1, mean = mean_3He_flux, sd = twosigma_3H3_flux/2)) # murphy et al 2010, pcc/cm^-2/kyr
  f = approxfun(x = c(-1000,0, 1000), y = c(0.5 * r, r, 2*r), rule = 2)
  return(f)
}
# decreasing flux in time domain
time_dec_gen_supp = function(){
  eps = 0.0001
  r = max(eps,rnorm(1, mean = mean_3He_flux, sd = twosigma_3H3_flux/2)) # murphy et al 2010, pcc/cm^-2/kyr
  f = approxfun(x = c(-1000,0, 1000), y = c(2 * r, r, 0.5*r), rule = 2)
  return(f)
}


# 3H3 flux observed in the stratigraphic domain assuming no error in 3He measurement
strat_cont_gen_det = function(){
  f = approxfun(x = data$Depth..mcd.,
                y = data$X3HeET..pcc.g...20. * data$DBD..g.cm3. * 100,
                rule = 2)
  return(f)
}

strat_cont_gen_rand = function(){
  # assuming the 10 % error mentioned in the ms are 1 sigma
  eps = 0.000001 # small number cutoff to prevent negative flux values
  fl_mean = data$X3HeET..pcc.g...20.
  fl_sigma = 0.1  * data$X3HeET..pcc.g...20.
  fl = pmax(eps,stats::rnorm(length(fl_sigma), mean = fl_mean, sd = fl_sigma))
  f = approxfun(x = data$Depth..mcd.,
                y = fl * data$DBD..g.cm3. * 100,
                rule = 2)
  return(f)
}

# construct adms for 6 cases: increasing, decreasing, and constant flux, with and without error in 3He flux in depth domain
build_adms = function(){
  adm_list = list()
  adm = strat_cont_to_multiadm(h_tp = h_tp,
                               t_tp = t_tp,
                               strat_cont_gen = strat_cont_gen_det,
                               time_cont_gen = time_const_gen_supp,
                               h = h_eval,
                               subdivisions = subdiv, stop.on.error = FALSE,
                               no_of_rep = 500)
  adm_list[['const_det']] = adm
  adm = strat_cont_to_multiadm(h_tp = h_tp,
                               t_tp = t_tp,
                               strat_cont_gen = strat_cont_gen_rand,
                               time_cont_gen = time_const_gen_supp,
                               h = h_eval,
                               subdivisions = subdiv, stop.on.error = FALSE,
                               no_of_rep = 500)
  adm_list[['const_rand']] = adm
  cat("done with constant flux \n")
  # increasing flux
  adm = strat_cont_to_multiadm(h_tp = h_tp,
                               t_tp = t_tp,
                               strat_cont_gen = strat_cont_gen_det,
                               time_cont_gen = time_inc_gen_supp,
                               h = h_eval,
                               subdivisions = subdiv, stop.on.error = FALSE,
                               no_of_rep = 500)
  adm_list[['inc_det']] = adm
  adm = strat_cont_to_multiadm(h_tp = h_tp,
                               t_tp = t_tp,
                               strat_cont_gen = strat_cont_gen_rand,
                               time_cont_gen = time_inc_gen_supp,
                               h = h_eval,
                               subdivisions = subdiv, stop.on.error = FALSE,
                               no_of_rep = 500)
  adm_list[['inc_rand']] = adm
  cat("done with increasing flux\n")
  
  # decreasing flux
  adm = strat_cont_to_multiadm(h_tp = h_tp,
                               t_tp = t_tp,
                               strat_cont_gen = strat_cont_gen_det,
                               time_cont_gen = time_dec_gen_supp,
                               h = h_eval,
                               subdivisions = subdiv, stop.on.error = FALSE,
                               no_of_rep = 500)
  adm_list[['dec_det']] = adm
  adm = strat_cont_to_multiadm(h_tp = h_tp,
                               t_tp = t_tp,
                               strat_cont_gen = strat_cont_gen_rand,
                               time_cont_gen = time_dec_gen_supp,
                               h = h_eval,
                               subdivisions = subdiv, stop.on.error = FALSE,
                               no_of_rep = 500)
  adm_list[['dec_rand']] = adm
  cat("done with decreasing flux")
  return(adm_list)
}
# change to TRUE to run analysis - this might take a bit
if (TRUE){
  adm_list = build_adms()
}


save.image(file = "data/res/site1266_data.RData")

clay_layer_top = 306.15
base_recovery = 306.4
recovery_1_top = 306.15
recovery_2_top = 304.7
recovery_3_top = 304.19

for (i in names(adm_list)){
  plot(adm_list[[i]])
}

iqr_dur = c()
med = c()
for (i in names(adm_list)){
  adm = adm_list[[i]]
  aa = get_time(adm, h = c(recovery_3_top, base_clay))
  hist(sapply(aa, diff), main = i)
  clay_dur[i] = IQR(sapply(aa, diff))
  med[i] = median(sapply(aa, diff))
}
# 
# plot(adm)
# 
# aa = get_time(adm, h = c(306.78, 306.15))
# 
# x = sapply(aa, diff)
# 
# quantile(x)
# 
# hist(x)
