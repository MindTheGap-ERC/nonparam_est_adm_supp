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


load(file = "data/res/site1266_data.RData")
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

petm_recovery_stats = function(){
  clay_int = c(base_clay, clay_layer_top)
  recovery_int = c(base_recovery, recovery_3_top)
  l = list()
  for (i in c("const_det", "inc_det", "dec_det")){
    adm = adm_list[[i]]
    l[[paste0(i, "_median_clay")]] = median(sapply(get_time(adm, h = rev(clay_int)), diff ))
    l[[paste0(i, "_iqr_clay")]] = IQR(sapply(get_time(adm, h = rev(clay_int)), diff ))
    l[[paste0(i, "_median_recovery")]] = median(sapply(get_time(adm, h = rev(recovery_int)), diff ))
    l[[paste0(i, "_iqr_recovery")]] = IQR(sapply(get_time(adm, h = rev(recovery_int)), diff ))
  }
  return(l)
}
petm_res = petm_recovery_stats()



dpi = 400
lab_size = 7
title_size = 8
fig_width_cm = 12
annot_size = 5
legend_size = 5
ax_size = 4

const_fl_col = "red"
inc_fl_col = "blue"
dec_fl_col = "black"
pre_col = "azure2"
main_col = "azure3"
rec_col = "azure4"
box_cols = c(pre_col, main_col, rec_col)
scenario_cols = c(const_fl_col, inc_fl_col, dec_fl_col)
lwd_env = 0.5
lwd_med = 1
med_lty = 1
env_lty = 6

## Plot: duration of clay layer and recovery interval
clay_int = c(base_clay, clay_layer_top)
recovery_int = c(base_recovery, recovery_3_top)
dur_clay_const = sapply(get_time(adm_list$const_det, h = rev(clay_int)), diff)
dur_clay_inc = sapply(get_time(adm_list$inc_det, h = rev(clay_int)), diff)
dur_clay_dec = sapply(get_time(adm_list$dec_det, h = rev(clay_int)), diff)
df = data.frame(duration = c(dur_clay_const, dur_clay_inc, dur_clay_dec),
                Scenario = c(rep("Constant flux", length(dur_clay_const)),
                             rep("Increasing flux", length(dur_clay_inc)),
                             rep("Decreasing flux", length(dur_clay_dec))))

clay_layer_dur = ggplot(df, aes(x = duration, fill = Scenario)) +
  geom_density(alpha = 0.5) +
  xlab("Duration [kyr[") +
  ylab("Density") +
  ggtitle("Clay duration") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.9),
        axis.text = element_text(size = ax_size),
        axis.title = element_text(size = title_size),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = legend_size),
        legend.title = element_blank(),
        plot.title = element_text(size = title_size))

dur_rec_const = sapply(get_time(adm_list$const_det, h = rev(recovery_int)), diff)
dur_rec_inc = sapply(get_time(adm_list$inc_det, h = rev(recovery_int)), diff)
dur_rec_dec = sapply(get_time(adm_list$dec_det, h = rev(recovery_int)), diff)
df = data.frame(duration = c(dur_rec_const, dur_rec_inc, dur_rec_dec),
                Scenario = c(rep("Constant flux", length(dur_rec_const)),
                             rep("Increasing flux", length(dur_rec_inc)),
                             rep("Decreasing flux", length(dur_rec_dec))))

petm_rec = ggplot(df, aes(x = duration, fill = Scenario)) +
  geom_density(alpha = 0.5) +
  ggtitle("PETM recovery duration") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.9),
        axis.text = element_text(size = ax_size),
        axis.title = element_text(size = title_size),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = legend_size),
        legend.title = element_blank(),
        plot.title = element_text(size = title_size))
  

plt = egg::ggarrange(clay_layer_dur, petm_rec, nrow = 1, ncol = 2, labels = LETTERS[1:2])


## Sed rate plot

a = median_adm(adm_list$const_det, h = h_eval)
plot(a)
sedr = sed_rate_l(a, h_eval)
plot(h_eval, sedr, type = "l")

sedr_const = adm_list$const_det |> median_adm(h = h_eval) |> sed_rate_l(h_eval)
sedr_dec = adm_list$dec_det |> median_adm(h = h_eval) |> sed_rate_l(h_eval)
sedr_inc = adm_list$inc_det |> median_adm(h = h_eval) |> sed_rate_l(h_eval)

df = data.frame(h = rep(h_eval, 3),
                sedr = c(sedr_const, sedr_dec, sedr_inc) * 100,
                Scenario = c(rep("Constant flux", length(sedr_const)),
                             rep("Decreasing flux", length(sedr_dec)),
                             rep("Increasing flux", length(sedr_inc))))
sedr_plot = ggplot(df, aes(x = h, y = sedr, col = Scenario)) +
  geom_line() +
  xlab("Meters below sea floor [m]") +
  ylab("Sedimentation rate [cm/kyr]") +
  ggtitle("Sedimentation Rate")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.9),
        axis.text = element_text(size = ax_size),
        axis.title = element_text(size = title_size),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = legend_size),
        legend.title = element_blank(),
        plot.title = element_text(size = title_size))

cond_const = adm_list$const_det |> median_adm(h = h_eval) |> condensation(h_eval)
cond_dec = adm_list$dec_det |> median_adm(h = h_eval) |> condensation(h_eval)
cond_inc = adm_list$inc_det |> median_adm(h = h_eval) |> condensation(h_eval)

df = data.frame(h = rep(h_eval, 3),
                cond = c(cond_const, cond_dec, cond_inc) / 100,
                Scenario = c(rep("Constant flux", length(cond_const)),
                             rep("Decreasing flux", length(cond_dec)),
                             rep("Increasing flux", length(cond_inc))))

cond_plot = ggplot(df, aes(x = h, y = cond, col = Scenario)) +
  geom_line() +
  xlab("Meters below sea floor [m]") +
  ylab("Condensation [kyr/cm]") +
  ggtitle("Condensation")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.9),
        axis.text = element_text(size = ax_size),
        axis.title = element_text(size = title_size),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = legend_size),
        legend.title = element_blank(),
        plot.title = element_text(size = title_size))

egg::ggarrange(sedr_plot, cond_plot, nrow = 1, ncol = 2, labels = LETTERS[1:2])

## age-depth model plots

const_med = quantile_adm(adm_list$const_det, h_eval, 0.5)
const_p1 = quantile_adm(adm_list$const_det, h_eval, 0.975)
const_p2 = quantile_adm(adm_list$const_det, h_eval, 0.025)

inc_med = quantile_adm(adm_list$inc_det, h_eval, 0.5)
inc_p1 = quantile_adm(adm_list$inc_det, h_eval, 0.975)
inc_p2 = quantile_adm(adm_list$inc_det, h_eval, 0.025)

dec_med = quantile_adm(adm_list$dec_det, h_eval, 0.5)
dec_p1 = quantile_adm(adm_list$dec_det, h_eval, 0.975)
dec_p2 = quantile_adm(adm_list$dec_det, h_eval, 0.025)

df = data.frame(he = rep(h_eval, 9),
                val = -c(const_med$t, const_p1$t, const_p2$t, inc_med$t, inc_p1$t, inc_p2$t, dec_med$t, dec_p1$t, dec_p2$t),
                Scenario = c(rep("Constant flux", 3 * length(const_med$t)), rep("Increasing flux", 3 * length(inc_med$t)), rep("Decreasing flux", 3 * length(dec_med$t))),
                Type = rep(rep(c("B","A", "A"), 3), each = length(h_eval)),
                ind = rep(LETTERS[1:9], each = length(h_eval)))

ggplot(df, aes(y = he, x = val, color = Scenario, group = ind)) +
  geom_line() +
  xlab("Time [kyr]") +
  ylab("Height [m]") +
  ggtitle("PETM age-depth models")
