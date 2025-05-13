#### Global plotting options ####

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
core_col = "azure1"
rec1_col = "azure2"
rec2_col = "azure3"
rec4_col = "azure4"
box_cols = c(core_col, rec1_col, rec2_col, rec4_col)
scenario_cols = c(const_fl_col, inc_fl_col, dec_fl_col)
lwd_env = 0.75
lwd_med = 1.5
med_lty = 1
env_lty = 6


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


mean_3He_flux = 0.48 # mu based on supplementary materials
twosigma_3H3_flux = 0.08 # 2 sigma based on supplementary materials

#### Three flux scenarios ####

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

# 3He flux observed in the strat domain with measurement error
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

#### construct adms for 6 cases: ####
# increasing, decreasing, and constant flux, with and without error in 3He flux in depth domain

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

#### Load and plot already estimated data ####

load(file = "data/res/site1266_data.RData")

##### Reference points and intervals #####
clay_layer_top = 306.15
base_recovery = 306.4
recovery_1_top = 306.15 # “Shoulder” δ13C inflection point F
recovery_2_top = 304.7 # δ13C inflection point G
recovery_3_top = 304.19 # End of anomalously high carbonate sedimentation

core_interval <- c(base_recovery, base_clay)
recovery1_interval <- c(recovery_1_top, base_recovery)
recovery2_interval <- c(recovery_2_top, recovery_1_top)
recovery3_interval <- c(recovery_3_top, recovery_2_top)


for (i in names(adm_list)){
  plot(adm_list[[i]])
}

iqr_dur = c()
clay_dur = c()
med = c()
for (i in names(adm_list)){
  adm = adm_list[[i]]
  aa = get_time(adm, h = c(recovery_3_top, base_clay))
  hist(sapply(aa, diff), main = i)
  clay_dur[i] = IQR(sapply(aa, diff))
  med[i] = median(sapply(aa, diff))
}

source("code/petm_recovery_stats.R")
petm_res = petm_recovery_stats()


##### Duration of clay layer #####

clay_int = c(base_clay, clay_layer_top)
dur_clay_const = sapply(get_time(adm_list$const_det, h = rev(clay_int)), diff)
dur_clay_inc = sapply(get_time(adm_list$inc_det, h = rev(clay_int)), diff)
dur_clay_dec = sapply(get_time(adm_list$dec_det, h = rev(clay_int)), diff)
df = data.frame(duration = c(dur_clay_const, dur_clay_inc, dur_clay_dec),
                Scenario = c(rep("Constant flux", length(dur_clay_const)),
                             rep("Increasing flux", length(dur_clay_inc)),
                             rep("Decreasing flux", length(dur_clay_dec))))

clay_layer_dur = ggplot(df, aes(x = duration, fill = Scenario)) +
  geom_density(alpha = 0.5) +
  xlab("Duration [kyr]") +
  ylab("Density") +
  ggtitle("Clay duration") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.85),
        axis.text = element_text(size = ax_size),
        axis.title = element_text(size = title_size),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = legend_size),
        legend.title = element_blank(),
        plot.title = element_text(size = title_size))

##### Duration of the recovery layer #####

recovery_int = c(base_recovery, recovery_3_top)
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
  xlab("Duration [kyr]") +
  ylab("Density") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.85),
        axis.text = element_text(size = ax_size),
        axis.title = element_text(size = title_size),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = legend_size),
        legend.title = element_blank(),
        plot.title = element_text(size = title_size))
  

plt = egg::ggarrange(clay_layer_dur, petm_rec, nrow = 1, ncol = 2, labels = LETTERS[1:2])

ggsave(filename = "figs/site1266_joint_duration.png",
       plot = plt,
       width = fig_width_cm,
       unit = "cm",
       height = 6)

##### Duration of PETM #####

PETM_int = c(base_clay, recovery_3_top)
dur_PETM_const = sapply(get_time(adm_list$const_det, h = rev(PETM_int)), diff)
dur_PETM_inc = sapply(get_time(adm_list$inc_det, h = rev(PETM_int)), diff)
dur_PETM_dec = sapply(get_time(adm_list$dec_det, h = rev(PETM_int)), diff)
df = data.frame(duration = c(dur_PETM_const, dur_PETM_inc, dur_PETM_dec),
                Scenario = c(rep("Constant flux", length(dur_PETM_const)),
                             rep("Increasing flux", length(dur_PETM_inc)),
                             rep("Decreasing flux", length(dur_PETM_dec))))

aggregate(df$duration, by=list(df$Scenario), FUN = median)

#### Determine sedimentation rate ####

# convert into cm/kyr
source("code/median_sed_rate_l.R")
sedr_const =   100 * median_sed_rate_l(adm_list$const_det, h_eval)
sedr_inc =  100 * median_sed_rate_l(adm_list$inc_det, h_eval)
sedr_dec =  100 * median_sed_rate_l(adm_list$dec_det, h_eval)

sedr_stats = list("sedr_range_const"= range(sedr_const),
                  "sedr_range_inc" = range(sedr_inc),
                  "sedr_range_dec" = range(sedr_dec),
                  "sedr_fac_const" = max(sedr_const)/ min(sedr_const),
                  "sedr_fac_ind" = max(sedr_inc)/ min(sedr_inc),
                  "sedr_fac_dec" = max(sedr_dec)/ min(sedr_dec),
                  "cond_range_const" = range(1/sedr_const),
                  "cond_range_inc" = range(1/sedr_inc),
                  "cond_range_dec" = range(1/sedr_dec))

##### sedimentation rate plot #####

source("code/sed_rate_plot.R")
sedr_plot = sed_rate_plot("site1266_sedrate")
ggsave("figs/site1266_sedrate.png",
       plot = sedr_plot,
       width = fig_width_cm,
       unit = "cm",
       height = 6)
#### Condensation plot ####

source("code/condensation_plot.R")
cond_plot = condensation_plot("site1266_condensation")
ggsave("figs/site1266_condensation.png",
       plot = cond_plot,
       width = fig_width_cm,
       unit = "cm",
       height = 6)

#### ADM plot ####

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
  annotate("rect", ymin = recovery_2_top, ymax = recovery_3_top, xmin = 0, xmax = 350, alpha=0.6, fill = "lightgray") +
  annotate("rect", ymin = recovery_1_top, ymax = recovery_2_top, xmin = 0, xmax = 350, alpha=0.6, fill = "gray") +
  annotate("rect", ymin = base_recovery, ymax = recovery_1_top, xmin = 0, xmax = 350, alpha=0.6, fill = "darkgray") +
  annotate("rect", ymin = base_clay, ymax = base_recovery, xmin = 0, xmax = 350, alpha=0.6, fill = "white") +
  annotate("text", y = recovery_3_top+(recovery_2_top-recovery_3_top)/2, x = 70, label = "Recovery III") +
  annotate("text", y = recovery_2_top+(recovery_1_top-recovery_2_top)/2, x = 70, label = "Recovery II") +
  annotate("text", y = recovery_1_top+(base_recovery-recovery_1_top)/2, x = 70, label = "Recovery I") +
  annotate("text", y = base_recovery+(base_clay-base_recovery)/2, x = 250, label = "Core") +
  geom_line(aes(linetype = Type, size = ind)) +
  scale_size_manual(values = c("A" = lwd_med, "B" = lwd_env, "C" = lwd_env, "D" = lwd_med, "E" = lwd_env, "F" = lwd_env, "G" = lwd_med, "H" =  lwd_env, "I" = lwd_env), guide = "none") +
  scale_linetype_manual(values = c("B" = med_lty, "A" = env_lty), guide = "none") +
  xlab("Relative time since the beginning of PETM [kyr]") +
  ylab("Depth [m composite depth]") +
  ggtitle("Age-depth models for PETM at ODP Site 1266") +
  scale_y_reverse()

ggsave("figs/site1266_adm.png",
       plot = last_plot(),
       width = fig_width_cm,
       dpi = dpi,
       unit = "cm",
       height = 10)

