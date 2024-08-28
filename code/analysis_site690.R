## load package ##
library(ggplot2)
library(admtools)


#### Set Seed ####
set.seed(42)

#### Constants ####
# taken from Farley and Eltgroth 2003, http://dx.doi.org/10.1016/S0012-821X(03)00017-7 , Fig 2
recovery_interval = c(-167, -169.65)
main_interval = c(-169.65, -170.65)
pre_interval = c(-170.65, -171.25)
max_obs = -166

obs_3He_flux = 0.69 # 3He flux
obs_3He_flux_error = 0.15 # % 2 sigma

eps = 0.001 # lowest possible 3He flux value


#### Constants for analysis ####
no_of_rep = 1000 # monte carlo iterations
subdivisions = 10000 # subdivisions of numeric integration

#### Load data & clean data ####
# from https://doi.org/10.1594/PANGAEA.723912 , supplement to Farley and Eltgroth 2003
data = read.csv(file = "data/raw/Farley_and_Eltgroth_2003_supp_data_1_site690.csv", header = TRUE, sep = ";")

# Exclude section 19H-5W-108-109 bc it has thickness 0
data = data[!data$core.section == "19H-5W-108-109",  ] 

# convert all length measurements to cm
data$mbsf_cm = data$mbsf * 100
data$mean.Depth..mbsf._cm = data$mean.Depth..mbsf. * 100
# determine start & end of bins
data["bin_start_cm"] =  (data$mean.Depth..mbsf._cm - 0.5 *data$Thickness..cm.)
data["bin_end_cm"] =  (data$mean.Depth..mbsf._cm + 0.5 *data$Thickness..cm.)
# convert mbsf into negative height
data$mbsf_cm = - data$mbsf_cm
data$mean.Depth..mbsf._cm = - data$mean.Depth..mbsf._cm
data$bin_end_cm = - data$bin_end_cm
data$bin_start_cm = -data$bin_start_cm

# tracer content per cm
data$ex_3He_pcc_p_cm = data$extraterrestrial.3He..pcc.g. * data$DBD..g.cm3.
data$ex_3He_pcc_p_cm_uncert = data$X3.He.... * data$DBD..g.cm3.

## convert everything to m
data$mbsf_m = data$mbsf_cm/100
data$bin_end_m = data$bin_end_cm/100
data$bin_start_m = data$bin_start_cm/100
data$ex_3He_pcc_p_m = data$ex_3He_pcc_p_cm * 100
data$ex_3He_pcc_p_m_uncert = data$ex_3He_pcc_p_cm_uncert * 100

binborders = sort(unique(c(data$bin_end_m, data$bin_start_m)), decreasing = TRUE)

#### Plot proxy ####
He = approxfun(x = binborders,
               y = c(data$ex_3He_pcc_p_m, data$ex_3He_pcc_p_m[length(data$ex_3He_pcc_p_m)]),
               rule = 2,
               f = 0,
               method = "constant")
He_uncert = approxfun(x = binborders,
                      y = c(data$ex_3He_pcc_p_m_uncert,data$ex_3He_pcc_p_m_uncert[length(data$ex_3He_pcc_p_m_uncert)]),
                      rule = 2,
                      f = 0,
                      method = "constant")
if (FALSE) { # plotting of tracer values
  plot(NULL,
       xlim = c(min(data$mbsf_m), max_obs),
       ylim = c(0, max(data$ex_3He_pcc_p_m) + max(data$ex_3He_pcc_p_m_uncert)),
       xlab = "Depth below sediment surface [m]",
       ylab = "ext. 3He [ppc / m]")
  h = seq(min(data$mbsf_m), max(data$mbsf_m), length.out = 10000)
  rect(ybottom = 0, ytop = 300, xleft = pre_interval[1], xright = pre_interval[2], col = "azure2")
  rect(ybottom = 0, ytop = 300, xleft = main_interval[1], xright = main_interval[2], col = "azure3")
  rect(ybottom = 0, ytop = 300, xleft = recovery_interval[1], xright = recovery_interval[2], col = "azure4")
  lines(h, He(h) + He_uncert(h), col = "red")
  lines(h, He(h) - He_uncert(h), col = "red")
  lines(h, He(h), col = "blue" )
  
  legend("topleft", lwd = 1, lty = 1, col = c("red", "blue"), legend = c("+/- 1 sd", "Mean"))
}


#### Estimate ADMs #### 

#### define adm helpers ####
mean_he = data$ex_3He_pcc_p_m
sd_he = data$ex_3He_pcc_p_m_uncert


#### Define proxy measurements in section ####
bin_borders = binborders
df = data.frame("sd" = sd_he, "mean" = mean_he)
st = strat_cont_gen_from_tracer(bin_borders = bin_borders, 
                               df = df, 
                               distribution = "normal")

#### Define Tie Points ####
t_tp = function(){ # fix time tie point at time 0
  return(0)
}

h_tp = function(){ # fix strat tie point at start of the main interval
  return(min(main_interval))
}

#### Define Proxy flux in time ####

time_cont_gen_const = function(){ # constant flux over 1 Myr based on the estimated tracer flux
  r = max(eps, rnorm(1, mean = obs_3He_flux  , sd = obs_3He_flux_error * obs_3He_flux / 2)) #pbb / kyr
  f = approxfun(x = c(-500, 500), y = c(r, r), rule = 2)
  return(f)
}

time_cont_gen_inc = function(){ # Flux increasing over 1 Myr
  r = max(eps, rnorm(1, mean = obs_3He_flux  , sd = obs_3He_flux_error * obs_3He_flux / 2)) #pbb / kyr
  f = approxfun(x = c(-500, 500), y = c(2/3  * r, 4/3  * r), rule = 2)
  return(f)
}
time_cont_gen_dec = function(){ # Flux decreasing over 1 Myr
  r = max(eps, rnorm(1, mean = obs_3He_flux  , sd = obs_3He_flux_error * obs_3He_flux / 2)) #pbb / kyr
  f = approxfun(x = c(-500, 500), y = c(4/3  * r, 2/3  * r), rule = 2)
  return(f)
}

#### Interpolation heights ####

h = seq(min(data$mbsf_m), max(data$mbsf_m), length.out = 100)
h = seq(min(data$mbsf_m), max_obs, length.out = 100)


#### run estimation ####
cat("determining age-depth models, please wait... \n")
cat("estimating constant flux adm\n")

my_adm_const = strat_cont_to_multiadm(h_tp = h_tp,
                                t_tp = t_tp,
                                strat_cont_gen = st,
                                time_cont_gen = time_cont_gen_const,
                                h = h,
                                no_of_rep = no_of_rep,
                                subdivisions = subdivisions,
                                T_unit = NULL,
                                L_unit = "m")
cat("estimating increasing flux adm\n")
my_adm_inc = strat_cont_to_multiadm(h_tp = h_tp,
                                      t_tp = t_tp,
                                      strat_cont_gen = st,
                                      time_cont_gen = time_cont_gen_inc,
                                      h = h,
                                      no_of_rep = no_of_rep,
                                      subdivisions = subdivisions,
                                      T_unit = NULL,
                                      L_unit = "m")
cat("estimating decreasing flux adm \n")

my_adm_dec = strat_cont_to_multiadm(h_tp = h_tp,
                                    t_tp = t_tp,
                                    strat_cont_gen = st,
                                    time_cont_gen = time_cont_gen_dec,
                                    h = h,
                                    no_of_rep = no_of_rep,
                                    subdivisions = subdivisions,
                                    T_unit = NULL,
                                    L_unit = "m")

cat("Creating plots for site 690 /n")

#### Condensation recovery interval to main interval ####
# main interval
dur_main_const =  - sapply(get_time(my_adm_const, h = main_interval), diff)
dur_main_dec =  - sapply(get_time(my_adm_dec, main_interval), diff)
dur_main_inc =  - sapply(get_time(my_adm_inc, main_interval), diff)
# recovery
dur_rec_const =  - sapply(get_time(my_adm_const, recovery_interval), diff)
dur_rec_inc = -  sapply(get_time(my_adm_inc, recovery_interval), diff)
dur_rec_dec =  - sapply(get_time(my_adm_dec, recovery_interval), diff)

# ratio main interval to recovery
time_ratio_const = dur_rec_const / dur_main_const
time_ratio_inc = dur_rec_inc / dur_main_inc
time_ratio_dec = dur_rec_dec / dur_main_dec

interval_cond_stats = list("thickness_ratio" = diff(recovery_interval)/ diff(main_interval),
                           "med_ratio_const" = 100 * (1 - median(time_ratio_const)),
                           "med_ratio_inc" = 100 * (1 - median(time_ratio_inc)),
                           "med_ratio_dec" = 100 * (1 - median(time_ratio_dec)))


#### Global plotting options ####

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

#### Plot ADMs ####

plot_adm_site_690 = function(file_name){
  const_med = quantile_adm(my_adm_const, h, 0.5)
  const_p1 = quantile_adm(my_adm_const, h, 0.975)
  const_p2 = quantile_adm(my_adm_const, h, 0.025)
  
  inc_med = quantile_adm(my_adm_inc, h, 0.5)
  inc_p1 = quantile_adm(my_adm_inc, h, 0.975)
  inc_p2 = quantile_adm(my_adm_inc, h, 0.025)
  
  dec_med = quantile_adm(my_adm_dec, h, 0.5)
  dec_p1 = quantile_adm(my_adm_dec, h, 0.975)
  dec_p2 = quantile_adm(my_adm_dec, h, 0.025)
  
  df = data.frame(he = rep(-h, 9),
                  val = c(const_med$t, const_p1$t, const_p2$t, inc_med$t, inc_p1$t, inc_p2$t, dec_med$t, dec_p1$t, dec_p2$t),
                  Scenario = c(rep("Constant flux", 3 * length(const_med$t)), rep("Increasing flux", 3 * length(inc_med$t)), rep("Decreasing flux", 3 * length(dec_med$t))),
                  Type = rep(rep(c("B","A", "A"), 3), each = 100),
                  ind = rep(LETTERS[1:9], each = 100))
  
  t_min = min(df$val)
  t_max = max(df$val)
  rect = data.frame(h_min = -c(pre_interval[1], main_interval[1], recovery_interval[1]),
                    h_max = -c(pre_interval[2], main_interval[2], recovery_interval[2]),
                    t_min = rep(t_min, 3),
                    t_max = rep(t_max, 3))
  
  
  ggplot(df, aes(y = he, x = val, color = Scenario,group = ind)) +
    geom_rect(data = rect, inherit.aes = FALSE, aes(ymin = h_min, ymax = h_max, xmin = t_min, xmax = t_max), fill = box_cols) +
    geom_line(aes(linetype = Type, size = ind)) +
    scale_color_manual(values = scenario_cols) +
    scale_size_manual(values = c("A" = lwd_med, "B" = lwd_env, "C" = lwd_env, "D" = lwd_med, "E" = lwd_env, "F" = lwd_env, "G" = lwd_med, "H" =  lwd_env, "I" = lwd_env), guide = "none") +
    scale_linetype_manual(values = c("B" = med_lty, "A" = env_lty), guide = "none") +
    xlab("Time since beginning of PETM main event [kyr]") +
    ylab("Meters below sea floor [m]") +
    scale_y_reverse() +
    annotate("text", x = mean(c(t_min, t_max)) - 150, y = mean(-recovery_interval), label = "recovery interval") +
    annotate("text", x = mean(c(t_min, t_max)) - 150, y = mean(-main_interval), label = "main interval") +
    annotate("text", x = mean(c(t_min, t_max)) - 150, y = mean(-pre_interval), label = "precursor interval") +
    ggtitle("Age-depth models for ODP Site 690") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.1, 0.9))
    ggsave(paste0("figs/", file_name, ".png"))
}

plot_adm_site_690("site690_adm")

#### Determine sedimentation rate ####

median_sed_rate_l = function(x, h){
  #' @title get median sed rate from multiadm
  #' 
  #' @param x a multiadm object
  #' @param h heights at which the sed rate is evaluated
  #' 
  #' @returns a vector of the same length as h, containing the sedimentation rates in the strat domain at the heights h
  adm_list = split_multiadm(x)
  sed_rate_list = list()
  for (i in seq_along(adm_list)){
    sed_rate_list[[i]] = sed_rate_l(adm_list[[i]], h)
  }
  l = list()
  for (j in seq_along(h)){
    l[[j]] = sapply(sed_rate_list, function(x) x[j])
  }
  sedr = sapply(l, median)
  return(sedr)
}
# convert into cm/kyr
sedr_const =   100 * median_sed_rate_l(my_adm_const, h)
sedr_inc =  100 * median_sed_rate_l(my_adm_inc, h)
sedr_dec =  100 * median_sed_rate_l(my_adm_dec ,h)

sedr_stats = list("sedr_range_const"= range(sedr_const),
                  "sedr_range_inc" = range(sedr_inc),
                  "sedr_range_dec" = range(sedr_dec),
                  "sedr_fac_const" = max(sedr_const)/ min(sedr_const),
                  "sedr_fac_ind" = max(sedr_inc)/ min(sedr_inc),
                  "sedr_fac_dec" = max(sedr_dec)/ min(sedr_dec),
                  "cond_range_const" = range(1/sedr_const),
                  "cond_range_inc" = range(1/sedr_inc),
                  "cond_range_dec" = range(1/sedr_dec))

#### sedimentation rate plot ####
sed_rate_plot = function(file_name){
  df = data.frame(h = rep(-h, 3), sedr = c(sedr_const, sedr_dec, sedr_inc),
                  Scenario = c(rep("Constant flux", length(sedr_const)),rep("Decreasing flux", length(sedr_inc)), rep("Increasing flux", length(sedr_dec))))
  m_sedr = max(df$sedr)
  
  rect = data.frame(h_min = -c(pre_interval[1], main_interval[1], recovery_interval[1]),
                    h_max = -c(pre_interval[2], main_interval[2], recovery_interval[2]),
                    ymax = rep(m_sedr, 3),
                    ymin = rep(0, 3))
  plt = ggplot(df, aes(x = h, y = sedr, group = Scenario, col = Scenario)) +
    ylim(c(0, m_sedr)) +
    geom_rect(data = rect, inherit.aes = FALSE, aes(xmin = h_min, xmax = h_max, ymin = ymin, ymax = ymax), fill = box_cols) +
    scale_color_manual(values = scenario_cols) +
    geom_line(size = 1, alpha = 0.7) +
    xlab("Meters below sea floor [m]") +
    ylab("Sedimentation Rate [cm/kyr]") +
    scale_x_reverse() +
    annotate("text", x = - mean(main_interval), y = 10, label = "main interval", angle = 90)  +
    annotate("text", x =  168, y = 3, label = "recovery interval", angle = 90) +
    annotate("text", x = - mean(pre_interval), y = 10, label = "precursor interval", angle = 90) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.1, 0.9)) +
    ggtitle("Sedimentation rate")
  ggsave(paste0("figs/",file_name, ".png"), plot = plt)
  return(plt)
}

sedr_plot = sed_rate_plot("site690_sedrate")

#### Condensation plot ####

condensation_plot = function(file_name){
  df = data.frame(h = rep(-h, 3), cond = c(1/sedr_const, 1/sedr_dec, 1/sedr_inc),
                  Scenario = c(rep("Constant flux", length(sedr_const)),rep("Decreasing flux", length(sedr_inc)), rep("Increasing flux", length(sedr_dec))))
  m_cond = max(df$cond)
  
  rect = data.frame(h_min = -c(pre_interval[1], main_interval[1], recovery_interval[1]),
                    h_max = -c(pre_interval[2], main_interval[2], recovery_interval[2]),
                    ymax = rep(m_cond, 3),
                    ymin = rep(0, 3))
  plt = ggplot(df, aes(x = h, y = cond, group = Scenario, col = Scenario)) +
    ylim(c(0, m_cond)) +
    geom_rect(data = rect, inherit.aes = FALSE, aes(xmin = h_min, xmax = h_max, ymin = ymin, ymax = ymax), fill = box_cols) +
    scale_color_manual(values = scenario_cols) +
    geom_line(size = 1, alpha = 0.7) +
    xlab("Meters below sea floor [m]") +
    ylab("Condensation [kyr/cm]") +
    scale_x_reverse() +
    annotate("text", x = - mean(main_interval), y = 0.2, label = "main interval", angle = 90)  +
    annotate("text", x =  -mean(recovery_interval), y = 1, label = "recovery interval", angle = 90) +
    annotate("text", x = - mean(pre_interval), y = 0.2, label = "precursor interval", angle = 90) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.8, 0.9)) +
    ggtitle("Condensation")
  ggsave(paste0("figs/",file_name, ".png"), plot = plt)
  return(plt)
}

cond_plot = condensation_plot("site690_condensation")

#### JOin plot of sedimentation and condensation

plt = egg::ggarrange(sedr_plot, cond_plot, nrow = 1, ncol = 2, labels = LETTERS[1:2])

ggsave("figs/site690_join_sedrate_cond.png", plot = plt)

#### Duration PETM main interval ####
petm_main_duration_plot = function(file_name){
  h_int = main_interval
  
  diff_const = sapply(admtools::get_time(my_adm_const, h_int), function(x) diff(x))
  diff_inc = sapply(admtools::get_time(my_adm_inc, h_int), function(x) diff(x))
  diff_dec = sapply(admtools::get_time(my_adm_dec, h_int), function(x) diff(x))
  
  df = data.frame(dur = -1 * c(diff_const, diff_inc, diff_dec),
                  Scenario = c(rep("Constant flux", length(diff_const)), rep("Increasing flux", length(diff_inc)), rep("Decreasing flux", length(diff_dec))) )
  plt = ggplot(df, aes(x = dur, fill = Scenario)) +
    geom_density(alpha = 0.5) +
    xlab("Duration [kyr]") +
    ylab("Density") +
    ggtitle("Duration of PETM main interval") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.9, 0.9))
  ggsave(paste0("figs/", file_name ,".png"), plot = plt)
  return(plt)
}
petm_dur_plot1 = petm_main_duration_plot("site690_PETM_main_interval_duration")

petm_main_stats = function(){
  h_int = main_interval
  diff_const = sapply(admtools::get_time(my_adm_const, h_int), function(x) diff(x))
  diff_inc = sapply(admtools::get_time(my_adm_inc, h_int), function(x) diff(x))
  diff_dec = sapply(admtools::get_time(my_adm_dec, h_int), function(x) diff(x))
  return(list("med_duration_const" = median(-diff_const),
              "med_duration_inc" = median(-diff_inc),
              "med_duration_dec" = median(-diff_dec),
              "iqr_duration_const" = IQR(-diff_const),
              "iqr_duration_inc" = IQR(-diff_inc),
              "iqr_duration_dec" = IQR(-diff_dec)
  ))
}
petm_main_stats = petm_main_stats()

#### Duration of PETM recovery interval ####

petm_recovery_duration_plot = function(file_name){
  
  h_int = recovery_interval
  
  diff_const = sapply(admtools::get_time(my_adm_const, h_int), function(x) diff(x))
  diff_inc = sapply(admtools::get_time(my_adm_inc, h_int), function(x) diff(x))
  diff_dec = sapply(admtools::get_time(my_adm_dec, h_int), function(x) diff(x))
  
  df = data.frame(dur = -1 * c(diff_const, diff_inc, diff_dec),
                  Scenario = c(rep("Constant flux", length(diff_const)), rep("Increasing flux", length(diff_inc)), rep("Decreasing flux", length(diff_dec))) )
  plt = ggplot(df, aes(x = dur, fill = Scenario)) +
    geom_density(alpha = 0.5) +
    xlab("Duration [kyr]") +
    ylab("Density") +
    ggtitle("Duration of PETM recovery interval") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.9, 0.9))
  ggsave(paste0("figs/", file_name ,".png"), plot = plt)
  return(plt)
}
petm_dur_plot2 = petm_recovery_duration_plot("site690_PETM_recovery_interval_duration")

petm_recovery_stats = function(){
  h_int = recovery_interval
  
  diff_const = sapply(admtools::get_time(my_adm_const, h_int), function(x) diff(x))
  diff_inc = sapply(admtools::get_time(my_adm_inc, h_int), function(x) diff(x))
  diff_dec = sapply(admtools::get_time(my_adm_dec, h_int), function(x) diff(x))
  return(list("med_duration_const" = median(-diff_const),
              "med_duration_inc" = median(-diff_inc),
              "med_duration_dec" = median(-diff_dec),
              "iqr_duration_const" = IQR(-diff_const),
              "iqr_duration_inc" = IQR(-diff_inc),
              "iqr_duration_dec" = IQR(-diff_dec)
  ))
}
petm_recovery_stats = petm_recovery_stats()


#### join plot of PETM duration

plt = egg::ggarrange(petm_dur_plot1, petm_dur_plot2, nrow = 1, ncol = 2, labels = LETTERS[1:2])

ggsave("figs/site690_joint_duration.png", plot = plt)
