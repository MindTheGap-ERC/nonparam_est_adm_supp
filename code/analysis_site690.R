#### install latest version of package ####
if (FALSE){
remotes::install_github(repo = "MindTheGap-ERC/admtools",
                        build_vignettes = TRUE,
                        ref = "dev",
                        force = TRUE)

}
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
no_of_rep = 1000
subdivisions = 10000

#### Load data & clean data ####
# from https://doi.org/10.1594/PANGAEA.723912 , supplement to Farley and Eltgroth 2003
data = read.csv(file = "data/Farley_and_Eltgroth_2003_supp_data_1_site690.csv", header = TRUE, sep = ";")

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

#### Duration PETM ####
# main interval
dur_main_const =  - sapply(get_time(my_adm_const, h = main_interval), diff)
dur_main_dec =  - sapply(get_time(my_adm_dec, main_interval), diff)
dur_main_inc =  - sapply(get_time(my_adm_inc, main_interval), diff)

quantile(dur_main_inc)
quantile(dur_main_const)
quantile(dur_main_dec)

IQR(dur_main_inc)
IQR(dur_main_const)
IQR(dur_main_dec)
# recovery
dur_rec_const =  - sapply(get_time(my_adm_const, recovery_interval), diff)
dur_rec_inc = -  sapply(get_time(my_adm_inc, recovery_interval), diff)
dur_rec_dec =  - sapply(get_time(my_adm_dec, recovery_interval), diff)

quantile(dur_rec_inc)
quantile(dur_rec_const)
quantile(dur_rec_dec)

IQR(dur_rec_inc)
IQR(dur_rec_const)
IQR(dur_rec_inc)
# pre interval
dur_pre_const = -sapply(get_time(my_adm_const, pre_interval), diff)
dur_pre_dec = -sapply(get_time(my_adm_dec, pre_interval), diff)
dur_pre_inc = -sapply(get_time(my_adm_inc, pre_interval), diff)

# ratio main interval to recovery
thickness_ratio = diff(recovery_interval)/ diff(main_interval)
time_ratio_const = dur_rec_const / dur_main_const
time_ratio_inc = dur_rec_inc / dur_main_inc
time_ratio_dec = dur_rec_dec / dur_main_dec

med_ratio_const = 100 * (1 - median(time_ratio_const))
med_ratio_inc = 100 * (1 - median(time_ratio_inc))
med_ratio_dec = 100 * (1 - median(time_ratio_dec))

dur_total_const =  sapply(get_time(my_adm_const, range(c(recovery_interval, main_interval, pre_interval))), diff)

#### Plotting ####
med_lty = 1
med_lwd = 3
env_lty = 4
env_lwd = 1

col_const = "black"
col_dec = "blue"
col_inc = "red"



#### make adm plot ####
adm_plot = function(){
const_med = quantile_adm(my_adm_const, h, 0.5)
const_p1 = quantile_adm(my_adm_const, h, 0.975)
const_p2 = quantile_adm(my_adm_const, h, 0.025)

inc_med = quantile_adm(my_adm_inc, h, 0.5)
inc_p1 = quantile_adm(my_adm_inc, h, 0.975)
inc_p2 = quantile_adm(my_adm_inc, h, 0.025)

dec_med = quantile_adm(my_adm_dec, h, 0.5)
dec_p1 = quantile_adm(my_adm_dec, h, 0.975)
dec_p2 = quantile_adm(my_adm_dec, h, 0.025)

plot(NULL,
     xlim = c(-200, 200),
     ylim = range(h),
     xlab = "time since beginning of PETM main event [kyr]",
     ylab = "stratigraphic position [m]")


rect(xleft = -1000, xright = 1000, ytop = pre_interval[1], ybottom = pre_interval[2], col = "azure2", border = NA)
rect(xleft = -1000, xright = 1000, ytop = main_interval[1], ybottom = main_interval[2], col = "azure3", border = NA)
rect(xleft = -1000, xright = 1000, ytop = recovery_interval[1], ybottom = recovery_interval[2], col = "azure4", border = NA)

lines(const_med$t, const_med$h, col = col_const, lwd = med_lwd, lty = med_lty)
lines(const_p1$t, const_p1$h, col = col_const, lwd = env_lwd, lty = env_lty)
lines(const_p2$t, const_p2$h, col = col_const, lwd = env_lwd, lty = env_lty)

lines(dec_med$t, dec_med$h, col = col_dec, lwd = med_lwd, lty = med_lty)
lines(dec_p1$t, dec_p1$h, col = col_dec, lwd = env_lwd, lty = env_lty)
lines(dec_p2$t, dec_p2$h, col = col_dec, lwd = env_lwd, lty = env_lty)

lines(inc_med$t, inc_med$h, col = col_inc, lwd = med_lwd, lty = med_lty)
lines(inc_p1$t, inc_p1$h, col = col_inc, lwd = env_lwd, lty = env_lty)
lines(inc_p2$t, inc_p2$h, col = col_inc, lwd = env_lwd, lty = env_lty)

points(0, min(main_interval), cex = 1.5, pch = 16)

legend("topleft",
       col = c(col_const, col_inc, col_dec),
       legend = c("Constant flux", "Increasing flux", "Decreasing flux"),
       lwd = med_lwd,
       lty = med_lty)
xpos = -200
text(x = xpos, y = mean(recovery_interval), "recovery interval", pos = 4)
text(x = xpos, y = mean(main_interval), "main interval", pos = 4)
text(x = xpos, y = mean(pre_interval), "pre interval", pos = 4)
}
png("figs/site690_adm.png")
adm_plot()
dev.off()
#### Sed rate plot ####

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

sed_rate_plot = function(){
plot(NULL, 
     xlim = range(h),
     ylim = c(0, 1.1 * max(c(sedr_const, sedr_dec, sedr_inc))),
     xlab = "Height [m]",
     ylab = "Sedimentation rate [cm/kyr]")
rect(ybottom = 0, ytop = 300, xleft = pre_interval[1], xright = pre_interval[2], col = "azure2")
rect(ybottom = 0, ytop = 300, xleft = main_interval[1], xright = main_interval[2], col = "azure3")
rect(ybottom = 0, ytop = 300, xleft = recovery_interval[1], xright = recovery_interval[2], col = "azure4")
lines(h, sedr_const, col = col_const, lwd = med_lwd)
lines(h, sedr_dec, col = col_dec, lwd = med_lwd)
lines(h, sedr_inc, col = col_inc, lwd = med_lwd)

legend("topleft",
       lwd = med_lwd,
       lty = med_lty,
       col = c(col_const, col_dec, col_inc),
       legend = c("constant flux", "decreasing flux", "increasing flux"))
}

png("figs/site690_sedrate.png")
sed_rate_plot()
dev.off()
#### codensation plot ####
make_cond_plot = function(){
  plot(NULL,
       xlim = range(h),
       ylim = range(c(0,(1.1 * max(c(1/sedr_const, 1/sedr_dec, 1/sedr_inc))))),
       xlab = "Height [m]",
       ylab = "Condensation [kyr/cm]")
  rect(ybottom = 0, ytop = 300, xleft = pre_interval[1], xright = pre_interval[2], col = "azure2")
  rect(ybottom = 0, ytop = 300, xleft = main_interval[1], xright = main_interval[2], col = "azure3")
  rect(ybottom = 0, ytop = 300, xleft = recovery_interval[1], xright = recovery_interval[2], col = "azure4")
  lines(h, 1/sedr_const, col = col_const, lwd = med_lwd)
  lines(h, 1/sedr_dec, col = col_dec, lwd = med_lwd)
  lines(h, 1/sedr_inc, col = col_inc, lwd = med_lwd)
  text(x = mean(pre_interval), y = 0.2, labels = "pre interval", srt = 90, pos = 3)
  text(x = mean(main_interval), y = 0.2, labels = "main interval", srt = 90, pos = 3)
  text(x = mean(recovery_interval),  y = 1, labels = "recovery interval", srt= 90, pos = 3)
  legend("topleft",
         lwd = med_lwd,
         lty = med_lty,
         col = c(col_const, col_dec, col_inc),
         legend = c("constant flux", "decreasing flux", "increasing flux"))
}

png("figs/site690_condensation.png")
make_cond_plot()
dev.off()

max_sedr = max(sedr_const, sedr_dec, sedr_inc)
min_sedr = min(sedr_const, sedr_dec, sedr_inc)
sedr_fact = max_sedr / min_sedr

sedr_range_const = range(sedr_const)
sedr_range_inc = range(sedr_inc)
sedr_range_dec = range(sedr_dec)

sedr_range_const[2]/ sedr_range_const[1]
sedr_range_dec[2]/ sedr_range_dec[1]
sedr_range_inc[2]/ sedr_range_inc[1]

range(1/sedr_inc)

#### Duration PETM ####
library(ggplot2)

h_int = unique(c(recovery_interval, main_interval)) # stratigraphic boundaries of main and recovery interval

l = get_time(my_adm_const, h_int)
diff_const = lapply(l, function(x) diff(x))
rat_const = sapply(diff_const, function(x) x[1]/x[2])

l = get_time(my_adm_inc, h_int)
diff_inc = lapply(l, function(x) diff(x))
rat_inc = sapply(diff_inc, function(x) x[1]/x[2])

l = get_time(my_adm_dec, h_int)
diff_dec = lapply(l, function(x) diff(x))
rat_dec = sapply(diff_dec, function(x) x[1]/x[2])

dur_const = -sapply(diff_const, function(x) x[2])
dur_inc = -sapply(diff_inc, function(x) x[2])
dur_dec = -sapply(diff_dec, function(x) x[2])

data = data.frame(type = c(rep("const", length(dur_const)), rep("dec", length(dur_dec)), rep("inc", length(dur_inc))),
                  value = c(dur_const, dur_dec, dur_inc))

ggplot(data, aes(x = value, fill = type )) +
  geom_density() + 
  xlab("Duration [kyr]") +
  ylab("Density") +
  ggtitle("Duration of PETM main event") 

quantile(dur_const)
quantile(dur_dec)
quantile(dur_inc)

median(dur_inc) / median(dur_const)
median(dur_dec) / median(dur_const)


quantile(rat_inc)
quantile(rat_dec)
quantile(rat_const)
diff(recovery_interval)/diff(main_interval)