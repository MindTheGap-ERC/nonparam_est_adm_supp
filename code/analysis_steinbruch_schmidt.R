#### Load packages ####
library(ggplot2)
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

save.image(file = "data/res/sbs.RData")


#### Plot ADMs ####
dpi = 400
lab_size = 7
title_size = 8
fig_width_cm = 12
annot_size = 5
legend_size = 3
ax_size = 4

box_col = grey(0.7)
col_med = "red"
col_env = "blue"
env_lwd = 0.5
med_lwd = 1
unc_lwd = 0.4


plot_sbs_adm_abs = function(adm, xlab){
  #adm = adm_prec_abs
  q2_adm = admtools::quantile_adm(adm, h, 0.975)
  q1_adm = admtools::quantile_adm(adm, h,  0.025)
  m_adm = admtools::quantile_adm(adm, h, 0.5)
  t_max = -  max(c(q1_adm$t/1000, q2_adm$t/1000, m_adm$t/1000))
  t_min =  - min(c(q1_adm$t/1000, q2_adm$t/1000, m_adm$t/1000))
  h_min = min(c(q1_adm$h, q2_adm$h, m_adm$h))
  
  df = data.frame(he =rep(h, 5),
                  t = -c(q2_adm$t/1000, q1_adm$t/1000, m_adm$t/1000, m_adm$t/1000-uncertainty_radiometric, m_adm$t/1000+uncertainty_radiometric),
                  type = c(rep("95 % Envelope", 2 * length(h)), rep("Median", length(h)), rep("Age uncertainty (2σ)", 2 * length(h))),
                  group = rep(LETTERS[1:5], each = length(h)))
  
  
  
  rect = data.frame(h_min = c(h_bottom_ukw, h_min), h_max = c(h_top_ukw, h_top_lkw),
                    t_min = rep(t_min, 2), t_max = rep(t_max, 2))
  
  plt = ggplot(df, aes(y = he, x = t, color = type, group = group, linetype = type)) + 
    geom_rect(rect, inherit.aes = FALSE, mapping = aes(xmin = t_min, xmax = t_max, ymin = h_min, ymax = h_max), fill = box_col) +
    geom_line(aes(size = group)) +
    scale_size_manual(values = c("A" = env_lwd, "B" = env_lwd, "C" = med_lwd, "D" = unc_lwd, "E" = unc_lwd), guide = "none") +
    scale_color_manual(values = c(col_env, col_med, col_med)) +
    scale_linetype_manual(values = c("solid", "dashed", "solid")) +
    xlab(xlab) +
    ylab("Stratigraphic position [m]") +
    annotate("text", x = mean(c(t_max, t_min)), y = mean(c(h_top_lkw, h_min)), label = "Lower Kellwasser Bed", size = annot_size/.pt) +
    annotate("text", x = mean(c(t_max, t_min)) - 0.2, y = mean(c(h_top_ukw, h_bottom_ukw)) - 0.1, label = "Upper Kellwasser Bed", size = annot_size/.pt) +
    geom_hline(yintercept = h_ashbed) +
    annotate("text", x = mean(c(t_max, t_min)), y = h_ashbed + 0.05, label = "Bentonite layer", col = grey(0.4) , size = annot_size/.pt) +
    theme(legend.title = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.1, 0.89),
          plot.title = element_text(size = title_size),
          axis.title = element_text(size = lab_size),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = legend_size),
          axis.text = element_text(size = ax_size)) +
    scale_x_reverse()

  return(plt)
}


## plot of anchored age-depth model without uncertainty of radiometric dates
plot_sbs_adm_float_abs = function(adm, xlab){
  q2_adm = admtools::quantile_adm(adm, h, 0.975)
  q1_adm = admtools::quantile_adm(adm, h,  0.025)
  m_adm = admtools::quantile_adm(adm, h, 0.5)
  t_max = -min(c(q1_adm$t/1000, q2_adm$t/1000, m_adm$t/1000)) 
  t_min = -max(c(q1_adm$t/1000, q2_adm$t/1000, m_adm$t/1000)) 
  h_min = min(c(q1_adm$h, q2_adm$h, m_adm$h))
  
  df = data.frame(he =rep(h, 3),
                  t = -c(q2_adm$t , q1_adm$t , m_adm$t )/1000,
                  type = c(rep("95 % Envelope", 2 * length(h)), rep("Median", length(h))),
                  group = rep(LETTERS[1:3], each = length(h)))
  
  
  rect = data.frame(h_min = c(h_bottom_ukw, h_min), h_max = c(h_top_ukw, h_top_lkw),
                    t_min = rep(t_min -0.1, 2), t_max = rep(t_max + 0.1, 2))
  
  plt = ggplot(df, aes(y = he, x = t, color = type, group = group)) + 
    geom_rect(rect, inherit.aes = FALSE, mapping = aes(xmin = t_min, xmax = t_max, ymin = h_min, ymax = h_max), fill = box_col) +
    geom_line(aes(size = group)) +
    scale_size_manual(values = c("A" = env_lwd, "B" = env_lwd, "C" = med_lwd), guide = "none") +
    scale_color_manual(values = c(col_env, col_med)) +
    xlab(xlab) +
    ylab("Stratigraphic position [m]")  +
    annotate("text", x = mean(c(t_max, t_min)), y = mean(c(h_top_lkw, h_min)), label = "Lower Kellwasser Bed", size = annot_size/.pt) +
    annotate("text", x = mean(c(t_max, t_min) - 0.2), y = mean(c(h_top_ukw, h_bottom_ukw)) - 0.1, label = "Upper Kellwasser Bed", size = annot_size/.pt) +
    geom_hline(yintercept = h_ashbed) +
    annotate("text", x = mean(c(t_max, t_min)), y = h_ashbed + 0.05, label = "Bentonite layer", col = grey(0.4) , size = annot_size/.pt) +
    theme(legend.title = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.1, 0.92),
          plot.title = element_text(size = title_size),
          axis.title = element_text(size = lab_size),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = legend_size),
          axis.text = element_text(size = ax_size)) +
    scale_x_reverse()
  plt
  return(plt)
}

#### JOint plot of adms ####
## precession
adm_anchor = plot_sbs_adm_abs(adm_prec_abs, "Age [Ma]")
adm_anchor_no_rad = plot_sbs_adm_float_abs(adm_prec_abs_no_error, "Age [Ma]")

# plot 2 age-depth models
plt = egg::ggarrange( adm_anchor_no_rad, adm_anchor, 
                      nrow = 1,
                      ncol = 2,
                      labels = LETTERS[1:2])
ggsave("figs/sbs_join_adm_prec.png",
       plot = plt, 
       width = fig_width_cm,
       height = 8,
       unit = "cm",
       dpi = dpi)

## short eccentricity
adm_anchor = plot_sbs_adm_abs(adm_secc_abs, "Age [Ma]")
adm_anchor_no_rad = plot_sbs_adm_float_abs(adm_secc_abs_no_error, "Age [Ma]")

plt = egg::ggarrange( adm_anchor_no_rad, adm_anchor,
                      nrow = 1, 
                      ncol = 2,
                      labels = LETTERS[1:2])
ggsave("figs/sbs_join_adm_secc.png",
       plot = plt,
       width = fig_width_cm,
       height = 8,
       unit = "cm",
       dpi = dpi)

#### filled contour plots for eTimeOpt ####
#### values of r^2_opt  ###
plot_r2_opt = function(output){
  v1 =  output$results[,which(output$heights == 1.56)]
  v2 =  output$results[,which(output$heights == 3.54)]
  s = output$sed_rate
  df = data.frame(s = rep(s, 2), val = c(v1, v2), Height = c(rep("1.56 m", length(v1)), rep("3.54 m", length(v2))))
  #df = data.frame(s, v1, v2)
  plt = ggplot(df, aes(x =s, y = val, group = Height, color = Height)) +
    geom_line(lwd = 1) +
    ggtitle(expression(r[opt]^2 ~ at ~ height)) +
    xlab("Sedimentation rate [cm/kyr]") +
    ylab(expression(r[opt]^2)) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.2, 0.9),
          plot.title = element_text(size = title_size),
          axis.title = element_text(size = lab_size),
          axis.text = element_text(size = ax_size),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = lab_size))
  return(plt)
}

## precession
r2_plot_prec = plot_r2_opt(fa_prec)
a = fa_prec$results
dimnames(a) = list("sedr" = fa_prec$sed_rate,
                   "height" = fa_prec$heights)
aa = reshape2::melt(a)

eTimeOpt_plot_prec = ggplot(aa, aes(x = sedr, y = height,  z = value)) +
  geom_contour_filled() +
  xlab("Sedimentation Rate [cm/kyr]") +
  ylab("Height [m]") +
  ggtitle(expression(r[opt]^2 ~ at ~ Steinbruch ~ Schmidt)) +
  theme( plot.title = element_text(size = title_size),
         axis.title = element_text(size = lab_size),
         axis.text = element_text(size = ax_size),
         legend.key.size = unit(0.25, "cm"),
         legend.text = element_text(size = 5),
         legend.title = element_blank())

plt = egg::ggarrange(r2_plot_prec, eTimeOpt_plot_prec,
                     ncol = 2,
                     nrow = 1, 
                     labels = LETTERS[1:2])
ggsave(filename = "figs/r2_plot_prec_join.png", 
       plot = plt, 
       width = fig_width_cm,
       height = 8,
       unit = "cm", 
       dpi = dpi)

### short eccentricity
r2_plot_secc = plot_r2_opt(fa_secc)
a = fa_secc$results
dimnames(a) = list("sedr" = fa_secc$sed_rate, 
                   "height" = fa_secc$heights)
aa = reshape2::melt(a)

eTimeOpt_plot_secc = ggplot(aa, aes(x = sedr, y = height,  z = value)) +
  geom_contour_filled() +
  xlab("Sedimentation Rate [cm/kyr]") +
  ylab("Height [m]") +
  ggtitle(expression(r[opt]^2 ~ at ~ height ~ and ~ sedimentation ~ rate)) +
  ggtitle(expression(r[opt]^2 ~ at ~ Steinbruch ~ Schmidt)) +
  theme( plot.title = element_text(size = title_size),
         axis.title = element_text(size = lab_size),
         axis.text = element_text(size = ax_size),
         legend.key.size = unit(0.25, "cm"),
         legend.text = element_text(size = 5),
         legend.title = element_blank())

plt = egg::ggarrange(r2_plot_secc, eTimeOpt_plot_secc,
                     ncol = 2, 
                     nrow = 1, 
                     labels = LETTERS[1:2])
ggsave(filename = "figs/r2_plot_secc_join.png", 
       plot = plt,
       width = fig_width_cm,
       height = 8, 
       unit = "cm", 
       dpi = dpi)


#### Density plot of durations ####
#### absolute age of ff boundary 

plot_age_ff_boundary = function(){
  l_secc = sapply(admtools::get_time(x = adm_secc_abs, h = h_f_f_bdry ), function(x) x[1])
  l_prec = sapply(admtools::get_time(x = adm_prec_abs, h = h_f_f_bdry ), function(x) x[1])
  df = data.frame(dur = -1 * c(l_prec, l_secc) / 1000, Test = c(rep("Precession", length(l_prec)), rep("Short eccentricity", length(l_secc))))

  plt = ggplot(df, aes(x = dur, fill = Test)) + 
    geom_density(alpha = 0.5) + 
    xlab("Age [Ma]") +
    ylab("Density") +
    ggtitle(paste0("Age of F-F boundary")) +
    scale_x_reverse() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.79, 0.95),
          legend.title = element_blank(),
          plot.title = element_text(size = title_size),
          axis.title = element_text(size = lab_size),
          axis.text = element_text(size = ax_size),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5))
  return(plt)
}

#### Duration of upper Kellwasser event 

plot_duration_uke = function(file_name, file_type){
  l_prec = sapply(admtools::get_time(x = adm_prec_abs, h = c(h_top_ukw, h_bottom_ukw )),  function(x) x[1] - x[2])
  l_secc = sapply(admtools::get_time(x = adm_secc_abs, h = c(h_top_ukw, h_bottom_ukw )),  function(x) x[1] - x[2])
  df = data.frame(dur = c(l_prec, l_secc), Test = c(rep("Precession", length(l_prec)), rep("Short eccentricity", length(l_secc))))
  
  plt = ggplot(data = df, aes(x = dur, fill = Test) ) + 
    geom_density(alpha = 0.5) + 
    xlab("Time [kyr]") + 
    ylab("Density") +
    ggtitle("Duration UKE") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.7, 0.95),
          legend.title = element_blank(),
          plot.title = element_text(size = title_size),
          axis.title = element_text(size = lab_size),
          axis.text = element_text(size = ax_size),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5))
  return(plt)
}

#### time elapsed between bentonite layer and the ff boundary 

plot_elapsed_time_ke = function(){
  l_prec = sapply(admtools::get_time(x = adm_prec_abs, h = c(h_bottom_ukw, h_top_lkw )),  function(x) x[1] - x[2])
  l_secc = sapply(admtools::get_time(x = adm_secc_abs, h = c(h_bottom_ukw, h_top_lkw )),  function(x) x[1] - x[2])
  df = data.frame(dur = c(l_prec, l_secc), Test = c(rep("Precession", length(l_prec)), rep("Short eccentricity", length(l_secc))))
  
  plt = ggplot(data = df, aes(x = dur, fill = Test) ) + 
    geom_density(alpha = 0.6) + 
    xlab("Time [kyr]") + 
    ylab("Density") +
    ggtitle("Time between KEs") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.7, 0.95),
          legend.title = element_blank(),
          plot.title = element_text(size = title_size),
          axis.title = element_text(size = lab_size),
          axis.text = element_text(size = ax_size),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5))
  return(plt)
  
}

#### Join plot of durations sbs ####
age_ff_plot = plot_age_ff_boundary()
duration_uke_plot = plot_duration_uke()
elapsed_time_plot = plot_elapsed_time_ke()
plt = egg::ggarrange(age_ff_plot, duration_uke_plot, elapsed_time_plot, 
                     ncol = 3,
                     nrow = 1,
                     labels = LETTERS[1:3])
ggsave("figs/sbs_durations_joint.png", 
       plot = plt, 
       width = fig_width_cm,
       height = 8,
       unit = "cm",
       dpi = dpi)

#### extract statistics ####
# get the numbers

ci_eqiv = function(x){
  #' get 95 % HD interval
  y = quantile(x, probs = c(0.025, 0.975))
  return(c(unname(diff(y))))
}

sbs_stats = list(
  "mean_dur_bentonite_ff" = mean(sapply(admtools::get_time(adm_prec_abs_no_error, h =  c(h_ashbed, h_f_f_bdry)), function(x) abs(x[1]- x[2]))),
  "95%_HDI_dur_bentonite_ff" = sapply(admtools::get_time(adm_prec_abs_no_error, h =  c(h_ashbed, h_f_f_bdry)), function(x) abs(x[1]- x[2])) |> ci_eqiv(),
  "sd_dur_bentonite_ff" = sapply(admtools::get_time(adm_prec_abs_no_error, h =  c(h_ashbed, h_f_f_bdry)), function(x) abs(x[1]- x[2])) |> sd(),
  "median_age_ff_bdry" = sapply(admtools::get_time(adm_prec_abs, h =  c(h_f_f_bdry)), function(x) abs(x[1]/1000)) |> median(),
  "HDI_age_ff_bdry" = sapply(admtools::get_time(adm_prec_abs, h =  c(h_f_f_bdry)), function(x) abs(x[1]/1000)) |> ci_eqiv(),
  "mean_age_ff_bdry" = sapply(admtools::get_time(adm_prec_abs, h =  c(h_f_f_bdry)), function(x) abs(x[1]/1000)) |> mean(),
  "twosigms_age_ff_bdry" = sapply(admtools::get_time(adm_prec_abs, h =  c(h_f_f_bdry)), function(x) abs(x[1]/1000)) |> sd() |> (\(x){x * 2})()
)

sbs_stats[["diff_this_study_daSilva"]] = sbs_stats$mean_age_ff_bdry - ff_mean_dasilva
sbs_stats[["quot_uncert_this_study_daSilva"]] = sbs_stats$twosigms_age_ff_bdry / ff_2sigma_dasilva
sbs_stats[["ff_mean_diff_da_silva"]] = ff_mean_gradstein - ff_mean_dasilva
sbs_stats[["ff_mean_diff_gradstein"]] = ff_mean_gradstein - sbs_stats$mean_age_ff_bdry
dur_uke = sapply(admtools::get_time(x = adm_prec_abs, h = c(h_top_ukw, h_bottom_ukw )),  function(x) x[1] - x[2])
sbs_stats[["duration_uke"]] = c("median" = dur_uke |> median(), 
                                dur_uke |> quantile(),
                                "IQR" = dur_uke |> IQR())
time_between_ke = sapply(admtools::get_time(x = adm_prec_abs, h = c(h_bottom_ukw, h_top_lkw )),  function(x) x[1] - x[2])
sbs_stats[["time_between_ke"]] = c("median" = time_between_ke |> median(),
                                   time_between_ke |> quantile(),
                                   "IQR" = time_between_ke |> IQR())
dur_lke_wichern = 250 # estimated duration of lower kellwasser event from Wichern et al. (2024)
sbs_stats[["time_onset_ke"]] = unname(sbs_stats$time_between_ke["median"]) + dur_lke_wichern 
sbs_stats[["duration_ke_crisis"]] = unname(sbs_stats$time_between_ke["median"]) + dur_lke_wichern + sbs_stats$duration_uke["median"] |> unname()