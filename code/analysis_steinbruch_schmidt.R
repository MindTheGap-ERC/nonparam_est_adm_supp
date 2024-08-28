#### install latest version ####
if(FALSE){
  remotes::install_github(repo = "MindTheGap-ERC/admtools",
                          build_vignettes = TRUE,
                          ref = "HEAD",
                          dependencies = TRUE)
  
}


#### Load packages ####
library(ggplot2)

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
t_ashbed_sd = 0.053 * 1000 / 2 # percival uses 2 sigma

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
etimeOptMS_prec=astrochron::eTimeOpt(MSi,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6, numsed=100,linLog=1,
                     limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                     detrend=T,ydir=1,output=1,genplot=F,check=T,verbose=1)
# test for short eccentricty amplitude modulation
etimeOptMS_secc=astrochron::eTimeOpt(MSi,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6, numsed=100,linLog=1,
                     limit=T,fit=2,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                     detrend=T,ydir=1,output=1,genplot=F,check=T,verbose=1)

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

# floating age tie point
t_tp_float = function(){
  return(0)
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
adm_prec_abs = admtools::sedrate_to_multiadm(h_tp , t_tp_absolute, sed_rate_gen = se_prec, h, no_of_rep = no_of_rep )
adm_prec_float = admtools::sedrate_to_multiadm(h_tp , t_tp_float, sed_rate_gen = se_prec, h, no_of_rep = no_of_rep )

adm_secc_abs = admtools::sedrate_to_multiadm(h_tp , t_tp_absolute, sed_rate_gen = se_secc, h = h, no_of_rep = no_of_rep )
adm_secc_float = admtools::sedrate_to_multiadm(h_tp , t_tp_float, sed_rate_gen = se_secc, h = h, no_of_rep = no_of_rep )

#### Plot ADMs ####
box_col = grey(0.7)
col_med = "red"
col_env = "blue"
env_lwd = 1
med_lwd = 2

plot_sbs_adm_float = function(adm, xlab, file_name){
  q2_adm = admtools::quantile_adm(adm, h, 0.925)
  q1_adm = admtools::quantile_adm(adm, h,  0.025)
  m_adm = admtools::quantile_adm(adm, h, 0.5)
  t_max = max(c(q1_adm$t, q2_adm$t, m_adm$t))
  t_min = min(c(q1_adm$t, q2_adm$t, m_adm$t)) 
  h_min = min(c(q1_adm$h, q2_adm$h, m_adm$h))
  
  df = data.frame(he =rep(h, 3),
                  t = c(q2_adm$t, q1_adm$t, m_adm$t),
                  type = c(rep("95 % Envelope", 2 * length(h)), rep("Median", length(h))),
                  group = rep(LETTERS[1:3], each = length(h)))
  
  
  
  rect = data.frame(h_min = c(h_bottom_ukw, h_min), h_max = c(h_top_ukw, h_top_lkw),
                    t_min = rep(t_min - 100, 2), t_max = rep(t_max + 100, 2))
  
  ggplot(df, aes(y = he, x = t, color = type, group = group)) + 
    geom_rect(rect, inherit.aes = FALSE, mapping = aes(xmin = t_min, xmax = t_max, ymin = h_min, ymax = h_max), fill = box_col) +
    geom_line(aes(size = group)) +
    scale_size_manual(values = c("A" = env_lwd, "B" = env_lwd, "C" = med_lwd), guide = "none") +
    scale_color_manual(values = c(col_env, col_med)) +
    xlab(xlab) +
    ylab("Stratigraphic position [m]") +
    ggtitle("Age-depth model") +
    annotate("text", x = mean(c(t_max, t_min)), y = mean(c(h_top_lkw, h_min)), label = "Lower Kellwasser Bed") +
    annotate("text", x = mean(c(t_max, t_min)), y = mean(c(h_top_ukw, h_bottom_ukw)), label = "Upper Kellwasser Bed") +
    geom_hline(yintercept = h_ashbed) +
    annotate("text", x = mean(c(t_max, t_min)), y = h_ashbed + 0.05, label = "Bentonite layer", col = grey(0.4) ) +
    theme(legend.title = element_blank())
  ggsave(paste0("figs/", file_name, ".png"))
}

plot_sbs_adm_float(adm_prec_float, "Time [kyr]", "sbs_floating_adm_prec")
plot_sbs_adm_float(adm_secc_float, "Time [kyr]", "sba_floating_adm_secc")

plot_sbs_adm_abs = function(adm, xlab, file_name){
  #adm = adm_prec_abs
  q2_adm = admtools::quantile_adm(adm, h, 0.925)
  q1_adm = admtools::quantile_adm(adm, h,  0.025)
  m_adm = admtools::quantile_adm(adm, h, 0.5)
  t_max = -  max(c(q1_adm$t/1000, q2_adm$t/1000, m_adm$t/1000))
  t_min =  - min(c(q1_adm$t/1000, q2_adm$t/1000, m_adm$t/1000))
  h_min = min(c(q1_adm$h, q2_adm$h, m_adm$h))
  
  df = data.frame(he =rep(h, 3),
                  t = -c(q2_adm$t/1000, q1_adm$t/1000, m_adm$t/1000),
                  type = c(rep("95 % Envelope", 2 * length(h)), rep("Median", length(h))),
                  group = rep(LETTERS[1:3], each = length(h)))
  
  
  
  rect = data.frame(h_min = c(h_bottom_ukw, h_min), h_max = c(h_top_ukw, h_top_lkw),
                    t_min = rep(t_min, 2), t_max = rep(t_max, 2))
  
  ggplot(df, aes(y = he, x = t, color = type, group = group)) + 
    geom_rect(rect, inherit.aes = FALSE, mapping = aes(xmin = t_min, xmax = t_max, ymin = h_min, ymax = h_max), fill = box_col) +
    geom_line(aes(size = group)) +
    scale_size_manual(values = c("A" = env_lwd, "B" = env_lwd, "C" = med_lwd), guide = "none") +
    scale_color_manual(values = c(col_env, col_med)) +
    xlab(xlab) +
    ylab("Stratigraphic position [m]") +
    ggtitle("Age-depth model") +
    annotate("text", x = mean(c(t_max, t_min)), y = mean(c(h_top_lkw, h_min)), label = "Lower Kellwasser Bed") +
    annotate("text", x = mean(c(t_max, t_min)), y = mean(c(h_top_ukw, h_bottom_ukw)), label = "Upper Kellwasser Bed") +
    geom_hline(yintercept = h_ashbed) +
    annotate("text", x = mean(c(t_max, t_min)), y = h_ashbed + 0.05, label = "Bentonite layer", col = grey(0.4) ) +
    theme(legend.title = element_blank()) +
    scale_x_reverse()
  ggsave(paste0("figs/", file_name, ".png"))
}

plot_sbs_adm_abs(adm_prec_abs, "Age [Ma]", "sbs_absolute_adm_prec")
plot_sbs_adm_abs(adm_secc_abs, "Age [Ma]", "sbs_absolute_adm_secc")

# auxiliary function

# plot_sbs = function(adm, file_name){
#   
#   png(file = paste0("figs/", file_name, ".png"))
#   
#   q2_adm = admtools::quantile_adm(adm, h, 0.925)
#   q1_adm = admtools::quantile_adm(adm, h,  0.025)
#   m_adm = admtools::quantile_adm(adm, h, 0.5)
#   
#   t_min = min(range(c(q2_adm$t, q1_adm$t, m_adm$t)))/1000
#   t_max = max(range(c(q2_adm$t, q1_adm$t, m_adm$t)))/1000
#   h_min = min(range(c(q1_adm$h, q2_adm$h, m_adm$h)))
#   h_max = max(range(c(q1_adm$h, q2_adm$h, m_adm$h)))
#   
#   plot(NULL,
#        xlim = c(t_min, t_max),
#        ylim = c(h_min, h_max),
#        xlab = "Time BP [Myr]",
#        ylab = "Stratigraphic Position [m]",
#        main = "ADM for Steinbruch Schmidt (Absolute Age)")
#   col_med = "red"
#   col_env = "blue"
#   lwd_med = 3
#   lwd_env = 3
#   rect(xleft = t_min, xright = t_max, ytop = h_top_lkw, ybottom = h_min - 1, col = grey(0.6), lty = 0)
#   text(x = mean(c(t_min, t_max)), y = mean(c(h_top_lkw, h_top_lkw - 0.1)), labels = "Lower Kellwasser Bed")
#   
#   rect(xleft = t_min, xright = t_max, ytop = h_top_ukw, ybottom = h_bottom_ukw, col = grey(0.6), lty = 0)
#   text(x = mean(c(t_min, t_max)), y = mean(c(h_top_ukw, h_bottom_ukw)), labels = "Upper Kellwasser Bed")
#   lines(m_adm$t/1000, m_adm$h, col = col_med, lwd = lwd_med)
#   lines(q1_adm$t/1000, q1_adm$h, col = col_env, lwd = lwd_env)
#   lines(q2_adm$t/1000, q2_adm$h, col = col_env, lwd = lwd_env)
#   
#   lines(c(t_min, t_max), rep(h_ashbed, 2), lty = 3, lwd = 3, col = grey(0.6))
#   text(x = t_max, y = h_ashbed + 0.1, labels = "Bentonite Layer", pos = 2)
#   
#   legend("topleft", col = c(col_med, col_env), lwd = c(lwd_med, lwd_env), lty = 1, legend = c("Median Age", "95 % Envelope"))
#   dev.off()
# }
# 
# plot_sbs(adm_prec_float, "sds_floating_adm_prec")
# 
# plot_sbs(adm_prec_abs, "sds_anchored_adm_prec")
# 
# plot_sbs(adm_secc_float, "sds_floating_adm_secc")
# 
# plot_sbs(adm_secc_abs, "sds_anchored_adm_secc")

#### values of r^2_opt  ####

plot_r2_opt = function(output, file_name, file_type){
  v1 =  output$results[,which(output$heights == 1.56)]
  v2 =  output$results[,which(output$heights == 3.54)]
  s = output$sed_rate
  df = data.frame(s = rep(s, 2), val = c(v1, v2), Height = c(rep("1.56 m", length(v1)), rep("3.54 m", length(v2))))
  #df = data.frame(s, v1, v2)
  ggplot(df, aes(x =s, y = val, group = Height, color = Height)) +
    geom_line(lwd = 2) +
    ggtitle(expression(r[opt]^2 ~ with~sedimentation~rate)) +
    xlab("Sedimentation rate [cm/kyr]") +
    ylab(expression(r[opt]^2))
  ggsave(paste0("figs/",file_name, file_type))
}
plot_r2_opt(fa_prec, "sbs_R2opt_vs_sedrate_prec", ".png")

#### absolute age of ff boundary ####

plot_age_ff_boundary = function(file_name, file_type){
  l_secc = sapply(admtools::get_time(x = adm_secc_abs, h = h_f_f_bdry ), function(x) x[1])
  l_prec = sapply(admtools::get_time(x = adm_prec_abs, h = h_f_f_bdry ), function(x) x[1])
  df = data.frame(dur = -1 * c(l_prec, l_secc) / 1000, Test = c(rep("Precession", length(l_prec)), rep("Short eccentricity", length(l_secc))))

  ggplot(df, aes(x = dur, fill = Test)) + 
    geom_density(alpha = 0.5) + 
    xlab("Age of F-F boundary [Ma]") +
    ylab("Density") +
    ggtitle(paste0("Age of F-F boundary")) +
    scale_x_reverse()
  ggsave(paste0("figs/", file_name, ".png"))
}

plot_age_ff_boundary("sbs_age_ff_bdry", ".png")

# get the numbers

get_ff_bdr_stats = function(adm){
  l = sapply(admtools::get_time(x = adm, h = h_f_f_bdry ), function(x) x[1])
  df = data.frame(dur = -l/1000)

  return(list("mean" = mean(df$dur),
              "2sig" = 2 * sd(df$dur),
              "median" = median(df$dur),
              "iqr" = IQR(df$dur),
              "quantile" = quantile(df$dur)))
}

prec_ages = get_ff_bdr_stats(adm_prec_abs)

ff_mean_diff = ff_mean_dasilva - prec_ages$mean
ff_2sigma_diff = ff_2sigma_dasilva - prec_ages$`2sig`
ff_2sigma_rat = 1- prec_ages$`2sig` / ff_2sigma_dasilva
ff_mean_diff_da_silva = ff_mean_gradstein - ff_mean_dasilva
ff_mean_diff_gradstein = ff_mean_gradstein - prec_ages$mean

#### Duration of upper Kellwasser event ####

plot_duration_uke = function(file_name, file_type){
  l_prec = sapply(admtools::get_time(x = adm_prec_abs, h = c(h_top_ukw, h_bottom_ukw )),  function(x) x[1] - x[2])
  l_secc = sapply(admtools::get_time(x = adm_secc_abs, h = c(h_top_ukw, h_bottom_ukw )),  function(x) x[1] - x[2])
  df = data.frame(dur = c(l_prec, l_secc), Test = c(rep("Precession", length(l_prec)), rep("Short eccentricity", length(l_secc))))
  
  ggplot(data = df, aes(x = dur, fill = Test) ) + 
    geom_density(alpha = 0.5) + 
    xlab("Duration of upper Kellwasser event [kyr]") + 
    ylab("Density") +
    ggtitle("Duration of upper Kellwasser event") 
  ggsave(paste0("figs/", file_name, file_type))
}
plot_duration_uke("sbs_duration_uke", ".png")

# some numbers
get_uke_duration_stats = function(adm){
  l = sapply(admtools::get_time(x = adm, h = c(h_top_ukw, h_bottom_ukw )),  function(x) x[1] - x[2])
  return(list("mean" = mean(l),
              "median" = median(l),
              "iqr" = IQR(l),
              "quantile" = quantile(l)))
}
uke_stats = get_uke_duration_stats(adm_prec_abs)

#### time elapsed between bentonite layer and the ff boundary ####

plot_elapsed_time_ke = function(file_name, file_type){
  l_prec = sapply(admtools::get_time(x = adm_prec_abs, h = c(h_bottom_ukw, h_top_lkw )),  function(x) x[1] - x[2])
  l_secc = sapply(admtools::get_time(x = adm_secc_abs, h = c(h_bottom_ukw, h_top_lkw )),  function(x) x[1] - x[2])
  df = data.frame(dur = c(l_prec, l_secc), Test = c(rep("Precession", length(l_prec)), rep("Short eccentricity", length(l_secc))))
  
  ggplot(data = df, aes(x = dur, fill = Test) ) + 
    geom_density(alpha = 0.6) + 
    xlab("Time elapsed between Kellwasser events [kyr]") + 
    ylab("Density") +
    ggtitle("Time elapsed between Kellwasser Events") 
  ggsave(paste0("figs/", file_name, file_type))
  
}
plot_elapsed_time_ke("sbs_time_elapsed_between_ke", ".png")


## some numbers on timing
time_between_ke_stats = function(adm){
  l = sapply(admtools::get_time(x = adm, h = c(h_bottom_ukw, h_top_lkw )),  function(x) x[1] - x[2])
  return(list("mean" = mean(l),
              "median" = median(l),
              "iqr" = IQR(l),
              "quantile" = quantile(l),
              "sd" = sd(l)))
}

time_ke_stats = time_between_ke_stats(adm_prec_abs)
