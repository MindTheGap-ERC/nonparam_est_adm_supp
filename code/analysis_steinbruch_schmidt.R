#### Load packages ####
library(astrochron)
library(admtools)
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
etimeOptMS_prec=eTimeOpt(MSi,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6, numsed=100,linLog=1,
                     limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                     detrend=T,ydir=1,output=1,genplot=F,check=T,verbose=1)
# test for short eccentricty amplitude modulation
etimeOptMS2_secc=eTimeOpt(MSi,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6, numsed=100,linLog=1,
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
fa1 = get_data_from_eTimeOpt(etimeOptMS_prec, index = 3)
# generate sed rate generator
se_prec = sed_rate_from_matrix(height = fa1$heights,
                               sedrate = fa1$sed_rate / 100, # convert cm/kyr to m/kyr
                               matrix = fa1$results,
                               rate = rate)
fa2 = get_data_from_eTimeOpt(etimeOptMS2_secc, index = 3)
se_secc = sed_rate_from_matrix(height = fa2$heights,
                               sedrate = fa2$sed_rate / 100, # convert cm/kyr to m/kyr
                               matrix = fa2$results,
                               rate = rate)

# heights of interest
h = seq(1.05, 4.25, by = 0.05)

#### estimate ADMs ####
adm_prec_abs = sedrate_to_multiadm(h_tp , t_tp_absolute, sed_rate_gen = se_prec, h, no_of_rep = no_of_rep )
adm_prec_float = sedrate_to_multiadm(h_tp , t_tp_float, sed_rate_gen = se_prec, h, no_of_rep = no_of_rep )

adm_secc_abs = sedrate_to_multiadm(h_tp , t_tp_absolute, sed_rate_gen = se_secc, h = h, no_of_rep = no_of_rep )
adm_secc_float = sedrate_to_multiadm(h_tp , t_tp_float, sed_rate_gen = se_secc, h = h, no_of_rep = no_of_rep )

#### Plotting ####
# auxiliary function

plot_sbs = function(adm){
  q2_adm = quantile_adm(adm, h, 0.925)
  q1_adm = quantile_adm(adm, h,  0.025)
  m_adm = quantile_adm(adm, h, 0.5)
  
  t_min = min(range(c(q2_adm$t, q1_adm$t, m_adm$t)))/1000
  t_max = max(range(c(q2_adm$t, q1_adm$t, m_adm$t)))/1000
  h_min = min(range(c(q1_adm$h, q2_adm$h, m_adm$h)))
  h_max = max(range(c(q1_adm$h, q2_adm$h, m_adm$h)))
  
  plot(NULL,
       xlim = c(t_min, t_max),
       ylim = c(h_min, h_max),
       xlab = "Time BP [Myr]",
       ylab = "Stratigraphic Position [m]",
       main = "ADM for Steinbruch Schmidt (Absolute Age)")
  col_med = "red"
  col_env = "blue"
  lwd_med = 3
  lwd_env = 3
  rect(xleft = t_min, xright = t_max, ytop = h_top_lkw, ybottom = h_min - 1, col = grey(0.6), lty = 0)
  text(x = mean(c(t_min, t_max)), y = mean(c(h_top_lkw, h_top_lkw - 0.1)), labels = "Lower Kellwasser Bed")
  
  rect(xleft = t_min, xright = t_max, ytop = h_top_ukw, ybottom = h_bottom_ukw, col = grey(0.6), lty = 0)
  text(x = mean(c(t_min, t_max)), y = mean(c(h_top_ukw, h_bottom_ukw)), labels = "Upper Kellwasser Bed")
  lines(m_adm$t/1000, m_adm$h, col = col_med, lwd = lwd_med)
  lines(q1_adm$t/1000, q1_adm$h, col = col_env, lwd = lwd_env)
  lines(q2_adm$t/1000, q2_adm$h, col = col_env, lwd = lwd_env)
  
  lines(c(t_min, t_max), rep(h_ashbed, 2), lty = 3, lwd = 3, col = grey(0.6))
  text(x = t_max, y = h_ashbed + 0.1, labels = "Bentonite Layer", pos = 2)
  
  legend("topleft", col = c(col_med, col_env), lwd = c(lwd_med, lwd_env), lty = 1, legend = c("Median Age", "95 % Envelope"))
  
}

png(file = "figs/sds_floating_adm_prec.png")
plot_sbs(adm_prec_float)
dev.off()

png(file = "figs/sds_anchored_adm_prec.png")
plot_sbs(adm_prec_abs)
dev.off()


png(file = "figs/sds_floating_adm_secc.png")
plot_sbs(adm_secc_float)
dev.off()

png(file = "figs/sds_anchored_adm_secc.png")
plot_sbs(adm_secc_abs)
dev.off()



#### timing of f-f boundary ####
# time elapsed between bentonite layer and the ff boundary
png("figs/time_elapsed_to_ff_secc.png")
dl = sapply(get_time(x = adm_secc_float, h = h_f_f_bdry), function(x) x[1])
hist(dl, main = "Time between bentonite and ff boundary \n using secc test",
     xlab = "kyr")
dev.off()

png("figs/time_elapsed_to_ff_prec.png")
dl = sapply(get_time(x = adm_prec_float, h = h_f_f_bdry), function(x) x[1])
hist(dl, main = "Time between bentonite and ff boundary \n using prec test",
     xlab = "kyr")
dev.off()

bent_ff_quantile = quantile(dl)
bent_ff_iqr = IQR(dl)
bent_ff_mean = mean(dl)
bent_ff_sd = sd(dl)

# absolute age of ff bdry
l = sapply(get_time(x = adm_secc_abs, h = h_f_f_bdry ), function(x) x[1])
df = data.frame(dur = -l/1000)

ggplot(df, aes(x = dur)) + 
  geom_density() + 
  xlab("Age of F-F boundary [Ma]") +
  ylab("Density") +
  ggtitle("Age of F-F boundary under secc test") +
  scale_x_reverse()
ggsave("figs/age of ff boundary secc.png")

l = sapply(get_time(x = adm_prec_abs, h = h_f_f_bdry ), function(x) x[1])
df = data.frame(dur = -l/1000)

ggplot(df, aes(x = dur)) + 
  geom_density() + 
  xlab("Age of F-F boundary [Ma]") +
  ylab("Density") +
  ggtitle("Age of F-F boundary under prec test") +
  scale_x_reverse()
ggsave("figs/age of ff boundary prec.png")

ff_mean = mean(df$dur)
ff_2sig = 2 * sd(df$dur)
ff_median = median(df$dur)
ff_iqr = IQR(df$dur)
ff_quantile = quantile(df$dur)


# difference to da silva age
ff_mean_diff = ff_mean_dasilva - ff_mean
ff_2sigma_diff = ff_2sigma_dasilva - ff_2sig
ff_2sigma_rat = 1- ff_2sig / ff_2sigma_dasilva
ff_mean_diff_da_silva = ff_mean_gradstein - ff_mean_dasilva
ff_mean_diff_gradstein = ff_mean_gradstein - ff_mean


#### Duration of upper Kellwasser event ####

l_prec = sapply(get_time(x = adm_prec_abs, h = c(h_top_ukw, h_bottom_ukw )),  function(x) x[1] - x[2])
l_secc = sapply(get_time(x = adm_secc_abs, h = c(h_top_ukw, h_bottom_ukw )),  function(x) x[1] - x[2])

df = data.frame(dur = c(l_prec, l_secc), 
                type = c(rep("prec", length(l_prec)), rep("secc", length(l_secc))))

ggplot(data = df, aes(x = dur, fill = type) ) + 
  geom_density() + 
  xlab("Duration of upper Kellwasser event [kyr]") + 
  ylab("Density") +
  ggtitle("Duration of upper Kellwasser event") + 
  xlim(70, 160)
ggsave("figs/durations_UKE.png")

quantile(l_prec)
quantile(l_secc)

l_prec = sapply(get_time(x = adm_prec_abs, h = c(h_top_ukw, h_bottom_ukw )),  function(x) x[1] - x[2])

df = data.frame(dur = c(l_prec), 
                type = c(rep("prec", length(l_prec))))

ggplot(data = df, aes(x = dur, fill = type) ) + 
  geom_density() + 
  xlab("Duration of upper Kellwasser event [kyr]") + 
  ylab("Density") +
  ggtitle("Duration of upper Kellwasser event") + 
  xlim(70, 160)
ggsave("figs/duration UKE.png")


uke_mean = mean(df$dur)
uke_sd = sd(df$dur)
uke_median = median(df$dur)
uke_quantile = quantile(df$dur)
uke_iqr = IQR(df$dur)


#### Time between lower and upper Kellwasser event ####
l_prec = sapply(get_time(x = adm_prec_abs, h = c(h_bottom_ukw, h_top_lkw )),  function(x) x[1] - x[2])
l_secc = sapply(get_time(x = adm_secc_abs, h = c( h_bottom_ukw , h_top_lkw)),  function(x) x[1] - x[2])

df = data.frame(dur = c(l_prec, l_secc), 
                type = c(rep("prec", length(l_prec)), rep("secc", length(l_secc))))

ggplot(data = df, aes(x = dur, fill = type) ) + 
  geom_density() + 
  xlab("Time elapsed between Kellwasser events [kyr]") + 
  ylab("Density") +
  ggtitle("Time elapsed between Kellwasser events") + 
  xlim(400, 800)



l_prec = sapply(get_time(x = adm_prec_abs, h = c(h_bottom_ukw, h_top_lkw )),  function(x) x[1] - x[2])

df = data.frame(dur = c(l_prec), 
                type = c(rep("prec", length(l_prec))))

ggplot(data = df, aes(x = dur, fill = type) ) + 
  geom_density() + 
  xlab("Time elapsed between Kellwasser events [kyr]") + 
  ylab("Density") +
  ggtitle("Time elapsed between Kellwasser events") + 
  xlim(400, 800)
ggsave("time between kellwasser events.png")
# Kellwasser elapsed time
ket_mean = mean(df$dur)
ket_median = median(df$dur)
ket_sd= sd(df$dur)
ket_quantile = quantile(df$dur)
ket_iqr = IQR(df$dur)


