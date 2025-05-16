#### Load and plot already estimated data ####

load(file = "data/res/site1266_data.RData")

library(ggplot2)

##### Reference points and intervals #####
# From Murphy et al. 2010, Table 1.
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
  plot(adm_list[[i]], main = i)
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
# no difference between results derived with and without incorporating measurement
# errors in 3He flux measured throughout the section
# -> focus on case with uncertainty from background flux

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
  ggtitle("Age-depth models for PETM at IODP Site 1266") +
  scale_y_reverse()

ggsave("figs/site1266_adm.png",
       plot = last_plot(),
       width = fig_width_cm,
       dpi = dpi,
       unit = "cm",
       height = 10)

