load(file = "data/res/sbs_data.RData")

library(ggplot2)
source("code/thesis_plot_settings.R")

#### Plot ADMs ####

box_col = grey(0.7)
col_med = "red"
col_env = "blue"
env_lwd = 0.5
med_lwd = 1
unc_lwd = 0.4


sbs_bed_annotations = function(rect){
  #' labels for the Kellwasser beds and bentonite layer, placed on the side of
  #' the panel that is free of age-depth model lines (x axis is reversed, so
  #' larger ages are on the left)
  x_old = max(c(rect$t_min, rect$t_max))
  x_young = min(c(rect$t_min, rect$t_max))
  pad = 0.02 * (x_old - x_young)
  list(
    annotate("text", x = x_young + pad, y = mean(c(h_top_lkw, min(rect$h_min))), label = "Lower Kellwasser Bed", hjust = 1, size = annot_size/.pt),
    annotate("text", x = x_old - pad, y = mean(c(h_top_ukw, h_bottom_ukw)), label = "Upper Kellwasser Bed", hjust = 0, size = annot_size/.pt),
    geom_hline(yintercept = h_ashbed),
    annotate("text", x = x_young + pad, y = h_ashbed, label = "Bentonite layer", hjust = 1, vjust = -0.4, col = grey(0.3), size = annot_size/.pt)
  )
}

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
    geom_line(aes(linewidth = group)) +
    scale_linewidth_manual(values = c("A" = env_lwd, "B" = env_lwd, "C" = med_lwd, "D" = unc_lwd, "E" = unc_lwd), guide = "none") +
    scale_color_manual(values = c(col_env, col_med, col_med)) +
    scale_linetype_manual(values = c("solid", "dashed", "solid")) +
    xlab(xlab) +
    ylab("Stratigraphic position [m]") +
    sbs_bed_annotations(rect) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
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
    geom_line(aes(linewidth = group)) +
    scale_linewidth_manual(values = c("A" = env_lwd, "B" = env_lwd, "C" = med_lwd), guide = "none") +
    scale_color_manual(values = c(col_env, col_med)) +
    xlab(xlab) +
    ylab("Stratigraphic position [m]")  +
    sbs_bed_annotations(rect) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
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
plt = egg::ggarrange(adm_anchor_no_rad, adm_anchor,
                     nrow = 2,
                     ncol = 1,
                     labels = LETTERS[1:2],
                     label.args = panel_label_args,
                     draw = FALSE)
save_thesis_fig("figs/sbs_join_adm_prec.png", plot = plt, height_cm = 17)

## short eccentricity
adm_anchor = plot_sbs_adm_abs(adm_secc_abs, "Age [Ma]")
adm_anchor_no_rad = plot_sbs_adm_float_abs(adm_secc_abs_no_error, "Age [Ma]")

plt = egg::ggarrange(adm_anchor_no_rad, adm_anchor,
                     nrow = 2,
                     ncol = 1,
                     labels = LETTERS[1:2],
                     label.args = panel_label_args,
                     draw = FALSE)
save_thesis_fig("figs/sbs_join_adm_secc.png", plot = plt, height_cm = 17)

#### filled contour plots for eTimeOpt ####
#### values of r^2_opt  ###
plot_r2_opt = function(output){
  v1 =  output$results[,which(output$heights == 1.56)]
  v2 =  output$results[,which(output$heights == 3.54)]
  s = output$sed_rate
  df = data.frame(s = rep(s, 2), val = c(v1, v2), Height = c(rep("1.56 m", length(v1)), rep("3.54 m", length(v2))))
  #df = data.frame(s, v1, v2)
  plt = ggplot(df, aes(x =s, y = val, group = Height, color = Height)) +
    geom_line(linewidth = 1) +
    ggtitle(expression(r[opt]^2 ~ at ~ height)) +
    xlab("Sedimentation rate [cm/kyr]") +
    ylab(expression(r[opt]^2)) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.2, 0.9),
          plot.title = element_text(size = title_size),
          axis.title = element_text(size = lab_size),
          axis.text = element_text(size = ax_size),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = legend_size),
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
         legend.text = element_text(size = legend_size),
         legend.title = element_blank())

plt = egg::ggarrange(r2_plot_prec, eTimeOpt_plot_prec,
                     nrow = 2,
                     ncol = 1,
                     labels = LETTERS[1:2],
                     label.args = panel_label_args,
                     draw = FALSE)
save_thesis_fig("figs/r2_plot_prec_join.png", plot = plt, height_cm = 15)

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
  ggtitle(expression(r[opt]^2 ~ at ~ Steinbruch ~ Schmidt)) +
  theme( plot.title = element_text(size = title_size),
         axis.title = element_text(size = lab_size),
         axis.text = element_text(size = ax_size),
         legend.key.size = unit(0.25, "cm"),
         legend.text = element_text(size = legend_size),
         legend.title = element_blank())

plt = egg::ggarrange(r2_plot_secc, eTimeOpt_plot_secc,
                     nrow = 2,
                     ncol = 1,
                     labels = LETTERS[1:2],
                     label.args = panel_label_args,
                     draw = FALSE)
save_thesis_fig("figs/r2_plot_secc_join.png", plot = plt, height_cm = 15)


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
          legend.text = element_text(size = legend_size))
  return(plt)
}

#### Duration of upper Kellwasser event 

plot_duration_uke = function(){
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
          legend.text = element_text(size = legend_size))
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
          legend.text = element_text(size = legend_size))
  return(plt)
  
}

#### Join plot of durations sbs ####
age_ff_plot = plot_age_ff_boundary()
duration_uke_plot = plot_duration_uke()
elapsed_time_plot = plot_elapsed_time_ke()
plt = egg::ggarrange(age_ff_plot, duration_uke_plot, elapsed_time_plot,
                     nrow = 3,
                     ncol = 1,
                     labels = LETTERS[1:3],
                     label.args = panel_label_args,
                     draw = FALSE)
save_thesis_fig("figs/sbs_durations_joint.png", plot = plt, height_cm = 18)

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