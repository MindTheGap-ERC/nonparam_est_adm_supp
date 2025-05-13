sed_rate_plot = function(file_name){
  df = data.frame(h = rep(h_eval, 3), sedr = c(sedr_const, sedr_dec, sedr_inc),
                  Scenario = c(rep("Constant flux", length(sedr_const)),rep("Decreasing flux", length(sedr_inc)), rep("Increasing flux", length(sedr_dec))))
  m_sedr = max(df$sedr)
  
  rect = data.frame(h_min = c(core_interval[1], recovery1_interval[1], recovery2_interval[1], recovery3_interval[1]),
                    h_max = c(core_interval[2], recovery1_interval[2], recovery2_interval[2], recovery3_interval[2]),
                    ymax = rep(m_sedr, 4),
                    ymin = rep(0, 4))
  plt = ggplot(df, aes(x = h, y = sedr, group = Scenario, col = Scenario)) +
    ylim(c(0, m_sedr)) +
    geom_rect(data = rect, inherit.aes = FALSE, aes(xmin = h_min, xmax = h_max, ymin = ymin, ymax = ymax), fill = box_cols) +
    scale_color_manual(values = scenario_cols) +
    geom_line(size = 0.5, alpha = 0.7) +
    xlab("Depth [m composite depth]") +
    ylab("Sedimentation rate [cm/kyr]") +
    annotate("text", x = mean(core_interval), y = 7.5, angle = 90, label = "Core interval", size = annot_size/.pt)  +
    annotate("text", x = mean(recovery1_interval), y = 7.5, label = "Recovery I", angle = 90, size = annot_size/.pt) +
    annotate("text", x = mean(recovery2_interval), y = 7.5, label = "Recovery II", size = annot_size/.pt) +
    annotate("text", x = mean(recovery3_interval), y = 0.5, label = "Recovery III", size = annot_size/.pt) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.1, 0.9),
          legend.title = element_blank(),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = legend_size),
          axis.title = element_text(size = lab_size),
          plot.title = element_text(size = title_size),
          axis.text = element_text(size = ax_size)) +
    ggtitle("Sedimentation rate")
  ggsave(paste0("figs/",file_name, ".png"), plot = plt)
  return(plt)
}