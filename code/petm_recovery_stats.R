petm_recovery_stats = function(){
  clay_int = c(base_clay, clay_layer_top)
  core_int = c(base_clay, base_recovery)
  recovery_int = c(base_recovery, recovery_3_top)
  PETM = c(base_clay, recovery_3_top)
  l = list()
  for (i in c("const_det", "inc_det", "dec_det")){
    adm = adm_list[[i]]
    l[[paste0(i, "_median_clay")]] = median(sapply(get_time(adm, h = rev(clay_int)), diff ))
    l[[paste0(i, "_iqr_clay")]] = IQR(sapply(get_time(adm, h = rev(clay_int)), diff ))
    l[[paste0(i, "_median_recovery")]] = median(sapply(get_time(adm, h = rev(recovery_int)), diff ))
    l[[paste0(i, "_iqr_recovery")]] = IQR(sapply(get_time(adm, h = rev(recovery_int)), diff ))
    l[[paste0(i, "_median_PETM")]] = median(sapply(get_time(adm, h = rev(PETM)), diff ))
    l[[paste0(i, "_2sigma_PETM")]] = 2 * sd(sapply(get_time(adm, h = rev(PETM)), diff ))
    l[[paste0(i, "_CI95_PETM")]] = quantile(sapply(get_time(adm, h = rev(clay_int)), diff ), probs=c(0.05, 0.95))
    l[[paste0(i, "_2sd_clay_int")]] = 2*sd(sapply(get_time(adm, h = rev(clay_int)), diff ))
    l[[paste0(i, "_median_core")]] = median(sapply(get_time(adm, h = rev(core_int)), diff ))
    l[[paste0(i, "_2sd_core")]] = 2*sd(sapply(get_time(adm, h = rev(core_int)), diff ))
  }
  return(l)
}