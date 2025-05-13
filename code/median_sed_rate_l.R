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
    sed_rate_list[[i]] = admtools::sed_rate_l(adm_list[[i]], h)
  }
  l = list()
  for (j in seq_along(h)){
    l[[j]] = sapply(sed_rate_list, function(x) x[j])
  }
  sedr = sapply(l, median)
  return(sedr)
}
