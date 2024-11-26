data = read.csv("data/raw/murphy_et_al_2010_1-s2.0-S0016703710003108-mmc3.csv", header = TRUE, sep = "\t")

data$Depth..mcd.

h = data$Depth..mcd.

plot(h, f(h))

base_clay = 306.78

t_tp = function(){
  return(0)
}

h_tp = function(){
  return(base_clay)
}

# based on the number in the main text
time_const_gen_text = function(){
  eps = 0.0001
  r = max(eps,rnorm(1, mean = 0.37, sd = 0.06/2)) # murphy et al 2010, pcc/cm^-2/kyr
  f = approxfun(x = c(-1000, 1000), y = rep(r, 2), rule = 2)
  return(f)
}
 # based on the number in the suppl. materials
time_const_gen_supp = function(){
  eps = 0.0001
  r = max(eps,rnorm(1, mean = 0.48, sd = 0.08/2)) # murphy et al 2010, pcc/cm^-2/kyr
  f = approxfun(x = c(-1000, 1000), y = rep(r, 2), rule = 2)
  return(f)
}

f = time_const_gen()

plot(1:300, time_const_gen()(1:300))


strat_cont_gen = function(){
  f = approxfun(x = data$Depth..mcd.,
                y = data$X3HeET..pcc.g...20. * data$DBD..g.cm3. * 100,
                rule = 2)
  return(f)
}

library(admtools)
adm = strat_cont_to_multiadm(h_tp = h_tp,
                             t_tp = t_tp,
                             strat_cont_gen = strat_cont_gen,
                             time_cont_gen = time_const_gen_supp,
                             h = seq(304, 307, by = 0.01),
                             subdivisions = 10000, stop.on.error = FALSE,
                             no_of_rep = 500)
plot(adm)

aa = get_time(adm, h = c(306.78, 306.15))

x = sapply(aa, diff)

quantile(x)

hist(x)
