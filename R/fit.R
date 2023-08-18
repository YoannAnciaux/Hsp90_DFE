#########################
# Script to fit the DFE #
#########################

library(here)
library(tidyverse)
library(fitdistrplus)
library(cowplot)
library(gridExtra)
source(here("R","dfe_functions.R"))

#### Import data and compute statistics for fit ####
data_DFE <- read_delim(here("Data", "database_clusters_reduced.csv"), delim = ";",
                       col_types = "_fffdddff",
                       locale=locale(decimal_mark = ","))

data_DFE <- data_DFE %>%
  filter(AA != "WILDTYPE") %>% # remove the wildtype as we want ony mutants for the DFE
  mutate(sd = abs(Upper_CI - Low_CI) / 4) #approximation of the standard deviation by the standard error.

# Count of mutations in each environment
stats_all <- data_DFE %>%
  group_by(Environment) %>%
  summarise(number = n())

# By environment, number of AA leading to codon stop and their mean and sd
stats_stop <-
  filter(data_DFE, AA == "*") %>%
  group_by(Environment) %>%
  summarise(mean_stop = mean(s), sd_stop = sd(s), number_stop = n())

# Regroup the parameter of the gaussian for the mixture by environement.
# The mean and the sd of the gaussian is the mean s of aa stop. Prob is the
# the proportion of aa no stop.
par_ini_norm <- left_join(stats_stop, stats_all) %>%
  transmute(Environment = Environment,
            mean = mean_stop,
            sd = sd_stop,
            prob = 1-number_stop/number)

envs <- data_DFE$Environment %>% unique() %>% as.vector()
names(envs) <- envs

#### MGE fit ####
# Fit DFE of all AA in each environment using a goodness of fit method based on
# the Cramer-von Mises distance. It gives results which are more plausible than the
# MLE method for the number of dimension in both 37C and Diamide environemnts. The
# absolute value of the parameters is a bit different between MLE and MGE for the
# other environments but the relative values of these parameters between environments
# stay similar.
# We fix the mean, and the sd of the gaussian to the mean and sd of the codon stops
# selection coefficients to reduce the number of parameters to fit.
# Here we use a mixture of the FGM DFE and a gaussian but lambda is not fitted, it
# is proportionnal to the n_dim and the mean of the data which are considered not
# in the gaussian.
mgefit <- sapply(envs,
                 function(e) {
                   print(e)
                   par_ini_norm_env <- filter(par_ini_norm, Environment == e);
                   data_env <- filter(data_DFE, Environment == e)$s
                   
                   fit_env <- fitdist(data_env, "normfgmsmix", method = "mge", optim.method = "default",
                                      start = list(n_dim = 2, so = max(data_env)*1.001, prob = par_ini_norm_env$prob),
                                      fix.arg = list(mean = par_ini_norm_env$mean, sd = par_ini_norm_env$sd, mean_data = mean(data_env)),
                                      lower = c(1., 0., (1-mean(data_env)/par_ini_norm_env$mean)*1.001), upper = c(10000., 10000., 1.))
                   
                   grid <- seq(min(data_env), max(data_env), length = 100)
                   dfe_fit_env <-  with(as.list(fit_env$estimate), {
                     par_ini_norm_env <- filter(par_ini_norm, Environment == e);
                     mean = par_ini_norm_env$mean;
                     sd = par_ini_norm_env$sd;
                     tibble(Environment = rep(e, length(grid)),
                            s = grid,
                            density = dnormfgmsmix(grid, mean, sd, n_dim, so, prob, mean(data_env)))
                   })
                   return(list(fit = fit_env, dfe = dfe_fit_env))
                 },
                 simplify = F,
                 USE.NAMES = T)


saveRDS(mgefit, here("Analysis", "mgefit_alldata_scaledlambda.rds"))


#### MLE fit ####
# mle fitting does no reach acceptable n_dim maybe because of lambda constrained.
# mlefit <- sapply(envs,
#                  function(e) {
#                    print(e)
#                    par_ini_norm_env <- filter(par_ini_norm, Environment == e);
#                    data_env <- filter(data_DFE, Environment == e)$s
#                    weights_env <- floor(1 / (filter(data_DFE, Environment == e)$sd)^2)
#                    
#                    fit_env <- fitdist(data_env, "normfgmsmix", method = "mle", optim.method = "default", #weights = weights_env,
#                                       start = list(n_dim = 2, so = max(data_env)*1.001, prob = par_ini_norm_env$prob),
#                                       fix.arg = list(mean = par_ini_norm_env$mean, sd = par_ini_norm_env$sd, mean_data = mean(data_env)),
#                                       lower = c(1., 0., (1-mean(data_env)/par_ini_norm_env$mean)*1.001), upper = c(10000., 10000., 1.))
#                    
#                    grid <- seq(min(data_env), max(data_env), length = 100)
#                    dfe_fit_env <-  with(as.list(fit_env$estimate), {
#                      par_ini_norm_env <- filter(par_ini_norm, Environment == e);
#                      mean = par_ini_norm_env$mean;
#                      sd = par_ini_norm_env$sd;
#                      tibble(Environment = rep(e, length(grid)),
#                             s = grid,
#                             density = dnormfgmsmix(grid, mean, sd, n_dim, so, prob, mean(data_env)))
#                    })
#                    return(list(fit = fit_env, dfe = dfe_fit_env))
#                  },
#                  simplify = F,
#                  USE.NAMES = T)
# 
# 
# saveRDS(mlefit, here("Analysis", "mlefit_alldata_scaledlambda.rds"))


