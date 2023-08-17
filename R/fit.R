#########################
# Script to fit the DFE #
#########################

library(here)
library(tidyverse)
source(here("R","dfe_functions.R"))

data_DFE <- read_delim(here("Data", "database_clusters_reduced.csv"), delim = ";",
                            col_types = "_fffdddff",
                            locale=locale(decimal_mark = ","))
                       
data_DFE <- data_DFE %>%
  filter(AA != "WILDTYPE") %>% # remove the wildtype as we want ony mutants for the DFE
  mutate(sd = abs(Upper_CI - Low_CI) / 4) #approximation of the standard deviation by the standard error.

# Count of mutationsin each environment
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

