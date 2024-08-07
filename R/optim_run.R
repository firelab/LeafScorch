library(tidyverse)
library(minpack.lm)
library(patchwork)
library(grid)
library(kimisc)
library(RColorBrewer)
library(car)
library(multcomp)
library(optimization)
library(viridis)
library(plotrix)

# universal constants ####
kB = 1.38*(10^-23)      ## Boltzman constant
h = 6.63*(10^-34)       ## Planck's constant
R = 8.31                ## universal gas constant 
a1 = 1/329              ## compensation law for the average protein, Rosenberg et al. (1971)
b1 = -64.9*4.184        ## empirically-derived intercept for compensation law for the average protein, Rosenberg et al. (1971)

# parameters from Dickinson and Johnson (2004) ####
df = data.frame(temperature = rep(seq(45, 65, by = 5), 2),
                species = c(rep("PSME", 5), rep("PIPO", 5)),
                deltaH = c(rep(285164, 5), rep(20195, 2), rep(478828, 3)))
df$cT = exp(c(-10.3, -8.7, -7.1, -5.5, -3.9, -10.625, -10.625, -8.75, -6.25, -3.7))

## all values fall between 1e-6 and 0.1
write.csv(df, "data/initial_parameters.csv", row.names = F)

# initial ####
label = "initial"
results_df = data.frame(species = c(rep("PIPO", 4), rep("PSME", 4)),
                        season = rep(c(rep("spring", 2), rep("autumn", 2)), 2),
                        age = rep(c("o", "n"), 4))
results_df$groupVar = paste0(results_df$species, "_", results_df$season, "_", results_df$age)
results_df$jackVar = paste0(results_df$species, "_", results_df$season, "_", results_df$age, "_1")

observed = read.csv("data/necrosis_data_scaled.csv") |>
  dplyr::select(groupVar, temperature, duration, survival)

# -------------- grid search ####
source('R/optim_functions.R')
grid.df = read.csv("data/initial_parameters.csv")
jackfun(observed, label)
select.jack = read.csv(paste0("data/necrosis_data_scaled_jackknifed_", label, ".csv"))
gridfun(grid.df, select.jack, label)

# -------------- select best values from grid search ####
source('R/optim_functions.R')
parameters = read.csv(paste0("outputs/gridsearch_", label, ".csv"))
gridwinfun(parameters, label)

# -------------- optimization ####
source('R/optim_functions.R')
observed_jk = read.csv(paste0("data/necrosis_data_scaled_jackknifed_", label, ".csv"))
parameters = read.csv(paste0("outputs/gridsearch_best_min_", label, ".csv"))
optimfun(parameters, observed_jk, label)

# -------------- survival curves ####
source('R/optim_functions.R')
parameters = readRDS(paste0("outputs/optimization_", label, ".rds"))
predictions.smooth = predictfun(parameters, observed, label)

# -------------- log-linear rate parameter models and significance ####
source('R/optim_functions.R')
parameters = readRDS(paste0("outputs/optimization_", label, ".rds"))
errors = readRDS(paste0("outputs/predictions.smoothed_", label, ".rds"))
loglinearfun(parameters, errors, label)

results_df = results_df %>% dplyr::select(-jackVar)

cT = read.csv(paste0("outputs/cT.loglinearmodelresults_", label, ".csv"))
pT = read.csv(paste0("outputs/pT.loglinearmodelresults_", label, ".csv"))
lplotsinitfun(cT, pT, parameters, label)

# -------------- grouping decision ####
# P. ponderosa
## merge old foliage categories

# P. menziesii
## merging autumn categories

# second- optimization run only for merged groups ####
label = "second"

results_df = data.frame(species = c(rep("PIPO", 3), rep("PSME", 3)),
                        season = c("spring", "autumn", "spring and autumn", "spring", "spring", "autumn"),
                        age = c("n", "n", "o", "n", "o", "n and o"),
                        groupVar = c("PIPO_spring_n", "PIPO_autumn_n", "PIPO_allold",
                                     "PSME_spring_n", "PSME_spring_o", "PSME_allautumn"))

observed = read.csv("data/necrosis_data_scaled.csv") |>
  dplyr::select(groupVar, temperature, duration, survival)
observed$groupVar = as.factor(observed$groupVar)
levels(observed$groupVar) = c("PIPO_autumn_n", "PIPO_allold", "PIPO_spring_n", "PIPO_allold",
                              "PSME_allautumn", "PSME_allautumn", "PSME_spring_n", "PSME_spring_o")

# -------------- grid search ####
source('R/optim_functions.R')
grid.df = read.csv("data/initial_parameters.csv")
jackfun(observed, label)
select.jack = read.csv(paste0("data/necrosis_data_scaled_jackknifed_", label, ".csv"))
gridfun(grid.df, select.jack, label)

# -------------- select best values from grid search ####
source('R/optim_functions.R')
parameters = read.csv(paste0("outputs/gridsearch_", label, ".csv"))
gridwinfun(parameters, label)

# -------------- optimization ####
source('R/optim_functions.R')
observed = read.csv(paste0("data/necrosis_data_scaled_jackknifed_", label, ".csv"))
parameters = read.csv(paste0("outputs/gridsearch_best_min_", label, ".csv"))
optimfun(parameters, observed, label)

# -------------- survival curves ####
source('R/optim_functions.R')
parameters = readRDS(paste0("outputs/optimization_", label, ".rds"))
predictions.smooth = predictfun(parameters, observed, label)

# -------------- log-linear rate parameter models and significance ####
source('R/optim_functions.R')
parameters = readRDS(paste0("outputs/optimization_", label, ".rds"))
errors = readRDS(paste0("outputs/predictions.smoothed_", label, ".rds"))
loglinearsecfun(parameters, errors, label)

cT = read.csv(paste0("outputs/cT.loglinearmodelresults_", label, ".csv"))
pT = read.csv(paste0("outputs/pT.loglinearmodelresults_", label, ".csv"))
lplotssecfun(cT, pT, parameters, label)

# -------------- grouping decision ####
# P. ponderosa
## merge autumn new foliage with all old foliage

# P. menziesii
## no changes
# final- optimization run only for merged groups ####
label = "final"
results_df = data.frame(species = "PIPO",
                        season = "spring and autumn",
                        age = "n and o",
                        groupVar = "PIPO_allothers")

results_fin = data.frame(species = c(rep("PIPO", 2), rep("PSME", 3)),
                         season = c("spring", "spring and autumn", "spring", "spring", "autumn"),
                         age = c("n", "n and o", "n", "o", "n and o"),
                         groupVar = c("PIPO_spring_n", "PIPO_allothers", "PSME_spring_n", "PSME_spring_o", "PSME_allautumn"))

observed = read.csv("data/necrosis_data_scaled.csv") |>
  dplyr::select(groupVar, temperature, duration, survival)
observed$groupVar = as.factor(observed$groupVar)
levels(observed$groupVar) = c("PIPO_allothers", "PIPO_allothers", "PIPO_spring_n", "PIPO_allothers",
                              "PSME_allautumn", "PSME_allautumn", "PSME_spring_n", "PSME_spring_o")
observed = observed |> 
  left_join(results_df) |> 
  filter(groupVar == "PIPO_allothers")

# -------------- grid search ####
source('R/optim_functions.R')
grid.df = read.csv("data/initial_parameters.csv")
jackfun(observed, label)
select.jack = read.csv(paste0("data/necrosis_data_scaled_jackknifed_", label, ".csv"))
gridfun(grid.df, select.jack, label)

# -------------- select best values from grid search ####
source('R/optim_functions.R')
parameters = read.csv(paste0("outputs/gridsearch_", label, ".csv"))
gridwinfun(parameters, label)

# -------------- optimization ####
source('R/optim_functions.R')
observed = read.csv(paste0("data/necrosis_data_scaled_jackknifed_", label, ".csv"))
parameters = read.csv(paste0("outputs/gridsearch_best_min_", label, ".csv"))
optimfun(parameters, observed, label)

# -------------- survival curves ####
source('R/optim_functions.R')
parameters = readRDS(paste0("outputs/optimization_", label, ".rds"))

predictions.smooth = predictfun(parameters, observed, label)
predictions.smooth1 = readRDS(paste0("outputs/predictions.smoothed_", label, ".rds"))
predictions.smooth2 = readRDS(paste0("outputs/predictions.smoothed_second.rds")) |>
  filter(groupVar != "PIPO_allold" & groupVar != "PIPO_autumn_n")
predictions.smooth = rbind(predictions.smooth1, predictions.smooth2)
saveRDS(predictions.smooth, paste0("outputs/predictions.smoothed_", label, "_all.rds"))
predictions.smooth = readRDS(paste0("outputs/predictions.smoothed_", label, "_all.rds"))

observed_fin = read.csv("data/necrosis_data_scaled.csv") |>
  dplyr::select(groupVar, temperature, duration, survival)

survivalplotsfinalfun(observed_fin, predictions.smooth)
options(scipen = 10000)
survivaldifffun(observed_fin, predictions.smooth)
SPPsurvivaldifffun(observed_fin, predictions.smooth)

# -------------- log-linear rate parameter models and significance ####
source('R/optim_functions.R')

parameters1 = readRDS(paste0("outputs/optimization_", label, ".rds"))
parameters2 = readRDS(paste0("outputs/optimization_second.rds")) |>
  filter(groupVar != "PIPO_allold" & groupVar != "PIPO_autumn_n")
parameters = rbind(parameters1, parameters2)
saveRDS(parameters, paste0("outputs/optimization_", label, "_all.rds"))

parameters = readRDS(paste0("outputs/optimization_", label, "_all.rds"))
errors = readRDS(paste0("outputs/predictions.smoothed_", label, "_all.rds"))

loglinearfinalfun(parameters, errors, label)

cT = read.csv(paste0("outputs/cT.loglinearmodelresults_", label, ".csv"))
pT = read.csv(paste0("outputs/pT.loglinearmodelresults_", label, ".csv"))
lplotsfinalfun(cT, pT, parameters, label)

# -------------- LD50 ####
source('R/optim_functions.R')
parameters = readRDS(paste0("outputs/optimization_", label, "_all.rds"))

observed_fin = read.csv("data/necrosis_data_scaled.csv") |>
  dplyr::select(groupVar, temperature, duration, survival)
observed_fin$groupVar = as.factor(observed_fin$groupVar)
levels(observed_fin$groupVar) = c("PIPO_allothers", "PIPO_allothers", "PIPO_spring_n", "PIPO_allothers",
                         "PSME_allautumn", "PSME_allautumn", "PSME_spring_n", "PSME_spring_o")
observed_fin = observed_fin |>
  left_join(results_fin)
write.csv(observed_fin, "outputs/necrosis_data_scaled_merged.csv")

observed_fin = read.csv("outputs/necrosis_data_scaled_merged.csv")
ld50fun(parameters, observed_fin)

# -------------- LD50 plots ####
source('R/optim_functions.R')
LD50 = read.csv("outputs/LD50.csv")
ld50plotsfun(LD50)