library(tidyverse)
library(lubridate)
library(car)
library(caret)

processfun = function(species, season){
  df = read.csv(paste0("data/necrosis_data_", season, ".csv")) %>% 
    filter(spp == species)
  names(df)[1] = "uniqueID"
  
  df$spc1 = df$EC1*(1 + 0.0191*(df$temp1 - 25))
  df$spc2 = df$EC2*(1 + 0.0191*(df$temp2 - 25))
  df$percentEC = as.numeric((df$spc1/df$spc2)*100)
  
  df.control = df |> 
    filter(temp == "25") |> 
    dplyr::select(spp, age, branch, percentEC)
  names(df.control)[names(df.control) == "percentEC"] = "EC.control"
  
  df.treatment = df |> 
    filter(temp != "25") |> 
    dplyr::select(spp, age, temp, dur1, branch, spc2, percentEC)
  names(df.treatment)[names(df.treatment) == "percentEC"] = "EC.treatment"
  names(df.treatment)[names(df.treatment) == "dur1"] = "duration"
  
  dfjoin = inner_join(df.control, df.treatment)
  
  dfjoin$survival.treatment = 100 - dfjoin$EC.treatment
  dfjoin$survival.control = 100 - dfjoin$EC.control
  dfjoin$survival_final = (dfjoin$survival.treatment/dfjoin$survival.control)*100
  
  dfjoin$treatment = paste0(dfjoin$temp, "_", dfjoin$duration, "_", dfjoin$age)
  dfjoin$temp = as.factor(dfjoin$temp)
  dfjoin$age = as.factor(dfjoin$age)
  dfjoin$duration = as.numeric(gsub(x = dfjoin$duration, pattern = ",", replacement = ""))
  
  df.agg.basic = (sapply(split(dfjoin, dfjoin$treatment), function(d) mean(d$survival_final)))
  df.agg.basic = as.data.frame(df.agg.basic)
  df.agg.basic$treatment = rownames(df.agg.basic)
  df.agg.basic$temp = unlist(map(str_split(df.agg.basic$treatment, "_", n = 3), 1))
  df.agg.basic$duration = unlist(map(str_split(df.agg.basic$treatment, "_", n = 3), 2))
  df.agg.basic$duration = as.numeric(gsub(x = df.agg.basic$duration, pattern = ",", replacement = ""))
  df.agg.basic$age = unlist(map(str_split(df.agg.basic$treatment, "_", n = 3), 3))
  names(df.agg.basic)[1] = "survival_final"
  df.agg.basic$temp = as.factor(df.agg.basic$temp)
  df.agg.basic$duration = as.numeric(df.agg.basic$duration)
  df.agg.basic$age = as.factor(df.agg.basic$age)
  
  df = df %>% 
    dplyr::select(uniqueID:branch, spc1, EC1, temp1, spc2, EC2, temp2, percentEC)
  write.csv(df, paste0("data/necrosis_data_", species, "_", season, ".csv"), row.names = F)
  
  dfjoin = dfjoin %>% 
    dplyr::select(spp, branch, age, temp, duration, treatment, EC.control, survival.control, EC.treatment, survival.treatment, survival_final)
  write.csv(dfjoin, paste0("data/survival_data_", species, "_", season, ".csv"), row.names = F)
}
processfun(species = "PIPO", season = "spring")
processfun(species = "PIPO", season = "autumn")
processfun(species = "PSME", season = "spring")
processfun(species = "PSME", season = "autumn")

scalefun = function(species, season){
  df = read.csv(paste0("data/survival_data_", species, "_", season, ".csv"))
  df = df |> 
    filter(spp == species)
  
  df$groupVar = paste0(df$temp, "_", df$age)
  
  df$treatment = paste0(df$temp, "_", df$duration, "_", df$age)
  df$temp = as.factor(df$temp)
  
  df.agg.basic = (sapply(split(df, df$treatment), function(d) mean(d$survival_final)))
  df.agg.basic = as.data.frame(df.agg.basic)
  df.agg.basic$treatment = rownames(df.agg.basic)
  df.agg.basic$temp = unlist(map(str_split(df.agg.basic$treatment, "_", n = 3), 1))
  df.agg.basic$duration = unlist(map(str_split(df.agg.basic$treatment, "_", n = 3), 2))
  df.agg.basic$age = unlist(map(str_split(df.agg.basic$treatment, "_", n = 3), 3))
  names(df.agg.basic)[1] = "survival_final"
  df.agg.basic$temp = as.factor(df.agg.basic$temp)
  df.agg.basic$duration = as.numeric(df.agg.basic$duration)
  df.agg.basic$age = as.factor(df.agg.basic$age)
  df.agg.basic$groupVar = paste0(df.agg.basic$temp, "_", df.agg.basic$age)
  
  # after scaling ####
  l = list()
  for(i in unique(df$groupVar)){
    temp = df |> 
      filter(groupVar == i)
    ref.temp = df.agg.basic |> 
      filter(groupVar == i)
    reftransform = preProcess(ref.temp |> dplyr::select(survival_final), method = c("range"), rangeBounds = c(0, 1))
    
    temp.transform = predict(reftransform, temp |> dplyr::select(survival_final))
    
    temp$survival_final = temp.transform$survival_final*100
    
    l[[i]] = temp
  }
  df_scaled = bind_rows(l)
  
  df_scaled$treatment = paste0(df_scaled$temp, "_", df_scaled$duration, "_", df_scaled$age)
  df_scaled$temp = as.factor(df_scaled$temp)
  df_scaled$season = season
  
  write.csv(df_scaled, paste0("data/necrosis_data_", species, "_", season, "_scaled.csv"), row.names = F)
}
scalefun("PIPO", "spring")
scalefun("PIPO", "autumn")
scalefun("PSME", "spring")
scalefun("PSME", "autumn")

df1 = read.csv(paste0("data/necrosis_data_PIPO_spring_scaled.csv"))
df2 = read.csv(paste0("data/necrosis_data_PIPO_autumn_scaled.csv"))
df3 = read.csv(paste0("data/necrosis_data_PSME_spring_scaled.csv"))
df4 = read.csv(paste0("data/necrosis_data_PSME_autumn_scaled.csv"))
df = rbind(df1, df2, df3, df4) |>
  dplyr::select(spp, season, age, temp, duration, survival_final)
names(df)[names(df) %in% c("spp", "temp", "survival_final")] = c("species", "temperature", "survival")
df$groupVar = paste0(df$species, "_", df$season, "_", df$age)
write.csv(df, "data/necrosis_data_scaled.csv", row.names = F)
