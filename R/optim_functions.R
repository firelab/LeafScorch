jackfun = function(dat, label){
  dat = dat |> 
    left_join(results_df)
  
  select.jack = list()
  for(n in 1:nrow(results_df)){
    group1 = results_df$groupVar[n]
    dat.temp = dat |>
      filter(groupVar == group1)
    
    selection.grid = data.frame(pt1 = combn(1:nrow(dat.temp), 3)[1,],
                                pt2 = combn(1:nrow(dat.temp), 3)[2,],
                                pt3 = combn(1:nrow(dat.temp), 3)[3,])
    
    nrow(selection.grid)
    set.seed(10)
    selection.grid = sample_n(selection.grid, size = 100)
    head(selection.grid)
    nrow(selection.grid)
    
    select.list = list()
    for(i in c(1:nrow(selection.grid))){
      select.pts = selection.grid[i,]
      dat.i = dat.temp[-c(select.pts$pt1, select.pts$pt2, select.pts$pt3),]
      dat.i$jackVar = paste0(group1, "_", i)
      select.list[[i]] = dat.i
    }
    select.jack[[n]] = bind_rows(select.list)
  }
  select.jack = bind_rows(select.jack)
  nrow(select.jack)
  ## 73,800
  write.csv(select.jack, paste0("data/necrosis_data_scaled_jackknifed_", label, ".csv"), row.names = F)
}

gridfun = function(dat, observed, label){
  jack_df = data.frame(jackVar = unique(observed$jackVar))
  params_list = list()
  for(n in 1:nrow(jack_df)){
    group1 = jack_df$jackVar[n]
    
    deltaHs = seq(min(dat$deltaH)*0.9, max(dat$deltaH)*1.1, length.out = 20) # expand potential values
    deltaSs = a1*deltaHs + b1
    pTs = list()
    for(i in c(1:20)){
      pT1 = ((kB*(seq(45, 65, by = 5) + 273.15))/h)*exp(deltaSs[i]/R)*exp(-deltaHs[i]/(R*(seq(45, 65, by = 5) + 273.15)))
      pTs[[i]] = pT1
    }
    pTs = sort(unlist(pTs))
    cTs = exp(seq(log(1e-15), log(10), length.out = 100))
    
    # expand grid to allow for all possible cT and pT values at all temperatures
    params.grid = expand.grid(pT = pTs,
                              cT = cTs,
                              temperature = c(45, 50, 55, 60, 65))
    params.grid$nrmse = NA
    
    df = observed |> 
      filter(jackVar == group1)
    
    params = data.frame()
    for(i in c(1:length(unique(df$temperature)))){
      temp = sort(unique(df$temperature))[i]
      df.temp = df %>% filter(temperature == temp)
      params.temp = params.grid |> 
        filter(temperature == temp)
      Twater = unique(as.numeric(as.character(df.temp$temperature)))
      T1 = Twater + 273.15
      
      for(a in c(1:nrow(params.temp))){
        cT = params.temp$cT[a]
        pT = params.temp$pT[a]
        
        df.temp$predicted = (exp((pT/cT)*(1 - cT*df.temp$duration - exp(-cT*df.temp$duration))))*100
        error = (df.temp$predicted - df.temp$survival)
        rmse <- sqrt(mean((error) ^ 2))
        nrmse <- rmse/sd(df.temp$survival)
        params.temp$nrmse[a] = nrmse
      }
      params = rbind(params, params.temp)
    }
    params$jackVar = group1
    params_list[[n]] = params
  }
  parameters = bind_rows(params_list)
  write.csv(parameters, paste0("outputs/gridsearch_", label, ".csv"), row.names = F)
}
gridwinfun = function(dat, label){
  jack_df = data.frame(jackVar = unique(dat$jackVar))
  winning_list = list()
  winning_pTlist = list()
  for(n in 1:nrow(jack_df)){
    group1 = jack_df$jackVar[n]
    
    params.temp = dat |> 
      filter(jackVar == group1)
    
    winning.pT = data.frame()
    winning.grid = data.frame()
    params.init = data.frame()
    for(i in c(45, 50, 55, 60, 65)){
      params.temp1 = params.temp |> 
        filter(temperature == i)
      min.limit = quantile(params.temp1$nrmse, probs = 0.01)
      params.temp1 = params.temp1 |> 
        filter(nrmse <= min.limit)
      pT1 = min(params.temp1$pT)
      params.temp1 = params.temp1 |> 
        filter(pT == pT1)
      params.temp1 = params.temp1 |>
        filter(cT == min(params.temp1$cT))
      winning.pT = rbind(winning.pT, params.temp1)
      params.init = rbind(params.init, params.temp1)
    }
    
    # use log-linear model to generate confidence limits
    mod1 = lm(data = params.init, log(cT) ~ temperature)
    cTs = as.data.frame(exp(predict(mod1, newdata = data.frame(temperature = seq(45, 65, by = 5)), interval = "confidence", level = 0.96))) ## alpha increased until confidence limits incorporate grid search value selection
    mod2 = lm(data = params.init, log(pT) ~ temperature)
    pTs = as.data.frame(exp(predict(mod2, newdata = data.frame(temperature = seq(45, 65, by = 5)), interval = "confidence", level = 0.96)))
    
    # select confidence limits to act as boundaries for optimization
    winning_list[[n]] = data.frame(jackVar = group1,
                                   temperature = seq(45, 65, by = 5),
                                   pT = params.init$pT,
                                   pT.min = pTs$lwr,
                                   pT.max = pTs$upr,
                                   cT = params.init$cT,
                                   cT.min = cTs$lwr,
                                   cT.max = cTs$upr,
                                   nrmse = params.init$nrmse)
    winning_pTlist[[n]] = winning.pT
  }
  winning.grid = bind_rows(winning_list)
  write.csv(winning.grid, paste0("outputs/gridsearch_best_min_", label, ".csv"), row.names = F)
  winning.pT = bind_rows(winning_pTlist)
  write.csv(winning.pT, paste0("outputs/gridsearch_best_", label, ".csv"), row.names = F)
}

optimfun = function(dat, observed, label){
  jack_df = data.frame(jackVar = unique(dat$jackVar))
  params_list = list()
  for(n in 1:nrow(jack_df)){
    group1 = jack_df$jackVar[n]
    print(group1)
    
    df = observed |> 
      filter(jackVar == group1)
    
    params.grid = dat |> 
      filter(jackVar == group1)
    
    init = c(params.grid$cT, 
             params.grid$pT)
    objective_function <- function(dat, pars){
      temperature = c(45, 50, 55, 60, 65)
      cTs = pars[1:5]
      pTs = pars[6:10]
      
      cT.pT.ss = list()
      curves.ss = list()
      for(a in c(1:5)){
        group.temp = df |>
          filter(temperature == temperature[a])
        cT.temp = cTs[a]
        pT.temp = pTs[a]
        group.temp$prediction = exp((pT.temp/cT.temp)*(1 - cT.temp*group.temp$duration - exp(-cT.temp*group.temp$duration)))*100
        error = (group.temp$prediction - group.temp$survival)
        rmse <- sqrt(mean((error)^2))
        nrmse <- rmse/sd(group.temp$survival)
        curves.ss[a] <- nrmse
      }
      nrmse1 = max(unlist(curves.ss))
      
      params.df = data.frame(cT = cTs,
                             pT = pTs,
                             temperature = c(45, 50, 55, 60, 65))
      mod1 = lm(data = params.df, log(cT) ~ temperature)
      error = fitted(mod1) - log(params.df$cT)
      rmse = sqrt(mean(error^2))
      nrmse2 = rmse/sd(log(params.df$cT))
      
      mod2 = lm(data = params.df, log(pT) ~ temperature)
      error = fitted(mod2) - log(params.df$pT)
      rmse = sqrt(mean(error^2))
      nrmse3 = rmse/sd(log(params.df$pT))
      
      stat.final = max(c(nrmse1, nrmse2, nrmse3))
      return(stat.final)
    }
    
    tryCatch({
      set.seed(12)
      optim_sa_result <- optim_sa(
        fun = function(pars) objective_function(dat = df, pars),
        start = init,
        maximization = FALSE,
        trace = TRUE,
        lower = c(params.grid$cT.min, params.grid$pT.min),
        upper = c(params.grid$cT.max, params.grid$pT.max),
        control = list(nlimit = 1000, stopac = 1)
      )
      params_list[[n]] = data.frame(jackVar = group1,
                                    groupVar = unique(df$groupVar),
                                    accuracy = optim_sa_result$function_value,
                                    cT = optim_sa_result$par[1:5],
                                    cT.init = init[1:5],
                                    pT = optim_sa_result$par[6:10],
                                    pT.init = init[6:10],
                                    temperature = c(45, 50, 55, 60, 65))
    }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    rm(optim_sa_result, init)
  }
  parameters = bind_rows(params_list)
  saveRDS(parameters, paste0("outputs/optimization_", label, ".rds"))
}

predictfun = function(jack, observed, label){
  df = observed
  df$modVar  = paste0(df$groupVar, "_", df$temperature)
  df = df |> 
    left_join(results_df)
  jack$temperature  = as.numeric(as.character(jack$temperature))
  jack$modVar  = paste0(jack$groupVar, "_", jack$temperature)
  jack.smooth = list()
  for(i in unique(df$modVar)){
    params.temp = jack[jack$modVar == i,]
    df.temp = df[df$modVar == i,]
    newdata = data.frame(duration = seq(0, max(df.temp$duration), by = 1))
    preds.temp = list()
    for(a in unique(params.temp$jackVar)){
      params.temp1 = params.temp |> 
        filter(jackVar == a)
      pT = unique(params.temp1$pT)
      cT = unique(params.temp1$cT)
      preds.temp.a = data.frame(groupVar = unique(params.temp1$groupVar),
                                   modVar = i,
                                   jackVar = a,
                                   species = unique(df.temp$species),
                                   season = unique(df.temp$season),
                                   age = unique(df.temp$age),
                                   temperature = unique(params.temp1$temperature),
                                   duration = newdata$duration,
                                   cT = cT,
                                   pT = pT,
                                   predicted = (exp((pT/cT)*(1 - cT*newdata$duration - exp(-cT*newdata$duration))))*100)
      preds.stat = preds.temp.a |> filter(preds.temp.a$duration %in% df.temp$duration) |> 
        inner_join(df.temp, by = c("groupVar", "modVar", "species", "season", "age", "temperature", "duration")) |> 
        filter(duration != 0)
      preds.stat$error = preds.stat$predicted - preds.stat$survival
      preds.temp.a$std.err = sd(preds.stat$error)
      preds.temp[[a]] = preds.temp.a
    }
    preds.temp = bind_rows(preds.temp)
    
    cT.lim = 1.96*std.error(log(preds.temp$cT))
    cT = mean(log(preds.temp$cT))
    pT.lim = 1.96*std.error(log(preds.temp$pT))
    pT = mean(log(preds.temp$pT))
    
    jack.temp = list()
    for(a in unique(preds.temp$duration)){
      params.temp1 = preds.temp |> 
        filter(duration == a)
      con.lim = 1.96*params.temp1$std.err
      pred.mean = mean(params.temp1$predicted)
      jack.temp[[a + 1]] = data.frame(groupVar = unique(params.temp1$groupVar),
                                      species = unique(params.temp1$species),
                                      season = unique(params.temp1$season),
                                      age = unique(params.temp1$age),
                                      modVar = i,
                                      temperature = unique(params.temp1$temperature),
                                      duration = a,
                                      cT.min = exp(cT - cT.lim),
                                      cT = exp(cT),
                                      cT.max = exp(cT + cT.lim),
                                      log.cT.se = std.error(log(preds.temp$cT)),
                                      pT.min = exp(pT - pT.lim),
                                      pT = exp(pT),
                                      pT.max = exp(pT + pT.lim),
                                      log.pT.se = std.error(log(preds.temp$pT)),
                                      predicted.min = pred.mean - con.lim,
                                      predicted = pred.mean,
                                      predicted.max = pred.mean + con.lim,
                                      predicted.se = params.temp1$std.err)
    }
    jack.smooth[[i]] = bind_rows(jack.temp)
  }
  jack.smooth = bind_rows(jack.smooth)
  saveRDS(jack.smooth, paste0("outputs/predictions.smoothed_", label, ".rds"))
  return(jack.smooth)
}

survivalplotsfinalfun = function(dat, jack){
  dat$groupVar = as.factor(dat$groupVar)
  levels(dat$groupVar) = c("PIPO_allothers", "PIPO_allothers", "PIPO_spring_n", "PIPO_allothers",
                                    "PSME_allautumn", "PSME_allautumn", "PSME_spring_n", "PSME_spring_o")
  dat = dat |> 
    left_join(results_fin)
  
  dfs = list(dat, jack)
  for(i in 1:2){
    df.temp = dfs[[i]]
    df.temp$temperature = as.factor(df.temp$temperature)
    df.temp$timestamp = as.POSIXct(seconds.to.hms(df.temp$duration), format = "%H:%M:%S")
    df.temp$season = factor(df.temp$season, levels = c("spring", "autumn", "spring and autumn"))
    df.temp$age = factor(df.temp$age, levels = c("n", "o", "n and o"))
    df.temp$plotVar = factor(paste0(df.temp$season, "_", df.temp$age), levels = c("spring_n",
                                                                                  "spring_o",
                                                                                  "autumn_n and o",
                                                                                  "spring and autumn_n and o"))
    dfs[[i]] = df.temp
  }
  dat = dfs[[1]]
  jack = dfs[[2]]
  
  gplots = list()
  for(n in 1:5){
    t = seq(45, 65, by = 5)[n]
    g1 =
      ggplot() +
      geom_ribbon(data = jack |> filter(temperature == t), aes(x = timestamp, ymin = predicted.min, ymax = predicted.max, fill = plotVar), alpha = 0.2) +
      geom_line(data = jack |> filter(temperature == t), aes(x = timestamp, y = predicted, col = plotVar)) +
      geom_point(data = dat |> filter(temperature == t), aes(x = timestamp, y = survival, col = plotVar, shape = plotVar)) +
      scale_x_datetime(date_labels = "%H:%M:%S") +
      xlab("") +
      coord_cartesian(ylim = c(0, 100)) +
      scale_shape_manual(name = "", labels = c("Spring - new",
                                               "Spring - 1yo",
                                               "All autumn",
                                               "Spring 1yo and all autumn"),
                         values = c(1, 19, 3, 3)) +
      scale_color_manual(name = "", labels = c("Spring - new",
                                               "Spring - 1yo",
                                               "All autumn",
                                               "Spring 1yo and all autumn"),
                         values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
      scale_fill_manual(name = "", labels = c("Spring - new",
                                              "Spring - 1yo",
                                              "All autumn",
                                              "Spring 1yo and all autumn"),
                        values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
      facet_wrap(facets = "species", nrow = 2) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
      ggtitle(paste0(t, " C"))
    if(t == 45){
      g1 = g1 +
        ylab("Survival (%)") +
        theme(legend.position = "none",
              strip.background = element_blank(),
              strip.text.x = element_blank())
    } else if(t > 45 & t < 65){
      g1 = g1 +
        ylab("") +
        theme(axis.text.y=element_blank(), 
              axis.ticks.y=element_blank(),
              legend.position = "none",
              strip.background = element_blank(),
              strip.text.x = element_blank())
    } else {
      g1 = g1 +
        ylab("") +
        theme(axis.text.y=element_blank(), 
              axis.ticks.y=element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_blank())
    }
    gplots[[n]] = g1
  }
  ggsave(paste0("outputs/survivalmodels_fit_smooth_final.jpg"), gplots[[1]]|gplots[[2]]|gplots[[3]]|gplots[[4]]|gplots[[5]], device = "jpeg", height = 4, width = 10)
  
  # g1 =
    ggplot() +
    geom_ribbon(data = jack, aes(x = timestamp, ymin = predicted.min, ymax = predicted.max, fill = plotVar), alpha = 0.2) +
    geom_line(data = jack, aes(x = timestamp, y = predicted, col = plotVar)) +
    geom_point(data = dat, aes(x = timestamp, y = survival, col = plotVar, shape = plotVar)) +
    scale_x_datetime(date_labels = "%H:%M:%S") +
    xlab("") +
    coord_cartesian(ylim = c(0, 100)) +
    scale_shape_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "All autumn",
                                             "Spring 1yo and all autumn"),
                       values = c(1, 19, 3, 3)) +
    scale_color_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "All autumn",
                                             "Spring 1yo and all autumn"),
                       values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
    scale_fill_manual(name = "", labels = c("Spring - new",
                                            "Spring - 1yo",
                                            "All autumn",
                                            "Spring 1yo and all autumn"),
                      values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
    facet_wrap(facets = c("species", "season", "age"), nrow = 2) +
    theme_bw()
  ggsave(paste0("outputs/survivalmodels_fit_smooth_final.jpg"), g1, device = "jpeg", height = 4, width = 10)
}

loglinearfun = function(dat, errors, label){
  cT_results_df = data.frame()
  pT_results_df = data.frame()
  for(n in c(1:nrow(results_df))){
    group1 = unique(results_df$groupVar)[n]
    params.temp = dat |> filter(groupVar == group1)
    errors.temp = errors |> filter(groupVar == group1)
    se.nonl = aggregate(data = errors.temp, predicted.se ~ temperature, FUN = mean)[,2]
    
    mod.temp = lm(data = params.temp, log(cT) ~ temperature)
    se.l = summary(mod.temp)$sigma
    
    results.temp = data.frame(groupVar = group1,
                              temperature = seq(45,65,by=5),
                              intercept = summary(mod.temp)$coefficients[1,1],
                              estimate = summary(mod.temp)$coefficients[2,1]*seq(45,65,by=5) + summary(mod.temp)$coefficients[1,1],
                              delta1 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2) + ((summary(mod.temp)$coefficients[2,1])^2) * (se.l)^2),
                              delta2 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2)),
                              std.err = summary(mod.temp)$coefficients[2,2],
                              AdjR2 = summary(mod.temp)$adj.r.squared,
                              f.stat = summary(mod.temp)$f[1])
    results.temp$est.min1 = results.temp$estimate - results.temp$delta1*1.96
    results.temp$est.max1 = results.temp$estimate + results.temp$delta1*1.96
    results.temp$est.min2 = results.temp$estimate - results.temp$delta2*1.96
    results.temp$est.max2 = results.temp$estimate + results.temp$delta2*1.96
    cT_results_df = rbind(cT_results_df, results.temp)
    
    mod.temp = lm(data = params.temp, log(pT) ~ temperature)
    
    results.temp = data.frame(groupVar = group1,
                              temperature = seq(45,65,by=5),
                              intercept = summary(mod.temp)$coefficients[1,1],
                              estimate = summary(mod.temp)$coefficients[2,1]*seq(45,65,by=5) + summary(mod.temp)$coefficients[1,1],
                              delta1 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2) + ((summary(mod.temp)$coefficients[2,1])^2) * (se.l)^2),
                              delta2 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2)),
                              std.err = summary(mod.temp)$coefficients[2,2],
                              AdjR2 = summary(mod.temp)$adj.r.squared,
                              f.stat = summary(mod.temp)$f[1])
    results.temp$est.min1 = results.temp$estimate - results.temp$delta1*1.96
    results.temp$est.max1 = results.temp$estimate + results.temp$delta1*1.96
    results.temp$est.min2 = results.temp$estimate - results.temp$delta2*1.96
    results.temp$est.max2 = results.temp$estimate + results.temp$delta2*1.96
    
    pT_results_df = rbind(pT_results_df, results.temp)
  }
  write.csv(cT_results_df, paste0("outputs/cT.loglinearmodelresults_", label, ".csv"), row.names = F)
  write.csv(pT_results_df, paste0("outputs/pT.loglinearmodelresults_", label, ".csv"), row.names = F)
  
  compare.grid1 = data.frame(par1 = combn(1:4, 2)[1,],
                             par2 = combn(1:4, 2)[2,])
  compare.grid2 = data.frame(par1 = combn(1:4, 2)[1,] + max(compare.grid1$par2),
                             par2 = combn(1:4, 2)[2,] + max(compare.grid1$par2))
  compare.grid = rbind(compare.grid1, compare.grid2)
  compare.grid$par1 = unique(results_df$groupVar)[compare.grid$par1]
  compare.grid$par2 = unique(results_df$groupVar)[compare.grid$par2]
  results.grid = data.frame()
  for(i in 1:nrow(compare.grid)){
    cT.temp1 = cT_results_df |> 
      filter(groupVar == compare.grid$par1[i])
    cT.temp2 = cT_results_df |> 
      filter(groupVar == compare.grid$par2[i])
    pT.temp1 = pT_results_df |> 
      filter(groupVar == compare.grid$par1[i])
    pT.temp2 = pT_results_df |> 
      filter(groupVar == compare.grid$par2[i])
    
    results.temp = data.frame(compare1 = compare.grid$par1[i],
                              compare2 = compare.grid$par2[i],
                              temperature = seq(45, 65, by = 5),
                              pT.overlap1 = ifelse(pT.temp1$est.max1 >= pT.temp2$est.min1 & pT.temp2$est.max1 >= pT.temp1$est.min1, TRUE, FALSE),
                              cT.overlap1 = ifelse(cT.temp1$est.max1 >= cT.temp2$est.min1 & cT.temp2$est.max1 >= cT.temp1$est.min1, TRUE, FALSE),
                              pT.overlap2 = ifelse(pT.temp1$est.max2 >= pT.temp2$est.min2 & pT.temp2$est.max2 >= pT.temp1$est.min2, TRUE, FALSE),
                              cT.overlap2 = ifelse(cT.temp1$est.max2 >= cT.temp2$est.min2 & cT.temp2$est.max2 >= cT.temp1$est.min2, TRUE, FALSE))
    results.grid = rbind(results.grid, results.temp)
  }
  write.csv(results.grid, paste0("outputs/significantdiff_parameters_", label, ".csv"), row.names = F)
  
  results.grid$temperature = paste0("pT_", results.grid$temperature)
  results.gridpT1 = results.grid |> 
    dplyr::select(compare1, compare2, pT.overlap1, temperature) |> 
    pivot_wider(names_from = temperature, values_from = pT.overlap1)
  results.gridpT2 = results.grid |> 
    dplyr::select(compare1, compare2, pT.overlap2, temperature) |> 
    pivot_wider(names_from = temperature, values_from = pT.overlap2)
  results.grid$temperature = paste0("cT_", substr(results.grid$temperature, 4, 5))
  results.gridcT1 = results.grid |> 
    dplyr::select(compare1, compare2, cT.overlap1, temperature) |> 
    pivot_wider(names_from = temperature, values_from = cT.overlap1)
  results.gridcT2 = results.grid |> 
    dplyr::select(compare1, compare2, cT.overlap2, temperature) |> 
    pivot_wider(names_from = temperature, values_from = cT.overlap2)
  results.gridall = cbind(results.gridpT1, results.gridpT2, results.gridcT1, results.gridcT2)
  
  write.csv(results.gridall, paste0("outputs/significantdiffTF_parameters_", label, ".csv"), row.names = F)
}
loglinearsecfun = function(dat, errors, label){
  cT_results_df = data.frame()
  pT_results_df = data.frame()
  for(n in c(1:nrow(results_df))){
    group1 = unique(results_df$groupVar)[n]
    params.temp = dat |> filter(groupVar == group1)
    errors.temp = errors |> filter(groupVar == group1)
    se.nonl = aggregate(data = errors.temp, predicted.se ~ temperature, FUN = mean)[,2]
    
    mod.temp = lm(data = params.temp, log(cT) ~ temperature)
    se.l = summary(mod.temp)$sigma
    
    results.temp = data.frame(groupVar = group1,
                              temperature = seq(45,65,by=5),
                              intercept = summary(mod.temp)$coefficients[1,1],
                              estimate = summary(mod.temp)$coefficients[2,1]*seq(45,65,by=5) + summary(mod.temp)$coefficients[1,1],
                              delta1 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2) + ((summary(mod.temp)$coefficients[2,1])^2) * (se.l)^2),
                              delta2 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2)),
                              std.err = summary(mod.temp)$coefficients[2,2],
                              AdjR2 = summary(mod.temp)$adj.r.squared,
                              f.stat = summary(mod.temp)$f[1])
    results.temp$est.min1 = results.temp$estimate - results.temp$delta1*1.96
    results.temp$est.max1 = results.temp$estimate + results.temp$delta1*1.96
    results.temp$est.min2 = results.temp$estimate - results.temp$delta2*1.96
    results.temp$est.max2 = results.temp$estimate + results.temp$delta2*1.96
    cT_results_df = rbind(cT_results_df, results.temp)
    
    mod.temp = lm(data = params.temp, log(pT) ~ temperature)
    
    results.temp = data.frame(groupVar = group1,
                              temperature = seq(45,65,by=5),
                              intercept = summary(mod.temp)$coefficients[1,1],
                              estimate = summary(mod.temp)$coefficients[2,1]*seq(45,65,by=5) + summary(mod.temp)$coefficients[1,1],
                              delta1 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2) + ((summary(mod.temp)$coefficients[2,1])^2) * (se.l)^2),
                              delta2 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2)),
                              std.err = summary(mod.temp)$coefficients[2,2],
                              AdjR2 = summary(mod.temp)$adj.r.squared,
                              f.stat = summary(mod.temp)$f[1])
    results.temp$est.min1 = results.temp$estimate - results.temp$delta1*1.96
    results.temp$est.max1 = results.temp$estimate + results.temp$delta1*1.96
    results.temp$est.min2 = results.temp$estimate - results.temp$delta2*1.96
    results.temp$est.max2 = results.temp$estimate + results.temp$delta2*1.96
    
    pT_results_df = rbind(pT_results_df, results.temp)
  }
  write.csv(cT_results_df, paste0("outputs/cT.loglinearmodelresults_", label, ".csv"), row.names = F)
  write.csv(pT_results_df, paste0("outputs/pT.loglinearmodelresults_", label, ".csv"), row.names = F)
  
  compare.grid1 = data.frame(par1 = combn(1:3, 2)[1,],
                             par2 = combn(1:3, 2)[2,])
  compare.grid2 = data.frame(par1 = combn(1:3, 2)[1,] + max(compare.grid1$par2),
                             par2 = combn(1:3, 2)[2,] + max(compare.grid1$par2))
  compare.grid = rbind(compare.grid1, compare.grid2)
  compare.grid$par1 = unique(results_df$groupVar)[compare.grid$par1]
  compare.grid$par2 = unique(results_df$groupVar)[compare.grid$par2]
  results.grid = data.frame()
  for(i in 1:nrow(compare.grid)){
    cT.temp1 = cT_results_df |> 
      filter(groupVar == compare.grid$par1[i])
    cT.temp2 = cT_results_df |> 
      filter(groupVar == compare.grid$par2[i])
    pT.temp1 = pT_results_df |> 
      filter(groupVar == compare.grid$par1[i])
    pT.temp2 = pT_results_df |> 
      filter(groupVar == compare.grid$par2[i])
    
    results.temp = data.frame(compare1 = compare.grid$par1[i],
                              compare2 = compare.grid$par2[i],
                              temperature = seq(45, 65, by = 5),
                              pT.overlap1 = ifelse(pT.temp1$est.max1 >= pT.temp2$est.min1 & pT.temp2$est.max1 >= pT.temp1$est.min1, TRUE, FALSE),
                              cT.overlap1 = ifelse(cT.temp1$est.max1 >= cT.temp2$est.min1 & cT.temp2$est.max1 >= cT.temp1$est.min1, TRUE, FALSE),
                              pT.overlap2 = ifelse(pT.temp1$est.max2 >= pT.temp2$est.min2 & pT.temp2$est.max2 >= pT.temp1$est.min2, TRUE, FALSE),
                              cT.overlap2 = ifelse(cT.temp1$est.max2 >= cT.temp2$est.min2 & cT.temp2$est.max2 >= cT.temp1$est.min2, TRUE, FALSE))
    results.grid = rbind(results.grid, results.temp)
  }
  write.csv(results.grid, paste0("outputs/significantdiff_parameters_", label, ".csv"), row.names = F)
  
  results.grid$temperature = paste0("pT_", results.grid$temperature)
  results.gridpT1 = results.grid |> 
    dplyr::select(compare1, compare2, pT.overlap1, temperature) |> 
    pivot_wider(names_from = temperature, values_from = pT.overlap1)
  results.gridpT2 = results.grid |> 
    dplyr::select(compare1, compare2, pT.overlap2, temperature) |> 
    pivot_wider(names_from = temperature, values_from = pT.overlap2)
  results.grid$temperature = paste0("cT_", substr(results.grid$temperature, 4, 5))
  results.gridcT1 = results.grid |> 
    dplyr::select(compare1, compare2, cT.overlap1, temperature) |> 
    pivot_wider(names_from = temperature, values_from = cT.overlap1)
  results.gridcT2 = results.grid |> 
    dplyr::select(compare1, compare2, cT.overlap2, temperature) |> 
    pivot_wider(names_from = temperature, values_from = cT.overlap2)
  results.gridall = cbind(results.gridpT1, results.gridpT2, results.gridcT1, results.gridcT2)
  
  write.csv(results.gridall, paste0("outputs/significantdiffTF_parameters_", label, ".csv"), row.names = F)
}
loglinearfinalfun = function(dat, errors, label){
  cT_results_fin = data.frame()
  pT_results_fin = data.frame()
  for(n in c(1:nrow(results_fin))){
    group1 = unique(results_fin$groupVar)[n]
    params.temp = dat |> filter(groupVar == group1)
    errors.temp = errors |> filter(groupVar == group1)
    se.nonl = aggregate(data = errors.temp, predicted.se ~ temperature, FUN = mean)[,2]
    
    mod.temp = lm(data = params.temp, log(cT) ~ temperature)
    se.l = summary(mod.temp)$sigma
    
    # se.delta sqrt(
    #   d(nonlinear)^2*log(se.nonl)^2 + # the slope of the nonlinear function at the value of x
    #         ## deriv(nonlinear) is substituted out for 1/T based on the delta method
    #         ## se.nonl = sd(x)/sqrt(n)
    #         ## so, se.nonl^2 = sd(x)^2/n
    #   d(linear)^2*(se.l)^2)
    #         ## d(linear) = slope(linear) = estimate
    
    results.temp = data.frame(groupVar = group1,
                              temperature = seq(45,65,by=5),
                              intercept = summary(mod.temp)$coefficients[1,1],
                              estimate = summary(mod.temp)$coefficients[2,1]*seq(45,65,by=5) + summary(mod.temp)$coefficients[1,1],
                              delta1 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2) + ((summary(mod.temp)$coefficients[2,1])^2) * (se.l)^2),
                              delta2 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2)),
                              std.err = summary(mod.temp)$coefficients[2,2],
                              AdjR2 = summary(mod.temp)$adj.r.squared,
                              f.stat = summary(mod.temp)$f[1])
    results.temp$est.min1 = results.temp$estimate - results.temp$delta1*1.96
    results.temp$est.max1 = results.temp$estimate + results.temp$delta1*1.96
    results.temp$est.min2 = results.temp$estimate - results.temp$delta2*1.96
    results.temp$est.max2 = results.temp$estimate + results.temp$delta2*1.96
    cT_results_fin = rbind(cT_results_fin, results.temp)
    
    mod.temp = lm(data = params.temp, log(pT) ~ temperature)
    
    results.temp = data.frame(groupVar = group1,
                              temperature = seq(45,65,by=5),
                              intercept = summary(mod.temp)$coefficients[1,1],
                              estimate = summary(mod.temp)$coefficients[2,1]*seq(45,65,by=5) + summary(mod.temp)$coefficients[1,1],
                              delta1 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2) + ((summary(mod.temp)$coefficients[2,1])^2) * (se.l)^2),
                              delta2 = sqrt((1/seq(45,65,by=5))^2 * ((log(se.nonl))^2)),
                              std.err = summary(mod.temp)$coefficients[2,2],
                              AdjR2 = summary(mod.temp)$adj.r.squared,
                              f.stat = summary(mod.temp)$f[1])
    results.temp$est.min1 = results.temp$estimate - results.temp$delta1*1.96
    results.temp$est.max1 = results.temp$estimate + results.temp$delta1*1.96
    results.temp$est.min2 = results.temp$estimate - results.temp$delta2*1.96
    results.temp$est.max2 = results.temp$estimate + results.temp$delta2*1.96
    
    pT_results_fin = rbind(pT_results_fin, results.temp)
  }
  write.csv(cT_results_fin, paste0("outputs/cT.loglinearmodelresults_", label, ".csv"), row.names = F)
  write.csv(pT_results_fin, paste0("outputs/pT.loglinearmodelresults_", label, ".csv"), row.names = F)
  
  compare.grid = data.frame(par1 = combn(1:length(unique(results_fin$groupVar)), 2)[1,],
                            par2 = combn(1:length(unique(results_fin$groupVar)), 2)[2,])
  compare.grid$par1 = unique(results_fin$groupVar)[compare.grid$par1]
  compare.grid$par2 = unique(results_fin$groupVar)[compare.grid$par2]
  results.grid = data.frame()
  for(i in 1:nrow(compare.grid)){
    cT.temp1 = cT_results_fin |> 
      filter(groupVar == compare.grid$par1[i])
    cT.temp2 = cT_results_fin |> 
      filter(groupVar == compare.grid$par2[i])
    pT.temp1 = pT_results_fin |> 
      filter(groupVar == compare.grid$par1[i])
    pT.temp2 = pT_results_fin |> 
      filter(groupVar == compare.grid$par2[i])
    
    results.temp = data.frame(compare1 = compare.grid$par1[i],
                              compare2 = compare.grid$par2[i],
                              temperature = seq(45, 65, by = 5),
                              pT.overlap1 = ifelse(pT.temp1$est.max1 >= pT.temp2$est.min1 & pT.temp2$est.max1 >= pT.temp1$est.min1, TRUE, FALSE),
                              cT.overlap1 = ifelse(cT.temp1$est.max1 >= cT.temp2$est.min1 & cT.temp2$est.max1 >= cT.temp1$est.min1, TRUE, FALSE),
                              pT.overlap2 = ifelse(pT.temp1$est.max2 >= pT.temp2$est.min2 & pT.temp2$est.max2 >= pT.temp1$est.min2, TRUE, FALSE),
                              cT.overlap2 = ifelse(cT.temp1$est.max2 >= cT.temp2$est.min2 & cT.temp2$est.max2 >= cT.temp1$est.min2, TRUE, FALSE))
    results.grid = rbind(results.grid, results.temp)
  }
  write.csv(results.grid, paste0("outputs/significantdiff_parameters_", label, ".csv"), row.names = F)
  
  results.grid$temperature = paste0("pT_", results.grid$temperature)
  results.gridpT1 = results.grid |> 
    dplyr::select(compare1, compare2, pT.overlap1, temperature) |> 
    pivot_wider(names_from = temperature, values_from = pT.overlap1)
  results.gridpT2 = results.grid |> 
    dplyr::select(compare1, compare2, pT.overlap2, temperature) |> 
    pivot_wider(names_from = temperature, values_from = pT.overlap2)
  results.grid$temperature = paste0("cT_", substr(results.grid$temperature, 4, 5))
  results.gridcT1 = results.grid |> 
    dplyr::select(compare1, compare2, cT.overlap1, temperature) |> 
    pivot_wider(names_from = temperature, values_from = cT.overlap1)
  results.gridcT2 = results.grid |> 
    dplyr::select(compare1, compare2, cT.overlap2, temperature) |> 
    pivot_wider(names_from = temperature, values_from = cT.overlap2)
  results.gridall = cbind(results.gridpT1, results.gridpT2, results.gridcT1, results.gridcT2)
  
  write.csv(results.gridall, paste0("outputs/significantdiffTF_parameters_", label, ".csv"), row.names = F)
}

lplotsinitfun = function(cT, pT, parameters, label){
  dfs = list(cT, pT, parameters)
  for(i in 1:3){
    df.temp = dfs[[i]]
    df.temp$groupVar = as.factor(df.temp$groupVar)
    df.temp = df.temp |> 
      left_join(results_df)
    df.temp$season = factor(df.temp$season, levels = c("spring", "autumn"))
    df.temp$age = factor(df.temp$age, levels = c("n", "o"))
    df.temp$plotVar = factor(paste0(df.temp$season, "_", df.temp$age), levels = c("spring_n",
                                                                                  "spring_o",
                                                                                  "autumn_n",
                                                                                  "autumn_o"))
    if(i %in% 1:2){
      df.temp$estimate = exp(df.temp$estimate)
      df.temp$est.min1 = exp(df.temp$est.min1)
      df.temp$est.max1 = exp(df.temp$est.max1)
      df.temp$est.min2 = exp(df.temp$est.min2)
      df.temp$est.max2 = exp(df.temp$est.max2)
    }
    dfs[[i]] = df.temp
  }
  cT = dfs[[1]]
  pT = dfs[[2]]
  parameters = dfs[[3]]
  
  g1 = ggplot() +
    geom_point(data = parameters, aes(y = cT, x = temperature, col = plotVar, shape = plotVar)) +
    geom_ribbon(data = cT, aes(x = temperature, ymin = est.min1, ymax = est.max1, fill = plotVar), alpha = 0.2) +
    geom_smooth(data = cT, aes(x = temperature, y = estimate, col = plotVar), method = "lm") +
    scale_y_continuous(trans = "log", label = function(x) formatC(x, format = "e", digits = 2)) +
    ylab(expression(italic("c"))) +
    xlab("Temperature") +
    scale_shape_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "Autumn - new",
                                             "Autumn - 1yo"),
                       values = c(1, 19, 1, 19)) +
    scale_color_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "Autumn - new",
                                             "Autumn - 1yo"),
                       values = c("#78B6D2", "#1E79B3", "#FDBF6F", "#FF7F00")) +
    scale_fill_manual(name = "", labels = c("Spring - new",
                                            "Spring - 1yo",
                                            "Autumn - new",
                                            "Autumn - 1yo"),
                      values = c("#78B6D2", "#1E79B3", "#FDBF6F", "#FF7F00")) +
    facet_wrap(facets = "species", nrow = 2) +
    theme_bw() +
    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())
  g2 = ggplot(pT, aes(x = temperature, y = estimate)) +
    geom_point(data = parameters, aes(y = pT, col = plotVar, shape = plotVar)) +
    geom_ribbon(aes(ymin = est.min1, ymax = est.max1, fill = plotVar), alpha = 0.2) +
    geom_smooth(method = "lm", aes(col = plotVar)) +
    scale_y_continuous(trans = "log", label = function(x) formatC(x, format = "e", digits = 2)) +
    ylab(expression(italic("p"))) +
    xlab("Temperature") +
    scale_shape_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "Autumn - new",
                                             "Autumn - 1yo"),
                       values = c(1, 19, 1, 19)) +
    scale_color_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "Autumn - new",
                                             "Autumn - 1yo"),
                       values = c("#78B6D2", "#1E79B3", "#FDBF6F", "#FF7F00")) +
    scale_fill_manual(name = "", labels = c("Spring - new",
                                            "Spring - 1yo",
                                            "Autumn - new",
                                            "Autumn - 1yo"),
                      values = c("#78B6D2", "#1E79B3", "#FDBF6F", "#FF7F00")) +
    facet_wrap(facets = "species", nrow = 2) +
    theme_bw() +
    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank())
  ggsave(paste0("outputs/parametermodels_cT_", label, ".jpg"), g1, device = "jpeg", height = 5, width = 4.5)
  ggsave(paste0("outputs/parametermodels_pT_", label, ".jpg"), g2, device = "jpeg", height = 5, width = 3)
}
lplotssecfun = function(cT, pT, parameters, label){
  dfs = list(cT, pT, parameters)
  for(i in 1:3){
    df.temp = dfs[[i]]
    df.temp$groupVar = as.factor(df.temp$groupVar)
    df.temp = df.temp |> 
      left_join(results_df)
    df.temp$season = factor(df.temp$season, levels = c("spring", "autumn", "spring and autumn"))
    df.temp$age = factor(df.temp$age, levels = c("n", "o", "n and o"))
    df.temp$plotVar = factor(paste0(df.temp$season, "_", df.temp$age), levels = c("spring_n",
                                                                                  "spring_o",
                                                                                  "autumn_n",
                                                                                  "autumn_n and o",
                                                                                  "spring and autumn_o"))
    if(i %in% 1:2){
      df.temp$estimate = exp(df.temp$estimate)
      df.temp$est.min1 = exp(df.temp$est.min1)
      df.temp$est.max1 = exp(df.temp$est.max1)
      df.temp$est.min2 = exp(df.temp$est.min2)
      df.temp$est.max2 = exp(df.temp$est.max2)
    }
    dfs[[i]] = df.temp
  }
  cT = dfs[[1]]
  pT = dfs[[2]]
  parameters = dfs[[3]]
  
  g1 = ggplot() +
    geom_point(data = parameters, aes(y = cT, x = temperature, col = plotVar, shape = plotVar)) +
    geom_ribbon(data = cT, aes(x = temperature, ymin = est.min1, ymax = est.max1, fill = plotVar), alpha = 0.2) +
    geom_smooth(data = cT, aes(x = temperature, y = estimate, col = plotVar), method = "lm") +
    scale_y_continuous(trans = "log", label = function(x) formatC(x, format = "e", digits = 2)) +
    ylab(expression(italic("c"))) +
    xlab("Temperature") +
    scale_shape_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "Autumn - new",
                                             "All autumn",
                                             "All 1yo"),
                       values = c(1, 19, 1, 3, 19)) +
    scale_color_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "Autumn - new",
                                             "All autumn",
                                             "All 1yo"),
                       values = c("#78B6D2", "#1E79B3", "#FDBF6F", "#FF7F00", "gray20")) +
    scale_fill_manual(name = "", labels = c("Spring - new",
                                            "Spring - 1yo",
                                            "Autumn - new",
                                            "All autumn",
                                            "All 1yo"),
                      values = c("#78B6D2", "#1E79B3", "#FDBF6F", "#FF7F00", "gray20")) +
    facet_wrap(facets = "species", nrow = 2) +
    theme_bw() +
    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())
  g2 = ggplot(pT, aes(x = temperature, y = estimate)) +
    geom_point(data = parameters, aes(y = pT, col = plotVar, shape = plotVar)) +
    geom_ribbon(aes(ymin = est.min1, ymax = est.max1, fill = plotVar), alpha = 0.2) +
    geom_smooth(method = "lm", aes(col = plotVar)) +
    scale_y_continuous(trans = "log", label = function(x) formatC(x, format = "e", digits = 2)) +
    ylab(expression(italic("p"))) +
    xlab("Temperature") +
    scale_shape_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "Autumn - new",
                                             "All autumn",
                                             "All 1yo"),
                       values = c(1, 19, 1, 3, 19)) +
    scale_color_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "Autumn - new",
                                             "All autumn",
                                             "All 1yo"),
                       values = c("#78B6D2", "#1E79B3", "#FDBF6F", "#FF7F00", "gray20")) +
    scale_fill_manual(name = "", labels = c("Spring - new",
                                            "Spring - 1yo",
                                            "Autumn - new",
                                            "All autumn",
                                            "All 1yo"),
                      values = c("#78B6D2", "#1E79B3", "#FDBF6F", "#FF7F00", "gray20")) +
    facet_wrap(facets = "species", nrow = 2) +
    theme_bw() +
    guides(linetype = guide_legend(override.aes = list(color = "black"))) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank())
  ggsave(paste0("outputs/parametermodels_cT_", label, ".jpg"), g1, device = "jpeg", height = 5, width = 4.5)
  ggsave(paste0("outputs/parametermodels_pT_", label, ".jpg"), g2, device = "jpeg", height = 5, width = 3)
}
lplotsfinalfun = function(cT, pT, parameters, label){
  dfs = list(cT, pT, parameters)
  for(i in 1:3){
    df.temp = dfs[[i]]
    df.temp$groupVar = as.factor(df.temp$groupVar)
    df.temp = df.temp |> 
      left_join(results_fin)
    df.temp$season = factor(df.temp$season, levels = c("spring", "autumn", "spring and autumn"))
    df.temp$age = factor(df.temp$age, levels = c("n", "o", "n and o"))
    df.temp$plotVar = factor(paste0(df.temp$season, "_", df.temp$age), levels = c("spring_n",
                                                                                  "spring_o",
                                                                                  "autumn_n and o",
                                                                                  "spring and autumn_n and o"))
    if(i %in% 1:2){
      df.temp$estimate = exp(df.temp$estimate)
      df.temp$est.min1 = exp(df.temp$est.min1)
      df.temp$est.max1 = exp(df.temp$est.max1)
      df.temp$est.min2 = exp(df.temp$est.min2)
      df.temp$est.max2 = exp(df.temp$est.max2)
    }
    dfs[[i]] = df.temp
  }
  cT = dfs[[1]]
  pT = dfs[[2]]
  parameters = dfs[[3]]
  
  g1 = ggplot() +
    geom_point(data = parameters, aes(y = cT, x = temperature, col = plotVar, shape = plotVar)) +
    geom_ribbon(data = cT, aes(x = temperature, ymin = est.min1, ymax = est.max1, fill = plotVar), alpha = 0.2) +
    geom_smooth(data = cT, aes(x = temperature, y = estimate, col = plotVar), method = "lm") +
    scale_y_continuous(trans = "log", label = function(x) formatC(x, format = "e", digits = 2)) +
    ylab(expression(italic("c"))) +
    xlab("Temperature") +
    scale_shape_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "All autumn",
                                             "Spring 1yo and all autumn"),
                       values = c(1, 19, 3, 3)) +
    scale_color_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "All autumn",
                                             "Spring 1yo and all autumn"),
                       values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
    scale_fill_manual(name = "", labels = c("Spring - new",
                                            "Spring - 1yo",
                                            "All autumn",
                                            "Spring 1yo and all autumn"),
                      values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
    facet_wrap(facets = "species", nrow = 2) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())
  g2 = ggplot(pT, aes(x = temperature, y = estimate)) +
    geom_point(data = parameters, aes(y = pT, col = plotVar, shape = plotVar)) +
    geom_ribbon(aes(ymin = est.min1, ymax = est.max1, fill = plotVar), alpha = 0.2) +
    geom_smooth(method = "lm", aes(col = plotVar)) +
    scale_y_continuous(trans = "log", label = function(x) formatC(x, format = "e", digits = 2)) +
    ylab(expression(italic("p"))) +
    xlab("Temperature") +
    scale_shape_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "All autumn",
                                             "Spring 1yo and all autumn"),
                       values = c(1, 19, 3, 3)) +
    scale_color_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "All autumn",
                                             "Spring 1yo and all autumn"),
                       values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
    scale_fill_manual(name = "", labels = c("Spring - new",
                                            "Spring - 1yo",
                                            "All autumn",
                                            "Spring 1yo and all autumn"),
                      values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
    facet_wrap(facets = "species", nrow = 2) +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank())
  ggsave(paste0("outputs/parametermodels_cT_", label, "_1.jpg"), g1, device = "jpeg", height = 5, width = 5)
  ggsave(paste0("outputs/parametermodels_pT_", label, "_1.jpg"), g2, device = "jpeg", height = 5, width = 3)
}

ld50fun = function(parameters, observed){
  df = observed
  
  parameters$temperature  = as.numeric(as.character(parameters$temperature))
  LD50 = list()
  for(i in unique(df$groupVar)){
    params.temp = parameters |> 
      filter(groupVar == i)
    df.temp = df |> 
      filter(groupVar == i)
    
    mod.temp1 = lm(data = params.temp, log(cT) ~ temperature)
    mod.temp2 = lm(data = params.temp, log(pT) ~ temperature)
    
    cT.slope = summary(mod.temp1)$coefficients[2,1]
    cT.intercept = summary(mod.temp1)$coefficients[1,1]
    pT.slope = summary(mod.temp2)$coefficients[2,1]
    pT.intercept = summary(mod.temp2)$coefficients[1,1]
    
    LD50.temp = list()
    for(dur in seq(1, 90, by = 1)){
      objective_function <- function(temperature, duration){
        pT = exp(pT.slope*temperature + pT.intercept)
        cT = exp(cT.slope*temperature + cT.intercept)
        x = (exp((pT/cT)*(1 - cT*duration - exp(-cT*duration))))*100
        return(abs(x - 50))
      }
      
      set.seed(12)
      optim_sa_result <- optim_sa(
        fun = function(temperature, duration) objective_function(temperature, duration = dur),
        start = 60,
        maximization = FALSE,
        trace = TRUE,
        lower = 45,
        upper = Inf,
        control = list(nlimit = 1000)
      )
      LD50.temp[[dur]] = data.frame(groupVar = i,
                                    duration = dur,
                                    LD50 = optim_sa_result$par,
                                    accuracy = optim_sa_result$function_value)
    }
    LD50[[i]] = bind_rows(LD50.temp)
  }
  LD50 = bind_rows(LD50)
  LD50 = left_join(LD50, results_fin)
  
  LD50$season = factor(LD50$season, levels = c("spring", "autumn", "spring and autumn"))
  LD50$age = factor(LD50$age, levels = c("n", "o", "n and o"))
  LD50$plotVar = factor(paste0(LD50$season, "_", LD50$age), levels = c("spring_n",
                                                                                "spring_o",
                                                                                "autumn_n and o",
                                                                                "spring and autumn_n and o"))
  write.csv(LD50, "outputs/LD50_final.csv", row.names = F)
}
ld50plotsfun = function(LD50){
  LD50$plotVar = factor(LD50$plotVar, levels = c("spring_n",
                                                 "spring_o",
                                                 "autumn_n and o",
                                                 "spring and autumn_n and o"))
  g1 =
    ggplot() +
    geom_ribbon(data = LD50, aes(x = duration, ymin = LD50 - accuracy, ymax = LD50 + accuracy, fill = plotVar), alpha = 0.5) +
    geom_line(data = LD50, aes(x = duration, y = LD50, col = plotVar)) +
    ylab("LD50 (C)") +
    xlab("Duration of heating (s)") +
    scale_fill_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "All autumn",
                                             "Spring 1yo and all autumn"),
                       values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
    scale_color_manual(name = "", labels = c("Spring - new",
                                             "Spring - 1yo",
                                             "All autumn",
                                             "Spring 1yo and all autumn"),
                       values = c("#78B6D2", "#1E79B3", "#FF7F00", "gray20")) +
    facet_wrap(facets = "species", nrow = 2) +
    scale_x_continuous(limits = c(0, 90), breaks = c(0, 30, 60, 90)) +
    facet_wrap(facets = "species", nrow = 2) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())
  ggsave("outputs/LD50_final.jpg", g1, device = "jpeg", height = 5, width = 5)
}