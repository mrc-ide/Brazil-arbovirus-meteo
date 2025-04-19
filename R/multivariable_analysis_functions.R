## multivariable model - functions for forward selection variable choice 

## test collinearity
collinarity_check<-function(data,all_variables, variable_chosen){
  
  m.cor = abs(cor(data[,all_variables],use="complete.obs",method="pearson"))
  var<-variable_chosen
  sub1 <- m.cor[,colnames(m.cor)==var] #subset just the chosen variable
  test_df1 <- data.frame(names = all_variables,
                         values = sub1)
  remove<-test_df1[which(test_df1$values>0.6 | test_df1$values < -0.6),1] # names of the variables colinear not to be included in next round
  
  return(remove)
}


#### run models
run_multivar_analysis <- function(data,GRAPH, variables_to_try, variables_chosen_already = NA, stage = NA,
                                  save_progress=TRUE){
  
  graph<-inla.read.graph(GRAPH)
  
  var_names <- variables_to_try
  
  models<-list()
  list<-list()
  
  for(i in 1:length(var_names)){
    
    if(length(variables_chosen_already)==1){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==2){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +data[,variables_chosen_already[2]] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==3){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +data[,variables_chosen_already[2]]+data[,variables_chosen_already[3]] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==4){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +data[,variables_chosen_already[2]]+
        data[,variables_chosen_already[3]] + data[,variables_chosen_already[4]] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==5){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +data[,variables_chosen_already[2]]+
        data[,variables_chosen_already[3]] + data[,variables_chosen_already[4]] + data[,variables_chosen_already[5]] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==6){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +data[,variables_chosen_already[2]]+
        data[,variables_chosen_already[3]] + data[,variables_chosen_already[4]] + data[,variables_chosen_already[5]] +
        data[,variables_chosen_already[6]] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==7){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +data[,variables_chosen_already[2]]+
        data[,variables_chosen_already[3]] + data[,variables_chosen_already[4]] + data[,variables_chosen_already[5]] +
        data[,variables_chosen_already[6]] + data[,variables_chosen_already[7]] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==8){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +data[,variables_chosen_already[2]]+
        data[,variables_chosen_already[3]] + data[,variables_chosen_already[4]] + data[,variables_chosen_already[5]] +
        data[,variables_chosen_already[6]] + data[,variables_chosen_already[7]] + data[,variables_chosen_already[8]] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==9){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +data[,variables_chosen_already[2]]+
        data[,variables_chosen_already[3]] + data[,variables_chosen_already[4]] + data[,variables_chosen_already[5]] +
        data[,variables_chosen_already[6]] + data[,variables_chosen_already[7]] + data[,variables_chosen_already[8]] +
        data[,variables_chosen_already[9]] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==10){
      formula <- CASES ~ 1 + lag_cases + data[,variables_chosen_already] +data[,variables_chosen_already[2]]+
        data[,variables_chosen_already[3]] + data[,variables_chosen_already[4]] + data[,variables_chosen_already[5]] +
        data[,variables_chosen_already[6]] + data[,variables_chosen_already[7]] + data[,variables_chosen_already[8]] +
        data[,variables_chosen_already[9]] + data[,variables_chosen_already[10]] +
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    
    models[[i]] <- inla(formula = formula,
                        offset=log(pop),
                        data = data, 
                        family="nbinomial",
                        control.family=list(link='log'),
                        control.compute = list(waic = TRUE),
                        control.inla=list(int.strategy = "eb"), 
                        verbose = F,
                        safe=T)
    
    waic<- models[[i]]$waic$waic
    beta<-models[[i]]$summary.fixed[3+stage,]
    
    list[[i]]<-data.frame(var = var_names[i],
                          waic = waic,
                          beta_mean = c(beta$mean),
                          beta_low = c(beta$`0.025quant`),
                          beta_upp = c(beta$`0.975quant`))
    
    
    models[[i]]<-NA
    print(paste0("finished model ", i, " out of ", length(var_names), " in stage ", stage))
    
  }
  
  result<- bind_rows(list)
  
  if(save_progress==TRUE){
    saveRDS(result, paste0("outputs/multivar/temp/",pathogen,"_stage_",stage,"_", i,".RDS")) 
  }
  
  # choose the best fitting
  result |> filter(beta_low < 0 & beta_upp <0 | beta_low > 0 & beta_upp >0)->sig
  sig$var[sig$waic==min(sig$waic)] ->var
  waic<-sig$waic[sig$var==var]
  
  return(c(var,waic))
}


#### run final model using variables from paper
run_final_model <-function(data,graph,pathogen){

## formulas
if(pathogen=="chikv"){
  formula <- CASES ~ 1 + lag_cases + ENSO_anom + spei6 + mean_rain_01 + ALPHA + cons_dryday_3 + min_hum_3 + mir_pop +
    f(WEEK_NO, replicate=STATE_NO, model="rw1", 
      cyclic=T) + 
    f(GEOM, model="bym2", 
      graph=graph, 
      replicate=YEAR_NO)
}
if(pathogen=="zikv"){
  formula <- CASES ~ 1 + lag_cases + ENSO_anom + spei12 + mean_rain_03 + ALPHA + min_min_hum_3 +
    f(WEEK_NO, replicate=STATE_NO, model="rw1", 
      cyclic=T) + 
    f(GEOM, model="bym2", 
      graph=graph, 
      replicate=YEAR_NO)
}
if(pathogen=="denv"){
  formula <- CASES ~ 1 + lag_cases + max_max_temp_02 + ENSO_anom + max_max_hum_02 + 
    spei3  + temp_range_day_03 + 
    ALPHA + min_min_temp_3 + mean_rain_1 +
    f(WEEK_NO, replicate=STATE_NO, model="rw1", 
      cyclic=T) + 
    f(GEOM, model="bym2", 
      graph=graph, 
      replicate=YEAR_NO)
}

## model
model <- inla(formula = formula,
              offset=log(pop),
              data = data, 
              family="nbinomial",
              control.family=list(link='log'),
              control.compute = list(waic = TRUE, config=TRUE),
              control.inla=list(int.strategy = "eb"), 
              verbose = F,
              safe=T)

return(model)
}


#### sampling
sampling_models <- function(model, index, nsamp=100){
  
  samples_x=inla.posterior.sample(nsamp,model)
  f2 <- matrix(NA, length(index), nsamp)
  xx.s <- inla.posterior.sample.eval(function(...) c(Predictor[1:(length(index))]), samples_x)
  for(i in 1:nsamp){
    f2[,i] = (xx.s[, i])
  }
  
  return(f2)
}
