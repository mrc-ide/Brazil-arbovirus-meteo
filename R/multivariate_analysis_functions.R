## multivariate - fitting models to all 3 timeseries together

## prepare multivariate data

prepare_multivariate_data <- function(years="all"){

## read in individual data
for(path in c("denv", "chikv", "zikv")){
  
  df <- readRDS(paste0("data/", path,"_timeseries",".RDS"))
  
  if(years=="all"){
  df <- df %>% filter(YEAR>=2016) %>% data_processing()
  }
  if(years=="single"){ # to run a shorter example
    df <- df %>% filter(YEAR==2016) %>% data_processing()  
  }
  
  if(path=="zikv"){
    df -> zikv
  }
  if(path=="chikv"){
    df -> chikv
  }
  if(path=="denv"){
    df -> denv
  }
 }

## bind - get response variable into suitable formula
tot<-(nrow(denv)+nrow(chikv)+nrow(zikv))
rowd <- nrow(denv)
rowz <- nrow(zikv)
rowc <- nrow(chikv)

df<-list(OBS = matrix(NA,nrow=tot, ncol = 3))
df$OBS[1:rowd, 1] <- denv$CASES
df$OBS[(rowd+1):(rowz+rowd), 2] <- zikv$CASES
df$OBS[(rowd+rowz+1):(tot), 3] <- chikv$CASES


## add random effect and baseline variables
df$pop <- c(denv$pop, zikv$pop, chikv$pop)
df$STATE_NO <- c(denv$STATE_NO, zikv$STATE_NO, chikv$STATE_NO)
df$YEAR <- c(denv$YEAR, zikv$YEAR, chikv$YEAR)

df$YEAR_NO <- c(denv$YEAR_NO-3, zikv$YEAR_NO, chikv$YEAR_NO-1)


df$lag <-matrix(NA,nrow=tot, ncol = 3)
df$lag[1:rowd, 1] <- denv$lag_cases
df$lag[(rowd+1):(rowz+rowd), 2] <- zikv$lag_cases
df$lag[(rowd+rowz+1):(tot), 3] <- chikv$lag_cases

df$intercept <-matrix(NA,nrow=tot, ncol = 3)
df$intercept[1:rowd, 1] <- rep(1, rowd)
df$intercept[(rowd+1):(rowz+rowd), 2] <-  rep(1, rowz)
df$intercept[(rowd+rowz+1):(tot), 3] <- rep(1, rowc)

df$spatial1 <- c(denv$GEOM, rep(NA, rowz), rep(NA, rowc))
df$spatial2 <- c(rep(NA, rowd), zikv$GEOM, rep(NA, rowc))
df$spatial3 <- c(rep(NA, rowd),rep(NA, rowz), chikv$GEOM)

df$temporal1 <- c(denv$WEEK_NO, rep(NA, rowz), rep(NA, rowc))
df$temporal2 <- c(rep(NA, rowd), zikv$WEEK_NO, rep(NA, rowc))
df$temporal3 <- c(rep(NA, rowd),rep(NA, rowz), chikv$WEEK_NO)


## add fixed effect variables
df$ENSO_anom <-c(denv$ENSO_anom, zikv$ENSO_anom, chikv$ENSO_anom) 
df$spei12 <-c(denv$spei12, zikv$spei12, chikv$spei12) 
df$mean_rain_03 <-c(denv$mean_rain_03, zikv$mean_rain_03, chikv$mean_rain_03) 
df$ALPHA <-c(denv$ALPHA, zikv$ALPHA, chikv$ALPHA) 
df$spei6 <-c(denv$spei6, zikv$spei6, chikv$spei6) 
df$max_max_temp_02 <-c(denv$max_max_temp_02, zikv$max_max_temp_02, chikv$max_max_temp_02) 
df$min_min_hum_3 <-c(denv$min_min_hum_3, zikv$min_min_hum_3, chikv$min_min_hum_3) 
df$mean_rain_01 <-c(denv$mean_rain_01, zikv$mean_rain_01, chikv$mean_rain_01) 
df$spei3 <-c(denv$spei3, zikv$spei3, chikv$spei3) 
df$cons_dryday_3 <-c(denv$cons_dryday_3, zikv$cons_dryday_3, chikv$cons_dryday_3) 
df$min_hum_3 <-c(denv$min_hum_3, zikv$min_hum_3, chikv$min_hum_3) 
df$max_max_hum_02 <-c(denv$max_max_hum_02, zikv$max_max_hum_02, chikv$max_max_hum_02) 
df$mir_pop <-c(denv$mir_pop, zikv$mir_pop, chikv$mir_pop) 
df$temp_range_day_03 <-c(denv$temp_range_day_03, zikv$temp_range_day_03, chikv$temp_range_day_03) 
df$min_min_temp_3 <-c(denv$min_min_temp_3, zikv$min_min_temp_3, chikv$min_min_temp_3) 
df$mean_rain_1 <-c(denv$mean_rain_1, zikv$mean_rain_1, chikv$mean_rain_1) 

return(df)
}


## run the final multivariate model
multivariate_final_run <- function(data, GRAPH){

graph <- inla.read.graph(GRAPH)
links<- list(link="log")

formula <-  OBS ~ -1 + intercept + lag   + ENSO_anom + min_hum_3  + spei3 + ALPHA + 
  cons_dryday_3 + mean_rain_03 +
  
  f(spatial1, model = "bym2", graph = graph, replicate=YEAR_NO) + 
  f(spatial2, copy = "spatial1", fixed=F) +
  f(spatial3, copy = "spatial1", fixed=F) +
  
  f(temporal1, model = "rw1", replicate=STATE_NO, cyclic=T) +
  f(temporal2, copy = "temporal1", fixed=F) +
  f(temporal3, copy = "temporal1", fixed=F) 


model <- inla(formula = formula,
              offset=log(pop),
              data = (data), 
              family = rep("nbinomial", 3),
              control.family=list(links, links,links),
              control.compute = list(waic = TRUE, 
                                     config=TRUE),
              safe=T,verbose=F)
return(model)
}


#### sampling
multivariate_sample <- function(model, nsamp=100, years="all"){
  
  ## data
  for(path in c("denv", "chikv", "zikv")){
    
    df <- readRDS(paste0("data/", path,"_timeseries",".RDS"))
    
    if(years=="all"){
      df <- df %>% filter(YEAR>=2016) %>% data_processing()
    }
    if(years=="single"){ # to run a shorter example
      df <- df %>% filter(YEAR==2016) %>% data_processing()  
    }
    
    if(path=="zikv"){
      df -> zikv
    }
    if(path=="chikv"){
      df -> chikv
    }
    if(path=="denv"){
      df -> denv
    }
  }
  samples_x=inla.posterior.sample(nsamp,model)
  
  #xx.d <- inla.posterior.sample.eval(function(...) c(Predictor[1:1034709]), samples)
  #xx.c <- inla.posterior.sample.eval(function(...) c(Predictor[1222584:1538496]), samples)
  #xx.z <- inla.posterior.sample.eval(function(...) c(Predictor[1034710:1222583]), samples)
  xx.d <- inla.posterior.sample.eval(function(...) c(Predictor[1:(nrow(denv))]), samples_x)
  xx.c <- inla.posterior.sample.eval(function(...) c(Predictor[(nrow(zikv)+nrow(denv)+1):(nrow(zikv)+nrow(denv)+nrow(chikv))]), samples_x)
  xx.z <- inla.posterior.sample.eval(function(...) c(Predictor[(nrow(denv)+1):(nrow(zikv)+nrow(denv))]), samples_x)
  
  
  fd <- matrix(NA, nrow(denv), nsamp)
  fc <- matrix(NA, nrow(chikv), nsamp)
  fz <- matrix(NA, nrow(zikv), nsamp)
  
  for(n in 1:nsamp){
    fd[,n] = exp(xx.d[, n])
    fc[,n] = exp(xx.c[, n])
    fz[,n] = exp(xx.z[, n])
  }
  
  return(list(fd,fc,fz))
}


## compare multivariate to univariate models