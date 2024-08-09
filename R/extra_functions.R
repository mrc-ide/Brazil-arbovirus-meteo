## extra functions


data_processing <- function(df){
  df<-as.data.frame(df)
  df$sanitation<-as.numeric(df$sanitation)
  df$illiteracy<-as.numeric(df$illiteracy)
  df$income<-as.numeric(df$income)
  df$spei1[is.infinite(df$spei1)==T]<-NA
  df$spei6[is.infinite(df$spei6)==T]<-NA
  df$spei12[is.infinite(df$spei12)==T]<-NA
  df$spei3[is.infinite(df$spei3)==T]<-NA
  
  # filter out municipalities with no local transmission
  df %>% 
    select(MUNICIP, CASES, YEAR)%>%
    group_by(MUNICIP, YEAR) %>%
    summarise(CASES=sum(CASES)) %>% 
    mutate(bin=ifelse(CASES>=5,1,0), # continuous low level transmission
           bin2=ifelse(CASES>20,1,0)) %>% # atleast one year with many cases
    group_by(MUNICIP) %>%
    summarise(years_n=sum(bin), years_n2=sum(bin2))%>%
    filter(years_n==5|years_n2>=1)-> df2
  df %>% 
    filter(MUNICIP %in% unique(df2$MUNICIP)==T)-> df
  
  return(df)
}

run_baseline <- function(data, GRAPH){
  
  graph<-inla.read.graph(GRAPH)
  
  formula <- CASES ~ 1 + lag_cases + 
    f(WEEK_NO, replicate=STATE_NO, model="rw1", 
      cyclic=T) + 
    f(GEOM, model="bym2", 
      graph=graph, 
      replicate=YEAR_NO)
  
  model <- inla(formula = formula,
                offset=log(pop),
                data = data, 
                family="nbinomial",
                control.family=list(link='log'),
                control.compute = list(waic = TRUE, cpo=T,
                                       config=TRUE, return.marginals.predictor=TRUE),
                verbose = F)
  
  return(model)
  
}


process_univar_names <- function(res){
  
  res %>% mutate(group=var) %>%
    mutate(group=case_when(
      var == "ENSO_anom" ~ "ENSO",
      var == "ALPHA" | var == "DEGREE" |var == "BETWEEN"  ~ "mobility",
      substring(var,1,4) =="spei" ~ "SPEI",
      substring(var,1,5) =="flood" ~ "flooding",
      substring(var,1,8)=="cum_rain"~"cumulative rain",
      substring(var,1,9)=="mean_rain"~"mean rain",
      substring(var,1,8)=="min_temp"~"min temp",
      substring(var,1,8)=="max_temp"~"max temp",
      substring(var,1,9)=="mean_temp"~"mean temp",
      substring(var,1,7)=="min_hum"~"min hum",
      substring(var,1,7)=="max_hum"~"max hum",
      substring(var,1,8)=="mean_hum"~"mean hum",
      substring(var,1,12)=="max_max_temp"~"Absolute maximum temp",
      substring(var,1,12)=="min_min_temp"~"Absolute minimum temp",
      substring(var,1,11)=="max_max_hum"~"Absolute maximum hum",
      substring(var,1,11)=="min_min_hum"~"Absolute minimum hum",
      var=="mir_pop" ~ "MIR",
      substring(var,1,8)=="heavy_raindays"~"heavy raindays",
      substring(var,1,16)=="v_heavy_raindays"~"very heavy raindays",
      substring(var,1,9)=="cons_rain"~"consecutive rain",
      substring(var,1,8)=="cons_dry"~"consecutive dry",
      substring(var,1,10)=="temp_range"~"temp range"
    ))->res
  
  return(res)
}

process_univar_results_for_plot <- function(){

a<-b<-c<-NULL

if(file.exists("outputs/univar/chikv_univar_all.RDS")==T){
a<-readRDS("outputs/univar/chikv_univar_all.RDS")
a$path<-"chikv"
a<-a[order(a$waic),]
a$ORDER<-1:nrow(a)
}

if(file.exists("outputs/univar/denv_univar_all.RDS")==T){
b<-readRDS("outputs/univar/denv_univar_all.RDS")
b$path<-"denv"
b<-b[order(b$waic),]
b$ORDER<-1:nrow(b)
  }
  
if(file.exists("outputs/univar/zikv_univar_all.RDS")==T){
c<-readRDS("outputs/univar/zikv_univar_all.RDS")
c$path<-"zikv"
c<-c[order(c$waic),]
c$ORDER<-1:nrow(c)
}

data<-rbind(a,b,c)

data$var <-sub("min_min_temp_03","Previous 0-3 weeks", data$var)
data$var <-sub("min_min_temp_02","Previous 0-2 weeks", data$var)
data$var <-sub("min_min_temp_01","Previous 0-1 weeks", data$var)
data$var <-sub("min_min_temp_3","Week t -3", data$var)
data$var <-sub("min_min_temp_2","Week t -2", data$var)
data$var <-sub("min_min_temp_1","Week t -1", data$var)
data$var <-sub("min_min_temp","Week t", data$var)
data$var <-sub("max_max_temp_03","Previous 0-3 weeks", data$var)
data$var <-sub("max_max_temp_02","Previous 0-2 weeks", data$var)
data$var <-sub("max_max_temp_01","Previous 0-1 weeks", data$var)
data$var <-sub("max_max_temp_3","Week t -3", data$var)
data$var <-sub("max_max_temp_2","Week t -2", data$var)
data$var <-sub("max_max_temp_1","Week t -1", data$var)
data$var <-sub("max_max_temp","Week t", data$var)
data$var <-sub("min_min_hum_03","Previous 0-3 weeks", data$var)
data$var <-sub("min_min_hum_02","Previous 0-2 weeks", data$var)
data$var <-sub("min_min_hum_01","Previous 0-1 weeks", data$var)
data$var <-sub("min_min_hum_3","Week t -3", data$var)
data$var <-sub("min_min_hum_2","Week t -2", data$var)
data$var <-sub("min_min_hum_1","Week t -1", data$var)
data$var <-sub("min_min_hum","Week t", data$var)
data$var <-sub("max_max_hum_03","Previous 0-3 weeks", data$var)
data$var <-sub("max_max_hum_02","Previous 0-2 weeks", data$var)
data$var <-sub("max_max_hum_01","Previous 0-1 weeks", data$var)
data$var <-sub("max_max_hum_3","Week t -3", data$var)
data$var <-sub("max_max_hum_2","Week t -2", data$var)
data$var <-sub("max_max_hum_1","Week t -1", data$var)
data$var <-sub("max_max_hum","Week t", data$var)
data$var <-sub("max_hum_03","Previous 0-3 weeks", data$var)
data$var <-sub("max_hum_02","Previous 0-2 weeks", data$var)
data$var <-sub("max_hum_01","Previous 0-1 weeks", data$var)
data$var <-sub("max_hum_3","Week t -3", data$var)
data$var <-sub("max_hum_2","Week t -2", data$var)
data$var <-sub("max_hum_1","Week t -1", data$var)
data$var <-sub("max_hum","Week t", data$var)
data$var <-sub("min_hum_03","Previous 0-3 weeks", data$var)
data$var <-sub("min_hum_02","Previous 0-2 weeks", data$var)
data$var <-sub("min_hum_01","Previous 0-1 weeks", data$var)
data$var <-sub("min_hum_3","Week t -3", data$var)
data$var <-sub("min_hum_2","Week t -2", data$var)
data$var <-sub("min_hum_1","Week t -1", data$var)
data$var <-sub("min_hum","Week t", data$var)
data$var <-sub("mean_hum_03","Previous 0-3 weeks", data$var)
data$var <-sub("mean_hum_02","Previous 0-2 weeks", data$var)
data$var <-sub("mean_hum_01","Previous 0-1 weeks", data$var)
data$var <-sub("mean_hum_3","Week t -3", data$var)
data$var <-sub("mean_hum_2","Week t -2", data$var)
data$var <-sub("mean_hum_1","Week t -1", data$var)
data$var <-sub("mean_hum","Week t", data$var)
data$var <-sub("max_temp_03","Previous 0-3 weeks", data$var)
data$var <-sub("max_temp_02","Previous 0-2 weeks", data$var)
data$var <-sub("max_temp_01","Previous 0-1 weeks", data$var)
data$var <-sub("max_temp_3","Week t -3", data$var)
data$var <-sub("max_temp_2","Week t -2", data$var)
data$var <-sub("max_temp_1","Week t -1", data$var)
data$var <-sub("max_temp","Week t", data$var)
data$var <-sub("min_temp_03","Previous 0-3 weeks", data$var)
data$var <-sub("min_temp_02","Previous 0-2 weeks", data$var)
data$var <-sub("min_temp_01","Previous 0-1 weeks", data$var)
data$var <-sub("min_temp_3","Week t -3", data$var)
data$var <-sub("min_temp_2","Week t -2", data$var)
data$var <-sub("min_temp_1","Week t -1", data$var)
data$var <-sub("min_temp","Week t", data$var)
data$var <-sub("mean_temp_03","Previous 0-3 weeks", data$var)
data$var <-sub("mean_temp_02","Previous 0-2 weeks", data$var)
data$var <-sub("mean_temp_01","Previous 0-1 weeks", data$var)
data$var <-sub("mean_temp_3","Week t -3", data$var)
data$var <-sub("mean_temp_2","Week t -2", data$var)
data$var <-sub("mean_temp_1","Week t -1", data$var)
data$var <-sub("mean_temp","Week t", data$var)
data$var <-sub("temp_range_day_03","Previous 0-3 weeks", data$var)
data$var <-sub("temp_range_day_02","Previous 0-2 weeks", data$var)
data$var <-sub("temp_range_day_01","Previous 0-1 weeks", data$var)
data$var <-sub("temp_range_day_3","Week t -3", data$var)
data$var <-sub("temp_range_day_2","Week t -2", data$var)
data$var <-sub("temp_range_day_1","Week t -1", data$var)
data$var <-sub("temp_range_day","Week t", data$var)
data$var <-sub("cum_rain_03","Previous 0-3 weeks", data$var)
data$var <-sub("cum_rain_02","Previous 0-2 weeks", data$var)
data$var <-sub("cum_rain_01","Previous 0-1 weeks", data$var)
data$var <-sub("cum_rain_3","Week t -3", data$var)
data$var <-sub("cum_rain_2","Week t -2", data$var)
data$var <-sub("cum_rain_1","Week t -1", data$var)
data$var <-sub("cum_rain","Week t", data$var)
data$var <-sub("mean_rain_03","Previous 0-3 weeks", data$var)
data$var <-sub("mean_rain_02","Previous 0-2 weeks", data$var)
data$var <-sub("mean_rain_01","Previous 0-1 weeks", data$var)
data$var <-sub("mean_rain_3","Week t -3", data$var)
data$var <-sub("mean_rain_2","Week t -2", data$var)
data$var <-sub("mean_rain_1","Week t -1", data$var)
data$var <-sub("mean_rain","Week t", data$var)
data$var <-sub("v_heavy_raindays_03","Previous 0-3 weeks", data$var)
data$var <-sub("v_heavy_raindays_02","Previous 0-2 weeks", data$var)
data$var <-sub("v_heavy_raindays_01","Previous 0-1 weeks", data$var)
data$var <-sub("v_heavy_raindays_3","Week t -3", data$var)
data$var <-sub("v_heavy_raindays_2","Week t -2", data$var)
data$var <-sub("v_heavy_raindays_1","Week t -1", data$var)
data$var <-sub("v_heavy_raindays","Week t", data$var)
data$var <-sub("heavy_raindays_03","Previous 0-3 weeks", data$var)
data$var <-sub("heavy_raindays_02","Previous 0-2 weeks", data$var)
data$var <-sub("heavy_raindays_01","Previous 0-1 weeks", data$var)
data$var <-sub("heavy_raindays_3","Week t -3", data$var)
data$var <-sub("heavy_raindays_2","Week t -2", data$var)
data$var <-sub("heavy_raindays_1","Week t -1", data$var)
data$var <-sub("heavy_raindays","Week t", data$var)
data$var <-sub("rain_days_03","Previous 0-3 weeks", data$var)
data$var <-sub("rain_days_02","Previous 0-2 weeks", data$var)
data$var <-sub("rain_days_01","Previous 0-1 weeks", data$var)
data$var <-sub("raindays_3","Week t -3", data$var)
data$var <-sub("raindays_2","Week t -2", data$var)
data$var <-sub("raindays_1","Week t -1", data$var)
data$var <-sub("raindays","Week t", data$var)
data$var <-sub("flood_03","Previous 0-3 weeks", data$var)
data$var <-sub("flood_02","Previous 0-2 weeks", data$var)
data$var <-sub("flood_01","Previous 0-1 weeks", data$var)
data$var <-sub("flood_3","Week t -3", data$var)
data$var <-sub("flood_2","Week t -2", data$var)
data$var <-sub("flood_1","Week t -1", data$var)
data$var <-sub("flood_0","Week t", data$var)
data$var <-sub("cons_dryday_3","Week t -3", data$var)
data$var <-sub("cons_dryday_2","Week t -2", data$var)
data$var <-sub("cons_dryday_1","Week t -1", data$var)
data$var <-sub("cons_dryday","Week t", data$var)
data$var <-sub("cons_rainday_3","Week t -3", data$var)
data$var <-sub("cons_rainday_2","Week t -2", data$var)
data$var <-sub("cons_rainday_1","Week t -1", data$var)
data$var <-sub("cons_rainday","Week t", data$var)
data$var <-sub("mir_pop","Mid infra-red index", data$var)
data$var <-sub("income","Median household\nincome", data$var)
data$var <-sub("illiteracy","Rate of adult\nilliteracy", data$var)
data$var <-sub("sanitation","Proportion of\nhouseholds with\naccess to adequate\nsanitation", data$var)
data$var <-sub("BETWEEN","Between centrality", data$var)
data$var <-sub("ALPHA","Alpha centrality", data$var)
data$var <-sub("DEGREE","Degree centrality", data$var)
data$var <-sub("spei12","SPEI 12 months", data$var)
data$var <-sub("spei1","SPEI 1 month", data$var)
data$var <-sub("spei3","SPEI 3 months", data$var)
data$var <-sub("spei6","SPEI 6 months", data$var)

data$path <-toupper(data$path)

return(data)
}



