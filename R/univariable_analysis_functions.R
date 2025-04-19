#### univariable analysis

run_univar_analysis <- function(data, variables, save_progress=TRUE, pathogen){
  
  options(mc.cores = parallel::detectCores())
  
  list<-list()
  models<-list()
  var_names <- variables
  
  graph<-inla.read.graph(GRAPH)
  
  for(i in 1:length(var_names)){
    
    if(var_names[i] %in% c("flood_0","flood_1","flood_2" ,"flood_3","flood_01","flood_02",           
                           "flood_03")){ # categorical variable
      formula <- CASES ~ 1 + lag_cases + 
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=T) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        f(data[,var_names[i]],model="iid")
    
    }else{
      formula <- CASES ~ 1 + lag_cases + 
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
    beta<-models[[i]]$summary.fixed[3,]
    list[[i]]<-data.frame(var = variables[i], 
                          waic = waic,
                          beta_mean = beta$mean, 
                          beta_low = beta$`0.025quant`, 
                          beta_upp = beta$`0.975quant`)
    models[[i]]<-NA
    
    print(paste0("finished model ", i, " out of ", length(variables)))
  }
  result<- bind_rows(list)
  
  if(save_progress==TRUE){
   saveRDS(result, paste0("outputs/univar/temp/",pathogen,"_", i,".RDS")) 
  }
  return(result)
  
  }
  

## bind univariable model results together
binding_univar_results <-function(pathogen, file_path="outputs/univar/temp/"){
  filenames <- list.files(file_path, pattern=pathogen, full.names=TRUE)
  ldf <- lapply(filenames, readRDS)
  my.fcn<-function(x){
    ifelse("cat" %in% colnames(x),x$cat<-as.character(x$cat), x$cat<-NA)
    return(x)
  }
  test<-lapply(ldf,my.fcn)
  res<- bind_rows(test, id=NULL)
  res<- res[duplicated(res$var)==F,]

  res <- process_univar_names(res)
  
  return(res)
}


## choose the best univariable model
best_univar <- function(univar_results_table){
  univar_results_table %>%filter(beta_low < 0 & beta_upp <0 | beta_low > 0 & beta_upp >0)->sig
  sig$var[sig$waic==min(sig$waic)] ->var
  waic<-min(sig$waic)
  return(c(var,waic))
}


## plot univariable beta results
plot_univar_betas <- function(choice = "temperature"){
  
  if(choice!="temperature" & choice !="non_temperature"){
    warning("choice must be temperature or non_temperature")
  }
  
  data <- process_univar_results_for_plot()
  
  if(choice=="temperature"){
  plots<-list()
  sub_tot<-data[data$group=="Absolute minimum temp" |data$group=="Absolute maximum temp"|data$group=="Absolute minimum temp"|
                  data$group=="min temp"|data$group=="mean temp"|data$group=="max temp" |
                  data$group=="temp range",]
  sub_tot$group <- sub("temp", "temperature", sub_tot$group)
  sub_tot$group <- sub("min temperature", "Minimum temperature", sub_tot$group)
  sub_tot$group <- sub("mean temperature", "Mean temperature", sub_tot$group)
  sub_tot$group <- sub("max temperature", "Maximum temperature", sub_tot$group)
  sub_tot$group <- sub("temperature range", "Temperature range", sub_tot$group)
  
  for(i in 1:length(unique(sub_tot$group))){
    sub<-sub_tot[sub_tot$group == unique(sub_tot$group)[i],]
    lim1<- ifelse(min(sub$beta_low)< -max(sub$beta_upp),
                  min(sub$beta_low), - max(sub$beta_upp))
    lim2<-ifelse(min(sub$beta_low)< -max(sub$beta_upp),
                 - min(sub$beta_low), max(sub$beta_upp))
    plots[[i]]<-ggplot(sub)+geom_hline(aes(yintercept = 0), linetype="dashed") +
      geom_point(aes(x=var, y=beta_mean, col=path), 
                 size=2,position=position_dodge(width=0.2)) +
      geom_pointrange(aes(x=var, y=beta_mean, 
                          ymin=beta_low, 
                          ymax=beta_upp, col=path), 
                      size=0.3,position=position_dodge(width=0.2)) + 
      labs(x = NULL, y = "Beta value", title = unique(sub_tot$group)[i]) +
      coord_flip(ylim=c(lim1,lim2)) + #makes horizontal
      theme_classic() +
      theme(legend.title=element_blank(),legend.position="none",
            plot.title = element_text(hjust=0.5),
            plot.margin = margin(10, 12, 10, 10))#, 
    
  }
  plot1<-sjPlot::plot_grid(plots)
  
  return(plot1)
  }
  
  if(choice=="non_temperature"){
    
    data$group2<-data$group
    new<-data[data$group2 !="min temp"&
        data$group2 !="mean temp"&
        data$group2 != "max temp"&
        data$group2 !="Absolute maximum temp"& 
        data$group2 !="Absolute minimum temp"&
        data$group2 !="temp range",]
    
    new$group2 <- sub("hum", "humidity", new$group2)
    new$group2 <- sub("min humidity", "Minimum humidity", new$group2)
    new$group2 <- sub("mean humidity", "Mean humidity", new$group2)
    new$group2 <- sub("max humidity", "Maximum humidity", new$group2)
    new$group2[new$group2=="mean rain"]<-"Mean rain"
    new$group2[new$group2=="cumulative rain"]<-"Cumulative rain"
    new$group2[new$group2=="income"|new$group2=="illiteracy"|new$group2=="sanitation"]<-"Socioeconomic"
    new$group2[new$group2=="consecutive rain"]<-"Consecutive days of rain"
    new$group2[new$group2=="consecutive dry"]<-"Consecutive days of no rain"
    new$group2[new$group2=="flooding"]<-"Flooding"
    new$group2[new$group2=="MIR"]<-"Land use"
    new$group2[new$group2=="raindays"]<-"Days of rain"
    new$group2[new$group2=="very heavy raindays"]<-"Days of very heavy rain"
    new$group2[new$group2=="heavy raindays"]<-"Days of heavy rain"
    new$group2[new$var=="Between centrality"]<-"Human connectivity:\nbetweenness centrality"
    new$group2[new$var=="Alpha centrality"]<-"Human connectivity:\nalpha centrality"
    new$group2[new$var=="Degree centrality"]<-"Human connectivity:\ndegree centrality"
    
    
    plots2<-list()
    for(i in 1:length(unique(new$group2))){
      sub<-new[new$group2 == unique(new$group2)[i],]
      lim1<- ifelse(min(sub$beta_low,na.rm=T)< -max(sub$beta_upp,na.rm=T),
                    min(sub$beta_low,na.rm=T), - max(sub$beta_upp,na.rm=T))
      lim2<-ifelse(min(sub$beta_low,na.rm=T)< -max(sub$beta_upp,na.rm=T),
                   - min(sub$beta_low,na.rm=T), max(sub$beta_upp,na.rm=T))
      
      plots2[[i]]<-
        ggplot(sub)+ 
        geom_hline(aes(yintercept = 0), linetype="dashed") +
        geom_point(aes(x=var, y=beta_mean, col=path), 
                   size=2,position=position_dodge(width=0.2)) +
        geom_pointrange(aes(x=var, y=beta_mean, 
                            ymin=beta_low, 
                            ymax=beta_upp, col=path), 
                        size=0.3,position=position_dodge(width=0.2)) + 
        labs(x = NULL,
             y = "Beta value", title = unique(new$group2)[i]) +
        coord_flip(ylim=c(lim1,lim2)) + #makes horizontal
        theme_classic() +
        theme(legend.title=element_blank(),#legend.position="none",
              plot.title = element_text(hjust=0.5))
    }
    
    
    plot2<- sjPlot::plot_grid(plots2)
    
    return(plot2)
  }
  }
