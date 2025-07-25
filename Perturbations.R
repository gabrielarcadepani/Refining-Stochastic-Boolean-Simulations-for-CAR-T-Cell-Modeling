############
# Packages #
############
PKG=c('BoolNet','tidyverse',
      'doSNOW','foreach')
for (pkg in PKG) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

#############
# Functions #
#############
getRates=function(net,specs,default=0){
  n=length(net$genes)
  rate=vector(mode='numeric',length = n)
  names(rate)=net$genes
  for (i in 1:n) {
    rate[i]=default
  }
  for(i in names(specs)){
    rate[which(names(rate)==i)]=specs[i]
  }
  return(rate)
}

getInitialState=function(net,scenario,rate_up,rate_down){
  State=generateState(net,specs = c('LCK'=1))
  State['TAex']=ifelse(scenario['TAex'],1,0)
  State['IL2ex']=ifelse(scenario['IL2ex'],1,0)
  State['PDL1L2ex']=ifelse(scenario['PDL1L2ex'],1,0)
  State['CD8086ex']=ifelse(scenario['CD8086ex'],1,0)
  return(State)
}

getTrajectory=function(net,Time,rate_up,rate_down,initialState,algorithm,num_points){
  if(algorithm=='SDDS'){
    #Initialization
    State=initialState
    n=length(net$genes)
    Trajectory=c()
    count=rep(0,n)
    names(count)=names(State)
    
    #Iterative process
    for (t in 1:Time) {
      s=sample(1:n,n)
      for (i in s) {
        count[i]=count[i]+1
        if (State[i]) { #1
          if (runif(1,0,1)<rate_down[i]&count[i]>1) {
            if (i==1|i==2|i==3|i==4) {
              State[i]=0
            } else {
              State[i]=stateTransition(net,State,'asynchronous',chosenGene = i)[i]
            }
            count[i]=0
          } else {
            #State[i]=State[i]
          }
        } else { #0
          if (runif(1,0,1)<rate_up[i]&count[i]>1) {
            if(i==1|i==2|i==3|i==4){
              State[i]=1
            } else {
              State[i]=stateTransition(net,State,'asynchronous',chosenGene = i)[i]
            }
            count[i]=0
          } else {
            #State[i]=State[i]
          }
        }
      }
      Trajectory=rbind(Trajectory,c(t,State))
    }
    
    #Consolidation
    Trajectory=as.data.frame(Trajectory)
    names(Trajectory)=c('Time',net$genes)
    
  }else if(algorithm=='BKMC'){
    #Initialization
    State=paste0(initState,collapse='')
    n=length(net$genes)
    t=0
    Trajectory=cbind(t,State)
    
    #Iterative process
    repeat{
      #get all the possible next states from the current state
      possible_states=c()
      for(i in 1:n){
        x=stateTransition(net,as.numeric(strsplit(State,'')[[1]]),
                          type = 'asynchronous',chosenGene = i)
        if(x[i]!=as.numeric(strsplit(State,'')[[1]][i])){
          possible_states=rbind(possible_states,
                                c(paste(x,collapse = ''),
                                  ifelse(x[i],rate_up[i],rate_down[i])))
        }
      }
      
      #Select the next state based on the ordered rates
      if(!is.null(nrow(possible_states))){
        possible_states=as.data.frame(possible_states)
        names(possible_states)=c('State','rate')
        possible_states=possible_states[order(possible_states$rate,decreasing = FALSE),]
        possible_states['upper_bound']=as.numeric(possible_states$rate)/sum(as.numeric(possible_states$rate))
        lower_bound=c(0)
        if(nrow(possible_states)>1){
          for (i in 2:nrow(possible_states)) {
            lower_bound=rbind(lower_bound,possible_states[i-1,3])
            possible_states[i,3]=possible_states[i,3]+possible_states[i-1,3]
          }
        }
        possible_states=cbind(possible_states,lower_bound)
        u=runif(1,0,1)
        for (j in 1:nrow(possible_states)) {
          if(u>possible_states[j,4]&u<possible_states[j,3]){
            State=possible_states$State[j]
            break
          }
        }
        tau=-log(runif(1,0,1))/sum(as.numeric(possible_states$rate))
        t=t+tau
        Trajectory=rbind(Trajectory,c(t,State))
        if(t>Time) break
      }else{
        break
      }
    }
    
    #Consolidation
    time_tick=Time/num_points
    Trajectory=as.data.frame(Trajectory)
    Trajectory[,1]=as.numeric(Trajectory[,1])
    names(Trajectory)=c('Time','State')
    n=nrow(Trajectory)
    bkmc=matrix(0,ncol=length(net$genes),
                nrow = (last(Trajectory$Time)%/%time_tick+1))
    for(i in 1:(n-1)){
      p=Trajectory$Time[i]%/%time_tick+1
      q=Trajectory$Time[i+1]%/%time_tick+1
      if(p==q){ #same time_tick
        delta_t=(Trajectory$Time[i+1]-Trajectory$Time[i])/time_tick
        bkmc[p,]=bkmc[p,]+delta_t*as.numeric(strsplit(Trajectory$State[i],'')[[1]])
      }else if(q-p==1){ #consecutive time_ticks
        delta_t=(time_tick*p-Trajectory$Time[i])/time_tick
        bkmc[p,]=bkmc[p,]+delta_t*as.numeric(strsplit(Trajectory$State[i],'')[[1]])
        delta_t=(Trajectory$Time[i+1]-time_tick*p)/time_tick
        bkmc[q,]=bkmc[q,]+delta_t*as.numeric(strsplit(Trajectory$State[i],'')[[1]])
      }else{ #spread time_ticks
        delta_t=(time_tick*p-Trajectory$Time[i])/time_tick
        bkmc[p,]=bkmc[p,]+delta_t*as.numeric(strsplit(Trajectory$State[i],'')[[1]])
        for(j in (p+1):(q-1)){
          bkmc[j,]=bkmc[j,]+as.numeric(strsplit(Trajectory$State[i],'')[[1]])
        }
        delta_t=(Trajectory$Time[i+1]-(time_tick*(q-1)))/time_tick
        bkmc[q,]=bkmc[q,]+delta_t*as.numeric(strsplit(Trajectory$State[i],'')[[1]])
      }
    }
    Trajectory=cbind('Time'=time_tick*c(1:num_points),
                     as.data.frame(bkmc[1:num_points,]))
    names(Trajectory)=c('Time',net$genes)
  }
  return(Trajectory)
}

runSDDS=function(net,scenario,Time,R,rate_up,rate_down,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages = c('BoolNet'),
            .combine = 'rbind',
            .export = c('getInitialState','getTrajectory')
  ) %dopar% {
    initState=getInitialState(net,scenario,rate_up,rate_down)
    getTrajectory(net,Time,rate_up,rate_down,initState,algorithm='SDDS')
  }
  stopCluster(cl)
  return(X)
}

runBKMC=function(net,scenario,Time,R,rate_up,rate_down,points,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages = c('BoolNet','dplyr'),
            .combine = 'rbind',
            .export = c('getInitialState','getTrajectory')
  ) %dopar% {
    initState=getInitialState(net,scenario,rate_up,rate_down)
    getTrajectory(net,Time,rate_up,rate_down,initState,
                  algorithm='BKMC',num_points=points)
  }
  stopCluster(cl)
  X=as.data.frame(X)
  return(X)
}

runExperiment=function(net,Time,R,p_a,p_d,scenario,algorithm,
                       points=100,threads=10){
  rate_up=getRates(net,
                   c('TAex'=ifelse(scenario['TAex'],p_a,0),
                     'IL2ex'=ifelse(scenario['IL2ex'],p_a,0),
                     'CD8086ex'=ifelse(scenario['CD8086ex'],p_a,0),
                     'PDL1L2ex'=ifelse(scenario['PDL1L2ex'],p_a,0)),
                   default = p_a)
  
  rate_down=getRates(net,
                     c('TAex'=ifelse(scenario['TAex'],p_d,0),
                       'IL2ex'=ifelse(scenario['IL2ex'],p_d,0),
                       'CD8086ex'=ifelse(scenario['CD8086ex'],p_d,0),
                       'PDL1L2ex'=ifelse(scenario['PDL1L2ex'],p_d,0)),
                     default = p_d)
  
  if(algorithm=='SDDS'){
    return(runSDDS(net,scenario,Time,R,rate_up,rate_down,threads))
  }else if(algorithm=='BKMC'){
    return(runBKMC(net,scenario,Time,R,rate_up,rate_down,points,threads))
  }
}

##############
# Parameters #
##############
R=10000
p_d=0.05
p_a=0.5 
t=1000
Scenarios=data.frame(TAex=c(1,1,1,1,1),
                     IL2ex=c(0,1,1,1,1),
                     PDL1L2ex=c(0,0,1,1,0),
                     CD8086ex=c(0,0,1,0,1))
Nets=c('net_shah_1.0','net_shah_3.1')
mutations=c('AKT','AP1','BATF','CD28',
            'CD3zeta','CFOS','CJUN','CREB',
            'CTLA4ex','Calcineurin','Calcium',
            'Calmodulin','DAG','ERK','GRB2SOS',
            'IL2ge','IP3','JAK13','LAT','LCK',
            'MTORC1','NFAT','NFKB','P38','PD1ex',
            'PDCD1ge','PDPK1','PI3K','PKCtheta',
            'PLCgamma1','RASGRP','SHP1','SHP2',
            'SLP76','STAT5','WASPCDC42','ZAP70',
            'dNFAT')

##############
# Experiment #
##############
CT=Data=c()

#Control
for(Net in Nets){
  t0=Sys.time()
  net=loadNetwork(paste0('Networks/',Net,'.txt'))
  print(paste0('Control: ',Net))
  ctrl=runExperiment(net,t,R,p_a,p_d,Scenarios[3,],'SDDS')
  CT=rbind(CT,cbind('Net'=Net,ctrl%>%filter(`Time`==t)%>%
                      select(Time,CFUNC,CINHI,TAPOP)%>%
                      pivot_longer(!Time)%>%
                      summarise(value=mean(value),
                                .by = c(Time,name))%>%
                      pivot_wider(values_from = value)%>%
                      select(!Time)))
  print(Sys.time()-t0)
}
CT$Net=factor(CT$Net,
              labels=c('Original',
                       'Alternative'))
#Perturbations
for (Net in Nets) {
  print(Net)
  net=loadNetwork(paste0('Networks/',Net,'.txt'))
  for (i in mutations) {
    t0=Sys.time()
    print(paste0('Current perturbation: ',i,' (',which(mutations==i),'/',length(mutations),')'))
    
    #KO
    net=fixGenes(net,i,c(0))
    data=runExperiment(net,t,R,p_a,p_d,Scenarios[3,],'SDDS')
    Data=rbind(Data,cbind('Net'=Net,'Perturbation'=i,type='KO',
                          data%>%filter(`Time`==t)%>%
                            select(Time,CFUNC,CINHI,TAPOP)%>%
                            pivot_longer(!Time)%>%
                            summarise(value=mean(value),
                                      .by = c(Time,name))%>%
                            pivot_wider(values_from = value)%>%
                            select(!Time)))
    
    #OE
    net=fixGenes(net,i,c(1))
    data=runExperiment(net,t,R,p_a,p_d,Scenarios[3,],'SDDS')
    Data=rbind(Data,cbind('Net'=Net,'Perturbation'=i,type='OE',
                          data%>%filter(`Time`==t)%>%
                            select(Time,CFUNC,CINHI,TAPOP)%>%
                            pivot_longer(!Time)%>%
                            summarise(value=mean(value),
                                      .by = c(Time,name))%>%
                            pivot_wider(values_from = value)%>%
                            select(!Time)))
    
    net=fixGenes(net,i,c(-1))
    print(Sys.time()-t0)
  }
}
Data$Net=factor(Data$Net,
              labels=c('Original',
                       'Alternative'))

###########################
# Processing and plotting #
###########################
setwd("Resultados")
if(!dir.exists('Perturbations')){
  dir.create('Perturbations')
  setwd('Perturbations')
  if(!dir.exists('Dados')){
    dir.create('Dados')
  }
  if(!dir.exists('Figuras')){
    dir.create('Figuras')
  }
}else{
  setwd('Perturbations')
  if(!dir.exists('Dados')){
    dir.create('Dados')
  }
  if(!dir.exists('Figuras')){
    dir.create('Figuras')
  }
}

setwd('Dados')
save(CT,file=paste0('Control_',R,'.RData'))
save(Data,file=paste0('Perturbations_',R,'.RData'))
setwd('..')
setwd('Figuras')

#Remapped
{
  {
    Figura=Data%>%
      group_by(Perturbation,type)%>%
      mutate(delta_CFUNC=ifelse(CFUNC>CT[1,2],
                                (100*(CFUNC-CT[1,2])/(1-CT[1,2])),
                                (100*(CFUNC/CT[1,2]-1))),
             delta_CINHI=ifelse(CINHI>CT[1,3],
                                (100*(CINHI-CT[1,3])/(1-CT[1,3])),
                                (100*(CINHI/CT[1,3]-1))),
             delta_TAPOP=ifelse(TAPOP>CT[1,4],
                                (100*(TAPOP-CT[1,4])/(1-CT[1,4])),
                                (100*(TAPOP/CT[1,4]-1))))
  }
  Fig=Figura%>%
    select(!c(CFUNC,CINHI,TAPOP))%>%
    pivot_longer(cols = c(delta_CFUNC,delta_CINHI,delta_TAPOP),
                 names_to = 'Output')%>%
    rowwise()%>%
    mutate(Output=strsplit(Output,'_')[[1]][2])%>%
    ggplot()+
    geom_point(aes(value,
                   factor(Perturbation,levels = rev(mutations)), 
                   colour = type, 
                   shape = Net), 
               alpha = 0.5,size=2)+
    xlab('% Change')+
    geom_vline(xintercept = 0,linetype = 2)+
    theme_bw()+
    facet_wrap(~Output)+
    ylab('')+
    theme(legend.position = 'right',
          legend.title = element_blank(),
          axis.text.x =element_text(size=7))
  ggsave('Perturbations_remapped.pdf',Fig,width = 6,height = 8)
}

#Normal
{
  {
    Figura=Data%>%
      group_by(Perturbation,type)%>%
      mutate(delta_CFUNC=(100*(CFUNC/CT[1,2]-1)),
             delta_CINHI=(100*(CINHI/CT[1,3]-1)),
             delta_TAPOP=(100*(TAPOP/CT[1,4]-1)))
  }

  range=Figura%>%
    select(!c(CFUNC,CINHI,TAPOP))%>%ungroup()%>%
    pivot_longer(cols = c(delta_CFUNC,delta_CINHI,delta_TAPOP),
                 names_to = 'Output')%>%
    summarise(max(abs(value)))
  
  Fig=Figura%>%
    select(Net,Perturbation,type,delta_TAPOP,delta_CINHI,delta_CFUNC)%>%
    pivot_longer(cols = c(delta_TAPOP,delta_CINHI,delta_CFUNC),
                 names_to = 'Output')%>%
    rowwise()%>%
    mutate(Output=strsplit(Output,'_')[[1]][2])%>%
    ggplot()+
    geom_point(aes(x=value,
                   y=factor(Perturbation,levels=rev(mutations)),
                   colour=type,
                   shape=Net),
               alpha=0.5,size=2)+
    xlab('% Change')+
    xlim(-range$`max(abs(value))`,range$`max(abs(value))`)+
    geom_vline(xintercept = 0,linetype=2)+
    theme_bw()+
    facet_wrap(~Output)+
    ylab('')+
    theme(legend.position = 'right',
          legend.title = element_blank(),
          axis.text.x =element_text(size=7))
  ggsave('Perturbations_normal.pdf',Fig,width = 6,height = 8)
}
