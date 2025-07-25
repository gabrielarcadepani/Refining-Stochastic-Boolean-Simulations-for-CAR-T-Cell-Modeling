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
    State=paste0(initialState,collapse='')
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
net=loadNetwork("Networks/net_shah_1.
                0.txt")
R=10000
p_d=0.05
p_a=0.5 
Time=1000
Scenarios=data.frame(TAex=c(1,1,1,1,1),
                     IL2ex=c(0,1,1,1,1),
                     PDL1L2ex=c(0,0,1,1,0),
                     CD8086ex=c(0,0,1,0,1))
algorithm='SDDS'

##############
# Experiment #
##############
Data=c()
for (l in 1:nrow(Scenarios)) {
  t0=Sys.time()
  print(paste0('Current scenario: ',l))
  data=runExperiment(net,Time,R,p_a,p_d,Scenarios[l,],algorithm)
  Data=rbind(Data,cbind('Scenario'=l,data))
  print(Sys.time()-t0)
}

Data$Scenario=factor(Data$Scenario,
                     levels=unique(Data$Scenario),
                     labels=c('TA+, IL2-, PDL1/L2-, CD80/86-',
                              'TA+, IL2+, PDL1/L2-, CD80/86-',
                              'TA+, IL2+, PDL1/L2+, CD80/86+',
                              'TA+, IL2+, PDL1/L2+, CD80/86-',
                              'TA+, IL2+, PDL1/L2-, CD80/86+'))
setwd('Resultados')
if(!dir.exists('Kinetic_BN')){
  dir.create('Kinetic_BN')
}
setwd('Kinetic_BN')
save(Data,file=paste0('KBN_',algorithm,'.RData'))

plot=Data%>%
  select(Scenario,Time,CFUNC,CINHI,TAPOP)%>%
  pivot_longer(!c(Scenario,Time))%>%
  summarise(`Activity level`=mean(value),
            .by = c(Scenario,Time,name))%>%
  ggplot()+geom_line(aes(Time,`Activity level`))+
  ylim(0,1)+facet_grid(name~Scenario)+
  theme_bw()+theme(legend.title = element_blank())
ggsave(paste0('KBN_',algorithm,'.pdf'),plot,
       width = 11,height = 6)

plot=Data%>%
  filter(Scenario=='TA+, IL2+, PDL1/L2+, CD80/86+')%>%
  select(Time,CFUNC)%>%
  pivot_longer(!c(Time))%>%
  summarise(`Activity level`=mean(value),
            .by = c(Time,name))%>%
  ggplot()+geom_line(aes(Time,`Activity level`))+
  ylim(0,1)+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'none')
ggsave(paste0('Stationarity',algorithm,'.jpg'),plot,
       width = 9,height = 3)
