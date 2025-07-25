############
# Packages #
############
PKG=c('BoolNet','doSNOW','foreach',
      'ggpubr','ggcorrplot','tidyverse',
      'latex2exp')
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

getInitialState=function(net,scenario,rate_up,rate_down,type){
  if(type=='random'){
    State=generateState(net,specs = c('LCK'=1))
    for (i in 1:length(net$genes)) {
      State[i]=rbinom(1,1,0.5)
    }
    State['LCK']=1
    State['TAex']=ifelse(scenario['TAex'],1,0)
    State['IL2ex']=ifelse(scenario['IL2ex'],1,0)
    State['PDL1L2ex']=ifelse(scenario['PDL1L2ex'],1,0)
    State['CD8086ex']=ifelse(scenario['CD8086ex'],1,0)
  }else if(type=='normal'){
    State=generateState(net,specs = c('LCK'=1))
    State['TAex']=ifelse(scenario['TAex'],1,0)
    State['IL2ex']=ifelse(scenario['IL2ex'],1,0)
    State['PDL1L2ex']=ifelse(scenario['PDL1L2ex'],1,0)
    State['CD8086ex']=ifelse(scenario['CD8086ex'],1,0)
  }
  return(State)
}

getTrajectory=function(net,Time,rate_up,rate_down,initialState,algorithm){
  if(algorithm=='SDDS'){
    State=initialState
    n=length(net$genes)
    Trajectory=c()
    count=rep(0,n)
    names(count)=names(State)
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
    Trajectory=as.data.frame(Trajectory)
    names(Trajectory)=c('Time',net$genes)
  }else if(algorithm=='BKMC'){
    State=paste0(initState,collapse='')
    n=length(net$genes)
    
    #Iterative process
    t=0
    Trajectory=c()
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
  }
  return(Trajectory)
}

runSDDS=function(net,scenario,Time,R,rate_up,rate_down,mode,threads,type){
  if(mode=='+'){
    cl=makeCluster(threads,type='SOCK')
    registerDoSNOW(cl)
    X=foreach(r=1:R,
              .packages=c('BoolNet'),
              .combine = '+',
              .export = c('getInitialState','getTrajectory')
    ) %dopar% {
      initState=getInitialState(net,scenario,rate_up,rate_down,type)
      as.matrix(getTrajectory(net,Time,rate_up,rate_down,initState,algorithm='SDDS'))
    }
    stopCluster(cl)
    X=as.data.frame(X/R)
  } else if(mode=='rbind'){
    cl=makeCluster(threads,type='SOCK')
    registerDoSNOW(cl)
    X=foreach(r=1:R,
              .packages=c('BoolNet'),
              .combine = 'rbind',
              .export = c('getInitialState','getTrajectory')
    ) %dopar% {
      initState=getInitialState(net,scenario,rate_up,rate_down,type)
      getTrajectory(net,Time,rate_up,rate_down,initState,algorithm='SDDS')
    }
    stopCluster(cl)
  }
  return(X)
}

runBKMC=function(net,scenario,Time,R,rate_up,rate_down,num_points,threads,type){
  time_tick=Time/num_points
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,.packages=c('BoolNet','dplyr'),.combine = '+',
            .export = c('getInitialState','getTrajectory')) %dopar% {
              initState=getInitialState(net,scenario,rate_up,rate_down,type)
              data=getTrajectory(net,Time,rate_up,rate_down,initState,algorithm='BKMC')
              data=as.data.frame(data)
              data[,1]=as.numeric(data[,1])
              names(data)=c('Time','State')
              n=nrow(data)
              bkmc=matrix(0,ncol=2*length(net$genes),nrow = (last(data$Time)%/%time_tick+1))
              for(i in 1:(n-1)){
                p=data$Time[i]%/%time_tick+1
                q=data$Time[i+1]%/%time_tick+1
                if(p==q){ #same time_tick
                  delta_t=(data$Time[i+1]-data$Time[i])/time_tick
                  bkmc[p,]=bkmc[p,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]),
                                      as.numeric(strsplit(data$State[i],'')[[1]]))
                }else if(q-p==1){ #consecutive time_ticks
                  delta_t=(time_tick*p-data$Time[i])/time_tick
                  bkmc[p,]=bkmc[p,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]),
                                      as.numeric(strsplit(data$State[i],'')[[1]]))
                  delta_t=(data$Time[i+1]-time_tick*p)/time_tick
                  bkmc[q,]=bkmc[q,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]),
                                      as.numeric(strsplit(data$State[i],'')[[1]]))
                }else{ #spread time_ticks
                  delta_t=(time_tick*p-data$Time[i])/time_tick
                  bkmc[p,]=bkmc[p,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]),
                                      as.numeric(strsplit(data$State[i],'')[[1]]))
                  for(j in (p+1):(q-1)){
                    bkmc[j,]=bkmc[j,]+c(as.numeric(strsplit(data$State[i],'')[[1]]),
                                        as.numeric(strsplit(data$State[i],'')[[1]]))
                  }
                  delta_t=(time_tick*q-data$Time[i+1])/time_tick
                  bkmc[q,]=bkmc[q,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]),
                                      as.numeric(strsplit(data$State[i],'')[[1]]))
                }
              }
              bkmc=bkmc[1:num_points,]
              bkmc
            }
  stopCluster(cl)
  X=X[,1:length(net$genes)]
  X=as.data.frame(X/R)
  names(X)=net$genes
  X=cbind('Time'=time_tick*c(1:num_points),X)
  return(X)
}

runExperiment=function(net,Time,R,p_a,p_d,scenario,algorithm,
                       mode='+',num_points=100,
                       threads=10,type='normal'){
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
    return(runSDDS(net,scenario,Time,R,rate_up,rate_down,mode,threads,type))
  }else if(algorithm=='BKMC'){
    return(runBKMC(net,scenario,Time,R,rate_up,rate_down,num_points,threads,type))
  }
}

##############
# Parameters # 
##############
p_d=0.05
p_a=0.5 
Time=1000
R=10000
Scenarios=data.frame(TAex=c(1,1,1,1,1),
                     IL2ex=c(0,1,1,1,1),
                     PDL1L2ex=c(0,0,1,1,0),
                     CD8086ex=c(0,0,1,0,1))

##############
# Experiment #
##############
Nets=c('net_shah_1.0','net_shah_3.1')
Data=c()
for(Net in Nets){
  print(Net)
  #Network uploading
  setwd('Networks')
  net=loadNetwork(paste0(Net,'.txt'))
  for (l in 1:nrow(Scenarios)) {
    t0=Sys.time()
    print(paste0('Current scenario: ',l))
    scenario=Scenarios[l,]
    data=runExperiment(net,Time,R,p_a,p_d,scenario,algorithm = 'BKMC',num_points = Time)
    data=cbind('Version'=Net,'Scenario'=l,data)
    Data=rbind(Data,data)
    print(Sys.time()-t0)
  } 
}
Data$Version=factor(Data$Version,
                    labels=c('Original',
                             'Alternative'))
Data$Scenario=factor(Data$Scenario,
                     levels=unique(Data$Scenario),
                     labels=c('TA+, IL2-, PDL1/L2-, CD80/86-',
                              'TA+, IL2+, PDL1/L2-, CD80/86-',
                              'TA+, IL2+, PDL1/L2+, CD80/86+',
                              'TA+, IL2+, PDL1/L2+, CD80/86-',
                              'TA+, IL2+, PDL1/L2-, CD80/86+'))


#########
# Plots #
#########
setwd('../Resultados/Versions')
save(Data,file='Comparação.RData')
plot=Data%>%
  select(Version,Scenario,Time,CFUNC,CINHI,TAPOP)%>%
  pivot_longer(!c(Version,Scenario,Time))%>%
  summarise(`% Activation` = mean(value),
            .by = c(Version,Scenario,Time,name))%>%
  ggplot()+
  geom_line(aes(Time,`% Activation`,
                linetype = Version))+
  ylim(0,1)+facet_grid(name~Scenario)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.title = element_blank())+
  labs(x='Time (a.u.)',y='Activity level')
ggsave('Comparison.pdf',plot,width =11,height =6)

if (!dir.exists('Simulações')) {
  dir.create('Simulações')
}
setwd('Simulações')
cl=makeCluster(10,type='SOCK')
registerDoSNOW(cl)
foreach(i=1:length(net$genes),.packages = 'tidyverse')%dopar%{
  X=Data%>%
    select(Version,Scenario,Time,i+3)%>%
    pivot_longer(!c(Version,Scenario,Time))%>%
    summarise(`% Activation` = mean(value),
              .by = c(Version,Scenario,Time,name))%>%
    ggplot()+
    geom_line(aes(Time,`% Activation`,
                  col = Version))+
    ylim(0,1.1)+facet_wrap(~Scenario,nrow=1)+
    theme_bw()+ggtitle(paste0(net$genes[i]))
  ggsave(paste(net$genes[i],'.jpg'),X,
         width = 10,height = 3)
}
stopCluster(cl)