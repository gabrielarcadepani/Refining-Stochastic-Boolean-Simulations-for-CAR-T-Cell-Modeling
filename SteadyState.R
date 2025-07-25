############
# Packages #
############
PKG=c('BoolNet','tidyverse','ggpubr',
      'doSNOW','foreach')
for (pkg in PKG) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

##############
# Parameters #
##############
net=loadNetwork("Networks/net_shah_3.1.txt")
f=tabela%>%
  summarise(m=round(mean(n)),
                     s=sd(n),
                     .by = c(l))%>%
  select(m)
R=f#tabela[,2]
p_d=0.05
p_a=0.5
Time_buffer=200
Time_window=80
wa=0.1
wb=0.5
eta=0.01
points=Time_buffer
Scenarios=data.frame(TAex=c(1,1,1,1,1),
                     IL2ex=c(0,1,1,1,1),
                     PDL1L2ex=c(0,0,1,1,0),
                     CD8086ex=c(0,0,1,0,1))

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

getFinalStates=function(Data,algorithm){
  if(algorithm=='SDDS'){
    x=Data%>%
      filter(Time==max(Time))%>%
      select(!Time)
  }else if(algorithm=='BKMC'){
    x=Data%>%
      filter(Time==-1)%>%
      select(!Time)
  }
  return(x)
}

getTrajectory=function(net,t0,Time_buffer,rate_up,rate_down,initialState,algorithm,points){
  if(algorithm=='SDDS'){
    #Initialization
    State=initialState
    n=length(net$genes)
    Trajectory=c()
    count=rep(0,n)
    names(count)=names(State)
    
    #Iterative process
    for (t in t0:(t0+Time_buffer)) {
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
        if(t>Time_buffer) break
      }else{
        break
      }
    }
    
    Final_State=as.vector(strsplit(Trajectory[(nrow(Trajectory)-1),2],'')[[1]],mode = 'integer')
    
    #Data consolidation
    time_tick=Time_buffer/points
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
    Trajectory=cbind('Time'=time_tick*c(1:points)+t0,
                     as.data.frame(bkmc[1:points,]))
    names(Trajectory)=c('Time',net$genes)
    Trajectory=rbind(c('Time'=(-1),Final_State),Trajectory)
  }
  return(Trajectory)
}

runSDDS=function(net,scenario,t0,Time_buffer,R,rate_up,rate_down,initialStates,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages = c('BoolNet'),
            .combine = 'rbind',
            .export = c('getInitialState','getTrajectory')
  ) %dopar% {
    if(is.null(initialStates)){
      initState=getInitialState(net,scenario,rate_up,rate_down)
    }else{
      initState=as.vector(initialStates[r,],mode='double')
    }
    getTrajectory(net,t0,Time_buffer,rate_up,rate_down,initState,
                  algorithm='SDDS')
  }
  stopCluster(cl)
  return(X)
}

runBKMC=function(net,scenario,t0,Time_buffer,R,rate_up,rate_down,initialStates,points,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages = c('BoolNet','dplyr'),
            .combine = 'rbind',
            .export = c('getInitialState','getTrajectory')
  ) %dopar% {
    if(is.null(initialStates)){
      initState=getInitialState(net,scenario,rate_up,rate_down)
    }else{
      initState=as.vector(initialStates[r,],mode='double')
    }    
    getTrajectory(net,t0,Time_buffer,rate_up,rate_down,initState,
                  algorithm='BKMC',points=points)
  }
  stopCluster(cl)
  X=as.data.frame(X)
  return(X)
}

runExperiment=function(net,t0,Time_buffer,R,p_a,p_d,scenario,algorithm,
                       initialStates=NULL,points=100,threads=10){
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
    return(runSDDS(net,scenario,t0,Time_buffer,R,rate_up,rate_down,initialStates,threads))
  }else if(algorithm=='BKMC'){
    return(runBKMC(net,scenario,t0,Time_buffer,R,rate_up,rate_down,initialStates,points,threads))
  }
}

###################
# Experiment SDDS #
###################
Data=FinalStates=Tabela=c()
for (l in 1:nrow(Scenarios)) {
  for(r in 1:10){
    t=Sys.time()
    print(paste0('Current scenario: ',l))
    t0=0
    data=c()
    condition=FALSE
    repeat{
      if(!t0){
        buffer=runExperiment(net,t0,Time_buffer,R[l,],p_a,p_d,Scenarios[l,],'SDDS')
      }else{
        buffer=runExperiment(net,t0,Time_buffer,R[l,],p_a,p_d,Scenarios[l,],'SDDS',initialStates = states)
      }
      data=rbind(data,buffer%>%
                   pivot_longer(!Time)%>%
                   summarise(N=n(),
                             p=sum(value)/N,
                             .by = c(Time,name))%>%
                   select(!N)%>%
                   pivot_wider(names_from = name,
                               values_from = p,
                               names_sort = FALSE)%>%
                   as.data.frame())
      p=data%>%pivot_longer(!Time)%>%
        filter(name=='CFUNC'|name=='CINHI'|name=='TAPOP')%>%
        ggplot()+geom_line(aes(Time,value,col=name))+ylim(0,1)+
        theme_bw()+theme(legend.position = 'none')
      #print(p)
      
      #Sliding window stop condition
      {
        for(i in (Time_window):(t0+Time_buffer)){
          X=data%>%filter(Time<i,Time>(i-Time_window))
          x=X%>%
            select(Time,CFUNC,CINHI,TAPOP)%>%
            pivot_longer(!Time)%>%
            summarise(p=mean(value),.by = c(Time,name))%>%
            summarise(q=mean(p),.by = name)
          a=X%>%
            select(Time,CFUNC,CINHI,TAPOP)%>%
            filter(Time<=(min(Time)+wa*Time_window))%>%
            pivot_longer(!Time)%>%
            summarise(p=mean(value),.by = c(Time,name))%>%
            summarise(q=mean(p),.by = name)
          b=X%>%
            select(Time,CFUNC,CINHI,TAPOP)%>%
            filter(Time>=(min(Time)+(1-wb)*Time_window))%>%
            pivot_longer(!Time)%>%
            summarise(p=mean(value),.by = c(Time,name))%>%
            summarise(q=mean(p),.by = name)
          
          Y=(b$q-a$q)/x$q
          
          for(j in 1:length(Y)){
            if(abs(Y[j])>eta)break
            if((j==3)&(abs(Y[j])<eta)) condition=!condition
          }
          
          print(i)
          print((b$q-a$q)/x$q)
          print(condition)
          if(condition) break
        }
      }
      if(condition) break
      states=getFinalStates(buffer,algorithm = 'SDDS')
      t0=t0+Time_buffer
    }
    #Data=rbind(Data,cbind('Scenario'=l,data%>%filter(Time<=i)))
    Tabela=rbind(Tabela,cbind(l,r,i))
    p=data%>%pivot_longer(!Time)%>%
      filter(name=='CFUNC'|name=='CINHI'|name=='TAPOP')%>%
      ggplot()+
      geom_line(aes(Time,value,col=name))+
      geom_vline(xintercept = i)+
      ylim(0,1)+
      theme_bw()+
      theme(legend.position = 'none')
    print(p)
    print(Sys.time()-t)
  }
}
#Results
Data%>%
  summarise(x=max(Time),.by=Scenario)
x11();Data%>%
  select(Scenario,Time,CFUNC,CINHI,TAPOP)%>%
  pivot_longer(!c(Scenario,Time))%>%
  ggplot()+geom_line(aes(Time,value,col=name))+
  ylim(0,1)+facet_grid(~Scenario)+
  theme_bw()+theme(legend.position = 'none')

###################
# Experiment BKMC #
###################
Data=c()
for (l in 1:nrow(Scenarios)) {
  for (r in 1:10) {
    t=Sys.time()
    print(paste0('Current scenario: ',l))
    t0=0
    data=c()
    condition=FALSE
    repeat{
      if(!t0){
        buffer=runExperiment(net,t0,Time_buffer,R[l,],p_a,p_d,Scenarios[l,],'BKMC',points = points)
      }else{
        buffer=runExperiment(net,t0,Time_buffer,R[l,],p_a,p_d,Scenarios[l,],'BKMC',initialStates = states,points = points)
      }
      data=rbind(data,buffer%>%
                   pivot_longer(!Time)%>%
                   summarise(N=n(),
                             p=sum(value)/N,
                             .by = c(Time,name))%>%
                   select(!N)%>%
                   pivot_wider(names_from = name,
                               values_from = p,
                               names_sort = FALSE)%>%
                   as.data.frame())
      
      p=data%>%pivot_longer(!Time)%>%
        filter(name=='CFUNC'|name=='CINHI'|name=='TAPOP')%>%
        ggplot()+geom_line(aes(Time,value,col=name))+ylim(0,1)+
        theme_bw()+theme(legend.position = 'none')
      print(p)
      
      #Sliding window stop condition
      {
        for(i in Time_window:(t0+Time_buffer)){
          X=data%>%filter(Time<i,Time>(i-Time_window))
          x=X%>%
            select(Time,CFUNC,CINHI,TAPOP)%>%
            pivot_longer(!Time)%>%
            summarise(p=mean(value),.by = c(Time,name))%>%
            summarise(q=mean(p),.by = name)
          a=X%>%
            select(Time,CFUNC,CINHI,TAPOP)%>%
            filter(Time<=(min(Time)+wa*Time_window))%>%
            pivot_longer(!Time)%>%
            summarise(p=mean(value),.by = c(Time,name))%>%
            summarise(q=mean(p),.by = name)
          b=X%>%
            select(Time,CFUNC,CINHI,TAPOP)%>%
            filter(Time>=(min(Time)+(1-wb)*Time_window))%>%
            pivot_longer(!Time)%>%
            summarise(p=mean(value),.by = c(Time,name))%>%
            summarise(q=mean(p),.by = name)
          
          Y=(b$q-a$q)/x$q
          
          for(j in 1:length(Y)){
            if(abs(Y[j])>eta)break
            if((j==3)&(abs(Y[j])<eta)) condition=!condition
          }
          
          print(i)
          print((b$q-a$q)/x$q)
          print(condition)
          if(condition) break
        }
      }
      if(condition) break
      states=getFinalStates(buffer,algorithm = 'BKMC') 
      t0=t0+Time_buffer
    }
    Data=rbind(Data,cbind('Scenario'=l,data%>%filter(Time<=i)))
    Tabela=rbind(Tabela,cbind(l,r,i))
    FinalStates=rbind(FinalStates,cbind('Scenario'=l,round(filter(buffer,Time>=(i-5),Time<=i),0)))
    print(Sys.time()-t)
  }
}
Data$Scenario=factor(Data$Scenario,
                     levels=unique(Data$Scenario),
                     labels=c('TA+, IL2-, PDL1/L2-, CD80/86-',
                              'TA+, IL2+, PDL1/L2-, CD80/86-',
                              'TA+, IL2+, PDL1/L2+, CD80/86+',
                              'TA+, IL2+, PDL1/L2+, CD80/86-',
                              'TA+, IL2+, PDL1/L2-, CD80/86+'))
#Results
setwd('Resultados/SteadyState')
if(!dir.exists('BKMC')){
  dir.create('BKMC')
}
setwd('BKMC')

save(Data,file='SteadyState.RData')
Tabela=Data%>%
  summarise(x=max(Time),.by=Scenario)
save(Tabela,file='Tabela.RData')

if(!dir.exists('Figuras')){
  dir.create('Figuras')
}
setwd('Figuras')

plot=Data%>%
  select(Scenario,Time,CFUNC,CINHI,TAPOP)%>%
  filter(Time>0)%>%
  pivot_longer(!c(Scenario,Time))%>%
  ggplot()+
  geom_line(aes(Time,
                value,
                col=name),
            linewidth=0.5)+
  ylim(0,1)+
  facet_wrap(~Scenario,
             scales='free_x',
             nrow=1,)+
  theme_bw()+
  theme(legend.position = 'none')+
  ylab('Activity level')+xlab('Time (a.u.)')
ggsave('SteadyState.pdf',plot,width =11,height = 4)
