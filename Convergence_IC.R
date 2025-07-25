############
# Packages #
############
PKG=c('BoolNet','tidyverse',
      'doSNOW','foreach','ggpubr')
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
                       points=100,threads=12){
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
net=loadNetwork("Networks/net_shah_3.1.txt")
p_d=0.05
p_a=0.5 
Time=500
max_iter=11000
buff_size=200
num_points=Time
Scenarios=data.frame(TAex=c(1,1,1,1,1),
                     IL2ex=c(0,1,1,1,1),
                     PDL1L2ex=c(0,0,1,1,0),
                     CD8086ex=c(0,0,1,0,1))

###############
# Experiments #
###############
# Granularity
TOTAL=c()
for (r in c(20,200,2000)) {
  t0=Sys.time()
  print(r)
  data=runExperiment(net,Time,r,p_a,p_d,Scenarios[3,],'BKMC',points = Time)
  TOTAL=rbind(TOTAL,cbind(r,data))
  print(Sys.time()-t0)
}
setwd('Resultados/Convergence/BKMC')
save(file = 'Granularity.RData',TOTAL)

plot=TOTAL%>%
  select(r,Time,CFUNC,CINHI,TAPOP)%>%
  pivot_longer(!c(r,Time))%>%
  summarise(value=mean(value),.by = c(r,name,Time))%>%
  ggplot()+geom_line(aes(Time,value,col=name))+
  facet_wrap(~r,nrow = 1)+
  theme_bw()+
  labs(x='Time (a.u.)',y='Activity level')+
  ylim(0,1)+
  theme(legend.position = 'none')
ggsave(filename='Granularity.pdf',
       plot=plot,
       path ='Figuras',
       width = 6,height = 3)

# Convergence SDDS 
# Control
Original=c()
for (l in 1:nrow(Scenarios)) {
  t0=Sys.time()
  print(paste0('Current scenario: ',l))
  original=runExperiment(net,Time,10000,p_a,p_d,Scenarios[l,],'SDDS')
  Original=rbind(Original,cbind('Scenario'=l,original))
  print(Sys.time()-t0)
}
Original$Scenario=factor(Original$Scenario,
                         levels=unique(Original$Scenario),
                         labels=c('TA+, IL2-, PDL1/L2-, CD80/86-',
                                  'TA+, IL2+, PDL1/L2-, CD80/86-',
                                  'TA+, IL2+, PDL1/L2+, CD80/86+',
                                  'TA+, IL2+, PDL1/L2+, CD80/86-',
                                  'TA+, IL2+, PDL1/L2-, CD80/86+'))

TOTAL=conv_media=conv_delta=tabela=c()
for(r in 1:10){
  for (l in 1:nrow(Scenarios)) {
    for (epsilon in c(0.1,0.01,0.0001,0.0001)) {
      t0=Sys.time()
      scenario=Scenarios[l,]
      print(r)
      print(scenario)
      print(epsilon)
      total=as.data.frame(matrix(0,nrow=Time,ncol=length(net$genes)+1))
      names(total)=c('Time',net$genes)
      i=0 
      condition=FALSE
      repeat{
        buffer=runExperiment(net,Time,buff_size,p_a,p_d,scenario,'SDDS')
        for (j in 1:buff_size) {
          #print(paste0('Scenario ',l,' / Iteration ',buff_size*i+j))
          
          # Add new Trajectories 
          total=total+as.matrix(buffer[(Time*(j-1)+1):(Time*j),])
          
          # Stop condition
          y=total/(buff_size*i+j)
          curr=y%>%
            filter(Time==125|Time==250|Time==375)%>%
            select(Time,CFUNC,CINHI,TAPOP)
          
          if((buff_size*i+j)>1){
            delta=(curr-prev)/prev
            for(a in 1:nrow(delta)){
              for (b in 2:ncol(delta)) {
                if(abs(delta[a,b])>epsilon|is.na(delta[a,b]))break
              }
              if(abs(delta[a,b])>epsilon|is.na(delta[a,b]))break
            }
            if (a==nrow(delta)&b==ncol(delta)) condition=!condition

            
          }
          
          # conv_media=rbind(conv_media,
          #                  cbind('Scenario'=l,
          #                        'Time'=c(125,250,375),
          #                        'Trajectories'=(buff_size*i+j),
          #                        curr[,-1]))
          
          if(condition | (buff_size*i+j)==max_iter) {
            break
          }
          prev=curr
        }
        if(condition | (buff_size*i+j)==max_iter) {
          break
        }
        i=i+1
      }
      #TOTAL=rbind(TOTAL,cbind('Scenario'=l,'N'=(buff_size*i+j),'epsilon'=epsilon,total))
      tabela=rbind(tabela,cbind(r,l,epsilon,(buff_size*i+j),Sys.time()-t0))
      print(Sys.time()-t0)
    }
  }
}
save(tabela,file = 'Tabela_conv_00001.RData')
tabela%>%summarise(m=mean(V4),
                   s=sd(V4),
                   .by = c(l,epsilon))

TOTAL$Scenario=factor(TOTAL$Scenario,
                      levels=unique(TOTAL$Scenario),
                      labels=c('TA+, IL2-, PDL1/L2-, CD80/86-',
                               'TA+, IL2+, PDL1/L2-, CD80/86-',
                               'TA+, IL2+, PDL1/L2+, CD80/86+',
                               'TA+, IL2+, PDL1/L2+, CD80/86-',
                               'TA+, IL2+, PDL1/L2-, CD80/86+'))

#Results
{
  setwd('C:/Users/Gabri/Desktop/Mestrado/Resultados/Convergence')
  if(!dir.exists('SDDS')){
    dir.create('SDDS')
  }
  setwd('SDDS')
  save(Original,file='Original.RData')
  save(TOTAL,file='Convergência.RData')
  save(conv_delta,file='Conv_delta.RData')
  save(conv_media,file='Conv_media.RData')
  save(tabela,file='Tabela.RData')
  print(tabela)
  
  #delta
  conv_delta$Scenario=factor(conv_delta$Scenario,
                             levels=unique(conv_delta$Scenario),
                             labels=c('TA+, IL2-, PDL1/L2-, CD80/86-',
                                      'TA+, IL2+, PDL1/L2-, CD80/86-',
                                      'TA+, IL2+, PDL1/L2+, CD80/86+',
                                      'TA+, IL2+, PDL1/L2+, CD80/86-',
                                      'TA+, IL2+, PDL1/L2-, CD80/86+'))
  conv_delta%>%
    pivot_longer(!c(Scenario,Time,Trajectories))%>%
    filter(name=='TAPOP',Time==125,Scenario==1)%>%
    ggplot()+
    geom_point(aes(Trajectories,abs(value),col=factor(Time)))+
    geom_line(data=cbind('x'=1:(1/epsilon),'y'=1/c(1:(1/epsilon))),aes(x,y))+
    facet_grid(Scenario~name)+
    theme_bw()+
    ylab('% Variation')
  
  #média
  conv_media$Scenario=factor(conv_media$Scenario,
                             levels=unique(conv_media$Scenario),
                             labels=c('TA+, IL2-, PDL1/L2-, CD80/86-',
                                      'TA+, IL2+, PDL1/L2-, CD80/86-',
                                      'TA+, IL2+, PDL1/L2+, CD80/86+',
                                      'TA+, IL2+, PDL1/L2+, CD80/86-',
                                      'TA+, IL2+, PDL1/L2-, CD80/86+'))
  plot=conv_media%>%
    pivot_longer(!c(Scenario,Time,Trajectories))%>%
    ggplot()+
    geom_line(aes(Trajectories,value,col=factor(Time)))+
    facet_grid(name~Scenario)+theme_bw()+
    theme(legend.position='none')+
    ylab('Mean value')
  ggsave('Conv_media.pdf',plot,width=11,height=6)
  
  #Comparação
  plot=TOTAL%>%
    select(c(Scenario,N,Time,CFUNC,CINHI,TAPOP))%>%
    pivot_longer(!c(Scenario,N,Time))%>%
    mutate(Time=Time/N,
           k=qt(0.975,(N-1)),
           p=value+(k^2)/2,
           N=N+k^2,
           p=value/N,
           h=k*sqrt(p*(1-p)/N),
           upper=p+h,
           lower=p-h)%>%
    ggplot()+
    geom_ribbon(aes(Time,p,ymin=lower,ymax = upper),alpha=0.15)+
    geom_line(aes(Time,p))+
    geom_line(data=Original%>%
                select(c(Scenario,Time,CFUNC,CINHI,TAPOP))%>%
                pivot_longer(!c(Scenario,Time))%>%
                summarise(value=mean(value),
                          .by=c(Scenario,Time,name)),
              aes(Time,value),col='red')+
    facet_grid(Scenario~name)+theme_bw()+ylim(0,1)+
    ylab('% Activation')+xlab('Time (a.u.)')
  ggsave('IC.pdf',plot,width=11,height=6)
}

################### 
# Experiment BKMC #
###################
#Control
Original=c()
for (l in 1:nrow(Scenarios)) {
  t0=Sys.time()
  print(paste0('Current scenario: ',l))
  original=runExperiment(net,Time,R,p_a,p_d,Scenarios[l,],'BKMC')
  Original=rbind(Original,cbind('Scenario'=l,original))
  print(Sys.time()-t0)
}
Original$Scenario=factor(Original$Scenario,
                         levels=unique(Original$Scenario),
                         labels=c('TA+, IL2-, PDL1/L2-, CD80/86-',
                                  'TA+, IL2+, PDL1/L2-, CD80/86-',
                                  'TA+, IL2+, PDL1/L2+, CD80/86+',
                                  'TA+, IL2+, PDL1/L2+, CD80/86-',
                                  'TA+, IL2+, PDL1/L2-, CD80/86+'))

TOTAL=conv_media=conv_delta=tabela=c()
for (l in 1:nrow(Scenarios)) {
  for(r in 1:10){
    t0=Sys.time()
    scenario=Scenarios[l,]
    total=as.data.frame(matrix(0,nrow=num_points,ncol=length(net$genes)+1))
    names(total)=c('Time',net$genes)
    i=0 
    condition=FALSE
    repeat{
      buffer=runExperiment(net,Time,buff_size,p_a,p_d,scenario,'BKMC',points = Time)
      for (j in 1:buff_size) {
        print(paste0('Scenario ',l,' / Iteration ',buff_size*i+j))
        
        # Add new Trajectories 
        total=total+as.matrix(buffer[((j-1)*num_points+1):((num_points)*j),])
        
        # Stop condition
        y=total/(buff_size*i+j)
        curr=y%>%
          filter(Time==125|Time==250|Time==375)%>%
          select(Time,CFUNC,CINHI,TAPOP)
        
        if((buff_size*i+j)>1){
          delta=(curr-prev)/prev
          for(a in 1:nrow(delta)){
            for (b in 2:ncol(delta)) {
              if(abs(delta[a,b])>epsilon|is.na(delta[a,b]))break
            }
            if(abs(delta[a,b])>epsilon|is.na(delta[a,b]))break
          }
          if (a==nrow(delta)&b==ncol(delta)) condition=!condition
          conv_delta=rbind(conv_delta,
                           cbind('Scenario'=l,
                                 'Time'=c(125,250,375),
                                 'Trajectories'=(buff_size*i+j),
                                 delta[,-1]))
        }
        
        conv_media=rbind(conv_media,
                         cbind('Scenario'=l,
                               'Time'=c(125,250,375),
                               'Trajectories'=(buff_size*i+j),
                               curr[,-1]))
        
        if(condition | (buff_size*i+j)==max_iter) {
          break
        }
        
        prev=curr
      }
      if(condition | (buff_size*i+j)==max_iter) {
        break
      }
      i=i+1
    }
    #TOTAL=rbind(TOTAL,cbind('Scenario'=l,'N'=(buff_size*i+j),total))
    tabela=rbind(tabela,cbind(l,r,(buff_size*i+j),Sys.time()-t0))
    print(Sys.time()-t0) 
  }
}
TOTAL$Scenario=factor(TOTAL$Scenario,
                      levels=unique(TOTAL$Scenario),
                      labels=c('TA+, IL2-, PDL1/L2-, CD80/86-',
                               'TA+, IL2+, PDL1/L2-, CD80/86-',
                               'TA+, IL2+, PDL1/L2+, CD80/86+',
                               'TA+, IL2+, PDL1/L2+, CD80/86-',
                               'TA+, IL2+, PDL1/L2-, CD80/86+'))

tabela[,3]=round(tabela[,3]*60)
#Results
{
  setwd('C:/Users/Gabri/Desktop/Resultados/Convergence')
  if(!dir.exists('BKMC')){
    dir.create('BKMC')
  }
  setwd('BKMC')
  save(Original,file='Original.RData')
  save(TOTAL,file='Convergência.RData')
  save(conv_delta,file='Conv_delta.RData')
  save(conv_media,file='Conv_media.RData')
  save(tabela,file='Tabela.RData')
  
  print(tabela)
  
  if(!dir.exists('Figuras')){
    dir.create('Figuras')
  }
  setwd('Figuras')
  
  #delta
  conv_delta%>%
    pivot_longer(!c(Scenario,Time,Trajectories))%>%
    filter(name=='TAPOP',Time==125,Scenario==1)%>%
    ggplot()+
    geom_point(aes(Trajectories,abs(value),col=factor(Time)))+
    geom_line(data=cbind('x'=1:(1/epsilon),'y'=1/c(1:(1/epsilon))),aes(x,y))+
    facet_grid(Scenario~name)+
    theme_bw()+
    ylab('% Variation')
  
  #média
  plot=conv_media%>%
    pivot_longer(!c(Scenario,Time,Trajectories))%>%
    ggplot()+
    geom_line(aes(Trajectories,value,col=factor(Time)),linewidth=0.15)+
    facet_grid(name~Scenario)+theme_bw()+
    theme(legend.position='none')+
    ylab('Mean value')
  ggsave('Conv_média.pdf',plot,width = 6,height = 6)
  
  #Comparação
  plot=TOTAL%>%
    select(c(Scenario,N,Time,CFUNC,CINHI,TAPOP))%>%
    pivot_longer(!c(Scenario,N,Time))%>%
    mutate(Time=Time/N,
           k=qt(0.975,(N-1)),
           p=value+(k^2)/2,
           N=N+k^2,
           p=value/N,
           h=k*sqrt(p*(1-p)/N),
           upper=p+h,
           lower=p-h)%>%
    ggplot()+
    geom_ribbon(aes(Time,p,ymin=lower,ymax = upper),alpha=0.3)+
    geom_line(aes(Time,p),linewidth=0.3)+
    geom_line(data=Original%>%
                select(c(Scenario,Time,CFUNC,CINHI,TAPOP))%>%
                pivot_longer(!c(Scenario,Time))%>%
                summarise(value=mean(value),
                          .by = c(Scenario,Time,name)),
              aes(Time,value),col='red',linewidth=0.3)+
    facet_grid(name~Scenario)+theme_bw()+ylim(0,1)+
    ylab('Activity level')+xlab('Time (a.u.)')
  ggsave('IC.pdf',plot,width = 11,height = 6)
}