############
# Packages #
############
PKG=c('BoolNet','tidyverse','ggpubr',
      'doSNOW','foreach','factoextra',
      'ape','dendextend','ggdendro')
for (pkg in PKG) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

##############
# Parameters #
##############
data("cellcycle")
x11();plotAttractors(getAttractors(cellcycle,'asynchronous'),
               allInOnePlot = TRUE,
               drawLegend = FALSE)

getAttractors(cellcycle,'asynchronous')
png("Basin2.png", width = 1000, height = 1000)
plotStateGraph(
  getBasinOfAttraction(getAttractors(cellcycle),2),
  drawLegend = FALSE,
  colorsAlpha = c(colorBasinsNodeAlpha    = .4,
                  colorBasinsEdgeAlpha    = .4,
                  colorAttractorNodeAlpha = 1,
                  colorAttractorEdgeAlpha = 1),
  layout = layout.kamada.kawai,
  basin.lty = 2,attractor.lty = 1)
dev.off()

p_d=0.05
p_a=0.5
Time=1000

#############
# Functions #
#############
getRates=function(net,RateSpecs=NULL,default=0){
  n=length(net$genes)
  rate=vector(mode='numeric',length = n)
  names(rate)=net$genes
  for (i in 1:n) {
    rate[i]=default
  }
  if(!is.null(RateSpecs)){
    for(i in names(RateSpecs)){
      rate[which(names(rate)==i)]=RateSpecs[i]
    } 
  }
  return(rate)
}

getInitialState=function(net,InitSpecs=NULL){
  state=generateState(net,specs = c('CycD'=1))
  for(i in 1:length(state)){
    state[i]=rbinom(1,1,0.5)
  }
  if(!is.null(InitSpecs)){
    for(i in names(InitSpecs)){
      state[i]=InitSpecs[i]
    }
  }
  return(state)
}

getTrajectory=function(net,Time,rate_up,rate_down,initialState,algorithm,num_points){
  if(algorithm=='SDDS'){
    #Initialization
    {
      State=initialState
      n=length(net$genes)
      Trajectory=c()
      count=rep(0,n)
      names(count)=names(State) 
    }
    
    #Iterative process
    for (t in 1:Time) {
      s=sample(1:n,n)
      for (i in s) {
        count[i]=count[i]+1
        if (State[i]) { #1
          if (runif(1,0,1)<rate_down[i]&count[i]>1) {
            State[i]=stateTransition(net,State,'asynchronous',chosenGene = i)[i]
            count[i]=0
          } else {
            #State[i]=State[i]
          }
        } else { #0
          if (runif(1,0,1)<rate_up[i]&count[i]>1) {
            State[i]=stateTransition(net,State,'asynchronous',chosenGene = i)[i]
            count[i]=0
          } else {
            #State[i]=State[i]
          }
        }
      }
      Trajectory=rbind(Trajectory,c(t,State))
    }
    
    #Consolidation
    {
      Trajectory=as.data.frame(Trajectory)
      names(Trajectory)=c('Time',net$genes) 
    }
  }else if(algorithm=='BKMC'){
    #Initialization
    {
      State=paste0(initialState,collapse='')
      n=length(net$genes)
      t=0
      Trajectory=cbind(t,State) 
    }
    
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
        Trajectory=rbind(Trajectory,c(Time,State))
        break
      }
    }
    
    #Consolidation
    {
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
  }
  return(Trajectory)
}

runSDDS=function(net,Time,R,rate_up,rate_down,InitSpecs,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages = 'BoolNet',
            .combine = 'rbind',
            .export = c('getInitialState','getTrajectory')) %dopar% {
              initState=getInitialState(net,InitSpecs=InitSpecs)
              x=getTrajectory(net,Time,rate_up,rate_down,initState,algorithm='SDDS')
              x
            }
  stopCluster(cl)
  return(X)
}

runBKMC=function(net,Time,R,rate_up,rate_down,InitSpecs,points,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages = c('BoolNet','dplyr'),
            .combine = 'rbind',
            .export = c('getInitialState','getTrajectory')
  ) %dopar% {
    initState=getInitialState(net,InitSpecs=InitSpecs)
    getTrajectory(net,Time,rate_up,rate_down,initState,algorithm='BKMC',num_points=points)
  }
  stopCluster(cl)
  X=as.data.frame(X)
  return(X)
}

runExperiment=function(net,Time,R,p_a,p_d,algorithm,
                       InitSpecs=NULL,RateSpecs=NULL,
                       points=100,threads=10){
  
  rate_up=getRates(net,RateSpecs=RateSpecs,default = p_a)
  
  rate_down=getRates(net,RateSpecs=RateSpecs,default = p_d)
  
  if(algorithm=='SDDS'){
    return(runSDDS(net,Time,R,rate_up,rate_down,InitSpecs,threads))
  }else if(algorithm=='BKMC'){
    return(runBKMC(net,Time,R,rate_up,rate_down,InitSpecs,points,threads))
  }
}
##############
# Experiment #
##############
Data=c()
for (i in c('SDDS','BKMC')) {
  for(j in c(0,1)){
    t0=Sys.time()
    print(c(i,j))
    data=runExperiment(net=cellcycle,
                       Time=Time,
                       R=10000,
                       p_a=p_a,
                       p_d=p_d,
                       algorithm = i,
                       points = Time,
                       InitSpecs = c('CycD'=j))
    Data=rbind(Data,cbind('Method'=i,'IC'=j,data))
    print(Sys.time()-t0)
  }
}
Data$IC=factor(Data$IC,
                     levels=unique(Data$IC),
                     labels=c('Absence of growth factors (CycD=0)',
                              'Presence of growth factors (CycD=1)'))
setwd('Resultados/Cellcycle')
save(Data,file = 'Cellcycle.RData')
fig=Data%>%
  filter(Method=='SDDS')%>%
  pivot_longer(!c(Time,Method,IC))%>%
  summarise(value=mean(value),.by=c(Time,name,Method,IC))%>%
  mutate(Node=factor(name,cellcycle$genes))%>%
  select(!c(name,Method))%>%
  ggplot()+geom_line(aes(Time,value,col=Node),linewidth=0.35)+
  theme_bw()+labs(y='Activity level',x='Time (a.u.)')+
  facet_wrap(~IC)+
  theme(legend.title = element_blank())
ggsave(filename = 'SDDS_cellcycle.pdf',
       plot = fig,
       width = 7,
       height = 3)

################
# Steady state #
################
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
            State[i]=stateTransition(net,State,'asynchronous',chosenGene = i)[i]
            count[i]=0
          } else {
            #State[i]=State[i]
          }
        } else { #0
          if (runif(1,0,1)<rate_up[i]&count[i]>1) {
            State[i]=stateTransition(net,State,'asynchronous',chosenGene = i)[i]
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
        Trajectory=rbind(Trajectory,c(t0+Time_buffer,State))
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

runSDDS=function(net,t0,Time_buffer,R,rate_up,rate_down,initialStates,InitSpecs,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages = c('BoolNet'),
            .combine = 'rbind',
            .export = c('getInitialState','getTrajectory')
  ) %dopar% {
    if(is.null(initialStates)){
      initState=getInitialState(net,InitSpecs=InitSpecs)
    }else{
      initState=as.vector(initialStates[r,],mode='double')
    }
    getTrajectory(net,t0,Time_buffer,rate_up,rate_down,initState,
                  algorithm='SDDS')
  }
  stopCluster(cl)
  return(X)
}

runBKMC=function(net,t0,Time_buffer,R,rate_up,rate_down,initialStates,InitSpecs,points,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages = c('BoolNet','dplyr'),
            .combine = 'rbind',
            .export = c('getInitialState','getTrajectory')
  ) %dopar% {
    if(is.null(initialStates)){
      initState=getInitialState(net,InitSpecs=InitSpecs)
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

runExperiment=function(net,t0,Time_buffer,R,p_a,p_d,algorithm,
                       initialStates=NULL,InitSpecs=NULL,RateSpecs=NULL,
                       points=100,threads=10){
  
  rate_up=getRates(net,RateSpecs=RateSpecs,default = p_a)
  
  rate_down=getRates(net,RateSpecs=RateSpecs,default = p_d)
  
  if(algorithm=='SDDS'){
    return(runSDDS(net,t0,Time_buffer,R,rate_up,rate_down,initialStates,InitSpecs,threads))
  }else if(algorithm=='BKMC'){
    return(runBKMC(net,t0,Time_buffer,R,rate_up,rate_down,initialStates,InitSpecs,points,threads))
  }
}

#Parameters
R=10000
Time_buffer=200
Time_window=50
w1=0.1
w2=0.5
eta=0.01

#Experiment
Data=FinalStates=c()
t=Sys.time()
t0=0
data=c()
condition=FALSE
repeat{
  if(!t0){
    buffer=runExperiment(net,t0,Time_buffer,R,p_a,p_d,'BKMC')
  }else{
    buffer=runExperiment(net,t0,Time_buffer,R,p_a,p_d,'BKMC',initialStates = states)
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
    ggplot()+geom_line(aes(Time,value,col=name))+ylim(0,1)+
    theme_bw()
  print(p)
  
  #Sliding window stop condition
  {
    for(i in Time_window:(t0+Time_buffer)){
      X=data%>%filter(Time<i,Time>(i-Time_window))
      x=X%>%
        select(Time,CycA,E2F,UbcH10)%>%
        pivot_longer(!Time)%>%
        summarise(p=mean(value),.by = c(Time,name))%>%
        summarise(q=mean(p),.by = name)
      a=X%>%
        select(Time,CycA,E2F,UbcH10)%>%
        filter(Time<=(min(Time)+w1*Time_window))%>%
        pivot_longer(!Time)%>%
        summarise(p=mean(value),.by = c(Time,name))%>%
        summarise(q=mean(p),.by = name)
      b=X%>%
        select(Time,CycA,E2F,UbcH10)%>%
        filter(Time>=(min(Time)+(1-w2)*Time_window))%>%
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
Data=rbind(Data,data%>%filter(Time<=i))
FinalStates=rbind(FinalStates,cbind(filter(buffer,Time>=(i-5),Time<=i)))
print(Sys.time()-t)
Time=i

#Clustering
getTWCSS=function(Data){
  TWCSS=Data%>%
    pivot_longer(!Cluster)%>%
    mutate(centroid=mean(value),
           .by = c(Cluster,name))%>%
    mutate(d=(value-centroid)^2)%>%
    summarise(TWCSS=sum(d))
  return(TWCSS$TWCSS)
}

#Analysis
cl=makeCluster(12,type='SOCK')
registerDoSNOW(cl)
X=foreach(i=1:nrow(FinalStates),
          .combine = 'rbind')%dopar%{
            cbind(FinalStates[i,1],'State'=paste0(round(FinalStates[i,-1],0),collapse = ''))
          }
stopCluster(cl)

#Piechart
{
  Y=X%>%as.data.frame()%>%select(State)%>%
    summarise(count=n(),.by = State)%>%
    mutate(percentage=100*count/sum(count))%>%
    arrange(desc(count),.by_group = TRUE)
  Z=c()
  for(i in 1:(nrow(Y))){
    d=nrow(getPathToAttractor(net,as.vector(strsplit(Y$State[i],'')[[1]],mode = 'double')))
    Z=rbind(Z,c('State'=i,'Value'=Y[i,2],'Percentage'=Y[i,3],'Dist'=d))
  }
  ggplot(Z)+
    geom_bar(aes(x='',y=Percentage,fill=factor(State)),stat = 'identity')+
    coord_polar("y", start=0)+
    theme_void()+
    theme(legend.position = 'none')
  head(Y,n = 10)
}

#Hierarchical clustering
{
  #Number of clusters
  {
    Y=X%>%as.data.frame()%>%select(State)%>%
      summarise(count=n(),.by = State)%>%
      mutate(percentage=100*count/sum(count))%>%
      arrange(desc(count),.by_group = TRUE)
    k_max=nrow(Y)
    x=c()
    for(i in 1:nrow(Y)){
      x=rbind(x,as.vector(strsplit(Y$State[i],'')[[1]],mode = 'double'))
    }
    x=as.data.frame(x)
    names(x)=net$genes
    d=dist(x,method = 'binary')
    tree=hclust(d)
    y=c()
    for(k in 1:k_max){
      x$Cluster=cutree(tree,k = k)
      y=rbind(y,c('k'=k,'TWCSS'=getTWCSS(x)))
    }
    plot(y,type='l')

    epsilon=0.025
    for(i in 2:nrow(y)){
      v=abs(y[i,2]-y[i-1,2])/y[i-1,2]
      print(i)
      if(v<epsilon) break
    }
    abline(v=i,lty='dashed',col='red')
  }
  
  #Dendrogram visualization
  {
    plot(color_branches(as.dendrogram(tree),k=i),
         ylim = c(0,1),leaflab = 'none',
         frame.plot = TRUE,edge.root = TRUE)
    
    #ggplot radial approach
    {
      # Convert to dendrogram object
      dend = as.dendrogram(tree)
      
      # Convert dendrogram to a format for ggplot2
      dendro_data = dendro_data(dend, type = "rectangle")
      segment_data=segment(dendro_data)
      
      # Create a circular plot using ggplot
      print(ggplot(segment_data) +
              geom_segment(aes(x=x, y=-y, xend=xend, yend=-yend),
                           linewidth = 0.1) +
              geom_hline(yintercept =-0.75,
                         col='red',alpha=0.5,
                         linetype='dashed')+
              coord_radial(theta="x",
                           expand = TRUE,
                           rotate_angle = TRUE)+
              theme_void())
    }
  }
  
  #Cluster's Information
  {
    #Size
    x$Cluster=cutree(tree,k = i)
    x$N=Y$count
    bar=x%>%ggplot()+
      geom_bar(aes(factor(Cluster)),
               position = 'identity',alpha=0.8)+
      theme_bw()+xlab('Cluster')+ylab('N° of states')+
      theme(axis.text.y = element_text(size=5),
            axis.text.x = element_text(size=7),
            axis.title.y = element_text(size=7),
            axis.title.x = element_text(size=7))
    
    #Composition
    x%>%summarise(S=sum(N),.by = c(Cluster))
    x%>%filter(Cluster==1)
    a=x%>%
      pivot_longer(!c(Cluster,N))%>%
      summarise(count=n(),
                value=mean(value),
                .by = c(Cluster,N,name))
    a$name=factor(a$name,levels = net$genes)
    heatmap=ggplot(a)+
      geom_tile(aes(factor(Cluster),name,
                    fill=value,alpha = N))+
      theme_bw()+xlab('Cluster')+
      guides(fill = guide_colourbar(barwidth = 0.5,
                                    barheight = 7,
                                    title = '',
                                    position = 'right',
                                    direction = 'vertical',
                                    cex=0.5))+
      theme(axis.text.y = element_text(size=5),
            axis.text.x = element_text(size=7),
            axis.title.y = element_text(size=7),
            axis.title.x = element_text(size=7),
            legend.text = element_text(size=5))
    
    print(ggarrange(bar,heatmap))
  }
  
  #Distance matrix visualization
  {
    # x11();fviz_dist(d,show_labels = FALSE,)+
    #   theme(legend.position = 'none')
  }
}

#K-means
{
  #set.seed(1)
  X=FinalStates%>%select(!Time)
  x=c()
  condition=FALSE
  for(i in 1:1000){
    #set.seed(1)
    print(i)
    modelo=kmeans(X,centers = i)
    x=rbind(x,c(i,modelo$tot.withinss))
    if(i>1){
      condition=(abs((x[i,2]-x[(i-1),2])/x[(i-1),2])<0.05)
    }
    if(condition) break
  }
  X=cbind(X,'Cluster'=modelo$cluster)
  plot(x,pch=20)
  y=c()
  for(i in 1:nrow(X)){
    y=rbind(y,cbind(X[i,11],paste0(X[i,-11],collapse = '')))
  }
  
  a=as.data.frame(unique(y))
  a$V1=as.factor(a$V1)
  names(a)=c('Cluster','States')
  
  a=modelo$centers
  a=cbind('Center'=1:nrow(a),a)
  a=as.data.frame(a)
  a=a%>%pivot_longer(!Center)
  a$name=factor(a$name,levels = net$genes)  
  print(ggplot(a)+geom_tile(aes(Center,name,fill=value)))
  
  y=c()
  for(i in 1:length(net$genes)){
    m=anova(lm(data=X,formula = Cluster~X[,i]))
    y=rbind(y,c('Gene'=net$genes[i],'F_stat'=m$`F value`[1]))
  }
  y=as.data.frame(y)
  y$F_stat=round(as.numeric(y$F_stat),2)
  y%>%arrange(desc(F_stat))
}

#PCA
{
  data.pca=prcomp(na.omit(x[,-c(ncol(x),ncol(x)-1)]))
  fviz_pca_biplot(data.pca,geom.ind = element_blank(),
                  axes = c(1,2),
                  repel = TRUE) 
}

x%>%filter(Cluster==4)%>%
  pivot_longer(!Cluster)%>%
  summarise(value=round(mean(value),0),.by = name)%>%
  filter(name!='N')%>%
  pivot_wider(names_from = name,values_from = value)
getAttractors(net)

###############
# Convergence #
###############
Original=runExperiment(net,Time,10000,p_a,p_d,'SDDS')

epsilon=0.001
max_iter=10000
buff_size=200
n_points=5
points=c()
for(i in 1:n_points){
  points=rbind(points,round(i*Time/n_points-Time/(2*n_points),0))
}

t0=Sys.time()
TOTAL=conv_media=conv_delta=tabela=c()
total=as.data.frame(matrix(0,nrow=Time,ncol=length(net$genes)+1))
names(total)=c('Time',net$genes)
i=0 
condition=FALSE
repeat{
  buffer=runExperiment(net,Time,buff_size,p_a,p_d,'SDDS')
  for (j in 1:buff_size) {
    print(paste0('Iteration ',buff_size*i+j))
    
    # Add new Trajectories 
    total=total+as.matrix(buffer[(Time*(j-1)+1):(Time*j),])
    
    # Stop condition
    y=total/(buff_size*i+j)
    curr=c()
    for(k in 1:nrow(points)){
      curr=rbind(curr,y%>%
                   filter(Time==points[k,])%>%
                   select(Time,E2F,CycA,UbcH10))
    }
    
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
                       cbind('Time'=unique(points),
                             'Trajectories'=(buff_size*i+j),
                             delta[,-1]))
    }
    
    conv_media=rbind(conv_media,
                     cbind('Time'=unique(points),
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
TOTAL=rbind(TOTAL,cbind('N'=(buff_size*i+j),total))
print(Sys.time()-t0)

conv_delta%>%
  pivot_longer(!c(Time,Trajectories))%>%
  ggplot()+
  geom_point(aes(Trajectories,abs(value),col=factor(Time)))+
  geom_line(data=cbind('x'=1:(1/epsilon),'y'=1/c(1:(1/epsilon))),aes(x,y))+
  facet_wrap(~name)+
  theme_bw()+
  ylab('% Variation')

conv_media%>%
  pivot_longer(!c(Time,Trajectories))%>%
  ggplot()+
  geom_line(aes(Trajectories,value,col=factor(Time)))+
  facet_grid(Time~name)+theme_bw()+
  theme(legend.position='none')+
  ylab('% Activation')

p=TOTAL%>%
  pivot_longer(!c(N,Time))%>%
  mutate(Time=Time/N,
         k=qt(0.975,(N-1)),
         p=value+(k^2)/2,
         N=N+k^2,
         p=value/N,
         h=k*sqrt(p*(1-p)/N),
         upper=p+h,
         lower=p-h)%>%
  ggplot()+
  geom_ribbon(aes(Time,p,ymin=lower,ymax=upper),alpha=0.15)+
  geom_line(aes(Time,p))+
  geom_line(data=Original%>%
              pivot_longer(!c(Time))%>%
              summarise(value=mean(value),.by = c(Time,name)),
            aes(Time,value),col='red')+
  facet_wrap(~name)+theme_bw()+ylim(0,1)+
  ylab('% Activation')+xlab('Time (a.u.)')
for(j in 1:nrow(points)){
  p=p+geom_vline(xintercept = points[j,],linetype='dotted')
}
p