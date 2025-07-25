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
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

##############
# Parameters #
##############
net=loadNetwork("C:/Users/Gabri/Desktop/Networks/net_shah_3.1.txt")
R=rep(1000,5)#tabela[,2]
p_d=0.05
p_a=0.5
Times=rep(500,5)
# Times=Data%>%
#   summarise(x=max(Time),
#             .by=Scenario)%>%
#   select(x)
window=10
points=Times
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

getTrajectory=function(net,Time,rate_up,rate_down,initialState,algorithm){
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
    Trajectory=as.data.frame(Trajectory[-nrow(Trajectory),])
    Trajectory=rbind(Trajectory,c(Time,Trajectory[nrow(Trajectory),2]))
    Trajectory$t=as.numeric(Trajectory$t)
    names(Trajectory)=c('Time','State')
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

runBKMC=function(net,scenario,Time,R,rate_up,rate_down,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages = c('BoolNet','dplyr'),
            .combine = 'rbind',
            .export = c('getInitialState','getTrajectory')
  ) %dopar% {
    initState=getInitialState(net,scenario,rate_up,rate_down)
    getTrajectory(net,Time,rate_up,rate_down,initState,algorithm='BKMC')
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
    return(runBKMC(net,scenario,Time,R,rate_up,rate_down,threads))
  }
}

getTWCSS=function(Data){
  TWCSS=Data%>%
    pivot_longer(!Cluster)%>%
    mutate(centroid=mean(value),
           .by = c(Cluster,name))%>%
    mutate(d=(value-centroid)^2)%>%
    summarise(TWCSS=sum(d))
  return(TWCSS$TWCSS)
}

###################
# Experiment SDDS #
###################
Data=c()
for (l in 1:nrow(Scenarios)) {
  t=Sys.time()
  print(paste0('Current scenario: ',l))
  data=runExperiment(net,(Times[l]+window),R[l],p_a,p_d,Scenarios[l,],'SDDS')
  Data=rbind(Data,cbind('Scenario'=l,data%>%filter(Time>=Times[l])%>%select(!Time)))
  print(Sys.time()-t)
}
cl=makeCluster(10,type='SOCK')
registerDoSNOW(cl)
X=foreach(i=1:nrow(Data),
          .combine = 'rbind')%dopar%{
            cbind(Data[i,1],'State'=paste0(Data[i,-1],collapse = ''))
          }
stopCluster(cl)
X=as.data.frame(X)
names(X)=c('Scenario','State')

setwd('C:/Users/Gabri/Desktop/Resultados')
if(!dir.exists('Clustering')){
  dir.create('Clustering')
}
setwd('Clustering')
if(!dir.exists('SDDS')){
  dir.create('SDDS')
}
setwd('SDDS')

#save('TOTAL.RData',X)
load('TOTAL.RData')

l=1
#Piechart
{
  Y=X%>%filter(Scenario==l)%>%
    summarise(value=n(),.by = State)%>%
    mutate(percentage=100*value/sum(value))%>%
    arrange(desc(value),.by_group = TRUE)
  
  Z=c()
  for(i in 1:nrow(Y)){
    Z=rbind(Z,c('State'=i,'Value'=Y[i,2]))
  }
  
  ggplot(Z)+
    geom_bar(aes(x='',y=Value,fill=factor(State)),stat='identity')+
    coord_polar("y", start=0)+
    theme_void()+
    theme(legend.position = 'none')
  head(Y,n=10)
}

#Hierarchical clustering
{
  k_max=100
  #Number of clusters
  {
    #Organizing the data
    Z=X%>%
      filter(Scenario==l)%>%
      select(State)
    Y=as.data.frame(unique(Z))
    x=c()
    for(i in 1:nrow(Y)){
      x=rbind(x,as.vector(strsplit(Y[i,],'')[[1]],mode = 'double'))
    }
    x=as.data.frame(x)
    names(x)=net$genes
    
    #Calculating the distance matrix
    d=dist(x,method = 'binary')
    #Distance matrix visualization
    {
      # x11();fviz_dist(d,show_labels = FALSE,)+
      #   theme(legend.position = 'none')
    }
    
    #Tree construction
    tree=hclust(d)
    y=c()
    for(k in 1:k_max){
      x$Cluster=cutree(tree,k = k)
      y=rbind(y,c('k'=k,'TWCSS'=getTWCSS(x)))
    }
    plot(y,type='l')
    epsilon=0.01
    for(i in 3:(nrow(y)-2)){
      print(i)
      x_a=mean(c(y[i-2,2],y[i-1,2]))
      x_b=mean(c(y[i+2,2],y[i+1,2]))
      x_m=mean(c(y[i-2,2],y[i-1,2],y[i,2],y[i+2,2],y[i+1,2]))
      v=abs(x_b-x_a)/x_m
      if(v<epsilon) break
    }
    abline(v=i,lty='dashed',col='red')
  }
  
  #Dendrogram visualization
  {
    plot(color_branches(as.dendrogram(tree),k=i),
         ylim = c(0,1),leaflab = 'none',
         frame.plot = TRUE,edge.root = TRUE)
  }
  
  #Cluster's Information
  {
    #Abundance
    {
      x$Cluster=cutree(tree,k = i)
      cluster=Z
      for(j in 1:nrow(x)){
        cluster[which(Z$State==paste0(x[j,-ncol(x)],collapse = '')),1]=x$Cluster[j]
      }
      Z=cbind(Z,'Cluster'=cluster)
      Z=as.data.frame(Z)
      names(Z)=c('State','Cluster')
      Z$Cluster=as.numeric(Z$Cluster)
      Z%>%select(Cluster)%>%
        summarise(N=n(),.by = Cluster)%>%
        mutate('%'=100*N/sum(N))%>%
        arrange(desc(N),.by_group = TRUE)
      abundance=Z%>%
        ggplot()+
        geom_bar(aes(Cluster),
                 position = 'identity',
                 alpha=0.8)+
        theme_bw()+xlab('Cluster')+
        theme(axis.text.y = element_text(size=5),
              axis.text.x = element_text(size=7),
              axis.title.y = element_text(size=7),
              axis.title.x = element_text(size=7))+
        ggtitle('Cluser abundance')
    }
    
    #Heat map
    {
      b=x$Cluster[tree$order]
      a=x%>%
        pivot_longer(!Cluster)%>%
        summarise(N=n(),
                  value=mean(value),
                  .by = c(Cluster,name))%>%
        mutate(N=N/sum(unique(N)))
      a$name=factor(a$name,levels = net$genes)
      heatmap=ggplot(a)+
        geom_tile(aes(factor(Cluster,levels = unique(b)),name,
                      fill=value#,alpha=N
        ))+
        theme_bw()+xlab('Cluster')+
        guides(fill = guide_colourbar(barwidth = 7,
                                      barheight = 0.5,
                                      title = '',
                                      position = 'bottom',
                                      direction = 'horizontal'))+
        theme(axis.text.y = element_text(size=5),
              axis.text.x = element_text(size=7),
              axis.title.y = element_text(size=7),
              axis.title.x = element_text(size=7),
              legend.text = element_text(size=5))+
        ggtitle('Cluster info')
    }
    
    print(ggarrange(heatmap,abundance))
  }
}

#PCA
{
  data.pca=prcomp(na.omit(x[,-ncol(x)]))
  fviz_pca_biplot(data.pca,
                  geom.ind = element_blank(),
                  axes = c(1,2),
                  repel = TRUE)
}

###################
# Experiment BKMC #
###################
Data=c()
for (l in 1:nrow(Scenarios)) {
  t=Sys.time()
  print(paste0('Current scenario: ',l))
  data=runExperiment(net,(Times[1]+window),R,p_a,p_d,Scenarios[l,],'BKMC')
  Data=rbind(Data,cbind('Scenario'=l,data%>%filter(Time>=Times[1])%>%select(!Time)))
  print(Sys.time()-t)
}
X=Data

setwd('C:/Users/Gabri/Desktop/Resultados')
if(!dir.exists('Clustering')){
  dir.create('Clustering')
}
setwd('Clustering')
if(!dir.exists('BKMC')){
  dir.create('BKMC')
}
setwd('BKMC')

#save('TOTAL.RData',X)
load('TOTAL.RData')

l=3
#Piechart
{
  Y=X%>%filter(Scenario==l)%>%
    summarise(value=n(),.by = State)%>%
    mutate(percentage=100*value/sum(value))%>%
    arrange(desc(value),.by_group = TRUE)
  
  Z=c()
  for(i in 1:nrow(Y)){
    Z=rbind(Z,c('State'=i,'Value'=Y[i,2]))
  }
  ggplot(Z)+
    geom_bar(aes(x='',y=Value,fill=factor(State)),stat = 'identity')+
    coord_polar("y", start=0)+
    theme_void()+
    theme(legend.position = 'none')
  head(Y,n = 10)
}

#Hierarchical clustering
{
  k_max=100
  #Number of clusters
  {
    Z=X%>%filter(Scenario==l)%>%
      select(State)
    Y=as.data.frame(unique(Z))
    x=c()
    for(i in 1:nrow(Y)){
      x=rbind(x,as.vector(strsplit(Y[i,],'')[[1]],mode = 'double'))
    }
    x=as.data.frame(x)
    names(x)=net$genes
    d=dist(x,method = 'binary')
    #Distance matrix visualization
    {
      # x11();fviz_dist(d,show_labels = FALSE,)+
      #   theme(legend.position = 'none')
    }
    tree=hclust(d)
    y=c()
    for(k in 1:k_max){
      x$Cluster=cutree(tree,k = k)
      y=rbind(y,c('k'=k,'TWCSS'=getTWCSS(x)))
    }
    plot(y,type='l')
    epsilon=0.01
    for(i in 3:(nrow(y)-2)){
      print(i)
      x_a=mean(c(y[i-2,2],y[i-1,2]))
      x_b=mean(c(y[i+2,2],y[i+1,2]))
      x_m=mean(c(y[i-2,2],y[i-1,2],y[i,2],y[i+2,2],y[i+1,2]))
      v=abs(x_b-x_a)/x_m
      if(v<epsilon) break
    }
    abline(v=i,lty='dashed',col='red')
  }
  
  #Dendrogram visualization
  {
    plot(color_branches(as.dendrogram(tree),k=i),
         ylim = c(0,1),leaflab = 'none',
         frame.plot = TRUE,edge.root = TRUE)
  }
  
  #Cluster's Information
  {
    #Abundance
    x$Cluster=cutree(tree,k = i)
    cluster=Z
    for(i in 1:nrow(x)){
      cluster[which(Z$State==paste0(x[i,-ncol(x)],collapse = '')),1]=x$Cluster[i]
    }
    Z=cbind(Z,'Cluster'=cluster)
    Z=as.data.frame(Z)
    names(Z)=c('State','Cluster')
    Z$Cluster=as.numeric(Z$Cluster)
    #Z$Cluster=as.factor(Z$Cluster)
    Z%>%select(Cluster)%>%
      summarise(N=n(),.by = Cluster)%>%
      mutate('%'=100*N/sum(N))%>%
      arrange(desc(N),.by_group = TRUE)
    abundance=Z%>%
      ggplot()+
      geom_bar(aes(as.factor(Cluster)),
               position = 'identity',
               alpha=0.8)+
      theme_bw()+xlab('Cluster')+
      theme(axis.text.y = element_text(size=5),
            axis.text.x = element_text(size=7),
            axis.title.y = element_text(size=7),
            axis.title.x = element_text(size=7))+
      ggtitle('Cluser abundance')
    
    #Heat map
    b=x$Cluster[tree$order]
    a=x%>%
      pivot_longer(!Cluster)%>%
      summarise(N=n(),
                value=mean(value),
                .by = c(Cluster,name))%>%
      mutate(N=N/sum(unique(N)))
    a$name=factor(a$name,levels = net$genes)
    heatmap=ggplot(a)+
      geom_tile(aes(factor(Cluster,levels = unique(b)),name,
                    fill=value#,alpha=N
      ))+
      theme_bw()+xlab('Cluster')+
      guides(fill = guide_colourbar(barwidth = 7,
                                    barheight = 0.5,
                                    title = '',
                                    position = 'bottom',
                                    direction = 'horizontal'))+
      theme(axis.text.y = element_text(size=5),
            axis.text.x = element_text(size=7),
            axis.title.y = element_text(size=7),
            axis.title.x = element_text(size=7),
            legend.text = element_text(size=5))+
      ggtitle('Cluster info')
    
    print(ggarrange(heatmap,abundance))
  }
}

#PCA
{
  data.pca=prcomp(na.omit(x[,-ncol(x)]))
  fviz_pca_biplot(data.pca,
                  geom.ind = element_blank(),
                  axes = c(1,2),
                  repel = TRUE)
}