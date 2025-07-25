############
# Packages #
############
PKG=c('BoolNet','tidyverse','gpuR',
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

getInitialState=function(net,scenario){
  State=generateState(net,specs = c('LCK'=1))
  State['TAex']=ifelse(scenario['TAex'],1,0)
  State['IL2ex']=ifelse(scenario['IL2ex'],1,0)
  State['PDL1L2ex']=ifelse(scenario['PDL1L2ex'],1,0)
  State['CD8086ex']=ifelse(scenario['CD8086ex'],1,0)
  return(State)
}

getTrajectory=function(net,Time,rate_up,rate_down,initialState){
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
  return(Trajectory)
}

runSDDS=function(net,scenario,Time,R,rate_up,rate_down,threads){
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages=c('BoolNet'),
            .combine = '+',
            .export = c('getInitialState','getTrajectory')
  ) %dopar% {
    as.matrix(getTrajectory(net,Time,rate_up,rate_down,getInitialState(net,scenario)))
  }
  stopCluster(cl)
  X=as.data.frame(X/R)
  return(X)
}

runExperiment=function(net,Time,R,p_a,p_d,scenario,rate_TA,rate_scFV,threads=12){
  rate_up=getRates(net,
                   c('TAex'=rate_TA,
                     'IL2ex'=ifelse(scenario['IL2ex'],p_a,0),
                     'CD8086ex'=ifelse(scenario['CD8086ex'],p_a,0),
                     'PDL1L2ex'=ifelse(scenario['PDL1L2ex'],p_a,0),
                     'SCFV'=rate_scFV),
                   default = p_a)
  rate_down=getRates(net,
                     c('TAex'=ifelse(scenario['TAex'],p_d,0),
                       'IL2ex'=ifelse(scenario['IL2ex'],p_d,0),
                       'CD8086ex'=ifelse(scenario['CD8086ex'],p_d,0),
                       'PDL1L2ex'=ifelse(scenario['PDL1L2ex'],p_d,0)),
                     default = p_d)
  return(runSDDS(net,scenario,Time,R,rate_up,rate_down,threads))
}

##############
# Parameters #
##############
R=1000
p_d=0.05
p_a=0.5 
Time=500
Scenarios=data.frame(TAex=c(1,1,1,1,1),
                     IL2ex=c(0,1,1,1,1),
                     PDL1L2ex=c(0,0,1,1,0),
                     CD8086ex=c(0,0,1,0,1))
Rates_TA=c(seq(0,0.1,0.0025),seq(0.105,0.2,0.005))
Rates_scFV=c(seq(0,0.1,0.0025),seq(0.105,0.2,0.005))

##############
# Experiment #
##############
Data=c()
net=loadNetwork(paste0('Networks/net_shah_3.1.txt'))
for(rate_scFV in Rates_scFV){
  for(rate_TA in Rates_TA){
    t0=Sys.time()
    print(paste0('Current Rate TA: ',rate_TA,' Current Rate scFV: ',rate_scFV))
    data=runExperiment(net,Time,R,p_a,p_d,Scenarios[1,],rate_TA,rate_scFV)
    data=cbind('Rate_TA'=rate_TA,'Rate_scFV'=rate_scFV,data)
    Data=rbind(Data,data)
    print(Sys.time()-t0)
  }
}
rio::export(as.data.frame(Data[,-2]),file='RateVariationSDDS.RData')
x11();Data[,-2]%>%
  pivot_longer(!c(Time,rate_scFV,Rate_TA))%>%
  filter(Time==max(Time),name=='CFUNC')%>%
  summarise(value=mean(value),
            .by=c(Rate_TA,rate_scFV,name))%>%
  ggplot()+
  geom_tile(aes(rate_scFV,Rate_TA,fill=value))+
  theme_bw()+
  theme(legend.position = 'none')
  

#########
# Plots #
#########
setwd("C:/Users/Gabri/Desktop/Resultados/RateVariation")
if (!dir.exists('SDDS')){
  dir.create('SDDS')
}
setwd('SDDS')
save(Data,file = 'Rate_Variation.RData')
if (!dir.exists('Simulações')){
  dir.create('Simulações')
}
setwd('Simulações')

plot=Data%>%
  select(Rate_TA,rate_scFV,Time,TAex,TA,CFUNC,CINHI,TAPOP)%>%
  filter(Time==max(Time))%>%
  select(!Time)%>%
  pivot_longer(!c(Rate_TA,rate_scFV))%>%
  ggplot()+
  geom_line(aes(Rate,value,col=name))+
  theme_bw()
ggsave('rate_vs_final.pdf',plot,width=4,height = 3)

Data$Rate=as.factor(Data$Rate)

scale_range = colorRampPalette(c("#e9e9e9", "#514e4c"))
color_scale = scale_range(length(Rates)) 

for(i in 1:length(net$genes)){
  print(net$genes[i])
  X=Data%>%
    select(Rate,Scenario,Time,i+3)%>%
    pivot_longer(!c(Rate,Scenario,Time))%>%
    summarise(`% Activation`=mean(value),
              .by = c(Rate,Scenario,Time,name))%>%
    ggplot()+geom_line(aes(Time,`% Activation`,col=Rate))+
    ylim(0,1.1)+scale_color_manual(values=setNames(color_scale, levels(as.factor(Rates))))+
    facet_wrap(~Scenario,nrow=1)+theme_classic()+ggtitle(paste0(net$genes[i]))
  ggsave(paste(net$genes[i],'.jpg'),X,width = 20,height = 7)
}

plot=Data%>%
  select(Rate,Time,
         TAex,TA,CFUNC,CINHI,TAPOP)%>%
  pivot_longer(!c(Rate,Time))%>%
  mutate(name=factor(name,
                     levels=c('TAex','TA','CFUNC','CINHI','TAPOP')))%>%
  ggplot()+geom_line(aes(Time,value,col=as.factor(Rate)))+
  ylim(0,1.1)+scale_color_manual(values=setNames(color_scale,levels(as.factor(Rates))))+
  facet_wrap(~name,nrow = 1)+theme_classic()+
  ylab('% Activation')+xlab('Time (a.u.)')+
  theme(legend.key.size = unit(0.4,'cm'),
        legend.text = element_text(size=7))
ggsave('Plot_rates.pdf',plot,width=12,height=7)
