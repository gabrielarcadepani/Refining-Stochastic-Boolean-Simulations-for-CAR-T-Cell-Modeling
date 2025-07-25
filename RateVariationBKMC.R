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

getInitialState=function(net,scenario){
  State=generateState(net,specs = c('LCK'=1))
  State['TAex']=ifelse(scenario['TAex'],1,0)
  State['IL2ex']=ifelse(scenario['IL2ex'],1,0)
  State['PDL1L2ex']=ifelse(scenario['PDL1L2ex'],1,0)
  State['CD8086ex']=ifelse(scenario['CD8086ex'],1,0)
  return(State)
}

getTrajectory=function(net,Time,rate_up,rate_down,initialState){
  State=paste0(initialState,collapse='')
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
    
    #manually include input activation/deactivation as possible transitions
    {
      x=as.numeric((strsplit(State,'')[[1]]))
      if(x[1]){
        x[1]=0
        possible_states=rbind(possible_states,c(paste0(x,collapse = ''),
                                                rate_down['TAex']))
      }else{
        x[1]=1
        possible_states=rbind(possible_states,c(paste0(x,collapse = ''),
                                                rate_up['TAex']))
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
  return(Trajectory)
}

runBKMC=function(net,scenario,Time,R,rate_up,rate_down,threads){
  time_tick=1
  cl=makeCluster(threads,type='SOCK')
  registerDoSNOW(cl)
  X=foreach(r=1:R,
            .packages=c('BoolNet','dplyr'),
            .combine = '+',
            .export = c('getInitialState','getTrajectory')) %dopar% {
              initState=getInitialState(net,scenario)
              data=getTrajectory(net,Time,rate_up,rate_down,initState)
              #Consolidation
              {
                data=rbind(c('0',paste0(initState,collapse = '')),data)
                data=as.data.frame(data)
                data[,1]=as.numeric(data[,1])
                names(data)=c('Time','State')
                n=nrow(data)
                bkmc=matrix(0,ncol = length(net$genes),nrow =Time)
                for(i in 1:(n-1)){
                  p=data$Time[i]%/%time_tick+1
                  q=data$Time[i+1]%/%time_tick+1
                  if(p==q){ #same time_tick
                    delta_t=(data$Time[i+1]-data$Time[i])
                    bkmc[p,]=bkmc[p,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]))
                  }else if(q-p==1){ #consecutive time_ticks
                    delta_t=(p-data$Time[i])
                    bkmc[p,]=bkmc[p,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]))
                    if(q>Time)break
                    delta_t=(data$Time[i+1]-p)
                    bkmc[q,]=bkmc[q,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]))
                  }else{ #spread time_ticks
                    delta_t=(p-data$Time[i])
                    bkmc[p,]=bkmc[p,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]))
                    for(j in (p+1):(q-1)){
                      if(j>Time)break
                      bkmc[j,]=bkmc[j,]+c(as.numeric(strsplit(data$State[i],'')[[1]]))
                    }
                    if(q>Time)break
                    delta_t=(data$Time[i+1]-j)
                    bkmc[q,]=bkmc[q,]+c(delta_t*as.numeric(strsplit(data$State[i],'')[[1]]))
                  }
                }
                bkmc[nrow(bkmc),]=as.numeric(strsplit(data$State[nrow(data)-1],split ='')[[1]])
              }
              bkmc
            }
  stopCluster(cl)
  X=as.data.frame(X/R)
  names(X)=net$genes
  X=cbind('Time'=1:Time,X)
  return(X)
}

runExperiment=function(net,Time,R,p_a,p_d,scenario,rate,threads=10){
  rate_up=getRates(net,
                   c('TAex'=rate,
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
  
  return(runBKMC(net,scenario,Time,R,rate_up,rate_down,threads))
}

##############  
# Parameters #
##############
load('Resultados/Convergence/BKMC/Tabela.RData')
R=tabela[,2]
p_d=0.05
p_a=0.5 
load('Resultados/SteadyState/BKMC/Tabela.RData')
Time=t(Tabela[2])
Scenarios=data.frame(TAex=c(1,1,1,1,1),
                     IL2ex=c(0,1,1,1,1),
                     PDL1L2ex=c(0,0,1,1,0),
                     CD8086ex=c(0,0,1,0,1))
Rates=c(seq(0,0.01,0.001),0.015,
        seq(0.02,0.1,0.01),0.15,
        seq(0.2,0.5,0.1))

##############
# Experiment #
##############
Data=c()
net=loadNetwork('Networks/net_shah_3.1.txt')
for(l in 1:nrow(Scenarios)){
  for(rate in Rates){
    t0=Sys.time()
    print(paste0('Current scenario: ',l,' / Current Rate: ',rate))
    data=runExperiment(net,Time[l],R[l],p_a,p_d,Scenarios[l,],rate)
    data=cbind('Scenario'=l,'Rate'=rate,'R'=R[l],'Max_t'=Time[l],data)
    Data=rbind(Data,data)
    print(Sys.time()-t0)
  }
}
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
setwd("Resultados/RateVariation")
if (!dir.exists('BKMC')){
  dir.create('BKMC')
}
setwd('BKMC')
save(Data,file = 'Rate_Variation.RData')
if (!dir.exists('Simulações')){
  dir.create('Simulações')
}
setwd('Simulações')

#Comparação
plot=Data%>%
  filter(Scenario=='TA+, IL2-, PDL1/L2-, CD80/86-')%>%
  filter(Time==max(Time),.by = Scenario)%>%
  dplyr::select(c(Scenario,Rate,R,TAex,TA,CFUNC,CINHI,TAPOP))%>%
  pivot_longer(!c(Scenario,Rate,R))%>%
  mutate(k=qt(0.975,(R-1)),
         p=value+(k^2)/2,
         N=R+k^2,
         p=value,
         h=k*sqrt(p*(1-p)/R),
         upper=p+h,
         lower=p-h,
         name=factor(name,levels=unique(name)))%>%
  ggplot()+
  geom_ribbon(aes(Rate,p,ymin=lower,ymax = upper,fill=name),alpha=0.2)+
  #geom_errorbar(aes(Rate,p,ymin=lower,ymax = upper,colour =name),alpha=0.4)+
  geom_line(aes(Rate,p,col=name),linewidth=0.5)+
  #facet_wrap(~Scenario,nrow = 1)+
  theme_bw()+ylim(0,1)+
  ylab('Activity level')+
  xlab('Rate')+
  theme(legend.position = 'right',
        legend.title = element_blank(),
        legend.key.size = unit(0.5,'cm'),
        legend.text = element_text(size=9))
ggsave('rate_vs_final.pdf',plot,width=7,height=3)

x=Data%>%
  filter(Scenario=='TA+, IL2+, PDL1/L2-, CD80/86-')%>%
  filter(Time==max(Time),.by = Scenario)%>%
  dplyr::select(c(Scenario,Rate,TAex,TA,CFUNC,CINHI,TAPOP))%>%
  pivot_longer(-c(Scenario,Rate))%>%
  dplyr::select(Scenario,name,Rate,value)

y=filter(x,name=='CINHI')%>%
  dplyr::select(-c(name,Scenario))%>%
  #spline()%>%
  as.data.frame()
names(y)=c('Rate','value')
model=drc::drm(value~Rate,data=y,fct = drc::MM.2())
summary(model)
y$pred=predict(model,as.data.frame(y$Rate))
ggplot(y)+
  geom_point(aes(Rate,value))+
  geom_line(aes(Rate,pred),col='red')+
  geom_hline(yintercept = coef(model)[1],
             linetype = 'dotted')+
  geom_vline(xintercept = coef(model)[2],
             linetype = 'dotted')+
  theme_bw()

Data%>%
  filter(Scenario=='TA+, IL2-, PDL1/L2-, CD80/86-',Rate>0)%>%
  dplyr::select(Rate,Time,Scenario,
         TAex,TA,CFUNC,CINHI,TAPOP)%>%
  pivot_longer(!c(Rate,Time,Scenario))%>%
  mutate(rate = Rate,
         Rate = factor(Rate),
         name = fct_inorder(name))%>%  
  ggplot()+
  geom_line(aes(Time,value,
                group=(log10(rate)+3)/6,
                colour=(log10(rate)+3)/6),
            linewidth=0.5,
            alpha=0.5)+
  ylim(0,1)+
  scale_color_gradient(name = 'Rate',
                       low ='#cccccc',
                       high = '#111111',
                       guide ='colourbar')+
  #facet_wrap(~name,scales = 'free_x',nrow=1)+
  facet_grid(Scenario~name,scales = 'free_x')+
  theme_bw()+
  ylab('Activity level')+
  xlab('Time (a.u.)')+
  theme(legend.position='right',
        legend.title.position = 'top',
        legend.title = element_text(size = 9,
                                    hjust = 0.5,
                                    vjust = 1,
                                    face = 'bold'),
        legend.key.height = unit(1,'cm'),
        legend.key.width = unit(0.3,'cm'),
        legend.text = element_text(size = 6))
ggsave('Plot_rates.pdf',plot,width=9,height=3)
