#load library
library(betapart)#for beta diversity partitioning (betapart)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

#set working directory
setwd(".../code_lis_210406")

#import data
all<-read.csv('all.csv')

R<-c('L1','L2','L3','L4','L5','R')
reactorname=c('L1'="L1",'L2'="L2",'L3'="L3",'L4'="L4",'L5'="L5",'R'="R")
P<-c('Closing-I','Recycle-10%','Recycle-50%','Recycle-80%','Closing-II')

#exclude undominant gates (never occurred at any sample with ra>0.0125)
Bll=all
for(i in 1:nrow(all)){
  for(j in 5:84){
    if(Bll[i,j]<=0.0125){
      Bll[i,j]=0
    }
  }
}#only keep the abundance of dominant subcommunities, the abundances of undominant subcommunities were set to zero

#transform abundance to 0/1 (undominance/dominance)
pre_Bll<-cbind(Bll[,1:4],1*(Bll[,5:84]>0.0125))[-1,]
pre_Bll$Phase<-factor(pre_Bll$Phase)
pre_Bll$Reactor<-factor(pre_Bll$Reactor)

#betapart based on all samples during balanced period per reactor per phase 
betapart_Bll_phase<-data.frame(Reactor=rep(reactorname,5),Phase=rep(P,each=6),SOR=rep(NA,30),SNE=rep(NA,30),SIM=rep(NA,30))
for (i in 1:nrow(betapart_Bll_phase)){
  dat<-subset(pre_Bll,Phase==betapart_Bll_phase[i,'Phase']&Reactor==betapart_Bll_phase[i,'Reactor']&!Time%in%c(0:7,26:30,47:51,64:70,89:93),5:84)
  bll_core<-betapart.core(dat)
  bll_multi<-beta.multi(bll_core,index.family="sor")
  betapart_Bll_phase[i,'SOR']<-bll_multi$beta.SOR
  betapart_Bll_phase[i,'SNE']<-bll_multi$beta.SNE
  betapart_Bll_phase[i,'SIM']<-bll_multi$beta.SIM
}

#betapart based on all samples from L1-L5 per day
betapart_Bll_reactor<-data.frame(subset(pre_Bll,Reactor=='L1')[,c(3,4)],SOR=rep(NA,76),SNE=rep(NA,76),SIM=rep(NA,76))
for (i in 1:nrow(betapart_Bll_reactor)){
  dat<-subset(pre_Bll,Time==betapart_Bll_reactor[i,'Time']&Reactor!='R',5:84)
  bll_core<-betapart.core(dat)
  bll_multi<-beta.multi(bll_core,index.family="sor")
  betapart_Bll_reactor[i,'SOR']<-bll_multi$beta.SOR
  betapart_Bll_reactor[i,'SNE']<-bll_multi$beta.SNE
  betapart_Bll_reactor[i,'SIM']<-bll_multi$beta.SIM
}
melt_betapart_phase<-melt(betapart_Bll_phase,id=c('Reactor','Phase'))
melt_betapart_phase$Phase<-factor(melt_betapart_phase$Phase,levels = P)
melt_betapart_phase$variable<-factor(melt_betapart_phase$variable,levels = c("SOR","SIM","SNE"))

figs7.6_a=ggplot(subset(melt_betapart_phase,variable!='SOR'),aes(x=Phase,y=value,fill=variable,group=variable,color=variable))+geom_bar(position = position_stack(reverse = TRUE),colour=NA,stat="identity")+
  facet_wrap(.~Reactor,ncol=2)+scale_x_discrete(breaks=P,labels = c('Insular I','RC10','RC50','RC80','Insular II'))+
  theme(axis.text.x = element_text(size = 8,angle=-30),axis.text.y = element_text(size = 10),axis.title = element_text(face='bold',size = 12),panel.background = element_rect(fill=NA,color='black'),strip.text =element_text(face='bold',size = 12),panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
  theme(legend.key.size = unit(0.4,"cm"),legend.text = element_text(size = 8),legend.title = element_text(size = 12),legend.position = 'bottom',legend.direction = "vertical",legend.margin=margin(0,0,5,0),legend.box.margin =margin(0,0,-5,0))+
  scale_fill_manual(values=brewer.pal(2,"Set2"),name=expression(paste(beta,"-diversity partitioning")),labels=c('turnover','nestedness'))+ylab(expression(bold(paste("multisite intra-community ",beta,"-diversity "))))+xlab('phases')+
  guides(fill = guide_legend(ncol = 2))+theme(strip.background = element_blank())

melt_betapart_reactor<-melt(betapart_Bll_reactor,id=c('Time','Phase'))
melt_betapart_reactor$Time<-factor(melt_betapart_reactor$Time)
melt_betapart_reactor$variable<-factor(melt_betapart_reactor$variable,levels = c("SOR","SIM","SNE"))

figs7.6_b=ggplot(subset(melt_betapart_reactor,variable!='SOR'),aes(x=Time,y=value,fill=variable,group=variable,color=variable,width=.65))+
  annotate("rect", xmin = 19, xmax = 33, ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = 33.1, xmax = 46, ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = 46.1, xmax = 61, ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_bar(position = position_stack(reverse = TRUE),colour=NA,stat="identity")+
  scale_x_discrete(breaks=c(1,26,47,64,89,110),expand = expand_scale(add = 2))+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title = element_text(face='bold',size = 12),panel.background = element_rect(fill=NA,color='black'),strip.text =element_text(face='bold',size = 12),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.key.size = unit(0.4,"cm"),legend.text = element_text(size = 8),legend.title = element_text(size = 12),legend.position = 'bottom',legend.direction = "vertical",legend.margin=margin(0,0,5,0),legend.box.margin =margin(0,0,-5,0))+
  scale_fill_manual(values=brewer.pal(2,"Set2"),name=expression(paste(beta,"-diversity partitioning")),labels=c('turnover','nestedness'))+ylab(expression(bold(paste("multisite inter-community ",beta,"-diversity "))))+xlab('time (d)')+
  guides(fill = guide_legend(ncol = 2))

figs7.6=plot_grid(figs7.6_a,figs7.6_b,labels=c('a','b'),ncol=2)

#comparing SIM and SNE between phase; melt_betapart_phase
data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[1],P[2])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[1],P[3])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[1],P[4])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[1],P[5])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[2],P[3])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[2],P[4])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[2],P[5])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[3],P[4])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[3],P[5])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_phase[which(melt_betapart_phase$Phase%in%c(P[4],P[5])&melt_betapart_phase$variable=='SNE'&melt_betapart_phase$Reactor!='R'),]
wilcox.test(data[,4]~data[,2])


#comparing SIM and SNE between phase; melt_betapart_reactor
data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[1],P[2])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[1],P[3])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[1],P[4])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[1],P[5])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[2],P[3])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[2],P[4])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[2],P[5])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[3],P[4])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[3],P[5])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])

data<-melt_betapart_reactor[which(melt_betapart_reactor$Phase%in%c(P[4],P[5])&melt_betapart_reactor$variable=='SIM'&!melt_betapart_reactor$Time%in%c(1:7,26:30,47:51,64:70,89:93)),]
wilcox.test(data[,4]~data[,2])
