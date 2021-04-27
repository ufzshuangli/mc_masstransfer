#load library
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

#set working directory
setwd(".../code_lis_210406")

#import data
abundance<-read.table("RA.txt",header=TRUE,na.strings = "NA",blank.lines.skip = FALSE)
Aabundance<-read.csv('Aabundance.csv')
all<-read.csv('all.csv')
CN_gr<-all[,c(1:4,95)]
CN_gr<-na.omit(CN_gr)
grAabundance<-Aabundance

R<-c('L1','L2','L3','L4','L5','R')
reactorname=c('L1'="L1",'L2'="L2",'L3'="L3",'L4'="L4",'L5'="L5",'R'="R")
P<-c('Closing-I','Recycle-10%','Recycle-50%','Recycle-80%','Closing-II')

#order subcommunities according to sum relative cell abundance
#abd:see 'code_lis.R'
abd=rbind(abundance[1,],abundance[which(abundance$Reactor=="L1"),],abundance[1,],abundance[which(abundance$Reactor=="L2"),],abundance[1,],abundance[which(abundance$Reactor=="L3"),],abundance[1,],abundance[which(abundance$Reactor=="L4"),],abundance[1,],abundance[which(abundance$Reactor=="L5"),],abundance[which(abundance$Reactor=="R"),])
abd$Reactor=c(rep("L1",77),rep("L2",77),rep("L3",77),rep("L4",77),rep("L5",77),rep("R",68))
abd=abd[,c(2:3,5:84)]
ab=abd[,c(3:82)]
absum<-apply(ab,2,sum)
abmean<-apply(ab,2,mean)
ab<-rbind(ab,absum)
ab<-ab[,order(-ab[454,])][1:453,]
abd=cbind(abd[,c(1:2)],ab)

#order subcommunities according to sum absolute cell abundance
Aabd=Aabundance[,c(2:3,5:84)]
Aab=Aabd[,c(3:82)]
Aabsum<-apply(Aab,2,sum)
Aab<-rbind(Aab,Aabsum)
Aab<-Aab[,order(-Aab[422,])][1:421,]
Aabd=cbind(Aabd[,c(1:2)],Aab)

#exclude undominant gates (never occurred at any sample with ra>0.0125)
Bll=all
for(i in 1:nrow(all)){
  for(j in 5:84){
    if(Bll[i,j]<=0.0125){
      Bll[i,j]=0
    }
  }
}#only keep the abundance of dominant subcommunities, the abundances of undominant subcommunities were set to zero

allmatrix<-Bll[1:449,5:84]#Bll: from all, ra<0.0125 is set to 0
j=0

for(i in 1:80){
  if(sum(allmatrix[,i])<=0.0125){
    j=c(j,i)
  }
}
j=c(j[-1])

ab_d<-ab[,-which(colnames(ab)%in%paste0('G',j))]
ab_d<-cbind(abd[,c(1:2)],ab_d)
ab_d1<-melt(ab_d,id=c("Reactor","Time"))

#plot mean absolute abundance of all 80 subcommunities
Aab1<-melt(Aabd,id=c("Reactor","Time"))
figs9.1_r=ggplot(Aab1,aes(x=rev(variable),y=value*100))+geom_boxplot(color='black',outlier.size=0.5)+ylab(bquote(bold(paste('cell number (mL'^-1,')'))))+
  scale_x_discrete(labels=rev(colnames(Aab)[-c(1:2)]))+scale_y_continuous(position = "right",breaks=c(0,100000000000,200000000000),expand = expand_scale(0,0))+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(face='bold',size=10),axis.text=element_text(size=6),panel.background = element_rect(colour ="black",fill = 'white'),axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+
  coord_flip(ylim = c(0,210000000000))
figs9.1_r

#1.calculate growth rate of total cell and each subcommunity
#total cell inflow of L1 from R(same for all L1-L5, cells/min)
dat<-subset(CN_gr, Reactor=='L1')
CN_sum_in<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_in=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_sum_in)){
  CN_sum_in$Time[i]=(dat$Time[i]+dat$Time[i+1])/2#time for interval=(t1+t2)/2
  CN_sum_in$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))#interval='t1 to t2'
 CN_sum_in[which(CN_sum_in$Phase%in%c(P[1],P[5])),'sum_in']<-0
    if (CN_sum_in$Phase[i]==P[2]){
      CN_sum_in[i,'sum_in']<-0.04*(subset(CN_gr,Reactor=='R'&Time==dat$Time[i],'Cell_number')+subset(CN_gr,Reactor=='R'&Time==dat$Time[i+1],'Cell_number'))/2}
    if (CN_sum_in$Phase[i]==P[3]){
      CN_sum_in[i,'sum_in']<-0.2*(subset(CN_gr,Reactor=='R'&Time==dat$Time[i],'Cell_number')+subset(CN_gr,Reactor=='R'&Time==dat$Time[i+1],'Cell_number'))/2}
    if (CN_sum_in$Phase[i]==P[4]){
      CN_sum_in[i,'sum_in']<-0.32*(subset(CN_gr,Reactor=='R'&Time==dat$Time[i],'Cell_number')+subset(CN_gr,Reactor=='R'&Time==dat$Time[i+1],'Cell_number'))/2}
}

#cells inflow per subcommunity of L1 from R (same for all L1-L5, cells/min)
dat<-subset(grAabundance, Reactor=='L1')
CN_in<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_in)){
  CN_in$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_in$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
  CN_in[which(CN_in$Phase%in%c(P[1],P[5])),j]<-0
  if (CN_in$Phase[i]==P[2]){
  CN_in[i,j]<-0.04*(subset(grAabundance,Reactor=='R'&Time==dat$Time[i],j+1)+subset(grAabundance,Reactor=='R'&Time==dat$Time[i+1],j+1))/2}
  if (CN_in$Phase[i]==P[3]){
  CN_in[i,j]<-0.2*(subset(grAabundance,Reactor=='R'&Time==dat$Time[i],j+1)+subset(grAabundance,Reactor=='R'&Time==dat$Time[i+1],j+1))/2}
  if (CN_in$Phase[i]==P[4]){
  CN_in[i,j]<-0.32*(subset(grAabundance,Reactor=='R'&Time==dat$Time[i],j+1)+subset(grAabundance,Reactor=='R'&Time==dat$Time[i+1],j+1))/2}
  }
}


#cells outflow per subcommunity of L1
dat<-subset(grAabundance, Reactor=='L1')
CN_out_L1<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_out_L1)){
  CN_out_L1$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_out_L1$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
      CN_out_L1[i,j]<-0.4*(subset(dat,Time==dat$Time[i],j+1)+subset(dat,Time==dat$Time[i+1],j+1))/2
  }
}

#cells outflow per subcommunity of L2
dat<-subset(grAabundance, Reactor=='L2')
CN_out_L2<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_out_L2)){
  CN_out_L2$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_out_L2$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_out_L2[i,j]<-0.4*(subset(dat,Time==dat$Time[i],j+1)+subset(dat,Time==dat$Time[i+1],j+1))/2
  }
}

#cells outflow per subcommunity of L3
dat<-subset(grAabundance, Reactor=='L3')
CN_out_L3<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_out_L3)){
  CN_out_L3$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_out_L3$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_out_L3[i,j]<-0.4*(subset(dat,Time==dat$Time[i],j+1)+subset(dat,Time==dat$Time[i+1],j+1))/2
  }
}

#cells outflow per subcommunity of L4
dat<-subset(grAabundance, Reactor=='L4')
CN_out_L4<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_out_L4)){
  CN_out_L4$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_out_L4$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_out_L4[i,j]<-0.4*(subset(dat,Time==dat$Time[i],j+1)+subset(dat,Time==dat$Time[i+1],j+1))/2
  }
}

#cells outflow per subcommunity of L5
dat<-subset(grAabundance, Reactor=='L5')
CN_out_L5<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_out_L5)){
  CN_out_L5$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_out_L5$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_out_L5[i,j]<-0.4*(subset(dat,Time==dat$Time[i],j+1)+subset(dat,Time==dat$Time[i+1],j+1))/2
  }
}

#cells inflow and outflow per subcommunity of R
dat<-subset(grAabundance, Reactor=='R')
CN_out_R<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
CN_in_R<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_out_R)){
  CN_out_R$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_out_R$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_in_R$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_in_R$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_out_R[i,j]<-0.4*5*(subset(dat,Time==dat$Time[i],j+1)+subset(dat,Time==dat$Time[i+1],j+1))/2
    CN_in_R[i,j]<-CN_out_L1[i+5,j]+CN_out_L2[i+5,j]+CN_out_L3[i+5,j]+CN_out_L4[i+5,j]+CN_out_L5[i+5,j]
  }
}


#changed cell number per subcommunity of L1 (cells/min)
dat<-subset(grAabundance, Reactor=='L1')
CN_d_L1<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_d_L1)){
  CN_d_L1$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_L1$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_d_L1[i,j]<-800*(dat[i+1,j+1]-dat[i,j+1])/(dat$Time[i+1]-dat$Time[i])/24/60
  }
}
#total changed cell number of L1
dat<-subset(CN_gr, Reactor=='L1')
CN_d_sum_L1<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_d=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_d_sum_L1)){
  CN_d_sum_L1$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_sum_L1$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_d_sum_L1[i,'sum_d']<-800*(dat[i+1,'Cell_number']-dat[i,'Cell_number'])/(dat$Time[i+1]-dat$Time[i])/24/60
}

#changed cell number per subcommunity of L2
dat<-subset(grAabundance, Reactor=='L2')
CN_d_L2<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_d_L2)){
  CN_d_L2$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_L2$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_d_L2[i,j]<-800*(dat[i+1,j+1]-dat[i,j+1])/(dat$Time[i+1]-dat$Time[i])/24/60
  }
}

#total changed cell number of L2
dat<-subset(CN_gr, Reactor=='L2')
CN_d_sum_L2<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_d=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_d_sum_L2)){
  CN_d_sum_L2$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_sum_L2$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_d_sum_L2[i,'sum_d']<-800*(dat[i+1,'Cell_number']-dat[i,'Cell_number'])/(dat$Time[i+1]-dat$Time[i])/24/60
}

#changed cell number per subcommunity of L3
dat<-subset(grAabundance, Reactor=='L3')
CN_d_L3<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_d_L3)){
  CN_d_L3$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_L3$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_d_L3[i,j]<-800*(dat[i+1,j+1]-dat[i,j+1])/(dat$Time[i+1]-dat$Time[i])/24/60
  }
}

#total changed cell number of L3
dat<-subset(CN_gr, Reactor=='L3')
CN_d_sum_L3<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_d=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_d_sum_L3)){
  CN_d_sum_L3$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_sum_L3$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_d_sum_L3[i,'sum_d']<-800*(dat[i+1,'Cell_number']-dat[i,'Cell_number'])/(dat$Time[i+1]-dat$Time[i])/24/60
}

#changed cell number per subcommunity of L4
dat<-subset(grAabundance, Reactor=='L4')
CN_d_L4<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_d_L4)){
  CN_d_L4$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_L4$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_d_L4[i,j]<-800*(dat[i+1,j+1]-dat[i,j+1])/(dat$Time[i+1]-dat$Time[i])/24/60
  }
}

#total changed cell number of L4
dat<-subset(CN_gr, Reactor=='L4')
CN_d_sum_L4<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_d=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_d_sum_L4)){
  CN_d_sum_L4$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_sum_L4$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_d_sum_L4[i,'sum_d']<-800*(dat[i+1,'Cell_number']-dat[i,'Cell_number'])/(dat$Time[i+1]-dat$Time[i])/24/60
}

#changed cell number per subcommunity of L5
dat<-subset(grAabundance, Reactor=='L5')
CN_d_L5<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_d_L5)){
  CN_d_L5$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_L5$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_d_L5[i,j]<-800*(dat[i+1,j+1]-dat[i,j+1])/(dat$Time[i+1]-dat$Time[i])/24/60
  }
}

#total changed cell number of L5
dat<-subset(CN_gr, Reactor=='L5')
CN_d_sum_L5<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_d=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_d_sum_L5)){
  CN_d_sum_L5$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_sum_L5$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_d_sum_L5[i,'sum_d']<-800*(dat[i+1,'Cell_number']-dat[i,'Cell_number'])/(dat$Time[i+1]-dat$Time[i])/24/60
}

#changed cell number per subcommunity of R
dat<-subset(grAabundance, Reactor=='R')
CN_d_R<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_d_R)){
  CN_d_R$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_d_R$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_d_R[i,j]<-800*(dat[i+1,j+1]-dat[i,j+1])/(dat$Time[i+1]-dat$Time[i])/24/60
  }
}

#growth rate of total cells in L1 (/d)
dat<-subset(CN_gr, Reactor=='L1')
CN_gr_sum_L1<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_gr=rep(NA,(nrow(dat)-1)),M=rep(NA,(nrow(dat)-1)),TGR=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_gr_sum_L1)){
  CN_gr_sum_L1$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_sum_L1$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_gr_sum_L1[i,'M']<-CN_sum_in[which(CN_sum_in$Time==CN_gr_sum_L1$Time[i]),'sum_in']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60#masstransfer rate
    CN_gr_sum_L1[i,'sum_gr']<-CN_d_sum_L1[i,'sum_d']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60+0.72-CN_gr_sum_L1[i,'M']#growth rate of total cells
    CN_gr_sum_L1[i,'TGR']<-CN_gr_sum_L1[i,'sum_gr']*(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/2*800#produced cells per day per reactor
}
CN_gr_sum_L1<-data.frame(Reactor='L1',CN_gr_sum_L1)

#growth rate per subcommunity in L1 (/d)
dat<-subset(grAabundance, Reactor=='L1')
CN_gr_L1<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_gr_L1)){
  CN_gr_L1$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_L1$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_gr_L1[i,j]<-(CN_d_L1[i,j]-CN_in[i,j]+CN_out_L1[i,j])/((dat[i+1,j+1]+dat[i,j+1]))/800*2
  }
}

#growth rate of total cells in L2 (/d)
dat<-subset(CN_gr, Reactor=='L2')
CN_gr_sum_L2<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_gr=rep(NA,(nrow(dat)-1)),M=rep(NA,(nrow(dat)-1)),TGR=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_gr_sum_L2)){
  CN_gr_sum_L2$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_sum_L2$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_gr_sum_L2[i,'M']<-CN_sum_in[which(CN_sum_in$Time==CN_gr_sum_L2$Time[i]),'sum_in']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60
  CN_gr_sum_L2[i,'sum_gr']<-CN_d_sum_L2[i,'sum_d']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60+0.72-CN_gr_sum_L2[i,'M']
  CN_gr_sum_L2[i,'TGR']<-CN_gr_sum_L2[i,'sum_gr']*(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/2*800
}
CN_gr_sum_L2<-data.frame(Reactor='L2',CN_gr_sum_L2)

#growth rate per subcommunity in L2 (/d)
dat<-subset(grAabundance, Reactor=='L2')
CN_gr_L2<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_gr_L2)){
  CN_gr_L2$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_L2$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_gr_L2[i,j]<-(CN_d_L2[i,j]-CN_in[i,j]+CN_out_L2[i,j])/((dat[i+1,j+1]+dat[i,j+1]))/800*2
  }
}

#growth rate of total cells in L3 (/d)
dat<-subset(CN_gr, Reactor=='L3')
CN_gr_sum_L3<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_gr=rep(NA,(nrow(dat)-1)),M=rep(NA,(nrow(dat)-1)),TGR=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_gr_sum_L3)){
  CN_gr_sum_L3$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_sum_L3$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_gr_sum_L3[i,'M']<-CN_sum_in[which(CN_sum_in$Time==CN_gr_sum_L3$Time[i]),'sum_in']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60
  CN_gr_sum_L3[i,'sum_gr']<-CN_d_sum_L3[i,'sum_d']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60+0.72-CN_gr_sum_L3[i,'M']
  CN_gr_sum_L3[i,'TGR']<-CN_gr_sum_L3[i,'sum_gr']*(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/2*800
}
CN_gr_sum_L3<-data.frame(Reactor='L3',CN_gr_sum_L3)

#growth rate per subcommunity in L3 (/d)
dat<-subset(grAabundance, Reactor=='L3')
CN_gr_L3<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_gr_L3)){
  CN_gr_L3$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_L3$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_gr_L3[i,j]<-(CN_d_L3[i,j]-CN_in[i,j]+CN_out_L3[i,j])/((dat[i+1,j+1]+dat[i,j+1]))/800*2
  }
}

#growth rate of total cells in L4 (/d)
dat<-subset(CN_gr, Reactor=='L4')
CN_gr_sum_L4<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_gr=rep(NA,(nrow(dat)-1)),M=rep(NA,(nrow(dat)-1)),TGR=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_gr_sum_L4)){
  CN_gr_sum_L4$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_sum_L4$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_gr_sum_L4[i,'M']<-CN_sum_in[which(CN_sum_in$Time==CN_gr_sum_L4$Time[i]),'sum_in']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60
  CN_gr_sum_L4[i,'sum_gr']<-CN_d_sum_L4[i,'sum_d']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60+0.72-CN_gr_sum_L4[i,'M']
  CN_gr_sum_L4[i,'TGR']<-CN_gr_sum_L4[i,'sum_gr']*(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/2*800
}
CN_gr_sum_L4<-data.frame(Reactor='L4',CN_gr_sum_L4)

#growth rate per subcommunity in L4 (/d)
dat<-subset(grAabundance, Reactor=='L4')
CN_gr_L4<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_gr_L4)){
  CN_gr_L4$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_L4$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_gr_L4[i,j]<-(CN_d_L4[i,j]-CN_in[i,j]+CN_out_L4[i,j])/((dat[i+1,j+1]+dat[i,j+1]))/800*2
  }
}

#growth rate of total cells in L5 (/d)
dat<-subset(CN_gr, Reactor=='L5')
CN_gr_sum_L5<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),Phase=dat[1:(nrow(dat)-1),4],sum_gr=rep(NA,(nrow(dat)-1)),M=rep(NA,(nrow(dat)-1)),TGR=rep(NA,(nrow(dat)-1)))
for (i in 1:nrow(CN_gr_sum_L5)){
  CN_gr_sum_L5$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_sum_L5$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  CN_gr_sum_L5[i,'M']<-CN_sum_in[which(CN_sum_in$Time==CN_gr_sum_L5$Time[i]),'sum_in']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60
  CN_gr_sum_L5[i,'sum_gr']<-CN_d_sum_L5[i,'sum_d']/(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/800*2*24*60+0.72-CN_gr_sum_L5[i,'M']
  CN_gr_sum_L5[i,'TGR']<-CN_gr_sum_L5[i,'sum_gr']*(dat[i+1,'Cell_number']+dat[i,'Cell_number'])/2*800
}
CN_gr_sum_L5<-data.frame(Reactor='L5',CN_gr_sum_L5)

#growth rate per subcommunity in L5 (/d)
dat<-subset(grAabundance, Reactor=='L5')
CN_gr_L5<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_gr_L5)){
  CN_gr_L5$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_L5$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_gr_L5[i,j]<-(CN_d_L5[i,j]-CN_in[i,j]+CN_out_L5[i,j])/((dat[i+1,j+1]+dat[i,j+1]))/800*2
  }
}

#growth rate per subcommunity in R (/d)
dat<-subset(grAabundance, Reactor=='R')
CN_gr_R<-data.frame(Time=rep(NA,(nrow(dat)-1)),Interval=rep(NA,(nrow(dat)-1)),dat[1:(nrow(dat)-1),4:84])
for (i in 1:nrow(CN_gr_R)){
  CN_gr_R$Time[i]=(dat$Time[i]+dat$Time[i+1])/2
  CN_gr_R$Interval[i]<-paste0(as.character(dat$Time[i]),"-",as.character(dat$Time[i+1]))
  for(j in 4:83){
    CN_gr_R[i,j]<-(CN_d_R[i,j]-CN_in_R[i,j]+CN_out_R[i,j])/((dat[i+1,j+1]+dat[i,j+1]))/800*2
  }
}

#2.plotting of total cell growth rate
#combine all growth rates of total cells 
CN_gr_sum_L<-rbind(CN_gr_sum_L1,CN_gr_sum_L2,CN_gr_sum_L3,CN_gr_sum_L4,CN_gr_sum_L5)
CN_gr_sum_L<-melt(CN_gr_sum_L,id=c('Reactor','Time','Interval','Phase'))

#theoratical value for masstransfer rate
ref_dm<-data.frame(x=c(0,26,26,47,47,64,64,89,89,110,0,26,26,47,47,64,64,89,89,110),y=c(0.72,0.72, 0.648,0.648,0.36,0.36,0.144,0.144,0.72,0.72,0,0,0.072,0.072,0.36,0.36,0.576,0.576,0,0),group=rep(1:10,each=2),group2=c(rep('ref.for gr',10),rep('ref. for M',10)))

#plot total grow rate and M
yname<-expression(bolditalic("µ' ")~bold(paste("(d"^-1,')'))~bold(" or ")~bolditalic("M ")~bold(paste("(d"^-1,')')))
fig3_a<-ggplot()+facet_wrap(.~Reactor, ncol = 1)+
  annotate("rect", xmin = 26, xmax = 47, ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = 47.1, xmax = 64, ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = 64.1, xmax = 89, ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_line(data=ref_dm[which(ref_dm$group2=='ref.for gr'),],aes(x=x,y=y,group=group),color='black')+
  geom_point(data=CN_gr_sum_L[which(CN_gr_sum_L$variable!='TGR'),], aes(x=Time,y=value, group=variable, color=Reactor,linetype=variable,shape=variable))+
  scale_x_continuous(breaks = c(0,26,47,64,89,110),expand = c(0,1))+
  ylab(yname)+xlab('time (d)')+
  scale_color_manual(values =c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'),name='reactors',guide = guide_legend(ncol = 5))+
  scale_shape_manual(values = c(16,17),labels=c(bquote(italic("?'")),bquote(italic("M"))))+
  guides(shape='none',color='none')+
  theme(panel.background = element_rect(color = 'black',fill='white'))+
  theme(axis.title = element_text(face="bold",size = 12),axis.title.y = element_text(size=12,margin = margin(t = 0, r = 0, b = 0, l = -10)),strip.text = element_text(face="bold",size = 12),axis.text = element_text(size = 10),panel.background = element_rect(fill='white',color = 'black'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.key.size = unit(0.4,"cm"),legend.text = element_text(size = 8),legend.position = 'bottom',legend.margin=margin(0,0,5,0),legend.box.margin =margin(-10,-10,-10,-10),legend.direction = "vertical")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())
fig3_a

#plot theoratical value for M
ann_text <- data.frame(x=c(12.5,12.5,36,36,55,55,76,76,99,99),y=c(0.85,-0.15,0.778,-0.078,0.49,0.21,0.706,-0.006,0.85,-0.15),group=c(1,2,1,2,1,2,2,1,1,2),
                       text=c("italic('D')","0","'90%'~italic('D')","'10%'~italic('D')","'50%'~italic('D')","'50%'~italic('D')","'80%'~italic('D')","'20%'~italic('D')","italic('D')","0"))
ann_text_m<-ann_text[which(ann_text$group==2),]
yname2<-expression(atop(textstyle(''),atop(textstyle(bold("hypothesized ")),paste(textstyle(bolditalic("µ' "))~textstyle(bold(paste("(d"^-1,')')))~textstyle(bold(" or "))~textstyle(bolditalic("M "))~textstyle(bold(paste("(d"^-1,')')))))))
fig3_b<-ggplot()+
  annotate("rect", xmin = 26, xmax = 47, ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = 47.1, xmax = 64, ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = 64.1, xmax = 89, ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_line(data=ref_dm[which(ref_dm$group2=='ref.for gr'),],aes(x=x,y=y,group=group),color='black')+geom_line(data=ref_dm[which(ref_dm$group2=='ref. for M'),],aes(x=x,y=y,group=group),color='black',linetype='dotted')+
  geom_point(data=CN_gr_sum_L[which(CN_gr_sum_L$variable!='TGR'),], aes(x=Time,y=value, group=variable, color=Reactor,linetype=variable,shape=variable))+#geom_line(data=CN_gr_sum_L[which(CN_gr_sum_L$variable!='TGR'),], aes(x=Time,y=value, group=variable, color=Reactor,linetype=variable,shape=variable))+
  scale_x_continuous(breaks = c(0,26,47,64,89,110),expand = c(0,1))+
  scale_shape_manual(values = c(NA,NA),labels=c(bquote(italic("?'")),bquote(italic("M"))))+
  guides(shape='none',color='none')+
  theme(panel.background = element_rect(color = 'black',fill='white'))+
  ylab(yname2)+xlab('time (d)')+
  theme(axis.title = element_text(face="bold",size = 12),axis.title.y = element_text(size=12,margin = margin(t = 0, r = 0, b = 0, l = -10)),strip.text = element_text(face="bold",size = 12),axis.text = element_text(size = 10),panel.background = element_rect(fill='white',color = 'black'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.key.size = unit(0.4,"cm"),legend.text = element_text(size = 8),legend.position = 'bottom',legend.margin=margin(0,0,5,0),legend.box.margin =margin(0,0,-5,0),legend.direction = "vertical")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  geom_text(data=ann_text_m,aes(x=x,y=y,label=text,color=as.factor(group)),parse = TRUE,inherit.aes = FALSE,size=3,color='black')
fig3_b

fig3_ab=plot_grid(fig3_a,fig3_b,rel_heights = c(4, 1),ncol = 1,align = "v")
fig3_ab

##3.growth rate per subcommunity in the whole system
#growth rate per subcommunity of whole system(L1-L5+R, /d)
dat1<-subset(grAabundance, Reactor=='L1')
dat2<-subset(grAabundance, Reactor=='L2')
dat3<-subset(grAabundance, Reactor=='L3')
dat4<-subset(grAabundance, Reactor=='L4')
dat5<-subset(grAabundance, Reactor=='L5')
dat6<-subset(grAabundance, Reactor=='R')
CN_gr_W<-data.frame(Time=rep(NA,(nrow(dat6)-1)),Interval=rep(NA,(nrow(dat6)-1)),dat6[1:(nrow(dat6)-1),4:84])
for (i in 1:nrow(CN_gr_W)){
  CN_gr_W$Time[i]=(dat6$Time[i]+dat6$Time[i+1])/2
  CN_gr_W$Interval[i]<-paste0(as.character(dat6$Time[i]),"-",as.character(dat6$Time[i+1]))
  for(j in 4:83){
    if (CN_gr_W[i,'Phase']%in%c(P[1],P[5])){
      CN_gr_W[i,j]<-(CN_d_L1[i+5,j]+CN_d_L2[i+5,j]+CN_d_L3[i+5,j]+CN_d_L4[i+5,j]+CN_d_L5[i+5,j]+CN_d_R[i,j]+CN_out_R[i,j])/((dat1[i+6,j+1]+dat1[i+5,j+1])+(dat2[i+6,j+1]+dat2[i+5,j+1])+(dat3[i+6,j+1]+dat3[i+5,j+1])+(dat4[i+6,j+1]+dat4[i+5,j+1])+(dat5[i+6,j+1]+dat5[i+5,j+1])+(dat6[i+1,j+1]+dat6[i,j+1]))/800*2
    }
    if (CN_gr_W[i,'Phase']==P[2]){
      CN_gr_W[i,j]<-(CN_d_L1[i+5,j]+CN_d_L2[i+5,j]+CN_d_L3[i+5,j]+CN_d_L4[i+5,j]+CN_d_L5[i+5,j]+CN_d_R[i,j]+0.9*CN_out_R[i,j])/((dat1[i+6,j+1]+dat1[i+5,j+1])+(dat2[i+6,j+1]+dat2[i+5,j+1])+(dat3[i+6,j+1]+dat3[i+5,j+1])+(dat4[i+6,j+1]+dat4[i+5,j+1])+(dat5[i+6,j+1]+dat5[i+5,j+1])+(dat6[i+1,j+1]+dat6[i,j+1]))/800*2
    }
    if (CN_gr_W[i,'Phase']==P[3]){
      CN_gr_W[i,j]<-(CN_d_L1[i+5,j]+CN_d_L2[i+5,j]+CN_d_L3[i+5,j]+CN_d_L4[i+5,j]+CN_d_L5[i+5,j]+CN_d_R[i,j]+0.5*CN_out_R[i,j])/((dat1[i+6,j+1]+dat1[i+5,j+1])+(dat2[i+6,j+1]+dat2[i+5,j+1])+(dat3[i+6,j+1]+dat3[i+5,j+1])+(dat4[i+6,j+1]+dat4[i+5,j+1])+(dat5[i+6,j+1]+dat5[i+5,j+1])+(dat6[i+1,j+1]+dat6[i,j+1]))/800*2
    }
    if (CN_gr_W[i,'Phase']==P[4]){
      CN_gr_W[i,j]<-(CN_d_L1[i+5,j]+CN_d_L2[i+5,j]+CN_d_L3[i+5,j]+CN_d_L4[i+5,j]+CN_d_L5[i+5,j]+CN_d_R[i,j]+0.2*CN_out_R[i,j])/((dat1[i+6,j+1]+dat1[i+5,j+1])+(dat2[i+6,j+1]+dat2[i+5,j+1])+(dat3[i+6,j+1]+dat3[i+5,j+1])+(dat4[i+6,j+1]+dat4[i+5,j+1])+(dat5[i+6,j+1]+dat5[i+5,j+1])+(dat6[i+1,j+1]+dat6[i,j+1]))/800*2
    }
  }
}
CN_gr_sc<-rbind(CN_gr_L1,CN_gr_L2,CN_gr_L3,CN_gr_L4,CN_gr_L5,CN_gr_R,CN_gr_W)
CN_gr_sc<-data.frame(Reactor=c(rep('L1',nrow(CN_gr_L1)),rep('L2',nrow(CN_gr_L2)),rep('L3',nrow(CN_gr_L3)),rep('L4',nrow(CN_gr_L4)),rep('L5',nrow(CN_gr_L5)),rep('R',nrow(CN_gr_R)),rep('W',nrow(CN_gr_W))),CN_gr_sc)

#calculate average relative cell abundance and absolute cell abundance of the system
abundance_W<-subset(abundance,Reactor=='R')
Aabundance_W<-subset(Aabundance,Reactor=='R')
Aabundance_W$Reactor<-'W'
abundance_W$Reactor<-'W'
for (k in 1:80){
  for (i in 1:nrow(Aabundance_W)){
    Aabundance_W[i,4+k]<-mean(Aabundance[which(Aabundance$Time==Aabundance_W[i,'Time']),4+k])#average of absolute cell abundance in each reactor
    abundance_W[i,4+k]<-sum(Aabundance[which(Aabundance$Time==Aabundance_W[i,'Time']),4+k])/sum(CN[which(CN$Time==Aabundance_W[i,'Time']),'Cell_number'])#sum(absolute cell abundance per subcommunity)/sum(total cell number)
  }
}
abundance1<-rbind(abundance,abundance_W)
Aabundance1<-rbind(Aabundance,Aabundance_W)

##4.average growth rate per subcommunity in L1-L5
#average the growth rates of 14day without adapation of each phase and each gate 
avg_gr_ra<-data.frame(Gate=rep(NA,80*35),Reactor=rep(c(R,'W'),5*80),Phase=rep(c(rep(P[1],7),rep(P[2],7),rep(P[3],7),rep(P[4],7),rep(P[5],7)),80),gr=rep(NA,2800),ra=rep(NA,2800),aa=rep(NA,2800))
for(k in 1:80){
  n<-35*k-34
  avg_gr_ra[c(n:(n+34)),'Gate']<-paste0('G',k)
 for(j in 1:5) {
   for(i in 1:7){
    #average growth rate
    avg_gr_ra[which(avg_gr_ra$Phase==P[j]&avg_gr_ra$Reactor==c(R,'W')[i]&avg_gr_ra$Gate==paste0('G',k)),'gr']<-mean(CN_gr_sc[which(CN_gr_sc$Reactor==c(R,'W')[i]&CN_gr_sc$Phase==P[j]&!CN_gr_sc$Time%in%c(1.5,3.5,6,26.5,27.5,28.5,29.5,31.5,47.5,48.5,49.5,50.5,52.5,64.5,65.5,66.5,67.5,68.5,69.5,70.5,89.5,90.5,91.5,92.5,94.5)),paste0('G',k)])
    #average relative cell abundance
    avg_gr_ra[which(avg_gr_ra$Phase==P[j]&avg_gr_ra$Reactor==c(R,'W')[i]&avg_gr_ra$Gate==paste0('G',k)),'ra']<-mean(abundance1[which(abundance1$Reactor==c(R,'W')[i]&abundance1$Phase==P[j]&!abundance1$Time%in%c(1:7,26:30,47:51,64:70,89:93)),paste0('G',k)])
    #average absolute cell abundance
    avg_gr_ra[which(avg_gr_ra$Phase==P[j]&avg_gr_ra$Reactor==c(R,'W')[i]&avg_gr_ra$Gate==paste0('G',k)),'aa']<-mean(Aabundance1[which(Aabundance1$Reactor==c(R,'W')[i]&Aabundance1$Phase==P[j]&!Aabundance1$Time%in%c(1:7,26:30,47:51,64:70,89:93)),paste0('G',k)])
   }
 } 
}

avg_gr_ra$Gate<-factor(avg_gr_ra$Gate,levels = paste0('G',1:80))
avg_gr_ra$Phase<-factor(avg_gr_ra$Phase,levels = P)

#set growth rate<=-1/24/60 to -1/24/60
avg_gr_ra2<-avg_gr_ra
avg_gr_ra2[which(avg_gr_ra2$gr<=-1/24/60),'gr']<--1/24/60

#map growth rate of each subcommunity in L with its relative cell abundance in R 
avg_gr_ra3<-avg_gr_ra2
avg_gr_ra3<-data.frame(subset(avg_gr_ra3,!Reactor%in%c('R','W')),ra_r=rep(avg_gr_ra3[which(avg_gr_ra3$Reactor=='R'),'ra'],each=5))
avg_gr_ra4<-avg_gr_ra3
avg_gr_ra4[,'ra_L']<-avg_gr_ra4[,'ra']

##5.calculating adaptation and balanced period growth rates per subcommunity
#growth rate calculated from beginning and ending data of the balanced period per phase
gr_period<-data.frame(Reactor=rep(reactorname,each=10),Period=rep(c('adaptation','balanced'),30),Phase=rep(rep(P,each=2),6),grAabundance[1:60,5:84])

#intercept for adaptation and balanced periods in phases
p_intercept<-c(1,7,26,33,47,54,64,71,89,96,107)

#inflow cells into reactors (cells/d)
gr_period_in_L<-data.frame(Period=rep(c('adaptation','balanced'),5),Phase=rep(P,each=2),grAabundance[1:10,5:84])
gr_period_in_R<-data.frame(Period=rep(c('adaptation','balanced'),5),Phase=rep(P,each=2),grAabundance[1:10,5:84])
gr_period_in_L[which(gr_period_in_L$Phase%in%c(P[1],P[5])),3:82]<-0

#L1-L5
for(j in 3:82){
  gr_period_in_L[which(gr_period_in_L$Period=='adaptation'&gr_period_in_L$Phase==P[2]),j]<-(grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[3]),j+2]+grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[4]),j+2])/2*0.04*60*24
  gr_period_in_L[which(gr_period_in_L$Period=='adaptation'&gr_period_in_L$Phase==P[3]),j]<-(grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[5]),j+2]+grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[6]),j+2])/2*0.2*60*24
  gr_period_in_L[which(gr_period_in_L$Period=='adaptation'&gr_period_in_L$Phase==P[4]),j]<-(grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[7]),j+2]+grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[8]),j+2])/2*0.32*60*24
  gr_period_in_L[which(gr_period_in_L$Period=='balanced'&gr_period_in_L$Phase==P[2]),j]<-(grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[4]),j+2]+grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[5]),j+2])/2*0.04*60*24
  gr_period_in_L[which(gr_period_in_L$Period=='balanced'&gr_period_in_L$Phase==P[3]),j]<-(grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[6]),j+2]+grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[7]),j+2])/2*0.2*60*24
  gr_period_in_L[which(gr_period_in_L$Period=='balanced'&gr_period_in_L$Phase==P[4]),j]<-(grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[8]),j+2]+grAabundance[which(grAabundance$Reactor=='R'&grAabundance$Time==p_intercept[9]),j+2])/2*0.32*60*24
  
  }

#R
for(i in 1:5){
for(j in 3:82){
  k=2*i
  gr_period_in_R[which(gr_period_in_R$Period=='adaptation'&gr_period_in_R$Phase==P[i]),j]<-(sum(grAabundance[which(grAabundance$Reactor!='R'&grAabundance$Time==p_intercept[k-1]),j+2])+sum(grAabundance[which(grAabundance$Reactor!='R'&grAabundance$Time==p_intercept[k]),j+2]))/2*0.4*60*24
  gr_period_in_R[which(gr_period_in_R$Period=='balanced'&gr_period_in_R$Phase==P[i]),j]<-(sum(grAabundance[which(grAabundance$Reactor!='R'&grAabundance$Time==p_intercept[k]),j+2])+sum(grAabundance[which(grAabundance$Reactor!='R'&grAabundance$Time==p_intercept[k+1]),j+2]))/2*0.4*60*24
  
  }
}
for(j in 3:82){
  gr_period_in_R[which(gr_period_in_R$Period=='adaptation'&gr_period_in_R$Phase==P[1]),j]<-NA#R absent from the adapation period of phase1
  gr_period_in_R[which(gr_period_in_R$Period=='balanced'&gr_period_in_R$Phase==P[1]),j]<-(sum(grAabundance[which(grAabundance$Reactor!='R'&grAabundance$Time==12),j+2])+sum(grAabundance[which(grAabundance$Reactor!='R'&grAabundance$Time==p_intercept[3]),j+2]))/2*0.4*60*24
  
}

#changed cells in one reactor
gr_period_d<-data.frame(Reactor=rep(reactorname,each=10),Period=rep(c('adaptation','balanced'),30),Phase=rep(rep(P,each=2),6),grAabundance[1:60,5:84])
gr_period_d[which(gr_period_d$Reactor=='R')[1],4:83]<-NA#R absent from the adapation period of phase1

#L1-L5
for(k in 1:5){
  dat<-subset(grAabundance,Reactor==reactorname[k])
    for(j in 4:83){
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[1],j]<-(dat[which(dat$Time==p_intercept[2]),j+1]-dat[which(dat$Time==p_intercept[1]),j+1])*800/(p_intercept[2]-p_intercept[1])
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[2],j]<-(dat[which(dat$Time==p_intercept[3]),j+1]-dat[which(dat$Time==p_intercept[2]),j+1])*800/(p_intercept[3]-p_intercept[2])
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[3],j]<-(dat[which(dat$Time==p_intercept[4]),j+1]-dat[which(dat$Time==p_intercept[3]),j+1])*800/(p_intercept[4]-p_intercept[3])
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[4],j]<-(dat[which(dat$Time==p_intercept[5]),j+1]-dat[which(dat$Time==p_intercept[4]),j+1])*800/(p_intercept[5]-p_intercept[4])
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[5],j]<-(dat[which(dat$Time==p_intercept[6]),j+1]-dat[which(dat$Time==p_intercept[5]),j+1])*800/(p_intercept[6]-p_intercept[5])
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[6],j]<-(dat[which(dat$Time==p_intercept[7]),j+1]-dat[which(dat$Time==p_intercept[6]),j+1])*800/(p_intercept[7]-p_intercept[6])
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[7],j]<-(dat[which(dat$Time==p_intercept[8]),j+1]-dat[which(dat$Time==p_intercept[7]),j+1])*800/(p_intercept[8]-p_intercept[7])
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[8],j]<-(dat[which(dat$Time==p_intercept[9]),j+1]-dat[which(dat$Time==p_intercept[8]),j+1])*800/(p_intercept[9]-p_intercept[8])
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[9],j]<-(dat[which(dat$Time==p_intercept[10]),j+1]-dat[which(dat$Time==p_intercept[9]),j+1])*800/(p_intercept[10]-p_intercept[9])
    gr_period_d[which(gr_period_d$Reactor==reactorname[k])[10],j]<-(dat[which(dat$Time==p_intercept[11]),j+1]-dat[which(dat$Time==p_intercept[10]),j+1])*800/(p_intercept[11]-p_intercept[10])
  }
}

#R
dat<-subset(grAabundance,Reactor=='R')
for(j in 4:83){
  gr_period_d[which(gr_period_d$Reactor=='R')[2],j]<-(dat[which(dat$Time==p_intercept[3]),j+1]-dat[which(dat$Time==12),j+1])*800/(p_intercept[3]-12)
  gr_period_d[which(gr_period_d$Reactor=='R')[3],j]<-(dat[which(dat$Time==p_intercept[4]),j+1]-dat[which(dat$Time==p_intercept[3]),j+1])*800/(p_intercept[4]-p_intercept[3])
  gr_period_d[which(gr_period_d$Reactor=='R')[4],j]<-(dat[which(dat$Time==p_intercept[5]),j+1]-dat[which(dat$Time==p_intercept[4]),j+1])*800/(p_intercept[5]-p_intercept[4])
  gr_period_d[which(gr_period_d$Reactor=='R')[5],j]<-(dat[which(dat$Time==p_intercept[6]),j+1]-dat[which(dat$Time==p_intercept[5]),j+1])*800/(p_intercept[6]-p_intercept[5])
  gr_period_d[which(gr_period_d$Reactor=='R')[6],j]<-(dat[which(dat$Time==p_intercept[7]),j+1]-dat[which(dat$Time==p_intercept[6]),j+1])*800/(p_intercept[7]-p_intercept[6])
  gr_period_d[which(gr_period_d$Reactor=='R')[7],j]<-(dat[which(dat$Time==p_intercept[8]),j+1]-dat[which(dat$Time==p_intercept[7]),j+1])*800/(p_intercept[8]-p_intercept[7])
  gr_period_d[which(gr_period_d$Reactor=='R')[8],j]<-(dat[which(dat$Time==p_intercept[9]),j+1]-dat[which(dat$Time==p_intercept[8]),j+1])*800/(p_intercept[9]-p_intercept[8])
  gr_period_d[which(gr_period_d$Reactor=='R')[9],j]<-(dat[which(dat$Time==p_intercept[10]),j+1]-dat[which(dat$Time==p_intercept[9]),j+1])*800/(p_intercept[10]-p_intercept[9])
  gr_period_d[which(gr_period_d$Reactor=='R')[10],j]<-(dat[which(dat$Time==p_intercept[11]),j+1]-dat[which(dat$Time==p_intercept[10]),j+1])*800/(p_intercept[11]-p_intercept[10])
}

#combine
gr_period_in<-data.frame(Reactor=rep(reactorname,each=10),rbind(gr_period_in_L,gr_period_in_L,gr_period_in_L,gr_period_in_L,gr_period_in_L,gr_period_in_R))

#calculte growth rate for L1-L5
for(k in 1:5){
  dat<-subset(grAabundance,Reactor==reactorname[k])
  for(j in 4:83){
    gr_period[which(gr_period$Reactor==reactorname[k])[1],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[1],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[1],j])/(dat[which(dat$Time==p_intercept[2]),j+1]+dat[which(dat$Time==p_intercept[1]),j+1])/800*2+0.72
    gr_period[which(gr_period$Reactor==reactorname[k])[2],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[2],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[2],j])/(dat[which(dat$Time==p_intercept[3]),j+1]+dat[which(dat$Time==p_intercept[2]),j+1])/800*2+0.72
    gr_period[which(gr_period$Reactor==reactorname[k])[3],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[3],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[3],j])/(dat[which(dat$Time==p_intercept[4]),j+1]+dat[which(dat$Time==p_intercept[3]),j+1])/800*2+0.72
    gr_period[which(gr_period$Reactor==reactorname[k])[4],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[4],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[4],j])/(dat[which(dat$Time==p_intercept[5]),j+1]+dat[which(dat$Time==p_intercept[4]),j+1])/800*2+0.72
    gr_period[which(gr_period$Reactor==reactorname[k])[5],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[5],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[5],j])/(dat[which(dat$Time==p_intercept[6]),j+1]+dat[which(dat$Time==p_intercept[5]),j+1])/800*2+0.72
    gr_period[which(gr_period$Reactor==reactorname[k])[6],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[6],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[6],j])/(dat[which(dat$Time==p_intercept[7]),j+1]+dat[which(dat$Time==p_intercept[6]),j+1])/800*2+0.72
    gr_period[which(gr_period$Reactor==reactorname[k])[7],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[7],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[7],j])/(dat[which(dat$Time==p_intercept[8]),j+1]+dat[which(dat$Time==p_intercept[7]),j+1])/800*2+0.72
    gr_period[which(gr_period$Reactor==reactorname[k])[8],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[8],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[8],j])/(dat[which(dat$Time==p_intercept[9]),j+1]+dat[which(dat$Time==p_intercept[8]),j+1])/800*2+0.72
    gr_period[which(gr_period$Reactor==reactorname[k])[9],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[9],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[9],j])/(dat[which(dat$Time==p_intercept[10]),j+1]+dat[which(dat$Time==p_intercept[9]),j+1])/800*2+0.72
    gr_period[which(gr_period$Reactor==reactorname[k])[10],j]<-(gr_period_d[which(gr_period$Reactor==reactorname[k])[10],j]-gr_period_in[which(gr_period$Reactor==reactorname[k])[10],j])/(dat[which(dat$Time==p_intercept[11]),j+1]+dat[which(dat$Time==p_intercept[10]),j+1])/800*2+0.72
    
      }
}

#calculte growth rate for R
gr_period[which(gr_period$Reactor=='R')[1],4:83]<-NA
dat<-subset(grAabundance,Reactor=='R')
for(j in 4:83){
  gr_period[which(gr_period$Reactor=='R')[2],j]<-(gr_period_d[which(gr_period$Reactor=='R')[2],j]-gr_period_in[which(gr_period$Reactor=='R')[2],j])/(dat[which(dat$Time==p_intercept[3]),j+1]+dat[which(dat$Time==12),j+1])/800*2+3.6
  gr_period[which(gr_period$Reactor=='R')[3],j]<-(gr_period_d[which(gr_period$Reactor=='R')[3],j]-gr_period_in[which(gr_period$Reactor=='R')[3],j])/(dat[which(dat$Time==p_intercept[4]),j+1]+dat[which(dat$Time==p_intercept[3]),j+1])/800*2+3.6
  gr_period[which(gr_period$Reactor=='R')[4],j]<-(gr_period_d[which(gr_period$Reactor=='R')[4],j]-gr_period_in[which(gr_period$Reactor=='R')[4],j])/(dat[which(dat$Time==p_intercept[5]),j+1]+dat[which(dat$Time==p_intercept[4]),j+1])/800*2+3.6
  gr_period[which(gr_period$Reactor=='R')[5],j]<-(gr_period_d[which(gr_period$Reactor=='R')[5],j]-gr_period_in[which(gr_period$Reactor=='R')[5],j])/(dat[which(dat$Time==p_intercept[6]),j+1]+dat[which(dat$Time==p_intercept[5]),j+1])/800*2+3.6
  gr_period[which(gr_period$Reactor=='R')[6],j]<-(gr_period_d[which(gr_period$Reactor=='R')[6],j]-gr_period_in[which(gr_period$Reactor=='R')[6],j])/(dat[which(dat$Time==p_intercept[7]),j+1]+dat[which(dat$Time==p_intercept[6]),j+1])/800*2+3.6
  gr_period[which(gr_period$Reactor=='R')[7],j]<-(gr_period_d[which(gr_period$Reactor=='R')[7],j]-gr_period_in[which(gr_period$Reactor=='R')[7],j])/(dat[which(dat$Time==p_intercept[8]),j+1]+dat[which(dat$Time==p_intercept[7]),j+1])/800*2+3.6
  gr_period[which(gr_period$Reactor=='R')[8],j]<-(gr_period_d[which(gr_period$Reactor=='R')[8],j]-gr_period_in[which(gr_period$Reactor=='R')[8],j])/(dat[which(dat$Time==p_intercept[9]),j+1]+dat[which(dat$Time==p_intercept[8]),j+1])/800*2+3.6
  gr_period[which(gr_period$Reactor=='R')[9],j]<-(gr_period_d[which(gr_period$Reactor=='R')[9],j]-gr_period_in[which(gr_period$Reactor=='R')[9],j])/(dat[which(dat$Time==p_intercept[10]),j+1]+dat[which(dat$Time==p_intercept[9]),j+1])/800*2+3.6
  gr_period[which(gr_period$Reactor=='R')[10],j]<-(gr_period_d[which(gr_period$Reactor=='R')[10],j]-gr_period_in[which(gr_period$Reactor=='R')[10],j])/(dat[which(dat$Time==p_intercept[11]),j+1]+dat[which(dat$Time==p_intercept[10]),j+1])/800*2+3.6
  
}

gr_period<-na.omit(gr_period)

gr_period2<-melt(gr_period,id=c('Reactor','Period','Phase'))
#order subcommunities according to sum absolute cell abundance
gr_period2$variable<-factor(gr_period2$variable,levels = rev(colnames(Aab)))
gr_period2$Phase<-factor(gr_period2$Phase,levels =P)
#exclude subcommunity never dominant in any sample
gr_period3<-gr_period2
gr_period3<-gr_period3[which(gr_period3$variable%in%colnames(ab_d)[-c(1:2)]),]
#set maximium abs(growth rate) to 1
gr_period3[which(gr_period3$value>=1),'value']<-1
gr_period3[which(gr_period3$value<=-1),'value']<--1


#set negative to zero in color, and uniform in size
nor_gr_period2<-gr_period2
nor_gr_period2[,'value2']<-nor_gr_period2[,'value']
#size-'value2'
#mark all negative value with -1
nor_gr_period2[which(nor_gr_period2$value2<=0),'value2']<--1
#set maximum value threshold to 1
nor_gr_period2[which(nor_gr_period2$value2>=1),'value2']<-1
#color-'value'
#set maximum value threshold to 1
nor_gr_period2[which(nor_gr_period2$value>=1),'value']<-1
#mark all negative value with 0
nor_gr_period2[which(nor_gr_period2$value<=0),'value']<-0

#assign the dominance of each subcommunity in each period (dominant in at least one sample during the period)
dominant<-rbind(apply(subset(Bll,Time%in%c(1:7))[,5:84],2,sum),apply(subset(Bll,Time%in%c(7:26))[,5:84],2,sum),apply(subset(Bll,Time%in%c(26:33))[,5:84],2,sum),apply(subset(Bll,Time%in%c(33:47))[,5:84],2,sum),apply(subset(Bll,Time%in%c(47:54))[,5:84],2,sum),apply(subset(Bll,Time%in%c(54:64))[,5:84],2,sum),apply(subset(Bll,Time%in%c(64:71))[,5:84],2,sum),apply(subset(Bll,Time%in%c(71:89))[,5:84],2,sum),apply(subset(Bll,Time%in%c(89:96))[,5:84],2,sum),apply(subset(Bll,Time%in%c(96:107))[,5:84],2,sum))
dominant<-data.frame(gr_period[1:10,c(2,3)],dominant)
dominant1<-melt(dominant,id=c('Phase','Period'))
#order subcommunities according to sum absolute cell abundance
dominant1$variable<-factor(dominant1$variable,levels=rev(colnames(Aab)))
dominant1$value<-1*(dominant1$value!=0)
#add dominance to growth rate data
nor_gr_period4<-nor_gr_period2

nor_gr_period4<-data.frame(nor_gr_period4,dominant=rep(NA,nrow(nor_gr_period4)))
for (i in 1:nrow(nor_gr_period4)){
  nor_gr_period4[i,'dominant']<-dominant1[which(dominant1$Phase==nor_gr_period4$Phase[i]&dominant1$Period==nor_gr_period4$Period[i]&dominant1$variable==nor_gr_period4$variable[i]),'value']
}

figs9.1_l=ggplot(nor_gr_period4[which(nor_gr_period4$Period=='balanced'),],aes(x=Reactor,y=variable, size=abs(value)*dominant, colour=value2*dominant))+geom_point(alpha=1.2)+facet_grid(.~ Phase,scales = 'free_x',labeller=as_labeller(c('Closing-I'="Insular I",'Recycle-10%'='RC 10%','Recycle-50%'='RC 50%','Recycle-80%'='RC 80%','Closing-II'='Insular II')))+ylab('subcommunities')+xlab('reactors')+
  scale_colour_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0,name=bquote(paste("netgrowth rate (d"^-1,')')),breaks=c(-1,0.25,0.5,0.75,1),labels=c(expression(""<="0"),'0.25','0.5','0.75',expression("">="1")))+scale_size(range = c(1, 2.5),name = bquote(paste("netgrowth rate (d"^-1,')')),breaks=c(0,0.25,0.5,0.75,1),labels=c(expression(""<="0"),'0.25','0.5','0.75',expression("">="1")))+
  guides(colour = guide_legend(title.position = "top",ncol=6,override.aes = list(size = c(1, 1.3,1.6, 1.9,2.2))),size = 'none')+
  theme(axis.text=element_text(size=8),strip.text = element_text(face='bold', size = 10),legend.title = element_text(size = 12),legend.direction = "horizontal", legend.text = element_text(size = 8),legend.key.size =unit(0.4,"cm") ,legend.position = 'bottom',panel.background = element_blank(),axis.title =element_text(face='bold',size=12))
figs9.1_l

plot_grid(figs9.1_l,figs9.1_r,rel_widths = c(2, 0.5),align = "h", axis = "bt")


#assign dominance together with abundance (balanced period)
gr_period_ra<-avg_gr_ra4
for(i in 1:nrow(gr_period_ra)){
  n<-which(nor_gr_period4$Period=='balanced'&nor_gr_period4$Reactor==as.character(gr_period_ra[i,'Reactor'])&nor_gr_period4$Phase==as.character(gr_period_ra[i,'Phase'])&nor_gr_period4$variable==as.character(gr_period_ra[i,'Gate']))
  gr_period_ra[i,'gr']<-nor_gr_period4[n,'value']
  gr_period_ra[i,'dominant']<-nor_gr_period4[n,'dominant']
}

gr_period_ra$ra_r<-gr_period_ra$ra_r*100
gr_period_ra<-gr_period_ra %>%
  mutate(ra_r=cut(ra_r,breaks = c(0,1.25,5,10,Inf), labels = c("<1.25%", "1.25-5%","5-10%",">10%"),
                  right = FALSE))

#keep dominant subcommunities during balanced period
nor_gr_period_balanced<-nor_gr_period4[which(nor_gr_period4$Period=='balanced'&nor_gr_period4$Reactor!='R'&nor_gr_period4$dominant==1),]
gr_period_ra$Reactor<-droplevels(gr_period_ra$Reactor)
#add abundance 
for(i in 1:nrow(nor_gr_period_balanced)){
  nor_gr_period_balanced[i,'ra']<-gr_period_ra[which(gr_period_ra$Gate==nor_gr_period_balanced[i,'variable']&gr_period_ra$Reactor==as.character(nor_gr_period_balanced[i,'Reactor'])&gr_period_ra$Phase==nor_gr_period_balanced[i,'Phase']),'ra']
  nor_gr_period_balanced[i,'ra_r']<-gr_period_ra[which(gr_period_ra$Gate==nor_gr_period_balanced[i,'variable']&gr_period_ra$Reactor==as.character(nor_gr_period_balanced[i,'Reactor'])&gr_period_ra$Phase==nor_gr_period_balanced[i,'Phase']),'ra_r']
  nor_gr_period_balanced[i,'ra_L']<-gr_period_ra[which(gr_period_ra$Gate==nor_gr_period_balanced[i,'variable']&gr_period_ra$Reactor==as.character(nor_gr_period_balanced[i,'Reactor'])&gr_period_ra$Phase==nor_gr_period_balanced[i,'Phase']),'ra_L']
}

#Figure.3d plotting relative cell abundance against growth rate
rect_dat<-data.frame(Phase=P,ymin=c(NA,-Inf,-Inf,-Inf,NA),ymax=c(NA,Inf,Inf,Inf,NA),xmin=c(NA,-Inf,-Inf,-Inf,NA),xmax=c(NA,Inf,Inf,Inf,NA),a=c(4,1,2,3,5))
rect_dat$Phase<-factor(rect_dat$Phase,levels = P)
rect_dat$a<-factor(rect_dat$a,levels = c(1,2,3,4,5,6))

ann_text <- data.frame(value=c(-0,-0.02,-0.02),ra_r=c(19,8.5,6),
                       Phase = factor(P[4],levels = P),gate=c('bold(G12)','bold(G5)','bold(G2)'))


nor_label_gr<-data.frame(Phase=P,
                         counts=c(nrow(subset(nor_gr_period_balanced,value>0&Phase==P[1])),nrow(subset(nor_gr_period_balanced,value>0&Phase==P[2])),nrow(subset(nor_gr_period_balanced,value>0&Phase==P[3])),nrow(subset(nor_gr_period_balanced,value>0&Phase==P[4])),nrow(subset(nor_gr_period_balanced,value>0&Phase==P[5]))),
                         counts_zero=c(nrow(subset(nor_gr_period_balanced,value==0&Phase==P[1])),nrow(subset(nor_gr_period_balanced,value==0&Phase==P[2])),nrow(subset(nor_gr_period_balanced,value==0&Phase==P[3])),nrow(subset(nor_gr_period_balanced,value==0&Phase==P[4])),nrow(subset(nor_gr_period_balanced,value==0&Phase==P[5]))))
nor_label_gr$counts<-as.numeric(as.character(nor_label_gr$counts))
nor_label_gr$counts_zero<-as.numeric(as.character(nor_label_gr$counts_zero))
nor_label_gr$x=c(1.2,1.2,1.2,1.2,1.2)
nor_label_gr$x2=c(0.3,0.3,0.3,0.3,0.3)
a=ggplot()+
  facet_wrap(.~ Phase,ncol=2,labeller = as_labeller(c('Closing-I'="Insular I",'Recycle-10%'='RC 10%','Recycle-50%'='RC 50%','Recycle-80%'='RC 80%','Closing-II'='Insular II')))+
  geom_rect(data = rect_dat, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax,alpha=a),fill="grey50")+guides(alpha='none')+
  geom_vline(xintercept=c(0),colour="grey50", linetype="dashed")+
  geom_point(data=nor_gr_period_balanced,aes(x=value,y=ra*100, colour=ra_r,shape=Reactor),size=2)+
  scale_color_manual(values=brewer.pal(6,"Blues")[3:6],name='relative cell\nabundance in\nregional pool',guide = guide_legend(ncol = 1))+
  scale_shape_manual(values = c('L1'=16,'L2'=16,'L3'=16,'L4'=16,'L5'=16,'R'=17,'W'=15),guide=FALSE)+
  xlab(bquote(bold(paste("netgrowth rate (d"^-1,')')))) + ylab(bquote(bold(paste("relative cell abundance in local community (%)"))))+
  theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),axis.text=element_text(size=10),strip.text = element_text(face="bold",size = 12),
        legend.text = element_text(size = 8), legend.key = element_rect(fill = NA, color=NA),legend.title = element_text(size = 12), legend.key.size =unit(0.4,"cm") ,legend.position = 'bottom',legend.margin=margin(0,0,5,0),legend.box.margin =margin(-10,-10,-10,-10),legend.direction = "vertical",panel.background = element_rect(fill='white',color = 'black'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
fig3_d=a+scale_x_continuous(limits = c(-0.1,1),breaks=c(0,0.5,1),labels = c(expression(""<="0.0"),'0.5','1.0'))+
  geom_text(data=ann_text,aes(x=value,y=ra_r,label=gate),parse = TRUE,inherit.aes = FALSE,size=3)+
  geom_text(data=nor_label_gr,aes(x = x*2/3,y=22.5,label=paste0(nor_label_gr$counts," SCs")),inherit.aes = FALSE,size = 3)+
  geom_text(data=nor_label_gr,aes(x = x2*2/3,y=21,label=paste0(nor_label_gr$counts_zero," SCs")),inherit.aes = FALSE,size = 3)+
  theme(panel.spacing = unit(0.3, "lines"))+ theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.background = element_blank())+
  theme(
    legend.position = c(0.93, .02),
    legend.justification = c("right", "bottom"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.background=element_blank()
  )
fig3_d

fig3_cd<-plot_grid(fig3_c,#see 'code_lis.R'
                   arrangeGrob(nullGrob(),fig3_d,widths=c(1,80)),rel_heights = c(1, 2),ncol = 1)+draw_plot_label(label = c( 'c', 'd'),x=0, y=c(0.97,0.67),size=12)
fig3<-plot_grid(arrangeGrob(nullGrob(),fig3_ab,nullGrob(),ncol=1,heights=c(0.93,17.8,0.21)),fig3_cd,rel_widths = c(1, 1),ncol = 2)+draw_plot_label(label = c( 'a', 'b'),x=0, y=c(0.97,0.25),size=12)


