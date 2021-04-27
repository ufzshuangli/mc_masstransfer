#1. Introduction===================================================================
#The R-script is for the data analysis and visualization of manuscript: 'Stabilizing microbial communities by looped mass transfer'
#Modified from the script published with the article: 
#'Liu et al. (2019). Neutral mechanisms and niche differentiation in steady-state insular microbial communities revealed by single cell analysis. Environmental microbiology, 21(1), 164-181.'

#----------------------------------------------------------------------------------
#Required packages
library(vegan)
library(reshape2)
library(ggplot2)
library(grid)
library(ggsignif)
library(data.table)
library(cowplot)
library(vegan)
library(Hmisc)        

#----------------------------------------------------------------------------------
#Data sets:
#RA.txt: relative cells abundance per subcommunity per sample
#Explination for each colunm
#"Sample": name of samples;
#"Reactor": reactors where samples picked from, 'L1','L2','L3','L4','L5','R';
#"Time": sampling time (day)
#"Phase": phases with recycle rate:0-10%-50%-80%:'Closing-I'(Insular I),'Recycle-10%'(RC10),'Recycle-50%'(RC50),'Recycle-80%'(RC80),'Closing-II'(Insular II)
#"G1-80": subcommunities

#----------------------------------------------------------------------------------
#set working directory
setwd(".../code_lis_210406")

#2. Overview of parameters and communities=========================================
#The R-script for Figure S3.1, Figure S5.1, Figure S7.1, Figure S7.2 & Figure S7.3
#----------------------------------------------------------------------------------
#Figure S3.1 Overview of the biotic and abiotic parameters
#read parameters of reactor L1, L2, L3, L4, L5 and R, from raw data
OD_1<-fread("parameters/OD.csv", select = 1:4)
OD_2<-fread("parameters/OD.csv", select = 5:8)
OD_3<-fread("parameters/OD.csv", select = 9:12)
OD_4<-fread("parameters/OD.csv", select = 13:16)
OD_5<-fread("parameters/OD.csv", select = 17:20)
OD_6<-fread("parameters/OD.csv", select = 21:24)
pH_1<-fread("parameters/pH.csv", select = 1:4)
pH_2<-fread("parameters/pH.csv", select = 5:8)
pH_3<-fread("parameters/pH.csv", select = 9:12)
pH_4<-fread("parameters/pH.csv", select = 13:16)
pH_5<-fread("parameters/pH.csv", select = 17:20)
pH_6<-fread("parameters/pH.csv", select = 21:24)
EC_1<-fread("parameters/EC.csv", select = 1:4)
EC_2<-fread("parameters/EC.csv", select = 5:8)
EC_3<-fread("parameters/EC.csv", select = 9:12)
EC_4<-fread("parameters/EC.csv", select = 13:16)
EC_5<-fread("parameters/EC.csv", select = 17:20)
EC_6<-fread("parameters/EC.csv", select = 21:24)
COD_t<-fread("parameters/COD.csv", select = 1:4)
COD_s<-fread("parameters/COD.csv", select = 5:8)
COD_b<-fread("parameters/COD.csv", select = 9:12)
NH4<-read.csv("parameters/NH4.csv")
PHO<-read.csv("parameters/PHO.csv")
DW<-read.csv("parameters/dry weight.csv")
CN<-read.table("parameters/live_cell_number.txt",header=TRUE,na.strings = "NA",blank.lines.skip = FALSE)

#combine all datasets
basic<-rbind(pH_1,pH_2,pH_3,pH_4,pH_5,pH_6
             ,EC_1,EC_2,EC_3,EC_4,EC_5,EC_6,NH4,PHO,
             COD_t,COD_s,COD_b,OD_1,OD_2,OD_3,OD_4,OD_5,OD_6,DW,cbind(rep('Cell_number',421),CN),use.names=FALSE)
basic<-na.omit(basic)

#plot
blank_data <- data.frame(parameter = c("CODt", "CODt", "CODs", "CODs", "CODb", "CODb","Cell_number","Cell_number","EC"), x = 0, y = c(0,3000,0,3000,0,3000,0,8000000000,5000))
ggplot(data = subset(basic,!parameter=="PHOs"), mapping= aes(x = time, y = value, color = parameter)) + 
  geom_line(size=0.5) + xlab('time') + 
  geom_point(size=0.5)+geom_blank(data = blank_data, aes(x = x, y = y))+
  facet_grid(parameter ~ reactor,scales = 'free_y',labeller = labeller(parameter=c("pH"="pH","EC"="EC","Ammonium"="NH4","PHOt"="PHOt","CODt"="CODt","CODs"="CODs","CODb"="CODb","OD"="OD","Dry_weight"="DW","Cell_number"="CN")))+
  xlab("time (d)") + ylab("biotic and abiotic parameters")+
  scale_colour_hue("parameters",breaks=c("pH","EC","Ammonium","PHOt","CODt","CODs","CODb","OD","Dry_weight","Cell_number"),labels=c("pH",bquote(paste('EC (µS cm'^-1,')')),bquote(paste("NH4 (mg N L"^-1,')')),bquote(paste("PHOt (mg P L"^-1,')')),bquote(paste("CODt (mg L"^-1,')')),bquote(paste("CODs (mg L"^-1,')')),bquote(paste("CODb (mg L"^-1,')')),"OD",bquote(paste("DW (mg mL"^-1,')')),bquote(paste("CN (mL"^-1,')'))))+
  geom_vline(xintercept=c(26,47,64,89),colour="#990000", linetype="dashed")+ theme(axis.title = element_text(face="bold",size = 12),strip.text = element_text(face="bold",size = 12),axis.text = element_text(size = 10))+
  theme(legend.key.size = unit(0.4,"cm"),legend.text = element_text(size = 8),legend.title = element_text(size = 12),legend.position = 'bottom',legend.direction = "vertical",legend.margin=margin(0,0,5,0),legend.box.margin =margin(0,0,-5,0))+
  guides(color = guide_legend(ncol = 5))

#----------------------------------------------------------------------------------
#Figure S5.1 overview of cell number
figs5.1<-ggplot(CN, aes(x=Time, y=Cell_number, color=Reactor,group=Reactor))+
  annotate("rect", xmin = 26, xmax = 47, ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = 47.1, xmax = 64, ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = 64.1, xmax = 89, ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_line()+geom_point()+
  scale_x_continuous(limit=c(0,110),breaks = c(0,26,47,64,89,110),expand = c(0,1))+
  scale_color_manual(values =c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'),name='reactors')+
  guides(color=guide_legend(ncol = 6))+
  theme(panel.background = element_rect(color = 'black',fill='white'))+
  ylab(bquote(bold(paste("cell number (mL"^-1,')'))))+xlab('time (d)')+
  theme(axis.title = element_text(face="bold",size = 12),axis.title.y = element_text(size=12),strip.text = element_text(face="bold",size = 12),axis.text = element_text(size = 10),legend.title = element_text(size = 12),panel.background = element_rect(fill='white',color = 'black'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.key.size = unit(0.4,"cm"),legend.text = element_text(size = 8),legend.position = 'bottom',legend.margin=margin(0,0,5,0),legend.box.margin =margin(0,0,-5,0),legend.direction = "vertical")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())
figs5.1

#----------------------------------------------------------------------------------
#Figure S7.1 relative cell abundance per dominant SC
#read relative cell abundance from raw data
abundance<-read.table("RA.txt",header=TRUE,na.strings = "NA",blank.lines.skip = FALSE)
ra_L1<-subset(abundance, Reactor=='L1')[,3:84]
ra_L2<-subset(abundance, Reactor=='L2')[,3:84]
ra_L3<-subset(abundance, Reactor=='L3')[,3:84]
ra_L4<-subset(abundance, Reactor=='L4')[,3:84]
ra_L5<-subset(abundance, Reactor=='L5')[,3:84]
ra_R<-subset(abundance, Reactor=='R')[,3:84]

#rank dominant subcommunities by their relative cell abundance
##L1
RA_L1<-ra_L1[,3:82]
rownames(RA_L1)<-ra_L1$Time
gates<-colnames(RA_L1)

for(i in 1:nrow(RA_L1)){
  for(j in 1:80){
    if(RA_L1[i,j]<0.0125){
      RA_L1[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(RA_L1)){
  c=sort(RA_L1[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(RA_L1)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L1melt=cbind(data2,reactor=rep("L1",80*nrow(RA_L1)))
L1melt<-na.omit(L1melt)

##L2
RA_L2<-ra_L2[,3:82]
rownames(RA_L2)<-ra_L2$Time
gates<-colnames(RA_L2)

for(i in 1:nrow(RA_L2)){
  for(j in 1:80){
    if(RA_L2[i,j]<0.0125){
      RA_L2[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(RA_L2)){
  c=sort(RA_L2[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(RA_L2)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L2melt=cbind(data2,reactor=rep("L2",80*nrow(RA_L2)))
L2melt<-na.omit(L2melt)

##L3
RA_L3<-ra_L3[,3:82]
rownames(RA_L3)<-ra_L3$Time
gates<-colnames(RA_L3)

for(i in 1:nrow(RA_L3)){
  for(j in 1:80){
    if(RA_L3[i,j]<0.0125){
      RA_L3[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(RA_L3)){
  c=sort(RA_L3[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(RA_L3)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L3melt=cbind(data2,reactor=rep("L3",80*nrow(RA_L3)))
L3melt<-na.omit(L3melt)

##L4
RA_L4<-ra_L4[,3:82]
rownames(RA_L4)<-ra_L4$Time
gates<-colnames(RA_L4)

for(i in 1:nrow(RA_L4)){
  for(j in 1:80){
    if(RA_L4[i,j]<0.0125){
      RA_L4[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(RA_L4)){
  c=sort(RA_L4[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(RA_L4)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L4melt=cbind(data2,reactor=rep("L4",80*nrow(RA_L4)))
L4melt<-na.omit(L4melt)

##L5
RA_L5<-ra_L5[,3:82]
rownames(RA_L5)<-ra_L5$Time
gates<-colnames(RA_L5)

for(i in 1:nrow(RA_L5)){
  for(j in 1:80){
    if(RA_L5[i,j]<0.0125){
      RA_L5[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(RA_L5)){
  c=sort(RA_L5[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(RA_L5)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L5melt=cbind(data2,reactor=rep("L5",80*nrow(RA_L5)))
L5melt<-na.omit(L5melt)

##R
RA_R<-ra_R[,3:82]
rownames(RA_R)<-ra_R$Time
gates<-colnames(RA_R)

for(i in 1:nrow(RA_R)){
  for(j in 1:80){
    if(RA_R[i,j]<0.0125){
      RA_R[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(RA_R)){
  c=sort(RA_R[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(RA_R)*(j-1)+i
    data2[n,"time"]<-i+8
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
Rmelt=cbind(data2,reactor=rep("R",80*nrow(RA_R)))
Rmelt<-na.omit(Rmelt)

##combined data sets and only keep dominant subcommunities for visulization
bind=rbind(L1melt,L2melt,L3melt,L4melt,L5melt,Rmelt)
meltdata<-na.omit(bind)
meltdata<-meltdata[order(meltdata$time),]
meltdata$value<-factor(meltdata$value, levels=paste0("G",1:80))

#plotting
bre=c(19,33,46,61)
colorcode=rainbow(80)
p=ggplot(data=meltdata,aes(x=time,y=abundance*100,fill=value,group=value))+#ggplot(data=meltdata,aes(x=time,y=abundance*100,fill=value,group=variable))+# if subcommunities stack in order of abundance
  annotate("rect", xmin = bre[1], xmax = bre[2], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[2]+0.1, xmax = bre[3], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[3]+0.1, xmax = bre[4], ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_bar(position = position_stack(reverse = TRUE),colour="black",size =0.01,alpha=0.7,stat="identity")+facet_wrap(.~reactor,scales="free",nrow=3,labeller = as_labeller(reactorname))+scale_fill_manual(values=colorcode[c(1:34,36:65,67:74,77,79,80)],name="subcommunities")+
  scale_x_continuous(limits = c(0,77),expand = c(0,1),breaks=c(1,bre,76),labels=c(timepoints[1],timepoints[bre],timepoints[76]),name = "time (d)")+scale_y_continuous(limits=c(0,100),expand=c(0,2),name = "relative cell abundance (%)")+
  theme(axis.title.x =element_text(face='bold',size=12), axis.title.y=element_text(face='bold',size=12),axis.text=element_text(size=10),strip.text = element_text(face='bold',size = 12),legend.title = element_text(size = 12), legend.text = element_text(size = 8),legend.key.size = unit(0.2,"cm"),panel.background = element_rect(color='black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom",legend.direction = "vertical")+ 
  guides(fill = guide_legend(override.aes = list(size = 0.01),ncol = 14))+
  theme(strip.background = element_blank())
p

#----------------------------------------------------------------------------------
#Figure S7.2 absolute cell abundance 
#calculate absolute cell abundance
Aabundance<-abundance
for (i in 1:449){
  for(j in 5:84){
    if (abundance[i,3]%in%CN$Time){
      k<-which(CN$Reactor==as.character(Aabundance[i,2]) & CN$Time==Aabundance[i,3])
      Aabundance[i,j]=abundance[i,j]*CN[k,3]  
    } else 
    {Aabundance[i,j]=NaN}
  }
}
Aabundance<-na.omit(Aabundance)
row.names(Aabundance)<-1:nrow(Aabundance)

#output absolute cell abundance for growth rate calculation
write.csv(Aabundance,'Aabundance.csv',row.names = FALSE)

aa_L1<-subset(Aabundance, Reactor=='L1')[,3:84]
aa_L2<-subset(Aabundance, Reactor=='L2')[,3:84]
aa_L3<-subset(Aabundance, Reactor=='L3')[,3:84]
aa_L4<-subset(Aabundance, Reactor=='L4')[,3:84]
aa_L5<-subset(Aabundance, Reactor=='L5')[,3:84]
aa_R<-subset(Aabundance, Reactor=='R')[,3:84]
CNL1<-subset(CN, Reactor=='L1')[,'Cell_number']
CNL2<-subset(CN, Reactor=='L2')[,'Cell_number']
CNL3<-subset(CN, Reactor=='L3')[,'Cell_number']
CNL4<-subset(CN, Reactor=='L4')[,'Cell_number']
CNL5<-subset(CN, Reactor=='L5')[,'Cell_number']
CNR<-subset(CN, Reactor=='R')[,'Cell_number']

##L1
AA_L1<-aa_L1[,3:82]
rownames(AA_L1)<-aa_L1$Time
gates<-colnames(AA_L1)

for(i in 1:nrow(AA_L1)){
  for(j in 1:80){
    if(AA_L1[i,j]<CNL1[i]/80){
      AA_L1[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(AA_L1)){
  c=sort(AA_L1[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(AA_L1)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L1melt_aa=cbind(data2,reactor=rep("L1",80*nrow(AA_L1)))
L1melt_aa<-na.omit(L1melt_aa)

##L2
AA_L2<-aa_L2[,3:82]
rownames(AA_L2)<-aa_L2$Time
gates<-colnames(AA_L2)

for(i in 1:nrow(AA_L2)){
  for(j in 1:80){
    if(AA_L2[i,j]<CNL2[i]/80){
      AA_L2[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(AA_L2)){
  c=sort(AA_L2[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(AA_L2)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L2melt_aa=cbind(data2,reactor=rep("L2",80*nrow(AA_L2)))
L2melt_aa<-na.omit(L2melt_aa)

##L3
AA_L3<-aa_L3[,3:82]
rownames(AA_L3)<-aa_L3$Time
gates<-colnames(AA_L3)

for(i in 1:nrow(AA_L3)){
  for(j in 1:80){
    if(AA_L3[i,j]<CNL3[i]/80){
      AA_L3[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(AA_L3)){
  c=sort(AA_L3[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(AA_L3)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L3melt_aa=cbind(data2,reactor=rep("L3",80*nrow(AA_L3)))
L3melt_aa<-na.omit(L3melt_aa)

##L4
AA_L4<-aa_L4[,3:82]
rownames(AA_L4)<-aa_L4$Time
gates<-colnames(AA_L4)

for(i in 1:nrow(AA_L4)){
  for(j in 1:80){
    if(AA_L4[i,j]<CNL4[i]/80){
      AA_L4[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(AA_L4)){
  c=sort(AA_L4[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(AA_L4)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L4melt_aa=cbind(data2,reactor=rep("L4",80*nrow(AA_L4)))
L4melt_aa<-na.omit(L4melt_aa)

##L5
AA_L5<-aa_L5[,3:82]
rownames(AA_L5)<-aa_L5$Time
gates<-colnames(AA_L5)

for(i in 1:nrow(AA_L5)){
  for(j in 1:80){
    if(AA_L5[i,j]<CNL5[i]/80){
      AA_L5[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(AA_L5)){
  c=sort(AA_L5[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(AA_L5)*(j-1)+i
    data2[n,"time"]<-i
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
L5melt_aa=cbind(data2,reactor=rep("L5",80*nrow(AA_L5)))
L5melt_aa<-na.omit(L5melt_aa)

##R
AA_R<-aa_R[,3:82]
rownames(AA_R)<-aa_R$Time
gates<-colnames(AA_R)

for(i in 1:nrow(AA_R)){
  for(j in 1:80){
    if(AA_R[i,j]<CNR[i]/80){
      AA_R[i,j]=0
    }
  }
}

data2<-data.frame()
for(i in 1:nrow(AA_R)){
  c=sort(AA_R[i,],decreasing = TRUE)
  for(j in 1:80){
    n<-nrow(AA_R)*(j-1)+i
    data2[n,"time"]<-i+8
    data2[n,"variable"]<-j
    data2[n,"abundance"]<-as.numeric(c[j])
    if(c[j]>0){
      c[j]=names(c[j])
      data2[n,"value"]<-c[j]
    }else{
      c[j]=NA
      data2[n,"value"]=NA
    }
  }
}
#'variable'-order of the subcommunity

##restructure the data set for plotting
Rmelt_aa=cbind(data2,reactor=rep("R",80*nrow(AA_R)))
Rmelt_aa<-na.omit(Rmelt_aa)

##combined data sets and only keep dominant subcommunities for visulization
bind=rbind(L1melt_aa,L2melt_aa,L3melt_aa,L4melt_aa,L5melt_aa,Rmelt_aa)
meltdata_aa<-na.omit(bind)
meltdata_aa<-meltdata_aa[order(meltdata_aa$time),]
meltdata_aa$value<-factor(meltdata_aa$value, levels=paste0("G",1:80))

bre=c(15,29,42,57)
p=ggplot(data=meltdata_aa,aes(x=time,y=abundance,fill=value,group=value))+#ggplot(data=meltdata_aa,aes(x=time,y=abundance,fill=value,group=variable))+ #if subcommunities stack in order of abundance
  annotate("rect", xmin = bre[1], xmax = bre[2], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[2]+0.1, xmax = bre[3], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[3]+0.1, xmax = bre[4], ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_bar(position = position_stack(reverse = TRUE),colour="black",size =0.01,alpha=0.7,stat="identity")+facet_wrap(.~reactor,scales="free",nrow=3,labeller = as_labeller(reactorname))+scale_fill_manual(values=colorcode[c(1:34,36:65,67:74,77,79,80)],name="subcommunities")+
  scale_x_continuous(limits = c(0,72),expand = c(0,1),breaks=c(1,bre,71),labels=c(timepoints[1],timepoints[bre],timepoints[71]),name = "time (d)")+scale_y_continuous(limits=c(0,5500000000),expand=c(0,2),name = bquote(bold(paste('absolute cell abundance (mL'^-1,')'))))+
  theme(axis.title.x =element_text(face='bold',size=12), axis.title.y=element_text(face='bold',size=12),axis.text=element_text(size=10),strip.text = element_text(face='bold',size = 12),legend.title = element_text(size = 12), legend.text = element_text(size = 8),legend.key.size = unit(0.2,"cm"),panel.background = element_rect(color='black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom",legend.direction = "vertical")+ 
  guides(fill = guide_legend(override.aes = list(size = 0.01),ncol = 14))+
  theme(strip.background = element_blank())
p

#----------------------------------------------------------------------------------
#Figure S7.3 fate of subcommunities (relative cell abundance)
#add preculture toeach reactor
abd=rbind(abundance[1,],abundance[which(abundance$Reactor=="L1"),],abundance[1,],abundance[which(abundance$Reactor=="L2"),],abundance[1,],abundance[which(abundance$Reactor=="L3"),],abundance[1,],abundance[which(abundance$Reactor=="L4"),],abundance[1,],abundance[which(abundance$Reactor=="L5"),],abundance[which(abundance$Reactor=="R"),])
abd$Reactor=c(rep("L1",77),rep("L2",77),rep("L3",77),rep("L4",77),rep("L5",77),rep("R",68))
abd=abd[,c(2:3,5:84)]
#order subcommunities according to sum abundance
ab=abd[,c(3:82)]
absum<-apply(ab,2,sum)
abmean<-apply(ab,2,mean)
ab<-rbind(ab,absum)
ab<-ab[,order(-ab[454,])][1:453,]
abd=cbind(abd[,c(1:2)],ab)

melt=melt(abd,id.vars = c("Reactor","Time"))

a=ggplot(data=melt[which(melt$variable%in%colnames(ab)[1:20]),],aes(x=Time,y=value*100))+facet_grid(variable ~ Reactor)+geom_line()+scale_x_continuous(name="time (d)",breaks=c(0,55,110))+scale_y_continuous(name="relative cell abundance per subcommunity (%)",breaks=c(0,50),limits=c(0,50))+theme(axis.title = element_text(face="bold",size=12),axis.text=element_text(size=8),strip.text = element_text(face='bold',size = 12),panel.background = element_blank())+geom_vline(xintercept=c(26,47,64,89),colour="grey", linetype="dashed",size=0.05)+theme(panel.spacing = unit(0.5, "lines"))
a
b=ggplot(data=melt[which(melt$variable%in%colnames(ab)[21:40]),],aes(x=Time,y=value*100))+facet_grid(variable ~ Reactor)+geom_line()+scale_x_continuous(name="time (d)",breaks=c(0,55,110))+scale_y_continuous(name="relative cell abundance per subcommunity (%)",breaks=c(0,50),limits=c(0,50))+theme(axis.title = element_text(face="bold",size=12),axis.text=element_text(size=8),strip.text = element_text(face='bold',size = 12),panel.background = element_blank())+geom_vline(xintercept=c(26,47,64,89),colour="grey", linetype="dashed",size=0.05)+theme(panel.spacing = unit(0.5, "lines"))
b
c=ggplot(data=melt[which(melt$variable%in%colnames(ab)[41:60]),],aes(x=Time,y=value*100))+facet_grid(variable ~ Reactor)+geom_line()+scale_x_continuous(name="time (d)",breaks=c(0,55,110))+scale_y_continuous(name="relative cell abundance per subcommunity (%)",breaks=c(0,50),limits=c(0,50))+theme(axis.title = element_text(face="bold",size=12),axis.text=element_text(size=8),strip.text = element_text(face='bold',size = 12),panel.background = element_blank())+geom_vline(xintercept=c(26,47,64,89),colour="grey", linetype="dashed",size=0.05)+theme(panel.spacing = unit(0.5, "lines"))
c
d=ggplot(data=melt[which(melt$variable%in%colnames(ab)[61:80]),],aes(x=Time,y=value*100))+facet_grid(variable ~ Reactor)+geom_line()+scale_x_continuous(name="time (d)",breaks=c(0,55,110))+scale_y_continuous(name="relative cell abundance per subcommunity (%)",breaks=c(0,50),limits=c(0,50))+theme(axis.title = element_text(face="bold",size=12),axis.text=element_text(size=8),strip.text = element_text(face='bold',size = 12),panel.background = element_blank())+geom_vline(xintercept=c(26,47,64,89),colour="grey", linetype="dashed",size=0.05)+theme(panel.spacing = unit(0.5, "lines"))
d

#result show in a pdf file which named "FateofSubcommunities.pdf"
pdf("FateofSubcommunities.pdf",width=6.8,height = 9.5)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(a, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(b, vp=viewport(layout.pos.row=1,layout.pos.col=2))

grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(c, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(d, vp=viewport(layout.pos.row=1,layout.pos.col=2))
dev.off()

#----------------------------------------------------------------------------------
#combine all parameters with relative cell abundance data, will be used for 5.correlation analysis
all<-abundance
all[c("OD","pH","EC","CODt","CODs","CODb","Ammonium","PHOs","PHOt","Dry_weight","Cell_number")]=NaN
for (i in 1:nrow(all)){
  if (all[i,3]%in%subset(basic,reactor==as.character(all[i,2]))$time)
  {k<-which(basic$reactor==as.character(all[i,2]) & basic$time==all[i,3])
  for (j in 1:length(k))
  {n=as.character(as.matrix(basic[k[j],1]))
  all[i,n]=basic[k[j],4] 
  } 
  }
}

write.csv(all,'all.csv',row.names = FALSE)


#3. diversity calculation (alpha, gamma, beta)=====================================
#The R-script is for diversity calculation and visualization Figure S7.5, Figure S7.4 and Figure 2
#for nestedness and turnover of beta diversity, see R-script 'code_betapart_lis.R' 

#----------------------------------------------------------------------------------
#3.1 alpha diversity
#The function SC.div was built to calculate gate based Hill numbers (Liu et al., 2019)
#In the function, requirements are q: the order(0, 1, 2) for calculation; tr: threshold for distinguish dominant subcommunities (i.e. 1/68=0.01471[1.47%] in this study); data: data set used for calculation.
SC.div<-function(q,tr,data){
  abundantdata=data.frame(data)
  nor.data=data.frame(data)
  D=matrix(nrow=nrow(data),ncol=1)
  rownames(D)=rownames(data)
  if(tr==0){
    nor.data=data
  }else
  {
    for(i in 1:ncol(data)) 
    {
      for(j in 1:nrow(data)) 
      {
        if(data[j,i]>=tr){
          abundantdata[j,i]=data[j,i] 
        }else{
          abundantdata[j,i]=0
        }
      }
    }
    
    for(j in 1:nrow(data)) 
    {
      rowsum=sum(abundantdata[j,])
      for(i in 1:ncol(data)) 
      {
        nor.data[j,i]=abundantdata[j,i]/rowsum
      }
    }
  }
  caldata=nor.data
  if (q==1){
    D[,1]=exp(diversity(caldata,index="shannon"))
  }else{
    for(i in 1:ncol(data)) 
    {
      for(j in 1:nrow(data)) 
      {
        if(nor.data[j,i]==0){
          caldata[j,i]=0
        }else{
          caldata[j,i]=nor.data[j,i]^q
        }
      }
    }
    
    for(k in 1:nrow(data))
    {
      D[k,1]=(sum(caldata[k,]))^(1/(1-q))
    }
  }
  
  return(D)  
}

#For example, to calculate the Dq=0 of dominant subcommunities in this study (relative cell abundance >1/number of total gates=0.0125)
D0dSc=SC.div(q=0,tr=0.0125,data=abundance[5:84])
rownames(D0dSc)=abundance$Sample
#Dq=0 of all subcommunities in this study
D0allSc=SC.div(q=0,tr=0,data=abundance[5:84])
#Dq=1 of all subcommunities in this study
D1allSc=SC.div(q=1,tr=0,data=abundance[5:84])
#Dq=2 of all subcommunities in this study
D2allSc=SC.div(q=2,tr=0,data=abundance[5:84])
#combine all alphadiversity
AlphaDiversity=data.frame('Sample'=abundance$Sample,'D0dSc'=D0dSc,'D0allSc'=D0allSc,'D1allSc'=D1allSc, 'D2allSc'=D2allSc)

#----------------------------------------------------------------------------------
#3.2 gamma diversity (number of subcommunities that were dominant in at least one local communities)
Bll=all
for(i in 1:nrow(all)){
  for(j in 5:84){
    if(Bll[i,j]<=0.0125){
      Bll[i,j]=0
    }
  }
}#only keep the abundance of dominant subcommunities, the abundances of undominant subcommunities were set to zero

timepoints<-subset(all,all$Reactor=='L1')$Time #extra timepoints
gamma<-data.frame(Time=timepoints,Phase=Bll$Phase[2:77],Gamma=rep(NA,length(timepoints)))
for(i in 1:nrow(gamma)){
  ga<-0
  dat<-subset(Bll,Time==gamma[i,'Time']&Reactor!='R',5:84)
  ga<-sum(1*(apply(dat,2,sum)>0.0125))
  gamma[i,'Gamma']<-ga
}#calculate gamma diversity

#----------------------------------------------------------------------------------
#Figure S7.5 gamma diversity
gamma2<-gamma
gamma2$Time<-1:76
bre=c(19,33,46,61)#breaks between phases
ggplot(gamma2,aes(x=Time,y=Gamma,group=1))+
  annotate("rect", xmin = bre[1], xmax = bre[2], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[2]+0.1, xmax = bre[3], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[3]+0.1, xmax = bre[4], ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_line(size=1)+geom_point(size=1.5,stroke =1.2)+
  scale_x_continuous(limits = c(0,77),expand = c(0,1),breaks=c(1,19,33,46,61,76),labels=c(timepoints[c(1,19,33,46,61,76)]))+scale_y_continuous(name = expression(bold(paste(gamma," diversity"))))+xlab("time (d)")+
  theme(axis.text = element_text(size=10),axis.title=element_text(face='bold',size=12),panel.background = element_rect(fill='white',color='black'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(strip.background = element_blank())

#----------------------------------------------------------------------------------
#compare gamma diversity among phases
gamma3<-gamma2
gamma3<-gamma3[which(gamma3$Time%in%c(8:18,24:32,38:45,51:60,68:76)),]#only keep data in balanced periods

#two-sided wilcox test between pairwise phases
data<-gamma3[which(gamma3$Phase%in%c(P[1],P[2])),]
wilcox.test(data[,3]~data[,2])

data<-gamma3[which(gamma3$Phase%in%c(P[1],P[3])),]
wilcox.test(data[,3]~data[,2])

data<-gamma3[which(gamma3$Phase%in%c(P[1],P[4])),]
wilcox.test(data[,3]~data[,2])

data<-gamma3[which(gamma3$Phase%in%c(P[1],P[5])),]
wilcox.test(data[,3]~data[,2])

data<-gamma3[which(gamma3$Phase%in%c(P[2],P[3])),]
wilcox.test(data[,3]~data[,2])

data<-gamma3[which(gamma3$Phase%in%c(P[2],P[4])),]
wilcox.test(data[,3]~data[,2])

data<-gamma3[which(gamma3$Phase%in%c(P[2],P[5])),]
wilcox.test(data[,3]~data[,2])

data<-gamma3[which(gamma3$Phase%in%c(P[3],P[4])),]
wilcox.test(data[,3]~data[,2])

data<-gamma3[which(gamma3$Phase%in%c(P[3],P[5])),]
wilcox.test(data[,3]~data[,2])

data<-gamma3[which(gamma3$Phase%in%c(P[4],P[5])),]
wilcox.test(data[,3]~data[,2])

#----------------------------------------------------------------------------------
#3.3 intra-community beta diversity
#L1
n<-length(ra_L1$Time[-1])
INAL1<-data.frame(Time=ra_L1$Time[-1], Values=rep(0,n), Group=rep(2,n), Reactor=rep("L1",n))
INAL1$Values[1]=0

for(i in 2:nrow(RA_L1)){
  for(j in 1:80){
    if(RA_L1[i-1,j]*RA_L1[i,j]==0&!RA_L1[i-1,j]+RA_L1[i,j]==0){
      INAL1[i-1,"Values"]=INAL1[i-1,"Values"]+1
    }
  }
}

INAL1$Time<-2:nrow(RA_L1)

#L2
n<-length(ra_L2$Time[-1])
INAL2<-data.frame(Time=ra_L2$Time[-1], Values=rep(0,n), Group=rep(2,n), Reactor=rep("L2",n))
INAL2$Values[1]=0

for(i in 2:nrow(RA_L2)){
  for(j in 1:80){
    if(RA_L2[i-1,j]*RA_L2[i,j]==0&!RA_L2[i-1,j]+RA_L2[i,j]==0){
      INAL2[i-1,"Values"]=INAL2[i-1,"Values"]+1
    }
  }
}

INAL2$Time<-2:nrow(RA_L2)

#L3
n<-length(ra_L3$Time[-1])
INAL3<-data.frame(Time=ra_L3$Time[-1], Values=rep(0,n), Group=rep(2,n), Reactor=rep("L3",n))
INAL3$Values[1]=0

for(i in 2:nrow(RA_L3)){
  for(j in 1:80){
    if(RA_L3[i-1,j]*RA_L3[i,j]==0&!RA_L3[i-1,j]+RA_L3[i,j]==0){
      INAL3[i-1,"Values"]=INAL3[i-1,"Values"]+1
    }
  }
}

INAL3$Time<-2:nrow(RA_L3)

#L4
n<-length(ra_L4$Time[-1])
INAL4<-data.frame(Time=ra_L4$Time[-1], Values=rep(0,n), Group=rep(2,n), Reactor=rep("L4",n))
INAL4$Values[1]=0

for(i in 2:nrow(RA_L4)){
  for(j in 1:80){
    if(RA_L4[i-1,j]*RA_L4[i,j]==0&!RA_L4[i-1,j]+RA_L4[i,j]==0){
      INAL4[i-1,"Values"]=INAL4[i-1,"Values"]+1
    }
  }
}

INAL4$Time<-2:nrow(RA_L4)

#L5
n<-length(ra_L5$Time[-1])
INAL5<-data.frame(Time=ra_L5$Time[-1], Values=rep(0,n), Group=rep(2,n), Reactor=rep("L5",n))
INAL5$Values[1]=0

for(i in 2:nrow(RA_L5)){
  for(j in 1:80){
    if(RA_L5[i-1,j]*RA_L5[i,j]==0&!RA_L5[i-1,j]+RA_L5[i,j]==0){
      INAL5[i-1,"Values"]=INAL5[i-1,"Values"]+1
    }
  }
}

INAL5$Time<-2:nrow(RA_L5)

#R
n<-length(ra_R$Time[-1])
INAR<-data.frame(Time=ra_R$Time[-1], Values=rep(0,n), Group=rep(2,n), Reactor=rep("R",n))
INAR$Values[1]=0

for(i in 2:nrow(RA_R)){
  for(j in 1:80){
    if(RA_R[i-1,j]*RA_R[i,j]==0&!RA_R[i-1,j]+RA_R[i,j]==0){
      INAR[i-1,"Values"]=INAR[i-1,"Values"]+1
    }
  }
}

INAR$Time<-10:(nrow(RA_R)+8)

#combine alpha diversity together with intra community beta diversity
divCs<-data.frame(Time=abundance$Time,Values=D0dSc,Group=1,Reactor=abundance$Reactor)
divCs<-divCs[-1,]
divCs$Time<-c(rep(1:76,5),9:76)
divCs<-rbind(divCs,INAL1,INAL2,INAL3,INAL4,INAL5,INAR)
divCs$Values<-as.numeric(divCs$Values)
divCs$Time<-as.numeric(divCs$Time)

#calculate reference line for intra-community beta diversity
div<-divCs
ref_L<-mean(div[which(div$Time%in%c(15:18)&div$Reactor!='R'&div$Group==2),'Values'])#reference for L using last five points
ref_L2<-mean(div[which(div$Time%in%c(8:18)&div$Reactor!='R'&div$Group==2),'Values'])#reference for L using last ten points
ref_R<-mean(div[which(div$Time%in%c(15:18)&div$Reactor=='R'&div$Group==2),'Values'])#reference for R using last five points
ref_R2<-mean(div[which(div$Time%in%c(8:18)&div$Reactor=='R'&div$Group==2),'Values'])#reference for R using last ten points
div<-rbind(div,c(0,ref_L,'ref_L','L1'),c(77,ref_L,'ref_L','L1'),c(0,ref_L2,'ref_L2','L1'),c(77,ref_L2,'ref_L2','L1'),c(0,ref_L,'ref_L','L2'),c(77,ref_L,'ref_L','L2'),c(0,ref_L2,'ref_L2','L2'),c(77,ref_L2,'ref_L2','L2'),c(0,ref_L,'ref_L','L3'),c(77,ref_L,'ref_L','L3'),c(0,ref_L2,'ref_L2','L3'),c(77,ref_L2,'ref_L2','L3'),c(0,ref_L,'ref_L','L4'),c(77,ref_L,'ref_L','L4'),c(0,ref_L2,'ref_L2','L4'),c(77,ref_L2,'ref_L2','L4'),c(0,ref_L,'ref_L','L5'),c(77,ref_L,'ref_L','L5'),c(0,ref_L2,'ref_L2','L5'),c(77,ref_L2,'ref_L2','L5'),c(0,ref_R,'ref_R','R'),c(77,ref_R,'ref_R','R'),c(0,ref_R2,'ref_R2','R'),c(77,ref_R2,'ref_R2','R'))
div$Time<-as.numeric(div$Time)
div$Values<-as.numeric(div$Values)

#----------------------------------------------------------------------------------
#Figure S7.4 alpha diversity
bre<-c(19,33,46,61,76)
p_alpha=ggplot(data=subset(div,Group==1),aes(x=Time,y=Values,color=Reactor,linetype = Group,group=Group,shape=Group))+
  annotate("rect", xmin = bre[1], xmax = bre[2], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[2]+0.1, xmax = bre[3], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[3]+0.1, xmax = bre[4], ymin = -Inf, ymax = Inf, alpha = .35)+
  facet_wrap(.~Reactor,scales="free",nrow=3,labeller = as_labeller(reactorname))+
  geom_line(size=1)+geom_point(size=1.5,stroke =1.2)+
  scale_color_manual(values =c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'))+
  scale_linetype_manual(values = c(1,2,4,4,4,4))+
  scale_shape_manual(values=c(19,21,NA,NA,NA,NA))+
  scale_x_continuous(limits = c(0,77),expand = c(0,1),breaks=c(1,19,33,46,61,76),labels=c(timepoints[c(1,19,33,46,61,76)]))+scale_y_continuous(limits=c(0,30),expand=c(0.03,0.025),breaks=c(0,5,10,15,20,25,30),labels=c(0,5,10,15,20,25,30),name = expression(bold(paste(alpha," diversity "))))+xlab("time (d)")+
  theme(axis.text = element_text(size = 10))+
  theme(panel.background  = element_rect(color = "black",fill = 'white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text = element_text(face="bold",size=12),axis.title = element_text(face="bold",size=12))+
  theme(legend.position = "none")+theme(strip.background = element_blank())
p_alpha

#----------------------------------------------------------------------------------
#3.4 inter community beta diversity
INTER<-data.frame(Pair1=rep(c('L1','L1','L1','L1','L2','L2','L2','L3','L3','L4'),nrow(subset(Bll,Reactor=='L1'))),Pair2=rep(c('L2','L3','L4','L5','L3','L4','L5','L4','L5','L5'),nrow(subset(Bll,Reactor=='L1'))), Time=rep(NA,10*nrow(subset(Bll,Reactor=='L1'))),Phase=rep(NA,10*nrow(subset(Bll,Reactor=='L1'))),Value=rep(0,10*nrow(subset(Bll,Reactor=='L1'))))
for (i in 1:nrow(subset(Bll,Reactor=='L1'))){
  INTER[c((10*i-9):(10*i)),'Time']<-subset(Bll,Reactor=='L1',Time)[i,]
  INTER[c((10*i-9):(10*i)),'Phase']<-as.character(subset(Bll,Reactor=='L1',Phase)[i,])
  for (j in (10*i-9):(10*i)){
    p1<-Bll[which(Bll$Reactor==as.character(INTER$Pair1[j])&Bll$Time==as.numeric(INTER$Time[j])),5:84]
    p2<-Bll[which(Bll$Reactor==as.character(INTER$Pair2[j])&Bll$Time==as.numeric(INTER$Time[j])),5:84]
    for(k in 1:80){
      if(p1[k]*p2[k]==0&!p1[k]+p2[k]==0){
        INTER[j,"Value"]=INTER[j,"Value"]+1
      }
    }
  }
}

INTERR<-data.frame(Pair1=rep(c('L1','L2','L3','L4','L5'),nrow(subset(Bll,Reactor=='R'))),Pair2=rep('R',5*nrow(subset(Bll,Reactor=='R'))), Time=rep(NA,5*nrow(subset(Bll,Reactor=='R'))),Phase=rep(NA,5*nrow(subset(Bll,Reactor=='R'))),Value=rep(0,5*nrow(subset(Bll,Reactor=='R'))))
for (i in 1:nrow(subset(Bll,Reactor=='R'))){
  INTERR[c((5*i-4):(5*i)),'Time']<-subset(Bll,Reactor=='R',Time)[i,]
  INTERR[c((5*i-4):(5*i)),'Phase']<-as.character(subset(Bll,Reactor=='R',Phase)[i,])
  for (j in (5*i-4):(5*i)){
    p1<-Bll[which(Bll$Reactor==as.character(INTERR$Pair1[j])&Bll$Time==as.numeric(INTERR$Time[j])),5:84]
    p2<-Bll[which(Bll$Reactor==as.character(INTERR$Pair2[j])&Bll$Time==as.numeric(INTERR$Time[j])),5:84]
    for(k in 1:80){
      if(p1[k]*p2[k]==0&!p1[k]+p2[k]==0){
        INTERR[j,"Value"]=INTERR[j,"Value"]+1
      }
    }
  }
}

INTER_c<-rbind(data.frame(INTER,Group=1),data.frame(INTERR,Group=2))
INTER_c2<-INTER_c
INTER_c2$Group<-as.factor(INTER_c2$Group)
INTER_c2$Time<-as.factor(INTER_c2$Time)

#----------------------------------------------------------------------------------
#Figure 2 intra- and inter-community beta diversity
#plot intra beta diversity, reference line using the last ten points of Insular I
p_intra<-ggplot(data=subset(div,Group%in%c(2,'ref_L2')),aes(x=Time,y=Values,color=Reactor,linetype = Group,group=Group,shape=Group))+
  annotate("rect", xmin = bre[1], xmax = bre[2], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[2]+0.1, xmax = bre[3], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[3]+0.1, xmax = bre[4], ymin = -Inf, ymax = Inf, alpha = .35)+
  facet_wrap(.~Reactor,scales="free",nrow=2,labeller = as_labeller(reactorname))+
  geom_line(size=0.5)+geom_point(size=1,stroke =1)+
  scale_color_manual(values =c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'))+
  scale_linetype_manual(values = c(1,2,4,4,4,4))+
  scale_shape_manual(values=c(21,NA,NA,NA,NA))+
  scale_x_continuous(limits = c(0,77),expand = c(0,1),breaks=c(1,19,33,46,61,76),labels=c(timepoints[c(1,19,33,46,61,76)]))+scale_y_continuous(limits=c(0,40),expand=c(0.03,0.025),breaks=c(0,10,20,30,40),labels=c(0,10,20,30,40),name = expression(bold(paste("intra-community ",beta,"-diversity "))) )+xlab("time (d)")+
  theme(axis.text = element_text(size = 10))+
  theme(panel.background  = element_rect(color = "black",fill = 'white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text = element_text(face="bold",size=12),axis.title = element_text(face="bold",size=12))+
  theme(legend.position = "none")+theme(strip.background = element_blank())
p_intra

#plot inter-community beta diversity
ann_p<-data.frame(x=c(9,26,39,53,68),y=c(35,35,35,35,35),
                  pha=c("bold('Insular I')","bold('RC 10')","bold('RC 50')","bold('RC 80')","bold('Insular II')"))
p_inter<-ggplot()+geom_boxplot(data=INTER_c2,aes(x=Time,y=Value,fill=Group),color='black',outlier.size=1,lwd=0.2)+
  annotate("rect", xmin = bre[1], xmax = bre[2], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[2]+0.1, xmax = bre[3], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[3]+0.1, xmax = bre[4], ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_text(data=ann_p,aes(x=x,y=y,label=pha),parse = TRUE,inherit.aes = FALSE,size=4)+
  xlab("time (d)")+ylab(expression(bold(paste("inter-community ",beta,"-diversity "))) )+theme(axis.title =element_text(face="bold",size=12), axis.text=element_text(size=10),panel.background =element_rect(colour ="black",fill = 'white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.text = element_text(size = 8),legend.position = 'bottom',legend.key.size = unit(0.4,"cm"),legend.direction = "vertical",legend.title = element_text(size=12),legend.margin=margin(0,0,5,0),legend.box.margin =margin(0,0,-5,0))+
  scale_fill_manual(name = 'pairs',values=c("white","black"),labels = c('L vs. L','L vs. R'))+scale_color_manual()+guides(color=FALSE,fill=guide_legend(ncol = 2))+
  scale_x_discrete(breaks=c(1,26,47,64,89,110),expand = expand_scale(add = 2))+scale_y_continuous(limits = c(0,40))
p_inter

plot_grid(p_intra,p_inter,labels=c('a','b'),ncol=1,vjust=1,rel_heights=c(1:1))

#----------------------------------------------------------------------------------
#compare alpha diversity among phases (balanced period)
div_alpha<-div[which(div$Group==1),]
div_alpha[which(div_alpha$Time%in%c(8:18)),'Phase']<-P[1]
div_alpha[which(div_alpha$Time%in%c(24:32)),'Phase']<-P[2]
div_alpha[which(div_alpha$Time%in%c(38:45)),'Phase']<-P[3]
div_alpha[which(div_alpha$Time%in%c(51:60)),'Phase']<-P[4]
div_alpha[which(div_alpha$Time%in%c(68:76)),'Phase']<-P[5]
div_alpha<-na.omit(div_alpha)

data<-div_alpha[which(div_alpha$Phase%in%c(P[1],P[2])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

data<-div_alpha[which(div_alpha$Phase%in%c(P[1],P[3])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

data<-div_alpha[which(div_alpha$Phase%in%c(P[1],P[4])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

data<-div_alpha[which(div_alpha$Phase%in%c(P[1],P[5])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

data<-div_alpha[which(div_alpha$Phase%in%c(P[2],P[3])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

data<-div_alpha[which(div_alpha$Phase%in%c(P[2],P[4])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

data<-div_alpha[which(div_alpha$Phase%in%c(P[2],P[5])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

data<-div_alpha[which(div_alpha$Phase%in%c(P[3],P[4])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

data<-div_alpha[which(div_alpha$Phase%in%c(P[3],P[5])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

data<-div_alpha[which(div_alpha$Phase%in%c(P[4],P[5])&div_alpha$Reactor!='R'),]
wilcox.test(data[,2]~data[,5])

#----------------------------------------------------------------------------------
#compare intra-community beta diversity among phases (balanced period)
INTRA_c1<-subset(div,Group==2)
INTRA_c1[which(INTRA_c1$Time%in%c(8:18)),'Phase']<-P[1]
INTRA_c1[which(INTRA_c1$Time%in%c(24:32)),'Phase']<-P[2]
INTRA_c1[which(INTRA_c1$Time%in%c(38:45)),'Phase']<-P[3]
INTRA_c1[which(INTRA_c1$Time%in%c(51:60)),'Phase']<-P[4]
INTRA_c1[which(INTRA_c1$Time%in%c(68:76)),'Phase']<-P[5]
INTRA_c2<-na.omit(INTRA_c1)

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[1],P[2])),]
wilcox.test(data[,2]~data[,5])

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[1],P[3])),]
wilcox.test(data[,2]~data[,5])

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[1],P[4])),]
wilcox.test(data[,2]~data[,5])

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[1],P[5])),]
wilcox.test(data[,2]~data[,5])

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[2],P[3])),]
wilcox.test(data[,2]~data[,5])

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[2],P[4])),]
wilcox.test(data[,2]~data[,5])

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[2],P[5])),]
wilcox.test(data[,2]~data[,5])

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[3],P[4])),]
wilcox.test(data[,2]~data[,5])

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[3],P[5])),]
wilcox.test(data[,2]~data[,5])

data<-INTRA_c2[which(INTRA_c2$Phase%in%c(P[4],P[5])),]
wilcox.test(data[,2]~data[,5])

#----------------------------------------------------------------------------------
#compare inter_community beta diversity among phases(balanced period)
INTER2<-INTER[which(!INTER$Time%in%c(0:7,26:30,47:51,64:70,89:93)),]
INTERR2<-INTERR[which(!INTERR$Time%in%c(0:7,26:30,47:51,64:70,89:93)),]
INTER_c3<-rbind(INTER2,INTERR2)
wilcox.test(INTER2[,5],INTERR2[,5]) 

#L vs. L
data<-INTER2[which(INTER2$Phase%in%c(P[1],P[2])),]
wilcox.test(data[,5]~data[,4])

data<-INTER2[which(INTER2$Phase%in%c(P[1],P[3])),]
wilcox.test(data[,5]~data[,4])

data<-INTER2[which(INTER2$Phase%in%c(P[1],P[4])),]
wilcox.test(data[,5]~data[,4])

data<-INTER2[which(INTER2$Phase%in%c(P[1],P[5])),]
wilcox.test(data[,5]~data[,4])

data<-INTER2[which(INTER2$Phase%in%c(P[2],P[3])),]
wilcox.test(data[,5]~data[,4])

data<-INTER2[which(INTER2$Phase%in%c(P[2],P[4])),]
wilcox.test(data[,5]~data[,4])

data<-INTER2[which(INTER2$Phase%in%c(P[2],P[5])),]
wilcox.test(data[,5]~data[,4])

data<-INTER2[which(INTER2$Phase%in%c(P[3],P[4])),]
wilcox.test(data[,5]~data[,4])

data<-INTER2[which(INTER2$Phase%in%c(P[3],P[5])),]
wilcox.test(data[,5]~data[,4])

data<-INTER2[which(INTER2$Phase%in%c(P[4],P[5])),]
wilcox.test(data[,5]~data[,4])
mean(INTER2[which(INTER2$Phase==P[1]),5])
mean(INTER2[which(INTER2$Phase==P[2]),5])
mean(INTER2[which(INTER2$Phase==P[3]),5])
mean(INTER2[which(INTER2$Phase==P[4]),5])
mean(INTER2[which(INTER2$Phase==P[5]),5])

sd(INTER2[which(INTER2$Phase==P[1]),5])
sd(INTER2[which(INTER2$Phase==P[2]),5])
sd(INTER2[which(INTER2$Phase==P[3]),5])
sd(INTER2[which(INTER2$Phase==P[4]),5])
sd(INTER2[which(INTER2$Phase==P[5]),5]) 

#L vs. R
data<-INTERR2[which(INTERR2$Phase%in%c(P[1],P[2])),]
wilcox.test(data[,5]~data[,4])

data<-INTERR2[which(INTERR2$Phase%in%c(P[1],P[3])),]
wilcox.test(data[,5]~data[,4])

data<-INTERR2[which(INTERR2$Phase%in%c(P[1],P[4])),]
wilcox.test(data[,5]~data[,4])

data<-INTERR2[which(INTERR2$Phase%in%c(P[1],P[5])),]
wilcox.test(data[,5]~data[,4])

data<-INTERR2[which(INTERR2$Phase%in%c(P[2],P[3])),]
wilcox.test(data[,5]~data[,4])

data<-INTERR2[which(INTERR2$Phase%in%c(P[2],P[4])),]
wilcox.test(data[,5]~data[,4])

data<-INTERR2[which(INTERR2$Phase%in%c(P[2],P[5])),]
wilcox.test(data[,5]~data[,4])

data<-INTERR2[which(INTERR2$Phase%in%c(P[3],P[4])),]
wilcox.test(data[,5]~data[,4])

data<-INTERR2[which(INTERR2$Phase%in%c(P[3],P[5])),]
wilcox.test(data[,5]~data[,4])

data<-INTERR2[which(INTERR2$Phase%in%c(P[4],P[5])),]
wilcox.test(data[,5]~data[,4])
#----------------------------------------------------------------------------------

#3.5 identify drifts according to intra-community beta diversity, >ref_L2 (8.76)
drift<-div[which(div$Group==2&div$Reactor!='R'&div$Values>8.76),]
drift<-drift[,c(1,2,4)]

#Table S7.2
sum_drift<-dcast(drift,Time~Reactor,value.var='Values')
sum_drift$Time<-timepoints[sum_drift$Time]
sum_drift<-sum_drift[which(!sum_drift$Time%in%c(0:7,27:33,48:54,65:71,90:96)),]
write.csv(sum_drift,'sum_drift.csv')

#4. netgrowth rate===================================================================
#For calculation of netgrowth rate, see 'code_growth rate_lis.R'
#The R script is for combined plotting of #Figure.3
#----------------------------------------------------------------------------------
#Figure.3
CN2<-CN
for(i in 1:nrow(CN2)){
  CN2[i,'Time']<-which(timepoints==CN2[i,'Time'])
}#rename timepoints according to continious time order
k<-setdiff(1:76,CN2[,'Time'])
CN_na<-data.frame(Reactor=rep('L1',5),Time=k, Cell_number=rep(NA,5))
CN2<-rbind(CN2,CN_na)
for(i in 1:nrow(CN2)){
  CN2[i,'Time']<-timepoints[CN2[i,'Time']]
}
CN2$Time<-factor(CN2$Time,levels=0:110)
fig3_c<-ggplot()+
  annotate("rect", xmin = bre[1], xmax = bre[2], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[2]+0.1, xmax = bre[3], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[3]+0.1, xmax = bre[4], ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_boxplot(data=CN2[which(CN2$Reactor!='R'),], aes(x=Time,y=Cell_number))+
  geom_point(data=CN2[which(CN2$Reactor!='R'),], aes(x=Time,y=Cell_number,color=Reactor))+
  scale_x_discrete(breaks=c(0,26,47,64,89,110),expand = c(0,1))+
  scale_y_continuous(name=bquote(bold(paste("cell number (mL"^-1,')'))))+xlab("time (d)")+
  scale_colour_manual(values=c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'),name='reactors',guide = guide_legend(ncol = 5))+
  theme(axis.text = element_text(size=10),axis.text.y =element_text(angle=90,hjust = 0.5),axis.title=element_text(face='bold',size=12),panel.background = element_rect(fill='white',color='black'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 8),legend.key.size =unit(0.4,"cm") ,legend.position = 'top',legend.margin=margin(0,0,5,0),legend.box.margin =margin(-5,-10,-10,-10),legend.direction = "vertical")+
  theme(legend.justification = c(0.85,1))
fig3_cd<-plot_grid(fig3_c,arrangeGrob(nullGrob(),fig3_d#see 'code_growth rate_lis.R'
                                      ,widths=c(1,80)),rel_heights = c(1, 2),ncol = 1)+draw_plot_label(label = c( 'c', 'd'),x=0, y=c(0.97,0.67),size=12)
fig3<-plot_grid(arrangeGrob(nullGrob(),fig3_ab#see 'code_growth rate_lis.R'
                            ,nullGrob(),ncol=1,heights=c(0.93,17.8,0.21)),fig3_cd,rel_widths = c(1, 1),ncol = 2)+draw_plot_label(label = c( 'a', 'b'),x=0, y=c(0.97,0.25),size=12)

#5.correlation analysis===================================================================
#The R-script is for correlation analysis and visualization Figure S10.1 & Figure S10.2 , except for: 
#analysis of subcommunity vs.subcommunity and subcommunity vs. parameter correlations (Table 2), for this see script 'code_correlation_lis.R'

#----------------------------------------------------------------------------------
#correlation of the relative cell abundance of each subcommunity with recycle rate (0-10%-50%-80%-0)
#exclude adaptation phase from relative cell abundance dataset
cell<-subset(abundance,!Time%in%c(0:7,26:30,47:51,64:70,89:93)&Reactor!='R')[,5:84]
Phase<-subset(abundance,!Time%in%c(0:7,26:30,47:51,64:70,89:93)&Reactor!='R')$Phase

#assign recycle rate (%) to each phase
Recycle<-as.character(Phase)
Recycle[which(Phase==P[1])]<-0
Recycle[which(Phase==P[2])]<-10
Recycle[which(Phase==P[3])]<-50
Recycle[which(Phase==P[4])]<-80
Recycle[which(Phase==P[5])]<-0
Recycle<-as.numeric(Recycle)

gate_phase_cor<-data.frame(r=rep(0,80),p=rep(0,80))
gate_phase_cor2<-data.frame(r=rep(0,80),p=rep(0,80))

#Spearson's test
for(i in 1:ncol(cell)){
  gcor<-rcorr(cell[,i],Recycle,type=c("spearman"))
  gate_phase_cor$r[i]<-gcor$r[1,2]
  gate_phase_cor$p[i]<-gcor$P[1,2]
  rownames(gate_phase_cor)[i]<-paste0('G',i)
}

#Pearson's test
for(i in 1:ncol(cell)){
  gcor<-rcorr(cell[,i],Recycle,type=c("pearson"))
  gate_phase_cor2$r[i]<-gcor$r[1,2]
  gate_phase_cor2$p[i]<-gcor$P[1,2]
  rownames(gate_phase_cor2)[i]<-paste0('G',i)
}

#keep significant correlations (p<0.05)
gate_phase_cor<-gate_phase_cor[which(gate_phase_cor$p<0.05),]
gate_phase_cor<-gate_phase_cor[order(-abs(gate_phase_cor$r)),]

gate_phase_cor2<-gate_phase_cor2[which(gate_phase_cor2$p<0.05),]
gate_phase_cor2<-gate_phase_cor2[order(-abs(gate_phase_cor2$r)),]

#----------------------------------------------------------------------------------
##correlation of the absolute cell abundance of each subcommunity with recycle rate (0-10%-50%-80%-0)
##exclude adaptation phase from absolute cell abundance dataset
Acell<-subset(Aabundance,!Time%in%c(0:7,26:30,47:51,64:70,89:93)&Reactor!='R')[,5:84]
APhase<-subset(Aabundance,!Time%in%c(0:7,26:30,47:51,64:70,89:93)&Reactor!='R')$Phase
ARecycle<-as.character(APhase)
ARecycle[which(APhase==P[1])]<-0
ARecycle[which(APhase==P[2])]<-10
ARecycle[which(APhase==P[3])]<-50
ARecycle[which(APhase==P[4])]<-80
ARecycle[which(APhase==P[5])]<-0
ARecycle<-as.numeric(ARecycle)
Agate_phase_cor<-data.frame(r=rep(0,80),p=rep(0,80))
Agate_phase_cor2<-data.frame(r=rep(0,80),p=rep(0,80))
Acell_sum<-apply(Acell,2,mean)
sum(Acell_sum[which(names(Acell_sum)%in%rownames(Agate_phase_cor2)[1:14])])/sum(Acell_sum)

#Spearson's test
for(i in 1:ncol(Acell)){
  gcor<-rcorr(Acell[,i],ARecycle,type=c("spearman"))
  Agate_phase_cor$r[i]<-gcor$r[1,2]
  Agate_phase_cor$p[i]<-gcor$P[1,2]
  rownames(Agate_phase_cor)[i]<-paste0('G',i)
}

#Pearson's test
for(i in 1:ncol(Acell)){
  gcor<-rcorr(Acell[,i],ARecycle,type=c("pearson"))
  Agate_phase_cor2$r[i]<-gcor$r[1,2]
  Agate_phase_cor2$p[i]<-gcor$P[1,2]
  rownames(Agate_phase_cor2)[i]<-paste0('G',i)
}

#keep significant correlations (p<0.05)
Agate_phase_cor<-Agate_phase_cor[which(Agate_phase_cor$p<0.05),]
Agate_phase_cor<-Agate_phase_cor[order(-abs(Agate_phase_cor$r)),]

Agate_phase_cor2<-Agate_phase_cor2[which(Agate_phase_cor2$p<0.05),]
Agate_phase_cor2<-Agate_phase_cor2[order(-abs(Agate_phase_cor2$r)),]

#----------------------------------------------------------------------------------
#Figure S10.1 correlation between absolute cell abundance of subcommunities and recycle rates
##combine all highly correlated subcommunities (Pearson's r>0.5)
dat<-data.frame(cell=Acell[,rownames(Agate_phase_cor2)[1]],Recycle=ARecycle,Gate=rep(rownames(Agate_phase_cor2)[1],nrow(Acell)))

##calculate linear coefficients and make labels
dat.lm <- lm(cell ~ Recycle, data = subset(dat,Gate==rownames(Agate_phase_cor2)[1]))
labels<-data.frame(Gate=rep(NA,14),x=rep(35,14),y1=rep(NA,14),y2=rep(NA,14),formula=rep(NA,14),r=rep(NA,14))
labels[1,'Gate']<-rownames(Agate_phase_cor2)[1]
labels[1,'formula'] <- sprintf("italic(y) == %.2e %+.2e * italic(x)",
                               round(coef(dat.lm)[1],digit=3),round(coef(dat.lm)[2]*100,digit=3))
labels[1,'r'] <- sprintf("italic(r) == %.2f",sqrt(summary(dat.lm)$r.squared))
labels[1,'y1']<-max(subset(dat,Gate==rownames(Agate_phase_cor2)[1],cell))*0.8
labels[1,'y2']<-max(subset(dat,Gate==rownames(Agate_phase_cor2)[1],cell))*0.64
for(i in 2:14){
  dat<-rbind(dat,data.frame(cell=Acell[,rownames(Agate_phase_cor2)[i]],Recycle=ARecycle,Gate=rep(rownames(Agate_phase_cor2)[i],nrow(Acell))))
  dat.lm <- lm(cell ~ Recycle, data = subset(dat,Gate==rownames(Agate_phase_cor2)[i]))
  labels[i,'Gate']<-rownames(Agate_phase_cor2)[i]
  labels[i,'formula'] <- sprintf("italic(y) == %.2e %+.2e * italic(x)",
                                 round(coef(dat.lm)[1],digit=3),round(coef(dat.lm)[2]*100,digit=3))
  labels[i,'r'] <- sprintf("italic(r) == %.2f",sqrt(summary(dat.lm)$r.squared))
  labels[i,'y1']<-max(subset(dat,Gate==rownames(Agate_phase_cor2)[i],cell))*0.8
  labels[i,'y2']<-max(subset(dat,Gate==rownames(Agate_phase_cor2)[i],cell))*0.64
}

figs10.1 <- ggplot(dat,aes(y=cell,x=Recycle,color=Gate)) +facet_wrap(Gate~.,ncol=4,scales="free")+
  geom_point(position = position_jitter(2),shape=19) +scale_color_manual(values=colorcode[c(12,5,14,36,80,33,75,13,16,23,20,3,31,37)],name='subcommunities')+
  geom_text(data=labels,mapping=aes(x = x,y=y2,label=r),parse = TRUE,inherit.aes = FALSE,size = 4) +
  xlab("recycle rate (%)") + ylab(bquote(bold(paste("absolute cell abaundance per SC (mL"^-1,')'))))+
  theme(axis.title =element_text(face='bold',size=12),axis.text.x=element_text(size=10),axis.text.y=element_text(size=6),strip.text = element_text(face='bold',size = 12),legend.title = element_text(size = 12), 
        legend.text = element_text(size = 8),legend.key.size = unit(0.4,"cm"),legend.position = "bottom",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.spacing = unit(0.2, "lines"))+
  guides(color=guide_legend(ncol = 7))+
  scale_x_continuous(breaks=c(0,10,50,80), labels = c(0,10,50,80))+geom_smooth(method = lm,colour='black')
figs10.1

#----------------------------------------------------------------------------------
#Figure S10.2 correlation between relative cell abundance in L and R
raLR<-subset(abundance, Time%in%subset(abundance,Reactor=='R')$Time&Reactor!='R')
raLR<-melt(raLR,id=c('Sample','Reactor','Time','Phase'))
raR<-melt(subset(abundance,Reactor=='R'),id=c('Sample','Reactor','Time','Phase'))
raLR<-data.frame(raLR,ra_r=NA)
for (i in 1:nrow(raLR)){
  raLR[i,'ra_r']<-raR[which(raR$Time==raLR[i,'Time']&raR$variable==raLR[i,'variable']),'value']
}

#calcualte Pearson's correlation coefficiency per phase
corLR_p1<-rcorr(subset(raLR,Phase==P[1])[,6],subset(raLR,Phase==P[1])[,7],type=c("pearson"))
corLR_p2<-rcorr(subset(raLR,Phase==P[2])[,6],subset(raLR,Phase==P[2])[,7],type=c("pearson"))
corLR_p3<-rcorr(subset(raLR,Phase==P[3])[,6],subset(raLR,Phase==P[3])[,7],type=c("pearson"))
corLR_p4<-rcorr(subset(raLR,Phase==P[4])[,6],subset(raLR,Phase==P[4])[,7],type=c("pearson"))
corLR_p5<-rcorr(subset(raLR,Phase==P[5])[,6],subset(raLR,Phase==P[5])[,7],type=c("pearson"))
labels<-data.frame(Phase=P,r=NA)
labels[1,'r']<-sprintf(" %.2f",corLR_p1$r[1,2])
labels[2,'r']<-sprintf(" %.2f",corLR_p2$r[1,2])
labels[3,'r']<-sprintf(" %.2f",corLR_p3$r[1,2])
labels[4,'r']<-sprintf(" %.2f",corLR_p4$r[1,2])
labels[5,'r']<-sprintf(" %.2f",corLR_p5$r[1,2])
raLR$Phase<-factor(raLR$Phase,levels = P)

rect_dat<-data.frame(Phase=P,ymin=c(NA,-Inf,-Inf,-Inf,NA),ymax=c(NA,Inf,Inf,Inf,NA),xmin=c(NA,-Inf,-Inf,-Inf,NA),xmax=c(NA,Inf,Inf,Inf,NA),a=c(4,1,2,3,5))
rect_dat$Phase<-factor(rect_dat$Phase,levels = P)
rect_dat$a<-factor(rect_dat$a,levels = c(1,2,3,4,5,6))

figs10.2=ggplot()+
  facet_wrap(.~Phase,ncol=3,labeller = as_labeller(c('Closing-I'='Insular I','Recycle-10%'='RC 10%','Recycle-50%'='RC 50%','Recycle-80%'='RC 80%','Closing-II'='Insular II')))+
  geom_rect(data = rect_dat, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax,alpha=a),fill="grey50")+guides(alpha='none')+
  geom_point(data=subset(raLR,Phase%in%P[1:5]),aes(x=ra_r*100,y=value*100,color=variable,group=1))+
  geom_smooth(data=subset(raLR,Phase%in%P[1:5]),aes(x=ra_r*100,y=value*100,color=variable,group=1),method = 'lm',show.legend = FALSE,inherit.aes = TRUE,color='black')+scale_color_manual(values=colorcode,name='subcommunities')+
  xlab("relative cell abundance in R (%)") + ylab("relative cell abundance in L1-L5(%)")+
  theme(axis.title.x =element_text(face='bold',size=12), axis.title.y=element_text(face='bold',size=12),axis.text=element_text(size=10),strip.text = element_text(face='bold',size = 12),legend.title = element_text(size = 12), 
        legend.text = element_text(size = 8),legend.key.size = unit(0.2,"cm"),legend.position = "bottom",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(shape=FALSE,color = guide_legend(override.aes = list(size = 1),ncol = 15))+
  geom_text(data=labels,aes(x = 20,y=37,label=paste0("Pearson's r = ",labels$r)),inherit.aes = FALSE,size = 4)+
  theme(strip.background=element_blank())
figs10.2

#6.stability test=================================================================== 
#The R-script is for visualization of stability properties, for calculation of stability properties, see R-script 'code_stability_lis.R'
#as a result of 'code_stability_lis.R', a stabilitytbl and a stabilitydata were got

#----------------------------------------------------------------------------------
#read the result 
stabilitytbl<-read.csv('stabilitytbl.csv')
stabilitydata<-read.csv('stabilitydata.csv')
#in this table
#Disturbance:
#1-Insular I to RC10
#2-RC10 to RC50
#3-RC50 to RC80
#4-RC80 to Insular II
#RS: resistance
#DS: displacement speed
#RL: resilience
#E: elasticity

#----------------------------------------------------------------------------------
#in addition, calculate r_Canberra of balanced period as constancy
#compare r_Canberra of different Phase (balanced period)
data<-data.frame(abundance,referencePhase = FALSE)
rownames(data)=data[,1]
R_Canberra<-data.frame()#blank table for saving r_Canberra
R_Canberra2<-data.frame()#r_Canberra for plotting

#calculate r_Canberra per reactor per phase
R<-c('L1','L2','L3','L4','L5','R')
P<-c('Closing-I','Recycle-10%','Recycle-50%','Recycle-80%','Closing-II')
numberOfGates<-80
for (i in 1:6){
  data<-data.frame(abundance,referencePhase = FALSE)
  rownames(data)=data[,1]
  data=subset(data,Reactor==R[i])
  for (j in 1:5){
    data$referencePhase = FALSE
    data[which(data$Phase==P[j]&!data$Time%in%c(1:7,26:30,47:51,64:70,89:93)),'referencePhase']=TRUE
    referenceState <- colMeans(data[data$referencePhase==TRUE, 5:(numberOfGates + 4)])
    data$Canberra <- apply(data[5:(numberOfGates+4)], 1, function(x) dist(rbind(referenceState, x), method = "canberra"))
    data$Canberra <- data$Canberra / numberOfGates
    r_Canberra <- max(data$Canberra[data$referencePhase==TRUE])
    R_Canberra[i,P[j]]<-r_Canberra
    R_Canberra2[(i-1)*6+j,'Reactor']<-R[i]
    R_Canberra2[(i-1)*6+j,'Phase']<-P[j]
    R_Canberra2[(i-1)*6+j,'value']<-r_Canberra
  }
}
R_Canberra<-cbind(Reactor=R,R_Canberra)
R_Canberra2<-na.omit(R_Canberra2)

#----------------------------------------------------------------------------------
#two-sided wilcox test of r_Canberra (constancy) between successive phases
data<-R_Canberra2[which(R_Canberra2$Reactor!='R'&R_Canberra2$Phase%in%c(P[1],P[2])),]
wilcox.test(data[,3]~data[,2])

data<-R_Canberra2[which(R_Canberra2$Reactor!='R'&R_Canberra2$Phase%in%c(P[2],P[3])),]
wilcox.test(data[,3]~data[,2])

data<-R_Canberra2[which(R_Canberra2$Reactor!='R'&R_Canberra2$Phase%in%c(P[3],P[4])),]
wilcox.test(data[,3]~data[,2])

data<-R_Canberra2[which(R_Canberra2$Reactor!='R'&R_Canberra2$Phase%in%c(P[4],P[5])),]
wilcox.test(data[,3]~data[,2])

#two-sided wilcox test of resistance between successive phases
data<-stabilitytbl[which(stabilitytbl$Reactor!='R'&stabilitytbl$Disturbance%in%c(1,2)),]
wilcox.test(data[,3]~data[,2])

data<-stabilitytbl[which(stabilitytbl$Reactor!='R'&stabilitytbl$Disturbance%in%c(2,3)),]
wilcox.test(data[,3]~data[,2])

data<-stabilitytbl[which(stabilitytbl$Reactor!='R'&stabilitytbl$Disturbance%in%c(3,4)),]
wilcox.test(data[,3]~data[,2])

#----------------------------------------------------------------------------------
#Figure.4 plot constancy together with resistance
#boxplot of r_Canberra (constancy), and pairwise significance comparision
compaired <- list(c("Recycle-10%", "Closing-I"),c("Recycle-50%","Recycle-10%"),c("Recycle-50%","Recycle-80%"),c("Recycle-80%","Closing-II"))
fig4_a<-ggplot(data = R_Canberra2, aes(x = factor(Phase, levels=c('Closing-I','Recycle-10%','Recycle-50%','Recycle-80%','Closing-II')), y = value, color=Reactor)) +
  geom_boxplot(color='black')+geom_point(position = position_jitter(0.1),size=3)+xlab("phases")+ylab("constancy space")+
  scale_x_discrete(labels=c('Insular I','RC10','RC50','RC80','Insular II'))+scale_y_continuous(limits = c(0.2,0.62))+
  theme(axis.text.x =element_text(angle=-30),axis.title =element_text(face='bold',size=12),axis.text=element_text(size=10),strip.text = element_text(face='bold',size = 12),legend.title = element_text(size = 12), 
        legend.text = element_text(size = 8),legend.key.size = unit(0.4,"cm"),legend.position = "none",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_colour_manual('reactors',values=c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'))+
  guides(color=guide_legend(ncol = 6))+
  geom_signif(data = R_Canberra2[which(R_Canberra2$Reactor!='R'),], comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = wilcox.test,color='black',textsize = 5)

#boxplot of resistance in stabilitytbl, and pairwise significance comparision
compaired <- list(c("1", "2"),c("2","3"),c("3","4"))
fig4_b<-ggplot(data = stabilitytbl, aes(x = factor(Disturbance, levels=c('1','2','3','4')), y = RS, color=Reactor)) +
  geom_boxplot(color='black')+geom_point(position = position_jitter(0.1),size=3)+xlab("disturbances")+ylab("resistance")+
  scale_x_discrete(labels=c('Insular I to RC10','RC10 to RC50','RC50 to RC80','RC80 to Insular II'))+scale_y_continuous(limits = c(0.35,0.77))+
  theme(axis.text.x =element_text(angle=-30),axis.title =element_text(face='bold',size=12),axis.text=element_text(size=10),strip.text = element_text(face='bold',size = 12),legend.title = element_text(size = 12), 
        legend.text = element_text(size = 8),legend.key.size = unit(0.4,"cm"),legend.position = "bottom",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_colour_manual('reactors',values=c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'))+
  guides(color=guide_legend(ncol = 6))+
  geom_signif(data=stabilitytbl[which(stabilitytbl$Reactor!='R'),],comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = wilcox.test,color='black',textsize = 5)

grobs <- ggplotGrob(p2)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

fig4<-plot_grid(fig4_a,fig_b+theme(legend.position="none"),labels=c('a','b'),vjust=1,ncol=2,align = "h")#rel_heights=c(1,1.25),
plot_grid(fig4, legend, ncol = 1, rel_heights = c(1, .15))


#7.dynamic summary=================================================================
#The R-script is for an overview of the community dynamic Figure.1

#----------------------------------------------------------------------------------
#Figure.1 NMDs and deviation
#NMDs based on relative cell abundance
BC <- metaMDS(abundance[,5:84], distance="bray", autotransform=FALSE, zerodist="add",try=100,trymax = 200)

#for deviation, calculate reference line as mean value of constancy spaces of L1-L5 per phase
R_Canberra3<-subset(R_Canberra2,Reactor=='R')

for(i in 1:5){
  R_Canberra3[5+i,'Reactor']<-'L'
  R_Canberra3[5+i,'Phase']<-P[i]
  R_Canberra3[5+i,'value']<-mean(as.numeric(as.matrix(subset(R_Canberra2,Reactor!='R'&Phase==P[i],value))))}

#deviation based on Canberra distance extracted from stabilitydata 
data_canberra<-subset(stabilitydata,Phase==P[1])[,c(1,3,85)]

#layout
colplate_reactor=c("1"='#d53e4f',"2"='#fc8d59',"3"='#fee08b',"4"='#e6f598',"5"='#99d594',"6"='black',"7"='#3288bd')
col_reactor=colplate_reactor[as.factor(abundance$Reactor)]
reactorname=c('L1'="L1",'L2'="L2",'L3'="L3",'L4'="L4",'L5'="L5",'R'="R")
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,11,11,11,11),3,5,byrow=TRUE),widths=c(1,1,1,1,1),height=c(2,2,1.2))
par(mai=c(0.3,0.2,0.25,0),omi=c(0,0.3,0.2,0.1),xpd=FALSE)
par(ps = 10, cex.main=1.8, cex.axis = 1,cex.lab=1.25)

#plot a.NMDs
#phase1-Insular I
plot(BC, type="n",xlab='',ylab='')
title("Insular I")
mtext("a", side = 3, line = -1,adj=-0.03,outer = TRUE,cex=1.25)

#other points in other phases-grey point
s=c(1:449)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#points in present phase-colored according to reactor
s=2:19#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+abundance[s,]$Time/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=78:95#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+abundance[s,]$Time/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=154:171#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+abundance[s,]$Time/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=230:247#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+abundance[s,]$Time/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=306:323#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+abundance[s,]$Time/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=382:391#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+abundance[s,]$Time/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

points(BC$points[1,1],BC$points[1,2],col="black",type = "p",pch=19)
arrows(BC$points[1,1],BC$points[1,2],BC$points[2,1],BC$points[2,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")
arrows(BC$points[1,1],BC$points[1,2],BC$points[78,1],BC$points[78,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")
arrows(BC$points[1,1],BC$points[1,2],BC$points[154,1],BC$points[154,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")
arrows(BC$points[1,1],BC$points[1,2],BC$points[230,1],BC$points[230,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")
arrows(BC$points[1,1],BC$points[1,2],BC$points[306,1],BC$points[306,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")

#phase2-RC10
plot(BC, type="n",xlab='',ylab='',yaxt="n")
title("RC 10")

#other points in other phases-grey points
s=c(1:449)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#points in present phase-colored according to reactor
s=20:33#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=96:109#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=172:185#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=248:261#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=324:337#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=392:405#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

#phase3-RC50
plot(BC, type="n",xlab='',ylab='',yaxt="n")
title("RC 50")

#other points in other phases-grey points
s=c(1:449)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#points in present phase-colored according to reactor
s=34:46#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=110:122#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=186:198#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=262:274#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=338:350#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=406:418#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

#phase4-RC80
plot(BC, type="n",xlab='',ylab='',yaxt="n")
title("RC 80")

#other points in other phases-grey points
s=c(1:449)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#points in present phase-colored according to reactor
s=47:61#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+a(abundance[s,]$Time-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=123:137#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=199:213#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=275:289#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=351:365#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=419:433#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

#phase5-Insular II
plot(BC, type="n",xlab='',ylab='',yaxt="n")
title("Insular II")

#other points in other phases-grey points
s=c(1:449)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#points in present phase-colored according to reactor
s=62:77#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=138:153#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=214:229#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=290:305#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=366:381#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=434:449#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(abundance[s,]$Time-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

#plot deviation
mtext("b", side = 3, line = -14,adj=-0.03,outer = TRUE,cex=1.25)
#canberra-p1
#deviation based on Canberra distance extracted from stabilitydata-p1 
data_canberra<-subset(stabilitydata,Phase==P[1])[,c(1,3,85)]
plot(data_canberra[,2:3], type="n",ylim=c(-0.025,1.025),xaxt="n",xlab="",ylab="")
axis(1,at=c(0,9,23),labels=c("0",'9',"23"))
abline(h=R_Canberra3[which(R_Canberra3$Reactor=='L'&R_Canberra3$Phase==P[1]),'value'],col="black",lty=2,lwd=1)
for (j in 1:6){
  s<-subset(data_canberra, Reactor==reactorname[j])
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],pch=19,cex=1)
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],type="l",lty=1)
  points(x=s[1,2],y=s[1,3],col='purple',pch=17,cex=1.5)
}

#canberra-p2
#deviation based on Canberra distance extracted from stabilitydata-p2
data_canberra<-subset(stabilitydata,Phase==P[2])[,c(1,3,85)]
plot(data_canberra[,2:3], type="n",ylim=c(-0.025,1.025),xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,at=c(26,44),labels=c("26",'44'))
abline(h=R_Canberra3[which(R_Canberra3$Reactor=='L'&R_Canberra3$Phase==P[2]),'value'],col="black",lty=2,lwd=1)
for (j in 1:6){
  s<-subset(data_canberra, Reactor==reactorname[j])
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],pch=19,cex=1)
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],type="l",lty=1)
  points(x=s[1,2],y=s[1,3],col='purple',pch=17,cex=1.5)
}

#canberra-p3
#deviation based on Canberra distance extracted from stabilitydata-p3
data_canberra<-subset(stabilitydata,Phase==P[3])[,c(1,3,85)]
plot(data_canberra[,2:3], type="n",ylim=c(-0.025,1.025),xaxt="n",yaxt="n",xlab="time (d)",ylab="")
axis(1,at=c(47,63),labels=c("47",'63'))
abline(h=R_Canberra3[which(R_Canberra3$Reactor=='L'&R_Canberra3$Phase==P[3]),'value'],col="black",lty=2,lwd=1)
for (j in 1:6){
  s<-subset(data_canberra, Reactor==reactorname[j])
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],pch=19,cex=1)
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],type="l",lty=1)
  points(x=s[1,2],y=s[1,3],col='purple',pch=17,cex=1.5)
}

#canberra-p4
#deviation based on Canberra distance extracted from stabilitydata-p4
data_canberra<-subset(stabilitydata,Phase==P[4])[,c(1,3,85)]
plot(data_canberra[,2:3], type="n",ylim=c(-0.025,1.025),xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,at=c(64,86),labels=c("64",'86'))
abline(h=R_Canberra3[which(R_Canberra3$Reactor=='L'&R_Canberra3$Phase==P[4]),'value'],col="black",lty=2,lwd=1)
for (j in 1:6){
  s<-subset(data_canberra, Reactor==reactorname[j])
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],pch=19,cex=1)
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],type="l",lty=1)
  points(x=s[1,2],y=s[1,3],col='purple',pch=17,cex=1.5)
}

#canberra-p5
#deviation based on Canberra distance extracted from stabilitydata-p5
data_canberra<-subset(stabilitydata,Phase==P[5])[,c(1,3,85)]
plot(data_canberra[,2:3], type="n",ylim=c(-0.025,1.025),xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,at=c(89,110),labels=c("89",'110'))
abline(h=R_Canberra3[which(R_Canberra3$Reactor=='L'&R_Canberra3$Phase==P[5]),'value'],col="black",lty=2,lwd=1)
for (j in 1:6){
  s<-subset(data_canberra, Reactor==reactorname[j])
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],pch=19,cex=1)
  points(x=s[,2],y=s[,3],col=colplate_reactor[c(1:5,7)][j],type="l",lty=1)
  points(x=s[1,2],y=s[1,3],col='purple',pch=17,cex=1.5)
}  

#legend
plot.new()
legend(title="",x="center",y="top",ncol=4,legend=c("L1","L2","L3","L4","L5","R","preculture"),col=c("1"='#d53e4f',"2"='#fc8d59',"3"='#fee08b',"4"='#e6f598',"5"='#99d594',"6"='#3288bd',"7"='black'),pch=c(19,19,19,19,19,19,19),lwd=1.5,title.adj = 0,box.col = "white",cex=1.25)
mtext("reactors", side = 3, line = -28,adj=0.25,outer = TRUE,cex=1.25)
mtext("NMDS2", side = 2, line = 1,adj=0.87,outer = TRUE,cex=1.25)
mtext("deviation", side = 2, line = 1,adj=0.40,outer = TRUE,cex=1.25)
mtext("time (d)", side = 3, line = -26.5,adj=0.53,outer = TRUE,cex=1.25)
mtext("NMDS1", side = 3, line = -14,adj=0.53,outer = TRUE,cex=1.25)

dev.off()
