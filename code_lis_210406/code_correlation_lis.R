##load library
library(vegan)
library(dplyr)
library(ggplot2)
library(ggcor)
library(data.table)
library(Hmisc)        
library(psych)  

#set working directory
setwd(".../code_lis_210406")

#input dataset
all<-read.csv('all.csv')

##function calculating subcommunity vs.subcommunity (bb) and subcommunity vs. abiotics (ab) correlations (Spearman's test)
#s:which rows of data are used
#parameters: which columns of data are used
#min: threshold of correlation's coefficieny of high correlation, subcommunity vs.subcommunity with r>min are highly correlated. 
COR=function(s,parameters,min){
  allmatrix<-do.call(cbind,all[s,parameters])
cell<-allmatrix[,1:80]
abio<-allmatrix[,81:ncol(allmatrix)]
corr <- rcorr(cell,type = "spearman")
dfbb <- as_cor_tbl(corr)
corr2<-rcorr(cell,abio,type ="spearman")
dfab <- as_cor_tbl(corr2)
dfab<-filter(dfab,.col.names%in%colnames(abio))
dfab<-filter(dfab,.row.names%in%colnames(cell))
dfbb<-filter(dfbb,.row.id+.col.id>ncol(cell)+1)
dfab$absr<-abs(dfab$r)
dfab$p.value<-sapply(as.numeric(dfab$p.value), p.adjust, method="BH")
dfbb$p.value<-sapply(as.numeric(dfbb$p.value), p.adjust, method="BH")
bb<-nrow(subset(dfbb,dfbb$p.value<0.05&abs(dfbb$r)>min))
dfab <- dfab %>% 
  mutate(r = cut(r, breaks = c(-Inf,0,Inf), 
                 labels = c("<0", ">0"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf,0.001,0.05, Inf),
                       labels = c("<0.001", "0.001-0.05",">=0.05"),
                       right = FALSE),
         absr = cut(absr, breaks = c(0,0.5,0.75,Inf), 
                    labels = c("0.00-0.5","0.5-0.75",">=0.75"),
                    right = FALSE))
dfab<-dfab[which((dfab$p.value == '<0.001') | (dfab$p.value == '0.001-0.05')&dfab$absr==">=0.75"), ]#keep significant ab correlation with coefficiency>=0.75
ab<-nrow(subset(dfab,dfab$absr==">=0.75"))
colnames(dfab)[1:2]<-c("env","spec")
extra.params <- extra_params(spec.label = text_params(size = 4,fontface = 2,hjust = 0,vjust = 1.2),link.params =link_params( env.point.hjust = -1.8,env.point.vjust = -0.5))
#plot correlation, bb correlations as heatmap, ab correlations as network 
p=quickcor(dfbb, type = "upper",cor.test = TRUE,title='SCs vs. SCs correlation r') + geom_color() +geom_cross(color='grey50')+
  add_link(dfab, mapping = aes(colour = r),diag.label = TRUE,extra.params = extra.params,size=2) +
  add_diag_label(size=1.2,hjust=1,vjust=0.4)+theme(plot.margin = margin(-1,-5,-1,-1))+
  remove_axis("x")+scale_x_continuous(limit=c(-45,80),expand = c(0,0))+guides(colour=guide_legend(title="r"))+
  theme(axis.text.y = element_text(size = 3),legend.position = 'top',legend.text = element_text(size = 8),legend.title = element_text(size = 12))

p
#output number of bb and ab correlations
out=c('bb'=bb,'ab'=ab)
print(p)
return(out)
}

#function output all correlations with r and p
#s:which rows of data are used
#parameters: which columns of data are used
#min: threshold of correlation's coefficieny of high correlation
cp_list=function(s,parameters,min){
   allmatrix<-do.call(cbind,all[s,parameters])
  n<-colnames(allmatrix)
  cor<-rcorr(x=allmatrix,type=c("spearman"))
  r<-cor$r
  p<-c(cor$P)
  colnames(r)<-n
  rownames(r)<-n
  bh<-sapply(as.numeric(p), p.adjust, method="BH")
  padj<-c(bh)
  padj[is.na(padj)]<-0
  
  cor1=rcorr(x=allmatrix,type=c("spearman"))
  r1<-cor1$r
  p1<-c(cor1$P)
  colnames(r1)<-n
  rownames(r1)<-n
  bh1<-sapply(as.numeric(p1), p.adjust, method="BH")
  padj1<-c(bh1)
  padj1[is.na(padj1)]<-0
  
  rr1=r
  for(i in 1:length(n)){
    for(j in 1:length(n)){
      if(is.na(r[i,j])|is.na(r1[i,j])){
        rr1[i,j]=0
      }else{
        if(i>j){
          if(r[i,j]>=min&&r1[i,j]>=min){
            rr1[i,j]=min
          }else{
            if(r[i,j]<=-min&&r1[i,j]<=-min){
              rr1[i,j]=-min
            }else{
              rr1[i,j]=0
            }
          }
        }else{rr1[i,j]=0}
      }
    }
  }
  for (i in 1:length(padj)) {
    if (padj[i]>0.05){rr1[[i]]<-0}
  }
  
  for (i in 1:length(padj1)) {
    if (padj1[i]>0.05){rr1[[i]]<-0}
  }
  BB=0
  correlation_list<-data.frame("node1"=NA,"node2"=NA,"group"=NA,"r"=NA,"p"=NA)
  for(i in 1:80){
    for(j in 1:80){
      if(i!=j){
        if(rr1[i,j]>=min|rr1[i,j]<=-min){
          BB=BB+1
          add<-c(n[i],n[j],"BB",r[i,j],padj[(i-1)*length(n)+j])
          correlation_list<-rbind(correlation_list,add)
          
        }
      }
    }
  }
  
  
  AB=0
  for(i in 81:length(n)){
    for(j in 1:80){
      if(rr1[i,j]>=min|rr1[i,j]<=-min)
      {
        AB=AB+1
        add<-c(n[i],n[j],"AB",r[i,j],padj[(i-1)*length(n)+j])
        correlation_list<-rbind(correlation_list,add)
      }
    }
  }
  AA=0
  for(i in 81:length(n)){
    for(j in 81:length(n)){
      if(rr1[i,j]>=min|rr1[i,j]<=-min)
      {
        AA=AA+1
        add<-c(n[i],n[j],"AA",r[i,j],padj[(i-1)*length(n)+j])
        correlation_list<-rbind(correlation_list,add)
      }
    }
  }
  correlation_list<-correlation_list[-1,]
  return(correlation_list)
}


#per phase per reactor
parameters<-c(5:91,93,95)#include G1-G80 and "OD","pH","EC","CODt","CODs","CODb","Ammonium","PHOt" and "Cell_number"
min<-0.75

##L1P1
tiff("L1P1.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-I"&all$Reactor=="L1")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum=data.frame('L1P1'=counts)
correlation_list<-data.frame(list,"ID"=rep("L1P1",nrow(list)))
dev.off()

##L2P1
tiff("L2P1.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-I"&all$Reactor=="L2")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L2P1'=counts
list<-data.frame(list,"ID"=rep("L2P1",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L3P1
tiff("L3P1.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-I"&all$Reactor=="L3")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L3P1'=counts
list<-data.frame(list,"ID"=rep("L3P1",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L4P1
tiff("L4P1.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-I"&all$Reactor=="L4")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L4P1'=counts
list<-data.frame(list,"ID"=rep("L4P1",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L5P1
tiff("L5P1.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-I"&all$Reactor=="L5")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L5P1'=counts
list<-data.frame(list,"ID"=rep("L5P1",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##RP1
tiff("RP1.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-I"&all$Reactor=="R")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'RP1'=counts
list<-data.frame(list,"ID"=rep("RP1",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L1P2
tiff("L1P2.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-10%"&all$Reactor=="L1")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L1P2'=counts
list<-data.frame(list,"ID"=rep("L1P2",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L2P2
tiff("L2P2.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-10%"&all$Reactor=="L2")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L2P2'=counts
list<-data.frame(list,"ID"=rep("L2P2",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L3P2
tiff("L3P2.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-10%"&all$Reactor=="L3")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L3P2'=counts
list<-data.frame(list,"ID"=rep("L3P2",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L4P2
tiff("L4P2.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-10%"&all$Reactor=="L4")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L4P2'=counts
list<-data.frame(list,"ID"=rep("L4P2",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L5P2
tiff("L5P2.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-10%"&all$Reactor=="L5")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L5P2'=counts
list<-data.frame(list,"ID"=rep("L5P2",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##RP2
tiff("RP2.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-10%"&all$Reactor=="R")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'RP2'=counts
list<-data.frame(list,"ID"=rep("RP2",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L1P3
tiff("L1P3.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-50%"&all$Reactor=="L1")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L1P3'=counts
list<-data.frame(list,"ID"=rep("L1P3",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L2P3
tiff("L2P3.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-50%"&all$Reactor=="L2")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L2P3'=counts
list<-data.frame(list,"ID"=rep("L2P3",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L3P3
tiff("L3P3.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-50%"&all$Reactor=="L3")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L3P3'=counts
list<-data.frame(list,"ID"=rep("L3P3",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L4P3
tiff("L4P3.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-50%"&all$Reactor=="L4")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L4P3'=counts
list<-data.frame(list,"ID"=rep("L4P3",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L5P3
tiff("L5P3.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-50%"&all$Reactor=="L5")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L5P3'=counts
list<-data.frame(list,"ID"=rep("L5P3",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##RP3
tiff("RP3.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-50%"&all$Reactor=="R")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'RP3'=counts
list<-data.frame(list,"ID"=rep("RP3",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L1P4
tiff("L1P4.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-80%"&all$Reactor=="L1")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L1P4'=counts
list<-data.frame(list,"ID"=rep("L1P4",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L2P4
tiff("L2P4.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-80%"&all$Reactor=="L2")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L2P4'=counts
list<-data.frame(list,"ID"=rep("L2P4",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L3P4
tiff("L3P4.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-80%"&all$Reactor=="L3")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L3P4'=counts
list<-data.frame(list,"ID"=rep("L3P4",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L4P4
tiff("L4P4.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-80%"&all$Reactor=="L4")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L4P4'=counts
list<-data.frame(list,"ID"=rep("L4P4",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L5P4
tiff("L5P4.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-80%"&all$Reactor=="L5")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L5P4'=counts
list<-data.frame(list,"ID"=rep("L5P4",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##RP4
tiff("RP4.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Recycle-80%"&all$Reactor=="R")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'RP4'=counts
list<-data.frame(list,"ID"=rep("RP4",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L1P5
tiff("L1P5.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-II"&all$Reactor=="L1")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L1P5'=counts
list<-data.frame(list,"ID"=rep("L1P5",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L2P5
tiff("L2P5.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-II"&all$Reactor=="L2")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L2P5'=counts
list<-data.frame(list,"ID"=rep("L2P5",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L3P5
tiff("L3P5.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-II"&all$Reactor=="L3")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L3P5'=counts
list<-data.frame(list,"ID"=rep("L3P5",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L4P5
tiff("L4P5.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-II"&all$Reactor=="L4")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L4P5'=counts
list<-data.frame(list,"ID"=rep("L4P5",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##L5P5
tiff("L5P5.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-II"&all$Reactor=="L5")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'L5P5'=counts
list<-data.frame(list,"ID"=rep("L5P5",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

##RP5
tiff("RP5.tif",height=9,width=8,units = "in",res = 600,compression = "lzw")
s=which(all$Phase=="Closing-II"&all$Reactor=="R")
counts=COR(s,parameters,min)
list=cp_list(s,parameters,min)
cornum$'RP5'=counts
list<-data.frame(list,"ID"=rep("RP5",nrow(list)))
correlation_list<-rbind(correlation_list,list)
dev.off()

write.csv(cornum,"cornum.csv")
write.csv(correlation_list,"correlation_list.csv")

#compare between phases
cornum<-as.data.frame(t(cornum))
colnames(cornum)<-c("bb","ab")
cornum<-cornum[-1,]

cornum[grep("P1",rownames(cornum)),'Phase']<-'Insular I'
cornum[grep("P2",rownames(cornum)),'Phase']<-'Recycle 10%'
cornum[grep("P3",rownames(cornum)),'Phase']<-'Recycle 50%'
cornum[grep("P4",rownames(cornum)),'Phase']<-'Recycle 80%'
cornum[grep("P5",rownames(cornum)),'Phase']<-'Insular II'

cornum[grep("L1",rownames(cornum)),'Reactor']<-'L1'
cornum[grep("L2",rownames(cornum)),'Reactor']<-'L2'
cornum[grep("L3",rownames(cornum)),'Reactor']<-'L3'
cornum[grep("L4",rownames(cornum)),'Reactor']<-'L4'
cornum[grep("L5",rownames(cornum)),'Reactor']<-'L5'
cornum[grep("R",rownames(cornum)),'Reactor']<-'R'

cornum$bb<-as.character(cornum$bb)
cornum$bb<-as.numeric(cornum$bb)
cornum$ab<-as.character(cornum$ab)
cornum$ab<-as.numeric(cornum$ab)
cornum$Phase<-factor(cornum$Phase,levels=c('Insular I','Recycle 10%','Recycle 50%','Recycle 80%','Insular II'))
PA<-c('Insular I','Recycle 10%','Recycle 50%','Recycle 80%','Insular II')
cornum<-melt(cornum, id=c('Phase','Reactor'))

kruskal.test(value~Phase,cornum[which(cornum$Reactor!='R'&cornum$variable=='bb'),])
wilcox.test(cornum[which(cornum$Phase=='Insular I'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'],cornum[which(cornum$Phase=='Recycle 10%'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Insular I'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'],cornum[which(cornum$Phase=='Recycle 50%'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Recycle 50%'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'],cornum[which(cornum$Phase=='Recycle 10%'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Recycle 50%'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'],cornum[which(cornum$Phase=='Insular II'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Recycle 50%'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'],cornum[which(cornum$Phase=='Recycle 80%'&cornum$Reactor!='R'&cornum$variable=='bb'),'value'])
wilcox.test(cornum[which(cornum$Phase%in%PA[2:4]&cornum$Reactor!='R'&cornum$variable=='bb'),'value'],cornum[which(!cornum$Phase%in%PA[2:4]&cornum$Reactor!='R'&cornum$variable=='bb'),'value'])
wilcox.test(cornum[which(cornum$Phase%in%PA[2:4]&cornum$Reactor!='R'&cornum$variable=='bb'),'value'],cornum[which(!cornum$Phase%in%PA[2:5]&cornum$Reactor!='R'&cornum$variable=='bb'),'value'])

kruskal.test(value~Phase,cornum[which(cornum$Reactor!='R'&cornum$variable=='ab'),])
wilcox.test(cornum[which(cornum$Phase%in%PA[2:4]&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(!cornum$Phase%in%PA[2:4]&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Insular I'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Recycle 10%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Insular I'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Recycle 50%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Insular I'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Recycle 80%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Insular I'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Insular II'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Recycle 50%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Recycle 10%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Recycle 80%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Recycle 10%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Recycle 10%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Insular II'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Recycle 50%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Insular II'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Recycle 50%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Recycle 80%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase=='Recycle 80%'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(cornum$Phase=='Insular II'&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
wilcox.test(cornum[which(cornum$Phase%in%PA[2:4]&cornum$Reactor!='R'&cornum$variable=='ab'),'value'],cornum[which(!cornum$Phase%in%PA[2:4]&cornum$Reactor!='R'&cornum$variable=='ab'),'value'])
