#load library
library(stringr)
library(reshape2)  
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(vegan)

#set working directory
setwd(".../code_lis_210406")

#1.bar plot of relative abundance of Class or Genus in whole communities and sorted gates
#set color for plotting
S2 <- brewer.pal(7, "Set2")
colors1 <- c("#50a33e", "#7eb38e", "#b1f41b", "#0e8121", "#5edbff", "#0097c3", "#0a3863", "#56B4E9", "#b2f5ff", "#0e5c63", "#1b8085", "#3aa6a9", "#ffb3b3", "#ff7b7b")
colors<- c( S2, colors1)

#1.1 sorted gates_class level
#read data, sorted gates_class level
table_sg_c<-read.delim('seq table/10_Asv_sg_after_curation_table_class.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)
tax_sg_c<-read.delim('seq table/11_Taxonomy_sg_after_curation_table_class.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)
meta_sg<-read.delim('seq table/02_Metadata_sorted_gates.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)

#leave out unuseful samples(sample 10,11,36),order Class according to abundance
table_sg_c_d<-table_sg_c
table_sg_c_d<-table_sg_c_d[,-c(10,11,46)]
rownames(table_sg_c_d)<-tax_sg_c[,'Class']
table_sg_c_d$sum <- rowSums(table_sg_c_d)
table_sg_c_d <- table_sg_c_d[order(table_sg_c_d$sum, decreasing = TRUE), ]
meta_sg_d<-meta_sg
meta_sg_d<-meta_sg_d[-c(10,11,46),]

#only keep Class ranking in top10, Class out of top10 sum up as 'Others'
table_sg_c_top10 <- table_sg_c_d[1:10, -ncol(table_sg_c_d)]
table_sg_c_top10['Others', ] <- 1 - colSums(table_sg_c_top10)  
table_sg_c_top10<-t(table_sg_c_top10)

#extract gates name
gate<-str_split(meta_sg_d$Type,"_")
gate<-sapply(gate,"[",2) 

#combine data for plotting
table_sg_c_top10<-cbind(gate,meta_sg_d[,c(3,4,7)],table_sg_c_top10)
table_sg_c_top10$Days <- factor(table_sg_c_top10$Days , levels = 1:107)
table_sg_c_top10$gate <- factor(table_sg_c_top10$gate , levels = paste0('G',1:80))
table_sg_c_top10 <- table_sg_c_top10[order(table_sg_c_top10$Days,table_sg_c_top10$gate, decreasing = FALSE), ]
table_sg_c_top10$Type <- factor(table_sg_c_top10$Type, levels = unique(table_sg_c_top10$Type))

data_sg_c<-melt(table_sg_c_top10, id=c('gate', 'Reactor','Days','Type'))
data_sg_c$variable <- factor(data_sg_c$variable, levels = rev(colnames(table_sg_c_top10)[5:15]))
data_sg_c$Type <- factor(data_sg_c$Type, levels = as.character(unique(table_sg_c_top10$Type)))

#figure S13.5 bar plot of Class relative abundance, sorted gates 
bre<-unique(table_sg_c_top10$Type)
color_sg_c<-c('blue', 'orange', 'green', 'yellow', 'red', 'hotpink', 'cyan','purple', 'burlywood1', 'skyblue','gray')
figs13.5 <- ggplot(data_sg_c, aes(Type, value*100, fill = variable)) +
 geom_col(position = 'stack', width = 0.6) +
  facet_wrap(~Reactor, ncol = 1) + 
  annotate("rect", xmin = 0, xmax = bre[6], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[6], xmax = bre[23], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[23], xmax = bre[35], ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_col(position = 'stack', width = 0.6)+
  scale_fill_manual(name='Class',values =  rev(color_sg_c),guide = guide_legend(reverse = TRUE,ncol = 4)) +
  labs(x = 'sorted subcommunities', y = 'relative abundance (%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(face='bold',size = 12),strip.background = element_blank()) +    
  theme(axis.text = element_text(size = 10), axis.text.x =element_text(angle=90,vjust = 0.5),axis.title =element_text(face='bold',size=12), legend.title = element_text(size = 12), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4,"cm"),legend.position = "bottom",legend.direction = "vertical")
figs13.5


#1.2 whole community_class level
#read data, whole community_class level
table_wc_c<-read.delim('seq table/10_Asv_wc_after_curation_class_table.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)
tax_wc_c<-read.delim('seq table/11_Taxonomy_after_curation_class_table.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)
meta_wc<-read.delim('seq table/02_Metadata_whole_community_ZS.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)

#only keep Class ranking in top10, Class out of top10 sum up as 'Others'
rownames(table_wc_c)<-tax_wc_c[,'Class']
table_wc_c$sum <- rowSums(table_wc_c)
table_wc_c <- table_wc_c[order(table_wc_c$sum, decreasing = TRUE), ]
table_wc_c_top10 <- table_wc_c[1:10, -ncol(table_wc_c)]
table_wc_c_top10['Others', ] <- 1 - colSums(table_wc_c_top10)  
table_wc_c_top10<-t(table_wc_c_top10)

#combine data for plotting
table_wc_c_top10<-cbind(meta_wc[,c(2,3)],table_wc_c_top10)
table_wc_c_top10$Days <- factor(table_wc_c_top10$Days , levels = 0:107)
table_wc_c_top10 <- table_wc_c_top10[order(table_wc_c_top10$Days, decreasing = FALSE), ]
table_wc_c_top10 <-rbind(table_wc_c_top10[1,], table_wc_c_top10[1,], table_wc_c_top10[1,],table_wc_c_top10[1,],table_wc_c_top10)
table_wc_c_top10$Reactor[1:5]<-c('L1','L2','L3','L4','L5')

data_wc_c<-melt(table_wc_c_top10, id=c('Reactor','Days'))
data_wc_c$variable <- factor(data_wc_c$variable, levels = rev(colnames(table_wc_c_top10)[3:13]))
data_wc_c$Days <- factor(data_wc_c$Days , levels = 0:107)

#figure s13.1 bar plot of Class relative abundance, whole community
#expand color plate for Class
add_c<-union(colnames(table_sg_c_top10)[5:15],colnames(table_wc_c_top10)[3:13])
color_c<-rbind(add_c,c(color_sg_c,colors[1:3]))
colnames(color_c)<-color_c[1,]
color_c<-color_c[-1,]

bre<-unique(table_wc_c_top10$Days)
figs13.1 <- ggplot(data_wc_c, aes(Days, value*100, fill = variable)) +
  geom_col(position = 'stack', width = 0.6) +
  facet_wrap(~Reactor, ncol = 1) + 
  annotate("rect", xmin = bre[5], xmax = bre[9], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[9], xmax = bre[13], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[13], xmax = bre[16], ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_col(position = 'stack', width = 0.6)+
  scale_fill_manual(name='Class',values =  color_c,guide = guide_legend(reverse = TRUE,ncol = 4)) +
  labs(x = 'time (d)', y = 'relative abundance (%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(face='bold',size = 12),strip.background = element_blank()) +    
  theme(axis.text = element_text(size = 10),axis.title =element_text(face='bold',size=12), legend.title = element_text(size = 12), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4,"cm"),legend.position = "bottom",legend.direction = "vertical")
figs13.1


#1.3 sorted gate_genus level
#read data, sorted gate_genus level
table_sg_g<-read.delim('seq table/04_Asv_sg_after_curation_table_genus.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)
tax_sg_g<-read.delim('seq table/05_Taxonomy_sg_after_curation_table_genus.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)

##leave out unuseful samples(sample 10,11,36),order Genera according to abundance
table_sg_g_d<-table_sg_g[,-c(10,11,46)]
table_sg_g_d$sum <- rowSums(table_sg_g_d)
tax_sg_g_d <- tax_sg_g[order(table_sg_g_d$sum, decreasing = TRUE), ]
table_sg_g_d <- table_sg_g_d[order(table_sg_g_d$sum, decreasing = TRUE), ]

#only keep Genera ranking in top20, Genera out of top20 sum up as 'Others'
table_sg_g_top20 <- table_sg_g_d[1:20, -ncol(table_sg_g_d)]
table_sg_g_top20['Others', ] <- 1 - colSums(table_sg_g_top20) 

#Genera absent for Genus assignment are assigned with Family or Order name
rownames(table_sg_g_top20)[1:2]<-tax_sg_g_d[1:2,'Genus']
rownames(table_sg_g_top20)[3]<-paste0(tax_sg_g_d[3,'Family'],' (Family)')
rownames(table_sg_g_top20)[4]<-paste0(tax_sg_g_d[4,'Order'],' (Order)')
rownames(table_sg_g_top20)[5:20]<-tax_sg_g_d[5:20,'Genus']
table_sg_g_top20<-t(table_sg_g_top20)

#extra gate name
gate<-str_split(meta_sg_d$Type,"_")
gate<-sapply(gate,"[",2)

#combine data for plotting
table_sg_g_top20<-cbind(gate,meta_sg_d[,c(3,4,7)],table_sg_g_top20)
table_sg_g_top20$Days <- factor(table_sg_g_top20$Days , levels = 1:107)
table_sg_g_top20$gate <- factor(table_sg_g_top20$gate , levels = paste0('G',1:80))
table_sg_g_top20 <- table_sg_g_top20[order(table_sg_g_top20$Days,table_sg_g_top20$gate, decreasing = FALSE), ]
table_sg_g_top20$Type <- factor(table_sg_g_top20$Type, levels = unique(table_sg_g_top20$Type))

data_sg_g<-melt(table_sg_g_top20, id=c('gate', 'Reactor','Days','Type'))
data_sg_g$variable <- factor(data_sg_g$variable, levels = rev(colnames(table_sg_g_top20)[5:25]))
data_sg_g$Type <- factor(data_sg_g$Type, levels = as.character(unique(table_sg_g_top20$Type)))

#figure S13.6 bar plot of Genus relative abundance, sorted gates
bre<-unique(table_sg_g_top20$Type)
figs13.6 <- ggplot(data_sg_g, aes(Type, value*100, fill = variable)) +
  geom_col(position = 'stack', width = 0.6) +
  facet_wrap(~Reactor, ncol = 1) + 
  annotate("rect", xmin = 0, xmax = bre[6], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[6], xmax = bre[23], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[23], xmax = bre[35], ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_col(position = 'stack', width = 0.6)+
  scale_fill_manual(name='Genus',values =  rev(c(colors[1:20],'grey')),guide = guide_legend(reverse = TRUE,ncol = 4)) +
  labs(x = 'sorted subcommunities', y = 'relative abundance (%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(face='bold',size = 12),strip.background = element_blank()) +    
  theme(axis.text = element_text(size = 10), axis.text.x =element_text(angle=90,vjust = 0.5),axis.title =element_text(face='bold',size=12), legend.title = element_text(size = 12), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4,"cm"),legend.position = "bottom",legend.direction = "vertical")
figs13.6


#1.4 whole community_genus level
#read data, whole community_genus level
table_wc_g<-read.delim('seq table/04_Asv_wc_after_curation_genus_table.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)
tax_wc_g<-read.delim('seq table/05_Taxonomy_after_curation_genus_table.txt', row.names = 1,sep = '\t',stringsAsFactors = FALSE, check.names = FALSE)

#order Genera according to abundance
table_wc_g$sum <- rowSums(table_wc_g)
tax_wc_g <- tax_wc_g[order(table_wc_g$sum, decreasing = TRUE), ]
table_wc_g <- table_wc_g[order(table_wc_g$sum, decreasing = TRUE), ]

#only keep Genera ranking in top20, Genera out of top20 sum up as 'Others'
table_wc_g_top20 <- table_wc_g[1:20, -ncol(table_wc_g)]
table_wc_g_top20['Others', ] <- 1 - colSums(table_wc_g_top20)  

#Genera absent for Genus assignment are assigned with Family 
rownames(table_wc_g_top20)[1]<-paste0(tax_wc_g[1,'Family'],' (Family)')
rownames(table_wc_g_top20)[2:20]<-tax_wc_g[2:20,'Genus']
table_wc_g_top20<-t(table_wc_g_top20)

#combine data for plotting
table_wc_g_top20<-cbind(meta_wc[,c(2,3)],table_wc_g_top20)
table_wc_g_top20$Days <- factor(table_wc_g_top20$Days , levels = 0:107)
table_wc_g_top20 <- table_wc_g_top20[order(table_wc_g_top20$Days, decreasing = FALSE), ]
table_wc_g_top20 <-rbind(table_wc_g_top20[1,], table_wc_g_top20[1,], table_wc_g_top20[1,],table_wc_g_top20[1,],table_wc_g_top20)
table_wc_g_top20$Reactor[1:5]<-c('L1','L2','L3','L4','L5')

data_wc_g<-melt(table_wc_g_top20, id=c('Reactor','Days'))
data_wc_g$variable <- factor(data_wc_g$variable, levels = rev(colnames(table_wc_g_top20)[3:23]))
data_wc_g$Days <- factor(data_wc_g$Days , levels = 0:107)

#expand color plate for Genus
add_g<-union(colnames(table_sg_g_top20)[5:25],colnames(table_wc_g_top20)[3:23])
color_g<-rbind(add_g,c(colors[1:20],'grey',colors[21], brewer.pal(7, "Set1")[c(1,4:7)],brewer.pal(5, "Pastel2")[c(1:3,5)]))
colnames(color_g)<-color_g[1,]
color_g<-color_g[-1,]

#figure S13.2 bar plot of Genus relative abundance, whole community
bre<-unique(table_wc_g_top20$Days)
figs13.2 <- ggplot(data_wc_g, aes(Days, value*100, fill = variable)) +
  geom_col(position = 'stack', width = 0.6) +
  facet_wrap(~Reactor, ncol = 1) + 
  annotate("rect", xmin = bre[5], xmax = bre[9], ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = bre[9], xmax = bre[13], ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = bre[13], xmax = bre[16], ymin = -Inf, ymax = Inf, alpha = .35)+
  geom_col(position = 'stack', width = 0.6)+
  scale_fill_manual(name='Genus',values =  color_g,guide = guide_legend(reverse = TRUE,ncol = 4)) +
  labs(x = 'time (d)', y = 'relative abundance (%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(face='bold',size = 12),strip.background = element_blank()) +    
  theme(axis.text = element_text(size = 10),axis.title =element_text(face='bold',size=12), legend.title = element_text(size = 12), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4,"cm"),legend.position = "bottom",legend.direction = "vertical")
figs13.2


#2.comparision analysis of whole community taxonomic and flow-cytometric compositions
#2.1 nmds, whole community at Genus level
#combine data for analysis
table_wc_g_t<-cbind(meta_wc[,c(2,3)],t(table_wc_g)[-ncol(table_wc_g),])
table_wc_g_t<-rbind(table_wc_g_t[96,],table_wc_g_t[96,],table_wc_g_t[96,],table_wc_g_t[96,],table_wc_g_t[96,],table_wc_g_t[1:95,])
table_wc_g_t$Reactor[1:5]<-c('L1','L2','L3','L4','L5')
table_wc_g_t<-table_wc_g_t[order( table_wc_g_t$Reactor,table_wc_g_t$Days,decreasing = FALSE), ]
table_wc_g_t$Reactor<-factor(table_wc_g_t$Reactor,levels = c('L1','L2','L3','L4','L5','R'))
rownames(table_wc_g_t)<-1:100

#nmds analysis
BC <- metaMDS(table_wc_g_t[,3:161], distance="bray", autotransform=FALSE, zerodist="add",try=100,trymax = 200)

#figure s13.3 nmds plot
colplate_reactor=c("1"='#d53e4f',"2"='#fc8d59',"3"='#fee08b',"4"='#e6f598',"5"='#99d594',"6"='#3288bd')
col_reactor=colplate_reactor[as.factor(table_wc_g_t$Reactor)]

layout(matrix(c(1,2,3,4,5,6,6,6,6,6),2,5,byrow=TRUE),widths=c(1,1,1,1,1),height=c(2,1.2))
par(mai=c(0.3,0.2,0.25,0),omi=c(0,0.3,0.2,0.1),xpd=FALSE)
par(ps = 10, cex.main=1.8, cex.axis = 1,cex.lab=1.25)

#phase1-Insular I
plot(BC, type="n",xlab='',ylab='')
title("Insular I")

#grey points in other phases
s=c(1:100)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#colored points in the present phase
s=1:4#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+table_wc_g_t[s,]$Days/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=18:21#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+table_wc_g_t[s,]$Days/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=35:38#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+table_wc_g_t[s,]$Days/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=53:56#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+table_wc_g_t[s,]$Days/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=69:72#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+table_wc_g_t[s,]$Days/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=87:88#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+table_wc_g_t[s,]$Days/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

arrows(BC$points[1,1],BC$points[1,2],BC$points[2,1],BC$points[2,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")
arrows(BC$points[1,1],BC$points[1,2],BC$points[19,1],BC$points[19,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")
arrows(BC$points[1,1],BC$points[1,2],BC$points[36,1],BC$points[36,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")
arrows(BC$points[1,1],BC$points[1,2],BC$points[54,1],BC$points[54,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")
arrows(BC$points[1,1],BC$points[1,2],BC$points[70,1],BC$points[70,2], length = 0.08, angle = 15, lwd=1.5, col ="grey70")
points(BC$points[1,1],BC$points[1,2],col="black",type = "p",pch=19)

#phase2-RC10
plot(BC, type="n",xlab='',ylab='',yaxt="n")
title("RC10")

#grey points in other phases
s=c(1:100)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#colored points in the present phase
s=4:7#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=21:25#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=38:41#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=56:59#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=72:76#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=88:91#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-25)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

#phase3-RC50
plot(BC, type="n",xlab='',ylab='',yaxt="n")
title("RC50")

#grey points in other phases
s=c(1:100)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#colored points in the present phase
s=7:10#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=25:28#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=41:45#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=59:62#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=76:80#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=91:94#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-46)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

#phase4-RC80
plot(BC, type="n",xlab='',ylab='',yaxt="n")
title("RC80")

#grey points in other phases
s=c(1:100)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#colored points in the present phase
s=10:13#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=28:31#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=45:48#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=62:65#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=80:83#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=94:97#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-63)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

#phase5-Insular II
plot(BC, type="n",xlab='',ylab='',yaxt="n")
title("Insular II")

#grey points in other phases
s=c(1:100)
points(BC$points[s,1],BC$points[s,2],col="grey70",type = "p",pch=19,cex=0.5)

#colored points in the present phase
s=13:17#L1
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=31:34#L2
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=48:52#L3
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=65:68#L4
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

s=83:86#L5
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

points(BC$points[83,1],BC$points[83,2],col=col_reactor[83],cex=0.5+(table_wc_g_t[83,]$Days-88)/110*3,pch=19)
points(BC$points[84,1],BC$points[84,2],col=col_reactor[84],cex=0.5+(table_wc_g_t[84,]$Days-88)/110*3,pch=19)
points(BC$points[85,1],BC$points[85,2],col=col_reactor[85],cex=0.5+(table_wc_g_t[85,]$Days-88)/110*3,pch=19)
points(BC$points[86,1],BC$points[86,2],col=col_reactor[86],cex=0.5+(table_wc_g_t[86,]$Days-88)/110*3,pch=19)

s=97:100#R
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],cex=0.5+(table_wc_g_t[s,]$Days-88)/110*3,pch=19)
points(BC$points[s,1],BC$points[s,2],col=col_reactor[s],type = "l")

#legend
plot.new()
legend(title="",x="center",y="top",ncol=4,legend=c("L1","L2","L3","L4","L5","R","preculture"),col=c("1"='#d53e4f',"2"='#fc8d59',"3"='#fee08b',"4"='#e6f598',"5"='#99d594',"6"='#3288bd',"7"='black'),pch=c(19,19,19,19,19,19,19),lwd=1.5,title.adj = 0,box.col = "white",cex=1.25)
mtext("reactors", side = 3, line = -15,adj=0.25,outer = TRUE,cex=1.25)
mtext("NMDS2", side = 2, line = 1,adj=0.72,outer = TRUE,cex=1.25)
mtext("NMDS1", side = 3, line = -14,adj=0.53,outer = TRUE,cex=1.25)

dev.off()

#2.2 whole community drift on Genus level 
#transform genus table to dominance (1) and undominance (0) table
tr<-1/(ncol(table_wc_g_t)-2)#treshold for dominant genra
table_wc_g_tr<-(table_wc_g_t[,3:ncol(table_wc_g_t)]>tr)*1
table_wc_g_tr<-cbind(table_wc_g_t[,1:2],table_wc_g_tr)

#calculate intra-community beta diversity based on Genus level
table_wc_g_tr$intra=0
for(i in 2:nrow(table_wc_g_tr)){
  for (j in 3:(ncol(table_wc_g_tr)-1)){
    if (table_wc_g_tr[i-1,j]*table_wc_g_tr[i,j]==0&!table_wc_g_tr[i-1,j]+table_wc_g_tr[i,j]==0){
      table_wc_g_tr$intra[i]<-table_wc_g_tr$intra[i]+1
    }
  }
}
table_wc_g_tr[which(table_wc_g_tr$Days==0),'intra']<-NA
intra_g<-table_wc_g_tr[,c(1,2,ncol(table_wc_g_tr))]

#for comparison, recalculate intra-community beta diversity based on subcommunity table, using the same pairwise samples as above
table_sc<-read.table("RA.txt",header=TRUE,na.strings = "NA",blank.lines.skip = FALSE)
table_sc<-table_sc[,c(2,3,5:84)]
table_sc<-rbind(table_sc[1,],table_sc[1,],table_sc[1,],table_sc[1,],table_sc)
table_sc$Reactor[1:5]<-reactorname[1:5]
table_sc_tr<-(table_sc[,3:ncol(table_sc)]>0.0125)*1
table_sc_tr<-cbind(table_sc[,1:2],table_sc_tr)
table_sc_tr$intra<-NA
table_sc_tr$Reactor<-droplevels(table_sc_tr$Reactor)

for (i in 1:nrow(intra_g)){
  table_sc_tr[which(table_sc_tr$Reactor==intra_g$Reactor[i]&table_sc_tr$Time==intra_g$Days[i]),'intra']<-0
}
table_sc_tr<-na.omit(table_sc_tr)
table_sc_tr<-table_sc_tr[order(table_sc_tr$Reactor),]

for(i in 2:nrow(table_sc_tr)){
  for (j in 3:(ncol(table_sc_tr)-1)){
    if (table_sc_tr[i-1,j]*table_sc_tr[i,j]==0&!table_sc_tr[i-1,j]+table_sc_tr[i,j]==0){
      table_sc_tr$intra[i]<-table_sc_tr$intra[i]+1
    }
  }
}

table_sc_tr[which(table_sc_tr$Time==0),'intra']<-NA
intra_sc<-table_sc_tr[,c(1,2,ncol(table_sc_tr))]
colnames(intra_sc)[2]<-'Days'
intra<-rbind(data.frame(intra_g, group=rep(1,nrow(intra_g))),data.frame(intra_sc, group=rep(2,nrow(intra_sc))))
intra$group<-as.factor(intra$group)

#figure S13.4 intra-community beta diversity comparison 
ann_text <- data.frame(x=c(47,47,61,61,100,100),y=c(35,35,25,25,25,25),
                       Reactor = factor(reactorname[c(2,5,3,5,1,3)],levels = reactorname),gate=c('drift1','drift1','drift2','drift2','drift3','drift3'))
ann_line<-data.frame(x1=c(47,47,61,61,100,100),x2=c(47,47,61,61,100,100),y1=c(32,32,22,22,22,22),y2=c(25,25,15,15,15,15),
                     Reactor = factor(reactorname[c(2,5,3,5,1,3)],levels = reactorname))


figs13.4=ggplot(data=intra,aes(x=Days,y=intra,color=Reactor,linetype = group,group=group,shape=group))+
  annotate("rect", xmin = 26, xmax = 47, ymin = -Inf, ymax = Inf, alpha = .1)+
  annotate("rect", xmin = 47, xmax = 64, ymin = -Inf, ymax = Inf, alpha = .25)+
  annotate("rect", xmin = 64, xmax = 89, ymin = -Inf, ymax = Inf, alpha = .35)+
  facet_wrap(.~Reactor,scales="free",nrow=2)+
  geom_hline(yintercept=8.76, linetype="dashed", color = "black")+
  geom_line(size=0.5)+geom_point(size=1,stroke =1)+
scale_color_manual(values =c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'))+
  scale_linetype_manual(values = c(1,2))+
  scale_shape_manual(values=c(19,1))+
  scale_x_continuous(limits = c(0,110),expand = c(0,1),breaks=c(1,26,47,64,89,107))+xlab("time (d)")+
  scale_y_continuous(limits=c(0,40),expand=c(0.03,0.025),breaks=c(0,10,20,30,40),labels=c(0,10,20,30,40),name = expression(bold(paste("intra-community ",beta,"-diversity "))) )+
  theme(axis.text = element_text(size = 10))+
  theme(panel.background  = element_rect(color = "black",fill = 'white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text = element_text(face="bold",size=12),axis.title = element_text(face="bold",size=12))+
  theme(legend.position = "none")+theme(strip.background = element_blank())
figs13.4=figs13.4+geom_segment(data=ann_line,aes(x=x1,xend=x2,y=y1,yend=y2),parse = TRUE,inherit.aes = FALSE,arrow=arrow(length=unit(0.2,"cm")),show_guide=F)+
  geom_text(data=ann_text,aes(x=x,y=y,label=gate,size=3),parse = TRUE,inherit.aes = FALSE,show_guide=F)
figs13.4


#figure S6.1 bar plot of Species relative abundance of Zymo community
#read data
table_zm_s<-read.csv( "seq table/16_asv_rel_abd_zymo_species.csv")
table_zm_g<-read.csv("seq table/17_asv_rel_abd_zymo_genus.csv")
meta_zm<-read.delim("02_metadata_zymo_and_negative_controls.tsv")
colnames(table_zm_s)[2:5]<-c('Zymo DNA','Zymo community1','Zymo community2','Zymo community3')
colnames(table_zm_g)[2:5]<-c('Zymo DNA','Zymo community1','Zymo community2','Zymo community3')

data_zm_s<-melt(table_zm_s, id=c('Species','Row.names'))
data_zm_g<-melt(table_zm_g, id=c('Genus','Row.names'))

figs6.1 <- ggplot(data_zm_s, aes(variable, value, fill = Species)) +
  geom_col(position = 'stack', width = 0.6) +
  geom_col(position = 'stack', width = 0.6)+
  scale_fill_manual(name='Species',values = brewer.pal(8,"Set1"),guide = guide_legend(ncol = 1)) +
  labs(x = 'positive controls', y = 'relative abundance (%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(face='bold',size = 12),strip.background = element_blank()) +    
  theme(axis.text = element_text(size = 10), axis.text.x =element_text(angle=40,vjust = 0.5),axis.title =element_text(face='bold',size=12), legend.title = element_text(size = 12), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4,"cm"),legend.direction = "vertical")
figs6.1

