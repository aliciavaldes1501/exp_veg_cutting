#####################################################################################
###############·####### Data exploration ants Remmene 2016 ##########################
#####################################################################################

library(ggplot2)
library(tidyr)
library(reshape2)
library(ggthemes)

ants_R<-read.table("./data/raw/ants_Remmene_2016.txt",header=T,sep="\t",dec=".")
ants_R$date<-NULL
head(ants_R)

ants_R_long<-gather(ants_R,species,number,Mrubra,Mscabrinodis,Mruginodis,Mschencki,others,factor_key=TRUE)
ants_R_long$plot_id<-as.factor(ants_R_long$plot_id)
ants_R_long$point_id<-as.factor(ants_R_long$point_id)
ants_R_long$pres<-ifelse(ants_R_long$number>0,1,0)
head(ants_R_long)

pdf("./results/figures/ant_counts_R.pdf", family="Times")
#total number of individuals of each species
ggplot(ants_R_long,aes(x=species,y=number))+geom_bar(stat="identity")+guides(fill=FALSE)+
  xlab("Species")+ylab("Total number of individuals")+
  theme_base()+theme(plot.background=element_rect(fill="white", colour=NA))
#number of points with presence of each species
ggplot(ants_R_long,aes(x=species,y=pres))+geom_bar(stat="identity")+guides(fill=FALSE)+
  xlab("Species")+ylab("Number of points with presence")+
  theme_base()+theme(plot.background=element_rect(fill="white", colour=NA))
dev.off()
