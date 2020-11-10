#Charlie's whaler sharks
library(data.table)
library(dplyr)

setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data")
Detections=fread("Detections.csv",data.table=FALSE)

South.Oz=read.csv(file="C:/Matias/Data/Tagging/Acoustic_tagging/Other researcher's tags/Charlie's tags.csv")
These.tags=c(South.Oz$Code2,23293,23294)
WA.tagged_detected.in.SA=c(23303,27698,29542,29454,29587,30870,30891,30894,30992,31000,31003)
These.tags=c(These.tags,WA.tagged_detected.in.SA)

  #tags
TAGS=read.csv("TAGS.csv",stringsAsFactors=F)%>%
  filter(Species2%in%c('bronze whaler','Dusky') & Project.rel=="SMN")

  #detections
SMN.tags.bronze.dusky=Detections%>%
  filter(Species%in%c('bronze whaler','Dusky') & Project=="SMN")%>%
  distinct(TagCode,.keep_all = T)
SMN=subset(Detections,Species%in%c('bronze whaler','Dusky'))
SMN$Sex=as.character(SMN$Sex)
SMN$Sex=ifelse(SMN$Sex=="m","M",SMN$Sex)
SMN=subset(SMN,TagCode%in%These.tags & Longitude<=129)

#Export stuff
setwd("C:/Matias/Analyses/Acoustic_tagging/For Charlie/Data")
write.csv(unique(c(TAGS$Code2,SMN.tags.bronze.dusky$TagCode)),"For.Charlie_WA.tags.csv",row.names=F)
write.csv(SMN,"For.Charlie.csv",row.names=F)
