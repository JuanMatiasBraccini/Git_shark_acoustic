#Charlie's whaler sharks
library(data.table)
library(dplyr)

setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data")
Detections=fread("Detections.csv",data.table=FALSE)


# setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data")
# 
# SMN.dusky=read.csv(file="SMN_donwloads_processed/dusky.2018-10-19.csv")
# SMN.bronzey=read.csv(file="SMN_donwloads_processed/copper.2018-10-19.csv")

South.Oz=read.csv(file="C:/Matias/Data/Tagging/Acoustic_tagging/Other researcher's tags/Charlie's tags.csv")
These.tags=c(South.Oz$Code2)
WA.tagged_detected.in.SA=c(23293,23294,23303,26438,27698,29542,29454,29587,30717,30870,30891,30894,30992,31000,31003)
These.tags=c(These.tags,WA.tagged_detected.in.SA)


  #tags
TAGS=read.csv("TAGS.csv",stringsAsFactors=F)%>%
  filter(Species2%in%c('bronze whaler','Dusky') & Project.rel=="SMN")
SMN.tags.bronze.dusky=Detections%>%
                      filter(Species%in%c('bronze whaler','Dusky') & Project=="SMN")%>%
                      distinct(TagCode,.keep_all = T)

write.csv(unique(c(TAGS$Code2,SMN.tags.bronze.dusky$TagCode)),"C:/Matias/Analyses/Acoustic_tagging/For Charlie/Data/For.Charlie_WA.tags.csv",row.names=F)


cc=c(23293,26438,30717)
subset(TAGS,Code2%in%cc)

  #detections
SMN=subset(Detections,Species%in%c('bronze whaler','Dusky'))
SMN$Sex=as.character(SMN$Sex)
SMN$Sex=ifelse(SMN$Sex=="m","M",SMN$Sex)
SMN=subset(SMN,TagCode%in%These.tags & Longitude<=129)
write.csv(SMN,"C:/Matias/Analyses/Acoustic_tagging/For Charlie/Data/For.Charlie.csv",row.names=F)
