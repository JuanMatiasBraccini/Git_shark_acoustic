#Charlie's whaler sharks
library(data.table)
library(dplyr)

send.only.shared=FALSE
send.all=TRUE
UTC.WA=8
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

setwd(handl_OneDrive("Data/Tagging/Acoustic_tagging/Acoustic_tagging_data"))
Detections=fread("Detections.csv",data.table=FALSE)

South.Oz=read.csv(file=handl_OneDrive("Data/Tagging/Acoustic_tagging/Other researcher's tags/Charlie's tags.csv"))
These.tags=c(South.Oz$Code2,23293,23294)
WA.tagged_detected.in.SA=c(23303,27698,29542,29454,29587,30870,30891,30894,30992,31000,31003)
These.tags=c(These.tags,WA.tagged_detected.in.SA)

  #tags
TAGS=read.csv("TAGS.csv",stringsAsFactors=F)%>%
  filter(Species2%in%c('bronze whaler','Dusky') & Project.rel=="SMN")
SMN.tags.bronze.dusky=Detections%>%
  filter(Species%in%c('bronze whaler','Dusky') & Project=="SMN")%>%
  distinct(TagCode,.keep_all = T)
tags=unique(c(TAGS$Code2,SMN.tags.bronze.dusky$TagCode))

  #Detections
SMN=Detections%>%
      filter(Species%in%c('bronze whaler','Dusky'))%>%
  mutate(Sex=as.character(Sex),
         Sex=ifelse(Sex=="m","M",Sex),
         DateTime.local=as.POSIXlt(DateTime.local,  format="%Y-%m-%d %H:%M",tz="GMT"),
         DateTime=DateTime.local-(UTC.WA*3600))
if(send.only.shared) SMN=subset(SMN,TagCode%in%These.tags & Longitude<=129)
if(send.all)   SMN=subset(SMN, Longitude<=129)

#Export stuff
setwd(handl_OneDrive("Analyses/Acoustic_tagging/For Charlie/Data"))
write.csv(tags,"For.Charlie_WA.tags.csv",row.names=F)
write.csv(TAGS,"For.Charlie_WA.tags_all.info.csv",row.names=F)
write.csv(SMN,"For.Charlie.csv",row.names=F)
