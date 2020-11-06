#Charlie's whaler sharks

setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data")
Detections=fread("Detections.csv",data.table=FALSE)


# setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data")
# 
# SMN.dusky=read.csv(file="SMN_donwloads_processed/dusky.2018-10-19.csv")
# SMN.bronzey=read.csv(file="SMN_donwloads_processed/copper.2018-10-19.csv")

South.Oz=read.csv(file="C:/Matias/Data/Tagging/Acoustic_tagging/Other researcher's tags/Charlie's tags.csv")
These.tags=c(South.Oz$Code2,23303, 27698, 31003,29542,
                 30992, 29587, 30891, 26438, 30894, 31003, 31000)

SMN=subset(Detections,Species%in%c('bronze whaler','Dusky'))
SMN$Sex=as.character(SMN$Sex)
SMN$Sex=ifelse(SMN$Sex=="m","M",SMN$Sex)
SMN=subset(SMN,TagCode%in%These.tags & Longitude<=129)

write.csv(SMN,"C:/Matias/Analyses/Acoustic_tagging/For Charlie/Data/For.Charlie.csv",row.names=F)
