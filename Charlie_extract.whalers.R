#Charlie's whaler sharks

setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data")

SMN.dusky=read.csv(file="SMN_donwloads_processed/dusky.csv")
SMN.bronzey=read.csv(file="SMN_donwloads_processed/bronzey.csv")

South.Oz=read.csv(file="C:/Matias/Data/Tagging/Acoustic_tagging/Other researcher's tags/Charlie's tags.csv")

SMN=rbind(SMN.dusky,SMN.bronzey)
SMN$Sex=as.character(SMN$Sex)
SMN$Sex=ifelse(SMN$Sex=="m","M",SMN$Sex)
SMN=SMN[,-match("X",names(SMN))]

SMN=subset(SMN,TagCode%in%South.Oz$Code2)
write.csv(SMN,"C:/Matias/Analyses/Acoustic_tagging/For Charlie/For.Charlie.csv",row.names=F)
