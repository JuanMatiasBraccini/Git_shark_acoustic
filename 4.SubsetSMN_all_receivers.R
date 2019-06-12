#SCRIPT FOR MANIPULATING SMN RECEIVERS INFOR

#1. lOAD DATA
setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/SMN_receiver_location")
SMN=read.csv(file="ReceiverStationDetailsReport.27_07_2015.csv")
#SMN=read.csv(file="ReceiverStationDetailsReport.10_04_2015.csv")
#SMN=read.csv(file="ReceiverStationDetailsReport.25_09_2012.csv")


#2. MANIPULATE DATA
THIS=match(c("Name1","Latitude1","Longitude1","SerialNumber"),colnames(SMN))
#THIS=match(c("Name1","Latitude1","Longitude1","ï..SerialNumber"),colnames(SMN))
colnames(SMN)[THIS]=c("station.name","latitude","longitude","receiver.ID")



SMN$Project="SMN"

SMN$EffectiveStartDate=as.POSIXlt(strptime(as.character(SMN$EffectiveStartDate1),format='%d-%B-%Y'))
SMN$EffectiveStartDate=format(SMN$EffectiveStartDate, "%d-%m-%y")

SMN$EffectiveEndDate=as.POSIXlt(strptime(as.character(SMN$EffectiveEndDate2),format='%d-%B-%Y'))
SMN$EffectiveEndDate=format(SMN$EffectiveEndDate, "%d-%m-%y")

keep.these.vars=match(c("receiver.ID","station.name","latitude","longitude",
                        "EffectiveStartDate","EffectiveEndDate","Project"),colnames(SMN))

SMN=SMN[,keep.these.vars]


#3. EXPORT FILE
write.csv(SMN,file="SMN_receivers.manipulated.csv")




