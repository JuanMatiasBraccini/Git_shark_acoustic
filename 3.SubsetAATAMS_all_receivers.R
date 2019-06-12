#SCRIPT FOR MANIPULATING AATAMS RECEIVERS INFOR

#1. lOAD DATA
setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/AATAMS_receiver_location")
#setwd("M:/Fisheries Research/FinFish/Shark/Braccini/Acoustic_tagging_data/AATAMS_receiver_location")

AATAMS=read.csv(file="AATAMS_receivers_25_09_2012.csv")

#2. MANIPULATE DATA
THIS=match(c("Station","code_name"),colnames(AATAMS))
colnames(AATAMS)[THIS]=c("station.name","receiver.ID")

AATAMS$EffectiveStartDate=as.POSIXlt(strptime(as.character(AATAMS$deployment_date),format='%d/%m/%Y'))
AATAMS$EffectiveStartDate=format(AATAMS$EffectiveStartDate, "%d-%m-%y")
AATAMS$EffectiveEndDate=NA


keep.these.vars=match(c("receiver.ID","station.name","latitude","longitude",
                        "EffectiveStartDate","EffectiveEndDate","Project"),colnames(AATAMS))

AATAMS=AATAMS[,keep.these.vars]



#3. EXPORT FILE
write.csv(AATAMS,file="AATAMS_receivers.manipulated_25_09_2012.csv")