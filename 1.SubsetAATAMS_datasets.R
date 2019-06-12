########################- SCRIPT FOR SUBSETTING ATTAMS DATASETS  -##########################

#notes: This script subset the various data sets downloaded from the AATAMS website
#       for creating smaller datasets for species of interest

library(chron) #for extracting time
library(lubridate)

setwd("C:/Users/myb/Desktop/AATAMS_donwloads")

#note: keep an eye on timestamp format when processing new data!!!!



  #Select batch of data to process
#Current.batch="Batch 21.09.2012"
#Current.batch="Batch 07.01.2014"
#Current.batch="Batch 05.08.2014"
#Current.batch="Batch 10.04.2015"
#Current.batch="Batch 1.10.2015"
#Current.batch="Batch 7.12.2016"
#Current.batch="Batch 22.2.2017"
Current.batch="Batch 1.11.2018"

#1. RECEIVER DATA

#Load receiver data

    #Batch 21.09.2012
if(Current.batch=="Batch 21.09.2012")
{
  A_Glen=read.csv(file="AATAMS_Glenelg_01_01_2007_to_21_09_2012.csv",stringsAsFactors =F)
  A_Row_S=read.csv(file="AATAMS_Rowley_S_01_01_2007_to_21_09_2012.csv",stringsAsFactors =F)
  A_Scott_R=read.csv(file="AATAMS_Scott_R_01_01_2007_to_21_09_2012.csv",stringsAsFactors =F)
  CSIRO_Shark=read.csv(file="CSIRO_Shark_01_01_2007_to_21_09_2012.csv",stringsAsFactors =F)
  G_St_Vin=read.csv(file="Gulf_St_Vincent_01_01_2007_to_21_09_2012.csv",stringsAsFactors =F)
  NRETA.1=read.csv(file="NRETA_01_01_2007_to_31_12_2009.csv",stringsAsFactors =F)
  NRETA.2=read.csv(file="NRETA_01_01_2010_to_31_07_2011.csv",stringsAsFactors =F)
  NRETA.3=read.csv(file="NRETA_01_08_2011_to_01_12_2011.csv",stringsAsFactors =F)
  NRETA.4=read.csv(file="NRETA_02_12_2011_to_21_05_2012.csv",stringsAsFactors =F)
  NRETA.5=read.csv(file="NRETA_22_05_2012_to_21_09_2012.csv",stringsAsFactors =F)
  
  List=list(A_Glen=A_Glen,A_Row_S=A_Row_S,A_Scott_R=A_Scott_R,CSIRO_Shark=CSIRO_Shark,G_St_Vin=G_St_Vin,
            NRETA.1=NRETA.1,NRETA.2=NRETA.2,NRETA.3=NRETA.3,NRETA.4=NRETA.4,NRETA.5=NRETA.5)
}

    #Batch 07.01.2014
if(Current.batch=="Batch 07.01.2014")
{
  Rowley_S=read.csv(file="Rowley_Shoals_R_01_01_2007_to_07_01_2014.csv",stringsAsFactors =F)
  Rowley_S$timestamp=as.POSIXlt(as.character(Rowley_S$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  CSIRO_Shark=read.csv(file="detection_CSIRO.animalTag. 8.12.2009.to. 22.04.2012.csv",stringsAsFactors =F)
  G_St_Vin.1=read.csv(file="detection_GSV. 6.8.2010.to6.18.2012.csv",stringsAsFactors =F)
  G_St_Vin.2=read.csv(file="detection_GSV. 6.8.2012.to17.5.2013.csv",stringsAsFactors =F)
  NRETA.1=read.csv(file="detection_NRETA.CB. 22.09.12.to. 18.08.13.csv",stringsAsFactors =F)
  NRETA.1$timestamp=as.POSIXlt(as.character(NRETA.1$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT") 
  NRETA.2=read.csv(file="detection_NRETA.MB. 22.09.12.to. 31.05.13.csv",stringsAsFactors =F)
  NRETA.3=read.csv(file="detection_NRETA.N.C.S.lines. 22.09.12.to. 17.08.13.csv",stringsAsFactors =F)
  A_Glen=read.csv(file="detection_AATAMS.Glenelg. 30.04.2011.to11.11.2012.csv",stringsAsFactors =F)
  A_Row_S=read.csv(file="detection_AATAMS.RS. 21.09.2012.to10.04.2013.csv",stringsAsFactors =F)
  A_Scott_R=read.csv(file="detection_AATAMS.SR. 21.09.2012.to21.11.2012.csv",stringsAsFactors =F)
  
  List=list(Rowley_S=Rowley_S,CSIRO_Shark=CSIRO_Shark,G_St_Vin.1=G_St_Vin.1,G_St_Vin.2=G_St_Vin.2,
            NRETA.1=NRETA.1,NRETA.2=NRETA.2,A_Glen=A_Glen,A_Row_S=A_Row_S,A_Scott_R=A_Scott_R,NRETA.3=NRETA.3)
  
}

    #Batch 05.08.2014
if(Current.batch=="Batch 05.08.2014")
{
  NRETA=read.csv(file="detection_NRETA. 1.09.12.to. 30.09.12.csv",stringsAsFactors =F)
  NRETA$timestamp=as.POSIXlt(as.character(NRETA$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  List=list(NRETA=NRETA)
}

  #Batch 10.04.2015
if(Current.batch=="Batch 10.04.2015")
{
  G_St_Vin=read.csv(file="Gulf_St_Vincent_01_01_2013_to_24_06_2014.csv",stringsAsFactors =F)
  G_St_Vin$timestamp=as.POSIXlt(as.character(G_St_Vin$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT") 
  
  NRETA=read.csv(file="NRETA_1_1_2014_to_13_12_2014.csv",stringsAsFactors =F)
  NRETA$timestamp=as.POSIXlt(as.character(NRETA$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  
  NRETA1=read.csv(file="NRETA_31_5_2013_to_1_1_2014.csv",stringsAsFactors =F)
  NRETA1$timestamp=as.POSIXlt(as.character(NRETA1$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  
  A_Glen=read.csv(file="AATAMS_Glenelg_20_04_2011_to_26_06_2014.csv",stringsAsFactors =F)
  A_Glen$timestamp=as.POSIXlt(as.character(A_Glen$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  
  A_Scott_R=read.csv(file="AATAMS_Scott_R_01_01_2012_to_26_02_2014.csv",stringsAsFactors =F)
  A_Scott_R$timestamp=as.POSIXlt(as.character(A_Scott_R$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  
  List=list(G_St_Vin=G_St_Vin,NRETA1=NRETA1,NRETA=NRETA,A_Glen=A_Glen,A_Scott_R=A_Scott_R)
  
}

#Batch 1.10.2015
if(Current.batch=="Batch 1.10.2015")
{
  A_Glen=read.csv(file="AATAMS_Glenelg_9_10_2014_to_08_08_2015.csv",stringsAsFactors =F)
  Location_Glen=read.csv("Glenelg Receiver lat and long.csv",stringsAsFactors =F)
  names(Location_Glen)[2]="Receiver"
  A_Glen$timestamp=as.POSIXlt(as.character(A_Glen$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  A_Glen$Receiver=as.numeric(substring(as.character(A_Glen$Receiver), first=6, last = 30))
  A_Glen=merge(A_Glen,Location_Glen[,-1],by="Receiver",all.x=T)
  names(A_Glen)[c(1,3)]=c("receiver.ID","transmitter.ID")
  List=list(A_Glen=A_Glen)
}


#Batch 7.12.2016
if(Current.batch=="Batch 7.12.2016")
{
  NRETA=read.csv(file="NRETA_14_12_2014_to_21_3_2016_.csv",stringsAsFactors =F)
  NRETA$timestamp=as.POSIXlt(as.character(NRETA$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  
  A_Glen=read.csv(file="AATAMS_Glenelg_09_08_2015_to_4_4_2016.csv",stringsAsFactors =F)
  A_Glen$timestamp=as.POSIXlt(as.character(A_Glen$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  
  A_Scott_R=read.csv(file="AATAMS_Scott_R_27_02_2014_to_6_09_2016.csv",stringsAsFactors =F)
  A_Scott_R$timestamp=as.POSIXlt(as.character(A_Scott_R$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  
  List=list(NRETA=NRETA,A_Glen=A_Glen,A_Scott_R=A_Scott_R)
  
}

if(Current.batch=="Batch 22.2.2017")
{
  NRETA=read.csv(file="NRETA_22_03_2016_to_11_12_2016.csv",stringsAsFactors =F)
  NRETA$timestamp=as.POSIXlt(as.character(NRETA$timestamp),  format="%d/%m/%Y %H:%M",tz="GMT")
  
  List=list(NRETA=NRETA)
  
}

if(Current.batch=="Batch 1.11.2018")
{
  NRETA=read.csv(file="IMOS-ATF_NRETA_detections_2007-2018.csv",stringsAsFactors =F)
  NRETA$timestamp=as_datetime(NRETA$timestamp)
  NRETA=subset(NRETA,timestamp>as_datetime("2016-12-11"))
  List=list(NRETA=NRETA)
  
}

#Convert factors to character
fn.fact.to.char=function(a)
{
  a$timestamp=as.character(a$timestamp)
  a$receiver.ID=as.character(a$receiver.ID)
  a$transmitter.ID=as.character(a$transmitter.ID)  
  return(a)
}
for(s in 1:length(List))List[[s]]=fn.fact.to.char(List[[s]])


#Combine receiver data
Data=do.call(rbind,List)
rm(List)


#Extract tag id
Data$ID=as.numeric(substring(as.character(Data$transmitter.ID), first=10, last = 30))


#Add tag information
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Source_acoustic_transmitters.R")  #Load tag ids. note: keep transmiters updated
cool.tags=as.numeric(unique(TAGS$Code2))


# Keep only tag IDs of interest
Data=subset(Data,ID%in%cool.tags)


#Extract Date and Time
Data$DateTime=as_datetime(as.character(Data$timestamp))


#4. Merge Data and TAGS
reused=table(TAGS$Code2)
reused=reused[reused>1]
TAGS.reused=subset(TAGS,Code2%in%names(reused))   #remove reused tags 
TAGS=subset(TAGS,!Code2%in%names(reused))

Detec.rango=range(Data$timestamp)

Uni.reused=unique(TAGS.reused$Code2)

dummy=NULL
for(u in 1:length(Uni.reused))
{
  a=subset(TAGS.reused,Code2==Uni.reused[u])
  within.obs.range=subset(a,as.Date(ReleaseDate2)>min(as.Date((Detec.rango))))
  if(nrow(within.obs.range)==0) within.obs.range=subset(a,ReleaseDate2==max(a$ReleaseDate2))
  dummy=rbind(dummy,within.obs.range)
}
TAGS=rbind(TAGS,dummy)

Data=merge(Data,TAGS,by.x="ID",by.y="Code2",all.x=T)



#6. Convert UTC detection time to local time in WA and SA and add date when needed
#note: daylight savings in SA was not consider as this is an artificial manipulation of time
#       must split Data by state because ifelse does not allow posixt!

UTC.WA=8
UTC.SA=9.5
UTC.SA.savn=10.5 
List.SA.savn=list(c("2007-10-28","2008-04-06"),c("2008-10-05","2009-04-05"),c("2009-10-04","2010-04-04"),
                  c("2010-10-10","2011-04-03"),c("2011-10-02","2012-04-01"),c("2011-10-07","2013-04-07"),
                  c("2013-10-06","2014-04-06"),c("2014-10-05","2015-04-05"))

Data.WA=subset(Data,longitude<=129)
Data.SA=subset(Data,longitude>129)

Data.WA$DateTime.local=Data.WA$DateTime+(UTC.WA*3600)
Data.SA$DateTime.local=Data.SA$DateTime+(UTC.SA*3600)

Data=rbind(Data.WA,Data.SA)
                  

#7. Separate date and time
#Data$Date.local=trunc(Data$DateTime.local,"days")
Data$Date.local=as.Date(as.POSIXlt(as.character(Data$DateTime.local)))
Data$Time.local=times(format(Data$DateTime.local, "%H:%M:%S"))
Data$Release.Date=Data$ReleaseDate2



#8. Some manipulations
these.cols=match(c("species","tag.ID","Name2","timestamp","DateTime","ReleaseDate2"),colnames(Data))#remove these columns
Data=Data[,-these.cols]



#9. flag white sharks tagged by shark dive operators
Data$Neptune.Chart=with(Data,ifelse(Species2=="White" & ReleaseLongitude2>=135 & ReleaseLongitude2<138 &
  Release.Date>'2011-07-01',"YES","NO"))


#10. Species manipulation 
Data$Species2=ifelse(Data$Species2=="OTN SENTINEL","SENTINEL",Data$Species2)




#8. Export data
setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/AATAMS_downloads_processed")

if(Current.batch=="Batch 21.09.2012") write.csv(Data,file="AATAMS.data.upto.21_09_2012.csv")
if(Current.batch=="Batch 07.01.2014") write.csv(Data,file="AATAMS.data.22_09_2012.to.7.1.2014.csv")
if(Current.batch=="Batch 05.08.2014") write.csv(Data,file="AATAMS.data.01_09_2012.to.30_09_2012.csv")
if(Current.batch=="Batch 10.04.2015") write.csv(Data,file="AATAMS.data.to.10_04_2015.csv")
if(Current.batch=="Batch 1.10.2015") write.csv(Data,file="Glenelg.data.to.08_08_2015.csv")
if(Current.batch=="Batch 7.12.2016") write.csv(Data,file="AATAMS.data.to.21_03_2016.csv")
if(Current.batch=="Batch 22.2.2017") write.csv(Data,file="AATAMS.data.to.11_12_2016.csv")
if(Current.batch=="Batch 1.11.2018") write.csv(Data,file="AATAMS.data.to.02_10_2018.csv")




