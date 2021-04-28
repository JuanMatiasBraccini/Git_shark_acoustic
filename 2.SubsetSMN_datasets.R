
########################- SCRIPT FOR SUBSETTING SMN DATASETS  -##########################

#notes: This script subset the data set downloaded from the SMN website
#steps:
        #1. For each species, create 'Event Details Report' in SMN/Reports 
#           from 12/04/2008 to date. Save in desktop (see setwd) using the 
#           same names specified below

library(chron) #for extracting time
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')
setwd("C:/Users/myb/Desktop/SMN_downloads")


#1. Load detection data

#Up to 19/10/2018
#   #select species
Data.copper=read.csv(file="EventDetailsReport_copper.csv",stringsAsFactors = F)
Data.dusky=read.csv(file="EventDetailsReport_dusky.csv",stringsAsFactors = F)
Data.greynurse=read.csv(file="EventDetailsReport_grey.nurse.csv",stringsAsFactors = F)
Data.gummy=read.csv(file="EventDetailsReport_gummy.csv",stringsAsFactors = F)
Data.sandbar=read.csv(file="EventDetailsReport_sandbar.csv",stringsAsFactors = F)
Data.tiger=read.csv(file="EventDetailsReport_tiger.csv",stringsAsFactors = F)
Data.whiskery=read.csv(file="EventDetailsReport_whiskery.csv",stringsAsFactors = F)
Data.white=read.csv(file="EventDetailsReport_white.csv",stringsAsFactors = F) 
DATA=list(copper=Data.copper,
          dusky=Data.dusky,
          gummy=Data.gummy,
          sandbar=Data.sandbar,
          whiskery=Data.whiskery,
          tiger=Data.tiger,
          greynurse=Data.greynurse,
          white=Data.white)

UTC.WA=8

fn.extract=function(Data)
{
  #Extract Date and Time
  #Data$DateTime=as.POSIXlt(as.character(Data$"?..Datetime"), format="%d-%b-%Y %r", tz="GMT") #use this for upto 11/10/2012
  Data$DateTime=as.POSIXlt(as.character(Data$Datetime),  format="%d/%m/%Y %H:%M",tz="GMT")#use this for 12/10/2012 to 27/12/2012
  Data$DateTime.local=Data$DateTime+(UTC.WA*3600)  
  Data$Date.local=trunc(Data$DateTime.local,"days")
  Data$Time.local=times(format(Data$DateTime.local, "%H:%M:%S"))
    
  
  #Some manipulations
  Data$Latitude=with(Data,ifelse(Latitude>0,-Latitude,Latitude))
  
  Data$Species=as.character(Data$Species)
  Data$Species=ifelse(Data$Species%in%c("OTN SENTINEL","Test"),"SENTINEL",Data$Species)
  these.cols=match(c("Datetime","StationCode","StationName","DateTime"),colnames(Data))#remove these columns
  Data=Data[,-these.cols]
  
  #Remove AATAMS data
 # Data=subset(Data,Latitude<(-30))
  
  
  
  return(Data)
}

#Export data
setwd(handl_OneDrive("Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/SMN_donwloads_processed"))
system.time(for(i in 1:length(DATA)) write.csv(fn.extract(DATA[[i]]),
                                    file=paste(names(DATA)[i],".",Sys.Date(),".csv",sep='')))
