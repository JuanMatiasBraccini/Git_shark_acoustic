  #----------------- SOURCE SCRIPT FOR ACOUSTIC TAGGING ANALYSIS ----------------#

#note: This script is used for any analysis (FRDC, white shark, etc) of acoustic tagging data for
#       loading the data from the different data bases


#note: data must be updated as receivers are downloaded

library(RODBC)  			#library for importing excel data
library(gridBase)			#for inset map
library(PBSmapping)   #for polygon
#library(ade4)				#to add subplots
library(plotrix)				#needed for graph legends
library(MASS)         #for fitting distributions
library(geosphere)   #for spatial statistics and trigonometry in a sphere
#library(VTrack)       #for handling acoustic data
library(mixtools)    #for fitting bimodal normal distributions
library(chron)      #for times
#library(epitools)   #for extracting hours from times object


#library(prada)   #for fitting bivariate normal distribution
# To install "prada"
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("prada")

# library(adehabitatLT)   #for different types of random walks
# library(VGAM)           #for fitting a Levy distribution


#source bubble plot functions
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Bubble.plot.detections.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Bubble.plot.R") 


#1. load data

#1.1 TRANSMITTERS
#note: make sure to have the most updated list of transmitters
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Source_acoustic_transmitters.R")


#1.2 RECEIVER DETECTIONS
setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data")

#ACA: improve this, read in everything in file!!!
  #1.2.1 AATAMS
AATAMS.depth=read.csv("NRETAstations.csv",stringsAsFactors =F)

filenames <- list.files("AATAMS_downloads_processed", pattern="*.csv", full.names=TRUE)
AATAMS <- lapply(filenames, read.csv,stringsAsFactors =F)
  
  #extract relevant columns
AATAMS[[1]]=AATAMS[[1]][,-match(c("sensor.value","sensor.unit"),names(AATAMS[[1]]))]
for(a in 2:length(AATAMS))
{
  id=match(names(AATAMS[[1]]),names(AATAMS[[a]]))
  idd=which(is.na(id))
  if(length(idd>0))
  {
    dummy=as.data.frame(matrix(nrow=nrow(AATAMS[[a]]),ncol=length(names(AATAMS[[1]])[idd])))
    colnames(dummy)=names(AATAMS[[1]])[idd]
    AATAMS[[a]]=cbind(AATAMS[[a]],dummy)
    id=match(names(AATAMS[[1]]),names(AATAMS[[a]]))
  }
  AATAMS[[a]]=AATAMS[[a]][,id]
}
AATAMS=do.call(rbind,AATAMS)


  #1.2.2 SMN
SMN.gummy=read.csv(file="SMN_donwloads_processed/gummy.2017-06-16.csv",stringsAsFactors =F)
SMN.whiskery=read.csv(file="SMN_donwloads_processed/whiskery.2017-06-16.csv",stringsAsFactors =F)
SMN.dusky=read.csv(file="SMN_donwloads_processed/dusky.2017-06-16.csv",stringsAsFactors =F)
SMN.sandbar=read.csv(file="SMN_donwloads_processed/sandbar.2017-06-16.csv",stringsAsFactors =F)

SMN.copper=read.csv(file="SMN_donwloads_processed/copper.2015-08-03.csv",stringsAsFactors =F)


if(use.all=="YES")
{
  SMN.white=read.csv(file="SMN_donwloads_processed/white.2017-08-14.csv",stringsAsFactors =F)
  SMN.tiger=read.csv(file="SMN_donwloads_processed/tiger.2015-08-03.csv",stringsAsFactors =F)
  SMN.greynurse=read.csv(file="SMN_donwloads_processed/greynurse.2015-08-03.csv",stringsAsFactors =F)
}

#combine all species
SMN=rbind(SMN.gummy,SMN.whiskery,SMN.dusky,SMN.sandbar,SMN.copper)
if(use.all=="YES") SMN=rbind(SMN,SMN.white,SMN.tiger,SMN.greynurse)
SMN$Sex=as.character(SMN$Sex)
SMN$Sex=ifelse(SMN$Sex=="m","M",ifelse(SMN$Sex=="?","U",SMN$Sex))


#1.3 ALL RECEIVER LOCATIONS
  #1.3.1 AATAMS
AATAMS.all=read.csv(file="AATAMS_receiver_location/AATAMS_receivers.manipulated_25_09_2012.csv")

  #1.3.2 SMN
SMN.all=read.csv(file="SMN_receiver_location/SMN_receivers.manipulated.csv")
SMN.depths=read.csv("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/Station depths and history.csv")

#drop Two Peoples bay and NRETA
SMN.all=subset(SMN.all,!station.name%in%c("SMN 2 Peoples Bay 01","SMN 2 Peoples Bay 02","SMN 2 Peoples Bay 03",       
 "SMN 2 Peoples Bay 04","SMN 2 Peoples Bay 05","SMN 2 Peoples Bay 06"))       
SMN.all=subset(SMN.all,latitude<(-25))




# #1.5 BATHYMETRY
# Bathymetry_120=read.table("C:/Matias/Data/Mapping/get_data112_120.cgi")
# Bathymetry_138=read.table("C:/Matias/Data/Mapping/get_data120.05_138.cgi")
# Bathymetry=rbind(Bathymetry_120,Bathymetry_138)


# #1.6. shapefile Perth islands
# PerthIs=read.table("C:/Matias/Data/Mapping/WAislandsPointsNew.txt", header=T)
# Rottnest.Is=subset(PerthIs,ID%in%c("ROTT1"))
# Garden.Is=subset(PerthIs,ID%in%c("ROTT3"))



#2. Manipulate data

  #2.1 remove duplicates (AATAMS records from SMN records & overlapping periods)
SMN$Unico=with(SMN,paste(TagCode,Latitude,Longitude,Date.local,Time.local))
duplis=SMN[duplicated(SMN$Unico),]
SMN=SMN[!duplicated(SMN$Unico),]

#deal with tags used multiple times
SMN$ReleaseDate=as.POSIXlt(as.character(SMN$ReleaseDate),format="%d-%b-%y")
AATAMS$Release.Date=as.POSIXlt(as.character(AATAMS$Release.Date))


TAGS$ReleaseDate2=as.POSIXlt(as.character(TAGS$ReleaseDate2),format='%d-%b-%y')
TAGS$ReleaseLatitude2=with(TAGS,ifelse(ReleaseLatitude2>0,-ReleaseLatitude2,ReleaseLatitude2))

a=subset(TAGS,select=c(Code2,Species2,Sex2,ReleaseDate2,ReleaseLatitude2,ReleaseLongitude2,ReleaseLength))
names(a)=paste(names(a),"TAG",sep='_')
b=a[duplicated(a$Code2_TAG),]

Multi.used.tags=unique(b$Code2_TAG)

AATAMS.multi=subset(AATAMS,ID%in%Multi.used.tags)
if(nrow(AATAMS.multi)>0)
{
  AATAMS=subset(AATAMS,!ID%in%Multi.used.tags)
  AATAMS.multi.tgs=unique(AATAMS.multi$ID)
  STR=vector('list',length(AATAMS.multi.tgs))
  for(m in 1:length(AATAMS.multi.tgs))
  {
    d=subset(AATAMS.multi,ID==AATAMS.multi.tgs[m])
    d.tag=subset(a,Code2_TAG==AATAMS.multi.tgs[m])
    d.tag$date.to=c(d.tag$ReleaseDate2_TAG[2:nrow(d.tag)],Sys.Date())
    Str=vector('list',nrow(d.tag))
    for(n in 1:nrow(d.tag))
    {
      x=d.tag[n,]
      xx=subset(d,Date.local>=x$ReleaseDate2_TAG & Date.local<x$date.to)
      if(nrow(xx)>0)
      {
        xx$Species2=x$Species2_TAG
        xx$Sex2=x$Sex2_TA
        xx$Release.Date=x$ReleaseDate2_TAG
        xx$ReleaseLatitude2=x$ReleaseLatitude2_TAG
        xx$ReleaseLongitude2=x$ReleaseLongitude2_TAG
        xx$ReleaseLength=x$ReleaseLength_TAG
        Str[[n]]=xx
      }
      
    }
    STR[[m]]=do.call(rbind,Str)
    
  }
  STR=do.call(rbind,STR)
  AATAMS=rbind(AATAMS,STR)
}



SMN.multi=subset(SMN,TagCode%in%Multi.used.tags)
if(nrow(SMN.multi)>0)
{
  SMN=subset(SMN,!TagCode%in%Multi.used.tags)
  SMN.multi.tgs=unique(SMN.multi$TagCode)
  STR=vector('list',length(SMN.multi.tgs))
  for(m in 1:length(SMN.multi.tgs))
  {
    d=subset(SMN.multi,TagCode==SMN.multi.tgs[m])
    d.tag=subset(a,Code2_TAG==SMN.multi.tgs[m])
    d.tag$date.to=c(d.tag$ReleaseDate2_TAG[2:nrow(d.tag)],Sys.Date())
    Str=vector('list',nrow(d.tag))
    for(n in 1:nrow(d.tag))
    {
      x=d.tag[n,]
      xx=subset(d,Date.local>=x$ReleaseDate2_TAG & Date.local<x$date.to)
      if(nrow(xx)>0)
      {
        xx$Species=x$Species2_TAG
        xx$Sex=x$Sex2_TA
        xx$ReleaseDate=x$ReleaseDate2_TAG
        xx$ReleaseLatitude=x$ReleaseLatitude2_TAG
        xx$ReleaseLongitude=x$ReleaseLongitude2_TAG
        Str[[n]]=xx
      }
      
    }
    STR[[m]]=do.call(rbind,Str)
  }
  STR=do.call(rbind,STR)
  SMN=rbind(SMN,STR)
}


#from AATAMS keep only species in SMN
AATAMS=subset(AATAMS,Species2%in%unique(SMN$Species)) 

SMN$Unico=with(SMN,paste(TagCode,round(Latitude,0),round(Longitude,0),Date.local))
AATAMS$Unico=with(AATAMS,paste(ID,round(latitude,0),round(longitude,0),Date.local))

#remove SMN records already available in AATAMS
AATAMS.in.SMN=subset(SMN,Unico%in%AATAMS$Unico)
SMN=subset(SMN,!Unico%in%AATAMS$Unico)

# AATAMS.in.SMN=subset(AATAMS,Unico%in%SMN$Unico)
# AATAMS=subset(AATAMS,!Unico%in%SMN$Unico)

AATAMS$Unico=with(AATAMS,paste(ID,latitude,longitude,Date.local,Time.local))
duplis=AATAMS[duplicated(AATAMS$Unico),]
AATAMS=AATAMS[!duplicated(AATAMS$Unico),]


SMN=SMN[,-match("Unico",names(SMN))]
if(match("StationName2",names(SMN))>0)SMN=SMN[,-match("StationName2",names(SMN))]


  #2.2. tidy up detection datasets
AATAMS$Name=NA
these.AATAMS.vars=c("ID","Species2","Sex2","receiver.ID","Release.Date","ReleaseLatitude2",
                    "ReleaseLongitude2","latitude","longitude","DateTime.local","Date.local","Time.local","Name")
matched.AATAMS.vars=match(these.AATAMS.vars,names(AATAMS))
AATAMS=AATAMS[,matched.AATAMS.vars]

not.this.SMN=match(c("X"),names(SMN))
not.this.SMN=not.this.SMN[!is.na(not.this.SMN)]
if(length(not.this.SMN)>0)SMN=SMN[,-not.this.SMN]


SMN$Project="SMN"
AATAMS$Project="AATAMS"
AATAMS$Depth=NA

re.arranged.cols=match(c("ID","Species2","Sex2","receiver.ID","Release.Date","ReleaseLatitude2","ReleaseLongitude2",
                         "latitude","longitude","Depth","DateTime.local","Date.local","Time.local","Project"),
                       names(AATAMS))
AATAMS=AATAMS[,re.arranged.cols]

colnames(AATAMS)=colnames(SMN)

AATAMS$SerialNumber=abs(as.numeric(gsub("^.*?-","-",as.character(AATAMS$SerialNumber))))



#Interpolate receiver's depth if missing
SMN.depths$Array=with(SMN.depths,ifelse(Latitude1<(-31)& Latitude1>(-32.09),"Perth",
            ifelse(Latitude1<(-32.09)& Latitude1>(-32.5),"Cock.Sound",
            ifelse(Latitude1<(-34)& Latitude1>(-34.5) & Longitude1<116,"SL1",
            ifelse(Latitude1<(-34.5) & Longitude1<117,"SL2",
            ifelse(Latitude1<(-34.5) & Longitude1>117,"SL3",NA))))))

id=which(is.na(SMN$Depth))
SMN.missing.z=SMN[id,]
missing.z=unique(SMN.missing.z$SerialNumber)

#see hits depth and missing depths
fn.plot.depth=function(a,factor,pt)
{
  plot(a$Longitude1,a$Latitude1,pch=19)
  points(pt$Longitude,pt$Latitude,pch=19,col=2)
  text(a$Longitude1*factor,a$Latitude1*1,a$Depth,cex=0.85)
  text(pt$Longitude*factor*0.9999,pt$Latitude,unique(pt$SerialNumber),cex=0.85,col=2)
}
fn.plot.depth(subset(SMN.depths,Array=="SL1"),1.00025,subset(SMN.missing.z,SerialNumber==119459))
fn.plot.depth(subset(SMN.depths,Array=="SL2"),1.000025,subset(SMN.missing.z,SerialNumber%in%c(119505,119508)))

SMN$Depth=with(SMN,ifelse(is.na(Depth) & SerialNumber==119459,140,
      ifelse(is.na(Depth) & SerialNumber%in%c(119505,119508),50,
      ifelse(is.na(Depth) & SerialNumber==105740,14.2,Depth))))
  
        
#add AATAMS depth
AATAMS.depth$Receivers=gsub("c", "", AATAMS.depth$Receivers)
AATAMS.depth$Receivers=with(AATAMS.depth,
               ifelse(substring(Receivers, first=1, last = 3)=="VR3",
               as.numeric(substring(Receivers, first=9, last = 15)),Receivers))

AATAMS.depth$Receivers=abs(as.numeric(gsub("^.*?-","-",AATAMS.depth$Receivers)))
names(AATAMS.depth)[match("Receivers",names(AATAMS.depth))]="SerialNumber"

AATAMS.depth$Start.date=as.POSIXlt(as.character(AATAMS.depth$Start.date),format='%d/%m/%Y')
Today=as.character(Sys.Date())
AATAMS.depth$End.Date=as.POSIXlt(as.character(AATAMS.depth$End.Date),format='%d/%m/%Y')
AATAMS.depth$End.Date=as.character(AATAMS.depth$End.Date)
AATAMS.depth$End.Date=with(AATAMS.depth,ifelse(is.na(End.Date),Today,End.Date))
AATAMS.depth$End.Date=as.POSIXlt(as.character(AATAMS.depth$End.Date))
AATAMS$Date.local=as.POSIXlt(as.character(AATAMS$Date.local))
AATAMS.depth$Latitude=with(AATAMS.depth,ifelse(Latitude>0,-Latitude,Latitude))

b=AATAMS[!duplicated(paste(AATAMS$SerialNumber,AATAMS$Date.local)),]
b=merge(b,AATAMS.depth[,match(c("SerialNumber","Start.date","End.Date","Depth"),names(AATAMS.depth))],
        by="SerialNumber",all=T)  
b$keep=with(b,ifelse(Date.local>=Start.date & Date.local<=End.Date,"Yes","No"))
b=subset(b,keep=="Yes")
b=b[,match(c("SerialNumber","Date.local","Depth.y"),names(b))]
AATAMS=merge(AATAMS,b,by=c("SerialNumber","Date.local"),all.x=T)
AATAMS$Depth=AATAMS$Depth.y

AATAMS$Unico=with(AATAMS,paste(TagCode,Latitude,Longitude,Date.local,Time.local))
duplis=AATAMS[duplicated(AATAMS$Unico),]
AATAMS=AATAMS[!duplicated(AATAMS$Unico),]
AATAMS=AATAMS[,match(names(SMN),names(AATAMS))]


#no depth for these receivers-dates
b=subset(AATAMS,is.na(Depth) & Longitude<118)
b$Unic=paste(b$SerialNumber,b$Date.local)
b=b[!duplicated(b$Unic),]
plot(b$Longitude,b$Latitude)
b$Date.local=as.character(b$Date.local)
b=b[order(b$SerialNumber,b$Date.local),c(4,12)]
Recs=unique(b$SerialNumber)

listas=vector("list",length(Recs))
names(listas)=Recs

for(i in 1:length(listas))
{
  x=subset(b,SerialNumber==Recs[i])
  listas[[i]]=c(Recs[i],x$Date.local[1],x$Date.local[nrow(x)])
}
Check=as.data.frame(do.call(rbind,listas))
names(Check)=c("SerialNumber","First.hit","Last.hit")

listas=vector("list",nrow(Check))
for (i in 1:nrow(Check))
{
  a=subset(AATAMS.depth,SerialNumber==Check[i,1],select=c(SerialNumber,Start.date,End.Date))
  listas[[i]]=merge(Check[i,],a,by="SerialNumber")
  
}
Recs=as.numeric(as.character(Check$SerialNumber))


#merge data files
SMN$Date.local=as.character(SMN$Date.local)
SMN$Time.local=as.character(SMN$Time.local)
SMN$Species=as.character(SMN$Species)

AATAMS$Date.local=as.character(AATAMS$Date.local)
AATAMS$Time.local=as.character(AATAMS$Time.local)
AATAMS$Species=as.character(AATAMS$Species)

Detections=rbind(SMN,AATAMS)
rm(list=c("AATAMS","SMN"))


#fix some sex info errors
Detections$Sex=with(Detections,ifelse(Sex %in% c("FALSE","?") & TagCode %in% c(31012,30934,30931,29599,
                        29523,29507,29474,29473,29470,29587),"F",
                        ifelse(Sex %in% c("FALSE","?") & TagCode %in% c(29566),"M",
                               
                               Sex)))
Detections$Sex=ifelse(Detections$Sex=="m","M",ifelse(Detections$Sex=="f","F",Detections$Sex))

Detections$Sex=ifelse(Detections$TagCode==29566,"M",Detections$Sex)



#Export inputs for analyses
setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data")
write.csv(Detections,"Detections.csv",row.names=F)
write.csv(AATAMS.all,"AATAMS.all.csv",row.names=F)
write.csv(SMN.all,"SMN.all.csv",row.names=F)
write.csv(TAGS,"TAGS.csv",row.names=F)
