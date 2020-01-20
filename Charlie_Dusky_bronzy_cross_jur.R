#--SCRIPT FOR ANALYSING CROSS-JURISDICTIONAL MOVEMENTS OF DUSKY & BRONZE WHALERS--

#note: must report the total number of WA and SA tags and detected to put those that 
#      move across in perspective...

#to do:
    # map release location and trajectories (showing where detected). Can do trajectories with progressive
    #       shades of grey to show initial and most recent...
    # seasonal patterns: GAM 

rm(list=ls(all=TRUE))
library(lubridate)
library(tidyverse)
library(geosphere)
library(chron) 
library(mgcv)
library(mgcViz)
library(vioplot)
library(beanplot)

options(stringsAsFactors = FALSE)

#1. Data Section---------------------------------------------------------
setwd('C:\\Matias\\Analyses\\Acoustic_tagging\\For Charlie\\Data')
Dat=read.csv('All detections_2018 01.csv')
SA.receivers=read.csv('GSV receiver location-IMOS.csv')
Blank_LatLong=read.csv('Blank_LatLong.csv')


#SMN individual tag download Nov 2019
setwd('C:\\Matias\\Analyses\\Acoustic_tagging\\For Charlie\\Data\\SMN_download_Nov_19')
file_list  <- list.files()
DATA.SMN=do.call(rbind,lapply(file_list,read.csv))
setwd('C:\\Matias\\Analyses\\Acoustic_tagging\\For Charlie')


#bring in conventional tagging data
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_conventional_data.R")


#2. Parameters Section---------------------------------------------------

  #define if doing exploratory analysis
do.expl="NO"

  #define arbitrary turning points
Cape.Leuwin=cbind(114.969,-34.459)
Mid.point=cbind(116.425,-35.043)
SA.Border=cbind(129.040964,-31.744653)
Eyre=c(135.694043,-34.883332)
Eyre_SA.Border=distGeo(Eyre,SA.Border)
SA.Border_Mid.point=distGeo(SA.Border,Mid.point)
Mid.point_Cape.Leuwin=distGeo(Mid.point,Cape.Leuwin)

#3. Procedure Section----------------------------------------------------

#3.1 initial data manipulation

#SMN
DATA.SMN$n=grepl("/",DATA.SMN$ï..Datetime)
DATA.SMN$n=with(DATA.SMN,ifelse(n=="TRUE","%d/%m/%Y %H:%M",ifelse(n=="FALSE","%d-%B-%Y %H:%M",NA)))
DATA.SMN$Datetime=as.POSIXlt(DATA.SMN$ï..Datetime,format=DATA.SMN$n,tz="UTC")

DATA.SMN$Date.local=date(DATA.SMN$Datetime)
DATA.SMN$Time.local=times(format(DATA.SMN$Datetime, "%H:%M:%S"))

DATA.SMN$n=grepl("20",DATA.SMN$ReleaseDate)
DATA.SMN$n=with(DATA.SMN,ifelse(n=="FALSE","%d-%B-%y",ifelse(n=="TRUE","%d-%B-%Y",NA)))
DATA.SMN$ReleaseDate=as.POSIXlt(DATA.SMN$ReleaseDate,format=DATA.SMN$n,tz="UTC")
DATA.SMN=DATA.SMN%>%select(-n)  

#Charlie's file  
Dat_no.name=subset(Dat,StationName=="") %>% select(-StationName)
Dat_no.name=left_join(Dat_no.name,Blank_LatLong,by="ReceiverSerialNumber")
Dat=subset(Dat,!StationName=="") %>% select(names(Dat_no.name))
Dat=rbind(Dat,Dat_no.name)
Dat = left_join(Dat,SA.receivers,by=c("StationName" = "Location.name")) %>%
      mutate(Datetime=as.POSIXct(Datetime..WST.,format="%d/%m/%Y %H:%M",tz="UTC"),
             Longitude=ifelse(is.na(Longitude.x),Longitude.y,Longitude.x),
             Latitude=ifelse(is.na(Latitude.x),Latitude.y,Latitude.x),
             Latitude=-abs(Latitude))

#manual addition of new shark
add.dumy=Dat[1:3,] 
add.dumy[,]=NA
add.dumy=add.dumy%>%mutate(TagCode=27698,
                           Species="bronze whaler",
                           State='SA',
                           Organisation='WA Fisheries',
                           Latitude.y=-35.23,
                           Longitude.y=136.07,
                           Datetime=as.POSIXct(c("2019-06-09 12:37:00",
                                                 "2019-06-09 13:37:00",
                                                 "2019-06-09 14:37:00")),
                           Month=6,
                           Year=2019,
                           Longitude=Longitude.y,
                           Latitude=Latitude.y,
                           Sex="M",
                           Station.type="IMOS")
Dat=rbind(Dat,add.dumy)


#Combine Charlie's file and SMN
  #remove duplicates
Dat$Unico=with(Dat,paste(TagCode,Latitude,Longitude,Datetime))
DATA.SMN$ReleaseLatitude=-abs(DATA.SMN$ReleaseLatitude)
DATA.SMN$Latitude=-abs(DATA.SMN$Latitude)
DATA.SMN$Unico=with(DATA.SMN,paste(TagCode,Latitude,Longitude,Datetime))
DATA.SMN$Station.type=with(DATA.SMN,ifelse(grepl("OTN",StationName),"IMOS","DoF"))
duplis=which(DATA.SMN$Unico%in%Dat$Unico)
if(length(duplis)>0)DATA.SMN=DATA.SMN[-duplis,]
dummy=Dat[1:nrow(DATA.SMN),]
dummy[,]=NA             

dummy=dummy%>%mutate(
  TagCode=DATA.SMN$TagCode,
  Organisation=ifelse(DATA.SMN$ReleaseLongitude>129,'Flinders/SARDI','WA Fisheries'),      
  State=ifelse(DATA.SMN$Longitude>129,'SA','WA'),
  Species=DATA.SMN$Species, 
  Longitude=DATA.SMN$Longitude,
  Latitude=DATA.SMN$Latitude,
  Datetime=as.POSIXct(DATA.SMN$Datetime),
  Unico=DATA.SMN$Unico)

Dat$Station.type=with(Dat,ifelse(grepl("OTN",StationName),"IMOS",
                          ifelse(is.na(Station.type) & Longitude<129,"DoF",Station.type)))
Dat=rbind(Dat,dummy) 

#Add release date and location as first location
DATA.SMN$ReleaseDate=paste(DATA.SMN$ReleaseDate,"00:01")
Rel.dat=DATA.SMN%>%
              distinct(TagCode,.keep_all = TRUE)%>%
              select(TagCode,ReleaseDate,ReleaseLatitude,ReleaseLongitude)%>%
              mutate(ReleaseState=ifelse(ReleaseLongitude>129,'SA','WA'))
                   
#add missing SA release data                          
#from Charlie:
# 33194:  SA tagged shark, Carcharhinus brachyurus 104 cm TL (86 cm FL) Female Tagged 22/02/2013

Rel.extra=data.frame(TagCode=c(30992, 23303),
                       ReleaseDate=c("2012-04-26 00:01","2015-02-14 00:01"),
                    ReleaseLatitude=c(-35.1058,-34.91),
                    ReleaseLongitude=c(117.85105,136.22)) %>%
                    mutate(ReleaseState=ifelse(ReleaseLongitude>129,'SA','WA'))
Rel.dat=rbind(Rel.dat,Rel.extra)%>%
  mutate(Relzone=ifelse(ReleaseLatitude>(-33) & ReleaseLatitude<=(-26) & ReleaseLongitude<116.5,"WC",
              ifelse(ReleaseLatitude<=(-33) & ReleaseLongitude<116.5,"Zone1",
              ifelse(ReleaseLongitude>=116.5 & ReleaseLongitude<129 & ReleaseLatitude<=(-26),"Zone2",
              ifelse(ReleaseLongitude>=129 & ReleaseLongitude<135.7,"SA.west",
              ifelse(ReleaseLongitude>=135.6,"SA.east",NA))))))

Dat=Dat%>%
  left_join(Rel.dat,by="TagCode")%>%
  arrange(TagCode, Datetime) %>%
  mutate(TagCode.prev=lag(TagCode,1),
         Datetime.prev=lag(Datetime,1),
         Longitude.prev=lag(Longitude,1),
         Latitude.prev=lag(Latitude,1),
         State.prev=lag(State,1),
         Datetime.prev=as.POSIXct(ifelse(is.na(TagCode.prev)|!TagCode==TagCode.prev,
                              as.character(ReleaseDate),as.character(Datetime.prev)),tz = "UTC"),
         Longitude.prev=ifelse(is.na(TagCode.prev)|!TagCode==TagCode.prev,
                               ReleaseLongitude,Longitude.prev),
         Latitude.prev=ifelse(is.na(TagCode.prev)|!TagCode==TagCode.prev,ReleaseLatitude,Latitude.prev),
         State.prev=ifelse(is.na(TagCode.prev)|!TagCode==TagCode.prev,ReleaseState,State.prev),
         TagCode.prev=ifelse(is.na(TagCode.prev),TagCode,TagCode.prev))

#Some manipulations
Dat = Dat %>%
        mutate(Mn=month(Datetime),
              Yr=year(Datetime),
              N=1,
              zone=ifelse(Latitude>(-33) & Latitude<=(-26) & Longitude<116.5,"WC",
                   ifelse(Latitude<=(-33) & Longitude<116.5,"Zone1",
                   ifelse(Longitude>=116.5 & Longitude<129 & Latitude<=(-26),"Zone2",
                   ifelse(Longitude>=129 & Longitude<135.7,"SA.west",
                   ifelse(Longitude>=135.6,"SA.east",NA))))),
              zone.prev=ifelse(Latitude.prev>(-33) & Latitude.prev<=(-26) & Longitude.prev<116.5,"WC",
                        ifelse(Latitude.prev<=(-33) & Longitude.prev<116.5,"Zone1",
                        ifelse(Longitude.prev>=116.5 & Longitude.prev<129 & Latitude.prev<=(-26),"Zone2",
                        ifelse(Longitude.prev>=129 & Longitude.prev<135.7,"SA.west",
                        ifelse(Longitude.prev>=135.6,"SA.east",NA))))),
              State.prev=ifelse(Longitude.prev>129,'SA','WA')) %>%
      arrange(TagCode, Datetime) %>%
      mutate(Same.station=ifelse(TagCode==TagCode.prev & 
                                 Longitude==Longitude.prev &
                                 Latitude==Latitude.prev  ,"YES","NO"),
             Time=ifelse(TagCode==TagCode.prev,difftime(Datetime,Datetime.prev,units="mins"),NA)) %>%   
      rename(Rel.state=ReleaseState)%>%
      select(c(TagCode,TagCode.prev,Species,Organisation,State,State.prev,
               Datetime,Datetime.prev,Mn,Yr,StationName,Longitude,Longitude.prev,
               Latitude,Latitude.prev,Same.station,zone,zone.prev,Time,N,Rel.state,Station.type))
      
if(sum(is.na(Dat$Latitude))>0) cat("------",sum(is.na(Dat$Latitude)),"RECORDS HAVE NO LATITUDE-----")


#Remove irrelevant TagCodes (whtie sharks)
Dat=subset(Dat,!TagCode%in%c(26438,33194))

#Check number of detections per state
A=as.data.frame(table(Dat$TagCode))
B=as.matrix(table(Dat$TagCode,Dat$State))
colnames(A)=c("TagCode","Total")
A$SA=B[,1]
A$WA=B[,2]


# Report table of tagcodes by species and state----------------------------------------------------------------
setwd('C:\\Matias\\Analyses\\Acoustic_tagging\\For Charlie\\Results')
Tab1= group_by(Dat, TagCode, Species,Organisation,State) %>%
  summarise(sum = sum(N)) %>%
  as.data.frame()
write.csv(Tab1,'Tab.tag.code_species_org_state.csv',row.names = F)



#Run data set scenarios-------------------------------------------------------
#Scen1: full data set
#Scen2: only state
#Scen3: only IMOS
Dat=Dat%>%mutate(Station.type=
                   ifelse(is.na(Station.type) & Longitude<115.626 & 
                            Latitude<(-31) & Latitude>(-32.28),
                          "IMOS",
                   ifelse(is.na(Station.type) & Longitude>116 & Longitude<=129,
                          "DoF",Station.type)))
Dat.scen1=Dat
Dat.scen2=subset(Dat, is.na(Station.type) | Station.type%in%c("DoF","Flinders"))
#Dat.scen3=subset(Dat, Station.type%in%'IMOS')

#proportion of time per jurisdiction (straight line movement assumption)
fn.prop.time.jur=function(d)  
{
  d=rbind(data.frame(Datetime=d$Datetime.prev[1],
                     Longitude=d$Longitude.prev[1],
                     Latitude=d$Latitude.prev[1],
                     zone=d$zone.prev[1],
                     Juris=d$Juris.prev[1],
                     Rel.state=d$Rel.state[1]),
          d%>%select(Datetime,Longitude,Latitude,zone,Juris,Rel.state))
  delta.t=max(d$Datetime)-min(d$Datetime)
  Rel.state=unique(d$Rel.state)
  d=d %>% mutate(date=date(Datetime)) %>%
    arrange(Datetime) %>%
    distinct(date,zone,.keep_all = TRUE) %>%
    mutate(Juris.prev=lag(Juris,1),
           Longitude.prev=lag(Longitude,1),
           Latitude.prev=lag(Latitude,1),
           date.prev=lag(date,1),
           delta.t=date-date.prev,
           Same.juris=ifelse(Juris==Juris.prev,"YES",
                             ifelse(!Juris==Juris.prev,"NO",NA)))
  d.diff.jur=subset(d,Same.juris=="NO")     
  d=subset(d,!Same.juris=="NO"|is.na(Same.juris)) 
  if(nrow(d.diff.jur)>0)
  {
    b1=vector('list',nrow(d.diff.jur))
    for(n in 1:length(b1))
    {
      dummy=with(d.diff.jur[n,],gcIntermediate(c(Longitude.prev,Latitude.prev),
                                               c(Longitude,Latitude),n=as.numeric(delta.t)-1, addStartEnd=T))
      dummy=data.frame(Longitude=dummy[,1],Latitude=dummy[,2])
      dummy$date=seq(d.diff.jur[n,'date.prev'],d.diff.jur[n,'date'],by='day')
      b1[[n]]=dummy
    }
    b1=do.call(rbind,b1)
    d.new=rbind(b1,subset(d,select=names(b1))) %>%
      arrange(date) %>%
      mutate(Juris=ifelse(Longitude>=129,"SA",
                          ifelse(Longitude<129,"WA",NA)),
             Juris.prev=lag(Juris,1),
             date.prev=lag(date,1),
             delta.t=date-date.prev)
    id=which(!d.new$Juris==d.new$Juris.prev)
    d.new$index=NA
    d.new$index[1:(id[1]-1)]=1
    if(length(id)>1)for(x in 2:length(id)) d.new$index[id[x-1]:(id[x]-1)]=x
    d.new$index[id[length(id)]:nrow(d.new)]=length(id)+1
    Tab1=group_by(subset(d.new,!is.na(delta.t)), index, Juris) %>% 
      summarise(Total.days = as.numeric(sum(delta.t))) %>%
      as.data.frame() 
    Tab1$Prop=with(Tab1,Total.days/sum(Total.days))
  }
  if(nrow(d.diff.jur)==0)
  {
    d=subset(d,!is.na(Juris))
    Tab1=data.frame(index=1,Juris=unique(d$Juris),Total.days=delta.t,Prop=1)
  }
  Tab1$Rel.state=Rel.state
  return(Tab1)
}

#proportion of time per jurisdiction
fn.plt.prop.time=function(d,CL1,CL2)
{
  plot(1:length(d),xlim=c(0,1),ylim=c(0,length(d)),col="transparent",ylab="",
       xlab="",yaxt='n',xaxt='n',bty='n')
  for(n in 1:length(d))
  {
    d[[n]]$col=with(d[[n]],ifelse(Juris=="WA",CL1,CL2))
    d[[n]]$cum.prop=cumsum(d[[n]]$Prop)
    if(nrow(d[[n]])>1)Total.days=sum(d[[n]]$Total.days)
    if(nrow(d[[n]])==1)Total.days=d[[n]]$Total.days
    
    for(r in 1:nrow(d[[n]]))
    {
      if(r==1) x1=0 else
        x1=d[[n]]$cum.prop[r-1]
      x2=d[[n]]$cum.prop[r]
      polygon(x=c(x1,rep(x2,2),x1),y=c(n-0.5,n-0.5,n,n),col=d[[n]]$col[r])
      text(mean(c(x1,x2)),mean(c(n-0.5,n)),
           paste(round(d[[n]]$Prop[r],2)*100,'%',sep=''),cex=.75)   
    }
    Clr=ifelse(d[[n]]$Rel.state[1]=="WA",CL1,CL2)
    text(0,n-.25,d[[n]]$Rel.state[1],pos=2,col=Clr,font=2)
    text(1,n-.25,round(Total.days),pos=4,cex=.85)
    text(1,n-.25,names(d)[n],pos=3,col='forestgreen',cex=.85)
    
  }
}

#residency
fn.plt.residency=function(d,CL1,CL2)
{
  plot(1:length(d),xlim=c(0,1),ylim=c(0,length(d)),col="transparent",ylab="",
       xlab="",yaxt='n',xaxt='n',bty='n')
  for(n in 1:length(d))
  {
    d[[n]]$col=with(d[[n]],ifelse(Juris=="WA",CL1,CL2))
    d[[n]]$cum.prop=cumsum(d[[n]]$Prop)
    
    for(r in 1:nrow(d[[n]]))
    {
      if(r==1) x1=0 else
        x1=d[[n]]$cum.prop[r-1]
      x2=d[[n]]$cum.prop[r]
      polygon(x=c(x1,rep(x2,2),x1),y=c(n-0.5,n-0.5,n,n),col=d[[n]]$col[r])
      text(mean(c(x1,x2)),mean(c(n-0.5,n)),
           paste(round(d[[n]]$Prop[r],2)*100,'%',sep=''),cex=.75)   
    }
    Clr=ifelse(d[[n]]$Rel.state[1]=="WA",CL1,CL2)
    text(0,n-.25,d[[n]]$Rel.state[1],pos=2,col=Clr,font=2)
    text(1,n-.4,names(d)[n],pos=3,col='forestgreen',cex=.85)
    
  }
}

#functions for plotting seasonality
fn.ggplot=function(dd,Y,X,Y.lab,X.lab)
{
  p=ggplot(dd,aes_string(y=Y,x=X,col='TagCode'))+
    geom_line()+geom_point()+scale_y_reverse()+
    labs(y = Y.lab,x = X.lab)
  if(X=='Mn')
  {
    p=p+scale_x_continuous(breaks = seq(1,12,2),
                         label = month.abb[seq(1,12,2)])
  }
  p
}

#ACA, missing  gam?? how to present, ask Charlie for his graph....
#function for wrapping scenarios
Adelaide=c(138.6,-34.928)
fun.run.scen=function(Dat,SCEN)
{
  #create useful objects
  state=unique(Dat$State)
  TAG=unique(Dat$TagCode)
  Sp=unique(Dat$Species)
  TAG.species=vector('list',length(Sp))
  names(TAG.species)=Sp
  for(t in 1:length(TAG.species)) TAG.species[[t]]=unique(subset(Dat,Species==Sp[t])$TagCode)
  
  #preliminary stuff
  Lon.range=range(Dat$Longitude,na.rm=T)
  Lon.range[2]=ifelse(Lon.range[2]<129,140,Lon.range[2])
  Lat.range=range(Dat$Latitude,na.rm=T)
  Lat.range[1]=ifelse(Lat.range[1]>(-35),-38,Lat.range[1])
  XLIM=c(Lon.range[1],Lon.range[2])
  YLIM=c(Lat.range[1],Lat.range[2])
  if(do.expl=="YES")
  {
    fn.plt1=function(TG)
    {
      a=subset(Dat,TagCode==TG)
      plot(a$Longitude,a$Latitude,ylim=YLIM,xlim=XLIM,las=1,ylab="Lat",xlab="Long",
           main=paste(unique(a$Species),"_",unique(a$TagCode)," (n=",nrow(a),
                      " detections)",sep=""))
    }
    pdf("Results/Exploratory.pdf") 
    sapply(TAG,fn.plt1)
    dev.off() 
  }
  
  #straight line distances (in km) between consecutive detections
  #       applying algorithm to avoid going over land
  Dat$Distance=distGeo(Dat[,c("Longitude.prev","Latitude.prev")],
                       Dat[,c("Longitude","Latitude")])
  Dat$Distance.c=Dat$Distance
  Dat$Distance.c=ifelse(Dat$zone.prev=='SA.east' & Dat$zone=='Zone2',
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Eyre)+
                          Eyre_SA.Border+
                          distGeo(SA.Border,Dat[,c("Longitude","Latitude")]),
                  ifelse(Dat$zone.prev=='SA.east' & Dat$zone=='Zone1' &
                          Dat$Longitude>= Cape.Leuwin[1] & Dat$Latitude<= Cape.Leuwin[2],
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Eyre)+
                          Eyre_SA.Border+SA.Border_Mid.point+
                          distGeo(Mid.point,Dat[,c("Longitude","Latitude")]),   
                  ifelse(Dat$zone.prev=='SA.east' & Dat$zone%in%c('Zone1','WC') & 
                          Dat$Latitude>Cape.Leuwin[2],
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Eyre)+
                          Eyre_SA.Border+SA.Border_Mid.point+Mid.point_Cape.Leuwin+
                          distGeo(Cape.Leuwin,Dat[,c("Longitude","Latitude")]), 
                                      
                  ifelse(Dat$zone.prev=='Zone2' & Dat$zone=='SA.east',
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],SA.Border)+
                          Eyre_SA.Border+
                          distGeo(Eyre,Dat[,c("Longitude","Latitude")]),   
                  ifelse(Dat$zone.prev=='Zone2' & Dat$zone=='Zone1' &
                          Dat$Longitude>= Cape.Leuwin[1] & Dat$Latitude<= Cape.Leuwin[2],
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Mid.point)+
                          distGeo(Mid.point,Dat[,c("Longitude","Latitude")]),
                  ifelse(Dat$zone.prev=='Zone2' & Dat$zone%in%c('Zone1','WC') & 
                          Dat$Latitude>Cape.Leuwin[2],
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Mid.point)+
                          Mid.point_Cape.Leuwin+
                          distGeo(Cape.Leuwin,Dat[,c("Longitude","Latitude")]),
                                                           
                  ifelse(Dat$zone.prev%in%c('Zone1','WC') & Dat$zone=='SA.east' & 
                          Dat$Latitude.prev> Cape.Leuwin[2],
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Cape.Leuwin)+
                          Mid.point_Cape.Leuwin+SA.Border_Mid.point+Eyre_SA.Border+
                          distGeo(Eyre,Dat[,c("Longitude","Latitude")]), 
                  ifelse(Dat$zone.prev%in%c('Zone1') & Dat$zone=='SA.east' & 
                          Dat$Latitude.prev<= Cape.Leuwin[2],
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Mid.point)+
                          SA.Border_Mid.point+Eyre_SA.Border+
                          distGeo(Eyre,Dat[,c("Longitude","Latitude")]), 
                                                                       
                  ifelse(Dat$zone.prev%in%c('Zone1','WC') & Dat$zone=='Zone2' &
                          Dat$Longitude< Cape.Leuwin[1] & Dat$Latitude> Cape.Leuwin[2],
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Cape.Leuwin)+
                          Mid.point_Cape.Leuwin+
                          distGeo(Mid.point,Dat[,c("Longitude","Latitude")]),
                  ifelse(Dat$zone.prev%in%c('Zone1') & Dat$zone=='Zone2' &
                          Dat$Longitude>= Cape.Leuwin[1] & Dat$Latitude<= Cape.Leuwin[2],
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Mid.point)+
                          distGeo(Mid.point,Dat[,c("Longitude","Latitude")]),
                  ifelse(Dat$zone.prev%in%c('Zone1','WC') & Dat$zone=='Zone1' &
                          Dat$Longitude< Cape.Leuwin[1] & Dat$Latitude> Cape.Leuwin[2],
                          distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Cape.Leuwin)+
                          distGeo(Cape.Leuwin,Dat[,c("Longitude","Latitude")]),
                          Dat$Distance.c)))))))))))
  
  Dat$Distance=with(Dat,ifelse(TagCode==TagCode.prev,Distance/1000,NA)) 
  Dat$Distance.c=with(Dat,ifelse(TagCode==TagCode.prev,Distance.c/1000,NA)) 
  
  
  #Distribution of displacements and ROM across jurisdictions
  Dat= Dat%>% 
    mutate(Juris=ifelse(zone%in%c("WC","Zone1","Zone2"),"WA",
                        ifelse(zone%in%c("SA.east","SA.west"),"SA",NA)),
           Juris.prev=ifelse(zone.prev%in%c("WC","Zone1","Zone2"),"WA",
                             ifelse(zone.prev%in%c("SA.east","SA.west"),"SA",NA)),
           Same.juris=ifelse(Juris==Juris.prev,"YES",
                             ifelse(!Juris==Juris.prev,"NO",NA)))
  
  Cros.jur=Dat%>%filter(Same.juris=="NO")%>%
    select(Species,TagCode,TagCode.prev,Juris,Juris.prev,
           Longitude,Latitude,Longitude.prev,Latitude.prev,
           Datetime,Datetime.prev,Distance.c,Time)%>%
    mutate(Time=ifelse(is.na(Time),difftime(Datetime,Datetime.prev,units='mins'),Time))
  
  Cros.jur$Distance.c=ifelse(is.na(Cros.jur$Distance.c),distGeo(Cros.jur[,c("Longitude.prev","Latitude.prev")],
                                  Cros.jur[,c("Longitude","Latitude")])/1000,Cros.jur$Distance.c)
  Cros.jur=Cros.jur%>%mutate(ROM=Distance.c/(Time/(24*60)))  # km/day
  
  tiff(file=paste('displacement and ROM/Cross.juris_hist_figure_',SCEN,'.tiff',sep=''),
       width=2400,height=2400,units="px",res=300,compression="lzw+p")
  par(mfrow=c(2,1),mar=c(1,2,1,.1),oma=c(.5,3,.1,.5),las=1,
      mgp=c(2,.5,0),cex.axis=1.25,cex.lab=1.5,xpd=T)
  
    #distance
   vioplot(Distance.c~Species,Cros.jur, horizontal=F, col="gray",ylab='')
  mtext("Distance (km)",2,line=3.25,las=3,cex=1.75)

    #ROM
  vioplot(ROM~Species,Cros.jur, horizontal=F, col="gray",ylab='')
  mtext("ROM (km/day)",2,line=3.25,las=3,cex=1.75)
  dev.off()
  
  ROM.sumery=tapply(Cros.jur$ROM, Cros.jur$Species, summary)
  Dist.sumery=tapply(Cros.jur$Distance.c, Cros.jur$Species, summary)
  Tab2<-rbind(do.call(rbind,Dist.sumery),do.call(rbind,ROM.sumery))
  rownames(Tab2)=paste(c(rep('Distance',2),rep('ROM',2)),rownames(Tab2),sep='_')
  write.csv(Tab2,paste('displacement and ROM/Cross.juris_summary_',SCEN,'.csv',sep=''))
  
  
  #cross jurisdictional displacements
  Prop.time=vector('list',length(TAG))
  names(Prop.time)=TAG
  for(i in 1:length(TAG))
  {
    Prop.time[[i]]=fn.prop.time.jur(d=subset(Dat,TagCode==TAG[i],
                                             select=c(Datetime,Longitude,Latitude,zone,Juris,Rel.state,
                                                      Datetime.prev,Longitude.prev,Latitude.prev,zone.prev,Juris.prev)))
  }
  
  tiff(file=paste('figure_prop.time_',SCEN,'.tiff',sep=''),width=1800,height=2400,units="px",res=300,
       compression="lzw+p")
  par(mfcol=c(2,1),mar=c(1,1,1,1.25),oma=c(.1,.1,.1,.5),las=1,
      mgp=c(1.25,.35,0),cex.axis=1.1,cex.lab=1.25,xpd=T)
  for(t in 1:length(TAG.species)) 
  {
    id=which(names(Prop.time)%in%TAG.species[[t]])
    fn.plt.prop.time(d=Prop.time[id],CL1=Cols[1],CL2=Cols[2])
    if(t==2)legend('bottom',c("WA","SA"),fill=Cols,horiz=T,bty='n')
    mtext(names(TAG.species)[t],3,cex=1.5)
  }
  mtext("Percent of time per zone",1,line=-2.5,cex=1.5)
  dev.off()
  
  
  #seasonal patterns plot distance from release thru time 
  Spics=Dat%>%distinct(TagCode,.keep_all = T)%>%select(TagCode,Species)
  Seasonal=Dat%>%
            left_join(subset(Rel.dat,select=c(TagCode,ReleaseDate,
                                              ReleaseLatitude,ReleaseLongitude)),by='TagCode')%>%
            arrange(TagCode,Datetime)%>%
            mutate(Delta.t=as.numeric(difftime(Datetime,ReleaseDate,units='days')),
                   AdelaideLon=Adelaide[1],
                   AdelaideLat=Adelaide[2])
  Seasonal$Dist.frm.rel=ifelse(Seasonal$TagCode==Seasonal$TagCode.prev,
                           distGeo(Seasonal[,c("ReleaseLongitude","ReleaseLatitude")],
                                   Seasonal[,c("Longitude","Latitude")])/1000,NA)
  Seasonal$Dist.frm.Adld=ifelse(Seasonal$TagCode==Seasonal$TagCode.prev,
                               distGeo(Seasonal[,c("AdelaideLon","AdelaideLat")],
                                       Seasonal[,c("Longitude","Latitude")])/1000,NA)
  
  Ril.dt=subset(Rel.dat,TagCode%in%unique(Seasonal$TagCode))
  Add.relis=Seasonal[1:length(unique(Seasonal$TagCode)),]
  Add.relis[,]=NA
  Add.relis=Add.relis%>%
                mutate(TagCode=Ril.dt$TagCode,
                               TagCode.prev=Ril.dt$TagCode,
                               Datetime=Ril.dt$ReleaseDate,
                               Dist.frm.rel=0,
                               ReleaseLatitude=Ril.dt$ReleaseLatitude,
                               ReleaseLongitude=Ril.dt$ReleaseLongitude,
                               AdelaideLon=Adelaide[1],
                               AdelaideLat=Adelaide[2],
                               Delta.t=0)%>%
                select(-Species)%>%
                left_join(Spics,by='TagCode')%>%
                select(names(Seasonal))
  Add.relis$Dist.frm.Adld=ifelse(Add.relis$TagCode==Add.relis$TagCode.prev,
                                distGeo(Add.relis[,c("AdelaideLon","AdelaideLat")],
                                        Add.relis[,c("ReleaseLongitude","ReleaseLatitude")])/1000,NA)
  
  Seasonal=rbind(Seasonal,Add.relis)%>%
        arrange(TagCode,Delta.t)%>%
        select(TagCode,TagCode.prev,Species,Datetime,
                             Dist.frm.rel,Dist.frm.Adld,Delta.t)%>%
        filter(!is.na(Dist.frm.rel))%>%
        mutate(Mn=month(Datetime),
               Yr=year(Datetime),
               ln.Dist.frm.rel=log(Dist.frm.rel+1e-4),
               ln.Dist.frm.Adld=log(Dist.frm.Adld),
               scaled.Dist.frm.Adld=scale(Dist.frm.Adld))
  
  #Bronzie    
  Seasonal.bronzie=Seasonal%>%
    filter(Species=='bronze whaler')
  Sisonls.bronzie=c(31003,31000,30894,29587,29542,27698)
  Seasonal.bronzie.Mod=Seasonal.bronzie%>%
                    filter(TagCode%in%Sisonls.bronzie)%>%
                    mutate(TagCode=as.factor(TagCode))%>%
                    arrange(TagCode,Datetime)

  fn.ggplot(dd=Seasonal.bronzie.Mod,
            Y='Dist.frm.Adld',
            X='Mn',
            Y.lab="Distance from Adelaide (km)",
            X.lab="Month")
  
  fn.ggplot(dd=Seasonal.bronzie.Mod,
            Y='Dist.frm.Adld',
            X='Delta.t',
            Y.lab="Distance from Adelaide (km)",
            X.lab="Delta.t")
  

  
  Mod.bronzie=gam(ln.Dist.frm.Adld~s(Delta.t,bs='cc')+s(TagCode,bs='re'),data=Seasonal.bronzie.Mod)
  Mod.bronzie=gam(ln.Dist.frm.Adld~s(Mn,bs='cc', k = 12)+s(TagCode,bs='re'),data=Seasonal.bronzie.Mod)
  
  newd=data.frame(Mn=seq(1,12,length.out=100),TagCode=)
  pred <- predict.gam(Mod.bronzie,newd,type = "response", se.fit = TRUE)
  plot(newd$Mn,exp(pred$fit))
  segments(newd$Mn,exp(pred$fit+1.96*pred$se.fit),
         newd$Mn,exp(pred$fit-1.96*pred$se.fit))
  
  
  #Dusky
  Seasonal.dusky=Seasonal%>%
    filter(Species=='Dusky')
  Sisonls.dusky=c(49144,49146)
  Seasonal.dusky.Mod=Seasonal.dusky%>%
    filter(TagCode%in%Sisonls.dusky)%>%
    mutate(TagCode=as.factor(TagCode))%>%
    arrange(TagCode,Datetime)
  
  fn.ggplot(dd=Seasonal.dusky.Mod,
            Y='Dist.frm.Adld',
            X='Mn',
            Y.lab="Distance from Adelaide (km)",
            X.lab="Month")
  
  fn.ggplot(dd=Seasonal.dusky.Mod,
            Y='Dist.frm.Adld',
            X='Delta.t',
            Y.lab="Distance from Adelaide (km)",
            X.lab="Delta.t")
  

  Mod.dusky=gam(ln.Dist.frm.Adld~s(Delta.t,bs='cc')+s(TagCode,bs='re'),data=Seasonal.dusky.Mod)
  ggplot(Seasonal.dusky.Mod,aes(x=Delta.t,y=ln.Dist.frm.Adld,col=TagCode))+
    geom_line()
 # Mod.dusky=gam(ln.Dist.frm.Adld~s(Mn,bs='cc')+s(TagCode,bs='re'),data=Seasonal.dusky.Mod)
  ggplot(Seasonal.dusky.Mod,aes(x=Delta.t,y=ln.Dist.frm.Adld,col=TagCode))+
    geom_line()
  

  #Residency
  Residency=Dat%>%
        arrange(TagCode,Datetime)%>%
        mutate(event=ifelse(TagCode==TagCode.prev & Same.juris=="NO",1,0),
               event=ifelse(TagCode==TagCode.prev,cumsum(event),NA),
               TagCode.Juris=paste(TagCode,Juris,Juris.prev))%>%
        group_by(TagCode) %>%
        mutate(Tot.time = difftime(max(Datetime),min(Datetime.prev),units='mins'))%>%   
        group_by(TagCode,event,TagCode.Juris)%>%
        mutate(Time.in.zn=ifelse(TagCode==TagCode.prev & Same.juris=="YES",
                difftime(max(Datetime),min(Datetime.prev),units='mins'),NA))%>%
        mutate(Residency=Time.in.zn/as.numeric(Tot.time))%>%
        distinct(TagCode,event,TagCode.Juris,.keep_all = T)%>%
        select(TagCode,Species,Rel.state,Juris,Residency,event,TagCode.Juris)%>%
        filter(!is.na(Residency))%>%
        data.frame()%>%mutate(Event.Juris=paste(event,Juris))
  
  Prop.time.Res=vector('list',length(unique(Residency$Species)))
  names(Prop.time.Res)=unique(Residency$Species)
  for(j in 1:length(Prop.time.Res))
  {
    d=Residency%>%filter(Species==names(Prop.time.Res)[j])
    this=unique(d$TagCode)
    dummy=vector('list',length(this))
    names(dummy)=this
    for(xx in 1:length(dummy))
    {
      dummy[[xx]]=d%>%filter(TagCode==names(dummy)[xx])%>%
        mutate(index=1:length(event),
               Prop=Residency,
               Total.days='')%>%
        select(index,Juris,Total.days,Prop,Rel.state)
    }
    Prop.time.Res[[j]]=dummy
  }

  tiff(file=paste('figure_residency_',SCEN,'.tiff',sep=''),width=1800,height=2400,units="px",res=300,
        compression="lzw+p")
  par(mfcol=c(2,1),mar=c(1,1,1,1.25),oma=c(.1,.1,.1,.5),las=1,
      mgp=c(1.25,.35,0),cex.axis=1.1,cex.lab=1.25,xpd=T)
  for(j in 1:length(Prop.time.Res))
  {
    fn.plt.residency(d=Prop.time.Res[[j]],CL1=Cols[1],CL2=Cols[2])
    if(j==2)legend('bottom',c("WA","SA"),fill=Cols,horiz=T,bty='n')
    mtext(names(Prop.time.Res)[j],3,cex=1.5)
  }
  mtext("Residency",1,line=-2.5,cex=1.5)
  dev.off()
  
}
Cols=c("steelblue","pink2")
names(Cols)=c("WA","SA")

fun.run.scen(Dat=Dat.scen1,SCEN="Scen1")
fun.run.scen(Dat=Dat.scen2,SCEN="Scen2")
#fun.run.scen(Dat=Dat.scen3,SCEN="Scen3")


# Conventional tagging----------------------------------------------------------------
Conv.Tagging=Tagging%>%
  filter(Species%in%c("BW","CP") & !is.na(Lat.rels))%>%
  mutate(Recaptured=ifelse(!is.na(Yr.rec)|!is.na(Lat.rec)|!is.na(Long.rec),"Yes","No"),
         Rel.juris=ifelse(Long.rels>129,"SA","WA"),
         Rec.juris=ifelse( Recaptured=="Yes"& Long.rec>129,"SA","WA"),
         Col.sp=ifelse(Species=="BW",2,3),
         Same.Juris=ifelse(Recaptured=="Yes"& Rel.juris==Rec.juris,"Yes",
                           ifelse(Recaptured=="Yes"& !Rel.juris==Rec.juris,"No",NA)))
TAB.conv.rel.rec=Conv.Tagging%>%
  group_by(Rel.juris,Rec.juris,COMMON_NAME,Recaptured)%>%
  summarise(n=n())%>%
  filter(!is.na(Rec.juris))%>%
  arrange(COMMON_NAME,Rel.juris,Recaptured)%>%
  data.frame%>%
  mutate(Rec.juris=ifelse(Recaptured=="No","N/A",Rec.juris))

write.csv(TAB.conv.rel.rec,"Summary.conv.tag.csv",row.names = F)

with(subset(Conv.Tagging,Recaptured=="Yes" & Same.Juris=="No"),{
  plot(1,1,ylim=c(-36,-30),xlim=c(113,140),ylab="",xlab="")
  arrows(Long.rels,Lat.rels,Long.rec,Lat.rec,col=Col.sp)
})
Siz.sex.cross.conv=Conv.Tagging%>%
  filter(Recaptured=="Yes" & Same.Juris=="No")%>%
  mutate(Date.rel=as.POSIXct(paste(Yr.rel,Mn.rel,Day.rel,sep="-")),
         Date.rec=as.POSIXct(paste(Yr.rec,Mn.rec,Day.rec,sep="-")),
         time.at.liberty=Date.rec-Date.rel)%>%
  select(Species,COMMON_NAME,Lat.rels,Long.rels,Lat.rec,Long.rec,
         Rel.juris,Rec.juris,Rel_FL,Sex,Date.rel,
         Date.rec,time.at.liberty)%>%
  arrange(Species,Date.rel,Sex)
write.csv(Siz.sex.cross.conv,"Summary.conv.tag_time.liberty.csv",row.names = F)







