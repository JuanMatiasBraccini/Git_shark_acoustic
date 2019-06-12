#--SCRIPT FOR ANALYSING CROSS-JURISDICTIONAL MOVEMENTS OF DUSKY & BRONZE WHALERS--

#note: must report the total number of WA and SA tags to put those that 
#      move across in perspective...

#to do for individuals that crossed jurisdictions:
   # Straight line distances between acoustic arrays (indicate travel distance)
   # Connectivity plots (The connectivity plot represented the frequency 
   #     and magnitude of daily shark movements, and therefore provided information
   #     about the incoming and outgoing movements of duskies and bronzies along  
   #     the southern coast of Australia
   #     use "circus.trackPlotRgiong()" from 'circlize' package)...this needs the total number of tagged sharks
   # Proportion of time per jurisdiction (WA-SA-WA-SA, etc)
   # Quantification of cross jurisdictional displacements and rates of movements
   # Any seasonal patterns??

   #Should also consider conventional tags of duskies that move to SA (any bronzies?)

#Clarify with Charlie: for each of these sharks, detections are the total number

rm(list=ls(all=TRUE))
library(lubridate)
library(dplyr)
library(geosphere)


#1. Data Section---------------------------------------------------------
setwd('C:\\Users\\myb\\Desktop\\new\\Charlie H')
Dat=read.csv('Data/All detections_2018 01.csv',stringsAsFactors = F)
SA.receivers=read.csv('Data/GSV receiver location.csv',stringsAsFactors = F)
Blank_LatLong=read.csv('Data/Blank_LatLong.csv',stringsAsFactors = F)

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
Dat_no.name=subset(Dat,StationName=="") %>% select(-StationName)
Dat_no.name=left_join(Dat_no.name,Blank_LatLong,by="ReceiverSerialNumber")
Dat=subset(Dat,!StationName=="") %>% select(names(Dat_no.name))
Dat=rbind(Dat,Dat_no.name)
Dat = left_join(Dat,SA.receivers,by=c("StationName" = "Location.name")) %>%
      mutate(Datetime=as.POSIXct(Datetime..WST.,format="%d/%m/%Y %H:%M",tz="UTC"),
              Mn=month(Datetime),
              Yr=year(Datetime),
              Longitude=ifelse(is.na(Longitude.x),Longitude.y,Longitude.x),
              Latitude=ifelse(is.na(Latitude.x),Latitude.y,Latitude.x),
              Latitude=-abs(Latitude),
              N=1,
              zone=ifelse(Latitude>(-33) & Latitude<=(-26) & Longitude<116.5,"WC",
                   ifelse(Latitude<=(-33) & Longitude<116.5,"Zone1",
                   ifelse(Longitude>=116.5 & Longitude<129 & Latitude<=(-26),"Zone2",
                   ifelse(Longitude>=129 & Longitude<135.7,"SA.west",
                   ifelse(Longitude>=135.6,"SA.east",NA)))))) %>%
      arrange(TagCode, Datetime) %>%
      mutate(TagCode.prev=lag(TagCode,1),
             Datetime.prev=lag(Datetime,1),
             Longitude.prev=lag(Longitude,1),
             Latitude.prev=lag(Latitude,1),
             State.prev=lag(State,1),
             zone.prev=lag(zone,1),
             Same.station=ifelse(TagCode==TagCode.prev & 
                                 Longitude==Longitude.prev &
                                 Latitude==Latitude.prev  ,"YES","NO"),
             Time=ifelse(TagCode==TagCode.prev,(Datetime-Datetime.prev)/60,NA)) %>%   #time in minutes
      select(c(TagCode,TagCode.prev,Species,Organisation,State,State.prev,
               Datetime,Datetime.prev,Mn,Yr,StationName,Longitude,Longitude.prev,
               Latitude,Latitude.prev,Same.station,zone,zone.prev,Time,N))
      
if(sum(is.na(Dat$Latitude))>0) cat("------",sum(is.na(Dat$Latitude)),"RECORDS HAVE NO LATITUDE-----")

#3.2 create useful objects
state=unique(Dat$State)
TAG=unique(Dat$TagCode)
Sp=unique(Dat$Species)
TAG.species=vector('list',length(Sp))
names(TAG.species)=Sp
for(t in 1:length(TAG.species)) TAG.species[[t]]=unique(subset(Dat,Species==Sp[t])$TagCode)

#3.3 preliminary stuff
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

#3.4 straight line distances (in km) between consecutive detections
#note: apply algorithm to avoid going over land
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


#3.5 cross jurisdictional displacements and rates of movements (km per day)
Dat= Dat%>% 
  mutate(ROM=Distance.c/(Time/(60*24)),
         Juris=ifelse(zone%in%c("WC","Zone1","Zone2"),"WA",
               ifelse(zone%in%c("SA.east","SA.west"),"SA",NA)),
         Juris.prev=ifelse(zone.prev%in%c("WC","Zone1","Zone2"),"WA",
                    ifelse(zone.prev%in%c("SA.east","SA.west"),"SA",NA)),
         Same.juris=ifelse(Juris==Juris.prev,"YES",
                    ifelse(!Juris==Juris.prev,"NO",NA))
        )


#3.6 proportion of time per jurisdiction (straight line movement assumption)
fn.prop.time.jur=function(d)
{
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
    Tab1=Tab1%>% select(index,Juris,Prop)
  }
  if(nrow(d.diff.jur)==0)
  {
    d=subset(d,!is.na(Juris))
    Tab1=data.frame(index=1,Juris=unique(d$Juris),Prop=1)
  }
   return(Tab1)
}
Prop.time=vector('list',length(TAG))
names(Prop.time)=TAG
for(i in 1:length(TAG))
{
  Prop.time[[i]]=fn.prop.time.jur(d=subset(Dat,TagCode==TAG[i],
                  select=c(Datetime,Longitude,Latitude,zone,Juris)))
}

#ACA
#3.7 connectivity plots 


#3.8 seasonal patterns??


#4. Report Section-------------------------------------------------------
setwd(paste(getwd(),"Results",sep="/"))

#4.1 table of tagcodes by species and state
Tab1= group_by(Dat, TagCode, Species,Organisation,State) %>%
      summarise(sum = sum(N)) %>%
      as.data.frame()
write.csv(Tab1,'Tab.tag.code_species_org_state.csv',row.names = F)

#4.2 proportion of time per jurisdiction
fn.plt.prop.time=function(d,CL1,CL2)
{
  plot(1:length(d),xlim=c(0,1),ylim=c(0,length(d)),col="transparent",ylab="",
       xlab="Proportion of time",yaxt='n')
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
          paste(round(d[[n]]$Prop[r],2)*100,'%',sep=''),cex=.85)                 
    }
  }
}

Cols=c("steelblue","pink")
names(Cols)=c("WA","SA")

jpeg(file='figure_prop.time.jpg',width=1800,height=2400,units="px",res=300)
par(mfcol=c(2,1),mar=c(2.5,1,1,.1),oma=c(.1,.1,.1,.5),las=1,mgp=c(1.25,.35,0),cex.axis=1.1,cex.lab=1.25)
for(t in 1:length(TAG.species)) 
{
  id=which(names(Prop.time)%in%TAG.species[[t]])
  fn.plt.prop.time(d=Prop.time[id],CL1=Cols[1],CL2=Cols[2])
  if(t==2)legend('bottom',c("WA","SA"),fill=Cols,horiz=T,bty='n')
  mtext(names(TAG.species)[t],3,cex=1.5)
}
dev.off()