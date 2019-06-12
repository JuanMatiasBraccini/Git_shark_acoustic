#GOOGLE MAPPING FOR FRDC SPATIAL AND TEMPORAL PATTERNS

#notes: this script maps the location of VR2s and hits received from tagged dusky, sandbar, whiskey and gummy sharks
#       remember to use latest version of input files

 #note: keep in mind faulty receivers must be flagged to avoid assuming that if no hits, no sharks!!!


setwd("H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Mapping/FRDC_tagging")

library(RgoogleMaps)  #library for plotting on google maps
library(RODBC)				#library for importing excel data
library(gridBase)			#for inset map
library(PBSmapping)   #for polygon
library(animation)    #for animations



windows(record=T)			#to keep al graphs

memory.limit(3900)   #set memory limit to maximum size


###### DATA SECTION ############

  #1. load data
    #1.1. listening stations                                                              #UPDATE
#note: these files need to be updated as receiver deployment progresses

path="M:/Fisheries Research/FinFish/Shark/Braccini/Ningaloo October 2011 Data/"


      #PERTH AND SOUTHWEST RECEIVERS
#SouthWest<- read.csv("H:/Matias WA Fisheries/Data/Mapping/AcousticLinesSWA.csv")
channel <- odbcConnectExcel("H:/Matias WA Fisheries/Data/Mapping/Receiver_locations_March_2012")
SMN_VR2W_Hamelin=sqlFetch(channel,"SMN_VR2W_Hamelin")
SMN_VR2W_Albany=sqlFetch(channel,"SMN_VR2W_Albany")
CSIRO_Chatham=sqlFetch(channel,"CSIRO_Chatham")
DoF_Chatham=sqlFetch(channel,"DoF_Chatham")

SMN_VR4G_metro=sqlFetch(channel,"SMN_VR4G_metro")
SMN_VR2W_metro=sqlFetch(channel,"SMN_VR2W_metro")
OTN_VR2W_metro=sqlFetch(channel,"OTN_VR2W_metro")
close(channel)

SouthWest=rbind(SMN_VR2W_Hamelin,SMN_VR2W_Albany,CSIRO_Chatham,DoF_Chatham)

      #Ningaloo
channel <- odbcConnectExcel(paste(path,"Field Sheets/Field Sheets October 2011",sep=""))
Ningaloo<- sqlFetch(channel,"Recoveries_Field_Sheets")
close(channel)
if (is.na(Ningaloo$"Receiver (S/N)"[1])) Ningaloo=Ningaloo[-1,]
Ningaloo=Ningaloo[!(is.na(Ningaloo$"Receiver (S/N)")),]
Ningaloo$"Recovery latitude"=-Ningaloo$"Recovery latitude"
colnames(Ningaloo)[match(c("Recovery latitude","Recovery longitude"),names(Ningaloo))]=c("Recovery_latitude" , "Recovery_longitude")

      #Perth
channel1 <- odbcConnectExcel("M:/Fisheries Research/FinFish/Shark/Braccini/Shark MonitoringOTN/SMN_OTN Running Sheet 2011 copy")
SMN_VR2_2010<- sqlFetch(channel1,"SMN VR2W 2010")
SMN_VR4G_2009<- sqlFetch(channel1,"VR4G Deployed 2009")
SMN_VR4G_2009$Long[2]=115.73338       #fill in R import stuff ups
SMN_VR4G_2009$Long[17]=115.47654
OTN_2009a<- sqlFetch(channel1,"OTN 2009a")
close(channel1)


    #1.2. shark detections                                                              #UPDATE
#note: these files need to be updated as tagging progresses

all.csv=list.files(path,pattern="csv")    #get all .csv files
Tagging.Data.Ningaloo=NULL
for(i in 1:length(all.csv))
{
  Tagging.Data.Ningaloo<- rbind(Tagging.Data.Ningaloo,read.csv(paste(path,all.csv[i],sep="")))
}

path.Perth="M:/Fisheries Research/FinFish/Shark/Braccini/PerthDummy/"
all.csv.Perth=list.files(path.Perth,pattern="csv")    #get all .csv files
Tagging.Data.Perth=NULL
for(i in 1:length(all.csv.Perth))
{
  Tagging.Data.Perth<- rbind(Tagging.Data.Perth,read.csv(paste(path.Perth,all.csv.Perth[i],sep="")))
}

path.SW="M:/Fisheries Research/FinFish/Shark/Braccini/SouthWestDummy/"
all.csv.SW=list.files(path.SW,pattern="csv")    #get all .csv files
Tagging.Data.SouthWest=NULL
for(i in 1:length(all.csv.SW))
{
  Tagging.Data.SouthWest<- rbind(Tagging.Data.SouthWest,read.csv(paste(path.SW,all.csv.SW[i],sep="")))
}


    #1.3. tag deployment data                                                              #UPDATE
#note: these files need to be updated as tagging progresses
Tag.Deployed=read.csv("H:/Matias WA Fisheries/Data/Tags_ID/AcousticsDeployed WA JunAug2011.csv")
Tag.Deployed$DateDeployed=as.Date(Tag.Deployed$DateDeployed,"%d/%m/%Y")      #convert to date
Tag.Deployed$MID.LAT=-Tag.Deployed$MID.LAT


  #2. specify mapping area
    #2.1. Australia
minLat=-38; maxLat=-20; minLong=112; maxLong=140;
lat=c(maxLat,minLat);long=c(minLong,maxLong)
center=c(mean(lat),mean(long))
#zoom=min(MaxZoom(range(lat),range(long)))
zoom=4
    #2.2. WA
lat.wa=c(-15,-45);long.wa=c(110,120)
#center.wa=c(mean(lat.wa),mean(long.wa))
center.wa=c(-26,112)
#zoom.wa=max(MaxZoom(range(lat.wa),range(long.wa)))
zoom.wa=5
 



###### PROCEDURE SECTION ############

#1. Manipulate Perth receivers
SMN_VR2_2010$Lat=ifelse(SMN_VR2_2010$Lat>0,-SMN_VR2_2010$Lat,SMN_VR2_2010$Lat)
thesecolumns=match(c("Station No","Lat","Long","Serial No"),names(SMN_VR2_2010))
SMN_VR2_2010=SMN_VR2_2010[,thesecolumns]
thesecolumnsVR4=match(c("Station No","Lat","Long","VR4G Serial No"),names(SMN_VR4G_2009))
thesecolumnsOTN09a=match(c("Station No","Lat","Long","Serial No"),names(OTN_2009a))
SMN_VR4G_2009=SMN_VR4G_2009[,thesecolumnsVR4]
OTN_2009a=OTN_2009a[,thesecolumnsOTN09a]
names(SMN_VR4G_2009)=names(OTN_2009a)=names(SMN_VR2_2010)
SMN_VR2_2010$RecType="VR2.SMN"
SMN_VR4G_2009$RecType="VR4.SMN"
OTN_2009a$RecType="VR2.OTN"
Perth=rbind(SMN_VR2_2010,OTN_2009a,SMN_VR4G_2009)
Perth=Perth[!(is.na(Perth$Lat)),]

#2. Manipulate tagging data and convert UTC detection time to Perth local time     (+8hours)
Tagging.Data.Ningaloo$DS="Ningaloo"
colnames(Tagging.Data.Ningaloo)[1]="Date.and.Time..UTC."
Tagging.Data.Ningaloo$Date.time=as.POSIXlt(as.POSIXct(as.character(Tagging.Data.Ningaloo$"Date.and.Time..UTC."),
tz="GMT"), "Australia/Perth")
Tagging.Data.Ningaloo$Date=as.Date(Tagging.Data.Ningaloo$Date.time)

Tagging.Data.Perth$DS="Perth"
Tagging.Data.Perth$Date.time=strptime(as.character(Tagging.Data.Perth$"Date.and.Time..UTC."),format='%d/%m/%Y %H:%M')
Tagging.Data.Perth$Date.time=as.POSIXlt(as.POSIXct(Tagging.Data.Perth$Date.time,tz="GMT"), "Australia/Perth")
Tagging.Data.Perth$Date=as.Date(Tagging.Data.Perth$Date.time)

Tagging.Data.SouthWest$DS="SouthWest"
Tagging.Data.SouthWest$Date.time=strptime(as.character(Tagging.Data.SouthWest$"Date.and.Time..UTC."),format='%d/%m/%Y %H:%M')
Tagging.Data.SouthWest$Date.time=as.POSIXlt(as.POSIXct(Tagging.Data.SouthWest$Date.time,tz="GMT"), "Australia/Perth")
Tagging.Data.SouthWest$Date=as.Date(Tagging.Data.SouthWest$Date.time)

Tagging.Data=rbind(Tagging.Data.Ningaloo,Tagging.Data.Perth,Tagging.Data.SouthWest)

      #extract tagged individuals by species
tagged.species=table(Tag.Deployed$Species)
Species=unique(as.character(Tag.Deployed$Species))


#3. Extract receiver and tag numbers
Tagging.Data$Receiver.Number=unlist(strsplit(as.character(Tagging.Data$Receiver), split="-"))[seq(2,2*nrow(Tagging.Data),2)]
Tagging.Data$Tag=unlist(strsplit(as.character(Tagging.Data$Transmitter), split="-"))[seq(3,3*nrow(Tagging.Data),3)]



#4. Merge receiver position and hits                            (FIX FOR REAL DATA FROM PERTH AND SOUTH WEST)
these.cols.Ning=match(c("Receiver (S/N)","Recovery_latitude","Recovery_longitude","Comments"),names(Ningaloo))
Ning.Rec=Ningaloo[,these.cols.Ning]
these.cols.Perth=match(c("Serial No","Lat","Long","RecType"),names(Perth))
Perth.Rec=Perth[,these.cols.Perth]
these.cols.SouthWest=match(c("Mooring.Number","Lat","Long","Line"),names(SouthWest))
SouthWest.Rec=SouthWest[,these.cols.SouthWest]
names(Perth.Rec)=names(SouthWest.Rec)=names(Ning.Rec)

All.Rec=rbind(Ning.Rec,Perth.Rec,SouthWest.Rec)


these.cols.Tagging=match(c("Date.time","Date","Receiver.Number","Tag","DS"),names(Tagging.Data))

Tagging.Data=merge(Tagging.Data[,these.cols.Tagging],All.Rec,by.y="Receiver (S/N)",
by.x="Receiver.Number",all.x=T)



#5. Select relevant tag numbers
these.tags=unique(Tag.Deployed$ATAG.NO)
Tagging.Data=subset(Tagging.Data,Tag%in%these.tags)

#6. Manipulate tagging data
             Tagging.Data$"Recovery_latitude"=ifelse(Tagging.Data$"Recovery_latitude">0,-Tagging.Data$"Recovery_latitude",Tagging.Data$"Recovery_latitude")
Tagging.Data=Tagging.Data[!(duplicated(paste(Tagging.Data$Receiver.Number,Tagging.Data$Date.time,Tagging.Data$Tag))),]  #remove duplicates

maxdate=strptime(max(Tagging.Data$Date.time), "%Y-%m-%d")    #min and max dates of data
mindate=strptime(min(Tagging.Data$Date.time), "%Y-%m-%d")

Tagging.Data=Tagging.Data[order(Tagging.Data$Date.time),]     #order by date time




#7. Expand hits to have all days within studied period
YEARS=unique(1900 + as.POSIXlt(Tagging.Data$Date.time)$year)
year=NULL
getDays <- function(year)
{
  seq(as.Date(paste(year, "-01-01", sep="")), as.Date(paste(year, "-12-31", sep="")), by="+1 day")
}
alldays=lapply(YEARS,getDays)
DATES=alldays[[1]][1]
for (i in 1:length(alldays))
{
  DATES=c(DATES,alldays[[i]])
}
DATES=DATES[-1]

matchDates=unique(match(Tagging.Data$Date,DATES))
DATES=DATES[-matchDates]
DATES=subset(DATES,DATES>=as.Date(mindate))
DATES=subset(DATES,DATES<=as.Date(maxdate))

times <- rep("00:00:01",length(DATES))
x <- paste(DATES, times)
Dateshours=strptime(x, "%Y-%m-%d %H:%M:%S")


dummyDATES=data.frame(Receiver.Number=NA,Date.time=Dateshours,Date=DATES,Tag="       ",DS=NA,"Recovery_latitude"=NA,
"Recovery_longitude"=NA,"Comments"=NA)
names(dummyDATES)=names(Tagging.Data)


Tagging.Data=rbind(Tagging.Data,dummyDATES)
Tagging.Data=Tagging.Data[order(Tagging.Data$Date.time),]   #order by date and time


#8. Create useful vectors  and data frames
Sharks=unique(Tagging.Data$Tag)     #unique sharks
Sharks=Sharks[-match("       ",Sharks)]
n.shark=length(Sharks)

Shark.Times=unique(Tagging.Data$Date.time)     #unique date and time

Shark.colors=rainbow(n.shark)       #colors for these sharks if detected simultaneously
Shark.colors=sample(Shark.colors, n.shark, replace = FALSE)       #make sequential colors different


#Tagging.Data.byShark=list()
#for (i in 1:n.shark)
#{
#    datos=subset(Tagging.Data,Tag==Sharks[i])
#    datos.Deploy=subset(Tag.Deployed,ATAG.NO==Sharks[i])
#
#    mindate.byShark=strptime(datos.Deploy$DateDeployed, "%Y-%m-%d")
#
#    maxdate.byShark=strptime(max(datos$Date.time), "%Y-%m-%d")    #min and max dates of data
#
#        #fill in missing dates in between hits
#    YEARS.byShark=unique(1900 + as.POSIXlt(datos$Date.time)$year)
#
#    year.byShark=NULL
#    getDays.byShark <- function(year.byShark)
#    {
#      seq(as.Date(paste(year.byShark, "-01-01", sep="")), as.Date(paste(year.byShark, "-12-31", sep="")), by="+1 day")
#    }
#    alldays.byShark=lapply(YEARS.byShark,getDays.byShark)
#    DATES.byShark=alldays.byShark[[1]][1]
#    for (j in 1:length(alldays.byShark))
#    {
#      DATES.byShark=c(DATES.byShark,alldays.byShark[[j]])
#    }
#    DATES.byShark=DATES.byShark[!(duplicated(DATES.byShark))]
#    datos$Date=as.Date(datos$Date)
#
#    matchDates.byShark=unique(match(datos$Date,DATES.byShark))
#    DATES.byShark=DATES.byShark[-matchDates.byShark]
#
#    if(as.Date(mindate.byShark)<as.Date(maxdate.byShark))
#    { DATES.byShark=subset(DATES.byShark,DATES.byShark>=as.Date(mindate.byShark))
#      DATES.byShark=subset(DATES.byShark,DATES.byShark<=as.Date(maxdate.byShark))
#      times <- rep("00:00:01",length(DATES.byShark))
#      x <- paste(DATES.byShark, times)
#      Dateshours.byShark=strptime(x, "%Y-%m-%d %H:%M:%S")
#
#      dummyDATES=data.frame(Receiver.Number=NA,Date.time=Dateshours,Date=DATES,Tag="       ",DS=NA,"Recovery_latitude"=NA,
#      "Recovery_longitude"=NA,"Comments"=NA)
#      names(dummyDATES)=names(Tagging.Data)
#
#      datos=rbind(datos,dummyDATES)
#    }
#
#    #if(as.Date(mindate.byShark)==as.Date(maxdate.byShark)) DATES.byShark=as.Date(maxdate.byShark)
#
#    datos=datos[order(datos$Date.time),]   #order by date and time
#
#    Tagging.Data.byShark[[i]]=datos
#
#}

  #add species name
these.Tag.Deployed=match(c("ATAG.NO","Species"),names(Tag.Deployed))
Tagging.Data=merge(Tagging.Data,Tag.Deployed[,these.Tag.Deployed],by.x="Tag", by.y="ATAG.NO",all.x=T)




  #create data frames by species
Tagging.Data.SandBar=subset(Tagging.Data,Species=="TK")
Tagging.Data.Dusky=subset(Tagging.Data,Species=="BW")
#Tagging.Data.Gummy=subset(Tagging.Data,Species=="xxxx")                            #UPDATE with gummy and whiskery code
#Tagging.Data.Whiskery=subset(Tagging.Data,Species=="xxxx")

    #extract hits by day, shark, receiver and data set
hits.day=function(whichspecies)
{
  whichspecies$Hits=1
  Tagging.Data.byday=aggregate(Hits ~ Date + DS + Receiver.Number + Tag + Recovery_latitude + Recovery_longitude, data = whichspecies, sum)
  maxdate=strptime(max(Tagging.Data.byday$Date), "%Y-%m-%d")    #min and max dates of data
  mindate=strptime(min(Tagging.Data.byday$Date), "%Y-%m-%d")

  #add continuous dates
  DATES=alldays[[1]][1]
  for (i in 1:length(alldays))
  {
    DATES=c(DATES,alldays[[i]])
  }
  DATES=DATES[-1]
  matchDates=unique(match(Tagging.Data.byday$Date,DATES))

  DATES=DATES[-matchDates]
  DATES=subset(DATES,DATES>=as.Date(mindate))
  DATES=subset(DATES,DATES<=as.Date(maxdate))

  dummyDATES=data.frame(Date=DATES,DS=NA,Receiver.Number=NA,Tag="       ","Recovery_latitude"=NA,
  "Recovery_longitude"=NA,"Hits"=NA)
  names(dummyDATES)=names(Tagging.Data.byday)

  Tagging.Data.byday=rbind(Tagging.Data.byday,dummyDATES)

  Tagging.Data.byday=Tagging.Data.byday[order(Tagging.Data.byday$Date),]
  
  return(Tagging.Data.byday)

}

Tagging.Data.byday.SandBar=hits.day(Tagging.Data.SandBar)
Tagging.Data.byday.Dusky=hits.day(Tagging.Data.Dusky)
#Tagging.Data.byday.Gummy=hits.day(Tagging.Data.Gummy)                                                     #UPDATE
#Tagging.Data.byday.Whiskery=hits.day(Tagging.Data.Whiskery)



#9. Mapping

#9.1 create Oz
MyMap <- GetMap(center=c(-29,134), zoom=zoom,destfile = "Oz.png",maptype = c("satellite"))
insetOz <- function()
{
  opar <- par(mai = c(.5,.5,.5,.5))
  on.exit(par(opar))
  par(mar=rep(.2,4),xaxt="n",yaxt="n",plt=par("plt"))

  	# add text and polygons
  TextOnStaticMap(MyMap, lat=-24,lon=130, "Australia", cex=1.1, col = "white")
  edgeX=c(112,127,127,112)
  edgeY=c(-12,-12,-36,-36)
  polys <- data.frame(PID=rep(1,4),POS=1:4,X=edgeX,Y=edgeY)
  polys <- as.PolySet(polys, projection=2)
  PlotPolysOnStaticMap(MyMap, polys, col="transparent", border = "orange", lwd = 2, verbose = 0, add=TRUE)
}


#9.2  Plot of WA with Oz inset
    # download WA
MyMap.wa <- GetMap(center=center.wa, zoom=zoom.wa,destfile = "WA.png",maptype = c("satellite"),frame=T)
Latcoor1= LatLon2XY.centered(MyMap.wa,-15,124)
Latcoor2= LatLon2XY.centered(MyMap.wa,-25,124)
Latcoor3= LatLon2XY.centered(MyMap.wa,-35,124)
Longcoor1= LatLon2XY.centered(MyMap.wa,-36.7,100)
Longcoor2= LatLon2XY.centered(MyMap.wa,-36.7,110)
Longcoor3= LatLon2XY.centered(MyMap.wa,-36.7,120)
Latvector=c("15ºS","25ºS","35ºS")
Longvector=c("100ºE","110ºE","120ºE")


insetWA <- function()
{
  opar <- par(mai = c(.5,.5,.5,.5))
  on.exit(par(opar))
  par(mar=rep(.2,4),xaxt="n",yaxt="n",plt=par("plt"))

  TextOnStaticMap(MyMap.wa, lat=-25,lon=121.25, "Western", cex=2, col = "cyan")
  TextOnStaticMap(MyMap.wa, lat=-26.25,lon=121.25, "Australia", cex=2, col = "cyan", add=TRUE)
  	# add Ningaloo polygon
  edgeX=c(113,115,115,113)
  edgeY=c(-21,-21,-23,-23)
  polys <- data.frame(PID=rep(1,4),POS=1:4,X=edgeX,Y=edgeY)
  polys <- as.PolySet(polys, projection=2)
  PlotPolysOnStaticMap(MyMap.wa, polys, col="transparent", border = "orange", lwd = 2, verbose = 0, add=TRUE)
  TextOnStaticMap(MyMap.wa, lat=-22,lon=117.5, "Ningaloo", cex=1.5, col = "cyan", add=TRUE)
  
  	# add Perth polygon
  edgeX=c(115,116,116,115)
  edgeY=c(-31.5,-31.5,-32.5,-32.5)
  polys <- data.frame(PID=rep(1,4),POS=1:4,X=edgeX,Y=edgeY)
  polys <- as.PolySet(polys, projection=2)
  PlotPolysOnStaticMap(MyMap.wa, polys, col="transparent", border = "orange", lwd = 2, verbose = 0, add=TRUE)
  TextOnStaticMap(MyMap.wa, lat=-32,lon=117.5, "Perth", cex=1.5, col = "cyan", add=TRUE)
  
  	# add SouthWest polygon
  edgeX=c(116,120,120,116)
  edgeY=c(-34,-34,-36,-36)
  polys <- data.frame(PID=rep(1,4),POS=1:4,X=edgeX,Y=edgeY)
  polys <- as.PolySet(polys, projection=2)
  PlotPolysOnStaticMap(MyMap.wa, polys, col="transparent", border = "orange", lwd = 2, verbose = 0, add=TRUE)
  TextOnStaticMap(MyMap.wa, lat=-34.85,lon=123, "South West", cex=1.5, col = "cyan", add=TRUE)

  # add lat and long
  legend(Latcoor1[1], Latcoor1[2],Latvector[1],bty="n",text.col="white",cex=1)
  legend(Latcoor2[1], Latcoor2[2],Latvector[2],bty="n",text.col="white",cex=1)
  legend(Latcoor3[1], Latcoor3[2],Latvector[3],bty="n",text.col="white",cex=1)
  legend(Longcoor1[1], Longcoor1[2],Longvector[1],bty="n",text.col="white",cex=1)
  legend(Longcoor2[1], Longcoor2[2],Longvector[2],bty="n",text.col="white",cex=1)
  legend(Longcoor3[1], Longcoor3[2],Longvector[3],bty="n",text.col="white",cex=1)
}

png(file="WA.with.OZ.png",width=700,height=700)
insetWA()
    # add Oz
vp <- baseViewports()
pushViewport(vp$inner,vp$figure,vp$plot)
pushViewport(viewport(x=0.8,y=1,width=.2,height=.2,just=c("center","center")))
par(fig=gridFIG(),new=T)
insetOz()
dev.off()


Colors.WA=c('yellow','orange','white')

#9.3  Animate movement in whole Western Australia by species

#################
Legcoor.date <- LatLon2XY.centered(MyMap.wa,-15,118)
Legcoor.hits=  LatLon2XY.centered(MyMap.wa,-19,117.35)
Legcoor.shark <- LatLon2XY.centered(MyMap.wa,-22.5,118)
LatcoorScaleB <- LatLon2XY.centered(MyMap.wa,-13,118)

WesternAustralia.movie=function(whichshark)
{
  datos=whichshark

  uniquedates=unique(datos$Date)

  ani.options(interval=0.2,nmax=length(uniquedates),outdir = getwd())
  LABEL=unique(datos$Tag)

  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date
  {
    datos1=subset(datos,Date==uniquedates[i])
    HITS=sum(datos1$Hits)
    PlotOnStaticMap(MyMap.wa, c(Ningaloo$"Recovery_latitude",Perth$Lat,SouthWest$Lat),
    c(Ningaloo$"Recovery_longitude",Perth$Long,SouthWest$Long),cex=1,pch=19,col=Colors.WA[1])
    ConvCoord <- LatLon2XY.centered(MyMap.wa,datos1$Recovery_latitude,datos1$Recovery_longitud)
    points(ConvCoord[[1]],ConvCoord[[2]],cex=2.5,col=2,lwd=2)
    legend(Legcoor.date[1],Legcoor.date[2],unique(datos1$Date),title="Date",title.adj=0.5,
    title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5,x.intersp=-0.5)
    if(is.na(datos$Receiver.Number[i]))
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],unique(datos1$Tag),title="ID code",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

      legend(Legcoor.hits[1],Legcoor.hits[2],HITS,title="Number of detections",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.3)

    }else
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],unique(datos1$Tag),title="ID code",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

      legend(Legcoor.hits[1],Legcoor.hits[2],HITS,title="Number of detections",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.3)

    }

      # add scale bar
    PlotOnStaticMap(MyMap.wa, c(-14.2,-14.2), c(118.,118.8999),cex=1,lwd=4,col="white",FUN=lines,add=T)
    legend(LatcoorScaleB[1], LatcoorScaleB[2],"100 km",bty="n",text.col="white",cex=1.5)

    ani.pause()         #pause between hits for the interval length of time
  }

}

Start.syst=date()
saveGIF({WesternAustralia.movie(Tagging.Data.byday.SandBar)}, interval=0.5,movie.name="WesternAustralia.Sandbar.movie.wmv",
 ani.width=600,ani.height=600)
saveGIF({WesternAustralia.movie(Tagging.Data.byday.Dusky)}, interval=0.5,movie.name="WesternAustralia.Dusky.movie.wmv",
 ani.width=600,ani.height=600)
#saveGIF({WesternAustralia.movie(Tagging.Data.byday.Gummy)}, interval=0.5,movie.name="WesternAustralia.movie.wmv",
# ani.width=600,ani.height=600)                                                                                #UPDATE
#saveGIF({WesternAustralia.movie(Tagging.Data.byday.Whiskery)}, interval=0.5,movie.name="WesternAustralia.movie.wmv",
# ani.width=600,ani.height=600)

End.syst=date()


###################

#9.4  Animate movement in Ningaloo by species
MyMap.Ningaloo <- GetMap(center=c(-22.5,114), zoom=9,destfile = "Ningaloo.png",maptype = c("satellite"),frame=T)
Legcoor.date <- LatLon2XY.centered(MyMap.Ningaloo,-22,114.5)
Legcoor.shark <- LatLon2XY.centered(MyMap.Ningaloo,-22.5,114.5)
Legcoor.hits=  LatLon2XY.centered(MyMap.Ningaloo,-22.25,114.35)
LatcoorScaleB <- LatLon2XY.centered(MyMap.Ningaloo,-23,114.5)

Ningaloo.movie=function(whichshark)
{
  datos=subset(whichshark,!(DS%in%c("Perth","SW"))& Date<=max(Tagging.Data.Ningaloo$Date))

  uniquedates=unique(datos$Date)

  ani.options(interval=0.2,nmax=length(uniquedates),outdir = getwd())
  LABEL=unique(datos$Tag)

  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date
  {
    datos1=subset(datos,Date==uniquedates[i])
    HITS=sum(datos1$Hits)
    PlotOnStaticMap(MyMap.Ningaloo, Ningaloo$"Recovery_latitude",Ningaloo$"Recovery_longitude",cex=1,pch=19,col=Colors.WA[1])
    ConvCoord <- LatLon2XY.centered(MyMap.Ningaloo,datos1$Recovery_latitude,datos1$Recovery_longitud)
    points(ConvCoord[[1]],ConvCoord[[2]],cex=2.5,col=2,lwd=2)
    legend(Legcoor.date[1],Legcoor.date[2],unique(datos1$Date),title="Date",title.adj=0.5,
    title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5,x.intersp=-0.5)
    if(is.na(datos$Receiver.Number[i]))
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],unique(datos1$Tag),title="ID code",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

      legend(Legcoor.hits[1],Legcoor.hits[2],HITS,title="Number of detections",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

    }else
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],unique(datos1$Tag),title="ID code",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

      legend(Legcoor.hits[1],Legcoor.hits[2],HITS,title="Number of detections",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

    }

      # add scale bar
    PlotOnStaticMap(MyMap.Ningaloo, c(-23.1,-23.1), c(114.60,114.693),cex=1,lwd=4,col="white",FUN=lines,add=T)
    legend(LatcoorScaleB[1], LatcoorScaleB[2],"10 km",bty="n",text.col="white",cex=1.5)

    ani.pause()         #pause between hits for the interval length of time
  }

}

Start.syst=date()
saveGIF({Ningaloo.movie(Tagging.Data.byday.SandBar)}, interval=0.5,movie.name="Ningaloo.Sandbar.movie.wmv",
 ani.width=600,ani.height=600)
saveGIF({Ningaloo.movie(Tagging.Data.byday.Dusky)}, interval=0.5,movie.name="Ningaloo.Dusky.movie.wmv",
 ani.width=600,ani.height=600)
#saveGIF({Ningaloo.movie(Tagging.Data.byday.Gummy)}, interval=0.5,movie.name="Ningaloo.movie.wmv",
# ani.width=600,ani.height=600)                                                                                #UPDATE
#saveGIF({Ningaloo.movie(Tagging.Data.byday.Whiskery)}, interval=0.5,movie.name="Ningaloo.movie.wmv",
# ani.width=600,ani.height=600)

End.syst=date()



#9.5  Animate movement in Perth by species
MyMap.Perth <- GetMap(center=c(-32.05,115.60), zoom=10,destfile = "Perth.png",maptype = c("satellite"),frame=T)
Legcoor.date <- LatLon2XY.centered(MyMap.Perth,-31.875,115.83)
Legcoor.shark <- LatLon2XY.centered(MyMap.Perth,-32.1,115.83)
Legcoor.hits=  LatLon2XY.centered(MyMap.Perth,-32,115.77)
LatcoorScaleB <- LatLon2XY.centered(MyMap.Perth,-32.30,115.8375)


Perth.movie=function(whichshark)
{
  datos=subset(whichshark,!(DS%in%c("Ningaloo","SW"))& Date<=max(Tagging.Data.Perth$Date)& Date>=min(Tagging.Data.Perth$Date))
  uniquedates=unique(datos$Date)

  ani.options(interval=0.2,nmax=length(uniquedates),outdir = getwd())
  LABEL=unique(datos$Tag)

  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date
  {
    datos1=subset(datos,Date==uniquedates[i])
    HITS=sum(datos1$Hits)
    PlotOnStaticMap(MyMap.Perth, Perth$Lat,Perth$Long,cex=1,pch=19,col=Colors.WA[1])
    ConvCoord <- LatLon2XY.centered(MyMap.Perth,datos1$Recovery_latitude,datos1$Recovery_longitud)
    points(ConvCoord[[1]],ConvCoord[[2]],cex=2.5,col=2,lwd=2)
    legend(Legcoor.date[1],Legcoor.date[2],unique(datos1$Date),title="Date",title.adj=0.5,
    title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5,x.intersp=-0.5)
    if(is.na(datos$Receiver.Number[i]))
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],unique(datos1$Tag),title="ID code",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

      legend(Legcoor.hits[1],Legcoor.hits[2],HITS,title="Number of detections",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.25)

    }else
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],unique(datos1$Tag),title="ID code",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

      legend(Legcoor.hits[1],Legcoor.hits[2],HITS,title="Number of detections",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.25)

    }

      # add scale bar
    PlotOnStaticMap(MyMap.Perth, c(-32.35,-32.35), c(115.90,115.993),cex=1,lwd=4,col="white",FUN=lines,add=T)
    legend(LatcoorScaleB[1], LatcoorScaleB[2],"10 km",bty="n",text.col="white",cex=1.5)

    ani.pause()         #pause between hits for the interval length of time
  }

}

Start.syst=date()
saveGIF({Perth.movie(Tagging.Data.byday.SandBar)}, interval=0.5,movie.name="Perth.Sandbar.movie.wmv",
 ani.width=600,ani.height=600)
saveGIF({Perth.movie(Tagging.Data.byday.Dusky)}, interval=0.5,movie.name="Perth.Dusky.movie.wmv",
 ani.width=600,ani.height=600)
#saveGIF({Perth.movie(Tagging.Data.byday.Gummy)}, interval=0.5,movie.name="Perth.movie.wmv",
# ani.width=600,ani.height=600)                                                                                #UPDATE
#saveGIF({Perth.movie(Tagging.Data.byday.Whiskery)}, interval=0.5,movie.name="Perth.movie.wmv",
# ani.width=600,ani.height=600)

End.syst=date()




#9.6  Animate movement in South West  by species
MyMap.SW <- GetMap(center=c(-34.5,118), zoom=8,destfile = "SW.png",maptype = c("satellite"),frame=T)

Legcoor.date <- LatLon2XY.centered(MyMap.SW,-33.1,116.4)
Legcoor.hits=  LatLon2XY.centered(MyMap.SW,-33.1,117.3)
Legcoor.shark <- LatLon2XY.centered(MyMap.SW,-33.1,118.4)
LatcoorScaleB <- LatLon2XY.centered(MyMap.SW,-33.1,118.9)


SW.movie=function(whichshark)
{
  datos=subset(whichshark,!(DS%in%c("Ningaloo","Perth"))& Date<=max(Tagging.Data.SouthWest$Date)& Date>=min(Tagging.Data.SouthWest$Date))
  uniquedates=unique(datos$Date)

  ani.options(interval=0.2,nmax=length(uniquedates),outdir = getwd())
  LABEL=unique(datos$Tag)

  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date
  {
    datos1=subset(datos,Date==uniquedates[i])
    HITS=sum(datos1$Hits)
    PlotOnStaticMap(MyMap.SW, SouthWest$Lat, SouthWest$Long,cex=1,pch=19,col=Colors.WA[1])
    ConvCoord <- LatLon2XY.centered(MyMap.SW,datos1$Recovery_latitude,datos1$Recovery_longitud)
    points(ConvCoord[[1]],ConvCoord[[2]],cex=2.5,col=2,lwd=2)
    legend(Legcoor.date[1],Legcoor.date[2],unique(datos1$Date),title="Date",title.adj=0.5,
    title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5,x.intersp=-0.5)
    if(is.na(datos$Receiver.Number[i]))
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],unique(datos1$Tag),title="ID code",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

      legend(Legcoor.hits[1],Legcoor.hits[2],HITS,title="Number of detections",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.25)

    }else
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],unique(datos1$Tag),title="ID code",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.5)

      legend(Legcoor.hits[1],Legcoor.hits[2],HITS,title="Number of detections",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.25)

    }

      # add scale bar
    PlotOnStaticMap(MyMap.SW, c(-33.35,-33.35), c(119.40,119.493),cex=1,lwd=4,col="white",FUN=lines,add=T)
    legend(LatcoorScaleB[1], LatcoorScaleB[2],"10 km",bty="n",text.col="white",cex=1.5)

    ani.pause()         #pause between hits for the interval length of time
  }

}

Start.syst=date()
saveGIF({SW.movie(Tagging.Data.byday.SandBar)}, interval=0.5,movie.name="SW.Sandbar.movie.wmv",
 ani.width=600,ani.height=600)
saveGIF({SW.movie(Tagging.Data.byday.Dusky)}, interval=0.5,movie.name="SW.Dusky.movie.wmv",
 ani.width=600,ani.height=600)
#saveGIF({SW.movie(Tagging.Data.byday.Gummy)}, interval=0.5,movie.name="SW.movie.wmv",
# ani.width=600,ani.height=600)                                                                                #UPDATE
#saveGIF({SW.movie(Tagging.Data.byday.Whiskery)}, interval=0.5,movie.name="SW.movie.wmv",
# ani.width=600,ani.height=600)

End.syst=date()
