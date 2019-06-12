#GOOGLE MAPPING OF VR4S AND WHITE POINTERS

#notes: this script maps the location of VR4s and hits received from tagged white pointers
#       remember to use latest version of input files


setwd("H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Mapping/WhitePointer_Vr4")

library(RgoogleMaps)  #library for plotting on google maps
library(RODBC)				#library for importing excel data
library(gridBase)			#for inset map
library(PBSmapping)   #for polygon
library(animation)    #for animations



windows(record=T)			#to keep al graphs

memory.limit(3900)   #set memory limit to maximum size


#DATA SECTION

  #1. load data
    #1.1. listening stations
# Metro area
channel1 <- odbcConnectExcel("M:/Fisheries Research/FinFish/Shark/Braccini/Shark MonitoringOTN/SMN_OTN Running Sheet 2011 copy")
SMN_VR2_2008<- sqlFetch(channel1,"SMN VR2W 2008")
SMN_VR2_2009<- sqlFetch(channel1,"SMN VR2W 2009")
SMN_VR2_2011<- sqlFetch(channel1,"SMN VR2W 2011")
SMN_VR2_2010<- sqlFetch(channel1,"SMN VR2W 2010")
SMN_VR4G_2011<- sqlFetch(channel1,"VR4G 2011")
SMN_VR4G_2009<- sqlFetch(channel1,"VR4G Deployed 2009")
SMN_VR4G_2009$Long[2]=115.73338       #fill in R import f... up
SMN_VR4G_2009$Long[17]=115.47654
SMN_VR4G_planned<- sqlFetch(channel1,"VR4G planned")
OTN_2009a<- sqlFetch(channel1,"OTN 2009a")
OTN_2009b<- sqlFetch(channel1,"OTN 2009b")
close(channel1)


    #1.2. shark detections
# Metro area
channel2 <- odbcConnectExcel("M:/Fisheries Research/FinFish/Shark/Braccini/Shark MonitoringOTN/2009 detections")
ALL_detections<- sqlFetch(channel2,"SMN")
close(channel2)

    #1.3. tag deployment date
Tag.Deployed=data.frame(Date=as.Date(c("2009-07-04","2009-05-24","2009-05-24","2009-07-14","2009-07-14")),
ID=c("58062","8500","8501","8510","8555"),Sex=c("male","male","female","unknown","unknown"),
Length=c(4.6,3.5,4.5,4,4),Lat=c(-35.211,-31.825,-31.825,-31.851,-31.851),
Long=c(136.833,115.69,115.69,115.587,115.587))

  #2. specify mapping area
    #2.1. Australia
minLat=-38; maxLat=-20; minLong=117; maxLong=147;
lat=c(maxLat,minLat);long=c(minLong,maxLong)
center=c(mean(lat),mean(long))
zoom=min(MaxZoom(range(lat),range(long)))

    #2.2. WA
lat.wa=c(-15,-45);long.wa=c(110,120)
#center.wa=c(mean(lat.wa),mean(long.wa))
center.wa=c(-26,119)
#zoom.wa=max(MaxZoom(range(lat.wa),range(long.wa)))
zoom.wa=5
 
    #2.3. Perth metro
lat.metro=c(-29,-32);long.metro=c(114,115)
#center.metro=c(mean(lat.metro),mean(long.metro))
center.metro=c(-32.05,115.60)
#zoom.metro=max(MaxZoom(range(lat.metro),range(long.metro)))
zoom.metro=10





#PROCEDURE SECTION

#1. Convert UTC detection time to Perth local time
#not needed, done already in data set



#2. Create dataframes with location of stations for each year
SMN_VR2_2008$year=2008
SMN_VR2_2009$year=2009
SMN_VR2_2010$year=2010
SMN_VR2_2010$Lat=ifelse(SMN_VR2_2010$Lat>0,-SMN_VR2_2010$Lat,SMN_VR2_2010$Lat) #correct missing -
SMN_VR2_2011$year=2011
SMN_VR2_2011$Lat=ifelse(SMN_VR2_2011$Lat>0,-SMN_VR2_2011$Lat,SMN_VR2_2011$Lat) #correct missing -
SMN_VR4G_2011$year=2011
SMN_VR4G_2009$year=2009
OTN_2009a$year=2009
OTN_2009b$year=2009

thesecolumns=match(c("Station No","Lat","Long","Serial No","year"),names(SMN_VR2_2008))
thesecolumns11=match(c("Station No","Lat","Long","Serial No","year"),names(SMN_VR2_2011))

SMN_VR2=list()
SMN_VR2[['2008']]=SMN_VR2_2008[,thesecolumns]
SMN_VR2[['2009']]=SMN_VR2_2009[,thesecolumns]
SMN_VR2[['2010']]=SMN_VR2_2010[,thesecolumns]
SMN_VR2[['2011']]=SMN_VR2_2011[,thesecolumns11]


thesecolumnsVR4=match(c("Station No","Lat","Long","VR4G Serial No","year"),names(SMN_VR4G_2011))
thesecolumnsOTN09a=match(c("Station No","Lat","Long","Serial No","year"),names(OTN_2009a))
thesecolumnsOTN09b=match(c("Station No","Lat","Long","Serial No","year"),names(OTN_2009b))

SMN_VR4G=list()
SMN_VR4G[['2009']]=SMN_VR4G_2009[,thesecolumnsVR4]
SMN_VR4G[['2011']]=SMN_VR4G_2011[,thesecolumnsVR4]

OTN_VR2=list()
OTN_VR2[['2009a']]=OTN_2009a[,thesecolumnsOTN09a]
OTN_VR2[['2009b']]=OTN_2009b[,thesecolumnsOTN09b]

# extract lat and long of receivers for plotting
LatLongVR2=SMN_VR2_2009[,match(c("Lat","Long"),names(SMN_VR2_2009))]
LatLongVR4=SMN_VR4G_2009[,match(c("Lat","Long"),names(SMN_VR4G_2009))]
LatLongOTN=OTN_2009a[,match(c("Lat","Long"),names(OTN_2009a))]



#3. Create dataframes with location of hits received from tagged white pointers
theseColsShkHits=c("Date","Year","Month","Day","Time","Event no#","Event ID","ID","Receiver S/N","Project",
"Lat","Long","Station","Location","Depth (m)")
 WhitePointers=ALL_detections[,match(theseColsShkHits,names(ALL_detections))]


WhitePointers=WhitePointers[!(duplicated(paste(WhitePointers$Time,WhitePointers$ID,WhitePointers$Lat))),]  #remove duplicates

maxdate=strptime(max(WhitePointers$Time), "%Y-%m-%d")    #min and max dates of data
mindate=strptime(min(WhitePointers$Time), "%Y-%m-%d")

WhitePointers=WhitePointers[order(WhitePointers$Date),]

WhitePointers$Date=strptime(WhitePointers$Time, "%Y-%m-%d")      #correct UTC date


#4. Expand hits to have all days within studied period
YEARS=unique(WhitePointers$Year)
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


DATES=as.character(DATES)                                #convert to character for merging
WhitePointers$Date=as.character(WhitePointers$Date)

matchDates=unique(match(WhitePointers$Date,DATES))
DATES=DATES[-matchDates]
DATES=subset(DATES,DATES>=mindate)
DATES=subset(DATES,DATES<=maxdate)

times <- rep("00:00:01",length(DATES))
x <- paste(DATES, times)
Dateshours=strptime(x, "%Y-%m-%d %H:%M:%S")


dummyDATES=data.frame(date=DATES,Year=NA,Month=NA,Day=NA,Time=Dateshours,"Event no#"=NA,"Event ID"=NA,ID="       ","Receiver S/N"=NA,
Project="OTN",Lat=40,Long=40,Station=NA,Location=NA,"Depth (m)"=NA)
names(dummyDATES)=names(WhitePointers)


WhitePointers=rbind(WhitePointers,dummyDATES)
WhitePointers=WhitePointers[order(WhitePointers$Time),]   #order by date and time


#5. Create useful vectors
Sharks=unique(WhitePointers$ID)     #unique sharks
Sharks=Sharks[-match("       ",Sharks)]
n.shark=length(Sharks)

Shark.Times=unique(WhitePointers$Time)     #unique date and time

Shark.colors=rainbow(n.shark)       #colors for these sharks if detected simultaneously
Shark.colors=sample(Shark.colors, n.shark, replace = FALSE)       #make sequential colors different



#6. Create dataframes with location of hits received by VR4s only
  #VR4 hits only
WhitePointers.VR4s=subset(WhitePointers,Lat%in%na.omit(LatLongVR4$Lat))

  #List of hits by Shark
#LastTime="2011-09-30 00:00:00 WST"
#maxdate.byShark=strptime(LastTime, "%Y-%m-%d")
    
WhitePointers.byShark=list()
for (i in 1:n.shark)
{
    datos=subset(WhitePointers,ID==Sharks[i])
#    datos$Date=strptime(datos$Time, "%Y-%m-%d")      #correct UTC date
    
    datos.Deploy=subset(Tag.Deployed,ID==Sharks[i])

    mindate.byShark=strptime(datos.Deploy$Date, "%Y-%m-%d")
    #mindate.byShark=strptime(min(datos$Time), "%Y-%m-%d")
    maxdate.byShark=strptime(max(datos$Time), "%Y-%m-%d")    #min and max dates of data

        #fill in missing dates in between hits
    YEARS.byShark=unique(datos$Year)
 #   First.day=as.POSIXlt(datos$Time[1])$mday       #extract days and months
 #   First.month=1+as.POSIXlt(datos$Time[1])$mon
 #   Last.day=as.POSIXlt(datos$Time[length(datos$Time)])$mday
 #   Last.month=1+as.POSIXlt(datos$Time[length(datos$Time)])$mon
 #    pastedFirst=paste("-",First.month,"-",First.day,sep="")
 #   pastedLast=paste("-",Last.month,"-",Last.day,sep="")
    
    year.byShark=NULL
    getDays.byShark <- function(year.byShark)
    {
      #seq(as.Date(paste(year.byShark,pastedFirst, sep="")), as.Date(paste(year.byShark,pastedLast, sep="")), by="+1 day")
      seq(as.Date(paste(year.byShark, "-01-01", sep="")), as.Date(paste(year.byShark, "-12-31", sep="")), by="+1 day")
    }
    alldays.byShark=lapply(YEARS.byShark,getDays.byShark)
    DATES.byShark=alldays.byShark[[1]][1]
    for (j in 1:length(alldays.byShark))
    {
      DATES.byShark=c(DATES.byShark,alldays.byShark[[j]])
    }
    DATES.byShark=DATES.byShark[!(duplicated(DATES.byShark))]
    DATES.byShark=as.character(DATES.byShark)                                #convert to character for merging
    datos$Date=as.character(datos$Date)

    matchDates.byShark=unique(match(datos$Date,DATES.byShark))
    DATES.byShark=DATES.byShark[-matchDates.byShark]

    DATES.byShark=subset(DATES.byShark,DATES.byShark>=mindate.byShark)
    DATES.byShark=subset(DATES.byShark,DATES.byShark<=maxdate.byShark)

    times <- rep("00:00:01",length(DATES.byShark))
    x <- paste(DATES.byShark, times)
    Dateshours.byShark=strptime(x, "%Y-%m-%d %H:%M:%S")

    dummyDATES=data.frame(date=DATES.byShark,Year=NA,Month=NA,Day=NA,Time=Dateshours.byShark,"Event no#"=NA,
    "Event ID"=NA,ID=unique(datos$ID),"Receiver S/N"=NA,Project="OTN",Lat=40,Long=40,Station=NA,Location=NA,"Depth (m)"=NA)
    names(dummyDATES)=names(datos)

    datos=rbind(datos,dummyDATES)
    datos=datos[order(datos$Time),]   #order by date and time

    WhitePointers.byShark[[i]]=datos

}





#7. Mapping

#7.1 create Oz
MyMap <- GetMap(center=center, zoom=zoom,destfile = "Oz.png",maptype = c("satellite"))
whatmap=NULL
insetOz <- function(whatmap)
{
  opar <- par(mai = c(.5,.5,.5,.5))
  on.exit(par(opar))
  par(mar=rep(.2,4),xaxt="n",yaxt="n",plt=par("plt"))

  	# add polygons
  if(whatmap=="WA")
  {
  TextOnStaticMap(MyMap, lat=-24,lon=130, "Australia", cex=1.1, col = "cyan")
  edgeX=c(112,128,128,112)
  edgeY=c(-13,-13,-36,-36)
  polys <- data.frame(PID=rep(1,4),POS=1:4,X=edgeX,Y=edgeY)
  polys <- as.PolySet(polys, projection=2)
  PlotPolysOnStaticMap(MyMap, polys, col="transparent", border = "orange", lwd = 2, verbose = 0, add=TRUE)
  }

  if(whatmap=="metro")
  {
  TextOnStaticMap(MyMap, lat=-24,lon=135, "Australia", cex=2.2, col = "cyan")
  edgeX=c(113.5,120,120,113.5)
  edgeY=c(-29,-29,-36,-36)
  polys <- data.frame(PID=rep(1,4),POS=1:4,X=edgeX,Y=edgeY)
  polys <- as.PolySet(polys, projection=2)
  PlotPolysOnStaticMap(MyMap, polys, col="transparent", border = "orange", lwd = 2.75, verbose = 0, add=TRUE)
  PlotOnStaticMap(MyMap,lat=-32,lon=116.2,cex=1.75,pch=20,col='cyan', add=TRUE)
  polysOut=data.frame(PID=rep(1,4),POS=1:4,X=c(104,160,160,104),Y=c(1,1,-50.5,-50.5))
  polysOut <- as.PolySet(polysOut, projection=2)
  PlotPolysOnStaticMap(MyMap, polysOut, col="transparent", border = "orange", lwd = 2.75, verbose = 0, add=TRUE)
  }
}


#7.2 Plot of Perth with Oz inset
    # download Perth
MyMap.metro <- GetMap(center=center.metro, zoom=zoom.metro,destfile = "Perth.png",maptype = c("satellite"),frame=T)

#png("Perth.png",640,640,res=100);

		# control width and height of graphic
#X11(width=14,height=12)
#par(mfcol=c(1,1),las=1,mai=c(.3, .3, .6, 0.6),omi=c(.3,.3,.6,0.6))

    # plot Perth
      #plot stations
Colors.metro=c('yellow','orange','white')
Legends.metro=c("VR2_SMN","VR4_SMN","VR2_OTN")
Legcoor <- LatLon2XY.centered(MyMap.metro,-32.30,115.2123)   #convert lat and longs to google coordinates

LatcoorScaleB <- LatLon2XY.centered(MyMap.metro,-32.370,115.879)

Latcoor1= LatLon2XY.centered(MyMap.metro,-31.75,115.92)
Latcoor2= LatLon2XY.centered(MyMap.metro,-32.00,115.92)
Latcoor3= LatLon2XY.centered(MyMap.metro,-32.25,115.92)
Longcoor1= LatLon2XY.centered(MyMap.metro,-31.675,115.20)
Longcoor2= LatLon2XY.centered(MyMap.metro,-31.675,115.45)
Longcoor3= LatLon2XY.centered(MyMap.metro,-31.675,115.70)

#degree=expression(K ~ group("[",degree,"]"))
Latvector=c("31º45'S","32º00'S","32º15'S")
Longvector=c("115º15'E","115º30'E","115º45'E")

metro=function(whatmap)
{
  PlotOnStaticMap(MyMap.metro, c(SMN_VR2[['2009']]$Lat,SMN_VR2[['2010']]$Lat),c(SMN_VR2[['2009']]$Long,SMN_VR2[['2010']]$Long),
  cex=1,pch=19,col=Colors.metro[1])
  PlotOnStaticMap(MyMap.metro, SMN_VR4G[['2009']]$Lat, SMN_VR4G[['2009']]$Long,cex=1,pch=17,col=Colors.metro[2],add=T)
  PlotOnStaticMap(MyMap.metro, OTN_VR2[['2009a']]$Lat, OTN_VR2[['2009a']]$Long,cex=1,pch=19,col=Colors.metro[3],add=T)

  if(whatmap=="metro")
  {
        #add text
    TextOnStaticMap(MyMap.metro, lat=-31.97,lon=115.8, "Perth", cex=1.5, col = "cyan",add=T)

      # add scale bar
    PlotOnStaticMap(MyMap.metro, c(-32.375,-32.375), c(115.90,115.993),cex=1,lwd=4,col="white",FUN=lines,add=T)
    legend(LatcoorScaleB[1], LatcoorScaleB[2],"10 km",bty="n",text.col="white",cex=1)

      # add lat and long
  legend(Latcoor1[1], Latcoor1[2],Latvector[1],bty="n",text.col="white",cex=1)
  legend(Latcoor2[1], Latcoor2[2],Latvector[2],bty="n",text.col="white",cex=1)
  legend(Latcoor3[1], Latcoor3[2],Latvector[3],bty="n",text.col="white",cex=1)
  legend(Longcoor1[1], Longcoor1[2],Longvector[1],bty="n",text.col="white",cex=1)
  legend(Longcoor2[1], Longcoor2[2],Longvector[2],bty="n",text.col="white",cex=1)
  legend(Longcoor3[1], Longcoor3[2],Longvector[3],bty="n",text.col="white",cex=1)
  }

      # add legends
  legend(Legcoor[1],Legcoor[2],Legends.metro,pch=c(19,17,19),bty="n",col=Colors.metro,text.col=Colors.metro,cex=1.1)
}



    # plot hits
whatshark=NULL
Interval=0.5    #interval between frames
Legcoor.shark <- LatLon2XY.centered(MyMap.metro,-32.01,115.8)   #convert lat and longs to google coordinates
Legcoor.date <- LatLon2XY.centered(MyMap.metro,-32.14,115.8)
Legcoor.hits <- LatLon2XY.centered(MyMap.metro,-32.25,115.8)


      # a)... Movie for displaying a subset of ID=8501 hits
Anim.movie.8501=function(whatshark)
{
  ani.options(interval=0.2,nmax=nrow(whatshark),outdir = getwd())
  LABEL=unique(whatshark$ID)
  LENGTH=paste(Tag.Deployed$Length[match(LABEL,Tag.Deployed$ID)],"m")
  SEX=as.character(Tag.Deployed$Sex[match(LABEL,Tag.Deployed$ID)])
  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date and time
  {
    metro("metro")
    ConvCoord <- LatLon2XY.centered(MyMap.metro,whatshark$Lat[i],whatshark$Long[i])
    points(ConvCoord[1],ConvCoord[2],cex=2,col=2,lwd=2)
    legend(Legcoor.date[1],Legcoor.date[2],whatshark$Time[i],title="Date",title.adj=0.5,
    title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1,x.intersp=-0.5)
    if(is.na(whatshark$Station[i])==T)
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],c(paste(whatshark$ID[i],"    ",sep=""),LENGTH,SEX),title="Shark",
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1)
    }else
    { legend(Legcoor.shark[1],Legcoor.shark[2],c(paste(whatshark$ID[i],"    ",sep=""),LENGTH,SEX),title="Shark",
    title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1)
    }
    ani.pause()         #pause between hits for the interval length of time
  }
}


 # save  as wmv
WhitePointer_8501=subset(WhitePointers.byShark[[2]],Date>="2009-07-23")
WhitePointer_8501=subset(WhitePointer_8501,Date<="2009-07-25")
Start.syst=date()
saveGIF({Anim.movie.8501(WhitePointer_8501)}, interval=Interval,
movie.name="Shark.8501subset.wmv", ani.width=600,ani.height=600)

End.syst=date()


      # b)... Movie of a summary of all sharks shown by shark
          # b.1) By blocks of days with hits and no hits

#datelist=list()
#for (i in 1:n.shark)
#{
#  datelist[[i]]=unique(WhitePointers.byShark[[i]]$Date)
#}


    #create list of ranges of hits and no hits for each shark
dates.ranges=list()
for (j in 1: n.shark)
{
  a=WhitePointers.byShark[[j]][,1:2]
  a=a[!duplicated(a$Date),]
  a$vector=with(a,ifelse(is.na(Year)==T,0,100))
  a$Cumsum=cumsum(a$vector)
  a$Cumsum=ifelse(a$vector==0,a$Cumsum+1,a$Cumsum)

  theseChunks=unique(a$Cumsum)

  lista=list()
  for(i in 1: length(theseChunks))
  {
    datos=subset(a$Date,a$Cumsum==theseChunks[i])
    if(length(datos)>1)
    {datos=range(datos)}else{datos=datos}
    lista[[i]]=datos
  }
  dates.ranges[[paste(Sharks[j])]]=lista
}


Anim.movie=function(whatshark,DatesRanges)
{
  Nmax=length(DatesRanges)
  whatshark$Date=as.Date(whatshark$Date)
  whatshark$ID=ifelse(whatshark$ID=="       ",NA,whatshark$ID)

  ani.options(interval=3,nmax=Nmax,outdir = getwd())
  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date and time
  {
    metro("metro")
    datos=subset(whatshark,Date>=min(DatesRanges[[i]]))
    datos=subset(datos,Date<=max(DatesRanges[[i]]))
    LABEL=unique(datos$ID)
    LENGTH=paste(Tag.Deployed$Length[match(LABEL,Tag.Deployed$ID)],"m")
    SEX=as.character(Tag.Deployed$Sex[match(LABEL,Tag.Deployed$ID)])
    hits=nrow(datos[is.na(datos$"Receiver S/N")==F,])
    ConvCoord <- LatLon2XY.centered(MyMap.metro,datos$Lat,datos$Long)
    points(ConvCoord[[1]],ConvCoord[[2]],cex=2,col=2,lwd=2)

    if(length(DatesRanges[[i]])>1)
    {
      legend(Legcoor.date[1],Legcoor.date[2],paste(DatesRanges[[i]]),title="Date",title.adj=0.5,
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1,x.intersp=-0.5)
     }else
     {
      legend(Legcoor.date[1],Legcoor.date[2],DatesRanges[[i]],title="Date",title.adj=0.5,
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1,x.intersp=2)
     }

    legend(Legcoor.shark[1],Legcoor.shark[2],c(paste(LABEL,"    ",sep=""),LENGTH,SEX),title="Shark",title.col="white",
    bg="black",box.col ="white",text.col="white",cex=1.1)

    legend(Legcoor.hits[1],Legcoor.hits[2],paste(hits,"    ",sep=""),title="# hits",title.col="white",
    bg="black",box.col ="white",text.col="white",cex=1.1)

    ani.pause()
  }
}

 # save as wmv
Start.syst=date()
for (k in 1: n.shark)
{
   saveGIF({Anim.movie(WhitePointers.byShark[[k]],dates.ranges[[k]])},
   movie.name=paste("Shark.",Sharks[k],".wmv",sep=''), ani.width=600,ani.height=600)
}
End.syst=date()



          # b.2) Varying speed of animation based on hit or no hit in a day
#note: not working yet, works fine in R but doesn't get through when exporting moving
#Anim.movie.byShark=function(whatshark)
#{
#  Nmax=length(unique(whatshark$Date))
#  Fecha=unique(whatshark$Date)

#  ani.options(interval=Interval,nmax=Nmax,outdir = getwd())
#  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date and time
#  {
#    metro("metro")
#    datos=subset(whatshark,Date==Fecha[i])
#    LABEL=unique(datos$ID)
#    hits=nrow(datos[is.na(datos$"Receiver S/N")==F,])

       #adjust frame speed for when there is no hits
#    if(hits==0)
#    {ani.options(interval=0.1)}else{ani.options(interval=1)}

#    ConvCoord <- LatLon2XY.centered(MyMap.metro,datos$Lat,datos$Long)
#    points(ConvCoord[[1]],ConvCoord[[2]],cex=2,col=2,lwd=2)

#    legend(Legcoor.shark[1],Legcoor.date[2],Fecha[i],title="Date",title.adj=0.5,
#    title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1,x.intersp=-0.5)

#    legend(Legcoor.shark[1],Legcoor.shark[2],paste(LABEL,"    ",sep=""),title="Shark",title.col="white",
#    bg="black",box.col ="white",text.col="white",cex=1.1)

#    legend(Legcoor.hits[1],Legcoor.hits[2],paste(hits,"    ",sep=""),title="# hits",title.col="white",
#    bg="black",box.col ="white",text.col="white",cex=1.1)

#    ani.pause()         #pause between hits for the interval length of time
#  }
#}

 # save as wmv
#Start.syst=date()
#for (i in 1: n.shark)
#{
#   saveGIF({Anim.movie.byShark(WhitePointers.byShark[[i]])},
#   movie.name=paste("Shark.",Sharks[i],".wmv",sep=''), ani.width=600,ani.height=600)
#}
#End.syst=date()




      # c)... Movie of a summary of all sharks together

    #create list of ranges of hits and no hits for each shark
All.dates.ranges=list()
# WhitePointers$Date=strptime(WhitePointers$Time, "%Y-%m-%d")      #correct UTC date
WhitePointers=WhitePointers[-1,]
a=WhitePointers[,1:2]
a=a[!duplicated(a$Date),]
a$vector=with(a,ifelse(is.na(Year)==T,0,100))
a$Cumsum=cumsum(a$vector)
a$Cumsum=ifelse(a$vector==0,a$Cumsum+1,a$Cumsum)
theseChunks=unique(a$Cumsum)

for(i in 1: length(theseChunks))
 {
   datos=subset(a$Date,a$Cumsum==theseChunks[i])
   if(length(datos)>1)
   {datos=range(datos)}else{datos=datos}
   All.dates.ranges[[i]]=as.Date(datos)
 }

 Legcoor.date.all= LatLon2XY.centered(MyMap.metro,-31.78,115.8)

Anim.movie.All=function(whatshark,DatesRanges)
{
  Nmax=length(DatesRanges)
  whatshark$Date=as.Date(whatshark$Date)

  ani.options(interval=3,nmax=Nmax,outdir = getwd())
  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date and time
  {
    metro("metro")
    datos=subset(whatshark,Date>=min(DatesRanges[[i]]))
    datos=subset(datos,Date<=max(DatesRanges[[i]]))
    LABEL=unique(datos$ID)
    LENGTH=paste(Tag.Deployed$Length[match(LABEL,Tag.Deployed$ID)],"m")
    SEX=as.character(Tag.Deployed$Sex[match(LABEL,Tag.Deployed$ID)])
    hits=nrow(datos[is.na(datos$"Receiver S/N")==F,])
    ConvCoord <- LatLon2XY.centered(MyMap.metro,datos$Lat,datos$Long)
    points(ConvCoord[[1]],ConvCoord[[2]],cex=2,col=2,lwd=2)

    if(length(DatesRanges[[i]])>1)
    {
      legend(Legcoor.date.all[1],Legcoor.date.all[2],paste(DatesRanges[[i]]),title="Date",title.adj=0.5,
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1,x.intersp=-0.5)
     }else
     {
      legend(Legcoor.date.all[1],Legcoor.date.all[2],DatesRanges[[i]],title="Date",title.adj=0.5,
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1,x.intersp=0)
     }
   if(is.na(SEX))
    {
      legend(Legcoor.shark[1],Legcoor.shark[2],c(paste(LABEL,"    ",sep=""),"   ","   "),title="Shark",title.col="white",
      bg="black",box.col ="white",text.col="white",cex=1.1)
      }else
     {
     legend(Legcoor.shark[1],Legcoor.shark[2],c(paste(LABEL,"    ",sep=""),LENGTH,SEX),title="Shark",title.col="white",
      bg="black",box.col ="white",text.col="white",cex=1.1)
     }
    legend(Legcoor.hits[1],Legcoor.hits[2],paste(hits,"    ",sep=""),title="# hits",title.col="white",
    bg="black",box.col ="white",text.col="white",cex=1.1)

    ani.pause()
  }
}

 # save as wmv
Start.syst=date()
saveGIF({Anim.movie.All(WhitePointers,All.dates.ranges)},
movie.name="All.Shark.wmv", ani.width=600,ani.height=600)
End.syst=date()






      # d)... Movie using all hits continuosly
Anim.movie.Cont=function(whatshark)
{
  ani.options(interval=Interval,nmax=nrow(whatshark),outdir = getwd())
  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date and time
  {
    metro("metro")
    ConvCoord <- LatLon2XY.centered(MyMap.metro,whatshark$Lat[i],whatshark$Long[i])
    points(ConvCoord[1],ConvCoord[2],cex=2,col=2,lwd=2)
    legend(Legcoor.date[1],Legcoor.date[2],whatshark$Time[i],title="Date",title.adj=0.5,
    title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1,x.intersp=-0.5)
    legend(Legcoor.shark[1],Legcoor.shark[2],paste(whatshark$ID[i],"    ",sep=""),title="Shark",title.col="white",
    bg="black",box.col ="white",text.col="white",cex=1.1)

    ani.pause()         #pause between hits for the interval length of time
  }
}

 # save as wmv
        #note: too long a file, R runs out of memory
#saveGIF({Anim.movie.Cont(WhitePointers)}, interval=Interval,movie.name="Allstations_demo.wmv", ani.width=600,ani.height=600)

 # save as html
#saveHTML({Anim.html(WhitePointers)}, imgname="test_plot", title="test",description="ccc",htmlfile = "please.html",
#single.opts = "'controls': ['first', 'previous', 'play', 'next', 'last', 'loop', 'speed'], 'delayMin': 0")



    
      #... Animation used within R
Animation=function(whatshark)
{
  ani.options(interval=Interval,nmax=nrow(whatshark),outdir = getwd())
  metro("metro")
  for( i in 1:ani.options("nmax"))  #loop for plotting all hits by date and time
  {
      #turn previous hit back to original color
      if(i>1)
      {
         if(whatshark$Project[i-1]=="OTN")
        {
          PlotOnStaticMap(MyMap.metro, whatshark$Lat[i-1], whatshark$Long[i-1],cex=1,pch=19,col=Colors.metro[3],add=T)
        }else if
        (whatshark$Lat[i-1]%in%LatLongVR2$Lat)
        {
          PlotOnStaticMap(MyMap.metro, whatshark$Lat[i-1], whatshark$Long[i-1],cex=1,pch=19,col=Colors.metro[1],add=T)
        }else
        {
          PlotOnStaticMap(MyMap.metro, whatshark$Lat[i-1], whatshark$Long[i-1],cex=1,pch=17,col=Colors.metro[2],add=T)
        }
      }

      # plot current hit
      if(whatshark$Project[i]=="OTN")
      {
        PlotOnStaticMap(MyMap.metro, whatshark$Lat[i], whatshark$Long[i],cex=1,pch=19,col="red",add=T)
      }else if
      (whatshark$Lat[i]%in%LatLongVR2$Lat)
      {
        PlotOnStaticMap(MyMap.metro, whatshark$Lat[i], whatshark$Long[i],cex=1,pch=19,col="red",add=T)
      }else
      {
        PlotOnStaticMap(MyMap.metro, whatshark$Lat[i], whatshark$Long[i],cex=1,pch=17,col="red",add=T)
      }

       legend(Legcoor.date[1],Legcoor.date[2],whatshark$Time[i],title="Date",title.adj=0.5,
      title.col="white",bg="black",box.col ="white",text.col="white",cex=1.1,x.intersp=-0.5)

      legend(Legcoor.shark[1],Legcoor.shark[2],paste(whatshark$ID[i],"    ",sep=""),title="Shark",title.col="white",
      bg="black",box.col ="white",text.col="white",cex=1.1)

      ani.pause()         #pause between hits for the interval length of time
  }

}


 # animate hits in VR4s only
 Animation(WhitePointers.VR4s)



      # inset Oz
 #   vp <- baseViewports()
 #   pushViewport(vp$inner,vp$figure,vp$plot)
 #   pushViewport(viewport(x=0.2,y=0.8,width=.3,height=.3,just=c("center","center")))
    #pushViewport(viewport(x=0.9,y=0.8,width=.2,height=.2,just=c("center","center")))
 #   par(fig=gridFIG(),new=T)
    png("H:/Matias WA Fisheries/Analyses/Mapping/WhitePointer_Vr4/OzMetro.png",640,640)
    insetOz("metro")
    dev.off()



#X11()	#reset with and height to defaut



#7.3 Plot of WA with Oz inset
    # download WA
MyMap.wa <- GetMap(center=center.wa, zoom=zoom.wa,destfile = "WA.png",maptype = c("satellite"),frame=T)

#png("Western_Australia.png",640,640,res=100);

		# control width and height of graphic
#X11(width=14,height=12)
#par(mfcol=c(1,1),las=1,mai=c(.3, .3, .6, 0.6),omi=c(.3,.3,.6,0.6))

    # plot WA
PlotOnStaticMap(MyMap.wa, c(SMN_VR2[['2009']]$Lat,SMN_VR2[['2010']]$Lat),c(SMN_VR2[['2009']]$Long,SMN_VR2[['2010']]$Long),
cex=0.5,pch=19,col=Colors.metro[1])
PlotOnStaticMap(MyMap.wa, SMN_VR4G[['2009']]$Lat, SMN_VR4G[['2009']]$Long,cex=0.5,pch=17,col=Colors.metro[2],add=T)
PlotOnStaticMap(MyMap.wa, OTN_VR2[['2009a']]$Lat, OTN_VR2[['2009a']]$Long,cex=0.5,pch=19,col=Colors.metro[3],add=T)

    # add Oz
vp <- baseViewports()
pushViewport(vp$inner,vp$figure,vp$plot)
pushViewport(viewport(x=0.9,y=0.8,width=.2,height=.2,just=c("center","center")))
par(fig=gridFIG(),new=T)
insetOz("WA")

    # add Perth zoom
#par(fig=c(0.125, 0.325, 0.175, 0.375), new = T)
#metro("WA")



#par(usr=c(-216, -66, 24, 144), new = T)
#par(usr=c(-236.7723, -134.0759, -86.53131, -203.49115), new = T)
#metro("WA")


# dev.off()

X11()	#reset with and height to defaut



   #ACA
#Figure 2. VR4 paper
library(rimage)
tagg <- read.jpeg("H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Outputs_VR4_paper/Tagging_Two Peoples Bay_21072010 001.jpg")
VR4 <- read.jpeg("H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Outputs_VR4_paper/_MG_8762.jpg")
VR41 <- read.jpeg("H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Outputs_VR4_paper/_MG_8792.jpg")

png("H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Outputs_VR4_paper/Figure2.png",640,640)
par(mar=c(0.2,0.2,0.2,0.2),oma=c(0.2,0.2,0.2,0.2))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(tagg)
text(70.96876,918.1875,"(a)",col="white",cex=2)
plot(VR4)
text(70.96876,584.26,"(b)",col="black",cex=2)
plot(VR41)
text(70.96876,602.1259,"(c)",col="black",cex=2)
dev.off()



#Figure 1. VR4 paper
Latcoor1= LatLon2XY.centered(MyMap.metro,-31.75,115.89)
Latcoor2= LatLon2XY.centered(MyMap.metro,-32.00,115.89)
Latcoor3= LatLon2XY.centered(MyMap.metro,-32.25,115.89)
Longcoor1= LatLon2XY.centered(MyMap.metro,-31.675,115.20)
Longcoor2= LatLon2XY.centered(MyMap.metro,-31.675,115.45)
Longcoor3= LatLon2XY.centered(MyMap.metro,-31.675,115.70)


map.VR4=function()
{
  PlotOnStaticMap(MyMap.metro, c(SMN_VR2[['2009']]$Lat,SMN_VR2[['2010']]$Lat),c(SMN_VR2[['2009']]$Long,SMN_VR2[['2010']]$Long),
  cex=1.2,pch=19,col=Colors.metro[1])
  PlotOnStaticMap(MyMap.metro, SMN_VR4G[['2009']]$Lat, SMN_VR4G[['2009']]$Long,cex=1.2,pch=17,col=Colors.metro[2],add=T)
  PlotOnStaticMap(MyMap.metro, OTN_VR2[['2009a']]$Lat, OTN_VR2[['2009a']]$Long,cex=1.2,pch=19,col=Colors.metro[3],add=T)

        #add text
    TextOnStaticMap(MyMap.metro, lat=-31.97,lon=115.8, "Perth", cex=2, col = "cyan",add=T)

      # add scale bar
    PlotOnStaticMap(MyMap.metro, c(-32.375,-32.375), c(115.90,115.993),cex=1.75,lwd=4,col="white",FUN=lines,add=T)
    legend(LatcoorScaleB[1], LatcoorScaleB[2],"10 km",bty="n",text.col="white",cex=1.75)

      # add lat and long
  legend(Latcoor1[1], Latcoor1[2],Latvector[1],bty="n",text.col="white",cex=1.5)
  legend(Latcoor2[1], Latcoor2[2],Latvector[2],bty="n",text.col="white",cex=1.5)
  legend(Latcoor3[1], Latcoor3[2],Latvector[3],bty="n",text.col="white",cex=1.5)
  legend(Longcoor1[1], Longcoor1[2],Longvector[1],bty="n",text.col="white",cex=1.5)
  legend(Longcoor2[1], Longcoor2[2],Longvector[2],bty="n",text.col="white",cex=1.5)
  legend(Longcoor3[1], Longcoor3[2],Longvector[3],bty="n",text.col="white",cex=1.5)


      # add legends
  legend(Legcoor[1],Legcoor[2],c("VR2W, PMN","VR4G","VR2W, OTN"),pch=c(19,17),bty="n",col=Colors.metro[1:3],text.col=Colors.metro[1:3],cex=1.5)
}



     # inset Oz
 #   vp <- baseViewports()
 #   pushViewport(vp$inner,vp$figure,vp$plot)
 #   pushViewport(viewport(x=0.2,y=0.8,width=.3,height=.3,just=c("center","center")))
    #pushViewport(viewport(x=0.9,y=0.8,width=.2,height=.2,just=c("center","center")))
 #   par(fig=gridFIG(),new=T)
png("H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Outputs_VR4_paper/Figure1.png",640,640)
map.VR4()
    # add Oz
vp <- baseViewports()
pushViewport(vp$inner,vp$figure,vp$plot)
pushViewport(viewport(x=0.2,y=0.75,width=.3,height=.3,just=c("center","center")))
par(fig=gridFIG(),new=T)

insetOz("metro")
dev.off()














#NOT USED
  #plot OZ
#MyMap <- GetMap(center=center, zoom=zoom,destfile = "Oz.png",maptype = c("satellite"))
#PlotOnStaticMap(MyMap, c(-10,-15,-20), lon = c(130,131,132),axes = TRUE, mar = rep(4, 4))
#PlotOnStaticMap(MyMap, lat = c(), lon = c(),lwd=1.5,col=c('red', 'blue', 'green'), FUN = lines, add=T)
#PlotOnStaticMap(MyMap, lat = c(-10,-15,-20), lon = c(130,131,132), cex=1.5,pch=20,col=c('red', 'blue', 'green'),add=T);
#
