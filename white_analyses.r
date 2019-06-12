####################################################################################
        # GENERAL ANALYSES OF VR4S AND VR2s FOR WHITE SHARKS #
####################################################################################
#notes: this script analyses hits received from tagged white sharks across all receiver lines
#       remember to update data files (AATAMS and SMN)

#       Outputs for VR4 paper are done at the end in "VR4 Outptus"

###################################
#TO DO: SMN all stations positions: use previous files sent by Rory or use unique lat and long
#       AATAMS: use unique lat and long
# Detections data frame: detection data from SMN and AATAMS
# Receivers data frame: location of all receivers
# Transmitters data frame: all sharks with tags

#compare US vs AATAMS webdatabase, is there missing info from webdatabase???


#note: . keep in mind faulty receivers must be flagged to avoid assuming that if no hits, no sharks!!!

# sentinel tags



##################################


#------DATA SECTION------

#source data
setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data")
Detections=read.csv("Detections.csv",stringsAsFactors=F)
AATAMS.all=read.csv("AATAMS.all.csv",stringsAsFactors=F)
SMN.all=read.csv("SMN.all.csv",stringsAsFactors=F)
TAGS=read.csv("TAGS.csv",stringsAsFactors=F)


# Add VR4s list
VR4s=read.csv(file="M:/Fisheries Research/FinFish/Shark/Braccini/Acoustic_tagging_data/SMN_receiver_location/VR4s.csv")


# Add transmitters from S.Africa and SA
South.Africa=read.csv(file="H:/Matias WA Fisheries/Data/Other researcher's tags/ATAP_South african tags_August2012.csv")
South.Oz=read.csv(file="H:/Matias WA Fisheries/Data/Other researcher's tags/Charlie's tags.csv")




# #------ PARAMETER SECTION ------
setwd("H:/Matias WA Fisheries/Analyses/Acoustic_tagging/White_Shark")

 
# min.ping.freq=50     #minimum and max pinging frequency of transmitters, in seconds, for detections <= Jan 2010
# max.ping.freq=130
 
# maxtime=30       #maximum time (in minutes) for succesive hits distribution               (check with Rory)
 
# max.succes.pings=max.ping.freq  #maximum time between succesive pings (in secs)
 



#------ PROCEDURE SECTION ------

#----1. Manipulate data and create data frames----

#Only SMN data
Detections=subset(Detections,Species=="white" & Project=="SMN")
Detections$Latitude=with(Detections,ifelse(Latitude>0,-Latitude,Latitude))

#ONly WA hits
Detections=subset(Detections,Longitude<=129)

#Fix latitude
Detections$ReleaseLatitude=with(Detections,ifelse(ReleaseLatitude=="",NA,ReleaseLatitude))
Detections$ReleaseLatitude=as.numeric(Detections$ReleaseLatitude)
Detections$ReleaseLatitude=with(Detections,ifelse(ReleaseLatitude>0,-ReleaseLatitude,ReleaseLatitude))



#Extract year and month
Detections$DateTime.local=as.POSIXlt(as.character(Detections$DateTime.local))
Detections$Date.local=as.POSIXlt(as.character(Detections$Date.local))
Detections$Time.local=times(as.character(Detections$Time.local))
Detections$ReleaseDate=as.POSIXlt(as.character(Detections$ReleaseDate))


Detections$Year=Detections$Date.local$year+1900
Detections$Month=Detections$Date.local$mon+1 
Detections$Day=Detections$Date.local$mday

Detections$Year.rel=Detections$ReleaseDate$year+1900
Detections$Month.rel=Detections$ReleaseDate$mon+1 
Detections$Day.rel=Detections$ReleaseDate$mday


#Id double-tagged sharks
TAG.tbl=with(Detections,table(TagCode,Name,useNA='ifany'))
TAG.tbl[TAG.tbl>0]=1

Dubl.tag.shks=colSums(TAG.tbl)
Dubl.tag.shks=Dubl.tag.shks[Dubl.tag.shks>1]
Dubl.tag.shks=names(Dubl.tag.shks)


Double.tg.list.Silas=c("Emmy","Lyn","Gabby","Colin","Tania","Cassie","Hamish","Bryn",
"Oscar","Elizabeth","Bethwyn","Nusu 22","Cupcake","Badger","Declan",
"Haley","Bruce","Dickie","Tiny","Chris","Juliet","Jemima","Byron","Luca",
"Jane","Domonic","Westy","Max","Arwyn","Smokey","Digger","Harry","Fran",
"Smudge","Margaux")




############################################################

  #1.2. Create Receivers file
#note: receivers data frame with location of all receivers
not.this.aatams.all=match("X",names(AATAMS.all))
not.this.aatams.all=not.this.aatams.all[!is.na(not.this.aatams.all)]
if(length(not.this.aatams.all)>0)AATAMS.all=AATAMS.all[,-not.this.aatams.all]

not.this.smn.all=match("X",names(SMN.all))
not.this.smn.all=not.this.smn.all[!is.na(not.this.smn.all)]
if(length(not.this.smn.all)>0)SMN.all=SMN.all[,-not.this.smn.all]

Receivers=rbind(SMN.all,AATAMS.all)
rm(list=c("AATAMS.all","SMN.all"))

Receivers$EffectiveStartDate=as.POSIXlt(strptime(as.character(Receivers$EffectiveStartDate),
                                                 format='%d-%m-%y'))
Receivers$EffectiveEndDate=as.POSIXlt(strptime(as.character(Receivers$EffectiveEndDate),
                                                 format='%d-%m-%y'))

  #1.3. Create Transmitters file
#note: transmitters data frame with all sharks with tags
WA.Fisheries$Project="SMN"
South.Africa$Project="South.Africa"
South.Oz$Project="South.Oz"

WA.Fisheries$Study="DoF"  #add study column
South.Africa$Study="SouthAfr"
South.Oz$Study="SARDI"

TAGS=rbind(WA.Fisheries,South.Africa,South.Oz)

TAGS=subset(TAGS,Species2%in%these.species) #keep white sharks and Sentinel tags only

TAGS$Sex2=as.character(TAGS$Sex2)
TAGS$Sex2=with(TAGS,ifelse(Sex2%in%c("","?","U"),NA,
                ifelse(Sex2=="f","F",ifelse(Sex2=="m","M",Sex2))))

TAGS.DoF=subset(TAGS,Study=="DoF")




#----2. Analyses----


  #2.1 All tags deployed

  #Figure 2. Size frequency distribution of all tagged sharks
#note: THIS IS DUMMY, USE REAL DATA. ALSO, NOTE THAT SMN IS FL WHEREAS NEPTUNE ONES IS TL??
TAGS.DoF$Size=with(TAGS.DoF,ifelse(Species2=="White",runif(nrow(TAGS.DoF),200,500),NA))    #dummy, remove!!!!!!!!!!!!!!!!

HISTO=table(TAGS.DoF$Sex2,10*floor(TAGS.DoF$Size/10)) #create histogram
rownames(HISTO)=c("Females","Males")


tiff(file="Outputs_movement/Figure2.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
barplot(HISTO, beside = TRUE,ylim=c(0,max(HISTO)),mgp = c(2, 0.6, 0),
        names.arg= as.numeric(colnames(HISTO)),xlab="Fork length (cm)",ylab="Frequency",cex.lab=1.5, 
        xpd=F,axis.lty=1, axes=T,col=rev(gray(1:nrow(HISTO)/nrow(HISTO))),cex.names=1.1,las=1,cex=1.1)
box()
legend("topright",rownames(HISTO),fill=rev(gray(1:nrow(HISTO)/nrow(HISTO))),yjust=0, horiz=F,bty="n",cex=1.1)
dev.off()



  #2.2 Detections

# Separate sentinel tags from white sharks
Detections.sentinel=subset(Detections,Species=="SENTINEL")
Detections=subset(Detections,Species=="White")


  # Add Julian day from start of release events
First.release=min(Detections$ReleaseDate,na.rm=T)
Detections$Julian=as.numeric(round(Detections$Date.local-First.release))

Start.of.Year.label=as.POSIXlt(as.character(c("2008-01-01","2009-01-01","2010-01-01",
                                        "2011-01-01","2012-01-01")))
Start.of.Year=as.numeric(round(Start.of.Year.label-First.release))

Months.for.plot.label=as.POSIXlt(as.character(c("2008-06-01","2009-06-01","2010-06-01",
                                          "2011-06-01","2012-06-01")))
Months.for.plot=as.numeric(round(Months.for.plot.label-First.release))


  # Summary statistics
Receiver_Tags.table=table(Detections$SerialNumber,Detections$TagCode)
Hits.by.tag=sort(colSums(Receiver_Tags.table))
Hits.by.receiver=sort(rowSums(Receiver_Tags.table))



# Drop false positives (i.e. tags detected in only 1 receiver only one time)              #REVIEW CRITERIA!!!!
#note: this keeps tags detected in more than 1 receiver at least once
drop.false.pos=colSums(Receiver_Tags.table)
drop.false.pos=names(drop.false.pos[drop.false.pos<=1])
Detections.single=subset(Detections,TagCode%in%drop.false.pos)
#Detections=subset(Detections,!(TagCode%in%drop.false.pos))              #REVIEW CRITERIA!!!!



# Drop tag drop-oups (i.e. tags detected in only 1 receiver continuously)

#manual testing, looking at detections of tags by time (Figure 3) and checking
# if tags 8563, 8519 and 8541 were detected by > 1 receiver thru time using
#test=subset(Detections,TagCode==8541);table(test$SerialNumber); plot(test$Julian,test$SerialNumber)


# Flag out sharks tagged east of Neptune is
Detections.TaggedEast=with(Detections,ifelse(ReleaseLongitude>146,"East.Bass.S","West.Bass.S"))


# Check if release < detected dates
TEST.2=(Detections$ReleaseDate-Detections$Date.local)/(24*3600)
id2=which(TEST.2==1)

Matcha=match("Date.local",colnames(Detections))
Detections[id2,Matcha]=Detections[id2,Matcha]+(3600*24)   # for detection day < release day, set detection day to release 
                                                          #day if only time difference of 1 day

TEST=Detections$ReleaseDate<=Detections$Date.local
id=which(TEST==FALSE)
TEST=cbind(Detections$TagCode,TEST)
Which.wrong=which(TEST[,2]==0)
This.wrong=unique(TEST[Which.wrong,1])


Detections=Detections[-id,] #drop sentinel tags used for tagging


# Double-tagged sharks (internal and external acoustic tag)
Double.tagged.internal=c(29584,29533,29546,29529,29469)
Double.tagged.external=c(64100,64115,64123,64125,64119)
Double=c(Double.tagged.internal,Double.tagged.external)

Detections$Double.tagged=with(Detections,ifelse(TagCode%in%Double,"Double","Single"))



  # Add geographical area for plotting
Detections$Area=with(Detections,ifelse(Latitude>-25,"North.WA",ifelse(Longitude>129,"SA",
                  ifelse(Latitude<-25&Longitude<=129,"South.WA",NA))))
unicas.Areas=unique(Detections$Area)
unicas.Areas=unicas.Areas[!is.na(unique(Detections$Area))]
unicas.Areas=unicas.Areas[match(c("North.WA","South.WA","SA"),unicas.Areas)]

Detections$ReleaseLatitude=-Detections$ReleaseLatitude
Detections$Area.release=with(Detections,ifelse(ReleaseLatitude>-25,"North.WA",ifelse(ReleaseLongitude>129,"SA",
                       ifelse(ReleaseLatitude<-25&ReleaseLongitude<=129,"South.WA",NA))))

Receivers$Area=with(Receivers,ifelse(latitude>-25,"North.WA",ifelse(longitude>129,"SA",
                      ifelse(latitude<-25&longitude<=129,"South.WA",NA))))


#Distinguish receiver type
VR4.receivers=unique(VR4s$SerialNumber)
Detections$Rec.Type=with(Detections,ifelse(SerialNumber%in%VR4.receivers,"VR4","VR2"))





  #Figure 1. All receiver locations 

data(worldLLhigh)

  #define coordinates of plots
# North.WA.lat=c(-24,-13); North.WA.long=c(112,122)
# South.WA.lat=c(-36,-30); South.WA.long=c(112,119)
# SA.lat=c(-36,-30); SA.long=c(134,139)
# OZ.lat=c(-44.5,-11);OZ.long=c(110,155)

North.WA.lat=c(-23.5,-21.5); North.WA.long=c(113,114.5)
South.WA.lat=c(-35.5,-31.5); South.WA.long=c(114,119)
SA.lat=c(-36,-33); SA.long=c(135,139)
OZ.lat=c(-44.5,-11);OZ.long=c(110,155)
  
plotlat=list(North.WA.lat,OZ.lat,South.WA.lat,SA.lat)
plotlong=list(North.WA.long,OZ.long,South.WA.long,SA.long)

fn.seq=function(Range)seq(Range[1]+1,Range[2]-1)
Lat.seq=list(fn.seq(North.WA.lat),fn.seq(OZ.lat),fn.seq(South.WA.lat),fn.seq(SA.lat))
Long.seq=list(fn.seq(North.WA.long),fn.seq(OZ.long),fn.seq(South.WA.long),fn.seq(SA.long))

  #define coordinates of polygons
N.WA.long=c(North.WA.long[2], North.WA.long[2], North.WA.long[1], North.WA.long[1])
N.WA.lat=c(North.WA.lat[2], North.WA.lat[1], North.WA.lat[1], North.WA.lat[2])
S.WA.long=c(South.WA.long[2], South.WA.long[2], South.WA.long[1], South.WA.long[1])
S.WA.lat=c(South.WA.lat[2], South.WA.lat[1], South.WA.lat[1], South.WA.lat[2])
SA.long=c(SA.long[2], SA.long[2], SA.long[1], SA.long[1])
SA.lat=c(SA.lat[2], SA.lat[1], SA.lat[1], SA.lat[2])

  #define Perth and Rotnest
Perth=c(115.866,-31.95)
Rotnest=c(115.50,-32.02)

  #define coordinates of receivers
Receiverlat=list(subset(Receivers,Area=="North.WA")$latitude,NA,
                 subset(Receivers,Area=="South.WA")$latitude,subset(Receivers,Area=="SA")$latitude)
Receiverlong=list(subset(Receivers,Area=="North.WA")$longitude,NA,
                  subset(Receivers,Area=="South.WA")$longitude,subset(Receivers,Area=="SA")$longitude)

#bathymetry
 Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
 xbat=sort(unique(Bathymetry$V1))
 ybat=sort(unique(Bathymetry$V2)) 
# reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))


#legends
Letter.leg=c("A",NA,"B","C")
Letter.leg.coor=cbind(c(113.17,NA,114.45,135.4),c(-21.7,NA,-31.85,-33.3))

tiff(file="Outputs_movement/Figure1.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
m <- rbind(c(0.05, 0.3, 0.5, .99),
           c(0.40, .99, 0.5, .99),
           c(0.001, 0.6, 0.01, 0.5),
           c(0.6, .99, 0.01, 0.5))
split.screen(m)

for(i in 1:4)
{
  screen(i)
  par(mar = c(0, 0, 0, 0),mgp=c(.1, 0.15, 0))
  plotMap(worldLLhigh, xlim=plotlong[[i]],ylim=plotlat[[i]],plt = c(.001, 1, 0.075, 1),
          col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  
  if(i==2)
  {
    polygon(x=N.WA.long,y=N.WA.lat,lwd=2)
    polygon(x=S.WA.long,y=S.WA.lat,lwd=2)
    polygon(x=SA.long,y=SA.lat,lwd=2)
    text(133,-25,("Australia"),col="black", cex=2)
    mtext("Latitude (ºS)",side=2,line=0,las=3,cex=2)
    mtext("Longitude (ºE)",side=1,line=0,cex=2)
    text(115.98,-22.6,("A"),col="black", cex=1.75)
    text(120.6,-31.1,("B"),col="black", cex=1.75)
    text(136.9,-31.4,("C"),col="black", cex=1.75)
  }
  if(i==3)
  {polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
   polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
  }
  if(!i==2)
    {
      points(Receiverlong[[i]],Receiverlat[[i]],col=1,pch=19)  #receiver location
      #points()  #hits  
      axis(side = 1, at =Long.seq[[i]], labels = Long.seq[[i]], tcl = .5,las=1,cex.axis=0.9)
      axis(side = 2, at = Lat.seq[[i]], labels = -Lat.seq[[i]],tcl = .5,las=2,cex.axis=0.9)
      box(lwd=2)
       contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=plotlat[[i]],xlim=plotlong[[i]], zlim=c(-1,-300),
               nlevels = 3,labcex=1.25,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
       text(Letter.leg.coor[i,1],Letter.leg.coor[i,2],Letter.leg[i],cex=1.75)
    }

}
close.screen(all = TRUE)
dev.off()





######################################################################
#DELETE ALL THIS ONCE SM NETWORK DATABASE IS RUNNING


# #2. Create dataframes with location of stations for each year
# SMN_VR2_2008$year=2008
# SMN_VR2_2009$year=2009
# SMN_VR2_2010$year=2010
# SMN_VR2_2010$Lat=ifelse(SMN_VR2_2010$Lat>0,-SMN_VR2_2010$Lat,SMN_VR2_2010$Lat) #correct missing -
# SMN_VR2_2011$year=2011
# SMN_VR2_2011$Lat=ifelse(SMN_VR2_2011$Lat>0,-SMN_VR2_2011$Lat,SMN_VR2_2011$Lat) #correct missing -
# SMN_VR4G_2011$year=2011
# SMN_VR4G_2009$year=2009
# OTN_2009a$year=2009
# OTN_2009b$year=2009
# 
# thesecolumns=match(c("Station No","Lat","Long","Serial No","year"),names(SMN_VR2_2008))
# thesecolumns11=match(c("Station No","Lat","Long","Serial No","year"),names(SMN_VR2_2011))
# 
# SMN_VR2=list()
# SMN_VR2[['2008']]=SMN_VR2_2008[,thesecolumns]
# SMN_VR2[['2009']]=SMN_VR2_2009[,thesecolumns]
# SMN_VR2[['2010']]=SMN_VR2_2010[,thesecolumns]
# SMN_VR2[['2011']]=SMN_VR2_2011[,thesecolumns11]
# 
# 
# thesecolumnsVR4=match(c("Station No","Lat","Long","VR4G Serial No","year"),names(SMN_VR4G_2011))
# thesecolumnsOTN09a=match(c("Station No","Lat","Long","Serial No","year"),names(OTN_2009a))
# thesecolumnsOTN09b=match(c("Station No","Lat","Long","Serial No","year"),names(OTN_2009b))
# 
# SMN_VR4G=list()
# SMN_VR4G[['2009']]=SMN_VR4G_2009[,thesecolumnsVR4]
# SMN_VR4G[['2011']]=SMN_VR4G_2011[,thesecolumnsVR4]
# 
# OTN_VR2=list()
# OTN_VR2[['2009a']]=OTN_2009a[,thesecolumnsOTN09a]
# OTN_VR2[['2009b']]=OTN_2009b[,thesecolumnsOTN09b]
# 
# # extract lat and long of receivers for plotting
# LatLongVR2=SMN_VR2_2009[,match(c("Lat","Long"),names(SMN_VR2_2009))]
# LatLongVR4=SMN_VR4G_2009[,match(c("Lat","Long"),names(SMN_VR4G_2009))]
# LatLongOTN=OTN_2009a[,match(c("Lat","Long"),names(OTN_2009a))]
######################################################################




#Figure 3. Timeline of hits by day by shark and by sentinel tag
    #create table
Table.hits.date1=with(Detections,table(TagCode,as.character(Date.local)))
Table.hits.date=with(Detections,table(TagCode,Julian))
First.Jul=colnames(Table.hits.date)[1] #this is first julian date
Hits.date=colnames(Table.hits.date1)

Table.hits=table(Detections$TagCode)
Range.hits=range(Table.hits)

    # add sex, size and release area
these.ones=match(c("TagCode","Sex","Area.release","ReleaseDate"),names(Detections))
noduplicates <- Detections[!duplicated(Detections$TagCode),these.ones]
noduplicates$Sex=ifelse(as.character(noduplicates$Sex)=="","U",as.character(noduplicates$Sex))
noduplicates=noduplicates[order(noduplicates$TagCode),]
Table.hits.date=cbind(noduplicates,as.data.frame.matrix(Table.hits.date))

Rel.date.Julian=as.numeric(round((Table.hits.date$ReleaseDate-First.release)/(3600*24)))
This.col=match(First.Jul,colnames(Table.hits.date))

    # create colors by area
Table.tag.day.area=with(Detections,table(TagCode,Julian,Area))
Table.tag.day.area[,,1]=ifelse(Table.tag.day.area[,,1]>0,1,0)   #convert observations to levels of areas (1:WA, 10:SA)
Table.tag.day.area[,,2]=ifelse(Table.tag.day.area[,,2]>0,10,0)
Table.tag.day.area[,,3]=ifelse(Table.tag.day.area[,,3]>0,1,0)
Table.tag.day.area=Table.tag.day.area[,,1]+Table.tag.day.area[,,2]+Table.tag.day.area[,,3]

Mat.colors=as.data.frame.matrix(Table.tag.day.area)
if(sum(rownames(Mat.colors)==Table.hits.date$TagCode)==nrow(Mat.colors))
  {
    Mat.colors$Area.release=Table.hits.date$Area.release
    Mat.colors$Area.release=ifelse(Mat.colors$Area.release%in%c("South.WA","North.WA"),2,
                            ifelse(Mat.colors$Area.release=="SA",3,NA))
  }
Add.This=Mat.colors$Area.release  
Mat.colors=Mat.colors[,-match("Area.release",colnames(Mat.colors))]
Mat.colors=Add.This*Mat.colors

Mat.colors=ifelse(Mat.colors==30,"brown3",ifelse(Mat.colors==3,"blue",ifelse(Mat.colors==2,"forestgreen",
                    ifelse(Mat.colors==20,"black",NA))))

Mat.colors.legends=c("Rel SA-Det SA","Rel SA-Det WA","Rel WA-Det WA","Rel WA-Det SA")
                                #   col brown3= Rel SA, Detect SA
                                      # blue= Rel SA, Detect WA
                                      # forestgreen= Rel WA, Detect WA
                                      # black= Rel WA, Detect SA
Mat.colors.legends.colors=c("brown3","blue","forestgreen","black")


    # plot
tiff(file="Outputs_movement/Figure3.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mai=c(.8,.8,.01,.01),las=1,mgp=c(2, 1, 0))

bubble.plot.detections(as.numeric(colnames(Table.hits.date[,This.col:ncol(Table.hits.date)])),
       1:nrow(Table.hits.date),
       Table.hits.date[,This.col:ncol(Table.hits.date)],3,"Date","Tag",Rel.date.Julian,
                       Mat.colors,Mat.colors.legends,Mat.colors.legends.colors,"proportion","topleft")
axis(side = 1, at = Months.for.plot,labels = Months.for.plot.label, tcl = -0.3,cex.axis=.75,padj=-2) # minor ticks
axis(side = 1, at = Start.of.Year, labels = Start.of.Year.label, tcl = -0.75) # major ticks & labels

axis(side = 2, at = 1:nrow(Table.hits.date), labels = Table.hits.date$TagCode, tcl = -0.5,cex.axis=.5,hadj=0.75) # major ticks & labels
#axis(side = 2, at = 1:nrow(Table.hits.date), labels = Table.hits.date$Sex, tcl = -0.5,cex.axis=.5,hadj=4) # major ticks & labels
#axis(side = 2, at = 1:nrow(Table.hits.date), labels = Table.hits.date$Area, tcl = -0.5,cex.axis=.5,hadj=3) # major ticks & labels
dev.off()






############# ACA, keep working.... #####################
############################################
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


dummyDATES=data.frame(date=DATES,Year=NA,Month=NA,Day=NA,Time=Dateshours,"Event no#"=NA,"Event ID"=NA,ID=NA,"Receiver S/N"=NA,
Project=NA,Lat=NA,Long=NA,Station=NA,Location=NA,"Depth (m)"=NA)
names(dummyDATES)=names(WhitePointers)


WhitePointers=rbind(WhitePointers,dummyDATES)
WhitePointers=WhitePointers[order(WhitePointers$Time),]   #order by date and time

WhitePointers$hit=ifelse(is.na(WhitePointers$"Receiver S/N")==T,0,1)    #create dummy hit column

WhitePointers$HrsMins=as.POSIXlt(WhitePointers$Time)$hour+(as.POSIXlt(WhitePointers$Time)$min/60)  #create hours an minutes
WhitePointers$HrsMins=ifelse(WhitePointers$HrsMins==0,NA,WhitePointers$HrsMins)

WhitePointers$Date=as.Date(WhitePointers$Date)  #convert back to date

#5. Create useful vectors
Sharks=unique(WhitePointers$ID)     #unique sharks
Sharks=subset(Sharks,is.na(Sharks)==F)
Sharks=sort(Sharks)
n.shark=length(Sharks)


Shark.Times=unique(WhitePointers$Time)     #unique date and time

Shark.colors=rainbow(n.shark)       #colors for these sharks if detected simultaneously
Shark.colors=sample(Shark.colors, n.shark, replace = FALSE)       #make sequential colors different

  #add dummy shark order for plotting
WhitePointers$ID.number=match(WhitePointers$ID,Sharks)

#6. Create dataframes with location of hits received by VR4s and by SMN receivers only
  #dummy NA data to keep all days
NA.Data=WhitePointers
KeepThesVars=match(c("Date","Time","Receiver S/N","HrsMins"),names(WhitePointers))
ALLnames=match(names(WhitePointers),names(WhitePointers))
ALLnames=ALLnames[-KeepThesVars]
NA.Data[,ALLnames]=NA

  #VR4 hits only
WhitePointers.VR4s=subset(WhitePointers,Lat%in%na.omit(LatLongVR4$Lat))
drop.VR4s=match(WhitePointers.VR4s$Time,NA.Data$Time)
WhitePointers.VR4s=rbind(WhitePointers.VR4s,NA.Data[-drop.VR4s,])
WhitePointers.VR4s=WhitePointers.VR4s[order(WhitePointers.VR4s$Date),]

  #SMN hits only
WhitePointers.SMN=subset(WhitePointers,Project=="SMN")
drop.SMN=match(WhitePointers.SMN$Time,NA.Data$Time)
WhitePointers.SMN=rbind(WhitePointers.SMN,NA.Data[-drop.SMN,])
WhitePointers.SMN=WhitePointers.SMN[order(WhitePointers.SMN$Date),]


  #matching deployed dataframes
Tag.Deployed.SMN=subset(Tag.Deployed,ID%in%unique(WhitePointers.SMN$ID[!is.na(WhitePointers.SMN$ID)]))
Tag.Deployed.VR4s=subset(Tag.Deployed,ID%in%unique(WhitePointers.VR4s$ID[!is.na(WhitePointers.VR4s$ID)]))



#7. General bar charts

    #--hits by day, receivers and time--
# get major and minor multiples for choosing labels
ntick.Days=20  #control the number of days between labels
ntick.Receiver=ntick.Time=2
MinorMajorTick=function(ntick)
{ if (ntick < 16) mult = c(2, 2)
    else if(ntick < 41) mult = c(5, 5)
    else if (ntick < 101) mult = c(10, 5)     else mult = c(20, 5)
    return(mult=mult)
}

#create range of hours
horas=1:23
horas.range=vector()
for (i in 0:length(horas-2)) horas.range=c(horas.range,paste(i,":","00","-",i+1,":","00",sep=""))

#plotting function
I.plot.You=function(database,ActiveFile)
{
  #hits by day
      #create table
  table.ID.by.Date=table(database$ID,database$Date)
  date.names=colnames(table.ID.by.Date)   #drop some dates for labelling
  thesecolors=Shark.colors[match(rownames(table.ID.by.Date),Sharks)]  #choose colors to keep it consistent among datasets
       #select ticks
  mult.Days=MinorMajorTick(ntick.Days)
  datalabels <- seq(from = 1, along.with = date.names)
  label.index <- which(datalabels %% mult.Days[1] == 0)
  minor.index = which(datalabels %% mult.Days[2] == 0) # Draw all ticks:
      #plot
  maxY= round(max(table.ID.by.Date)*1.2)
  png(file=paste("hits by day.",ActiveFile,".png",sep=""),width=900,height=600)
  par(las=1,mar=c(6,4,1,1), oma=c(1,1,1,1))
  bp=barplot(table.ID.by.Date,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="Detection numbers",args.legend = list(bty="n",cex=1.1),
  main=paste("2009 Shark Monitoring Network Detections by date (",ActiveFile,")",sep=""),cex.main=1.0,font.main=1,legend=T,
  ylim=c(0,maxY),las=2,col=thesecolors,cex.names=1.0,cex.lab=1.1,cex.axis=1.1,axis.lty=4)
  axis(side = 1, at = bp, labels = FALSE, tcl = -0.2) # Draw minor ticks:
  axis(side = 1, at = bp[minor.index], labels = FALSE, tcl = -0.5) # Draw major ticks & labels:
  axis(side = 1, at = bp[label.index], labels = date.names[label.index], tcl = -0.7,las=2,cex.axis=1)
  mtext("Date",side=1,outer=T,line=-0.2,font=1,cex=1)
  box()
  dev.off()

 #hits by station
      #create table
  table.ID.by.Receiver=table(database$ID,database$"Receiver S/N")
  Receiver.names=colnames(table.ID.by.Receiver)   #drop some dates for labelling
       #select ticks
  mult.Receiver=MinorMajorTick(ntick.Receiver)
  datalabels <- seq(from = 1, along.with = Receiver.names)
  label.index <- which(datalabels %% mult.Receiver[1] == 0)
  minor.index = which(datalabels %% mult.Receiver[2] == 0) # Draw all ticks:
      #plot
  maxY= round(max(table.ID.by.Receiver)*1.2)
  png(file=paste("hits by receiver.",ActiveFile,".png",sep=""),width=900,height=600)
  par(las=1,mar=c(6,4,1,1), oma=c(1,1,1,1))
  bp=barplot(table.ID.by.Receiver,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="Detection numbers",args.legend = list(bty="n",cex=1.1),
  main=paste("2009 Shark Monitoring Network Detections by receiver (",ActiveFile,")",sep=""),cex.main=1.0,font.main=1,legend=T,
  ylim=c(0,maxY),las=2,col=thesecolors,cex.names=1.0,cex.lab=1.1,cex.axis=1.1,axis.lty=4)
  axis(side = 1, at = bp, labels = FALSE, tcl = -0.2) # Draw minor ticks:
  axis(side = 1, at = bp[minor.index], labels = FALSE, tcl = -0.5) # Draw major ticks & labels:
  axis(side = 1, at = bp[label.index], labels = Receiver.names[label.index], tcl = -0.7,las=2,cex.axis=1)
  mtext("Receiver",side=1,outer=T,line=-0.2,font=1,cex=1)
  box()
  dev.off()

 #daily hits
      #create table
  table.ID.by.time=table(database$ID,floor(database$HrsMins))
  if(is.na(match(0,floor(database$HrsMins)))==T) table.ID.by.time=cbind("0"=rep(0,nrow(table.ID.by.time)),table.ID.by.time) #add dummy if no records at 0 hours
  Time.names=horas.range
      #plot
  maxY= round(max(table.ID.by.time)*1.2)
  png(file=paste("hits by time.",ActiveFile,".png",sep=""),width=900,height=600)
  par(las=1,mar=c(6,4,1,1), oma=c(1,1,1,1))
  bp=barplot(table.ID.by.time,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="Detection numbers",args.legend = list(bty="n",cex=1.1),
  main=paste("2009 Shark Monitoring Network Detections by time of day (",ActiveFile,")",sep=""),cex.main=1.0,font.main=1,legend=T,
  ylim=c(0,maxY),las=2,col=thesecolors,cex.names=1.0,cex.lab=1.1,cex.axis=1.1,axis.lty=4)
  axis(side = 1, at = bp, labels = Time.names, tcl = -0.7,las=2,cex.axis=1)
  mtext("Time",side=1,outer=T,line=-0.2,font=1,cex=1)
  box()
  dev.off()
}

# submit functions to create grapsh
I.plot.You(WhitePointers,"OTN & SMN")       #all receivers
I.plot.You(WhitePointers.SMN,"SMN")         #SMN only
I.plot.You(WhitePointers.VR4s,"VR4s")       #VR4s only


    #--hits by receiver type--
    
WhitePointers$Rec.Type=with(WhitePointers,ifelse(Lat%in%na.omit(LatLongVR4$Lat),"VR4",ifelse(Lat%in%na.omit(LatLongVR2$Lat),
"VR2.SMN",ifelse(Lat%in%na.omit(LatLongOTN$Lat),"VR2.OTN",NA))))
table.ID.by.Type=table(WhitePointers$ID,WhitePointers$Rec.Type)
Receiver.names=colnames(table.ID.by.Type)   #drop some dates for labelling
thesecolors=Shark.colors[match(rownames(table.ID.by.Type),Sharks)]  #choose colors to keep it consistent among datasets
# plot
maxY= round(max(table.ID.by.Type)*1.2)
png(file="hits by receiver type.png",width=900,height=600)
par(las=1,mar=c(6,4,1,1), oma=c(1,1,1,1))
bp=barplot(table.ID.by.Type,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="Detection numbers",args.legend = list(bty="n",cex=1.1),
main="2009 Shark Monitoring Network Detections by receiver type",cex.main=1.0,font.main=1,legend=T,
ylim=c(0,maxY),las=2,col=thesecolors,cex.names=1.0,cex.lab=1.1,cex.axis=1.1,axis.lty=4)
axis(side = 1, at = bp, labels = Receiver.names, tcl = -0.7,las=1,cex.axis=1)
mtext("Receiver type",side=1,outer=T,line=-0.2,font=1,cex=1)
box()
dev.off()






#9.1 Create presence of tagged sharks thru time

#plotting function
PresenceTagTime=function(tagshark,database,ActiveFile)
{

    #date range
  daterange=c(min(as.numeric(tagshark$Date)),max(as.numeric(database$Date)))
  datadailyticks <- seq(from = min(daterange), max(daterange))

    #get major and minor multiples for choosing labels
  date.names=unique(sort(database$Date))
  date.names=c(seq(min(tagshark$Date),date.names[1]-1, "days"),date.names)   #add tagging date

  mult.Days=7       #days in between labels
  datalabels <- seq(from = min(as.numeric(date.names)), along.with = date.names)
  label.index = which(datalabels %% mult.Days == 0)
  Atshark=sort(unique(database$ID[!is.na(database$ID)]))

  png(file=paste("PresenceTagThruTime.",ActiveFile,".png",sep=""),width=900,height=600)
  par(las=1,mar=c(6,4,1,1), oma=c(1,1,1,1),mgp=c(3,0.75,0))
  with(database,{
    plot(Date,as.factor(ID),type="n", ylab="Tag ID code",xlab="",yaxt='n',xaxt='n',cex.main=1,
    main=paste("Presence through time (",ActiveFile,")",sep=""),xlim=daterange)
    points(Date,as.factor(ID),pch=19,col="black",cex=1.25)
    axis(side = 2, at =as.factor(Atshark) , labels = Atshark, tcl = -0.5,las=2,cex.axis=1)
    axis(side = 1, at = datadailyticks, labels = FALSE, tcl = -0.2) # Draw minor ticks
    axis(side = 1, at = datadailyticks[label.index], labels = date.names[label.index],
    tcl = -0.5,las=2,cex.axis=0.85) # Draw major ticks & labels
    mtext("Date",side=1,outer=T,line=-1.5,font=1,cex=1)
  })
  if(length(Atshark)==length(Sharks))points(tagshark$Date,as.factor(tagshark$ID),pch=17,col="gray42",cex=1.25)
  if(length(Atshark)<length(Sharks))points(tagshark$Date,as.factor(as.character(tagshark$ID)),pch=17,col="gray42",cex=1.25)
  dev.off()

  TableHitsDate=table(database$Date,database$ID)
  TableHitsDate=ifelse(TableHitsDate>0,1,0)
  TableHitsDate=colSums(TableHitsDate)
  return(TableHitsDate=TableHitsDate)
}

# submit functions to create grapsh
PresenceTagTime(Tag.Deployed,WhitePointers,"OTN & SMN")       #all receivers
PresenceTagTime(Tag.Deployed.SMN,WhitePointers.SMN,"SMN")         #SMN only
PresenceTagTime(Tag.Deployed.VR4s,WhitePointers.VR4s,"VR4s")       #VR4s only



#9.2 Create Table 1 and Figure 3 for publication VR4 use
 hits.by.receiver.type.ID=print(table.ID.by.Type)
 hits.by.receiver.type=print(colSums(table.ID.by.Type))
 Table1=rbind(hits.by.receiver.type.ID,hits.by.receiver.type)
 rownames(Table1)[nrow(Table1)]="Total"
 write.csv(Table1,file="H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Outputs_VR4_paper/Table1.csv", rownames=T)
 uniqueVR4.receivers=unique(VR4$"Receiver S/N")
 
#plotting function
Figure_3=function(tagshark,database,ActiveFile)
{

    #date range
  daterange=c(min(as.numeric(tagshark$Date)),max(as.numeric(database$Date)))
  datadailyticks <- seq(from = min(daterange), max(daterange))

    #get major and minor multiples for choosing labels
  date.names=unique(sort(database$Date))
  date.names=c(seq(min(tagshark$Date),date.names[1]-1, "days"),date.names)   #add tagging date

  mult.Days=7       #days in between labels
  datalabels <- seq(from = min(as.numeric(date.names)), along.with = date.names)
  label.index = which(datalabels %% mult.Days == 0)
  Atshark=sort(unique(database$ID[!is.na(database$ID)]))

  png(file="H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Outputs_VR4_paper/Figure_3.png",width=900,height=600)
  par(las=1,mar=c(6,4,1,1), oma=c(1,1,1,1),mgp=c(3,0.75,0))
  OTN=subset(database,Rec.Type=="VR2.OTN")
  VR4=subset(database,Rec.Type=="VR4")
  VR2.SMN=subset(database,Rec.Type=="VR2.SMN")
  with(database,{
    plot(Date,ID.number,type="n", ylab="Tag ID code",xlab="",yaxt='n',xaxt='n',cex.main=1,
    xlim=daterange)
    points(Date,as.factor(ID),pch=19,col="white",cex=1.25)
    axis(side = 2, at =as.factor(Atshark) , labels = Atshark, tcl = -0.5,las=2,cex.axis=1)
    axis(side = 1, at = datadailyticks, labels = FALSE, tcl = -0.2) # Draw minor ticks
    axis(side = 1, at = datadailyticks[label.index], labels = date.names[label.index],
    tcl = -0.5,las=2,cex.axis=0.85) # Draw major ticks & labels
    mtext("Date",side=1,outer=T,line=-1.5,font=1,cex=1)
  })
  with(OTN,points(Date,ID.number,pch=19,col="blue",cex=1.25))
  with(VR2.SMN,points(Date,ID.number+0.1,pch=19,col="red",cex=1.25))
  with(VR4,points(Date,ID.number-0.1,pch=19,col="gold",cex=1.25))
  if(length(Atshark)==length(Sharks))points(tagshark$Date,as.factor(tagshark$ID),pch=17,col="chartreuse4",cex=1.25)
  if(length(Atshark)<length(Sharks))points(tagshark$Date,as.factor(as.character(tagshark$ID)),pch=17,col="chartreuse4",cex=1.25)
  legend('bottomright',c("tagging date","VR2W, OTN","VR2W, PMN","VR4G"),pch=c(17,19,19,19),bty="n",
  col=c("chartreuse4","blue","red","gold"),text.col="black",cex=1.1)
  dev.off()
}

# submit functions to create grapsh
Figure_3(Tag.Deployed,WhitePointers,"OTN & SMN")       #all receivers





#10.  Create table and figures of movement summaries

table.function=function(tagshark,database,ActiveFile)
{
  Atshark=sort(unique(database$ID[!is.na(database$ID)]))

  Tabla=data.frame(Transmitter.number=tagshark$ID,Total.length.cm=tagshark$Length,Sex=tagshark$Sex,Release.date=tagshark$Date,
  Release.location=tagshark$Location)

  hits=table(database$ID)
  receivers=rowSums(ifelse(table(database$ID,database$"Receiver S/N")>0,1,0))
  days.monitored=as.numeric(max(max(database$Date),as.Date("2011-11-01"))-tagshark$Date)
  days.detected=rowSums(ifelse(table(database$ID,database$Date)>0,1,0))
  Min.cons.det.days=Max.cons.det.days=days.start.end.det=Total.cons.det.days=NULL
  for(i in 1:length(Atshark))
  {
    datos=subset(database,ID==Atshark[i])
    datos=datos[!(duplicated(datos$Date)),]
    Breaks <- c(0, which(diff(datos$Date) != 1), length(datos$Date))
    consdays=sapply(seq(length(Breaks) - 1),function(i) length(datos$Date[(Breaks[i] + 1):Breaks[i+1]]))
    Min.cons.det.days=rbind(Min.cons.det.days,min(consdays))
    Max.cons.det.days=rbind(Max.cons.det.days,max(consdays))
    Total.cons.det.days=c(Total.cons.det.days,consdays)
    if(length(consdays)==1) days.start.end.det=rbind(days.start.end.det,NA)
    if(length(consdays)>1) days.start.end.det=rbind(days.start.end.det,max(datos$Date)-min(datos$Date))
  }

  Residency=days.detected/days.start.end.det

  Tabla$Number.of.hits=hits
  Tabla$Number.of.receivers=receivers
  Tabla$Number.of.days.monitored=days.monitored
  Tabla$Number.of.days.detected=days.detected
  #Tabla$Min.Consecutive.days.present=Min.cons.det.days
  #Tabla$Max.Consecutive.days.present=Max.cons.det.days
  Tabla$Temp.Residence.time.percent=round(100*Residency,0)

  write.table(Tabla,paste("Table_1.",ActiveFile,".csv",sep=""), col.names=T,row.names=F, sep=",")
  
  
      #frequency of continuous detections
  missingdays=match(1:max(Total.cons.det.days),Total.cons.det.days)
  names(missingdays)=1:max(Total.cons.det.days)
  missingdays=missingdays[is.na(missingdays)]
  TableConsDy=table(Total.cons.det.days)
  TableConsDy=c(TableConsDy,missingdays)
  TableConsDy=TableConsDy[match(1:max(Total.cons.det.days),names(TableConsDy))]
  #maxY= round(max(table.ID.by.Date)*1.2)
  png(file=paste("frequency.cont.det.",ActiveFile,".png",sep=""),width=600,height=600)
  bp=barplot(TableConsDy, space = 0,xlab="Days continuously detected",ylab="Frequency",args.legend = list(bty="n",cex=1.1),
  main=paste("Continuous detections (",ActiveFile,")",sep=""),cex.main=1.0,font.main=1,ylim=c(0,max(TableConsDy,na.rm=T)+10),
  las=1,col="gray42",cex.names=1.0,cex.lab=1.1,cex.axis=1.1,axis.lty=1)
  box()
  dev.off()
}

table.function(Tag.Deployed,WhitePointers,"OTN & SMN")       #all receivers
table.function(Tag.Deployed.SMN,WhitePointers.SMN,"SMN")         #SMN only
table.function(Tag.Deployed.VR4s,WhitePointers.VR4s,"VR4s")       #VR4s only



#11. Calculate time intervals between hits
 time.intervals=function(database,ActiveFile)
{
  Atshark=sort(unique(database$ID[!is.na(database$ID)]))
  Time.Intervals=Succesive.Time.Intervals=NULL
  for(i in 1:length(Atshark))
  {
    datos=subset(database,ID==Atshark[i])
    datos=datos[order(datos$Time),]
    time.int=diff.rec=time.int.rec=distance.int=station_t=Station_t_1=NULL
    for (j in 2:nrow(datos))
    {
      time.int=c(time.int,(as.numeric(datos$Time[j])- as.numeric(datos$Time[j-1]))/60)  #time differences in minutes
      diff.rec=c(diff.rec,ifelse(datos$"Receiver S/N"[j]==datos$"Receiver S/N"[j-1],"same","diff"))       #check if succesive hits in different receivers
      delta.Lat=abs(datos$Lat[j]- datos$Lat[j-1])  #lat and long differences
      delta.Long=abs(datos$Long[j]- datos$Long[j-1])
      distance.int=c(distance.int,round((((delta.Lat)^2+(delta.Long)^2)^0.5)*60*1.8*1000,2))  #distance in metres
      Station_t_1=c(Station_t_1,datos$"Receiver S/N"[j])
      station_t=c(station_t,datos$"Receiver S/N"[j-1])
      time.int.rec=data.frame(time.int,diff.rec,ID=unique(datos$ID),dist.int_metres=distance.int,station_t,Station_t_1)

    }
    Time.Intervals=rbind(Time.Intervals,time.int.rec)
   }

    #remove succesive hits < min ping freq and in same receiver
  Time.Intervals.1=Time.Intervals[Time.Intervals$time.int>=min.ping.freq/60,]
  Time.Intervals.2=subset(Time.Intervals,time.int<min.ping.freq/60 & diff.rec=="diff")

  Time.Intervals=rbind(Time.Intervals.1,Time.Intervals.2)
  Time.Intervals=Time.Intervals[order(Time.Intervals$time.int),]
  
  Succesive.Time.Intervals<<-subset(Time.Intervals,time.int<=max.succes.pings/60)

  Quantiles=quantile(Time.Intervals$time.int, probs = c(seq(0, 0.9, 0.1),seq(0.9, 1, 0.01)))
  Succesive.Quantiles=quantile(Succesive.Time.Intervals$time.int, probs = c(seq(0, 0.9, 0.1),seq(0.9, 1, 0.01)))
   
   #all hits
  png(file=paste("Time.intervals.between.hits",ActiveFile,".png",sep=""),width=600,height=600)
  uptohere=which.min(abs(Quantiles - maxtime))
  m=hist(Time.Intervals$time.int[Time.Intervals$time.int<Quantiles[uptohere]],breaks=seq(0,Quantiles[uptohere],1),xaxt='n',
  freq=FALSE,ylim=c(0,0.6),main=paste("Time intervals between hits (",ActiveFile,")",sep=""),ylab="density",xlab="",las=1)
  lines(density(Time.Intervals$time.int[Time.Intervals$time.int<Quantiles[uptohere]],adjust=2,kernel="gaussian"),col=2,lwd=1.5)
  axis(side = 1, at = seq(0.5,Quantiles[uptohere]+0.5,1), labels = FALSE, tcl = -0.2)
  axis(side = 1, at = seq(0.5,Quantiles[uptohere]+0.5,2), labels =seq(1,Quantiles[uptohere]+1,2), tcl = -0.7,las=1,cex.axis=1)
  mtext("Minutes",side=1,outer=T,line=-2,font=1,cex=1)
  x=m$mids[length(m$mids)];y=m$density[length(m$density)]
  arrows(x, y+y*50, x,y,  col= 4)
  text(x-x*.05, y+y*60,paste(names(Quantiles[uptohere]),"of data"),cex = 1,col=4)
  box()
  dev.off()


#   #succesive hits
#  png(file=paste("Succesive.Time.intervals.between.hits",ActiveFile,".png",sep=""),width=600,height=600)
#  uptohere=which.min(abs(Succesive.Quantiles - maxtime))
#  m=hist(Succesive.Time.Intervals$time.int[Succesive.Time.Intervals$time.int<Succesive.Quantiles[uptohere]],
#  breaks=seq(0,Succesive.Quantiles[uptohere],1),xaxt='n',
#  freq=FALSE,ylim=c(0,0.6),main=paste("Time intervals between hits (",ActiveFile,")",sep=""),ylab="density",xlab="",las=1)
#  lines(density(Succesive.Time.Intervals$time.int[Succesive.Time.Intervals$time.int<Succesive.Quantiles[uptohere]],
#  adjust=2,kernel="gaussian"),col=2,lwd=1.5)
#  axis(side = 1, at = seq(0.5,Succesive.Quantiles[uptohere]+0.5,1), labels = FALSE, tcl = -0.2)
#  axis(side = 1, at = seq(0.5,Succesive.Quantiles[uptohere]+0.5,2), labels =seq(1,Succesive.Quantiles[uptohere]+1,2), tcl = -0.7,las=1,cex.axis=1)
#  mtext("Minutes",side=1,outer=T,line=-2,font=1,cex=1)
#  x=m$mids[length(m$mids)];y=m$density[length(m$density)]+0.0025
#  arrows(x, y+y*50, x,y,  col= 4)
#  text(x-x*.05, y+y*60,paste(names(Succesive.Quantiles[uptohere]),"of data"),cex = 1,col=4)
#  box()
#  dev.off()
}
time.intervals(WhitePointers,"OTN & SMN")       #all receivers


  ##
Table.succ.hits.sameOrNot=table(Succesive.Time.Intervals$diff.rec)
 #ACA!!!!!


#XX. Create velocity vector and time-step movement

  #obtain velocity vector (split vector in x and y components)
tagging.Data$delta.Long=tagging.Data$recaptured.Long-tagging.Data$tagged.Long
tagging.Data$delta.Lat=tagging.Data$recaptured.Lat-tagging.Data$tagged.Lat
   WhitePointers
   
   
   
   if(sum(a==d)==2)print("SIII")








############################
## VR4 Outptus ##
############################
Detections.VR4.paper=subset(Detections,Latitude>(-32.5) & Latitude <(-31.5) & Longitude <116)





#Figure 1
insetOz <- function()
  {
  opar <- par(mai = c(1.5,2,1,1))
  on.exit(par(opar))
  par(mar=rep(.1,4),xaxt="n",yaxt="n",plt=par("plt"))
  plotMap(worldLLhigh, xlim=c(110,155), ylim=c(-44.5,-11),col="dark grey", 
          plt = c(.001, 1, 0.075, 1),tck = 0.025, tckMinor = 0.0125,axes=F, xlab="", ylab="")
  text(133,-25,("Australia"),col="black", cex=1.35)
  edgeX=c(114,117,117,114) 
  edgeY=c(-30.25,-30.25,-33.25,-33.25)
  polygon(x=edgeX,y=edgeY,lwd=2)
  box(,lwd=2)
}

  
tiff(file="Outputs_movement/Figure1.VR4.paper.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mar = c(0, 0, 0, 0),mgp=c(.1, 0.15, 0))
plotMap(worldLLhigh, xlim=c(115.18,116),ylim=c(-32.41,-31.7),plt = c(.1, 1, 0.075, 1),
          col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)


polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
points(Receiverlong[[3]],Receiverlat[[3]],col=1,pch=21,bg=1)  #receiver location

points(VR4s$Longitude1,VR4s$Latitude1,col=1,pch=24,cex=1.5,bg="gray80")  #receiver location

axis(side = 1, at =seq(115,115.75,by=.25), labels = seq(115,115.75,by=.25), tcl = .5,las=1,cex.axis=1.25)
axis(side = 2, at = seq(-32.4,-31.6,by=.2), labels = -seq(-32.4,-31.6,by=.2),tcl = .5,las=2,cex.axis=1.25)
box(lwd=2)
contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=plotlat[[i]],xlim=plotlong[[i]], zlim=c(-1,-300),
        nlevels = 3,labcex=1.25,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
text(115.80,-31.96656,"Perth",cex=1.75)
mtext("Latitude (ºS)",side=2,line=2,las=3,cex=2)
mtext("Longitude (ºE)",side=1,line=1.75,cex=2)
legend('bottomright',c("VR2W","VR4G"),pch=c(21,24),bty="n",
       col=1,pt.bg=c("black","gray80"),text.col="black",cex=1.25)


vp <- baseViewports()
pushViewport(vp$inner,vp$figure,vp$plot)
pushViewport(viewport(x=0.05,y=0.4,width=.35,height=.35,just=c("left","top")))
par(fig=gridFIG(),new=T)  
insetOz()


dev.off()




#Table 1
Table1=with(Detections.VR4.paper,table(TagCode,Rec.Type))

write.table(Table1,"Outputs_movement/Table1.VR4.paper.csv", col.names=T,row.names=T, sep=",")



#Figure 3

#create table
Table.hits.date1=with(Detections.VR4.paper,table(TagCode,as.character(Date.local),Rec.Type))
Table.hits.date=with(Detections.VR4.paper,table(TagCode,Julian,Rec.Type))
Table.hits.date=ifelse(Table.hits.date>0,1,0)
First.Jul=as.numeric(colnames(Table.hits.date)[1]) 
Last.Jul=as.numeric(colnames(Table.hits.date)[NNN])
NNN=ncol(Table.hits.date)  
This.Juls=as.numeric(colnames(Table.hits.date))

First.Date=colnames(Table.hits.date1)[1]

b=matrix(ncol=NNN,nrow=nrow(Table.hits.date))
for(i in 1:ncol(b))b[,i]=This.Juls[i]

Table.hits.date.VR2=Table.hits.date[,,1]*b
Table.hits.date.VR4=Table.hits.date[,,2]*b

VR4.paper.TAGS=rownames(Table.hits.date)
Tag.order=1:length(VR4.paper.TAGS)
Tag.order.VR4=Tag.order+.35

Start.of.Year.label=as.POSIXlt(as.character(First.Date))
More.dates=as.POSIXlt(as.character(c("2009-12-01","2010-06-01","2010-12-01",
                                     "2011-06-01","2011-12-01","2012-06-01","2012-12-01")))
Diff.dates=as.numeric(round(More.dates-Start.of.Year.label))
AT.here=c(First.Jul,Diff.dates+First.Jul)
AT.here.labels=c("06-06-2009","01-12-2009","01-06-2010", "01-12-2010"
                 ,"01-06-2011" ,"01-12-2011","01-06-2012", "01-12-2012")

par(mar = c(3.5, 4, 0, 0))
tiff(file="Outputs_movement/Figure3.VR4.paper.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
plot(Table.hits.date.VR2[1,],rep(Tag.order[1],length(Table.hits.date.VR2[1,])),
     xlim=c(First.Jul-5,Last.Jul+28),ylim=c(1,length(VR4.paper.TAGS)),pch=19,col=1,xlab="",ylab="",xaxt='n',,yaxt='n')

for (i in 1:length(Tag.order)) 
{
  points(Table.hits.date.VR2[i,],rep(Tag.order[i],length(Table.hits.date.VR2[i,])),pch=19,col=1)
  points(Table.hits.date.VR4[i,],rep(Tag.order.VR4[i],length(Table.hits.date.VR4[i,])),pch=19,col="grey50")
}  

axis(side = 1, at = AT.here,labels = AT.here.labels, tcl = -0.3,cex.axis=.75,padj=-1.5) # minor ticks
axis(side = 2, at = Tag.order, labels = rownames(Table.hits.date), tcl = -0.5,cex.axis=0.9,hadj=0.8,las=2) # major ticks & labels
mtext("Date",side=1,line=1.85,cex=1.5)
mtext("Tag code",side=2,line=3,cex=1.5)
legend("topleft",c("VR2W", "VR4G"),col=c("black","grey50"),pch=19,bty='n')

dev.off()