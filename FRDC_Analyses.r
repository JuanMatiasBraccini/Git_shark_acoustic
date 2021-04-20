####################################################################################
# ANALYSES FOR FRDC SPATIAL AND TEMPORAL PATTERNS
####################################################################################

#notes: .this script analyses hits received from tagged dusky, copper, sandbar, whiskery and gummy sharks
#         across all receiver lines
#       . remember to update data files (AATAMS and SMN)
#       . See Willis and Hobday 2007 for neat analysis


#PROCEDURE INDEX
#1. --Reported recaptured acoustically tagged sharks ---
#2. --Create Receivers file ---
#3. --Create Transmitters file ---
#4. --Detection manipulations ---
#5. -- Summary table ---
#6. -- Proportion of time within zones ---
#6.1-- Movement rates of gummy and whiskery among zones ---
#7. -- Data set for pop dyn modelling ---
#8. -- Home range (Kernel Density) for Ningaloo ---      #this is incomplete
#9. -- Natal migration ---
#10. -- Distance correction for corners ---
#11. -- Calculate days straight-line movement ---
#12. -- Drop tags detected in only 1 receiver only one time ---
#13. -- Create detections data sets for each species ---
#14. -- Table 2 FRDC milestone reports ---
#15. -- Table 3 FRDC milestone reports ---
#16. -- Histogram of all tags deployed ---
#17. -- Map of study area ---
#18. -- Figure Timeline of hits by day by shark ---
#18.1 -- Behavioural polymorphism ---
#19. -- Dusky natal migration analysis ---
#20. --  Compare sex ratio, condition and size of individuals tagged in WA vs detected  ---
#21. -- Proportion per array ---
#22. -- Bubble plots of proportion of hits by station for each array ---
#23. -- Daily patterns --- 
#24. -- Co-detection of individuals ---
#25. -- Speed and distance travelled ---
#26. -- Cumulative Distance ---
#27. -- Patterns in residency, area use and speed ---
#28. -- MEPS figures ---
#29. -- Fisheries Oceanography figures ---


rm(list=ls(all=TRUE))

library(sp)
library(rgdal)
library(RODBC)
library(plotrix)
library(PBSmapping)
library(chron)
library(geosphere)
library(adehabitat)
library(circlize)
library(diagram)  #for curved arrows
library(lme4)
library(lubridate)
library(oce)   #for day length
library(ReporteRs)
library(randomForest)
library(data.table)

options(stringsAsFactors = FALSE)

#source bubble plot functions
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Bubble.plot.R")) 
source(handl_OneDrive("Analyses/Acoustic_tagging/Git_shark_acoustic/Bubble.plot.detections.R"))


source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Deviance.explained.R"))

###### DATA SECTION ############

#source data
#note: the Detections file is created in "5.Source_acoustic_data.R"
#use.all="NO"   #set to "YES" if also anlysing white, tiger and grey nurse sharks
#source("C:/Matias/Analyses/Acoustic_tagging/5.Source_acoustic_data.R")

setwd(handl_OneDrive("Data/Tagging/Acoustic_tagging/Acoustic_tagging_data"))
Detections=fread("Detections.csv",data.table=FALSE)
AATAMS.all=fread("AATAMS.all.csv",data.table=FALSE)
SMN.all=fread("SMN.all.csv",data.table=FALSE)
TAGS=read.csv("TAGS.csv",stringsAsFactors=F)
#ACA

PerthIs=read.table(handl_OneDrive("Data/Mapping/WAislandsPointsNew.txt", header=T))
Rottnest.Is=subset(PerthIs,ID%in%c("ROTT1"))
Garden.Is=subset(PerthIs,ID%in%c("ROTT3"))
  
setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC"))


#Sharks data base for size of all tagged sharks
Source.FL.externally="NO"
if(Source.FL.externally=="YES")  #if FL not extracted from SMN web
{
  DAT=2  #if reading data from Shark data base
  if(DAT==2)
  {
    setwd("M:/Fisheries Research/Production Databases/Shark")  # working directory
    channel <- odbcConnectAccess("Sharks.mdb")      
    Boat_bio=sqlFetch(channel, "Boat_bio", colnames = F)   
    Boat_hdr=sqlFetch(channel, "Boat_hdr", colnames = F)   
    close(channel)
    DATA=merge(Boat_bio,Boat_hdr,by="SHEET_NO",all.x=T)
    names(DATA)[match(c("SPECIES","ATAG NO"),names(DATA))]=c("Species","TagCode")
    DATA=subset(DATA,Species%in%c("BW","TK","GM","WH","CP")& !is.na(TagCode),select=c(Species,DATE,TagCode,FL,TL))
    DATA=subset(DATA,!(TagCode==29597 & as.character(DATE)=="2013-07-26"))  #remove duplication
  }
  if(DAT==1) DATA=read.csv(handl_OneDrive("Analyses/Acoustic_tagging/DATA.csv"))
  
  DATA$Project.rel='SMN'
  
  #Size and Release info from all South Australian sharks
  #note: update this file if new SA sharks are tagged
  DATA.SA=read.csv(handl_OneDrive("Data/Tagging/Acoustic_tagging/Other researcher's tags/Charlie's tags.csv"))
  DATA.SA$DATE=as.POSIXlt(as.character(DATA.SA$ReleaseDate),format='%d/%m/%Y')
  DATA.SA$Project.rel='South.Australia'
  DATA.SA$FL=100*DATA.SA$FL
  DATA.SA$TL=round(with(DATA.SA,FL*1.1849+2.9835))
  
  DATA=rbind(DATA,DATA.SA[,match(names(DATA),names(DATA.SA))])  #combine SMN and South Australia
  
}

#Release condition
channel <- odbcConnectAccess2007("U:/Shark/Sharks.mdb")   #use for 64 bit
#channel <- odbcConnectAccess("M:/Fisheries Research/Production Databases/Shark/Sharks.mdb")      #32 bit
Condition=sqlFetch(channel, "Tag data", colnames = F)   
close(channel)
names(Condition)[match("ATAG NO",names(Condition))]="ATAG"
Condition=subset(Condition,!is.na(ATAG),select=c(ATAG,SPECIES,FL,CONDITION))


#Shark zones
JA_Northern_Shark=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/JA_Northern_Shark.shp", layer="JA_Northern_Shark")) 
WA_Northern_Shark=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/NorthCoastShark_s43.shp", layer="NorthCoastShark_s43")) 
WA_Northern_Shark_2=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/NorthWestCoastShark_s43.shp", layer="NorthWestCoastShark_s43")) 
SDGDLL_zone1=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/SDGDLL_zone1.shp", layer="SDGDLL_zone1")) 
SDGDLL_zone2=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/SDGDLL_zone2.shp", layer="SDGDLL_zone2")) 
WCDGDLL=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/WCDGDLL.shp", layer="WCDGDLL")) 


#Reported recaptures
setwd("U:/Shark")  # working directory
channel <- odbcConnectAccess2007("Sharks v20200323.mdb")      
Rep.Recap=sqlFetch(channel, "Tag data", colnames = F)   
close(channel)


#Bathymetry
Bathymetry_120=read.table(handl_OneDrive("Data/Mapping/get_data112_120.cgi"))
Bathymetry_138=read.table(handl_OneDrive("Data/Mapping/get_data120.05_138.cgi"))
Bathymetry=rbind(Bathymetry_120,Bathymetry_138)
rm(Bathymetry_138);rm(Bathymetry_120)

Temp=read.csv(handl_OneDrive("Data/SST.nice.format.csv"))
#Temp=read.csv("C:/Matias/Data/Reynolds SST/Coast_temperatures.csv")


###### PARAMETER SECTION ############
START=as.POSIXlt("2011-01-01")   #update start and end of data
END=as.POSIXlt(as.character(sort(unique(Detections$Date.local))
        [length(sort(unique(Detections$Date.local)))]))     

Migration.threshold=c(100,100,100,100) #minimum distance (in km) for movement to be consider migration for Dusky, sandbar, gummy and whiskery
#Migration.threshold=c(200,200,100,100)
detec.range=400       #receiver detection range (in m)
Delta.t=1             #minimum time (in minutes) for succesive hits in different receivers. Parameter for dropping too-close hits 
                      # at adjacent receivers inflating speed!!

# speed.threshold=5        #maximum possible speed (m/s). 
#                         # White shark average speed: 3 km/h (0.84 m/s); hammerhead: 0.65m/s;
#                         # whale shark max speed: 3.9 km/h (1.084 m/s);tiger shark: 3.85 km/h (1.07 m/s);
#                         # misc ranges: 18 m/s,26 m/s,9.7-16 m/s,14 m/s,10-27 m/s; mako shark: 15 m/s
#                         # dusky: 1km/h (0.28 m/s)

#growth parameters
Gr=list(BW=c(K.f=.0367,Linf.f=374.4,to.f=-3.3,K.m=0.045,Linf.m=337,to.m=-3),
        TK=c(K.f=.040,Linf.f=244.2,to.f=-4.8,K.m=0.044,Linf.m=226,to.m=-4),
        GM=c(K.f=0.123,Linf.f=201.9,to.f=-1.55,K.m=0.266,Linf.m=137,to.m=-0.8),
        WH=c(K.f=0.369,Linf.f=120.7,to.f=-0.6,K.m=.423,Linf.m=121.5,to.m=-0.472))
Mx.age=list(BW=45,TK=33,GM=16,WH=15)

#FL at 50% maturity (minimum value of McAuley et al 2007)
FL_0.5.dusky=2.54  #in m
FL_0.5.dusky.male=1.91
FL_0.5.sandbar=1.30
FL_0.5.sandbar.male=1.26
FL_0.5.gummy=(125-4.6424)/1.0837
FL_0.5.gummy=FL_0.5.gummy/100
FL_0.5.gummy.male=(93-4.6424)/1.0837
FL_0.5.gummy.male=FL_0.5.gummy.male/100
FL_0.5.whiskery=1.12
FL_0.5.whiskery.male=1.07
FL_0.5=list(Dusky=c(FL_0.5.dusky,FL_0.5.dusky.male),Thickskin=c(FL_0.5.sandbar,FL_0.5.sandbar.male),
            Gummy=c(FL_0.5.gummy,FL_0.5.gummy.male),Whiskery=c(FL_0.5.whiskery,FL_0.5.whiskery.male))

#Corners for moving around land
Exmouth=cbind(113.843,-21.81416)
Shark.bay=cbind(112.921,-25.497)
Cape.Leuwin=cbind(114.969,-34.459)
Mid.point=cbind(116.425,-35.043)

#minimum time between detections
minHours=1 #set to one hour to be consistent with speed units. Hits must be at least
#           1 hour apart to avoid nonsense high speeds due hits in the outskirt of detection range
#           of the two receivers and not really apart the distance between the receivers (800-1000 m)


#Control colors

  #Arrays
#CLS=c("red","blue","forestgreen")
CLS=c("red","blue","black")
names(CLS)=c("Ningaloo","Perth","Southern.lines")

  #Zones
Zns=c("North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2")
Zns.leg=c("WANCSF","Ningaloo","WCDGDLF","Metro","Zone1","Zone2")
COL.prop=c("lightseagreen","seagreen4","lightgreen","olivedrab4","olivedrab3","mediumseagreen")
names(COL.prop)=Zns
#colfunc <- colorRampPalette(c("black", "white"))
#colfunc <- colorRampPalette(c("cadetblue", "deepskyblue4"))
#COL.prop=colfunc(6)
funky.cols=c(North="navy",Closed.ningaloo="yellowgreen",West="turquoise4",
             Closed.metro="violetred",Zone1="honeydew4",Zone2="slateblue1",SA="red")


  #Single color
Singl.col="grey60"


#Control where movement rates are calculated
do_mov_rates="elsewhere"  #movement rates calculated in "Matias/Analyses/Movement rate estimation/"  
#do_mov_rates="here"



###### PROCEDURE SECTION ############

#1. --Reported recaptured acoustically tagged sharks ---

No.recap.pos=c(30926,30950,29436,30929,29463)   #Reported by J. Tindal but whithout recaptured info
No.recap.date=c(30926,30929,30928,29499,29463)   #Reported by J. Tindal but whithout recaptured info
names(Rep.Recap)[match(c("Captured?","ATAG NO"),names(Rep.Recap))]=c("Recaptured","ATAG_NO")
Rep.Recap=subset(Rep.Recap,Recaptured=="Y" & !is.na(ATAG_NO))
Rep.Recap$RELLATDECDEG=-with(Rep.Recap,REL_LATD+(REL_LATM/60))
Rep.Recap$RECLATDECDEG=-with(Rep.Recap,CAP_LATD+(CAP_LATM/60))
Rep.Recap$RELLNGDECDEG=with(Rep.Recap,REL_LNGD+(REL_LNGM/60))
Rep.Recap$RECLNGDECDEG=with(Rep.Recap,CAP_LNGD+(CAP_LNGM/60))

  #assign Zone 2 recapture position for some unknowns that were recaptured in that zone
Rep.Recap$RECLATDECDEG=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.pos,-35,RECLATDECDEG))
Rep.Recap$RECLNGDECDEG=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.pos,119,RECLNGDECDEG))
Rep.Recap$DATE_CAPTR1=as.character(Rep.Recap$DATE_CAPTR)
Rep.Recap$DATE_CAPTR1=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.date 
                        & is.na(DATE_CAPTR1),"2013-06-01",DATE_CAPTR1))
Rep.Recap$DATE_CAPTR1=as.POSIXlt(Rep.Recap$DATE_CAPTR1)
Rep.Recap$DATE_CAPTR=Rep.Recap$DATE_CAPTR1
Rep.Recap=Rep.Recap[,-match("DATE_CAPTR1",names(Rep.Recap))]
Rep.Recap$Rep.pos=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.pos,"No rec. pos","OK"))
Rep.Recap$Rep.date=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.date,"No rec. date","OK"))

#add recaptured gummy by J. Cooke but with no recapture information
Add.Jeff=Rep.Recap[1,]
Add.Jeff[,]=NA
Add.Jeff$SPECIES="GM"
Add.Jeff$Recaptured="YES"
Add.Jeff$Rep.pos="No rec. pos"
Add.Jeff$Rep.date="OK"
Add.Jeff$DATE_CAPTR=as.POSIXlt("2015-05-15")
Rep.Recap=rbind(Rep.Recap,Add.Jeff)

Tab.rep.rec=table(as.character(Rep.Recap$SPECIES))

setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement"))
this=match(c("SPECIES","ATAG_NO","SEX","Rep.pos","Rep.date","RELEASE DATE","DATE_CAPTR","RELLATDECDEG",
            "RELLNGDECDEG","RECLATDECDEG","RECLNGDECDEG"),names(Rep.Recap))
Re.Rec.Tab=Rep.Recap[order(Rep.Recap$SPECIES,Rep.Recap$ATAG_NO),this]
write.csv(Re.Rec.Tab,"Tab.rep.rec.csv",row.names=F)

#Export recaptures for population dynamics modelling
Rep.Recap$year.rel=year(Rep.Recap$"RELEASE DATE")
Rep.Recap$year.rec=year(Rep.Recap$"DATE_CAPTR")
Rep.Recap$Recaptured=as.character(Rep.Recap$Recaptured)

Pop.din.sp=c("BW",'WH','GM','TK')
for(i in 1:length(Pop.din.sp))      
{
  a=subset(Rep.Recap,SPECIES==Pop.din.sp[i])
  NmS=ifelse(Pop.din.sp[i]=="BW",'Dusky shark',
             ifelse(Pop.din.sp[i]=='WH','Whiskery shark',
                    ifelse(Pop.din.sp[i]=='GM','Gummy shark',
                           ifelse(Pop.din.sp[i]=='TK','Sandbar shark',NA))))
  write.csv(a,paste(handl_OneDrive('Analyses/Data_outs/'),NmS,'/',NmS,"_Acous.Tag_Rep.Recap.csv",sep=""),row.names=F)
}

#Export as word table
Scenarios.tbl=function(WD,Tbl,Doc.nm,caption,paragph,HdR.col,HdR.bg,Hdr.fnt.sze,Hdr.bld,
                       body.fnt.sze,Zebra,Zebra.col,Grid.col,Fnt.hdr,Fnt.body,HEDR,HEDR.cols,HEDR2,HEDR3)
{
  mydoc = docx(Doc.nm)  #create r object
  mydoc = addSection( mydoc, landscape = T )   #landscape table
  # add title
  if(!is.na(caption))mydoc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
  
  # add a paragraph
  if(!is.na(paragph))mydoc = addParagraph(mydoc , paragph, stylename="Citationintense")
  
  #add table
  MyFTable=FlexTable(Tbl,header.column=F,add.rownames =F,
                     header.cell.props = cellProperties(background.color=HdR.bg), 
                     header.text.props = textProperties(color=HdR.col,font.size=Hdr.fnt.sze,
                                                        font.weight="bold",font.family =Fnt.hdr), 
                     body.text.props = textProperties(font.size=body.fnt.sze,font.family =Fnt.body))
  
  #Add header
  MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value = HEDR, colspan = HEDR.cols)
  #Add second header
  if(length(HEDR2)>1)MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value = HEDR2)
  #Add third header
  if(length(HEDR3)>1)MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value = HEDR3)
  
  
  # zebra stripes - alternate colored backgrounds on table rows
  if(Zebra=="YES") MyFTable = setZebraStyle(MyFTable, odd = Zebra.col, even = "white" )
  
  # table borders
  MyFTable = setFlexTableBorders(MyFTable,
                                 inner.vertical = borderNone(),inner.horizontal = borderNone(),
                                 outer.vertical = borderNone(),
                                 outer.horizontal = borderProperties(color=Grid.col, style="solid", width=4))
  
  # set columns widths (in inches)
  #MyFTable = setFlexTableWidths( MyFTable, widths = Col.width)
  
  mydoc = addFlexTable( mydoc, MyFTable) 
  
  mydoc = addSection( mydoc, landscape = F ) 
  
  # write the doc 
  writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
}

idds=match(c("SPECIES","SEX","RELEASE DATE","DATE_CAPTR","RELLATDECDEG","RELLNGDECDEG",
             "RECLATDECDEG","RECLNGDECDEG"),names(Re.Rec.Tab))
Re.Rec.Word=Re.Rec.Tab[,idds]
Re.Rec.Word$SPECIES=with(Re.Rec.Word,ifelse(SPECIES=="BW","Dusky shark",
                 ifelse(SPECIES=="TK","Sandbar shark",
                 ifelse(SPECIES=="WH","Whiskery shark",
                 ifelse(SPECIES=="GM","Gummy shark",
                 ifelse(SPECIES=="CP","Copper shark",NA))))))
Re.Rec.Word[,5:8]=round(Re.Rec.Word[,5:8],2)
Scenarios.tbl(WD=getwd(),Tbl=Re.Rec.Word,Doc.nm="Re.Rec.Tab",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
              HEDR=c('Species','Sex','Date','Release','Recapture'),
              HEDR.cols=c(1,1,2,2,2),HEDR2=c("","","Release","Recapture","Latitude",
              "Longitude","Latitude","Longitude"),HEDR3=NA)



Recap.sp=sort(as.character(unique(Rep.Recap$SPECIES)))
NAMES.rec=c("Dusky shark","Copper shark","Gummy shark","Sandbar shark","Whiskery shark")
data(worldLLhigh)
fn.plt.recap=function(spec,NME)    
{
  dat=subset(Rep.Recap,SPECIES==spec)
  plotMap(worldLLhigh, xlim=xlm,ylim=ylm,plt = c(.001, 1, 0.075, 1),
          col="grey88",xlab="",ylab="",axes=F,main="")
  arrows(dat$RELLNGDECDEG,dat$RELLATDECDEG,dat$RECLNGDECDEG,dat$RECLATDECDEG,
         col=2,lwd=1.7,length=0.1,angle=35)
  mtext(NME,side=3,line=0.25,cex=1.25) 
  axis(side = 1, at =round(xlm[1]):xlm[2], labels = F, tcl = .5,las=1,cex.axis=1.2)
  axis(side = 2, at = round(ylm[1]):ylm[2], labels = F,tcl = .5,las=2,cex.axis=1.2)
  box()
}
ylm=c(-36,-13)
xlm=c(112,129)
#tiff(file="Rep.recap.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
#par(mfcol=c(2,3),mai=c(.5,.5,.1,.1),oma=c(2.5,4,2,.1),mgp=c(1,.7,0))
#for (i in 1:length(Recap.sp))
tiff(file="Rep.recap.tiff",width = 1600, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mai=c(.1,.1,.1,.4),oma=c(4,4,2,.1),mgp=c(1,.8,0))
for (i in c(1,3:5))
{
  fn.plt.recap(Recap.sp[i],NAMES.rec[i])
#   if(i%in%c(2,4))axis(side = 1, at =seq(xlm[1],xlm[2],4), labels = seq(xlm[1],xlm[2],4), tcl = .5,las=1,cex.axis=1.2)
#   if(i%in%1:2)axis(side = 2, at = seq(ylm[1],ylm[2],4), labels = -seq(ylm[1],ylm[2],4),tcl = .5,las=2,cex.axis=1.2)
  if(i%in%c(3,5))axis(side = 1, at =seq(xlm[1],xlm[2],4), labels = seq(xlm[1],xlm[2],4), tcl = .5,las=1,cex.axis=1.2)
  if(i%in%1:3)axis(side = 2, at = seq(ylm[1],ylm[2],4), labels = -seq(ylm[1],ylm[2],4),tcl = .5,las=2,cex.axis=1.2)
}
mtext("    Latitude (ºS)",side=2,line=2,las=3,cex=2,outer=T)
mtext("Longitude (ºE)",side=1,line=2,cex=2,outer=T)
dev.off()



#2. --Create Receivers file ---

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
Receivers$Area=with(Receivers,ifelse(latitude>-29,"North.WA",ifelse(longitude>129,"SA",
                ifelse(latitude<-29 & longitude<=129,"South.WA",NA))))


#3. --Create Transmitters file ---
these.species=c("SENTINEL","Dusky","Thickskin","Gummy","Whiskery")
Other.TAGS=subset(TAGS,!Project.rel=="SMN")
Show.South.Oz="YES" 
if(Show.South.Oz=="YES") TAGS=subset(TAGS,Project.rel%in%c("SMN","South.Australia")) 
if(Show.South.Oz=="NO") TAGS=subset(TAGS,Project.rel=="SMN")    
TAGS$Sex2=as.character(TAGS$Sex2)
TAGS$Sex2=with(TAGS,ifelse(Sex2%in%c("","?","U"),NA,ifelse(Sex2=="f","F",ifelse(Sex2=="m","M",Sex2))))
TAGS$Species2=as.character(TAGS$Species2)

TAGS$ReleaseLongitude2=with(TAGS,
            ifelse(ReleaseLatitude2<=(-23.74)& ReleaseLatitude2>(-24.5)& ReleaseLongitude2>113.8,113.2038,
                   ReleaseLongitude2))
TAGS$ReleaseLatitude2=with(TAGS,  
            ifelse(ReleaseLatitude2<=(-34.1)& ReleaseLatitude2>(-34.73)& 
                     ReleaseLongitude2>117.94 & ReleaseLongitude2<118.68,-34.83,
                   ReleaseLatitude2))


#4. --Detection manipulations ---

#Remove duplications
Detections$DUPLI=with(Detections,paste(TagCode,DateTime.local,Species,Sex))
Detections=Detections[!duplicated(Detections$DUPLI),]
Detections=Detections[,-match("DUPLI",colnames(Detections))]


  #4.1. fix dates and times and other variables
Detections$DateTime.local=as.POSIXlt(paste(Detections$Date.local,Detections$Time.local))
Detections$Date.local=as.POSIXlt(Detections$Date.local)
Detections$Time.local=times(Detections$Time.local)
Detections$ReleaseDate=as.POSIXlt(Detections$ReleaseDate)

fn=function(x)as.POSIXlt(x) #& Date.local< fn("2013-05-04")

# add missing depths
Detections$Depth=with(Detections,
            ifelse(is.na(Depth) & SerialNumber==120751 ,14,
            ifelse(is.na(Depth) & SerialNumber==101838 ,14.4,
            ifelse(is.na(Depth) & SerialNumber==105724 ,12,
            ifelse(is.na(Depth) & SerialNumber==101707 ,14,
            ifelse(is.na(Depth) & SerialNumber==101807 ,30,
            ifelse(is.na(Depth) & SerialNumber==101861 ,11,
            ifelse(is.na(Depth) & SerialNumber==101824 ,19,
            ifelse(is.na(Depth) & SerialNumber==105729 ,15,
            ifelse(is.na(Depth) & SerialNumber==105730 ,7,
            ifelse(is.na(Depth) & SerialNumber==101781 ,12.5,
            ifelse(is.na(Depth) & SerialNumber==106860 ,13,
            ifelse(is.na(Depth) & SerialNumber==101861 ,11,
            ifelse(is.na(Depth) & SerialNumber==106881 ,9,
            ifelse(is.na(Depth) & SerialNumber==106894 ,89,
            ifelse(is.na(Depth) & SerialNumber==101771 ,2,
            ifelse(is.na(Depth) & SerialNumber==120741 ,30,
            ifelse(is.na(Depth) & SerialNumber==101684 ,30,
            ifelse(is.na(Depth) & SerialNumber==113945 ,15.8,
            ifelse(is.na(Depth) & SerialNumber==106652 ,15,
            ifelse(is.na(Depth) & SerialNumber==106663 ,69,
            ifelse(is.na(Depth) & SerialNumber== 101796,16,
            ifelse(is.na(Depth) & SerialNumber== 105722,10,
            ifelse(is.na(Depth) & SerialNumber== 106897,14,
            ifelse(is.na(Depth) & SerialNumber== 120758,11,
            ifelse(is.na(Depth) & SerialNumber== 104490,148,
            ifelse(is.na(Depth) & SerialNumber== 107452,10,
            ifelse(is.na(Depth) & SerialNumber== 119865,14,
            Depth))))))))))))))))))))))))))))

# remove dates prior to first tagging event
dodgy.date.time=subset(Detections,DateTime.local<as.POSIXlt("2010-01-01 14:55:52"))
Detections=subset(Detections,DateTime.local>as.POSIXlt("2010-01-01 14:55:52"))


# add FL to Detections
if(Source.FL.externally=="NO")    
{
  TAGS=subset(TAGS,Species2%in%unique(Detections$Species))
  TAGS$DATE=strptime(TAGS$ReleaseDate2,format="%Y-%m-%d")
  #TAGS$DATE=strptime(TAGS$ReleaseDate2,format="%d-%b-%y")
  #if(class(TAGS$ReleaseDate2)=="factor") TAGS$DATE=as.POSIXlt(as.character(TAGS$ReleaseDate2),format="%d-%b-%y")
  TAGS$ReleaseLength=as.numeric(TAGS$ReleaseLength)
  TAGS$ReleaseLength=with(TAGS,ifelse(ReleaseLength==189,1.89,
            ifelse(ReleaseLength>5,ReleaseLength/100,ReleaseLength)))
  TAGS$ReleaseLongitude2=with(TAGS,ifelse(ReleaseLongitude2>160 & Code2==31019,113.7134,ReleaseLongitude2)) #fix dodgy release long
  Detections=merge(Detections,TAGS[,match(c("Species2","Code2","ReleaseLength","Project.rel"),names(TAGS))],
          by.x=c("TagCode","Species"),by.y=c("Code2","Species2"),all.x=T)
  Detections$FL=Detections$ReleaseLength
  Detections=Detections[,-match("ReleaseLength",names(Detections))]
}
  
if(Source.FL.externally=="YES")
{
  DATA$Species=as.character(DATA$Species)
  DATA$Species=with(DATA,ifelse(Species=="BW","Dusky",ifelse(Species=="TK","Thickskin",
       ifelse(Species=="GM","Gummy",ifelse(Species=="WH","Whiskery",
       ifelse(Species=="CP","bronze whaler",Species))))))
  DATA$FL=with(DATA,ifelse(Species=="BW" & TagCode==29531,229,FL))
  DATA$DATE=as.POSIXlt(as.character(DATA$DATE))
  DATA$Year.rel=DATA$DATE$year+1900
  DATA$Month.rel=DATA$DATE$mon+1
  DATA$Month.rel=with(DATA,ifelse(Species=="BW" & TagCode%in%c(30998,29534),8,Month.rel))
  DATA$Month.rel=with(DATA,ifelse(Species=="TK" & TagCode%in%c(31030,31029,
                                                               31027,31037,31033,31031,31028),7,Month.rel))
  DATA=DATA[,-match("DATE",names(DATA))]
  Detections=merge(Detections,DATA[,match(c("Species","TagCode","FL","Project.rel"),names(DATA))],by=c("TagCode","Species"),all.x=T)
  
}


  #separate copper from other species
TAGS.copper=subset(TAGS,Species2=="bronze whaler")
TAGS=subset(TAGS,Species2%in%these.species) 


  #4.2 add needed variables
Detections$Species=as.character(Detections$Species)

  #4.3 Unify different receivers in same station       
Detections$Station=paste(Detections$Latitude,Detections$Longitude)


  #4.4 Separate sentinel tags 
Detections.sentinel=subset(Detections,Species=="SENTINEL")
Detections=subset(Detections,!(Species=="SENTINEL"))


  #4.5 Fix dodgy longitude and latitudes
Detections$Longitude=with(Detections,
                        ifelse(is.na(Longitude) & SerialNumber==104490 &
                          Date.local>fn("2013-11-21"),115.26013,
                        ifelse(is.na(Longitude) & SerialNumber==113126 &
                          Date.local<fn("2013-11-21"),113.9107,
                        ifelse(is.na(Longitude) & SerialNumber==113140 &
                          Date.local<fn("2013-11-21"),113.9047,
                        Longitude))))
Detections$Latitude=with(Detections,
                        ifelse(is.na(Latitude) & SerialNumber==104490 &
                          Date.local>fn("2013-11-21"),-32.02237,
                        ifelse(is.na(Latitude) & SerialNumber==113126 &
                          Date.local<fn("2013-11-21"),-21.8848,
                        ifelse(is.na(Latitude) & SerialNumber==113140 &
                          Date.local<fn("2013-11-21"),-21.8813,
                        Latitude))))

Detections$ReleaseLongitude=with(Detections,ifelse(ReleaseLongitude=="",NA,ReleaseLongitude))
Detections$ReleaseLatitude=with(Detections,ifelse(ReleaseLatitude=="",NA,ReleaseLatitude))

Detections$ReleaseLongitude=as.numeric(Detections$ReleaseLongitude)
Detections$ReleaseLatitude=as.numeric(Detections$ReleaseLatitude)

Detections$ReleaseLatitude=with(Detections,ifelse(ReleaseLatitude>0,-ReleaseLatitude,ReleaseLatitude))

Detections$ReleaseLongitude=with(Detections, 
        ifelse(ReleaseLongitude>114 & ReleaseLongitude<114.5 & ReleaseLatitude>(-25) & ReleaseLatitude<(-24),113.15,ReleaseLongitude))


#tags fixed here
Detections$ReleaseLongitude=with(Detections,ifelse(TagCode==31019,113.7134,ifelse(TagCode==49146,138.24,
                        ifelse(TagCode==29531,118.342,ReleaseLongitude))))
Detections$ReleaseLatitude=with(Detections,ifelse(TagCode==49146,-34.56,
                        ifelse(TagCode==29531,-34.843,ReleaseLatitude)))

Dodgy.release.position=subset(Detections,is.na(ReleaseLatitude)|is.na(ReleaseLongitude) |ReleaseLongitude>150)
Dodgy.release.position=subset(Dodgy.release.position,!Project.rel=="South.Australia")
if(nrow(Dodgy.release.position)>0)
{
  n=table(Dodgy.release.position$Project.rel)
  if(n[1]>0)cat("check release locations for these tags",unique(Dodgy.release.position$TagCode))
}
  
Dodgy.receiver.position=subset(Detections,Latitude>0 |Longitude>150|is.na(Latitude)|is.na(Longitude))
if(nrow(Dodgy.receiver.position)>0) cat("check receiver locations for these codes",unique(Dodgy.receiver.position$TagCode))

  


  #4.6 Add geographical area for plotting
Detections=subset(Detections,!is.na(Latitude))
Detections$Area=with(Detections,ifelse(Latitude>-29,"North.WA",ifelse(Longitude>129,"SA",
                      ifelse(Latitude<-29 & Longitude<=129,"South.WA",NA))))
unicas.Areas=unique(Detections$Area)
unicas.Areas=unicas.Areas[!is.na(unique(Detections$Area))]
unicas.Areas=unicas.Areas[match(c("North.WA","South.WA","SA"),unicas.Areas)]
Detections$Area.release=with(Detections,ifelse(ReleaseLatitude>-29,"North.WA",
                      ifelse(ReleaseLongitude>129,"SA",ifelse(ReleaseLatitude<-29 & ReleaseLongitude<=129,"South.WA",NA))))
Detections$Array=with(Detections,ifelse(Latitude>(-24),"Ningaloo",
                  ifelse(Latitude<=(-24) & Latitude>(-33),"Perth",
                  ifelse(Latitude<=(-33),"Southern.lines",NA))))

 
  #4.7 select species of interest

write.csv(table(Detections$Species,useNA='ifany'),"total.hits.species.csv",row.names=F)

  # list of bronze and dusky South Australian sharks detected in WA
South.Oz.shks.detect=subset(Detections,Project.rel=='South.Australia' & Longitude<=129) 
South.Oz.shks.detect$TagCode.original=South.Oz.shks.detect$TagCode
South.Oz.shks.list=sort(unique(South.Oz.shks.detect$TagCode.original))
write.csv(South.Oz.shks.list,"South.Oz.shks.list.csv",row.names=F)

Other.detections=subset(Detections,!Species%in%these.species)
Detections.copper=subset(Other.detections,Species=="bronze whaler")
Detections=subset(Detections,Species%in%these.species)


  #4.8 Testing and removing strange hits
table(Detections$Area.release,Detections$Area)

#sandbar released in Ningaloo and one hit in SA
STRANGE=subset(Detections,TagCode==31038)
STRANGE=STRANGE[order(STRANGE$Date.local),]

#remove strange hit
id=which(Detections$TagCode==31038 & Detections$Area =="SA")
Detections=Detections[-id,]

#reset species name for reused tags
Rep.Recap$Species.rec=with(Rep.Recap,ifelse(SPECIES=="BW","Dusky",
                      ifelse(SPECIES=="GM","Gummy",
                      ifelse(SPECIES=="WH","Whiskery",
                      ifelse(SPECIES=="TK","Thickskin",as.character(SPECIES))))))
Rep.Recap$TagCode=Rep.Recap$ATAG_NO

Detections=merge(Detections,subset(Rep.Recap,select=c(TagCode,Species.rec,DATE_CAPTR)),by=c("TagCode"),all.x=T)
Detections$Species.old=Detections$Species

#correct Gummy= 29525 post recapture dates
Detections$Species=with(Detections,ifelse(TagCode==29525 & Date.local>DATE_CAPTR,"Unknwon",Species))

#add recapture information as last detections for recaptured sharks
b=Rep.Recap[,match(c("TagCode","Species.rec","SEX","DATE_CAPTR","RECLATDECDEG",
            "RECLNGDECDEG","RELLATDECDEG","RELLNGDECDEG","RELEASE DATE","FL","DATE_CAPTR"),names(Rep.Recap))]
b=subset(b,Species.rec%in%unique(Detections$Species))
names(b)=c("TagCode","Species","Sex","Date.local","Latitude","Longitude",
           "ReleaseLatitude","ReleaseLongitude","ReleaseDate","FL","DATE_CAPTR")
b$SerialNumber=b$Depth=b$DateTime.local=b$Time.local=b$Station=b$Array=b$Project.rel=NA
b$Project="SMN"
b$Area=with(b,ifelse(Latitude>-29,"North.WA",ifelse(Longitude>129,"SA",
      ifelse(Latitude<-29 & Longitude<=129,"South.WA",NA))))

x=subset(Detections,TagCode%in%unique(b$TagCode),select=c(Species,TagCode,Area.release,
        Species.rec,Species.old))
x=x[!duplicated(paste(x$Species,x$TagCode)),]
b=merge(b,x,by=c("Species","TagCode"),all.x=T)
b=b[,match(names(Detections),names(b))]
b$Area.release=with(b,ifelse(ReleaseLatitude>-29,"North.WA",
    ifelse(ReleaseLongitude>129,"SA",ifelse(ReleaseLatitude<-29 & ReleaseLongitude<=129,"South.WA",NA))))
b$Area.release=with(b,ifelse(is.na(Area.release)&TagCode==49146,"SA",Area.release))

Detections$Recapture.hit="NO"
b$Recapture.hit="YES"
Detections=rbind(Detections,b)

Detections$Array=with(Detections,ifelse(Latitude>(-24),"Ningaloo",
    ifelse(Latitude<=(-24) & Latitude>(-33),"Perth",
    ifelse(Latitude<=(-33),"Southern.lines",NA))))

Detections=subset(Detections,!is.na(TagCode))


  #4.9 Create dummy TagCode for confidentiallity issues  
TAGS$TagCode.original=TAGS$Code2
Detections$TagCode.original=Detections$TagCode
Detections=Detections[,-match("TagCode",names(Detections))]

TAGS=TAGS[order(TAGS$Species2,TAGS$DATE),]
unk=c("Dusky","Gummy","Thickskin","Whiskery")
names(unk)=c("DS","GS","SS","WS")
STR=vector('list',length(unk))
for(i in 1:length(unk))
{
  b=subset(TAGS,Species2==unk[i])
  b$TagCode=paste(names(unk)[i],'.',1:nrow(b),sep='')
  STR[[i]]=b
}
TAGS=do.call(rbind,STR)
TAGS$Species=TAGS$Species2
dummy=subset(TAGS,TagCode.original%in%Detections$TagCode.original,
             select=c(TagCode.original,TagCode,Species))
Detections=merge(Detections,dummy,by=c("TagCode.original","Species"),all.x=T)

Detections=subset(Detections,!is.na(TagCode.original))


#Add datetime for recaptures
Detections$DateTime.local=as.character(Detections$DateTime.local)
Detections$DateTime.local=with(Detections,ifelse(is.na(DateTime.local) & DATE_CAPTR==Date.local,
                               paste(Date.local,"00:00:00"),DateTime.local))
Detections$DateTime.local=as.POSIXlt(Detections$DateTime.local)


#Export data for Dave Jacoby
DaveJ=subset(Detections,Species%in%c("Dusky","Thickskin") & Project.rel=="SMN",
        select=c(TagCode,TagCode.original,SerialNumber,Species,Sex,FL,ReleaseDate,ReleaseLatitude,
                 ReleaseLongitude,Latitude,Longitude,Depth,DateTime.local,
                 Date.local,Time.local))
write.csv(DaveJ,handl_OneDrive("Data/Tagging/Acoustic_tagging/Dave Jacoby/Data.csv"),row.names=F)
write.csv(Receivers,handl_OneDrive("Data/Tagging/Acoustic_tagging/Dave Jacoby/Receivers.csv"),row.names=F)



  #4.10 Calculate straight line distance movement

    #4.10.1. movement (in km) between consecutive detections 
Detections=subset(Detections, !(is.na(Latitude) | is.na(Longitude)))
N.det=nrow(Detections)
Detections=Detections[order(Detections$Species,Detections$TagCode,Detections$DateTime.local),]

Detections$Latitude.prev=NA
Detections$Latitude.prev[1]=Detections$Latitude[1]
Detections$Latitude.prev[2:N.det]=Detections$Latitude[1:(N.det-1)]

Detections$Longitude.prev=NA
Detections$Longitude.prev[1]=Detections$Longitude[1]
Detections$Longitude.prev[2:N.det]=Detections$Longitude[1:(N.det-1)]



Detections$TagCode.prev=NA
Detections$TagCode.prev[1]=Detections$TagCode[1]
Detections$TagCode.prev[2:N.det]=Detections$TagCode[1:(N.det-1)]   

Detections$SerialNumber.prev=NA
Detections$SerialNumber.prev[1]=Detections$SerialNumber[1]
Detections$SerialNumber.prev[2:N.det]=Detections$SerialNumber[1:(N.det-1)]

Detections=Detections[order(Detections$Species,Detections$TagCode,Detections$DateTime.local),]
Detections$Dist.moved.conseq.det=with(Detections,
      ifelse(TagCode==TagCode.prev & !(Longitude.prev==Longitude & Latitude.prev==Latitude),
                   distCosine(cbind(Longitude.prev,Latitude.prev),cbind(Longitude,Latitude))/1000,
      ifelse(TagCode==TagCode.prev & Longitude.prev==Longitude & Latitude.prev==Latitude,0,NA)))
Detections$Dist.moved.conseq.det=with(Detections,ifelse(is.na(Dist.moved.conseq.det) & TagCode==TagCode.prev & 
                    SerialNumber.prev==SerialNumber,0,Dist.moved.conseq.det))


Detections$Latitude.prev=with(Detections,ifelse(TagCode==TagCode.prev,Latitude.prev,NA))
Detections$Longitude.prev=with(Detections,ifelse(TagCode==TagCode.prev,Longitude.prev,NA))

#set first .prev to release location
all.tags=unique(Detections$TagCode)
for(i in 1:length(all.tags))
{
  id=which(Detections$TagCode==all.tags[i])[1]
  Detections$Latitude.prev[id]=Detections$ReleaseLatitude[id]
  Detections$Longitude.prev[id]=Detections$ReleaseLongitude[id]
}


    #4.10.2. movement (in km) and days between release and first detection 
Detections=Detections[order(Detections$Species,Detections$TagCode,Detections$DateTime.local),]
Detections$Dup=duplicated(paste(Detections$TagCode,Detections$Species))

Detections$ReleaseLatitude=with(Detections, ifelse(is.na(ReleaseLatitude)& TagCode.original==49136,
                   -34.51,ReleaseLatitude))
Detections$ReleaseLongitude=with(Detections, ifelse(is.na(ReleaseLongitude)& TagCode.original==49136,
                    138.07,ReleaseLongitude))

Detections$Dist.moved.rel.det=with(Detections,ifelse(Dup=="FALSE",
      distCosine(cbind(ReleaseLongitude,ReleaseLatitude),cbind(Longitude,Latitude)),NA))/1000

Detections$days.rel.det=with(Detections,ifelse(Dup=="FALSE",
                  as.numeric(Date.local-ReleaseDate),NA))/(24*3600)
ID.1=which(Detections$days.rel.det<0)  #re set the release day of records with negative day.rel.det
Detections[ID.1,"ReleaseDate"]=Detections[ID.1,"ReleaseDate"]-(24*3600)
Detections$days.rel.det=with(Detections,ifelse(Dup=="FALSE",
                  as.numeric(Date.local-ReleaseDate),NA))/(24*3600)

Detections=Detections[order(Detections$Species,Detections$TagCode,Detections$DateTime.local),]
Detections$TagCode.prev[1]=NA
Detections$Latitude.prev=with(Detections,ifelse(TagCode==TagCode.prev,Latitude.prev,NA))
Detections$Longitude.prev=with(Detections,ifelse(TagCode==TagCode.prev,Longitude.prev,NA))
Detections$SerialNumber.prev=with(Detections,ifelse(TagCode==TagCode.prev,SerialNumber.prev,NA))



#Extract year and month
Detections$Year=Detections$Date.local$year+1900
Detections$Month=Detections$Date.local$mon+1 
Detections$Day=Detections$Date.local$mday

Detections$Year.rel=Detections$ReleaseDate$year+1900
Detections$Month.rel=Detections$ReleaseDate$mon+1 
Detections$Day.rel=Detections$ReleaseDate$mday

#Remove faulty date stamps (AATAMS receivers...)
table((Detections$Year))
Detections=subset(Detections,Year>=2011)


#Only keep WA sharks detected in WA and SA and SA sharks detected in WA
Detections=subset(Detections,!(Project=="AATAMS" & Project.rel=="South.Australia"))
Detections.copper=subset(Detections.copper,!(Project=="AATAMS" & Project.rel=="South.Australia"))

#remove unknown species
Detections=subset(Detections,!Species=="Unknwon")

#create list of original and new tagIDs
TAg.list=Detections[!duplicated(Detections$TagCode.original),match(c("TagCode","TagCode.original"),names(Detections))]


#5. -- Summary table ---
Detections$FL=with(Detections,ifelse(FL>10,FL/100,FL))  #make sure all FL are in m

fn.summary.tbl=function(tgs,dtect,ORDR,what)
{
  dtect=subset(dtect,Longitude<=129)  #report only WA detections
  tgs.SA=subset(tgs,Project.rel=="South.Australia" & Code2%in%unique(dtect$TagCode))
  tgs=subset(tgs,Project.rel=="SMN")
  if(nrow(tgs.SA)>0) tgs=rbind(tgs,tgs.SA)
  
  #Released individuals
  N.tag=table(tgs$Species2)
  N.tag=N.tag[match(ORDR,names(N.tag))]
  mean.FL.tagged=aggregate(ReleaseLength~Species2,tgs,mean,na.rm=T)
  SD.FL.tagged=aggregate(ReleaseLength~Species2,tgs,sd,na.rm=T)
  Sx.ratio=with(subset(tgs,!Species2=="SENTINEL"),table(Sex2,Species2))
  Sx.ratio.tagged=paste("1:",round(Sx.ratio[1,]/Sx.ratio[2,],1),sep="")
  names(Sx.ratio.tagged)=colnames(Sx.ratio)
  
  #detected individuals (including those released in SA)
  only.recaptured=subset(dtect, Species%in% ORDR & Recapture.hit=="YES")
  only.recaptured=subset(dtect,TagCode.original%in%only.recaptured$TagCode.original)
  only.recaptured=table(only.recaptured$TagCode.original)
  only.recaptured=as.numeric(names(only.recaptured[only.recaptured==1]))
  N.det=subset(dtect, Species%in% ORDR & !TagCode.original%in%only.recaptured)
  N.det$Sex=as.character(N.det$Sex)
  N.det$Sex=with(N.det,ifelse(Sex%in%c("","?","U"),NA,ifelse(Sex=="f","F",ifelse(Sex=="m","M",Sex))))
  N.det$Dupli.Station=with(N.det,paste(Species,Station)) #number of receivers detecting sharks
  Rec.det.sp=with(N.det[!duplicated(N.det$Dupli.Station),],table(Station,Species))
  Rec.det.sp[Rec.det.sp>0]=1
  Rec.det.sp=colSums(Rec.det.sp)
  N.det.array=table(N.det$Array,N.det$Species) #number of detections by array
  N.det.total=colSums(N.det.array)
  
  #number of days monitored and detected per species
  shks=unique(N.det$TagCode)
  Total.time.monitored=vector('list',length(shks))
  Tim.mon=function(what)
  {
    Rel=what$ReleaseDate[1]
    fst=sort(what$Date.local)[1]
    lst=sort(what$Date.local)[length(what$Date.local)]
    return(data.frame(Species=unique(what$Species),TagCode=unique(what$TagCode),Release=Rel,first.detec=fst,last.detect=lst))
  }
  for (i in 1:length(shks))Total.time.monitored[[i]]=Tim.mon(subset(N.det,TagCode==shks[i]))
  Total.time.monitored=do.call(rbind,Total.time.monitored)  
  
    #number of days between release and last detection 
  if(what=="bronze")Total.time.monitored$days.mon=with(Total.time.monitored,as.numeric(last.detect-Release))
  if(!what=="bronze")Total.time.monitored$days.mon=with(Total.time.monitored,as.numeric(last.detect-Release)/(24*3600))
  Total.time.monitored$days.mon=with(Total.time.monitored,ifelse(Release==last.detect,1,days.mon))
  Max.day.mon=aggregate(days.mon~Species,Total.time.monitored,max)
  Min.day.mon=aggregate(days.mon~Species,Total.time.monitored,min)
  Mean.day.mon=aggregate(days.mon~Species,Total.time.monitored,mean)
  SD.day.mon=aggregate(days.mon~Species,Total.time.monitored,sd)
  
    #number of days between release and first detection
  #if(what=="bronze")Total.time.monitored$days.first.det=with(Total.time.monitored,as.numeric(first.detec-Release))
  #if(!what=="bronze")Total.time.monitored$days.first.det=with(Total.time.monitored,as.numeric(first.detec-Release)/(24*3600))
  Total.time.monitored$days.first.det=with(Total.time.monitored,as.numeric(first.detec-Release)/(24*3600))
  Total.time.monitored$days.first.det=with(Total.time.monitored,ifelse(Release==last.detect,1,days.first.det))
  Max.day.frst=aggregate(days.first.det~Species,Total.time.monitored,max)
  Min.day.frst=aggregate(days.first.det~Species,Total.time.monitored,min)
  Mean.day.frst=aggregate(days.first.det~Species,Total.time.monitored,mean)
  SD.day.frst=aggregate(days.first.det~Species,Total.time.monitored,sd)
 
    #number of days between first and last detection 
  if(what=="bronze")Total.time.monitored$days.det=with(Total.time.monitored,as.numeric(last.detect-first.detec))
  if(!what=="bronze")Total.time.monitored$days.det=with(Total.time.monitored,as.numeric(last.detect-first.detec)/(24*3600))
  Total.time.monitored$days.det=with(Total.time.monitored,ifelse(first.detec==last.detect,1,days.det))
  Max.day.det=aggregate(days.det~Species,Total.time.monitored,max)
  Min.day.det =aggregate(days.det~Species,Total.time.monitored,min)
  Mean.day.det=aggregate(days.det~Species,Total.time.monitored,mean)
  SD.day.det=aggregate(days.det~Species,Total.time.monitored,sd)
  
  
  N.det$dupli=with(N.det,paste(TagCode,Species))
  N.det=N.det[!duplicated(N.det$dupli),]
  mean.FL.det=aggregate(FL~Species,N.det,mean,na.rm=T)
  SD.FL.det=aggregate(FL~Species,N.det,sd,na.rm=T)
  Sx.ratio=with(N.det,table(Sex,Species))
  Sx.ratio.det=paste("1:",round(Sx.ratio[1,]/Sx.ratio[2,],1),sep="")
  names(Sx.ratio.det)=colnames(Sx.ratio)
  N.det=table(N.det$Species)
  
  
  #combine in table
  fn.mean.sd=function(a,b) paste(round(a,0),"(",round(b,0),")",sep="")
  fn.sort=function(a,b) b[match(names(a),names(b))]
  
  #release
  FL.tagged=mean.FL.tagged[,2]
  FL.tagged=fn.mean.sd(FL.tagged*100,SD.FL.tagged[,2]*100)
  names(FL.tagged)=mean.FL.tagged[,1]
  FL.tagged=fn.sort(N.tag,FL.tagged)
  Sx.ratio.tagged=fn.sort(N.tag,Sx.ratio.tagged)
  dummy=rbind(N.tag,FL.tagged,Sx.ratio.tagged)
  
  #detections
  FL.det=mean.FL.det[,2]
  FL.det=fn.mean.sd(FL.det*100,SD.FL.det[,2]*100)
  names(FL.det)=mean.FL.det[,1]
  FL.det=fn.sort(N.tag,FL.det)   
  Sx.ratio.det=fn.sort(N.tag,Sx.ratio.det)    
  N.det=fn.sort(N.tag,N.det)    
  Mn.det=Mean.day.det[,2]
  Mn.det=fn.mean.sd(Mn.det,SD.day.det[,2])
  names(Mn.det)=Mean.day.det[,1]
  Mn.det=fn.sort(N.tag,Mn.det)
  Max.det=Max.day.det[,2]
  names(Max.det)=Max.day.det[,1]
  Max.det=fn.sort(N.tag,Max.det)
  Min.det=Min.day.det[,2]
  names(Min.det)=Min.day.det[,1]
  Min.det=fn.sort(N.tag,Min.det)
  Mn.mon=Mean.day.mon[,2]
  Mn.mon=fn.mean.sd(Mn.mon,SD.day.mon[,2])
  names(Mn.mon)=Mean.day.mon[,1]
  Mn.mon=fn.sort(N.tag,Mn.mon)
  Max.mon=Max.day.mon[,2]
  names(Max.mon)=Max.day.mon[,1]
  Min.mon=Min.day.mon[,2]
  names(Min.mon)=Min.day.mon[,1]
  Max.mon=fn.sort(N.tag,Max.mon)
  Min.mon=fn.sort(N.tag,Min.mon)
  Mn.frst=Mean.day.frst[,2]
  Mn.frst=fn.mean.sd(Mn.frst,SD.day.frst[,2])
  names(Mn.frst)=Mean.day.frst[,1]
  Mn.frst=fn.sort(N.tag,Mn.frst) 
  Max.frst=Max.day.frst[,2]
  names(Max.frst)=Max.day.frst[,1]
  Min.frst=Min.day.frst[,2]
  names(Min.frst)=Min.day.frst[,1]
  Max.frst=fn.sort(N.tag,Max.frst)
  Min.frst=fn.sort(N.tag,Min.frst)
  if(ncol(N.det.array)==1)rownames(N.det.array)=paste("n.det.per.array_",rownames(N.det.array),sep="") 
  if(ncol(N.det.array)>1)
  {
    N.det.array=N.det.array[,match(names(N.tag),colnames(N.det.array))]
    rownames(N.det.array)=paste("n.det.per.array_",rownames(N.det.array),sep="")
  }
  
  N.det.total=fn.sort(N.tag,N.det.total)
  Rec.det.sp =fn.sort(N.tag,Rec.det.sp)  
  
  dummy.det=rbind(N.det,FL.det,Sx.ratio.det,Mn.frst,Max.frst,Min.frst,Mn.mon,Max.mon,Min.mon,
                 Mn.det,Max.det,Min.det,N.det.total,N.det.array,Rec.det.sp)
  dummy=rbind(dummy,dummy.det)
  
  return(list(dummy=dummy,Total.time.monitored=Total.time.monitored,only.recaptured=only.recaptured))  
  
}
a=fn.summary.tbl(TAGS,Detections,c("Dusky","Thickskin","Gummy","Whiskery"),"notbronze")
Table1=a$dummy
Dist.days.monitored=a$Total.time.monitored
only.recaptured=a$only.recaptured

Detections.copper$Recapture.hit="YES"
Detections.copper$TagCode.original=Detections.copper$TagCode
if(nrow(TAGS.copper)>0)
{
  a=fn.summary.tbl(TAGS.copper,Detections.copper,"bronze whaler","bronze")
  Table1.copper=a$dummy
  Table1=cbind(Table1,Table1.copper)
}

write.csv(Table1,"Table1.data_summary.csv")

#export word table
Table1=cbind(Variable=c("N° of individuals.tg","Mean FL (SD).tg","Sex ratio (male:female).tg",
    "N° of individuals.det","Mean FL (SD).det","Sex ratio (male:female).det","Mean (SD).fst.de",
    "Max.fst.de","Min.fst.de",  "Mean (SD).monitor","Max.monitor","Min.monitor",                
    "Mean (SD).last.de","Max.last.de","Min.last.de",
    "Total N° of detections","N° of detections in the Ningaloo array",
    "N° of detections in the Perth array","N° of detections in the Southern Lines array",
    "Total N° of receivers detecting individuals"),Table1)
if(nrow(TAGS.copper)>0)Scenarios.tbl(WD=getwd(),Tbl=Table1,Doc.nm="Table1.data_summary",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
              HEDR=c('Variable','Dusky','Thickskin','Gummy','Whiskery','bronze whaler'),
              HEDR.cols=c(1,1,1,1,1,1),HEDR2=NA,HEDR3=NA)


#Numbers released by area and time  
a=subset(TAGS,Project.rel=="SMN")
a.sp=c("Dusky","Thickskin","Gummy","Whiskery")
nm.a.sp=c("Dusky","Sandbar","Gummy","Whiskery")
fn.plot.rel=function(spec,nm)
{
  dat=subset(a,Species2==spec)
  plotMap(worldLLhigh, xlim=xlm,ylim=ylm,plt = c(.001, 1, 0.075, 1),
          col="grey95",xlab="",ylab="",axes=F,main="")
  dat$ReleaseLatitude2=ifelse(dat$ReleaseLatitude2>0,-dat$ReleaseLatitude2,dat$ReleaseLatitude2)
  dat$ReleaseLongitude2=with(dat,ifelse(ReleaseLongitude2>113.65 & ReleaseLatitude2<(-23.7)
                                        & ReleaseLatitude2>(-24.7),113.36,ReleaseLongitude2))
  dat$ReleaseLatitude2=with(dat,ifelse(ReleaseLongitude2>117.86 & ReleaseLongitude2<118.73 &
                                         ReleaseLatitude2<(-33.96) & ReleaseLatitude2>(-34.49),-35.02,ReleaseLatitude2))
  dat$ReleaseLatitude2=with(dat,ifelse(ReleaseLongitude2>114.66 &
                                         ReleaseLatitude2<(-25.43243) & ReleaseLatitude2>(-26.62),-21.332,ReleaseLatitude2))
  
  dat$Yr=substr(dat$DATE,1,4)
  Tab=table(dat$Yr)
  Yrs=names(Tab)
  CLl=c("black","grey99","grey30","grey75","grey55")
  for(x in 1:length(Yrs))
  {
    dd=subset(dat,Yr==Yrs[x])
    with(dd,points(ReleaseLongitude2,ReleaseLatitude2,bg=CLl[x],pch=21,cex=1.4))
  }
  legend('right',paste(Yrs," (n=",Tab,")",sep=""),bty='n',pch=rep(21,length(Yrs)),cex=1.25,pt.cex=1.75,pt.bg=CLl)
  legend('topleft',nm,bty='n',cex=1.8)
  box()
  axis(side = 1, at =round(xlm[1]):xlm[2], labels = F, tcl = -.25)
  axis(side = 2, at = round(ylm[1]):ylm[2], labels = F,tcl = -.25)
  axis(side = 1, at =seq(xlm[1],xlm[2],5), labels = F,tcl = -.5)
  axis(side = 2, at = seq(ylm[1],ylm[2],5), labels = F,tcl = -.5)
  
}
tiff(file="Map_releases.tiff",width = 1600, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mai=c(.1,.1,.3,1),oma=c(3,4,2,.1),mgp=c(1,.6,0))
for(i in 1:length(a.sp))
{
  fn.plot.rel(spec=a.sp[i],nm=nm.a.sp[i])
  if(i%in%c(2:4))  axis(side = 1, at =seq(xlm[1],xlm[2],5), labels = seq(xlm[1],xlm[2],5), tcl = -.5,las=1,cex.axis=1.2)
  if(i%in%c(1:2))  axis(side = 2, at = seq(ylm[1],ylm[2],5), labels = seq(ylm[1],ylm[2],5),tcl = -.5,las=2,cex.axis=1.2)
}
mtext("    Latitude (ºS)",side=2,line=2,las=3,cex=2,outer=T)
mtext("Longitude (ºE)",side=1,line=1,cex=2,outer=T)
dev.off()
rm(a)


#6. -- Proportion of time within zones ---
#note: this is used for Risk Assessment of dusky and sandbar and for estimating movement rates of gummy and whiskery
Prop.time=subset(Detections,select=c(TagCode.original,TagCode,TagCode.prev,ReleaseDate,
            DateTime.local,Date.local,Time.local,ReleaseLatitude,ReleaseLongitude,Latitude,
            Longitude,Latitude.prev,Longitude.prev,Recapture.hit,Species,Sex,Depth,Project.rel))
                 
    #Add fishing zones
Prop.time$zone=as.character(with(Prop.time,
   ifelse(Longitude>=116.5 & Latitude<=(-26),"Zone2",
   ifelse(Longitude<116.5 & Latitude<=(-33),"Zone1",
   ifelse(Latitude>(-33) & Latitude<=(-26) & Longitude<116.5,"West",
   ifelse(Latitude>(-26) & Longitude<114.833,"Closed.ningaloo",
   ifelse(Latitude>(-22) & Longitude>=114.833 & Longitude<123.75,"North",
   ifelse(Latitude>(-22) & Longitude>=123.75,"Joint",NA))))))))
Prop.time$zone=with(Prop.time,ifelse(Latitude>(-33) & 
    Latitude<=(-31) & Longitude>=114.8476 & Longitude<116,"Closed.metro",zone))


#set first .prev to release location
all.tags=unique(Prop.time$TagCode)
for(i in 1:length(all.tags))
{
  id=which(Prop.time$TagCode==all.tags[i])[1]
  Prop.time$Latitude.prev[id]=Prop.time$ReleaseLatitude[id]
  Prop.time$Longitude.prev[id]=Prop.time$ReleaseLongitude[id]
}


# Prop.time$Dist.moved.rel.det=with(Prop.time,ifelse(Recapture.hit=="YES" & is.na(Dist.moved.rel.det),
#       distCosine(cbind(ReleaseLongitude,ReleaseLatitude),cbind(Longitude,Latitude))/1000,Dist.moved.rel.det))
# Prop.time$Dist.moved.conseq.det=with(Prop.time,ifelse(!is.na(Dist.moved.rel.det),
#                                 Dist.moved.rel.det,Dist.moved.conseq.det))
Prop.time$zone.prev=as.character(with(Prop.time,
       ifelse(Longitude.prev>=116.5 & Latitude.prev<=(-26),"Zone2",
       ifelse(Longitude.prev<116.5 & Latitude.prev<=(-33),"Zone1",
       ifelse(Latitude.prev>(-33) & Latitude.prev<=(-26) & Longitude.prev<116.5,"West",
       ifelse(Latitude.prev>(-26) & Longitude.prev<114.833,"Closed.ningaloo",
       ifelse(Latitude.prev>(-22) & Longitude.prev>=114.833 & Longitude.prev<123.75,"North",
       ifelse(Latitude.prev>(-22) & Longitude.prev>=123.75,"Joint",NA))))))))
Prop.time$zone.prev=with(Prop.time,ifelse(Latitude.prev>(-33) & 
       Latitude.prev<=(-31) & Longitude.prev>=114.8476 & Longitude.prev<116,"Closed.metro",zone.prev))

Prop.time=Prop.time[order(Prop.time$TagCode,Prop.time$Date.local,Prop.time$Time.local,Prop.time$zone),]
Prop.time$Date.local.prev=c(NA,as.character(Prop.time$Date.local)[1:(nrow(Prop.time)-1)])
Prop.time$Date.local.prev=with(Prop.time,ifelse(TagCode==TagCode.prev,Date.local.prev,NA))

#set first .prev to release location
all.tags=unique(Prop.time$TagCode)
for(i in 1:length(all.tags))
{
  id=which(Prop.time$TagCode==all.tags[i])[1]
  Prop.time$Date.local.prev[id]=as.character(Prop.time$ReleaseDate[id])
}

Prop.time$Date.local.prev=as.POSIXlt(Prop.time$Date.local.prev)
#Prop.time$days.conseq.det=with(Prop.time,ifelse(TagCode==TagCode.prev,as.numeric(Date.local-Date.local.prev),NA))/(24*3600)
Prop.time$same.zone=with(Prop.time,ifelse(TagCode==TagCode.prev & !(zone==zone.prev),"N","Y"))
Prop.time$same.zone=with(Prop.time,ifelse(ReleaseLatitude==Latitude.prev & ReleaseLongitude==Longitude.prev
                                          & !(zone==zone.prev),"N",same.zone))

#remove duplicated dates
Prop.time$DupliDateZone=with(Prop.time,paste(zone,TagCode,Species,Date.local))
Prop.time=Prop.time[!duplicated(Prop.time$DupliDateZone),]

#Total time monitored
shks=unique(Detections$TagCode)
Total.time.monitored=vector('list',length(shks))
Tim.mon=function(what)
{
  Rel=what$ReleaseDate[1]
  fst=sort(what$Date.local)[1]
  lst=sort(what$Date.local)[length(what$Date.local)]
  return(data.frame(TagCode=unique(what$TagCode),Release=Rel,first.detec=fst,last.detect=lst))
}
for (i in 1:length(shks))Total.time.monitored[[i]]=Tim.mon(subset(Detections,TagCode==shks[i]))
Total.time.monitored=do.call(rbind,Total.time.monitored)  
Total.time.monitored$days.mon=with(Total.time.monitored,as.numeric(last.detect-Release)/(24*3600))
Total.time.monitored$days.mon=with(Total.time.monitored,ifelse(Release==last.detect,1,days.mon))

# time monitored by yr
Time.monitored.yr=vector('list',length(shks))
Tim.mon.by.yr=function(what)
{
  Yr.rel=as.numeric(substr(what$ReleaseDate[1],1,4))
  Yr.lst=max(what$Year)
  Each.yr=Yr.rel:Yr.lst
  dd=vector('list',length(Each.yr))
  TG=unique(what$TagCode)
  for(y in 1:length(Each.yr))
  {
    d=subset(what,Year==Each.yr[y])
    
    if(nrow(d)==0)
    {
      fst=what$ReleaseDate[1]
      lst=as.POSIXlt(paste(Each.yr[y],"-12-31",sep=""))
    }
    
    if(nrow(d)>0)
    {
      if(Each.yr[y]==Yr.rel)fst=what$ReleaseDate[1]
      if(Each.yr[y]>Yr.rel) fst=as.POSIXlt(paste(Each.yr[y],"-01-01",sep=""))
      
      if(Each.yr[y]<Yr.lst)   lst=as.POSIXlt(paste(Each.yr[y],"-12-31",sep=""))
      if(Each.yr[y]==Yr.lst) lst=sort(what$Date.local)[length(what$Date.local)]
    }
    
    
    dd[[y]]=data.frame(TagCode=TG,Year=Each.yr[y],first.detec=fst,last.detect=lst)
  }
  return(do.call(rbind,dd))
}

system.time(for (i in 1:length(shks))Time.monitored.yr[[i]]=Tim.mon.by.yr(subset(Detections,TagCode==shks[i])))

Time.monitored.yr=do.call(rbind,Time.monitored.yr)  
Time.monitored.yr$days.mon=with(Time.monitored.yr,as.numeric(last.detect-first.detec)/(24*3600))



shks=unique(Prop.time$TagCode)
Time.mon.zone=vector('list',length(shks))
Tim.mon.zn=function(SHK)
{  
  THIS=c("TagCode","Species","Sex","ReleaseDate","ReleaseLatitude","ReleaseLongitude","Latitude","Longitude",
         "Latitude.prev","Longitude.prev","Date.local","Time.local","Date.local.prev","zone","zone.prev","same.zone")
  dat1=subset(Prop.time,TagCode==SHK,select=THIS)
  dat1$same.location=with(dat1,ifelse(Latitude==Latitude.prev & Longitude==Longitude.prev,"Y","N"))
  dat1$Interpolate=with(dat1,ifelse(same.location=="Y","N",NA))
  
  Leave.zone=dat1
  
  if(nrow(Leave.zone)>0)
  {
    Store=vector('list',nrow(Leave.zone))
    for(s in 1:nrow(Leave.zone))
    {
      a=Leave.zone[s,]
      
      n=as.numeric(a$Date.local-a$Date.local.prev)
      if(n==0)n=1
      b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),c(Longitude,Latitude),n=n, addStartEnd=F))
      
      #replicate location the number of days if same location
      if(nrow(b)<n & a$same.location=="Y") b= matrix(rep(b,n),ncol=2,byrow=T)
      
      #Add corners if not same location
      if(a$same.location=="N")
      {
        if(a$zone.prev=="North" & a$zone%in%c("Closed.metro","Zone1","West"))
        {
          if(a$Longitude<=Cape.Leuwin[,1] & a$Latitude>=Cape.Leuwin[,2])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
            n2=distCosine(Exmouth,Shark.bay)/1000
            n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
            b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
            b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
          if(a$Longitude>Cape.Leuwin[,1] & a$Latitude<Cape.Leuwin[,2])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
            n2=distCosine(Exmouth,Shark.bay)/1000
            n3=distCosine(Shark.bay,Cape.Leuwin)/1000
            n4=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3+n4
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            n4=round(n*(n4/TN))
            
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
            b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
            b2=with(a,gcIntermediate(Shark.bay,Cape.Leuwin,n=n3, addStartEnd=F))
            b3=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n4, addStartEnd=F))
            b=rbind(b,b1,b2,b3)
          }
          if(a$Latitude<Shark.bay[,2] & a$Latitude>Cape.Leuwin[,2])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
            n2=distCosine(Exmouth,Shark.bay)/1000
            n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
            b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
            b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
        }
        if(a$zone.prev=="North" & a$zone=="Zone2")
        {
          if(a$Longitude<=Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
            n2=distCosine(Exmouth,Shark.bay)/1000
            n3=distCosine(Shark.bay,Cape.Leuwin)/1000
            n4=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3+n4
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            n4=round(n*(n4/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
            b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
            b2=gcIntermediate(Shark.bay,Cape.Leuwin,n=n3, addStartEnd=F)
            b3=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n4, addStartEnd=F))
            b=rbind(b,b1,b2,b3)
          }
          if(a$Longitude>Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
            n2=distCosine(Exmouth,Shark.bay)/1000
            n3=distCosine(Shark.bay,Cape.Leuwin)/1000
            n4=distCosine(Cape.Leuwin,Mid.point)/1000
            n5=distCosine(Mid.point,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3+n4+n5
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            n4=round(n*(n4/TN))
            n5=round(n*(n5/TN))
            
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
            b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
            b2=gcIntermediate(Shark.bay,Cape.Leuwin,n=n3, addStartEnd=F)
            
            b3=with(a,gcIntermediate(Cape.Leuwin,Mid.point,n=n4, addStartEnd=F))
            b4=with(a,gcIntermediate(Mid.point,c(Longitude,Latitude),n=n5, addStartEnd=F))
            b=rbind(b,b1,b2,b3,b4)
          }
        }
        if(a$zone.prev%in%c("Closed.ningaloo") & a$zone=="Zone2")
        {
          if(a$Longitude<=Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
            n2=distCosine(Shark.bay,Cape.Leuwin)/1000    
            n3=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
            b1=gcIntermediate(Shark.bay,Cape.Leuwin,n=n2, addStartEnd=F)
            b2=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
          if(a$Longitude>Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
            n2=distCosine(Shark.bay,Cape.Leuwin)/1000    
            n3=distCosine(Cape.Leuwin,Mid.point)/1000
            n4=distCosine(Mid.point,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3+n4
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            n4=round(n*(n4/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
            b1=gcIntermediate(Shark.bay,Cape.Leuwin,n=n2, addStartEnd=F)
            b2=with(a,gcIntermediate(Cape.Leuwin,Mid.point,n=n4, addStartEnd=F))
            b3=with(a,gcIntermediate(Mid.point,c(Longitude,Latitude),n=n4, addStartEnd=F))
            b=rbind(b,b1,b2,b3)
          }
        }
        if(a$zone.prev%in%c("Closed.ningaloo") & a$zone%in%c("Zone1","Closed.metro"))
        {
          if(a$Latitude>=Cape.Leuwin[,2])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
            n2=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          if(a$Latitude<Cape.Leuwin[,2])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
            n2=distCosine(Shark.bay,Cape.Leuwin)/1000
            n3=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Shark.bay,Cape.Leuwin,n=n2, addStartEnd=F))
            b2=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
         }
        if(a$zone.prev%in%c("Closed.metro") & a$zone=="Zone2")   
        {
          if(a$Longitude<=Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          if(a$Longitude>Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,Mid.point)/1000
            n3=distCosine(Mid.point,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Cape.Leuwin,Mid.point,n=n2, addStartEnd=F))
            b2=with(a,gcIntermediate(Mid.point,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
        }
        if(a$zone.prev%in%c("Closed.metro","West") & a$zone=="Zone1")   
        {
          if(a$Latitude<Cape.Leuwin[,2])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Cape.Leuwin,c(a$Longitude,a$Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
        }
        if(a$zone.prev%in%c("Closed.metro","Zone1") & a$zone%in%c("Closed.ningaloo"))
        {
          n=as.numeric(a$Date.local-a$Date.local.prev)
          n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
          n2=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
          TN=n1+n2
          n1=round(n*(n1/TN))
          n2=round(n*(n2/TN))
          b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
          b1=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n2, addStartEnd=F))
          b=rbind(b,b1)
        }
        if(a$zone.prev=="Zone1" & a$zone%in%c("Closed.ningaloo"))
        {
          if(a$Latitude.prev>=Cape.Leuwin[,2])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
            n2=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          if(a$Latitude.prev<Cape.Leuwin[,2])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,Shark.bay)/1000
            n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Cape.Leuwin,Shark.bay,n=n2, addStartEnd=F))
            b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }

        }
        if(a$zone.prev%in%c("Zone1") & a$zone=="Closed.metro")   
        {
          if(a$Latitude>Cape.Leuwin[,2])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Cape.Leuwin,c(a$Longitude,a$Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
        }
        if(a$zone.prev=="Zone1" & a$zone%in%c("Zone2"))     
        {
          if(a$Longitude.prev<Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Mid.point)/1000    #divide time proportional to distance
            n2=distCosine(Mid.point,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Mid.point,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Mid.point,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          
        }
        if(a$zone.prev=="Zone2" & a$zone%in%c("Closed.ningaloo"))     
        {
          if(a$Longitude.prev<=Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,Shark.bay)/1000
            n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=gcIntermediate(Cape.Leuwin,Shark.bay,n=n2, addStartEnd=F)
            b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
          if(a$Longitude.prev>Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Mid.point)/1000
            n2=distCosine(Mid.point,Cape.Leuwin)/1000    #divide time proportional to distance
            n3=distCosine(Cape.Leuwin,Shark.bay)/1000
            n4=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3+n4
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            n4=round(n*(n4/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Mid.point,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Mid.point,Cape.Leuwin,n=n2, addStartEnd=F))
            b2=gcIntermediate(Cape.Leuwin,Shark.bay,n=n3, addStartEnd=F)
            b3=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n4, addStartEnd=F))
            b=rbind(b,b1,b2,b3)
          }
        }
        if(a$zone.prev=="Zone2" & a$zone%in%c("Closed.metro"))     
        {
          if(a$Longitude.prev<=Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          if(a$Longitude.prev>Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Mid.point)/1000    #divide time proportional to distance
            n2=distCosine(Mid.point,Cape.Leuwin)/1000
            n3=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Mid.point,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Mid.point,Cape.Leuwin,n=n2, addStartEnd=F))
            b2=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
          
        }
        if(a$zone.prev=="Zone2" & a$zone%in%c("Zone1"))     
        {
          if(a$Longitude.prev>Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Mid.point)/1000    #divide time proportional to distance
            n2=distCosine(Mid.point,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Mid.point,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Mid.point,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          
        }
        if(a$zone.prev=="Zone2" & a$zone%in%c("North"))      
        {
          if(a$Longitude.prev<=Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,Shark.bay)/1000
            n3=distCosine(Shark.bay,Exmouth)/1000
            n4=distCosine(Exmouth,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3+n4
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            n4=round(n*(n4/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=gcIntermediate(Cape.Leuwin,Shark.bay,n=n2, addStartEnd=F)
            b2=gcIntermediate(Shark.bay,Exmouth,n=n3, addStartEnd=F)
            b3=with(a,gcIntermediate(Exmouth,c(Longitude,Latitude),n=n4, addStartEnd=F))
            b=rbind(b,b1,b2,b3)
          }
          if(a$Longitude.prev>Mid.point[,1])
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Mid.point)/1000    #divide time proportional to distance
            n2=distCosine(Mid.point,Cape.Leuwin)/1000
            n3=distCosine(Cape.Leuwin,Shark.bay)/1000
            n4=distCosine(Shark.bay,Exmouth)/1000
            n5=distCosine(Exmouth,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3+n4+n5
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            n4=round(n*(n4/TN))
            n5=round(n*(n5/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Mid.point,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Mid.point,Cape.Leuwin,n=n2, addStartEnd=F))
            b2=gcIntermediate(Cape.Leuwin,Shark.bay,n=n3, addStartEnd=F)
            b3=gcIntermediate(Shark.bay,Exmouth,n=n4, addStartEnd=F)
            b4=with(a,gcIntermediate(Exmouth,c(Longitude,Latitude),n=n5, addStartEnd=F))
            b=rbind(b,b1,b2,b3,b4)
          }
        }

      }
      
      x=a[rep(seq_len(nrow(a)), n),]
      this.na=match(c("Latitude.prev","Longitude.prev","same.zone"),names(x))
      x[,this.na]=NA
      if(!nrow(x)==nrow(b)) x=x[1:nrow(b),]
      x$Longitude=b[,1]
      x$Latitude=b[,2]
      Day.seq=seq(1,n,1)*24*60*60
      if(!length(Day.seq)==nrow(b)) Day.seq=Day.seq[1:nrow(b)]
      x$Date.local=a$Date.local.prev+Day.seq    
      Store[[s]]=x
    }    
    Store=do.call(rbind,Store)
    Store$zone=as.character(with(Store,
                                 ifelse(Longitude>=116.5 & Latitude<=(-26),"Zone2",
                                        ifelse(Longitude<116.5 & Latitude<=(-33),"Zone1",
                                               ifelse(Latitude>(-33) & Latitude<=(-26) & Longitude<116.5,"West",
                                                      ifelse(Latitude>(-26) & Longitude<114.833,"Closed.ningaloo",
                                                             ifelse(Latitude>(-22) & Longitude>=114.833 & Longitude<123.75,"North",
                                                                    ifelse(Latitude>(-22) & Longitude>=123.75,"Joint",NA))))))))
    Store$zone=with(Store,ifelse(Latitude>(-33) & 
                                   Latitude<=(-31) & Longitude>=114.8476 & Longitude<116,"Closed.metro",zone))
    Store$Interpolate="Y"
  }
  
  
  All=Store
  All$Dupl=with(All,paste(Date.local,zone))
  All=All[!duplicated(All$Dupl),]
  All=All[order(All$TagCode,All$Date.local),]
  
  if(nrow(All)>1)
  {
    All$Date.local.prev[2:nrow(All)]=as.character(All$Date.local)[1:(nrow(All)-1)]
    All$Date.local.prev=as.POSIXlt(All$Date.local.prev)
    All$zone.prev[2:nrow(All)]=All$zone[1:(nrow(All)-1)]
    All=subset(All,!is.na(Date.local))    
  }
  
  return(All)
}
system.time(for (i in 1:length(shks))Time.mon.zone[[i]]=Tim.mon.zn(shks[i]))
Time.mon.zone=do.call(rbind,Time.mon.zone)  

#Plot example of track
fn.plot.track=function(what,ylm,xlm)
{
  d=subset(Time.mon.zone,TagCode==what)
  e=subset(Prop.time,TagCode==what)  
  
  plotMap(worldLLhigh, xlim=xlm,ylim=ylm,plt = c(.001, 1, 0.075, 1),
          col="grey88",xlab="",ylab="",axes=F)
  
  plot(JA_Northern_Shark,ylim=c(-39,-11),xlim=c(108,130),add=T,col="grey75")
  plot(WA_Northern_Shark,add=T,col="grey45")
  plot(WA_Northern_Shark_2,add=T,col="grey90")
  plot(WCDGDLL,add=T,col="white")
  plot(SDGDLL_zone1,add=T,col="grey70")
  plot(SDGDLL_zone2,add=T,col="grey55")
  
  text(117.7227,-17.95964,"WANCSF")
  text(113,-20.5,"Ningaloo", srt=45)
  text(113,-21.5,"closure", srt=45)
  text(113.5,-31.45,"WCDGDLF")
  text(113.4807,-33.5,"JASDGDLF")
  text(113.4807,-34.25,"(Zone 1)")
  text(123.7437,-34.5,"JASDGDLF")
  text(123.7437,-35.25,"(Zone 2)")
  
  axis(side = 1, at =round(xlm[1]):xlm[2], labels = round(xlm[1]):xlm[2], tcl = .5,las=1,cex.axis=1.2)
  axis(side = 2, at = round(ylm[1]):ylm[2], labels = -(round(ylm[1]):ylm[2]),tcl = .5,las=2,cex.axis=1.2)
 
  points(d$Longitude,d$Latitude,pch=19,cex=0.8)
  points(e$Longitude,e$Latitude,col=3,pch=19,cex=2.5)
  points(e$ReleaseLongitude[1],e$ReleaseLatitude[1],col=2,pch=19,cex=2.5)
  if(what=="DS.87")
  {
    NNN=length(d$Longitude)
    ss=gcIntermediate(cbind(d$Longitude[NNN-1],d$Latitude[NNN-1]),
                      cbind(d$Longitude[NNN],d$Latitude[NNN]),n=8, addStartEnd=F)
    points(ss[,1],ss[,2],pch=19,cex=0.8)
  }
  
  mtext("    Latitude (ºS)",side=2,line=3.5,las=3,cex=1.8)
  mtext("Longitude (ºE)",side=1,line=3,cex=1.8)
  legend("center",paste(unique(d$Species)," shark"," (tag: ", what,")",sep=""),bty='n',cex=1.5)
  
}

id=TAg.list[which(TAg.list$TagCode.original==29583),]$TagCode
tiff(file="Track.example.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(oma = c(2, 1, 1, 0))
fn.plot.track(id,c(-36,-16),c(112,129))
box()
dev.off()

  #Animate Dusky movement
Animate.Dusky="NO"
if(Animate.Dusky=="YES")
{
  library(geosphere) 
  library(splancs)        #for testing if location is within boundaries
  library(sp)
  library(adehabitatLT)
  library(PBSmapping)   #for polygon
  data(worldLLhigh)
  library(CircStats)
  
  library(Rcpp)
  library(ggplot2)
  library(ggmap)
  require(animation) # NB, must install ImageMagick
  require(ggsn)   #for scale baer
  
  Animate.dusky=subset(Detections,Species=='Dusky' & ReleaseLongitude<124)
  Animate.dusky.traj=subset(Time.mon.zone,Species=='Dusky' & ReleaseLongitude<124)
  
  xrng=c(112, 124)
  yrng=c(-36, -19)
  
  p <- ggmap(get_map(c(116,-28),maptype="satellite", zoom = 5))
  p<-p+labs(x="Longitude",y="Latitude")+
    scale_x_continuous(limits = xrng, expand = c(0, 0))+
    scale_y_continuous(limits = yrng, expand = c(0, 0))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))+
     scalebar(dist = 250, st.size=4, height=0.04, dd2km = TRUE, model = 'WGS84',
             x.min=xrng[1],x.max=xrng[2]*.999,y.min=yrng[1]*.725,y.max=yrng[2])

  Show.what="arrows"
  
  Animate.dusky$Julian=as.numeric(Animate.dusky$Date.local-min(Animate.dusky$Date.local))/(60*60*24)
  Animate.dusky.traj$Julian=as.numeric(Animate.dusky.traj$Date.local-min(Animate.dusky.traj$Date.local))/(60*60*24)
  
  Animate.dusky$dummy=with(Animate.dusky,paste(Julian,TagCode))
  Animate.dusky.traj$dummy=with(Animate.dusky.traj,paste(Julian,TagCode))
  
  Animate.dusky=Animate.dusky[!duplicated(Animate.dusky$dummy),]
  Animate.dusky.traj=Animate.dusky.traj[!duplicated(Animate.dusky.traj$dummy),]
  
  
  N.traj.per.frame=10
  
  Animate.dusky$Julian.frame=N.traj.per.frame*floor(Animate.dusky$Julian/N.traj.per.frame)
  Animate.dusky.traj$Julian.frame=N.traj.per.frame*floor(Animate.dusky.traj$Julian/N.traj.per.frame)
  
  a=subset(Animate.dusky,select=c(TagCode,FL))
  a=a[!duplicated(a$TagCode),]
  Animate.dusky.traj=merge(Animate.dusky.traj,a,by="TagCode",all.x=T)
  Animate.dusky.traj$SIZE.cat=with(Animate.dusky.traj,ifelse(FL>2.15,"Large (>2 m)","Small (<2 m)"))
  Animate.dusky.traj$Sex=with(Animate.dusky.traj,ifelse(Sex=='U','F',Sex))
  Animate.dusky.traj$Sex=with(Animate.dusky.traj,ifelse(Sex=='F','Female','Male'))
  #Animate.dusky.traj$LGND=with(Animate.dusky.traj,paste(SIZE.cat,Sex,TagCode))
  Animate.dusky.traj$LGND=Animate.dusky.traj$TagCode
  colfunc.F <- colorRampPalette(c("hotpink4", "deeppink"))
  colfunc.M <- colorRampPalette(c("dodgerblue", "dodgerblue4"))
  
  #Full trajectories
  time=sort(unique(Animate.dusky.traj$Julian.frame))
  
      #add increasing time line
  Time.spn.lktn=data.frame(lon=c(115,121),lat=c(-24,-24))
  X.SEQ=seq(Time.spn.lktn$lon[1],Time.spn.lktn$lon[2],l=length(time))
  Time.spn.seq=time/max(time)
  fn.Time.spn=function(x)
  {
    path=data.frame(x= X.SEQ[1:match(x,time)],y=mean(Time.spn.lktn$lat))
    return(path)
  }
  p=p+geom_path(data = fn.Time.spn(x=time[length(time)]), aes(x, y),size=7,col="grey",alpha=0.5)+
    annotate('text', x = Time.spn.lktn$lon[1]*1.001, y=Time.spn.lktn$lat[1]*1.025,label=as.character(min(Animate.dusky.traj$Date.local)), col="white",size=4.5)+
    annotate('text', x = Time.spn.lktn$lon[2]*.9925, y=Time.spn.lktn$lat[2]*1.025,label=as.character(max(Animate.dusky.traj$Date.local)),col="white",size=4.5)

  system.time({saveGIF({
    N.shk=NULL
    suppressWarnings(for(x in time)
    {
      DF <- subset(Animate.dusky.traj,Julian.frame==x)
      N.shk=unique(c(N.shk,unique(DF$TagCode)))
      
      if(Show.what=="arrows")
      {
        print(p+geom_path(data = DF,mapping=aes(x = Longitude, y = Latitude,color=LGND),size=.75,
                          arrow=arrow(angle = 15,length=unit(0.1, "inches")),alpha=0.9,show.legend=F)+
              geom_path(data = fn.Time.spn(x=x), aes(x, y),size=7,col="white")+
                annotate('text', x = 116,y=-19.5,label=paste(length(N.shk),"sharks detected in total"), col="white",size=4.5))
      }
    })
  },movie.name=handl_OneDrive("Analyses/Acoustic_tagging/Dusky_migration/Animation_data_Trajectories.gif"),interval=0.25,loop =1)})   

}

  #aggregate number of days within each zone
Time.mon.zone$Days=with(Time.mon.zone,ifelse(zone==zone.prev,
              as.numeric(Date.local-Date.local.prev),NA))
#Time.mon.zone$Days=with(Time.mon.zone,ifelse(zone==zone.prev,
#       as.numeric(Date.local-Date.local.prev)/(24*3600),NA))
Time.mon.zone$Days=with(Time.mon.zone,ifelse(Date.local==Date.local.prev,1,Days))
Species.time.zone=aggregate(Days~zone+TagCode+Species,Time.mon.zone,sum,na.rm=T)
Species.time.zone=merge(Species.time.zone,
    Total.time.monitored[,match(c("TagCode","days.mon"),names(Total.time.monitored))],
    by="TagCode",all.x=T)
Species.time.zone$prop.time.in.zn=Species.time.zone$Days/Species.time.zone$days.mon
Species.time.zone=Species.time.zone[,-match("Days",names(Species.time.zone))]  

Species.time.zone=reshape(Species.time.zone,v.names = "prop.time.in.zn", idvar = c("TagCode","Species","days.mon"),
          timevar = "zone", direction = "wide")
ID=c("days.mon","prop.time.in.zn.Closed.metro","prop.time.in.zn.Closed.ningaloo",
  "prop.time.in.zn.West","prop.time.in.zn.Zone1","prop.time.in.zn.Zone2","prop.time.in.zn.North")
names(Species.time.zone)[match(ID,names(Species.time.zone))]=c("Tot.days.mon.",
          "Closed.metro","Closed.ningaloo","West","Zone1","Zone2","North")
Species.time.zone=Species.time.zone[order(Species.time.zone$Species,Species.time.zone$TagCode),
            match(c("Species","TagCode","North","Closed.ningaloo","Closed.metro",
                    "West","Zone1","Zone2","Tot.days.mon."),names(Species.time.zone))]  
Species.time.zone[is.na(Species.time.zone)]=0
IJ=match(c("North","Closed.ningaloo","Closed.metro","West","Zone1","Zone2"),names(Species.time.zone))
Species.time.zone[,IJ]=round(Species.time.zone[,IJ],2)

#add zone release
Time.mon.zone$Zone.rel=as.character(with(Time.mon.zone,
      ifelse(ReleaseLongitude>=116.5 & ReleaseLatitude<=(-26),"Zone2",
      ifelse(ReleaseLongitude<116.5 & ReleaseLatitude<=(-33),"Zone1",
      ifelse(ReleaseLatitude>(-33) & ReleaseLatitude<=(-26) & ReleaseLongitude<116.5,"West",
      ifelse(ReleaseLatitude>(-26) & ReleaseLongitude<114.833,"Closed.ningaloo",
      ifelse(ReleaseLatitude>(-22) & ReleaseLongitude>=114.833 & ReleaseLongitude<123.75,"North",
      ifelse(ReleaseLatitude>(-22) & ReleaseLongitude>=123.75,"Joint",NA))))))))
Time.mon.zone$Zone.rel=with(Time.mon.zone,ifelse(ReleaseLatitude>(-33) & 
      ReleaseLatitude<=(-31) & ReleaseLongitude>=114.8476 & ReleaseLongitude<116,"Closed.metro",Zone.rel))

dummy=subset(Time.mon.zone,select=c(TagCode,Species,Zone.rel))
dummy=dummy[!duplicated(paste(dummy$TagCode,dummy$Species,dummy$Zone.rel)),]
Species.time.zone=merge(Species.time.zone,dummy,by=c("TagCode","Species"),all.x=T)

#add dropped tagcodes for staying all time in same zone
aa=unique(Prop.time$TagCode)
bb=unique(Species.time.zone$TagCode)
a=which(!aa%in%bb)
b=subset(Prop.time,TagCode%in%aa[a])
if(nrow(b)>0)
{
  b=merge(b,dummy,by=c("TagCode","Species"),all.x=T)
  b1=Species.time.zone[1:nrow(b),]
  b1[,]=NA
  b1$TagCode=b$TagCode
  b1$Species=b$Species
  b1$North=with(b,ifelse(zone=="North" & Zone.rel=="North",1,0))
  b1$Closed.ningaloo=with(b,ifelse(zone=="Closed.ningaloo" & Zone.rel=="Closed.ningaloo",1,0))
  b1$Closed.metro=with(b,ifelse(zone=="Closed.metro" & Zone.rel=="Closed.metro",1,0))
  b1$West=with(b,ifelse(zone=="West" & Zone.rel=="West",1,0))
  b1$Zone1=with(b,ifelse(zone=="Zone1" & Zone.rel=="Zone1",1,0))
  b1$Zone2=with(b,ifelse(zone=="Zone2" & Zone.rel=="Zone2",1,0))
  b1$Tot.days.mon.=1
  b1$Zone.rel=b$Zone.rel
  b1=b1[!duplicated(paste(b1$TagCode,b1$Species)),]
  Species.time.zone=rbind(Species.time.zone,b1)
}
Species.time.zone=Species.time.zone[order(Species.time.zone$Species,Species.time.zone$Zone.rel,
    -Species.time.zone$Tot.days.mon.),match(c("Species","Tot.days.mon.","TagCode","Zone.rel","North",
    "Closed.ningaloo","West" ,"Closed.metro","Zone1","Zone2"),colnames(Species.time.zone))]

#export
write.csv(Species.time.zone,"Prop.time.by.zone.csv",row.names=F)

#plot proportions
fn.plot.prop=function(what,sort.criteria)
{
  a=subset(Species.time.zone,Species==what)
  id=match(Zns,names(a))
  a[,id]=a[,id]/rowSums(a[,id])
  crap=as.data.frame(as.matrix(COL.prop))
  names(crap)="Col"
  crap$Zone.rel=rownames(crap)
  a=merge(a,crap,by="Zone.rel")
  a$Sort=with(a,ifelse(Zone.rel=="North",6,
                       ifelse(Zone.rel=="Closed.ningaloo",5,
                              ifelse(Zone.rel=="West",4,
                                     ifelse(Zone.rel=="Closed.metro",3,
                                            ifelse(Zone.rel=="Zone1",2,1))))))
  find.zn.rl=match(a$Zone.rel,colnames(a))
  Dummy=find.zn.rl
  for(dd in 1:length(Dummy))Dummy[dd]=a[dd,find.zn.rl[dd]]
  a$p.same=Dummy
  if(sort.criteria=="Zone_only")a=a[order(a$Sort),]
  if(sort.criteria=="Zone_prop")a=a[order(a$Sort,a$p.same),]
  par(mai=c(.55,.5,.25,.65),oma = c(.4, .5, 1, 0),xpd=T,mgp=c(1,.55,0))
  r=barplot(t(a[,id]), col = COL.prop,horiz=T,beside=F,yaxt='n',cex.axis=1.25)
  legend("top",Zns.leg,bty='n',pt.cex=2,pch=15,col=COL.prop,horiz=T,inset=c(0,-.07),
         cex=1)
  box()
  points(rep(-0.045,length(r)),r,pch=15,cex=1.5,col=as.character(a$Col))
  #axis(2,r,a$Zone.rel,las=1,cex.axis=1.25)
  mtext("Proportion of time",1,line=2,cex=1.75)
  mtext("Release zone",2,line=1.6,cex=1.75)
  axis(4,r,(a$Tot.days.mon.),las=1,cex.axis=0.56)
  mtext("Days monitored",4,line=2.1,cex=1.75)
}

tiff(file="Proportion.zone.dusky.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
#fn.plot.prop("Dusky",sort.criteria="Zone_only")
fn.plot.prop("Dusky",sort.criteria="Zone_prop")
dev.off()

tiff(file="Proportion.zone.whiskery.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
fn.plot.prop("Whiskery",sort.criteria="Zone_prop")
dev.off()

tiff(file="Proportion.zone.gummy.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
fn.plot.prop("Gummy",sort.criteria="Zone_prop")
dev.off()

tiff(file="Proportion.zone.sandbar.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
fn.plot.prop("Thickskin",sort.criteria="Zone_prop")
dev.off()


#Dummy tag list of original tag codes
TAG.LST=subset(Detections,select=c(TagCode.original,TagCode,Species))
TAG.LST=TAG.LST[!duplicated(TAG.LST$TagCode.original),]



#6.1-- Movement rates of gummy and whiskery among zones ----
if (do_mov_rates=="here")
{
  #Mov.rates="annual"
  Mov.rates="daily.matrix.expand"     #switch movement rates on/off
  
  
  if(Mov.rates=="daily.matrix.expand")
  {
    #Drop sharks monitored for less than one month
    drop.tag.cod.mvov.rate=subset(Species.time.zone,Tot.days.mon.>=30,select=c(Species,TagCode)) 
    drop.tag.cod.mvov.rate=subset(drop.tag.cod.mvov.rate,Species%in%c("Gummy","Whiskery"))
    
    shks=unique(drop.tag.cod.mvov.rate$TagCode)
    Move.rate=subset(Prop.time,TagCode%in%shks)
    Move.rate$Block=with(Move.rate,paste(floor(-Latitude),floor(Longitude),sep=""))
    Move.rate$Block.prev=with(Move.rate,paste(floor(-Latitude.prev),floor(Longitude.prev),sep=""))
    Move.rate$Block.rel=with(Move.rate,paste(floor(-ReleaseLatitude),floor(ReleaseLongitude),sep=""))
    
    #Calculation shark position per day
    #note: time step= 1 day as sharks cover less than 2 degrees in one day (Figure Speed distributions)
    #      assume straight line movement between detections and correct for corners
    Time.mon.grid=vector('list',length(shks))
    Tim.mon.grid=function(SHK)
    {  
      dat1=subset(Move.rate,TagCode==SHK)
      dat1$Latitude.prev[1]=dat1$ReleaseLatitude[1]
      dat1$Longitude.prev[1]=dat1$ReleaseLongitude[1]
      dat1$Block.prev[1]=dat1$Block.rel[1]
      dat1$Date.local.prev[1]=dat1$ReleaseDate[1]
      
      THIS=c("TagCode","Species","Sex","ReleaseDate","Latitude","Longitude","Latitude.prev",
             "Longitude.prev","Date.local","Date.local.prev","Block","Block.prev","same.blk","zone.prev","zone")
      dat1$same.blk=with(dat1,ifelse(Block==Block.prev,"Y","N"))
      
      dat1$Dummy=with(dat1,paste(Block,Date.local))
      dat1=dat1[!duplicated(dat1$Dummy),]
      dat1=dat1[order(dat1$Date.local),]
      
      Stay.in.zone=subset(dat1,same.blk=="Y"| is.na(same.blk),select=THIS)
      Leave.zone=subset(dat1,same.blk=="N",select=THIS)
      
      #add intermediate days for those that stay in grid
      if(nrow(Stay.in.zone)>0)
      {
        Store.stay=vector('list',nrow(Stay.in.zone))
        for(s in 1:nrow(Stay.in.zone))
        {
          a=Stay.in.zone[s,]
          n=as.numeric(a$Date.local-a$Date.local.prev)
          if(n>1)
          {
            a<- a[rep(row.names(a), n), ]
            Day.seq=seq(1,n,1)*24*60*60
            if(!length(Day.seq)==nrow(a)) Day.seq=Day.seq[1:nrow(a)]
            a$Date.local=a$Date.local.prev+Day.seq  
            a$Block.prev=a$Block
            a$Date.local.prev=as.POSIXlt(c(NA,as.character(a$Date.local[1:(nrow(a)-1)])))
          }
          Store.stay[[s]]=a
          
        }
        Store.stay=do.call(rbind,Store.stay)
        Stay.in.zone=Store.stay
      }
      
      #interpolate position and day for those that left grid
      if(nrow(Leave.zone)>0)
      {
        Store=vector('list',nrow(Leave.zone))
        for(s in 1:nrow(Leave.zone))
        {
          a=Leave.zone[s,]
          n=as.numeric(a$Date.local-a$Date.local.prev)
          b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),c(Longitude,Latitude),n=n, addStartEnd=F))
          x=a
          if(nrow(b)>1)
          {
            #Add corners
            if(a$zone.prev=="North" & a$zone%in%c("Closed.metro","Zone1","West"))
            {
              n=as.numeric(a$Date.local-a$Date.local.prev)
              n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
              n2=distCosine(Exmouth,Shark.bay)/1000
              n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
              TN=n1+n2+n3
              n1=round(n*(n1/TN))
              n2=round(n*(n2/TN))
              n3=round(n*(n3/TN))
              
              b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
              b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
              b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
              b=rbind(b,b1,b2)
            }
            if(a$zone.prev=="North" & a$zone=="Zone2")
            {
              n=as.numeric(a$Date.local-a$Date.local.prev)
              n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
              n2=distCosine(Exmouth,Shark.bay)/1000
              n3=distCosine(Shark.bay,Cape.Leuwin)/1000
              n4=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
              TN=n1+n2+n3+n4
              n1=round(n*(n1/TN))
              n2=round(n*(n2/TN))
              n3=round(n*(n3/TN))
              n4=round(n*(n4/TN))
              b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
              b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
              b2=gcIntermediate(Shark.bay,Cape.Leuwin,n=n3, addStartEnd=F)
              b3=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n4, addStartEnd=F))
              b=rbind(b,b1,b2,b3)
            }
            if(a$zone.prev%in%c("Closed.ningaloo") & a$zone=="Zone2")
            {
              n=as.numeric(a$Date.local-a$Date.local.prev)
              n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
              n2=distCosine(Shark.bay,Cape.Leuwin)/1000    
              n3=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
              TN=n1+n2+n3
              n1=round(n*(n1/TN))
              n2=round(n*(n2/TN))
              n3=round(n*(n3/TN))
              b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
              b1=gcIntermediate(Shark.bay,Cape.Leuwin,n=n2, addStartEnd=F)
              b2=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n3, addStartEnd=F))
              b=rbind(b,b1,b2)
            }
            if(a$zone.prev%in%c("Closed.ningaloo") & a$zone%in%c("Zone1","Closed.metro"))
            {
              n=as.numeric(a$Date.local-a$Date.local.prev)
              n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
              n2=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
              TN=n1+n2
              n1=round(n*(n1/TN))
              n2=round(n*(n2/TN))
              b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
              b1=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n2, addStartEnd=F))
              b=rbind(b,b1)
            }
            if(a$zone.prev%in%c("Closed.metro") & a$zone=="Zone2")
            {
              n=as.numeric(a$Date.local-a$Date.local.prev)
              n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
              n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
              TN=n1+n2
              n1=round(n*(n1/TN))
              n2=round(n*(n2/TN))
              b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
              b1=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n2, addStartEnd=F))
              b=rbind(b,b1)
            }
            if(a$zone.prev=="Zone2" & a$zone%in%c("Closed.ningaloo"))
            {
              n=as.numeric(a$Date.local-a$Date.local.prev)
              n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
              n2=distCosine(Cape.Leuwin,Shark.bay)/1000
              n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
              TN=n1+n2+n3
              n1=round(n*(n1/TN))
              n2=round(n*(n2/TN))
              n3=round(n*(n3/TN))
              b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
              b1=gcIntermediate(Cape.Leuwin,Shark.bay,n=n2, addStartEnd=F)
              b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
              b=rbind(b,b1,b2)
            }
            if(a$zone.prev=="Zone2" & a$zone%in%c("Closed.metro"))
            {
              n=as.numeric(a$Date.local-a$Date.local.prev)
              n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
              n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
              TN=n1+n2
              n1=round(n*(n1/TN))
              n2=round(n*(n2/TN))
              b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
              b1=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n2, addStartEnd=F))
              b=rbind(b,b1)
            }
            if(a$zone.prev=="Zone2" & a$zone%in%c("North"))
            {
              n=as.numeric(a$Date.local-a$Date.local.prev)
              n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
              n2=distCosine(Cape.Leuwin,Shark.bay)/1000
              n3=distCosine(Shark.bay,Exmouth)/1000
              n4=distCosine(Exmouth,c(a$Longitude,a$Latitude))/1000
              TN=n1+n2+n3+n4
              n1=round(n*(n1/TN))
              n2=round(n*(n2/TN))
              n3=round(n*(n3/TN))
              n4=round(n*(n4/TN))
              b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
              b1=gcIntermediate(Cape.Leuwin,Shark.bay,n=n2, addStartEnd=F)
              b2=gcIntermediate(Shark.bay,Exmouth,n=n3, addStartEnd=F)
              b3=with(a,gcIntermediate(Exmouth,c(Longitude,Latitude),n=n4, addStartEnd=F))
              b=rbind(b,b1,b2,b3)
            }
            if(a$zone.prev=="Zone1" & a$zone%in%c("Closed.ningaloo"))
            {
              n=as.numeric(a$Date.local-a$Date.local.prev)
              n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
              n2=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
              TN=n1+n2
              n1=round(n*(n1/TN))
              n2=round(n*(n2/TN))
              b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
              b1=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n2, addStartEnd=F))
              b=rbind(b,b1)
            }
            x=a[rep(seq_len(nrow(a)), n),]
            this.na=match(c("Latitude.prev","Longitude.prev","same.blk"),names(x))
            x[,this.na]=NA
            if(!nrow(x)==nrow(b)) x=x[1:nrow(b),]
            x$Longitude=b[,1]
            x$Latitude=b[,2]
            Day.seq=seq(1,n,1)*24*60*60
            if(!length(Day.seq)==nrow(b)) Day.seq=Day.seq[1:nrow(b)]
            x$Date.local=a$Date.local.prev+Day.seq       
            x$Latitude.prev= c(NA,x$Latitude[1:(nrow(x)-1)])
            x$Longitude.prev=c(NA,x$Longitude[1:(nrow(x)-1)])
            x$Date.local.prev=as.POSIXlt(c(NA,as.character(x$Date.local[1:(nrow(x)-1)])))
            x$Block=with(x,paste(floor(-Latitude),floor(Longitude),sep=""))
            x$Block.prev=with(x,paste(floor(-Latitude.prev),floor(Longitude.prev),sep=""))
            
          }
          
          Store[[s]]=x
        }  
        Store=do.call(rbind,Store)
      }
      All=NULL
      if(nrow(Leave.zone)==0) All=Stay.in.zone
      if(nrow(Leave.zone)>0) All=rbind(Store,Stay.in.zone)
      All=All[order(All$Date.local),match(c("TagCode","Species","Sex","Date.local","Date.local.prev","Block","Block.prev",
                                            "Latitude","Longitude"),names(All))]
      All$dummy= c(NA,All$Block[1:(nrow(All)-1)])
      All$Block.prev=with(All,ifelse(Block.prev=="NANA",dummy,Block.prev))  
      All=subset(All,!is.na(Date.local.prev))
      
      return(All[,-match('dummy',names(All))])
    }
    system.time(for (i in 1:length(shks))Time.mon.grid[[i]]=Tim.mon.grid(shks[i]))
    Time.mon.grid=do.call(rbind,Time.mon.grid)  
    
    
    source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Deviance.explained.R"))
    move.rate.fn=function(dat,MOD,zns,set.non.adj.0)
    {
      #Put data in right format (add pseudo-0)
      dat$Block=as.numeric(dat$Block)
      dat$Block.prev=as.numeric(dat$Block.prev)
      dat=subset(dat,!Block%in%c(33116,34117))
      dat=subset(dat,!Block.prev%in%c(33116,34117))  #drop land blocks
      dat$Indx=paste(dat$Block.prev,dat$TagCode)
      uni.blk.rel=unique(dat$Indx)
      dat$TagCode=as.character(dat$TagCode)
      dat$Species=as.character(dat$Species)
      dat$Sex=as.character(dat$Sex)    
      fn.add.adj=function(blk)
      {
        NW=blk-1001;N=blk-1000;NE=N+1;W=blk-1;SW=W+1000;Self=blk;S=blk+1000;SE=S+1;E=blk+1
        blks=c(NW,N,NE,W,SW,Self,S,SE,E)
        names(blks)=c("NW","N","NE","W","SW","Self","S","SE","E")
        STore.add.stuff=vector('list',length(a$Block))
        for(n in 1:length(a$Block))
        {
          add.0=blks[which(!blks%in%a$Block[n])]
          add.stuff= a[rep(1, length(add.0)), ]
          add.stuff$Coordinate=names(add.0)
          add.stuff$Recapture=0
          STore.add.stuff[[n]]=add.stuff
        }
        add.stuff=do.call(rbind,STore.add.stuff)
        
        blks=data.frame(Block=blks,Coordinate=names(blks))
        a=merge(a,blks,by="Block",all.x=T)
        a$Recapture=1
        a=rbind(a,add.stuff)
        return(a[,match(c("TagCode","Block.prev","Coordinate","Recapture"),names(a))])
      }  
      Store=vector('list',length(uni.blk.rel))
      for(i in 1:length(uni.blk.rel))
      {
        a=subset(dat,Indx==uni.blk.rel[i])
        BLK=unique(a$Block.prev)
        Store[[i]]=fn.add.adj(BLK)
      }
      dat1=do.call(rbind,Store)
      
      #Observations
      TT=unique(dat1$Block.prev)
      Collect=vector('list',length(TT))
      for( s in 1:length(TT))
      {
        a=subset(dat1,Block.prev==TT[s])
        NewData=a[!duplicated(paste(a$Block.prev,a$Coordinate)),]
        TAB=table(a$Recapture,a$Coordinate)
        test=round(TAB[2,]/sum(TAB[2,]),3)
        test=data.frame(OBS=test,Coordinate=names(test))
        NewData=merge(NewData,test,by="Coordinate")
        Collect[[s]]=NewData
      }
      Observed.props=do.call(rbind,Collect)
      
      #Model predictions
      dat1$Coordinate=as.factor(dat1$Coordinate)
      dat1$Block.prev=as.factor(dat1$Block.prev)
      dat1$TagCode=as.factor(dat1$TagCode)
      
      dat1$Lat=-as.numeric(substr(dat1$Block.prev,1,2))
      dat1$Lon=as.numeric(substr(dat1$Block.prev,3,6))
      
      dat1$Zone=with(dat1,ifelse(Lat<=(-26) & Lat>(-33) &Lon<117,"West",
                                 ifelse(Lat<=(-33) &Lon<117,"Zone1",
                                        ifelse(Lat<(-30) & Lon>=117,"Zone2",NA))))
      dat1$Zone=as.factor(dat1$Zone)
      
      
      #if(MOD=="GLM")prob <- glm(Recapture ~ Coordinate+Coordinate:Zone/Block.prev,data =dat1,family=binomial)
      if(MOD=="GLM")prob <- glm(Recapture ~ Coordinate+Coordinate:Block.prev,data =dat1,family=binomial)
      if(MOD=="GLMM") prob <- glmer(Recapture ~ Coordinate+Coordinate:Block.prev+(1 | TagCode),
                                    data =dat1, family = binomial, control = glmerControl(optimizer = "bobyqa"),nAGQ = 1)
      
      dat1$Preds=predict(prob,type='response')
      if(class(prob)[1]=="glm")Dev.exp=Dsquared(prob,adjust = T)
      
      Predicted.props=dat1[!duplicated(paste(dat1$Block.prev,dat1$Coordinate)),]
      
      #Compare observations with prediction
      #   par(mfcol=c(6,4))
      #   for( s in 1:length(TT))
      #   {
      #     a=subset(Observed.props,Block.prev==TT[s])
      #     b=subset(Predicted.props,Block.prev==TT[s])
      #     plot(a$Coordinate,a$OBS,ylim=c(0,1))
      #     points(b$Coordinate,b$Preds,col=2,cex=2)
      #   }
      #   
      #Create square matrix
      #1st allocate block to coordinates
      names(Predicted.props)[match("Block.prev",names(Predicted.props))]="Block.from"
      Predicted.props$Block.from=as.numeric(as.character(Predicted.props$Block.from))
      Predicted.props$Block.to=with(Predicted.props,
                                    ifelse(Coordinate=="NW",Block.from-1001,
                                           ifelse(Coordinate=="N",Block.from-1000, 
                                                  ifelse(Coordinate=="W",Block.from-1,
                                                         ifelse(Coordinate=="Self",Block.from,
                                                                ifelse(Coordinate=="S",Block.from+1000,
                                                                       ifelse(Coordinate=="E",Block.from+1,
                                                                              ifelse(Coordinate=="NE",Block.from-999,  
                                                                                     ifelse(Coordinate=="SW",Block.from+999,  
                                                                                            ifelse(Coordinate=="SE",Block.from+1001,NA))))))))))
      
      #2nd remove land blocks and rescale props
      Land.blk=c(as.numeric(paste(30,116:128,sep='')),as.numeric(paste(31,116:127,sep='')),
                 as.numeric(paste(32,116:123,sep='')),as.numeric(paste(33,116:119,sep='')),34117)
      Too.deep.blk=c(29113,29114,29115,30113,31113,32113,33113,34113,35113,35114,3415,34125,34126,35119,
                     35120,35121,35122,35123,35124,paste(36,113:128,sep=''))
      Predicted.props=subset(Predicted.props,!Block.to%in%c(Land.blk,Too.deep.blk))
      uni=unique(Predicted.props$Block.from)
      dummy=vector('list',length(uni))
      for(u in 1:length(uni))
      {
        aa=subset(Predicted.props,Block.from==uni[u])
        aa$Preds=aa$Preds/sum(aa$Preds)
        dummy[[u]]=aa
      }
      Predicted.props=do.call(rbind,dummy)
      
      #3rd predict for Block.from with no data for square matrix
      BLKS.to=unique(Predicted.props$Block.to)
      BLKS.frm=unique(Predicted.props$Block.from)
      no.dat=BLKS.to[which(!BLKS.to%in%BLKS.frm)]
      Add.this=vector('list',length(no.dat))
      for(bb in 1:length(no.dat))
      {
        ww=BLKS.frm[which(abs(BLKS.frm-no.dat[bb])==min(abs(BLKS.frm-no.dat[bb])))]
        ww=subset(Predicted.props,Block.from==ww)
        ww$Block.from=no.dat[bb]
        ww$Block.to=with(ww,
                         ifelse(Coordinate=="NW",Block.from-1001,
                                ifelse(Coordinate=="N",Block.from-1000, 
                                       ifelse(Coordinate=="W",Block.from-1,
                                              ifelse(Coordinate=="Self",Block.from,
                                                     ifelse(Coordinate=="S",Block.from+1000,
                                                            ifelse(Coordinate=="E",Block.from+1,
                                                                   ifelse(Coordinate=="NE",Block.from-999,  
                                                                          ifelse(Coordinate=="SW",Block.from+999,  
                                                                                 ifelse(Coordinate=="SE",Block.from+1001,NA))))))))))
        ww=subset(ww,!Block.to%in%c(Land.blk,Too.deep.blk))
        ww$Preds=ww$Preds/sum(ww$Preds)
        Add.this[[bb]]=ww
      }
      Add.this=do.call(rbind,Add.this)
      
      Predicted.props=rbind(Predicted.props,Add.this)
      Predicted.props=subset(Predicted.props,!Block.to%in%c(32127,33127))
      Predicted.props=Predicted.props[order(Predicted.props$Block.from,Predicted.props$Block.to),]
      Reshaped=reshape(subset(Predicted.props,select=c(Block.from,Preds,Block.to)),
                       v.names = "Preds", idvar = "Block.from",timevar = "Block.to", direction = "wide")
      names(Reshaped)[2:length(names(Reshaped))]=substr(names(Reshaped)[2:length(names(Reshaped))],7,20)
      
      Reshaped[is.na(Reshaped)]=0
      MTM.day=as.matrix(Reshaped[,2:ncol(Reshaped)])
      rownames(MTM.day)=Reshaped$Block.from
      MTM.day=MTM.day[,match(rownames(MTM.day),colnames(MTM.day))]
      
      
      ##...expand one day movement transition matrix
      #weekly
      MAT7=MTM.day%*%MTM.day%*%MTM.day%*%MTM.day%*%MTM.day%*%MTM.day%*%MTM.day
      MAT7=MAT7/rowSums(MAT7)
      
      #monthly
      MAT30=MAT7%*%MAT7%*%MAT7%*%MAT7
      MAT30=MAT30/rowSums(MAT30)
      
      #annual
      MAT365=MAT30%*%MAT30%*%MAT30%*%MAT30%*%MAT30%*%MAT30%*%MAT30%*%MAT30%*%MAT30%*%MAT30%*%MAT30%*%MAT30
      MAT365=MAT365/rowSums(MAT365)
      
      
      ##...aggregate from grid to regions  
      
      #1. define grids per zone
      RegnMat=data.frame(Grid=rownames(MAT365))
      RegnMat$Lat=-as.numeric(substr(RegnMat$Grid,1,2))
      RegnMat$Lon=as.numeric(substr(RegnMat$Grid,3,6))
      RegnMat$Zone=with(RegnMat,ifelse(Lat<=(-26) & Lat>(-33) &Lon<117,"West",
                                       ifelse(Lat<=(-33) &Lon<117,"Zone1",
                                              ifelse(Lat<(-30) & Lon>=117,"Zone2",NA))))
      RegnMat$dummy=1
      RegnMat=RegnMat[,-match(c("Lat","Lon"),names(RegnMat))]
      RegnMat2=reshape(RegnMat,v.names="dummy",idvar="Zone",timevar="Grid",direction="wide") 
      RegnMat=reshape(RegnMat,v.names="dummy",idvar="Grid",timevar="Zone",direction="wide")
      
      RegnMat[is.na(RegnMat)]=0
      names(RegnMat)[2:ncol(RegnMat)]=substr(names(RegnMat)[2:ncol(RegnMat)],7,15) 
      rownames(RegnMat)=RegnMat$Grid
      RegnMat=RegnMat[,match(zns,names(RegnMat))]
      colnames(RegnMat)=names(RegnMat)
      RegnMat=as.matrix(RegnMat)
      
      RegnMat2=RegnMat2[order(RegnMat2$Zone),]
      RegnMat2[is.na(RegnMat2)]=0
      RegnMat2[,2:ncol(RegnMat2)]=RegnMat2[,2:ncol(RegnMat2)]/rowSums(RegnMat2[,2:ncol(RegnMat2)])
      names(RegnMat2)[2:ncol(RegnMat2)]=substr(names(RegnMat2)[2:ncol(RegnMat2)],7,15)  
      rownames(RegnMat2)=RegnMat2$Zone
      RegnMat2=RegnMat2[,-1]
      colnames(RegnMat2)=names(RegnMat2)
      RegnMat2=as.matrix(RegnMat2)
      
      #2. aggregate columns and rows into regions
      fn.ag=function(dat)
      {
        MAT.zn=dat
        MAT.zn=MAT.zn%*%RegnMat
        MAT.zn=RegnMat2%*%MAT.zn  
        if(set.non.adj.0=="YES") MAT.zn[1,3]=MAT.zn[3,1]=0  #rescale to 0 prob of moving to non-adjacent zone
        MAT.zn=MAT.zn/rowSums(MAT.zn)  
        return(MAT.zn)
      }
      MAT30.zn=fn.ag(MAT30)    #monthly
      MAT365.zn=fn.ag(MAT365) #annual
      
      return(list(MTM.day=MTM.day,MAT7=MAT7,MAT30=MAT30,MAT365=MAT365,
                  MAT30.zn=MAT30.zn,MAT365.zn=MAT365.zn))
    }
    
    MOV.RT=list(Gummy=NA,Whiskery=NA)
    sp.vr.rt=c('Gummy','Whiskery')
    Zn.species=list(c("West","Zone1","Zone2"),c("Zone1","Zone2"))   #no detections in West for Whiskery
    Set.0s=list("NO","NO")
    
    #note: GLMM takes >2 days... still undone...
    for(y in 1:length(sp.vr.rt)) MOV.RT[[y]]=move.rate.fn(dat=subset(Time.mon.grid,Species==sp.vr.rt[y]),
                                                          MOD="GLM",zns=Zn.species[[y]],set.non.adj.0=Set.0s[[y]])
    
    #display matrix (origin is the rows, destination is the columns)
    fn.show.MTM=function(MATRIZ,TITLE,Scale)
    {
      N.int=50
      colfunc <- colorRampPalette(c("navy", "cadetblue","white"))
      couleurs=rev(colfunc(N.int))
      BREAKS=seq(0,1,length.out=N.int+1)
      xx=1:nrow(MATRIZ)
      yy=1:ncol(MATRIZ)
      image(xx,yy,MATRIZ,ylab="",xlab="",xaxt='n',yaxt='n',col =couleurs,breaks=BREAKS)
      axis(1,xx,F,tck=0.025)
      axis(2,yy,F,tck=0.025)
      mtext(TITLE,3,.75,cex=2)
      box()
      SQ=rev(seq(BREAKS[1],BREAKS[length(BREAKS)],.25))
      if(Scale=='block')
      {
        axis(1,seq(1,nrow(MATRIZ),1),colnames(MATRIZ)[seq(1,nrow(MATRIZ),1)],las=2,cex.axis=.75,hadj=0.6,tck=0.025)
        axis(2,seq(1,nrow(MATRIZ),1),(colnames(MATRIZ)[seq(1,nrow(MATRIZ),1)]),cex.axis=.75,tck=0.025,hadj=.6,las=1)
        color.legend(xx[length(xx)*.94],yy[length(yy)*.1],xx[length(xx)*.99],yy[length(yy)*.5],
                     SQ,rect.col=rev(couleurs),gradient="y",col=1,cex=.8)
        mtext("Origin grid",1,2.5,cex=1.5)
        mtext("Destination grid",2,3,cex=1.5)
      }
      if(Scale=='zone')
      {
        color.legend(3,1,3.35,2,SQ,rect.col=rev(couleurs),gradient="y",col=1,cex=1.25)
        axis(1,1:nrow(MATRIZ),colnames(MATRIZ)[1:nrow(MATRIZ)],las=2,cex.axis=1.45,hadj=0.8,tck=0.025)
        axis(2,1:nrow(MATRIZ),colnames(MATRIZ)[1:nrow(MATRIZ)],cex.axis=1.45,tck=0.025,hadj=.8,las=1)
        mtext("Origin grid",1,5,cex=1.95)
        mtext("Destination grid",2,5,cex=1.95)
      }
    }
    
    for(y in 1:length(sp.vr.rt))
    {
      hndl=paste("Movement_rates/",sp.vr.rt[y],sep='')
      This=MOV.RT[[y]]
      fn.tif=function(FILE) tiff(file=FILE,width=2400,height=2400,units="px",res=300,compression="lzw")
      
      fn.tif(paste(hndl,".MTM.day.blk.tiff",sep=''))
      fn.show.MTM(This$MTM.day,"Daily exchange rate among blocks",'block')
      dev.off()
      
      fn.tif(paste(hndl,".MTM.week.blk.tiff",sep=''))
      fn.show.MTM(This$MAT7,"Weekly exchange rate among blocks",'block')
      dev.off()
      
      fn.tif(paste(hndl,".MTM.mon.blk.tiff",sep=''))
      fn.show.MTM(This$MAT30,"Monthly exchange rate among blocks",'block')
      dev.off()
      
      fn.tif(paste(hndl,".MTM.year.blk.tiff",sep=''))
      fn.show.MTM(This$MAT365,"Annual exchange rate among blocks",'block')
      dev.off()
      
      fn.tif(paste(hndl,".MTM.year.zn.tiff",sep=''))
      par(mfcol=c(1,1),mai=c(1.3,1.3,.5,.1))
      fn.show.MTM(This$MAT365.zn,paste(names(MOV.RT)[y],"shark"),'zone')
      dev.off()
      
      write.csv(This$MAT365.zn,paste(hndl,"shark.Annual.zone.MTM.csv"),row.names=T)
    }
    
    
    
    #Plot management zones and receiver lines for TDGDLF
    source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/TDGLDF.zones.receivers.R"))
    tiff(file="TDGDLF_zone_receiver.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    fn.map.zones.receivers(a=112:129,PLATE=c(.01,.9,.075,.9),OZ.lat=c(-44.5,-11),OZ.long=c(113,155),
                           South.WA.lat=c(-36,-25), South.WA.long=c(112,129))
    dev.off()
    
  }
  
  
  #Simon's approach
  if(Mov.rates=="annual")
  {
    #Drop sharks monitored for less than one month
    drop.tag.cod.mvov.rate=subset(Species.time.zone,Tot.days.mon.>=30,select=c(Species,TagCode)) 
    drop.tag.cod.mvov.rate=subset(drop.tag.cod.mvov.rate,Species%in%c("Gummy","Whiskery"))
    
    shks=unique(drop.tag.cod.mvov.rate$TagCode)
    Move.rate=subset(Prop.time,TagCode%in%shks)
    
    #Drop first 30 days to reduce release effect
    n.filter=30
    Move.rate$Delta.t=with(Move.rate,ifelse(Recapture.hit=="NO",as.numeric(Date.local-ReleaseDate),NA))/(24*3600)
    Move.rate=subset(Move.rate,Delta.t>n.filter | is.na(Delta.t))
    shks=unique(Move.rate$TagCode)
    
    Move.rate$Block=with(Move.rate,paste(floor(-Latitude),floor(Longitude),sep=""))
    Move.rate$Block.prev=with(Move.rate,paste(floor(-Latitude.prev),floor(Longitude.prev),sep=""))
    Move.rate$Block.rel=with(Move.rate,paste(floor(-ReleaseLatitude),floor(ReleaseLongitude),sep=""))
    
    Move.rate$Lat=-as.numeric(substr(Move.rate$Block,1,2))
    Move.rate$Lon=as.numeric(substr(Move.rate$Block,3,6))
    Move.rate$zone=with(Move.rate,ifelse(Lat<=(-26) & Lat>(-33) &Lon<117,"West",
                                         ifelse(Lat<=(-33) &Lon<117,"Zone1",
                                                ifelse(Lat<(-30) & Lon>=117,"Zone2",NA))))
    
    Move.rate$N=1
    Move.rate$Year=as.POSIXlt(as.character(Move.rate$Date.local))$year+1900
    Move.rate$Month=as.POSIXlt(as.character(Move.rate$Date.local))$mon+1 
    Move.rate$Day=as.POSIXlt(as.character(Move.rate$Date.local))$mday
    
    #aggregate by month
    agg.month.year=aggregate(N~Year+Month+TagCode+Species,Move.rate,sum)
    names(agg.month.year)[match("N",names(agg.month.year))]="Total"  
    agg.zone.month.year=aggregate(N~Year+Month+zone+TagCode+Species,Move.rate,sum)
    agg.zone.month.year=merge(agg.zone.month.year,agg.month.year,
                              by=c("Year","Month","TagCode","Species"),all.x=T)
    agg.zone.month.year$prop=agg.zone.month.year$N/agg.zone.month.year$Total
    
    #aggregate by year
    agg.year=aggregate(N~Year+TagCode+Species,Move.rate,sum)
    names(agg.year)[match("N",names(agg.year))]="Total"  
    agg.zone.year=aggregate(N~Year+zone+TagCode+Species,Move.rate,sum)
    agg.zone.year=merge(agg.zone.year,agg.year,
                        by=c("Year","TagCode","Species"),all.x=T)
    agg.zone.year$prop=agg.zone.year$N/agg.zone.year$Total
    
    
    #Detection area by zone for weighting
    Number.Stations.Line=data.frame(N=c(30,25,28,48,44,33),Line=c("Rot.out","Rot.in","Perth.coast",
                                                                  "Hamelin","South2","Albany"),zone=c("West","West","West","Zone1","Zone1","Zone2"))
    Number.Stations.Line=aggregate(N~zone,Number.Stations.Line,sum)
    Number.Stations.Line$Det.area=detec.range*Number.Stations.Line$N
    Number.Stations.Line$Rel.det.area=Number.Stations.Line$Det.area/min(Number.Stations.Line$Det.area)
    
    agg.zone.month.year=merge(agg.zone.month.year,subset(Number.Stations.Line,select=c(zone,Rel.det.area)),by="zone",all.x=T)
    agg.zone.year=merge(agg.zone.year,subset(Number.Stations.Line,select=c(zone,Rel.det.area)),by="zone",all.x=T)
    
    setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement/Movement_rates/Simon"))
    write.csv(Move.rate,"rawdata.csv",row.names=F)
    write.csv(agg.zone.year,"agg.zone.year.csv",row.names=F)
    
    interpolate.between.detections="NO"
    if(interpolate.between.detections=="YES")
    {
      #Calculation shark position per day
      #note: time step= 1 day as sharks cover less than 2 degrees in one day (Figure Speed distributions)
      #      assume straight line movement between detections and correct for corners
      Time.mon.grid=vector('list',length(shks))
      Tim.mon.grid=function(SHK)
      {  
        dat1=subset(Move.rate,TagCode==SHK)
        dat1$Latitude.prev[1]=dat1$ReleaseLatitude[1]
        dat1$Longitude.prev[1]=dat1$ReleaseLongitude[1]
        dat1$Block.prev[1]=dat1$Block.rel[1]
        dat1$Date.local.prev[1]=dat1$ReleaseDate[1]+(n.filter)*(24*3600)  #add 30 days to filter out
        
        THIS=c("TagCode","Species","Sex","ReleaseDate","Latitude","Longitude","Latitude.prev",
               "Longitude.prev","Date.local","Date.local.prev","Block","Block.prev","same.blk",
               "zone.prev","zone")
        dat1$same.blk=with(dat1,ifelse(Block==Block.prev,"Y","N"))
        
        dat1$Dummy=with(dat1,paste(Block,Date.local))
        dat1=dat1[!duplicated(dat1$Dummy),]
        dat1=dat1[order(dat1$Date.local),]
        
        Stay.in.zone=subset(dat1,same.blk=="Y"| is.na(same.blk),select=THIS)
        Leave.zone=subset(dat1,same.blk=="N",select=THIS)
        
        #add intermediate days for those that stay in grid
        if(nrow(Stay.in.zone)>0)
        {
          Store.stay=vector('list',nrow(Stay.in.zone))
          for(s in 1:nrow(Stay.in.zone))
          {
            a=Stay.in.zone[s,]
            n=as.numeric(a$Date.local-a$Date.local.prev)
            if(n>1)
            {
              a<- a[rep(row.names(a), n), ]
              Day.seq=seq(1,n,1)*24*60*60
              if(!length(Day.seq)==nrow(a)) Day.seq=Day.seq[1:nrow(a)]
              a$Date.local=a$Date.local.prev+Day.seq  
              a$Block.prev=a$Block
              a$Date.local.prev=as.POSIXlt(c(NA,as.character(a$Date.local[1:(nrow(a)-1)])))
            }
            Store.stay[[s]]=a
            
          }
          Store.stay=do.call(rbind,Store.stay)
          Stay.in.zone=Store.stay
        }
        
        #interpolate position and day for those that left grid
        if(nrow(Leave.zone)>0)
        {
          Store=vector('list',nrow(Leave.zone))
          for(s in 1:nrow(Leave.zone))
          {
            a=Leave.zone[s,]
            n=as.numeric(a$Date.local-a$Date.local.prev)
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),c(Longitude,Latitude),n=n, addStartEnd=F))
            x=a
            if(nrow(b)>1)
            {
              #Add corners
              if(a$zone.prev=="North" & a$zone%in%c("Closed.metro","Zone1","West"))
              {
                n=as.numeric(a$Date.local-a$Date.local.prev)
                n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
                n2=distCosine(Exmouth,Shark.bay)/1000
                n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
                TN=n1+n2+n3
                n1=round(n*(n1/TN))
                n2=round(n*(n2/TN))
                n3=round(n*(n3/TN))
                
                b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
                b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
                b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
                b=rbind(b,b1,b2)
              }
              if(a$zone.prev=="North" & a$zone=="Zone2")
              {
                n=as.numeric(a$Date.local-a$Date.local.prev)
                n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
                n2=distCosine(Exmouth,Shark.bay)/1000
                n3=distCosine(Shark.bay,Cape.Leuwin)/1000
                n4=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
                TN=n1+n2+n3+n4
                n1=round(n*(n1/TN))
                n2=round(n*(n2/TN))
                n3=round(n*(n3/TN))
                n4=round(n*(n4/TN))
                b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
                b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
                b2=gcIntermediate(Shark.bay,Cape.Leuwin,n=n3, addStartEnd=F)
                b3=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n4, addStartEnd=F))
                b=rbind(b,b1,b2,b3)
              }
              if(a$zone.prev%in%c("Closed.ningaloo") & a$zone=="Zone2")
              {
                n=as.numeric(a$Date.local-a$Date.local.prev)
                n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
                n2=distCosine(Shark.bay,Cape.Leuwin)/1000    
                n3=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
                TN=n1+n2+n3
                n1=round(n*(n1/TN))
                n2=round(n*(n2/TN))
                n3=round(n*(n3/TN))
                b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
                b1=gcIntermediate(Shark.bay,Cape.Leuwin,n=n2, addStartEnd=F)
                b2=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n3, addStartEnd=F))
                b=rbind(b,b1,b2)
              }
              if(a$zone.prev%in%c("Closed.ningaloo") & a$zone%in%c("Zone1","Closed.metro"))
              {
                n=as.numeric(a$Date.local-a$Date.local.prev)
                n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
                n2=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
                TN=n1+n2
                n1=round(n*(n1/TN))
                n2=round(n*(n2/TN))
                b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
                b1=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n2, addStartEnd=F))
                b=rbind(b,b1)
              }
              if(a$zone.prev%in%c("Closed.metro") & a$zone=="Zone2")
              {
                n=as.numeric(a$Date.local-a$Date.local.prev)
                n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
                n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
                TN=n1+n2
                n1=round(n*(n1/TN))
                n2=round(n*(n2/TN))
                b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
                b1=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n2, addStartEnd=F))
                b=rbind(b,b1)
              }
              if(a$zone.prev=="Zone2" & a$zone%in%c("Closed.ningaloo"))
              {
                n=as.numeric(a$Date.local-a$Date.local.prev)
                n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
                n2=distCosine(Cape.Leuwin,Shark.bay)/1000
                n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
                TN=n1+n2+n3
                n1=round(n*(n1/TN))
                n2=round(n*(n2/TN))
                n3=round(n*(n3/TN))
                b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
                b1=gcIntermediate(Cape.Leuwin,Shark.bay,n=n2, addStartEnd=F)
                b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
                b=rbind(b,b1,b2)
              }
              if(a$zone.prev=="Zone2" & a$zone%in%c("Closed.metro"))
              {
                n=as.numeric(a$Date.local-a$Date.local.prev)
                n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
                n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
                TN=n1+n2
                n1=round(n*(n1/TN))
                n2=round(n*(n2/TN))
                b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
                b1=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n2, addStartEnd=F))
                b=rbind(b,b1)
              }
              if(a$zone.prev=="Zone2" & a$zone%in%c("North"))
              {
                n=as.numeric(a$Date.local-a$Date.local.prev)
                n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
                n2=distCosine(Cape.Leuwin,Shark.bay)/1000
                n3=distCosine(Shark.bay,Exmouth)/1000
                n4=distCosine(Exmouth,c(a$Longitude,a$Latitude))/1000
                TN=n1+n2+n3+n4
                n1=round(n*(n1/TN))
                n2=round(n*(n2/TN))
                n3=round(n*(n3/TN))
                n4=round(n*(n4/TN))
                b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
                b1=gcIntermediate(Cape.Leuwin,Shark.bay,n=n2, addStartEnd=F)
                b2=gcIntermediate(Shark.bay,Exmouth,n=n3, addStartEnd=F)
                b3=with(a,gcIntermediate(Exmouth,c(Longitude,Latitude),n=n4, addStartEnd=F))
                b=rbind(b,b1,b2,b3)
              }
              if(a$zone.prev=="Zone1" & a$zone%in%c("Closed.ningaloo"))
              {
                n=as.numeric(a$Date.local-a$Date.local.prev)
                n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
                n2=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
                TN=n1+n2
                n1=round(n*(n1/TN))
                n2=round(n*(n2/TN))
                b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
                b1=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n2, addStartEnd=F))
                b=rbind(b,b1)
              }
              x=a[rep(seq_len(nrow(a)), n),]
              this.na=match(c("Latitude.prev","Longitude.prev","same.blk"),names(x))
              x[,this.na]=NA
              if(!nrow(x)==nrow(b)) x=x[1:nrow(b),]
              x$Longitude=b[,1]
              x$Latitude=b[,2]
              Day.seq=seq(1,n,1)*24*60*60
              if(!length(Day.seq)==nrow(b)) Day.seq=Day.seq[1:nrow(b)]
              x$Date.local=a$Date.local.prev+Day.seq       
              x$Latitude.prev= c(NA,x$Latitude[1:(nrow(x)-1)])
              x$Longitude.prev=c(NA,x$Longitude[1:(nrow(x)-1)])
              x$Date.local.prev=as.POSIXlt(c(NA,as.character(x$Date.local[1:(nrow(x)-1)])))
              x$Block=with(x,paste(floor(-Latitude),floor(Longitude),sep=""))
              x$Block.prev=with(x,paste(floor(-Latitude.prev),floor(Longitude.prev),sep=""))
              
            }
            
            Store[[s]]=x
          }  
          Store=do.call(rbind,Store)
        }
        All=NULL
        if(nrow(Leave.zone)==0) All=Stay.in.zone
        if(nrow(Leave.zone)>0) All=rbind(Store,Stay.in.zone)
        
        All=All[order(All$Date.local),match(c("TagCode","Species","Sex","Date.local",
                                              "Date.local.prev","Block","Block.prev","zone"),names(All))]
        All$dummy= c(NA,All$Block[1:(nrow(All)-1)])
        All$Block.prev=with(All,ifelse(Block.prev=="NANA",dummy,Block.prev))  
        All=subset(All,!is.na(Date.local.prev))
        
        return(All[,-match('dummy',names(All))])
      }
      system.time(for (i in 1:length(shks))Time.mon.grid[[i]]=Tim.mon.grid(shks[i]))
      Time.mon.grid=do.call(rbind,Time.mon.grid)  
      
      Time.mon.grid$Year=as.POSIXlt(as.character(Time.mon.grid$Date.local))$year+1900
      Time.mon.grid$Month=as.POSIXlt(as.character(Time.mon.grid$Date.local))$mon+1 
      Time.mon.grid$Day=as.POSIXlt(as.character(Time.mon.grid$Date.local))$mday
      
      
      #Aggregate number of days monitored by zone
      Time.mon.grid$N=1
      
      #monthly analysis
      fn.prob=function(dat)
      {
        agg.month.year=aggregate(N~Year+Month+TagCode,dat,sum)
        names(agg.month.year)[match("N",names(agg.month.year))]="Total"  
        agg.zone.month.year=aggregate(N~Year+Month+zone+TagCode,dat,sum)
        agg.zone.month.year=merge(agg.zone.month.year,agg.month.year,
                                  by=c("Year","Month","TagCode"),all.x=T)
        agg.zone.month.year$N.out=agg.zone.month.year$Total-agg.zone.month.year$N
        
        #by year
        agg.zone.year=aggregate(cbind(N,Total,N.out)~Year+zone+TagCode,agg.zone.month.year,sum)
        
        
        
        Zns=unique(agg.zone.month.year$zone)
        Store.dat=vector('list',length(Zns))
        names(Store.dat)=Zns
        Store.mod=Store.dat
        for(x in 1:length(Zns))
        {
          d1=subset(agg.zone.month.year,zone==Zns[x])
          MOD=glm(cbind(N, N.out) ~ Year*Month,data =d1,family=binomial)
          #MOD=glmer(cbind(N, N.out) ~ Year+Month+(1 | TagCode),data=d1, family = binomial)
          d1$prop=d1$N/d1$Total
          d1$pred=predict(MOD,type='response')
          Store.dat[[x]]=d1
          Store.mod[[x]]=MOD
        }
        Preds=do.call(rbind,Store.dat)
        
        Preds.rshp=reshape(subset(Preds,select=c(pred,Year,Month,TagCode,zone)),
                           v.names = "pred", idvar = c("Year","Month","TagCode"),timevar = "zone", direction = "wide")
      }
      fn.prob(subset(Time.mon.grid,Species=="Gummy")) 
      
    }
    
  }
}





#Connectivity plot of overall movement
  #Image option
Image.mig.fn=function(what,SP)
{
  a=subset(Species.time.zone,Species==what)
  
  id=match(Zns,names(a))
  a[,id]=a[,id]/rowSums(a[,id])
  crap=as.data.frame(as.matrix(COL.prop))
  names(crap)="Col"
  crap$Zone.rel=rownames(crap)
  a=merge(a,crap,by="Zone.rel")
  
  a$Sort=with(a,ifelse(Zone.rel=="North",6,
                       ifelse(Zone.rel=="Closed.ningaloo",5,
                              ifelse(Zone.rel=="West",4,
                                     ifelse(Zone.rel=="Closed.metro",3,
                                            ifelse(Zone.rel=="Zone1",2,1))))))
  
  a$Zone.rel=as.character(a$Zone.rel)
  
  a$Zone.rel=with(a,ifelse(Zone.rel=="Closed.ningaloo","5.Ningaloo",
                           ifelse(Zone.rel=="Closed.metro","3.Metro",
                                  ifelse(Zone.rel=="Zone1","2.Zn1",
                                         ifelse(Zone.rel=="Zone2","1.Zn2",
                                                ifelse(Zone.rel=="North","6.North","4.West"))))))
  a$Zone.rel=as.factor(a$Zone.rel)
  MATRX=aggregate(cbind(North,Closed.ningaloo,West,Closed.metro,Zone1,Zone2)~Zone.rel,a,mean)

  colnames(MATRX)[match(c("North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2"),
                        colnames(MATRX))]=c("6.North","5.Ningaloo","4.West","3.Metro","2.Zn1","1.Zn2")
  MATRX=MATRX[,match(c("Zone.rel","1.Zn2","2.Zn1","3.Metro","4.West","5.Ningaloo","6.North"),names(MATRX))]
  MATRX$Zone.rel=as.character(MATRX$Zone.rel)
  
  ALL.zns=c("1.Zn2","2.Zn1","3.Metro","4.West","5.Ningaloo","6.North")
  ID=ALL.zns[which(!ALL.zns%in%MATRX$Zone.rel)]
  if(length(ID)>0)
  {
    ss=as.data.frame(matrix(nrow=length(ID),ncol=ncol(MATRX)))
    names(ss)=names(MATRX)
    ss[,]=0
    ss$Zone.rel=ID
    MATRX=rbind(MATRX,ss)
    MATRX=MATRX[order(MATRX$Zone.rel),]
  }
  MTRX=t(MATRX[,-1])  
  image(1:nrow(MTRX),1:ncol(MTRX),as.matrix(MTRX),col =couleurs,breaks=BREAKS,xaxt='n',yaxt='n',ylab="",xlab="")
  
  ZN=rev(c("North","Ning.","West","Metro","Zn.1","Zn.2"))
  ZN.x=substr(as.character(MATRX$Zone.rel),3,20)
  ZN.x=ifelse(ZN.x=="Zn1","Zn.1",ifelse(ZN.x=="Zn2","Zn.2",ifelse(ZN.x=="Ningaloo","Ning.",ZN.x)))
  
  axis(2,1:ncol(MTRX),F,las=1)
  axis(1,1:nrow(MTRX),F)
  
  if(what%in%c("Dusky","Thickskin")) axis(2,1:ncol(MTRX),ZN.x,las=1,cex.axis=1)
  if(what%in%c("Thickskin","Whiskery")) axis(1,1:nrow(MTRX),ZN,cex.axis=.9)
  box()
  mtext(SP,3)
  lines(.5:6.5,.5:6.5,lwd=1.5)
  return(MTRX)
}

N.int=1000
colfunc <- colorRampPalette(c("navy", "cadetblue","white"))
couleurs=rev(colfunc(N.int))
#couleurs=rev(gray(seq(0,0.9,length=N.int)))
BREAKS=seq(0,1,length.out=N.int+1)

tiff(file="Connectivity/Image.All.together.tiff",width = 2000, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mar = c(1, 1, 1, 1),oma=c(3,4,0,0),xpd=T,mgp=c(1,.7,0))
Image.mig.fn("Dusky","Dusky shark")
Image.mig.fn("Thickskin","Sandbar shark")
Image.mig.fn("Gummy","Gummy shark")
Image.mig.fn("Whiskery","Whiskery shark")
#SQ=BREAKS
SQ=seq(BREAKS[1],BREAKS[length(BREAKS)],.1)
color.legend(6.5,0.5,7,6.5,SQ,rect.col=couleurs,gradient="y",col=1,cex=1)
mtext("Zone released",2,outer=T,line=2.5,cex=1.5)
mtext("Zone detected",1,outer=T,line=1.5,cex=1.5)
dev.off()



#Linear option
linear.mig.fn=function(what,Inc)
{
  a=subset(Species.time.zone,Species==what)
  
  id=match(Zns,names(a))
  a[,id]=a[,id]/rowSums(a[,id])
  crap=as.data.frame(as.matrix(COL.prop))
  names(crap)="Col"
  crap$Zone.rel=rownames(crap)
  a=merge(a,crap,by="Zone.rel")
  
  a$Sort=with(a,ifelse(Zone.rel=="North",1,
                       ifelse(Zone.rel=="Closed.ningaloo",2,
                              ifelse(Zone.rel=="West",3,
                                     ifelse(Zone.rel=="Closed.metro",4,
                                            ifelse(Zone.rel=="Zone1",5,6))))))
  ZN=c("North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2")
  n.zn=length(ZN)
  
  add=ZN[which(!ZN%in%unique(a$Zone.rel))]
  if(length(add)>0)     #add dummy to plot zones with no releases
  {
    add1=a[1:length(add),]
    add1$Zone.rel=add
    add1[match(ZN,names(add1))]=0
    a=rbind(a,add1,add1)
  }
  a=a[order(a$Sort),]
  a$Zone.rel=as.character(a$Zone.rel)
  
  a$Zone.rel=with(a,ifelse(Zone.rel=="Closed.ningaloo","2.Ningaloo",
                           ifelse(Zone.rel=="Closed.metro","4.Metro",
                                  ifelse(Zone.rel=="Zone1","5.Zn1",
                                         ifelse(Zone.rel=="Zone2","6.Zn2",
                                                ifelse(Zone.rel=="North","1.North","3.West"))))))
  factors=a$Zone.rel
  Tab=table(factors)
  Prop=Tab/sum(Tab)
  CLs=COL.prop
  names(CLs)=names(Tab)
  BG.col=COL.prop[match(names(Tab),names(CLs))]
  
  plot(1:n.zn,1:n.zn,xaxt='n',yaxt='n',col="transparent",ylab="",xlab="",axes=F,ylim=c(0,n.zn))
  NN.zn=c(n.zn,n.zn-cumsum(Prop*n.zn))
  ZN.lab=c("North","Ningaloo","West","Metro","Zone1","Zone2")
  x1=n.zn*.4
  x2=n.zn*.6
  YYs=matrix(ncol=2,nrow=n.zn)
  for(i in 1:n.zn) YYs[i,]=c(NN.zn[i+1],NN.zn[i])
  
  YYs=as.data.frame(YYs)
  colnames(YYs)=c("FROM","TO")
  YYs$ZonE=c("1.North","2.Ningaloo","3.West","4.Metro","5.Zn1","6.Zn2")
  
  #Add connections
  a$Zone.rel=as.factor(a$Zone.rel)
  MATRX=aggregate(cbind(North,Closed.ningaloo,West,Closed.metro,Zone1,Zone2)~Zone.rel,a,sum)
  MATRX[,2:ncol(MATRX)]=MATRX[,2:ncol(MATRX)]/rowSums(MATRX[,2:ncol(MATRX)])
  
  colnames(MATRX)[match(c("North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2"),
                        colnames(MATRX))]=c("1.North","2.Ningaloo","3.West","4.Metro","5.Zn1","6.Zn2")
  MATRX$Zone.rel=as.character(MATRX$Zone.rel)
  
  PAIRS=expand.grid(MATRX$Zone.rel,MATRX$Zone.rel)
  names(PAIRS)=c("Start","End")
  CLs=BG.col
  names(CLs)[match(c("North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2"),
                   names(CLs))]=c("1.North","2.Ningaloo","3.West","4.Metro","5.Zn1","6.Zn2")
  CLs=data.frame(col=CLs,Zone.rel=names(CLs))
  PAIRS=merge(PAIRS,CLs,by.x="Start",by.y="Zone.rel",all.x=T)
  
  LWD=reshape(MATRX,idvar="Zone.rel",varying=list(2:7),direction="long")
  LWD=LWD[order(LWD$Zone.rel),]
  n=length(2:ncol(MATRX))
  LWD$time=rep(names(MATRX)[2:ncol(MATRX)],n)
  names(LWD)=c("Start","End","LWD")
  PAIRS=merge(PAIRS,LWD,by=c("Start","End"))
  PAIRS=subset(PAIRS,LWD>0)
  
  PAIRS=merge(PAIRS,YYs,by.x="Start",by.y="ZonE",all.x=T)
  names(PAIRS)[match(c("FROM","TO"),names(PAIRS))]=c("From.st","To.st")
  
  PAIRS=merge(PAIRS,YYs,by.x="End",by.y="ZonE",all.x=T)
  names(PAIRS)[match(c("FROM","TO"),names(PAIRS))]=c("From.end","To.end")
  
  
  
  for(s in 1:nrow(PAIRS))
  {
    Lcol=Arcol=as.character(PAIRS$col[s])
    Lwd= PAIRS$LWD[s]*Inc
    Ar.l=Lwd/10
    Ar.w=Ar.l/2
    if(!PAIRS$Start[s]==PAIRS$End[s])
    {
      Y1=sample(seq(PAIRS$From.st[s],PAIRS$To.st[s],.01),1)
      Y2=sample(seq(PAIRS$From.end[s],PAIRS$To.end[s],.01),1)
      Cur=Y1/Y2*Ar.w
      if(Cur==0) Cur=.5
      curvedarrow(from=c(x1*0.975,Y1),to=c(x1*0.975,Y2),arr.pos=1,arr.width=Ar.w,arr.length= Ar.l,lcol=Lcol,arr.col=Arcol,
                  lwd=Lwd,arr.type='triangle',curve = Cur,col=1)
    }
    
    if(PAIRS$Start[s]==PAIRS$End[s])
    {
      Y1=PAIRS$From.st[s]
      Y2=PAIRS$To.end[s]
      Cur=Y1/Y2
      if(Cur==0) Cur=.5
      curvedarrow(from=c(x2*1.025,Y1),to=c(x2*1.025,Y2),arr.pos=1,arr.width=Ar.w,arr.length= Ar.l,lcol=Lcol,arr.col=Arcol,
                  lwd=Lwd,arr.type='triangle',curve = Cur,col=1)
    }
    
  }
  
  
  #Add zones
  for(i in 1:n.zn) 
  {
    polygon(x=c(x1,x2,x2,x1),y=c(NN.zn[i+1],NN.zn[i+1],NN.zn[i],NN.zn[i]),col=BG.col[i])
    text(n.zn/1.225,mean(c(NN.zn[i],NN.zn[i+1])),ZN.lab[i],cex=1.5,pos=4,col=BG.col[i])
  }
  
  
}



  #Circular option
#note: with of arrow is proportion to movement proportion
#      with or each bar is proportion to number of sharks detected
cicle.mig.fn=function(what,Inc,CEX)
{
  a=subset(Species.time.zone,Species==what)
  a=subset(Species.time.zone,Species==what)
  id=match(Zns,names(a))
  a[,id]=a[,id]/rowSums(a[,id])
  crap=as.data.frame(as.matrix(COL.prop))
  names(crap)="Col"
  crap$Zone.rel=rownames(crap)
  a=merge(a,crap,by="Zone.rel")
  
   a$Sort=with(a,ifelse(Zone.rel=="North",1,
         ifelse(Zone.rel=="Closed.ningaloo",2,
        ifelse(Zone.rel=="West",3,
       ifelse(Zone.rel=="Closed.metro",4,
      ifelse(Zone.rel=="Zone1",5,6))))))
  ZN=c("North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2")
  
  add=ZN[which(!ZN%in%unique(a$Zone.rel))]
  if(length(add)>0)     #add dummy to plot zones with no releases
  {
    add1=a[1:length(add),]
    add1$Zone.rel=add
    add1[match(ZN,names(add1))]=0
    a=rbind(a,add1,add1,add1)
  }
  a=a[order(a$Sort),]
  a$Zone.rel=as.character(a$Zone.rel)
    
  a$Zone.rel=with(a,ifelse(Zone.rel=="Closed.ningaloo","2.Ningaloo",
        ifelse(Zone.rel=="Closed.metro","4.Metro",
        ifelse(Zone.rel=="Zone1","5.Zn1",
        ifelse(Zone.rel=="Zone2","6.Zn2",
        ifelse(Zone.rel=="North","1.North","3.West"))))))
  
  
  factors=a$Zone.rel
  Tab=table(factors)
  Prop=Tab/sum(Tab)
  CLs=COL.prop
  names(CLs)=names(Tab)
  BG.col=COL.prop[match(names(Tab),names(CLs))]
  circos.initialize(factors, xlim = c(0,1),sector.width = Prop)
   circos.trackPlotRegion(factors, ylim = c(0, 1), track.height = 0.15,
     bg.border = 1, bg.col = BG.col, 
     panel.fun = function(x, y)circos.text(.5, 1.5, get.cell.meta.data("sector.index"),cex=CEX))
  
  #Add connections
  a$Zone.rel=as.factor(a$Zone.rel)
  MATRX=aggregate(cbind(North,Closed.ningaloo,West,Closed.metro,Zone1,Zone2)~Zone.rel,a,sum)
  MATRX[,2:ncol(MATRX)]=MATRX[,2:ncol(MATRX)]/rowSums(MATRX[,2:ncol(MATRX)])
  
  colnames(MATRX)[match(c("North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2"),
                        colnames(MATRX))]=c("1.North","2.Ningaloo","3.West","4.Metro","5.Zn1","6.Zn2")
  MATRX$Zone.rel=as.character(MATRX$Zone.rel)
  
  PAIRS=expand.grid(MATRX$Zone.rel,MATRX$Zone.rel)
  names(PAIRS)=c("Start","End")
  CLs=BG.col
  names(CLs)[match(c("North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2"),
                   names(CLs))]=c("1.North","2.Ningaloo","3.West","4.Metro","5.Zn1","6.Zn2")
  CLs=data.frame(col=CLs,Zone.rel=names(CLs))
  PAIRS=merge(PAIRS,CLs,by.x="Start",by.y="Zone.rel",all.x=T)
  
  LWD=reshape(MATRX,idvar="Zone.rel",varying=list(2:7),direction="long")
  LWD=LWD[order(LWD$Zone.rel),]
  n=length(2:ncol(MATRX))
  LWD$time=rep(names(MATRX)[2:ncol(MATRX)],n)
  names(LWD)=c("Start","End","LWD")
  PAIRS=merge(PAIRS,LWD,by=c("Start","End"))
  PAIRS=subset(PAIRS,LWD>0)
  Seq=seq(0.1,.9,.1)
  Seq=Seq[!Seq%in%c(.1,.5)]
  PAIRS$Pos.end=with(PAIRS,ifelse(Start==End,.5,sample(Seq,n-1)))
  PAIRS$Pos.strt=with(PAIRS,ifelse(Start==End,.1,sample(Seq,n-1)))
  
  for(s in 1:nrow(PAIRS))
  {
    circos.link(as.character(PAIRS$Start[s]),PAIRS$Pos.strt[s],
                as.character(PAIRS$End[s]),PAIRS$Pos.end[s],
                lwd=PAIRS$LWD[s]*Inc,col = as.character(PAIRS$col[s]),
                rou1 = 0.6, rou2 = 0.6)
   #directional = -1,arr.length = PAIRS$LWD[s]*2,arr.width = PAIRS$LWD[s]*2/2, 
  }

}

LEG=function(INSET,CX) 
{
  legend("bottomleft", paste(c(10,50),"%",sep=""),lwd=c(.1,.5)*25,title="Per. movement",bty='n',cex=CX,
         horiz=F,inset=INSET)
}


tiff(file="Connectivity/Connectivity.plot.All.together.tiff",width = 2000, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mar = c(1, 1, 1, 1),lwd=1.5)
cicle.mig.fn("Dusky",Inc=25,1.15)
mtext("Dusky",3,line=-0.5,cex=1.5)

cicle.mig.fn("Thickskin",Inc=25,0.8)
mtext("Sandbar",3,line=-0.5,cex=1.5)
LEG(c(-0.01,-.03),0.95)

cicle.mig.fn("Gummy",Inc=25,0.7)
mtext("Gummy",3,line=-0.5,cex=1.5)

cicle.mig.fn("Whiskery",Inc=25,1.25)
mtext("Whiskery",3,line=-0.5,cex=1.5)

dev.off()

tiff(file="Connectivity/Connectivity.plot.Dusky.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mar = c(1, 1, 1, 1),lwd=1.5)
cicle.mig.fn("Dusky",Inc=25,1.75)
LEG(c(0,0),1.5)
dev.off()

tiff(file="Connectivity/Connectivity.plot.Sandbar.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mar = c(1, 1, 1, 1),lwd=1.5)
cicle.mig.fn("Thickskin",Inc=25,1.75)
LEG(c(0,0),1.5)
dev.off()

tiff(file="Connectivity/Connectivity.plot.Gummy.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
cicle.mig.fn("Gummy",Inc=25,1.5)
LEG(c(0,0),1.5)
dev.off()

tiff(file="Connectivity/Connectivity.plot.Whiskery.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
cicle.mig.fn("Whiskery",Inc=25,1.75)
LEG(c(0,0),1.35)
dev.off()



#7. -- Data set for pop dyn modelling ---

Tagging.pop.dyn=Detections[,-match(c("SerialNumber","Time.local",
              "Station","Area","Area.release","Array","SerialNumber.prev"),names(Detections))]
#Tagging.pop.dyn=subset(Detections, Dist.moved.conseq.det>0 | Dist.moved.rel.det>0 )
names(Tagging.pop.dyn)[match(c("Latitude","Longitude"),names(Tagging.pop.dyn))]=c("Lat.rec","Long.rec")
Tagging.pop.dyn$FL=with(Tagging.pop.dyn,ifelse(FL<10,FL*100,FL))
AVrg.FL=aggregate(FL~Species,Tagging.pop.dyn,mean)
names(AVrg.FL)[2]="Mean.FL"
Tagging.pop.dyn=merge(Tagging.pop.dyn,AVrg.FL,by="Species")
Tagging.pop.dyn$FL=with(Tagging.pop.dyn,ifelse(is.na(FL),Mean.FL,FL))

#add release and recapture BLOCK 
Tagging.pop.dyn$Block.rel=with(Tagging.pop.dyn,-(ceiling(ReleaseLatitude)) * 100 +(floor(ReleaseLongitude)-100))
Tagging.pop.dyn$Block.rec=with(Tagging.pop.dyn,-(ceiling(Lat.rec)) * 100 +(floor(Long.rec)-100))

#add zone released
Tagging.pop.dyn$Rel.zone=as.character(with(Tagging.pop.dyn,
        ifelse(ReleaseLongitude>=116.5 & ReleaseLatitude<=(-26),"Zone2",
        ifelse(ReleaseLongitude<116.5 & ReleaseLatitude<=(-33),"Zone1",
        ifelse(ReleaseLatitude>(-33) & ReleaseLatitude<=(-26) & ReleaseLongitude<116.5,"West",
        ifelse(ReleaseLatitude>(-26) & ReleaseLongitude<114,"Closed",
        ifelse(ReleaseLatitude>(-26) & ReleaseLongitude>=114 & ReleaseLongitude<123.75,"North",
        ifelse(ReleaseLatitude>(-26) & ReleaseLongitude>=123.75,"Joint",NA))))))))
Tagging.pop.dyn$Rel.zone=with(Tagging.pop.dyn,ifelse(ReleaseLongitude>=129 & ReleaseLatitude<=(-26),"SA",Rel.zone))

#add zone rec
Tagging.pop.dyn$Rec.zone=as.character(with(Tagging.pop.dyn,
        ifelse(Long.rec>=116.5 & Lat.rec<=(-26),"Zone2",
        ifelse(Long.rec<116.5 & Lat.rec<=(-33),"Zone1",
        ifelse(Lat.rec>(-33) & Lat.rec<=(-26) & Long.rec<116.5,"West",
        ifelse(Lat.rec>(-26) & Long.rec<114,"Closed",
        ifelse(Lat.rec>(-26) & Long.rec>=114 & Long.rec<123.75,"North",
        ifelse(Lat.rec>(-26) & Long.rec>=123.75,"Joint",NA))))))))

Tagging.pop.dyn$Rec.zone=with(Tagging.pop.dyn,ifelse(Long.rec>=129 & Lat.rec<=(-26),"SA",Rec.zone))

Tagging.pop.dyn$Week=with(Tagging.pop.dyn,ifelse(Day%in%1:8,1,ifelse(Day%in%9:16,2,ifelse(Day%in%17:24,3,4))))
Tagging.pop.dyn$Week.rel=with(Tagging.pop.dyn,ifelse(Day.rel%in%1:8,1,ifelse(Day.rel%in%9:16,2,ifelse(Day.rel%in%17:24,3,4))))

# Minimum time step for modelling movement rates
#check if weekly time step is ok (output is paste of TagCode, Yr, Month, Week). Any value >1 means two different zone
fn.what.time.step=function(DAT,TIMSTP)
{
  a=sort(unique(DAT$TagCode))
  b=table(DAT$TagCode,DAT$Rec.zone)
  b[b>0]=1
  this=rowSums(b)
  this=this[this>1]
  id=names(this)
  d=subset(DAT,TagCode%in%id)
  
  if(TIMSTP=='WEEKLY')
  {
    d$dummy=with(d,paste(TagCode,Year,Month,Week,Rec.zone))
    d=d[!duplicated(d$dummy),]
    d=d[order(d$dummy),]
    d$dummy2=with(d,paste(TagCode,Year,Month,Week))
    x=table(d$dummy2)
  }
  
  if(TIMSTP=='MONTHLY')
  {
    d$dummy=with(d,paste(TagCode,Year,Month,Rec.zone))
    d=d[!duplicated(d$dummy),]
    d=d[order(d$dummy),]
    d$dummy2=with(d,paste(TagCode,Year,Month))
    x=table(d$dummy2)
  }
  
  return(x[x>1])
}
Pop.din.sp=c("Dusky","Thickskin","Gummy","Whiskery")
Time.step=vector('list',length(Pop.din.sp))
names(Time.step)=Pop.din.sp
#for(i in 1:length(Time.step))Time.step[[i]]=fn.what.time.step(DAT=subset(Tagging.pop.dyn,Species==Pop.din.sp[i]),'WEEKLY')
for(i in 1:length(Time.step))Time.step[[i]]=fn.what.time.step(DAT=subset(Tagging.pop.dyn,Species==Pop.din.sp[i]),'MONTHLY')

#model time step= week (several indivdiuals of most species detected in more than 1 zone within a month)
fn.group=function(DAT,K.f,Linf.f,to.f,K.m,Linf.m,to.m,MX.AGE,do.age,Resp.var)
{
  DAT$Sex=with(DAT,ifelse(Sex=="U","F",Sex))
  
  if(do.age=="YES")
  {
    #assign age
    if(unique(DAT$Species)=="Gummy")DAT$FL=DAT$FL*1.0837+4.6424   #gummy growth pars are in TL
    
    
    #Calculate age from length
    DAT$Age=with(DAT,ifelse(Sex=="F" & FL<Linf.f,to.f-(1/K.f)*log(1-(FL/Linf.f)),
                            ifelse(Sex=="M" & FL<Linf.m,to.m-(1/K.m)*log(1-(FL/Linf.m)),NA))) 
    #fix age for length > Linf 
    fem=subset(DAT,Sex=="F");mal=subset(DAT,Sex=="M")
    max.A.f=max(subset(fem$Age,!fem$Age=='Inf'),na.rm=T);max.A.m=max(subset(mal$Age,!mal$Age=='Inf'),na.rm=T)
    
    
    DAT$Age=with(DAT,ifelse(Sex=="F" & FL>=Linf.f,runif(1,max.A.f,MX.AGE),
                            ifelse(Sex=="M" & FL>=Linf.m,runif(1,max.A.m,MX.AGE),Age)))
    DAT$Age=with(DAT,ifelse(Age<0,0,Age))
    DAT$Age=round(DAT$Age)
  }
  
  DAT$Number=1
  
  #get single hit per day per spatial strata
  DAT$Dummy.blk=with(DAT,paste(TagCode,Year,Month,Day,Block.rec))   
  Dat.blk=DAT[!duplicated(DAT$Dummy.blk),]
  DAT$Dummy.zn=with(DAT,paste(TagCode,Year,Month,Day,Rec.zone))
  Dat.zn=DAT[!duplicated(DAT$Dummy.zn),]
  
  #get proportion of time (days over total days detected within the time step)
  Ag.blk=aggregate(Number~Block.rec+Week+Month+Year+TagCode,Dat.blk,sum)
  Ag.blk$dummy=with(Ag.blk,paste(Week,Month,Year,TagCode))
  
  Ag.zn=aggregate(Number~Rec.zone+Week+Month+Year+TagCode,Dat.zn,sum)
  Ag.zn$dummy=with(Ag.zn,paste(Week,Month,Year,TagCode))
  
  dummy.blk=aggregate(Number~dummy,Ag.blk,sum)
  dummy.zn=aggregate(Number~dummy,Ag.zn,sum)
  names(dummy.blk)[2]=names(dummy.zn)[2]="Total"
  
  Ag.blk=merge(Ag.blk,dummy.blk,by="dummy",all.x=T)
  Ag.blk$Prop=with(Ag.blk,Number/Total)
  
  Ag.zn=merge(Ag.zn,dummy.zn,by="dummy",all.x=T)
  Ag.zn$Prop=with(Ag.zn,Number/Total)
  
  if(Resp.var=="Number")
  {
    Ag.blk.sub=subset(Ag.blk,Prop<1)
    Ag.blk=subset(Ag.blk,Prop==1)    
    if(nrow(Ag.blk.sub)>0)
    {
      esto=unique(Ag.blk.sub$dummy)
      store=Ag.blk.sub[1:length(esto),]
      store[,]=NA
      for(x in 1:length(esto))
      {
        a=subset(Ag.blk.sub,dummy==esto[x])      
        id=which(a$Prop==max(a$Prop))
        if(length(id>1)) id=sample(id,1)
        store[x,]=a[id,]    
      }    
      Ag.blk=rbind(Ag.blk,store)
    }    
    Ag.blk$Prop=1
    
    Ag.zn.sub=subset(Ag.zn,Prop<1)
    Ag.zn=subset(Ag.zn,Prop==1)    
    if(nrow(Ag.zn.sub)>0)
    {
      esto=unique(Ag.zn.sub$dummy)
      store=Ag.zn.sub[1:length(esto),]
      store[,]=NA
      for(x in 1:length(esto))
      {
        a=subset(Ag.zn.sub,dummy==esto[x])      
        id=which(a$Prop==max(a$Prop))
        if(length(id>1)) id=sample(id,1)
        store[x,]=a[id,]    
      }
      Ag.zn=rbind(Ag.zn,store)
    }    
    Ag.zn$Prop=1    
  }
  
  if(Resp.var=="Proportion")
  {
    Ag.blk=Ag.blk
    Ag.zn=Ag.zn
  }
  
  BLK.rel=subset(DAT,select=c(TagCode,Sex,FL,Block.rel,Year.rel,Month.rel,Week.rel))
  if(do.age=="YES") BLK.rel=subset(DAT,select=c(TagCode,Sex,Age,FL,Block.rel,Year.rel,Month.rel,Week.rel))
  BLK.rec=subset(Ag.blk,select=c(TagCode,Block.rec,Year,Month,Week,Prop))    
  
  Zn.rel=subset(DAT,select=c(TagCode,Sex,FL,Rel.zone,Year.rel,Month.rel,Week.rel))
  if(do.age=="YES") Zn.rel=subset(DAT,select=c(TagCode,Sex,Age,FL,Rel.zone,Year.rel,Month.rel,Week.rel))
  Zn.rec=subset(Ag.zn,select=c(TagCode,Rec.zone,Year,Month,Week,Prop)) 
  
  return(list(BLK.rel=BLK.rel,Zn.rel=Zn.rel,BLK.rec=BLK.rec,Zn.rec=Zn.rec,DAT.all=DAT))
  
}
Store.group=vector('list',length(Pop.din.sp))
names(Store.group)=Pop.din.sp
Store.group.prop=Store.group
for(i in 1:length(Store.group))
{
  Store.group[[i]]=fn.group(subset(Tagging.pop.dyn,Species==Pop.din.sp[i]),
    K.f=Gr[[i]][1],Linf.f=Gr[[i]][2],to.f=Gr[[i]][3],K.m=Gr[[i]][4],Linf.m=Gr[[i]][5],
    to.m=Gr[[i]][6],MX.AGE=Mx.age[[i]],do.age="NO",Resp.var="Number")
  
  Store.group.prop[[i]]=fn.group(subset(Tagging.pop.dyn,Species==Pop.din.sp[i]),
        K.f=Gr[[i]][1],Linf.f=Gr[[i]][2],to.f=Gr[[i]][3],K.m=Gr[[i]][4],Linf.m=Gr[[i]][5],
        to.m=Gr[[i]][6],MX.AGE=Mx.age[[i]],do.age="NO",Resp.var="Proportion")
}

#data for individual-based model  
Tagging.pop.dyn$Index=1:nrow(Tagging.pop.dyn)
fn.dat.for.ind.base.mdl=function(dat,n.min.days)    
{
  Tgs=unique(dat$TagCode)
  dummy=vector('list',length(Tgs))
  for(t in 1:length(Tgs))
  {
    x=subset(dat,TagCode==Tgs[t],select=c(TagCode,Index,
        ReleaseLatitude,ReleaseLongitude,Rel.zone,Year.rel,Month.rel,Day.rel,ReleaseDate,
        Lat.rec,Long.rec,Rec.zone,Year,Month,Day,Date.local))
    x$ReleaseDate=min(x$ReleaseDate)
    x=x[order(x$ReleaseDate),]
    x$Dup=with(x,paste(TagCode,Date.local))
    x=x[!duplicated(x$Dup),]
    x$dummy=NA
    x$Delta.t=as.numeric(with(x,difftime(Date.local,ReleaseDate,units="secs"))/(24*3600))
    x=x[order(x$Date.local),]
    x$Jul.day.rel=0
    x$Jul.day.rec=x$Delta.t
    x=subset(x,Delta.t>=n.min.days)
    N=nrow(x)
    if(N>1)
    {
 #     x$Jul.day.rec=ifelse(x$Jul.day.rec>=86400,x$Jul.day.rec/(24*3600),x$Jul.day.rec)
      x$Rel.zone=c(x$Rel.zone[1],x$Rec.zone[1:(N-1)])
 #     x$Jul.day.rel=c(x$Jul.day.rel[1],x$Jul.day.rec[1:(N-1)])
      x$dummy=c(x$Jul.day.rel[1],x$Jul.day.rec[1:(N-1)])
      x$Jul.day.rec=x$Jul.day.rec-x$dummy
    }
     dummy[[t]]=x
  }
  x=do.call(rbind,dummy)
  x=x[order(x$TagCode,x$Date.local),]
  x=x[,match(c("TagCode","Rel.zone","Rec.zone","Jul.day.rec","Index"),names(x))]
  names(x)=c("TagID","Rel.zn","Rec.zn","DaysAtLarge","Index")
  return(x)
}
Store.group_ind.base=vector('list',length(Pop.din.sp))
names(Store.group_ind.base)=Pop.din.sp
for(q in 1:length(Store.group))Store.group_ind.base[[q]]=fn.dat.for.ind.base.mdl(subset(Tagging.pop.dyn,Species==Pop.din.sp[q]),
                            n.min.days=30)


#export
setwd(handl_OneDrive("Analyses/Data_outs"))
for(i in 1:length(Store.group))
{
  #SS3 approach
  a=Store.group[[i]][1:4]
  NmS=ifelse(names(Store.group)[i]=="Dusky",'Dusky shark',
      ifelse(names(Store.group)[i]=='Whiskery','Whiskery shark',
      ifelse(names(Store.group)[i]=='Gummy','Gummy shark',
      ifelse(names(Store.group)[i]=='Thickskin','Sandbar shark',NA))))
  for(p in 1:length(a)) write.csv(a[[p]],paste(getwd(),'/',NmS,'/',NmS,"_Acous.Tag_",names(a)[p],"_","Acous.Tag.csv",sep=""),row.names=F)

  b=Store.group.prop[[i]][1:4]
  for(p in 1:length(b)) write.csv(b[[p]],paste(getwd(),'/',NmS,'/',NmS,"_Acous.Tag_",names(b)[p],"_","Acous.Tag.prop.csv",sep=""),row.names=F)
  
  #Individual-based model approach
  a=Store.group_ind.base[[i]]
  a=a[,-match("Index",names(a))]
  write.csv(a,paste(getwd(),'/',NmS,'/',NmS,'_Acous.Tag_Acous.Tag.Ind_based.csv',sep=""),row.names=F)
  
  #export data for mapping 
  if(names(Store.group)[i]%in%c("Gummy" ,"Whiskery"))
  {
    HnDDl=handl_OneDrive("Analyses/Movement rate estimation/Joint.estim_ind.base.mod/Show Gummy and whiskery outputs/")
    a=subset(Tagging.pop.dyn,Index%in%Store.group_ind.base[[i]]$Index,select=c(TagCode,
            Lat.rec, Long.rec,ReleaseLatitude,ReleaseLongitude,Latitude.prev,Longitude.prev))
    a$Latitude.prev=with(a,ifelse(is.na(Latitude.prev),ReleaseLatitude,Latitude.prev))
    a$Longitude.prev=with(a,ifelse(is.na(Longitude.prev),ReleaseLongitude,Longitude.prev))
    a=a[,-match(c("ReleaseLatitude","ReleaseLongitude"),names(a))]
    write.csv(a,paste(HnDDl,names(Store.group.prop)[i],"_","Raw.Acous.Tag.csv",sep=""),row.names=F)
    
  }
}

#explore time for moving to adjancent and non-adjacent zones
fn.time.zone=function(dat,ADJ)
{
  #Release- to first detection
  different.zone=subset(dat,Dist.moved.rel.det>0)
  different.zone$Rel.rec.zone=different.zone$Rec.zone
  different.zone=subset(different.zone,!Rel.zone==Rel.rec.zone)
  different.zone$zone.number=with(different.zone,ifelse(Rel.zone=="North",1,
        ifelse(Rel.zone=="Closed",2,ifelse(Rel.zone=="West",3,
        ifelse(Rel.zone=="Zone1",4,ifelse(Rel.zone=="Zone2",5,6))))))
  different.zone$zone.rec.number=with(different.zone,ifelse(Rel.rec.zone=="North",1,
        ifelse(Rel.rec.zone=="Closed",2,ifelse(Rel.rec.zone=="West",3,
        ifelse(Rel.rec.zone=="Zone1",4,ifelse(Rel.rec.zone=="Zone2",5,6))))))
  different.zone$delta.zone=abs(different.zone$zone.number-different.zone$zone.rec.number)
  if(ADJ=="NO")diff=subset(different.zone,delta.zone>1)
  if(ADJ=="YES")diff=subset(different.zone,delta.zone==1)
  
  diff$days.rel.det=as.numeric(diff$Date.local-diff$ReleaseDate)
  
  STORE=vector('list',4)
  names(STORE)=c('# indiv. moving to a non-adjancent zone in < 1 week',
                 '# indiv. moving to a non-adjancent zone in < 1 month',
                 '# indiv. moving to a non-adjancent zone in < 1 year',
                 'Min # of days for moving to a non-adjancent zone')
  STORE[1:4]=0
  if(nrow(diff)>0)
  {
    this=subset(diff,days.rel.det<=7)
    n.week=length(unique(this$TagCode))
    STORE[[1]]=n.week
    
    
    this=subset(diff,days.rel.det<=30)
    n.month=length(unique(this$TagCode))
    STORE[[2]]=n.month
    
    this=subset(diff,days.rel.det<=365)
    n.yr=length(unique(this$TagCode))
    STORE[[3]]=n.yr
    
    STORE[[4]]=min(diff$days.rel.det)
  }
    
  #Conseq recaptures thereafter
  different.zone=dat
  different.zone=different.zone[order(different.zone$TagCode,different.zone$DateTime.local),]
  different.zone$Rel.rec.zone=c(NA,different.zone$Rec.zone[1:(nrow(different.zone)-1)])
  
  different.zone$Year.prev=c(NA,different.zone$Year[1:(nrow(different.zone)-1)])
  different.zone$Month.prev=c(NA,different.zone$Month[1:(nrow(different.zone)-1)])
  different.zone$Day.prev=c(NA,different.zone$Day[1:(nrow(different.zone)-1)])
  
  different.zone$Rec.rec.zone=different.zone$Rec.zone
  different.zone=subset(different.zone, TagCode==TagCode.prev & !Rec.rec.zone==Rel.rec.zone)
  
  different.zone$zone.number=with(different.zone,ifelse(Rec.rec.zone=="North",1,
      ifelse(Rec.rec.zone=="Closed",2,ifelse(Rec.rec.zone=="West",3,
      ifelse(Rec.rec.zone=="Zone1",4,ifelse(Rec.rec.zone=="Zone2",5,6))))))
  
  different.zone$zone.rec.number=with(different.zone,ifelse(Rel.rec.zone=="North",1,
      ifelse(Rel.rec.zone=="Closed",2,ifelse(Rel.rec.zone=="West",3,
      ifelse(Rel.rec.zone=="Zone1",4,ifelse(Rel.rec.zone=="Zone2",5,6))))))
  different.zone$delta.zone=abs(different.zone$zone.number-different.zone$zone.rec.number)
  if(ADJ=="NO")diff=subset(different.zone,delta.zone>1)
  if(ADJ=="YES")diff=subset(different.zone,delta.zone==1)
  
  diff$days.rel.det=round((diff$Year-diff$Year.prev)*365+(diff$Month-diff$Month.prev)*30.5+(diff$Day-diff$Day.prev))
  
  if(nrow(diff)>0)
  {
    this=subset(diff,days.rel.det<=7)
    n.week=length(unique(this$TagCode))
    STORE[[1]]=STORE[[1]]+n.week
    
    this=subset(diff,days.rel.det<=30)
    n.month=length(unique(this$TagCode))
    STORE[[2]]=STORE[[2]]+n.month
    
    this=subset(diff,days.rel.det<=365)
    n.yr=length(unique(this$TagCode))
    STORE[[3]]=STORE[[3]]+n.yr 
  
    STORE[[4]]=min(c(STORE[[4]],min(diff$days.rel.det)))
  }

    return(STORE)
}
fn.plt=function(what,xLAB)
{
  a=barplot(what,horiz=T,xlab="",main=xLAB,ylab="",col=Singl.col,cex.main=1.75,cex.axis=1.5)
  axis(2,a,c("Dusky","Sandbar","Gummy","Whiskery"),las=1,cex.axis=1.5)
  box()
}

sps=c("Dusky","Thickskin","Gummy","Whiskery")
Time.adj.zn=vector('list',length(sps))
names(Time.adj.zn)=sps

  #Movement among non-adjacent zones
for(s in 1:length(sps))Time.adj.zn[[s]]=fn.time.zone(subset(Tagging.pop.dyn,Species==sps[s]),ADJ="NO")

tiff(file=handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement/Move_among_non_adjacent_zones.tiff"),
     width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(3,1),mai=c(.5,.8,.25,.1),mgp=c(1,.7,0))
fn.plt(cbind(Time.adj.zn[[1]][2],Time.adj.zn[[2]][2],Time.adj.zn[[3]][2],Time.adj.zn[[4]][2]),
       "Number of individuals that moved to non-adjacent areas in one month")
fn.plt(cbind(Time.adj.zn[[1]][3],Time.adj.zn[[2]][3],Time.adj.zn[[3]][3],Time.adj.zn[[4]][3]),
       "Number of individuals that moved to non-adjacent areas in one year")
fn.plt(cbind(Time.adj.zn[[1]][4],Time.adj.zn[[2]][4],Time.adj.zn[[3]][4],Time.adj.zn[[4]][4]),
       "Minimum number of days to move to a non-adjacent area")
dev.off()

#Movement among adjacent zones
for(s in 1:length(sps))Time.adj.zn[[s]]=fn.time.zone(subset(Tagging.pop.dyn,Species==sps[s]),ADJ="YES")
tiff(file=handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement/Move_among_adjacent_zones.tiff"),
     width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(3,1),mai=c(.5,.8,.25,.1),mgp=c(1,.7,0))
fn.plt(cbind(Time.adj.zn[[1]][2],Time.adj.zn[[2]][2],Time.adj.zn[[3]][2],Time.adj.zn[[4]][2]),
       "Number of individuals that moved to an adjacent area in one month")
fn.plt(cbind(Time.adj.zn[[1]][3],Time.adj.zn[[2]][3],Time.adj.zn[[3]][3],Time.adj.zn[[4]][3]),
       "Number of individuals that moved to an adjacent area in one year")
fn.plt(cbind(Time.adj.zn[[1]][4],Time.adj.zn[[2]][4],Time.adj.zn[[3]][4],Time.adj.zn[[4]][4]),
       "Minimum number of days to move to an adjacent area")
dev.off()

  #define coordinates of receivers
      #1 unify different receivers in same stations
Receivers$Station=paste(Receivers$latitude,Receivers$longitude) 
Receivers=Receivers[order(Receivers$Station),]

      #2 extract location of each station, removing duplicates
STATIONS=Receivers[!duplicated(Receivers$Station),]
STATIONS=subset(STATIONS,longitude<=150)

Receiverlat=list(subset(STATIONS,Area=="North.WA")$latitude,NA,
                 subset(STATIONS,Area=="South.WA")$latitude)
Receiverlong=list(subset(STATIONS,Area=="North.WA")$longitude,NA,
                  subset(STATIONS,Area=="South.WA")$longitude)

    #define coordinates of plots
North.WA.lat=c(-23.5,-21.5); North.WA.long=c(113.5,114.25)
South.WA.lat=c(-35.5,-31.5);# South.WA.long=c(114,119)
South.WA.long=c(114,124)
OZ.lat=c(-44.5,-11);OZ.long=c(110,155)

plotlat=list(North.WA.lat,OZ.lat,South.WA.lat)
plotlong=list(North.WA.long,OZ.long,South.WA.long)

fn.seq=function(Range)seq(Range[1]+1,Range[2]-1)
Lat.seq=list(c(-23,-22.5,-22),fn.seq(OZ.lat),fn.seq(South.WA.lat))
Long.seq=list(c(113.5,113.75,114,114.25),fn.seq(OZ.long),fn.seq(South.WA.long))

LONG.lst=list(c(113.5,114.5),c(114,119),c(112.5,119))
LAT.lst=list(c(-23.25,-21.5),c(-35.5,-31.5),c(-36,-21.5))

fn.table=function(y)
{
  id=match(c("Mn.REL","Yr.REL","ZN.REL","Mn.REC","Yr.REC","ZN.REC"),names(y))
  return(y[,id])
}

colfunc <- colorRampPalette(c("chartreuse4","cyan2","firebrick1"))  #color gradient
#colfunc <- colorRampPalette(c("gold","orangered", "red2"))  
fn.see=function(x,y,XLIM,YLIM)
{    
 m <- rbind(c(0, 0.55, 0.5, 1),
             c(0, 0.55, 0, 0.5),
             c(0.6, 1, 0, 1))
  split.screen(m)
   
  for(i in 1:3)
  {
    screen(i)
    par(mar = c(0, 0, 0.1, 0),mgp=c(.1, 0.15, 0))
    xlm=XLIM[[i]]
    ylm=YLIM[[i]]
    plotMap(worldLLhigh, xlim=xlm,ylim=ylm,plt = c(.001, 1, 0.075, 1),
            col="grey88",xlab="",ylab="",axes=F)

    polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
    polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
    
    if(i==3)
    {
      axis(side = 1, at =round(xlm[1]):xlm[2], labels = round(xlm[1]):xlm[2], tcl = .5,las=1,cex.axis=1.2)
      axis(side = 2, at = round(ylm[1]):ylm[2], labels = -(round(ylm[1]):ylm[2]),tcl = .5,las=2,cex.axis=1.2)
      plot(JA_Northern_Shark,ylim=c(-39,-11),xlim=c(108,130),add=T,col="grey75")
      plot(WA_Northern_Shark,add=T,col="grey50")
      plot(WA_Northern_Shark_2,add=T,col="grey85")
      plot(WCDGDLL,add=T,col="white")
      plot(SDGDLL_zone1,add=T,col="grey70")
      plot(SDGDLL_zone2,add=T,col="grey50")
    }
    for(r in 1:length(Receiverlong))points(Receiverlong[[r]],Receiverlat[[r]],col="grey30",pch=20,cex=.75)  #receiver location
    
    if(!i==3)
    {
      #Julian=1+x$Date.local$yday+365*(x$Year-x$Year.rel[1])-x$Date.local[1]$yday
      Julian=as.numeric(x$Date.local-x$ReleaseDate)
      CoLs=colfunc(max(Julian))
      
      points(x$ReleaseLongitude,x$ReleaseLatitude,pch=4,cex=1.25,col="chartreuse4",lwd=3)
      points(x$Long.rec,x$Lat.rec,pch=20,col=CoLs[Julian],cex=2)
    }
    if(i==3)
    {
      #Julian=1+y$Date.local$yday+365*(y$Year-x$Year.rel[1])-y$Date.local[1]$yday
      Julian=as.numeric(y$Date.local-y$ReleaseDate)
      CoLs=colfunc(max(Julian))
      
      arrows(y$LONG.REL,y$LAT.REL,y$LONG.REC,y$LAT.REC,col=CoLs[Julian],lwd=2,length=0.1,angle=20)
      text(117,-24.5,"Days from release",cex=1.25)
      color.legend(117,-25,118,-30,round(seq(1,max(Julian),length.out=5)),CoLs,gradient="y")
      mtext("    Latitude (ºS)",side=2,line=2,las=3,cex=2)
      mtext("Longitude (ºE)",side=1,line=1.5,cex=2)
      points(y$ReleaseLongitude,y$ReleaseLatitude,pch=4,cex=1.25,col="chartreuse4",lwd=3)
    }
      
    box(lwd=1)
   }
  close.screen(all = TRUE)
}

setwd(handl_OneDrive("Analyses/Acoustic_tagging/Acoustic_outputs_pop.dyn"))

# tag=29599
# tiff(file=paste("Sandbar shark.tag code.",tag,".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
# fn.see(x=subset(Store.group$Thickskin$DAT.all,TagCode==tag),y=subset(Store.group$Thickskin$DAT.zn,TagCode==tag),
#        XLIM=LONG.lst,YLIM=LAT.lst)
# dev.off()
# write.csv(fn.table(y=subset(Store.group$Thickskin$DAT.zn,TagCode==tag)),paste("Sandbar shark.tag code.",tag,".csv",sep=""),row.names=F)
 
# tag=31062
# tiff(file=paste("Dusky shark.tag code.",tag,".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
# fn.see(x=subset(Store.group$Dusky$DAT.all,TagCode==tag),y=subset(Store.group$Dusky$DAT.zn,TagCode==tag),
#        XLIM=LONG.lst,YLIM=LAT.lst)
# dev.off()
# write.csv(fn.table(y=subset(Store.group$Dusky$DAT.zn,TagCode==tag)),paste("Dusky shark.tag code.",tag,".csv",sep=""),row.names=F)

#Explore patterns for pop dyn model data
fn.explr.pop=function(dat)
{
  dat=dat[order(dat$TagCode,dat$DateTime.local),]
  dat.f=subset(dat,Sex=="F")
  dat.m=subset(dat,Sex=="M")
  Y=range(c(dat$ReleaseLatitude,dat$Lat.rec),na.rm=T)
  X=range(c(dat$ReleaseLongitude,dat$Long.rec),na.rm=T)
  
  par(mfcol=c(2,1),mai=c(1,1.1,.175,.1),las=1)
  plot(X,Y,col='transparent',ylab="Lat",xlab="Long",cex.lab=1.75,cex.axis=1.25)
  with(subset(dat.f,!is.na(Latitude.prev)),arrows(Long.rec,Lat.rec,
              Longitude.prev,Latitude.prev,col="pink"))
  with(subset(dat.f,is.na(Latitude.prev)),arrows(Long.rec,Lat.rec,
              ReleaseLongitude,ReleaseLatitude,col="pink"))
  
  plot(X,Y,col='transparent',ylab="Lat",xlab="Long",cex.lab=1.75,cex.axis=1.25)
  with(subset(dat.m,!is.na(Latitude.prev)),arrows(Long.rec,Lat.rec,
             Longitude.prev,Latitude.prev,col="blue"))
  with(subset(dat.m,is.na(Latitude.prev)),arrows(Long.rec,Lat.rec,
             ReleaseLongitude,ReleaseLatitude,col="blue"))
  mtext(paste(unique(dat$Species)),3,outer=T,line=-3,cex=2)
}
for(i in 1:length(Store.group))
{
  tiff(file=paste(Pop.din.sp[i],".All.tiff"),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  fn.explr.pop(subset(Tagging.pop.dyn,Species==Pop.din.sp[i]))
  dev.off()
  
  tiff(file=paste(Pop.din.sp[i],".Pop.dyn.selected.tiff"),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  fn.explr.pop(Store.group[[i]][5]$DAT)
  dev.off()
}


#8. -- Home range (Kernel Density) for Ningaloo ---
setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement"))

#number of days detected for 4 species out of total number of days monitored
fn.n.days.monitor=function(SPec,SPEC.nm,txt.CX)
{
  Dat=subset(Detections,Species==SPec & !(TagCode.original%in%only.recaptured))
  Tgs=unique(Dat$TagCode)
  n.Tgs=length(Tgs)
  b=subset(Total.time.monitored,TagCode%in%Tgs,select=c(TagCode,days.mon))
  a=matrix(nrow=length(Tgs),ncol=2)
  for(e in 1:n.Tgs)
  {
    qq=subset(Dat,TagCode==Tgs[e])
    a[e,]=cbind(Tgs[e],length(unique(qq$Date.local)))
  }
  a=as.data.frame(a) 
  names(a)=c("TagCode","Days.detected")
  
  a=merge(a,b,by="TagCode")
  a$Days.detected=as.numeric(as.character(a$Days.detected))
  a$Prop=a$Days.detected/a$days.mon
  a=a[order(a$days.mon),]
  plot(1:n.Tgs,a$Prop,type='h',ylab="",xlab="",xaxt="n",lwd=2,col=Singl.col,ylim=c(0,max(a$Prop,na.rm=T)*1.1))
  axis(1,1:n.Tgs,a$TagCode,las=3,cex.axis=.65)
  mtext(SPEC.nm,3,cex=0.9)
  text(1:n.Tgs,a$Prop,a$days.mon,pos = 3,cex=txt.CX,srt=30)
}
All.spec=these.species[2:5]
SPEC.nms=c("Dusky","Sandbar","Gummy","Whiskery")

tiff(file="Proportion.days.detect.Gum.Whi.tiff",width = 2400, height = 1800,units = "px", res = 300,compression = "lzw")
#par(mfcol=c(2,2),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.75, 0))
par(mfcol=c(2,1),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.7, 0),las=1)
for(p in 3:length(All.spec)) fn.n.days.monitor(All.spec[p],SPEC.nms[p],txt.CX=.65)
mtext("Porportion of days detected",2,outer=T,cex=1.25,line=0.65,las=3)
dev.off()

tiff(file="Proportion.days.detect.Dus.San.tiff",width = 2400, height = 1800,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,1),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.7, 0),las=1)
for(p in 1:2) fn.n.days.monitor(All.spec[p],SPEC.nms[p],txt.CX=.425)
mtext("Porportion of days detected",2,outer=T,cex=1.25,line=0.65,las=3)
dev.off()


  #criteria for calculating: at least detected 5 days within the month
Criteria=5 #at least 5 days detected per month to be selected for that month's analysis
  
fn.plot=function(YLM,XLIM)plot(STATIONS$longitude,STATIONS$latitude,ylim=YLM,xlim=XLIM,pch=19,ylab="",xlab="")
Tgt.spec=these.species[2:3]     #species to consider
Y.range=c(-23.14389,-21.81604)   #spatial range of observations
X.range=c(113.5427,113.9585)
  
#Number of days per month per tag and plotting
do.this="NO"
fn.pre.krnl=function(SPec)
{
    Dat=subset(Detections,Species==SPec & Latitude>=Y.range[1] & Latitude<=Y.range[2])
    Tgs=unique(Dat$TagCode)
    n.Tgs=length(Tgs)
    
    #1. plot detections in Ningaloo for each shark
    fn.see.detec=function(TG,SUBS,start,end)
    {
      dat=subset(Dat,TagCode==TG)
      if(SUBS=="Y") dat=subset(dat,Date.local>=start & Date.local<=end)
      fn.plot(Y.range,X.range)
      points(dat$Longitude,dat$Latitude,col=2,cex=3)
      legend("topleft",paste(TG),bty='n',cex=1.5,text.col=4)
    }
    par(mfcol=c(6,5),mai=c(.1,.1,.1,.1))  
    for(t in 1:n.Tgs)fn.see.detec(TG=Tgs[t],"N",START,END)
    
    #2. table of detections by month and SerialNumber for each tagcode
    Dat$hits=1
    Table1=aggregate(hits~Month+SerialNumber+TagCode,Dat,sum)
    Table1=reshape(Table1,v.names = "hits", idvar = c("TagCode","SerialNumber"),
                   timevar = "Month", direction = "wide")
    names(Table1)[match(c("hits.11","hits.10","hits.8","hits.9","hits.7","hits.6","hits.1",
                          "hits.2","hits.5","hits.3","hits.4","hits.12"),names(Table1))]=
      c("Nov","Oct","Aug","Sep","Jul","Jun","Jan","Feb","May","Mar","Apr","Dec")
    Table1=Table1[,match(c("SerialNumber","TagCode","Jan","Feb","Mar","Apr","May","Jun","Jul",
                           "Aug","Sep","Oct","Nov","Dec"),names(Table1))]
    Table1[is.na(Table1)]=0
    
    #3. table of detections by month and tagcode
    Table2=table(Dat$Month,Dat$TagCode)
    
    
    #4. Kernel density per month
    Dat$Month1=factor(Dat$Month,levels=1:12)
    Dat$Day1=factor(Dat$Day,levels=1:31)
    b=vector('list',n.Tgs)
    names(b)=Tgs
    
    for(l in 1:n.Tgs)
    {
      q=subset(Dat,TagCode==Tgs[l])
      q1=aggregate(hits~Year+Month1+Day1,q,sum)
      q1$hits[q1$hits>0]=1 
      q1=aggregate(hits~Month1,q1,sum)
      q1$TagCode=Tgs[l]
      names(q1)[1:2]=c("Month","Days")
      b[[l]]=q1
    }
    
    b=do.call(rbind,b)
    b=reshape(b,v.names = "Days", idvar ="TagCode",timevar = "Month", direction = "wide")
    Mns=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    names(b)[2:ncol(b)]=Mns
    b[is.na(b)]=0
    b$Total.days=rowSums(b[2:13])
    
    return(list(Detections.month.receiver=Table1,Detections.month=Table2,Days.per.month=b))
  }
Store.stuff=vector("list",length(Tgt.spec))
names(Store.stuff)=Tgt.spec
if(do.this=="YES")
{
  for(q in 1:length(Tgt.spec))Store.stuff[[q]]=fn.pre.krnl(Tgt.spec[q])
  
  #export number of days per month (note that 4 years combined so month could be upto 120 days max)
  Store.stuff$Dusky$Days.per.month$Species="Dusky"
  Store.stuff$Thickskin$Days.per.month$Species="Sandbar"
  Ningaloo.days.monitored=rbind(Store.stuff$Dusky$Days.per.month,Store.stuff$Thickskin$Days.per.month)
  Ningaloo.days.monitored=Ningaloo.days.monitored[order(Ningaloo.days.monitored$Species,Ningaloo.days.monitored$Total.days),]
  write.csv(Ningaloo.days.monitored,"Ningaloo.days.monitored.csv",row.names=F)
  
  #do word table
  Ning.d.mon.word=Ningaloo.days.monitored
  THIS.nms=c("Species","TagCode","Jan","Feb","Mar","Apr","May","Jun",
             "Jul","Aug","Sep","Oct","Nov","Dec","Total.days")
  Ning.d.mon.word=Ning.d.mon.word[,match(THIS.nms,names(Ning.d.mon.word))]
  Scenarios.tbl(WD=getwd(),Tbl=Ning.d.mon.word,Doc.nm="Ningaloo.days.monitored",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
                HEDR=THIS.nms,HEDR.cols=rep(1,length(THIS.nms)),HEDR2=NA,HEDR3=NA)
  
}



#9 -- Natal migration ---
par(new=T)
fn.natal.mig=function(dat) unique(subset(dat,!Rel.zone==Rec.zone)$TagCode)
Migrating.dusky=fn.natal.mig(Store.group$Dusky$DAT.all)
Migrating.sandbar=fn.natal.mig(Store.group$Thickskin$DAT)
colfunc <- colorRampPalette(c("chartreuse4","cyan2","firebrick1"))
fn.natal.mig=function(shks)
{
  dat=subset(Detections,TagCode%in%shks)
  for(i in 1:length(shks))
  {
    d=subset(dat,TagCode==shks[i])
    if(nrow(d)>1)
    {
      d$Mn.prev=c(NA,d$Month[1:(length(d$Month)-1)])
      d$Mn.prev[1]=d$Month.rel[1]
      d$Yr.prev=c(NA,d$Year[1:(length(d$Year)-1)])
      d$Yr.prev[1]=d$Year.rel[1]
      
      d$Label=with(d,ifelse(!Mn.prev==Month,paste(Month,"/",Year),NA))
      d$Label[1]=paste(d$Month[1],"/",d$Year[1])
      
      yli=c(range(d$ReleaseLatitude,d$Latitude))
      xli=c(range(d$ReleaseLongitude,d$Longitude))
      plot(d$ReleaseLongitude,d$ReleaseLatitude,ylim=yli,main=shks[i],
           xlim=xli,ylab="Latitude",xlab="Longitude",pch=4,cex=2,lwd=2)
      
      Julian=as.numeric(d$Date.local-d$ReleaseDate)
      CoLs=colfunc(max(Julian))
      
      arrows(d$ReleaseLongitude,d$ReleaseLatitude,d$Longitude[1],d$Latitude[1],col=CoLs[Julian][1],lwd=2,length=0.1,angle=20)
      arrows(d$Longitude.prev,d$Latitude.prev,d$Longitude,d$Latitude,col=CoLs[Julian],lwd=2,length=0.1,angle=20)
      here=c(yli[2],yli[2]*1.05)
      there= c(xli[2]*.999,xli[2]*0.9999) 
      color.legend(there[1],here[1],there[2],here[2],round(seq(1,max(Julian),length.out=5)),CoLs,gradient="y")
      text(d$Longitude*1.0005,d$Latitude,d$Label)
      text(d$ReleaseLongitude[1],d$ReleaseLatitude[1],paste(d$Month.rel[1],"/",d$Year.rel[1]))
      
    }
    
  }
}
fn.natal.mig(Migrating.dusky)


#time line with seasons as polygons
Solstice=seq(as.POSIXlt("2011-06-01"),as.POSIXlt("2016-12-09"),by="quarter")

names(Solstice)=c(rep(c("Win","Spr","Su","Au"),5),c("Win","Spr","Su"))   

function.timeline=function(shks,Spec,SPEC,FL_50,FL_50.male)
{
  DAT=subset(Detections,TagCode%in%shks & Species== Spec & !Area.release=="SA")
  
  DAT$zone=as.character(with(DAT,
     ifelse(Longitude>=116.5 & Latitude<=(-26),"Zone2",
     ifelse(Longitude<116.5 & Latitude<=(-33),"Zone1",
     ifelse(Latitude>(-33) & Latitude<=(-26) & Longitude<116.5,"West",
     ifelse(Latitude>(-26) & Longitude<114.833,"Closed.ningaloo",
     ifelse(Latitude>(-22) & Longitude>=114.833 & Longitude<123.75,"North",
     ifelse(Latitude>(-22) & Longitude>=123.75,"Joint",NA))))))))
  
  DAT$zone=with(DAT,ifelse(Latitude>(-33) & 
     Latitude<=(-31) & Longitude>=114.8476 & Longitude<116,"Closed.metro",zone))
  
  
  DAT$zone=with(DAT,ifelse(Longitude>=129 & Latitude<=(-26),"SA",zone))
  First.rel=min(DAT$ReleaseDate)
  DAT$Julian=as.numeric(DAT$Date.local-First.rel)
  if(max(DAT$Julian)>10000) DAT$Julian=as.numeric(DAT$Date.local-First.rel)/(24*3600)
  DAT$Julian.release=as.numeric(DAT$ReleaseDate-First.rel)/(24*3600)
  DAT$COLS=with(DAT,ifelse(zone=="Closed",1,ifelse(zone=="West",2,
                ifelse(zone=="Zone1",3,ifelse(zone=="Zone2",4,
                ifelse(zone=="North",5,ifelse(zone=="Joint",6,NA)))))))
   
  DAT$zone.rel=as.character(with(DAT,
    ifelse(ReleaseLongitude>=116.5 & ReleaseLatitude<=(-26),"Zone2",
    ifelse(ReleaseLongitude<116.5 & ReleaseLatitude<=(-33),"Zone1",
    ifelse(ReleaseLatitude>(-33) & ReleaseLatitude<=(-26) & ReleaseLongitude<116.5,"West",
    ifelse(ReleaseLatitude>(-26) & ReleaseLongitude<114,"Closed",
    ifelse(ReleaseLatitude>(-26) & ReleaseLongitude>=114 & ReleaseLongitude<123.75,"North",
    ifelse(ReleaseLatitude>(-26) & ReleaseLongitude>=123.75,"Joint",NA))))))))
  
  DAT$COLS.rel=with(DAT,ifelse(zone.rel=="Closed",1,ifelse(zone.rel=="West",2,
    ifelse(zone.rel=="Zone1",3,ifelse(zone.rel=="Zone2",4,
    ifelse(zone.rel=="North",5,ifelse(zone.rel=="Joint",6,NA)))))))
  
  
  shks=unique(DAT$TagCode)
  N=length(shks)
  names(shks)=1:N
  
  #add size 
  size=subset(DAT,Species==SPEC) 
  size=size[,-match("Species",names(size))]
  #DAT=merge(DAT,size,by=c("TagCode","Year.rel","Month.rel"),all.x=T)

    #order by release date and separate by sex
  fems=subset(DAT,Sex=='F')    
  #a=fems[order(fems$ReleaseDate),]
  a=fems[order(fems$FL),]
  a=a[!duplicated(a$TagCode),]
  fem.shks=a$TagCode
  n.fems=length(fem.shks)
  
  mals=subset(DAT,Sex=='M')
  #a=mals[order(mals$ReleaseDate),]
  a=mals[order(mals$FL),]
  a=a[!duplicated(a$TagCode),]  
  mal.shks=a$TagCode
  n.mals=length(mal.shks)
  
    #add maturity
  fems$Mature=with(fems,ifelse(FL>=FL_50,1,0))
  mals$Mature=with(mals,ifelse(FL>=FL_50.male,1,0))
  
  #a=fems[order(fems$ReleaseDate),]
  a=fems[order(fems$FL),]
  a=a[!duplicated(a$TagCode),]
  fem.mat=a$Mature
  names(fem.mat)=a$TagCode
  fem.mat=ifelse(fem.mat==1,"*","")
  
  #a=mals[order(mals$ReleaseDate),]
  a=mals[order(mals$FL),]
  a=a[!duplicated(a$TagCode),]
  mal.mat=a$Mature
  names(mal.mat)=a$TagCode
  mal.mat=ifelse(mal.mat==1,"*","")
  
    #setup plot
  par(xpd=T,mgp=c(1,.7,0))
  plot(rep(min(DAT$Julian.release):max(DAT$Julian),N),
  rep(1:N,length(min(DAT$Julian.release):max(DAT$Julian))),ylab="",xlab="",
  col='transparent',xlim=c(0,max(DAT$Julian)),ylim=c(1,N),yaxt='n',xaxt='n')
  text(1,-10,First.rel)    #add first date
  arrows(1,-8,1,-1,lwd=2,length=0.1,angle=20)
  
  yr.u=sort(unique(DAT$Year))
  yr.u=yr.u[-1]
  YRS=paste(yr.u,"01-01",sep="-")
  yrs=as.numeric(as.POSIXlt(as.character(YRS))-First.rel)
  ids=yrs<=max(DAT$Julian)
  yrs=yrs[ids]
  YRS=YRS[ids]
  axis(1,yrs,YRS)
  
  #add seasons
  Winter=rgb(0, 0, 1, 0.2);Spring=rgb(0, 1, 0, 0.2);
  Summer=rgb(0.9, 0.9, 0.1, 0.2);Autumn=rgb(0.9, 0.3, 0.1, 0.2)
  COLS=names(Solstice)
  COLS=ifelse(COLS=="Win",Winter,ifelse(names(Solstice)=="Spr",Spring,
            ifelse(names(Solstice)=="Su",Summer,ifelse(names(Solstice)=="Au",Autumn,NA))))
  Polys=as.numeric(Solstice-First.rel)
  names(Polys)=Solstice
  Polys=subset(Polys,Polys<=max(DAT$Julian))
  id=match(names(Polys),as.character(Solstice))
  COLS=COLS[id]  
  for(t in 1:(length(Polys)-1))
  {
    polygon(c(Polys[t],Polys[t+1],Polys[t+1],Polys[t]),c(1,1,N,N),col=COLS[t],border=COLS[t])
  }
  if(Polys[1]>0)
  {
    Cl=names(Solstice)[match(names(Polys[1]),Solstice)-1]
    Cl=ifelse(Cl=="Win",rgb(0, 0, 1, 0.2),
        ifelse(Cl=="Spr",rgb(0, 1, 0, 0.2),
        ifelse(Cl=="Su",rgb(0.9, 0.9, 0.1, 0.2),rgb(0.9, 0.3, 0.1, 0.2))))
    polygon(c(0,Polys[1],Polys[1],0),c(1,1,N,N),col=Cl,border=Cl)
  }
  if(Polys[length(Polys)]<max(DAT$Julian))
  {
    Cl=names(Solstice)[match(names(Polys[length(Polys)]),Solstice)]
    Cl=ifelse(Cl=="Win",rgb(0, 0, 1, 0.2),
      ifelse(Cl=="Spr",rgb(0, 1, 0, 0.2),
      ifelse(Cl=="Su",rgb(0.9, 0.9, 0.1, 0.2),rgb(0.9, 0.3, 0.1, 0.2))))
    polygon(c(Polys[length(Polys)],max(DAT$Julian),max(DAT$Julian),
              Polys[length(Polys)]),c(1,1,N,N),col=Cl,border=Cl)
   }
      
  for(i in 1:n.fems)
  {
    a=subset(fems,TagCode==fem.shks[i])
    points(a$Julian,rep(i,nrow(a)),pch=19,col=a$COLS,cex=1.25)
    points(a$Julian.release[1]*0.95,i,pch="*",col=a$COLS.rel[1],cex=2)
  }

  axis(2,at=1:n.fems,paste(fem.mat,fem.shks),tcl = .3,cex.axis=.75,hadj=0.8,
       col.axis="hotpink1",las=2)
  
  if(nrow(mals)>0)
  {
    lines(c(0,max(DAT$Julian)),rep(n.fems+0.5,2),lty=2,col="grey60",lwd=1.5)
    
    for(i in 1:n.mals)
    {
      a=subset(mals,TagCode==mal.shks[i])
      points(a$Julian,rep(i+n.fems,nrow(a)),pch=19,col=a$COLS,cex=1.25)
      points(a$Julian.release[1]*0.95,i+n.fems,pch="*",col=a$COLS.rel[1],cex=2)
    }
    axis(2,at=(n.fems+1):(n.mals+n.fems),paste(mal.mat,mal.shks),tcl = 0.3,cex.axis=.8,hadj=0.75,
         col.axis="darkblue",las=2)
    
  }

  legend("top", inset=c(0.1,-0.175),c("Closed","West","Zone1","Zone2","North",
      "Joint"),col=1:6,bty='n',pch=19,horiz=T,cex=1.25,title="Fishing zone")

  legend("topright", inset=c(-0.05,-0.175),c("Winter","Spring","Summer","Autumn"),
         fill=c(Winter,Spring,Summer,Autumn),bty='n',horiz=F,cex=0.75,title="Season")
  
  legend("topleft", inset=c(-0.1,-0.15),c("Release"),bty='n',pch="*",
         col="darkorange1",pt.cex=3,cex=1.25)  
  
  mtext("Tag id",2,outer=T,line=-1.4,cex=2)
   mtext("Date",1,outer=T,line=-2,cex=2) 
}

setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC"))

plot.seasonal.color.code.time.lime="NO"
if(plot.seasonal.color.code.time.lime=="YES")
{
  #migrating duskies
  tiff(file="Outputs_movement/Natal_migration/Dusky_timeline.tiff",width = 2500, height = 2000,units = "px", res = 300,compression = "lzw")
  function.timeline(shks=Migrating.dusky,Spec="Dusky",SPEC="BW",FL_50=FL_0.5.dusky,FL_50.male=FL_0.5.dusky.male)
  dev.off()
  
  #non-migrating duskies
  Dusky.tag=subset(Detections,Species=="Dusky")
  Non.mig.dusky=unique(Dusky.tag$TagCode)[which(!unique(Dusky.tag$TagCode)%in%Migrating.dusky)]
  
  tiff(file="Outputs_movement/Natal_migration/Dusky_timeline_non.migratory.tiff",width = 2500, height = 2000,units = "px", res = 300,compression = "lzw")
  function.timeline(Non.mig.dusky,"Dusky","BW",FL_0.5.dusky,FL_0.5.dusky.male)
  dev.off()
  rm(Dusky.tag)
  
  #migrating sandbars
  tiff(file="Outputs_movement/Natal_migration/Sandbar_timeline.tiff",width = 2500, height = 2000,units = "px", res = 300,compression = "lzw")
  function.timeline(Migrating.sandbar,"Thickskin","TK",FL_0.5.sandbar,FL_0.5.sandbar.male)
  dev.off()
  
  #non-migrating sandbars
  Sandbar.tag=subset(Detections,Species=="Thickskin")
  Non.mig.sandbar=unique(Sandbar.tag$TagCode)[which(!unique(Sandbar.tag$TagCode)%in%Migrating.sandbar)]
  
  tiff(file="Outputs_movement/Natal_migration/Sandbar_timeline_non.migratory.tiff",width = 2500, height = 2000,units = "px", res = 300,compression = "lzw")
  function.timeline(Non.mig.sandbar,"Thickskin","TK",FL_0.5.sandbar,FL_0.5.sandbar.male)
  dev.off()
  rm(Sandbar.tag)
  
}



#10. -- distance correction for corners ---
Detections=Detections[order(Detections$Species,Detections$TagCode,Detections$DateTime.local),]

#update Cons detections for release-first detection
Detections$Dist.moved.conseq.det=with(Detections,ifelse(Dup=="FALSE",Dist.moved.rel.det,Dist.moved.conseq.det))

#Corrections for release-first hit
Detections$Dist.moved.conseq.det=with(Detections,
      #north to south movement
    ifelse(Dup=="FALSE" & ReleaseLatitude>Exmouth[2] & Latitude>(-34.65) & Latitude<(-31.5),
   (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Exmouth)+
    distCosine(Exmouth,Shark.bay)+distCosine(Shark.bay,cbind(Longitude,Latitude)))/1000,
                           
   ifelse(Dup=="FALSE" & ReleaseLatitude>Exmouth[2] & Latitude<(-34.65) & Longitude<(116.8),
    (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Exmouth)+distCosine(Exmouth,Shark.bay)+
     distCosine(Shark.bay,Cape.Leuwin)+distCosine(Cape.Leuwin,cbind(Longitude,Latitude)))/1000,
                                                    
   ifelse(Dup=="FALSE" & ReleaseLatitude>Exmouth[2] & Latitude<(-32) & Longitude>=(116.8),         
    (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Exmouth)+distCosine(Exmouth,Shark.bay)+
     distCosine(Shark.bay,Cape.Leuwin)+distCosine(Cape.Leuwin,Mid.point)+
     distCosine(Mid.point,cbind(Longitude,Latitude)))/1000,
                                                        
   ifelse(Dup=="FALSE" & ReleaseLatitude<=Exmouth[2] & ReleaseLatitude>Shark.bay[2] & Latitude>(-34.65) & Latitude<(-31.5),
    (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Shark.bay)+
    distCosine(Shark.bay,cbind(Longitude,Latitude)))/1000,
                                                               
   ifelse(Dup=="FALSE" & ReleaseLatitude<=Exmouth[2] & ReleaseLatitude>Shark.bay[2] & Latitude<(-34.65) & Longitude<(116.8),
    (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Shark.bay)+
    distCosine(Shark.bay,Cape.Leuwin)+distCosine(Cape.Leuwin,cbind(Longitude,Latitude)))/1000,
                                                                      
    ifelse(Dup=="FALSE" & ReleaseLatitude<=Exmouth[2] & ReleaseLatitude>Shark.bay[2] & Latitude<(-32) & Longitude>=(116.8),
       (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Shark.bay)+distCosine(Shark.bay,Cape.Leuwin)+
       distCosine(Cape.Leuwin,Mid.point)+distCosine(Mid.point,cbind(Longitude,Latitude)))/1000,
                                                                             
    ifelse(Dup=="FALSE" & ReleaseLatitude<=Shark.bay[2] & ReleaseLatitude>(-34.65) & ReleaseLongitude<Cape.Leuwin[1] & Latitude<(-34.65) & Longitude<(116.8),
        (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Cape.Leuwin)+
        distCosine(Cape.Leuwin,cbind(Longitude,Latitude)))/1000,
                                                                                    
    ifelse(Dup=="FALSE" & ReleaseLatitude<=Shark.bay[2] & ReleaseLatitude>(-34.65) & ReleaseLongitude<Cape.Leuwin[1] & Latitude<(-32) & Longitude>=(116.8),         
       (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Cape.Leuwin)+
       distCosine(Cape.Leuwin,Mid.point)+distCosine(Mid.point,cbind(Longitude,Latitude)))/1000,
                                                                                           
    ifelse(Dup=="FALSE" & ReleaseLatitude<=(-34.65) & ReleaseLongitude< 116.8 & Latitude<(-32) & Longitude>=(116.8),
      (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Mid.point)+
      distCosine(Mid.point,cbind(Longitude,Latitude)))/1000,
                                                                                                  
    #south to north  movement        
    ifelse(Dup=="FALSE" & ReleaseLatitude<=(-33) & ReleaseLongitude>= 116.8 & Latitude<=(-34) & Longitude<(116.8),
       (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Mid.point)+
       distCosine(Mid.point,cbind(Longitude,Latitude)))/1000,
                                                           
    ifelse(Dup=="FALSE" & ReleaseLatitude<=(-33) & ReleaseLongitude>= 116.8 & Latitude>(-34) & Latitude<(-31.5) & Longitude<115,
      (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Mid.point)+distCosine(Mid.point,Cape.Leuwin)+
      distCosine(Cape.Leuwin,cbind(Longitude,Latitude)))/1000,
                                                                                                                
    ifelse(Dup=="FALSE" & ReleaseLatitude<=(-33) & ReleaseLongitude>= 116.8 & Latitude>Shark.bay[2],
       (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Mid.point)+distCosine(Mid.point,Cape.Leuwin)+
       distCosine(Cape.Leuwin,Shark.bay)+distCosine(Shark.bay,cbind(Longitude,Latitude)))/1000,
                                                                                                                       
    ifelse(Dup=="FALSE" & ReleaseLatitude>(-34.65) & ReleaseLongitude< Cape.Leuwin[1] & ReleaseLatitude< Shark.bay[2]& Latitude>Shark.bay[2],
      (distCosine(cbind(ReleaseLongitude,ReleaseLatitude),Shark.bay)+
      distCosine(Shark.bay,cbind(Longitude,Latitude)))/1000,
                                                                                                                          
    Dist.moved.conseq.det))))))))))))))       


      #Corrections for conseq. hits other than release-first hit
a=subset(Detections,Dup=="FALSE")
Detections=subset(Detections,Dup=="TRUE")
Detections$Dist.moved.conseq.det=with(Detections,
  #north to south movement                                      
  ifelse(Dup=="TRUE" & Latitude.prev>Shark.bay[2] & Latitude>(-34.65) & Latitude<(-31.5),
    (distCosine(cbind(Longitude.prev,Latitude.prev),Shark.bay)+
    distCosine(Shark.bay,cbind(Longitude,Latitude)))/1000,
                                             
  ifelse(Dup=="TRUE" & Latitude.prev>Shark.bay[2] & Latitude<(-34.65) & Longitude<(116.8),
    (distCosine(cbind(Longitude.prev,Latitude.prev),Shark.bay)+
     distCosine(Shark.bay,Cape.Leuwin)+distCosine(Cape.Leuwin,cbind(Longitude,Latitude)))/1000,
                                                    
   ifelse(Dup=="TRUE" & Latitude.prev>Shark.bay[2] & Latitude<(-32) & Longitude>=(116.8),
     (distCosine(cbind(Longitude.prev,Latitude.prev),Shark.bay)+distCosine(Shark.bay,Cape.Leuwin)+
     distCosine(Cape.Leuwin,Mid.point)+distCosine(Mid.point,cbind(Longitude,Latitude)))/1000,
                                                           
   ifelse(Dup=="TRUE" & Latitude.prev<=Shark.bay[2] & Latitude.prev>(-34.65) & Longitude.prev<Cape.Leuwin[1] & Latitude<(-34.65) & Longitude<(116.8),
     (distCosine(cbind(Longitude.prev,Latitude.prev),Cape.Leuwin)+
     distCosine(Cape.Leuwin,cbind(Longitude,Latitude)))/1000,
                                                                  
   ifelse(Dup=="TRUE" & Latitude.prev<=Shark.bay[2] & Latitude.prev>(-34.65) & Longitude.prev<Cape.Leuwin[1] & Latitude<(-32) & Longitude>=(116.8),         
     (distCosine(cbind(Longitude.prev,Latitude.prev),Cape.Leuwin)+
      distCosine(Cape.Leuwin,Mid.point)+distCosine(Mid.point,cbind(Longitude,Latitude)))/1000,
                                                                         
   ifelse(Dup=="TRUE" & Latitude.prev<=(-34.65) & Longitude.prev< 116.8 & Latitude<(-32) & Longitude>=(116.8),
     (distCosine(cbind(Longitude.prev,Latitude.prev),Mid.point)+
     distCosine(Mid.point,cbind(Longitude,Latitude)))/1000,
                                                                                
  #south to north  movement        
   ifelse(Dup=="TRUE" & Latitude.prev<=(-33) & Longitude.prev>= 116.8 & Latitude<=(-34) & Longitude<(116.8),
     (distCosine(cbind(Longitude.prev,Latitude.prev),Mid.point)+
     distCosine(Mid.point,cbind(Longitude,Latitude)))/1000,
                                                                                    
  ifelse(Dup=="TRUE" & Latitude.prev<=(-33) & Longitude.prev>= 116.8 & Latitude>(-34) & Latitude<(-31.5) & Longitude<115,
    (distCosine(cbind(Longitude.prev,Latitude.prev),Mid.point)+distCosine(Mid.point,Cape.Leuwin)+
    distCosine(Cape.Leuwin,cbind(Longitude,Latitude)))/1000,
                                                                                              
  ifelse(Dup=="TRUE" & Latitude.prev<=(-33) & Longitude.prev>= 116.8 & Latitude>Shark.bay[2],
    (distCosine(cbind(Longitude.prev,Latitude.prev),Mid.point)+distCosine(Mid.point,Cape.Leuwin)+
    distCosine(Cape.Leuwin,Shark.bay)+distCosine(Shark.bay,cbind(Longitude,Latitude)))/1000,
                                                                                                     
  ifelse(Dup=="TRUE" & Latitude.prev>(-34.65) & Longitude.prev< Cape.Leuwin[1]& Latitude.prev< Shark.bay[2]& Latitude>Shark.bay[2],
    (distCosine(cbind(Longitude.prev,Latitude.prev),Shark.bay)+
    distCosine(Shark.bay,cbind(Longitude,Latitude)))/1000,
                                                                                                            
  Dist.moved.conseq.det)))))))))))       

Detections=rbind(Detections,a)


#11. -- Calculate days straight-line movement ---

    #11.1. days between consecutive detections
Detections=Detections[order(Detections$Species,Detections$TagCode,Detections$DateTime.local),]

Detections$DateTime.local.prev=c(NA,as.character(Detections$DateTime.local)[1:(nrow(Detections)-1)])
Detections$DateTime.local.prev=as.POSIXlt(Detections$DateTime.local.prev)
Detections$days.conseq.det=as.numeric(with(Detections,
                        ifelse(TagCode==TagCode.prev,DateTime.local-DateTime.local.prev,NA)))/(24*3600)

Detections$hours.conseq.det=as.numeric(with(Detections,
                    ifelse(TagCode==TagCode.prev,DateTime.local-DateTime.local.prev,NA)))/(3600)


    #11.2. add rel-first detection to movement between detections
Detections$days.conseq.det=with(Detections,ifelse(Dup=="FALSE",days.rel.det,days.conseq.det))
Detections$hours.conseq.det=with(Detections,ifelse(Dup=="FALSE",days.rel.det*24,hours.conseq.det))



#12. -- Drop tags detected in only 1 receiver only one time ---          
#note: all hits from single detection sharks done days after release
keep.all="YES"
if(keep.all=="NO")
{
  Station_Tags.table=table(Detections$Station,Detections$TagCode)
  Hits.by.tag=sort(colSums(Station_Tags.table))
  Hits.by.Station=sort(rowSums(Station_Tags.table))
  Hits.by.species.by.area=table(Detections$Area,Detections$Species)
  
  drop.false.pos=colSums(Station_Tags.table)
  drop.false.pos=names(drop.false.pos[drop.false.pos<=1])
  Detections.single=subset(Detections,TagCode%in%drop.false.pos)
  Detections=subset(Detections,!(TagCode%in%drop.false.pos))                        
}


  # Check if release < detected dates 
a=subset(Detections,ReleaseDate>Date.local)
Check.These.TagCodes=a[!duplicated(a$TagCode),]


  #Add dummy
Detections$hit=1


#Add fishing zone of release
Detections$Zone.rel=as.character(with(Detections,
  ifelse(ReleaseLongitude>=116.5 & ReleaseLatitude<=(-26),"Zone2",
  ifelse(ReleaseLongitude<116.5 & ReleaseLatitude<=(-33),"Zone1",
  ifelse(ReleaseLatitude>(-33) & ReleaseLatitude<=(-26) & ReleaseLongitude<116.5,"West",
  ifelse(ReleaseLatitude>(-26) & ReleaseLongitude<114.833,"Closed.ningaloo",
  ifelse(ReleaseLatitude>(-22) & ReleaseLongitude>=114.833 & ReleaseLongitude<123.75,"North",
  ifelse(ReleaseLatitude>(-22) & ReleaseLongitude>=123.75,"Joint",NA))))))))
Detections$Zone.rel=with(Detections,ifelse(ReleaseLatitude>(-33) & 
  ReleaseLatitude<=(-31) & ReleaseLongitude>=114.8476 & ReleaseLongitude<116,"Closed.metro",Zone.rel))
Detections$Zone.rel=with(Detections,ifelse(ReleaseLongitude>129,"SA",Zone.rel))
                         
#create dummy var to separate hits in same receiver in less than minHours
Detections$Rec.hours.conseq.det=with(Detections,ifelse(SerialNumber.prev==SerialNumber &                                              
          hours.conseq.det<minHours,"same","diff"))


#13. -- Create detections data sets for each species ---
SPECIES=these.species[-(match("SENTINEL",these.species))]
N.sp=length(SPECIES)

NUM.TAGGED=as.numeric(Table1[which(rownames(Table1)=="N.tag"),])
names(NUM.TAGGED)=colnames(Table1)
NUM.TAGGED=NUM.TAGGED[which(names(NUM.TAGGED)%in%SPECIES)]

Detections.species=vector("list",length=length(SPECIES))
names(Detections.species)=SPECIES

Months.for.plot=Months.for.plot.label=Start.of.Year=First.Release=Last.Hit=Detections.species

Start.project.date=START    #Date start and end of Analysis
End.project.date=END

#create time line data
for(i in 1:N.sp)
{
  dummy=subset(Detections,Species==SPECIES[i])  
  
  # Add Julian day from start of release events
  First.release=max(Start.project.date,min(dummy$ReleaseDate,na.rm=T))
  Last.hit=max(dummy$Date.local)
  
  dummy$Julian=as.numeric(round(dummy$Date.local-First.release))
  dummy$Julian.release=as.numeric(round(dummy$ReleaseDate-First.release))
  
  if(max(dummy$Julian)>10000) dummy$Julian=dummy$Julian/(3600*24)     #convert back to days if in seconds
  if(max(dummy$Julian.release)>10000) dummy$Julian.release=dummy$Julian.release/(3600*24)
  
  Labels=seq(First.release, Last.hit, by = "6 month")
  Label.value=as.numeric(round(Labels-First.release))
  Labels=format(Labels, format = "%Y-%m-%d")
  if(max(Label.value)>10000) Label.value=Label.value/(3600*24)
  
  Detections.species[[i]]=dummy
  Months.for.plot[[i]]=Label.value
  Months.for.plot.label[[i]]=Labels
  First.Release[[i]]=First.release
  Last.Hit[[i]]=Last.hit
}



#14. -- Table 2 FRDC milestone reports ---
  #Number of individuals detected of each species in each area
Detections$Uni.Sp.Arr.Tag=with(Detections,paste(Species,Array,TagCode))
table2=aggregate(hit~Species+Array,Detections[!duplicated(Detections$Uni.Sp.Arr.Tag),],sum)
names(table2)[3]="No. individuals"
table2=reshape(table2, v.names = "No. individuals", idvar = "Species",
                timevar = "Array", direction = "wide")
write.csv(table2,"Outputs_movement/table2.Numb.ind.array.csv",row.names=F)


#15. -- Table 3 FRDC milestone reports ---
  #Number of individuals detected in >1 array
more.than.one=function(SPECIES)
{
  dat=subset(Detections,Species==SPECIES)
  unik=unique(dat$TagCode)
  one=subset(dat, Array=="Ningaloo")
  one=sort(unique(one$TagCode))
  two=subset(dat, Array=="Perth")
  two=sort(unique(two$TagCode))
  three=subset(dat, Array=="Southern.lines")
  three=sort(unique(three$TagCode))
  
  Ning.Perth.South=length(which(one%in%two & one%in%three))
  Ning.Perth=length(which(one%in%two))
  Ning.South=length(which(one%in%three))
  Perth.South=length(which(two%in%three))
  Perth.Ning=length(which(two%in%one))
  South.Perth=length(which(three%in%two))
  South.Ning=length(which(three%in%one))
  return(list(Ning.Perth.South=Ning.Perth.South,Ning.Perth=Ning.Perth,Ning.South=Ning.South,
              Perth.South=Perth.South,Perth.Ning=Perth.Ning,South.Perth=South.Perth,South.Ning=South.Ning))
}
Duksy.multi.detec=unlist(more.than.one("Dusky"))
Sandbar.multi.detec=unlist(more.than.one("Thickskin"))
Gummy.multi.detec=unlist(more.than.one("Gummy"))
Whiskery.multi.detec=unlist(more.than.one("Whiskery"))
Table3=rbind(Duksy.multi.detec,Sandbar.multi.detec,Gummy.multi.detec,Whiskery.multi.detec)
write.csv(Table3,"Outputs_movement/table3.Multi.array.detections.csv",row.names=T)


#16. -- Histogram of all tags deployed ---
SPECIES.san=c("Dusky","Sandbar","Gummy","Whiskery")
  #2.1 Size frequency distribution of all tagged sharks
TAGS$Size=as.numeric(as.character(TAGS$ReleaseLength))*100 # length in cm
HISTO=vector('list',length=N.sp)
names(HISTO)=SPECIES
for(i in 1:N.sp)
{
  data=subset(TAGS,Species2==SPECIES[i])
  histogram=table(data$Sex2,10*floor(data$Size/10)) #create histogram
  rownames(histogram)=c("Females","Males")
  
    #add missing columns
  COLNS=as.numeric(colnames(histogram))
  Rang=range(COLNS)
  SEC=seq(Rang[1],Rang[2],10)
  ID.mis=match(COLNS,SEC)
  SEC.mis=SEC[-ID.mis]
  MIS.len=matrix(0,nrow=2,ncol=length(SEC.mis))
  colnames(MIS.len)=SEC.mis
  histogram=cbind(as.matrix(histogram),MIS.len)
  histogram=histogram[,match(SEC,colnames(histogram))]
  
  HISTO[[i]]=histogram
}
tiff(file="Outputs_movement/Figure2.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),oma=c(2,4,0.001,5),mar=c(2,2,2,2),mgp=c(2, 0.75, 0))
for(i in 1:N.sp)
{
  #COLS=rev(gray(1:nrow(HISTO[[i]])/nrow(HISTO[[i]])))
  COLS=c("hotpink1","blue4")
  numbers=NUM.TAGGED[i]
  barplot(HISTO[[i]], beside = TRUE,ylim=c(0,max(HISTO[[i]])+1),mgp = c(2, 0.6, 0),
        names.arg= as.numeric(colnames(HISTO[[i]])),xlab="",ylab="",cex.lab=1.5, 
        xpd=F,axis.lty=1, axes=T,col=COLS,
          cex.names=1.1,las=1,cex=1.1,main=SPECIES.san[i])
  box()
  if(!i==1) legend("topright",paste("N=",numbers),cex=1.1,bty='n')
  
  if(i==1)
  {
    legend("topleft",c(rownames(HISTO[[i]]),paste("N=",numbers)),fill=c(COLS,"transparent"),
         border=c(1,1,"transparent"),yjust=0, horiz=F,bty="n",cex=1.1)
  }
}
mtext("Numbers",side=2,outer=T,line=.5,font=1,las=0,cex=1.5)
mtext("Fork length (cm)",side=1,outer=T,line=.5,font=1,las=0,cex=1.5)
dev.off()


#16.1 Check if systematic bias in non-detected duskies
aa=subset(TAGS,Species2=="Dusky" & Project.rel=="SMN")
bb=subset(Detections,Species=="Dusky")
bb=bb[!duplicated(bb$TagCode.original),]
aa=aa[which(!aa$Code2%in%bb$TagCode.original),]
plot(aa$ReleaseLongitude2,aa$ReleaseLatitude2)
barplot(table(aa$Sex2,10*round(aa$ReleaseLength*100/10)),legend.text =c("F","M"),beside=T)

#17. -- Map of study area ---

do.map="NO"   #control if doing map or not (takes time to reshape depth)

if(do.map=="YES")
{
  #define coordinates of polygons
  N.WA.long=c(North.WA.long[2], North.WA.long[2], North.WA.long[1], North.WA.long[1])
  N.WA.lat=c(North.WA.lat[2], North.WA.lat[1], North.WA.lat[1], North.WA.lat[2])
  S.WA.long=c(South.WA.long[2], South.WA.long[2], South.WA.long[1], South.WA.long[1])
  S.WA.lat=c(South.WA.lat[2], South.WA.lat[1], South.WA.lat[1], South.WA.lat[2])
  
  #define Perth and Rotnest
  Perth=c(115.866,-31.95)
  Rotnest=c(115.50,-32.02)
  
  #bathymetry
  Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
  xbat=sort(unique(Bathymetry$V1))
  ybat=sort(unique(Bathymetry$V2)) 
  if(!exists("reshaped"))reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
  rm(Bathymetry)
  #legends
  Letter.leg=c("A",NA,"B")
  Letter.leg.coor=cbind(c(113.61,NA,114.45),c(-21.61,NA,-31.85))
  tiff(file="Outputs_movement/Figure1.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  m <- rbind(c(0.0385, 0.3, 0, 1),
             c(0.35, 1, 0.5, 1),
             c(0.3, 1, 0, 0.5))
  split.screen(m)
  for(i in 1:3)
  {
    screen(i)
    par(mar = c(0, 0, 0, 0),mgp=c(.1, 0.15, 0))
    plotMap(worldLLhigh, xlim=plotlong[[i]],ylim=plotlat[[i]],plt = c(.001, 1, 0.075, 1),
            col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    if(i==2)
    {
      polygon(x=N.WA.long,y=N.WA.lat,lwd=2)
      polygon(x=S.WA.long,y=S.WA.lat,lwd=2)
      text(133,-25,("Australia"),col="black", cex=2)
      mtext("Latitude (ºS)",side=2,line=0,las=3,cex=2)
      mtext("Longitude (ºE)",side=1,line=0,cex=2)
      text(115.98,-22.6,("A"),col="black", cex=1.75)
      text(120.6,-30,("B"),col="black", cex=1.75)
      
    }
    if(i==3)
    { 
      text(116.4,Perth[2],"Perth",col="black", cex=1.5)
      polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
      polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
    }
    if(!i==2)
    {
      points(Receiverlong[[i]],Receiverlat[[i]],col='firebrick3',pch=19)  #receiver location
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
}


#18. -- Figure Timeline of hits by day by shark ---
What.COL.proP="Greens"
#What.COL.proP="Blue_reds"
function.timeline=function(DATA,ARRAY,whereLeg,FL_50,FL_50.male,THRESH,CEX.AX,CEX.TOP)
{
  DATA=subset(DATA,Array%in%ARRAY & !TagCode.original%in%only.recaptured)
  DATA$TagCode1=DATA$TagCode
  DATA$Sex=as.character(DATA$Sex)
  n.detected=DATA[!(duplicated(DATA$TagCode)),]
  n.detected=table(n.detected$Sex)
  
  #SA="paleturquoise4"
  SA="grey45"
  names(SA)="SA"
  if(What.COL.proP=="Greens")COL.proP=c(COL.prop,SA)
  if(What.COL.proP=="Blue_reds")COL.proP=funky.cols
  CL=data.frame(rilis.col=COL.proP,Zone.rel=names(COL.proP))
  DATA=merge(DATA,CL,by='Zone.rel',all.x=T)
  names(COL.proP)[c(2,4)]=c("Ningaloo","Metro")
  
  #add maturity
  DATA$Mature=with(DATA,ifelse(Sex=="F" & FL>=FL_50,1,ifelse(Sex=="F" & FL<FL_50,0,
        ifelse(Sex=="M" & FL>=FL_50.male,1,ifelse(Sex=="M" & FL<FL_50.male,0,NA)))))
  
  #dummy to sort by maturity, date release and sex
  DATA=DATA[order(DATA$ReleaseDate,DATA$Sex,DATA$Mature),]
  Mat.or.not=DATA[!duplicated(DATA$TagCode),match(c("TagCode","Mature"),names(DATA))]
  
  #create table
  Table.hits.date1=with(DATA,table(TagCode,as.character(Date.local)))
  Table.hits.date=with(DATA,table(TagCode,Julian))
  
  First.Jul=colnames(Table.hits.date)[1] #this is the first julian date
  Hits.date=colnames(Table.hits.date1)
  
  Table.hits=table(DATA$TagCode)
  Range.hits=range(Table.hits)
  
  # add sex, size and release area
  these.ones=match(c("TagCode","Sex","Zone.rel","Julian.release"),names(DATA))
  noduplicates <- DATA[!duplicated(DATA$TagCode),these.ones]
  noduplicates$Sex=ifelse(as.character(noduplicates$Sex)=="","U",as.character(noduplicates$Sex))
  noduplicates=noduplicates[order(noduplicates$TagCode),]
  if(sum(row.names(Table.hits.date)==noduplicates$TagCode)==nrow(noduplicates))
  {
    Table.hits.date=cbind(noduplicates,as.data.frame.matrix(Table.hits.date))  
  }
  
  #sort detections
  Table.hits.date=merge(Table.hits.date,Mat.or.not,by="TagCode",all.x=T) 
  Table.hits.date=Table.hits.date[order(Table.hits.date$Sex,Table.hits.date$Zone.rel,
                                        Table.hits.date$Mature,Table.hits.date$Julian.release),]
  Table.hits.date$Mature=with(Table.hits.date,ifelse(Mature==1,"*",""))
  Names=Table.hits.date[,match(c("TagCode","Mature","Sex"),colnames(Table.hits.date))]
  Table.hits.date=Table.hits.date[,-match("Mature",colnames(Table.hits.date))]
  
  
  Rel.date.Julian=Table.hits.date$Julian.release
  This.col=match(First.Jul,colnames(Table.hits.date)) 
  
  # create colors by array
  Table.tag.day.area=with(DATA,table(TagCode,Julian,as.character(Array)))
  DIM.area=dim(Table.tag.day.area)[3]
  
  #convert observations to levels (1=Ningaloo, 10=Perth, 5=Southern.lines)
  Replace.number=dimnames(Table.tag.day.area)[[3]]
  Replace.number=ifelse(Replace.number=="Ningaloo",1,ifelse(Replace.number=="Perth",10,
                                                            ifelse(Replace.number=="Southern.lines",5,NA)))        
  for(j in 1:DIM.area) Table.tag.day.area[,,j]=ifelse(Table.tag.day.area[,,j]>0,Replace.number[j],0)
  if(DIM.area==1) Table.tag.day.area=Table.tag.day.area[,,1]
  if(DIM.area==2)Table.tag.day.area=Table.tag.day.area[,,1]+Table.tag.day.area[,,2]
  if(DIM.area==3)Table.tag.day.area=Table.tag.day.area[,,1]+Table.tag.day.area[,,2]+Table.tag.day.area[,,3]
  
  #reorder matrix
  ID=match(Table.hits.date$TagCode,row.names(Table.tag.day.area))
  Table.tag.day.area=Table.tag.day.area[ID,]
  
  #specify colors
  Mat.colors=as.data.frame.matrix(Table.tag.day.area)
  Mat.colors=ifelse(Mat.colors==1,CLS[1],ifelse(Mat.colors==10,CLS[2],ifelse(Mat.colors==5,CLS[3],NA)))
  Mat.colors.legends=c("Ningaloo","Perth","Southern Lines")
  Mat.colors.legends.colors=CLS
  
  # create pch by array
  Mat.pch=as.data.frame.matrix(Table.tag.day.area)
  Mat.pch=ifelse(Mat.pch==1,19,ifelse(Mat.pch==10,19,ifelse(Mat.pch==5,19,NA)))
  Mat.colors.legends.pch=c(19,19,19)
  
  colores=DATA[!duplicated(DATA$TagCode),match(c("TagCode","rilis.col"),names(DATA))]
  colores=merge(subset(Table.hits.date,select=TagCode),colores,by="TagCode")  
  ID=match(Table.hits.date$TagCode,colores$TagCode)
  colores=colores[ID,]
  Y.leg=paste("Tag ID (n=",n.detected[2] ,"males &",n.detected[1],"females)")
  
  bubble.plot.detections(as.numeric(colnames(Table.hits.date[,This.col:ncol(Table.hits.date)])),
                         1:nrow(Table.hits.date),
                         Table.hits.date[,This.col:ncol(Table.hits.date)],3,"Date",Y.leg,Rel.date.Julian,
                         Mat.colors,Mat.pch,Mat.colors.legends,Mat.colors.legends.colors,
                         Mat.colors.legends.pch,"presence",whereLeg,as.character(colores$rilis.col))
  a=subset(COL.proP,COL.proP%in%as.character(unique(DATA$rilis.col)))
  names(a)=ifelse(names(a)=="North","WANCSF",ifelse(names(a)=="West","WCDGDLF",
      ifelse(names(a)=="Metro","Metro closure",names(a))))
  legend("top",names(a),bty='n',pt.cex=1.5,pch="R",col=a,horiz=T,inset=c(0,-.07),cex=CEX.TOP)
  
  Names$Mature=with(Names,ifelse(is.na(Mature),"",Mature))
  FEM=subset(Names,Sex=="F")
  Legnd=paste(FEM$Mature,FEM$TagCode,sep="")
  axis(side = 2, at = 1:nrow(FEM), labels = Legnd, tcl = 0.5,
       cex.axis=CEX.AX,hadj=0.75,col.axis="hotpink1",col.ticks="hotpink1",las=1) # major ticks & labels
  MAL= subset(Names,Sex=="M") 
  Legnd=paste(MAL$Mature,MAL$TagCode,sep="")
  axis(side = 2, at = (nrow(FEM)+1):(nrow(FEM)+nrow(MAL)), labels = Legnd, tcl = 0.5,
       cex.axis=CEX.AX,hadj=0.75,col.axis="blue4",col.ticks="blue4",las=1) 
  
  U= subset(Names,Sex=="U") 
  if(nrow(U)>0)
  {
    Legnd=paste(U$Mature,U$TagCode,sep="")
    axis(side = 2, at = (nrow(FEM)+nrow(MAL)+1):(nrow(FEM)+nrow(MAL)+nrow(U)), labels = Legnd, tcl = 0.5,
         cex.axis=CEX.AX,hadj=0.75,col.axis=1,col.ticks=1,las=1) 
  }
  
  
  #Extract individuals moving to different areas
  DATA$Array.prev=with(DATA,ifelse(Latitude.prev>(-24),"Ningaloo",
                                   ifelse(Latitude.prev<=(-24) & Latitude.prev>(-33),"Perth",
                                          ifelse(Latitude.prev<=(-33),"Southern.lines",NA))))
  
  
  nn=nrow(DATA)
  DATA=DATA[order(DATA$Species,DATA$TagCode,DATA$DateTime.local),]
  DATA$Month.prev=NA
  DATA$Month.prev[1]=DATA$Month.prev[1]
  DATA$Month.prev[2:nn]=DATA$Month[1:(nn-1)]
  DATA$Month.prev=with(DATA,ifelse(TagCode1==TagCode.prev,Month.prev,NA)) 
  DATA$Year.prev=NA
  DATA$Year.prev[1]=DATA$Year.prev[1]
  DATA$Year.prev[2:nn]=DATA$Year[1:(nn-1)]
  DATA$Year.prev=with(DATA,ifelse(TagCode1==TagCode.prev,Year.prev,NA))  
  
  different.areas=subset(DATA,!Area==Area.release & !is.na(Dist.moved.rel.det))  
  different.areas$Array.prev=with(different.areas,ifelse(ReleaseLatitude>(-24),"Ningaloo",
                                                         ifelse(ReleaseLatitude<=(-24) & ReleaseLatitude>(-33),"Perth",
                                                                ifelse(ReleaseLatitude<=(-33),"Southern.lines",NA))))
  different.areas$Month.prev=different.areas$Month.rel
  different.areas$Year.prev=different.areas$Year.rel
  different.areas1=subset(DATA,!Array==Array.prev)
  different.areas1=subset(different.areas1,Dist.moved.conseq.det>THRESH)
  different.areas=rbind(different.areas,different.areas1)
  different.areas$Moved.to=with(different.areas,
                                ifelse(Array.prev=="Ningaloo" & Array=="Southern.lines","Moved.Ningaloo-Southern",           
                                       ifelse(Array.prev=="Ningaloo" & Array=="Perth","Moved.Ningaloo-Perth",     
                                              ifelse(Array.prev=="Southern.lines" & Array=="Ningaloo","Moved.Southern-Ningaloo",           
                                                     ifelse(Array.prev=="Perth" & Array=="Ningaloo","Moved.Perth-Ningaloo",
                                                            ifelse(Array.prev=="Southern.lines" & Array=="Perth","Moved.Southern-Perth",
                                                                   ifelse(Array.prev=="Perth" & Array=="Southern.lines","Moved.Perth-Southern",
                                                                          NA)))))))
  different.areas=subset(different.areas,!is.na(Moved.to),
                         select=c(TagCode1,Sex,Month.prev,Year.prev,Month,Year,Array.prev,Array,Moved.to))
  different.areas=different.areas[order(different.areas$TagCode1,different.areas$Sex,
                                        different.areas$Moved.to,different.areas$Month),]
  
  return(different.areas)
}       

LEGND=c("left","left","topright","bottomright")
CX=c(.65,.65,.8,.8)
CEX.TOp=c(0.85,1.1,1.1,1.1)
NsPs=as.numeric(Table1[rownames(Table1)=="N.det"])
names(NsPs)=colnames(Table1)
POLYS=list(c(-1,-1,NsPs[1]+2,NsPs[1]+2),
           c(-1,-1,NsPs[2]+2,NsPs[2]+2),
           c(0,0,NsPs[3]+1.25,NsPs[3]+1.25),
           c(0.65,0.65,NsPs[4]+.5,NsPs[4]+.5))

#Check if Sanbars constantly detected in Ningaloo are dead
check=c("SS.12","SS.26","SS.24","SS.23","SS.22")
for(i in 1:length(check))
{
  test1=subset(Detections.species[[2]],TagCode==check[i])
  test1$LONLAT=with(test1,paste(Longitude,Latitude))
  test1=test1[!duplicated(test1$LONLAT),]
  test1$colr=substr(test1$Year,4,4)
  tiff(file=paste("Outputs_movement/Time.lines/Sandbar.not.dead/",check[i],".tiff",sep=""),width = 2400, height = 2400,
       units = "px", res = 300,compression = "lzw")
  XLIMsand=c(min(c(test1$Longitude,test1$ReleaseLongitude)),max(c(test1$Longitude,test1$ReleaseLongitude)))
  YLIMsand=c(min(c(test1$Latitude,test1$ReleaseLatitude)),max(c(test1$Latitude,test1$ReleaseLatitude)))
  with(test1,plot(Longitude,Latitude,main=check[i],xlim=XLIMsand,ylim=YLIMsand,col=colr))
  with(test1,points(ReleaseLongitude,ReleaseLatitude,pch=19,col="orange",cex=2)) 
  legend('topleft',paste(unique(test1$Year)),bty='n',text.col=unique(test1$colr))
  dev.off()
}

  #all species timeline
for (s in 1:length(Detections.species))             
{ 
  # plot presence absence time line
  tiff(file=paste("Outputs_movement/Time.lines/Figure3.",SPECIES[s],".pres_abs.tiff",sep=""),width = 2400, height = 2400,
       units = "px", res = 300,compression = "lzw")
  #par(mai=c(.8,.8,.01,.01),las=1,mgp=c(2.3, 1, 0))
  par(mai=c(.05,.5,.4,.01),oma=c(2,1,1,1),las=1,xpd=T,mgp=c(2,.8,0))
  store=function.timeline(Detections.species[[s]],c("Ningaloo","Perth","Southern.lines"),
      LEGND[s],FL_0.5[[s]][1],FL_0.5[[s]][2],Migration.threshold[s],CX[s],CEX.TOp[s])
   axis(side = 1, at = Months.for.plot[[s]],labels = Months.for.plot.label[[s]], 
        tcl = 0.5,cex.axis=.85,padj=-2) # minor ticks
  hand=Months.for.plot[[s]]
  COL.t=rgb(.1,.1,.1,alpha=0.1)
  ADD.polys=function(H1,H2)polygon(x=c(H1,H2,H2,H1),y=POLYS[[s]],lwd=2,col=COL.t,border=COL.t)
  is.even <- function(x) x %% 2 == 0 #get even numbers
  Even=1:length(hand)
  Even=Even[is.even(1:length(hand))]
 #  for (p in 1:length(Even))ADD.polys(hand[Even[p]],hand[Even[p]+1])
  dev.off()

  #report migrations
  if(nrow(store)>0)
  {
    store$time.diff=round(365*with(store,Year+(Month/12)-(Year.prev+(Month.prev/12))),0)
      
    Traj=c("Moved.Ningaloo-Perth","Moved.Ningaloo-Southern","Moved.Perth-Southern",
           "Moved.Perth-Ningaloo","Moved.Southern-Ningaloo","Moved.Southern-Perth")
    names(Traj)=c("From Ningaloo to Perth","From Ningaloo to Southern Lines","From Perth to Southern Lines",
                  "From Perth to Ningaloo","From Southern Lines to Ningaloo","From Southern Lines to Perth")
    if(!s==1)Traj=Traj[match(unique(store$Moved.to),Traj)]
    
    n.graf=length(Traj)
    if(n.graf>2)Mfcol=c(ceiling(n.graf/2),2)
    if(n.graf==2)Mfcol=c(2,ceiling(n.graf/2))
    tiff(file=paste("Outputs_movement/Time.lines/Figure3.migrations",SPECIES[s],".tiff",sep=""),width = 2000, height = 2400,
         units = "px", res = 300,compression = "lzw")
    if(s==1)par(mfcol=Mfcol,mai=c(.3,.3,.1,.01),oma=c(2,2,1,1),las=1,mgp=c(2.3, 0.8, 0))
    if(!s==1)par(mfcol=Mfcol,mai=c(.4,.4,.25,.01),oma=c(2,2,1,1),las=1,mgp=c(2.3, 0.8, 0))
    for (t in 1:length(Traj))
    {
      a=subset(store,Moved.to==Traj[t])
      a$Month=factor(a$Month,levels=1:12)
      TAB=table(a$Sex,a$Month)
      Cols=c("hotpink1","blue4")
      names(Cols)=c("F","M")
      Cols=Cols[which(names(Cols)%in%rownames(TAB))]
#       if(t==1)barplot(TAB,beside=T,col=Cols,cex.axis=1.15,cex.names=1.15,
#                       legend.text=c("Female","Male"),args.legend=list("topleft",
#                       cex=1.75,bty='n',fill=c("hotpink1","blue4")))
#       if(!t==1)barplot(TAB,beside=T,col=Cols,cex.axis=1.15,cex.names=1.15)
      barplot(TAB,beside=T,col=Cols,cex.axis=1.15,cex.names=1.15)
      box()
      if(s==1)mtext(names(Traj)[t],3,line=0,cex=1.05)
      if(!s==1)mtext(names(Traj)[t],3,line=0,cex=1.25)
    }
    mtext("Month arriving",1,outer=T,line=0.2,cex=1.4)
    mtext("Frequency",2,outer=T,las=3,line=0.3,cex=1.4)
    dev.off()
    
  }
  write.csv(store,paste("Outputs_movement/Time.lines/Migrations.table",SPECIES[s],".csv",sep=""),row.names=F)
}                                                      



#18.1 -- Behavioural polymorphism ---   

BP_sp=c("Dusky","Thickskin")
hndl.BP=handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement/Behavioural_Polymorphism/")

#18.1.1 Proportion of days monitored
fn.n.days.monitor_BP=function(SPec,SPEC.nm,txt.CX)
{
  Dat=subset(Detections,Array=="Ningaloo" & Species==SPec & !(TagCode.original%in%only.recaptured))
  Tgs=unique(Dat$TagCode)
  n.Tgs=length(Tgs)
  
  a=vector('list',length=n.Tgs)
  for(n in 1:n.Tgs)
  {
    d=subset(Dat,TagCode==Tgs[n],select=c(TagCode,Sex,FL,ReleaseDate,Date.local))
    d=d[order(d$Date.local),]
    Monitored=d[nrow(d),]
    Monitored=as.numeric(with(Monitored,Date.local-ReleaseDate))
    Detected=length(unique(d$Date.local))
    a[[n]]=data.frame(TagCode=unique(d$TagCode),Sex=unique(d$Sex),FL=unique(d$FL),
                      Prop=Detected/Monitored,days.mon=Monitored)
  }
  a=do.call(rbind,a)
  a=a[order(a$days.mon),]
  
  Y=1.05
  plot(1:n.Tgs,a$Prop,type='h',ylab="",xlab="",xaxt="n",lwd=2,col="grey20",ylim=c(0,Y))
  axis(1,1:n.Tgs,a$TagCode,las=3,cex.axis=.65)
  mtext(SPEC.nm,3,cex=1.25)
  text(1:n.Tgs,a$Prop,a$days.mon,pos = 3,cex=txt.CX,srt=55)
}
tiff(file=paste(hndl.BP,"Appendix1.tiff",sep=""),width = 2400, height = 1800,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,1),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.7, 0),las=1)
for(p in 1:length(BP_sp)) fn.n.days.monitor_BP(SPec=BP_sp[p],SPEC.nm=paste(SPEC.nms[p],"shark"),txt.CX=.6)
mtext("Porportion of days detected",2,outer=T,cex=1.25,line=0.65,las=3)
dev.off()


#Beta regression testing effect of sex and size
library("betareg")

fn.n.beta.reg=function(SPec)
{
  Dat=subset(Detections,Array=="Ningaloo" & Species==SPec & !(TagCode.original%in%only.recaptured))
  Tgs=unique(Dat$TagCode)
  n.Tgs=length(Tgs)
  
  a=vector('list',length=n.Tgs)
  for(n in 1:n.Tgs)
  {
    d=subset(Dat,TagCode==Tgs[n],select=c(TagCode,Sex,FL,ReleaseDate,Date.local,Station))
    d=d[order(d$Date.local),]
    Monitored=d[nrow(d),]
    Monitored=as.numeric(with(Monitored,Date.local-ReleaseDate))
    Detected=length(unique(d$Date.local))
    N.receiver=length(unique(d$Station))
    a[[n]]=data.frame(TagCode=unique(d$TagCode),Sex=unique(d$Sex),FL=unique(d$FL),
                      Prop.det=Detected/Monitored,Monitored=Monitored,N.receiver=N.receiver)
  }
  DAT=do.call(rbind,a)
  DAT=subset(DAT,Monitored>MinDays & !Sex=="U")
  DAT$Sex=as.factor(DAT$Sex)
  
  #beta regression
  BetaReg <- betareg(Prop.det ~ FL + Sex, data = DAT)
  
  return(list(BetaReg=BetaReg,DAT=DAT))
}
MinDays=30
STore.Beta=vector("list",length(BP_sp))
for(p in 1:length(BP_sp)) STore.Beta[[p]]=fn.n.beta.reg(BP_sp[p])

DuskyBetaReg=as.data.frame(summary(STore.Beta[[1]]$BetaReg)[[1]]$mean)
SandbarBetaReg=as.data.frame(summary(STore.Beta[[2]]$BetaReg)[[1]]$mean)

write.csv(round(DuskyBetaReg,2),paste(hndl.BP,"BetaReg_Dusky.csv",sep=""))
write.csv(round(SandbarBetaReg,2),paste(hndl.BP,"BetaReg_Sandbar.csv",sep=""))

with(STore.Beta[[1]]$DAT,plot(log(Prop.det),N.receiver))
with(STore.Beta[[2]]$DAT,plot(log(Prop.det),N.receiver))

#18.1.2 Timelines
Detections.species_BP=Detections.species[match(c("Dusky","Thickskin"),names(Detections.species))]
function.timeline_BP=function(DATA,CEX.AX,CEX.TOP,Example_BP)
{
  DATA$FL=round(DATA$FL*100,0)
  DATA=subset(DATA,Array=="Ningaloo" & !TagCode.original%in%only.recaptured)
  DATA$TagCode1=DATA$TagCode
  DATA$Sex=as.character(DATA$Sex)
  n.detected=DATA[!(duplicated(DATA$TagCode)),]
  n.detected=table(n.detected$Sex)
  
  #SA="paleturquoise4"
  SA="grey45"
  names(SA)="SA"
  if(What.COL.proP=="Greens")COL.proP=c(COL.prop,SA)
  if(What.COL.proP=="Blue_reds")COL.proP=funky.cols
  CL=data.frame(rilis.col=COL.proP,Zone.rel=names(COL.proP))
  DATA=merge(DATA,CL,by='Zone.rel',all.x=T)
  names(COL.proP)[c(2,4)]=c("Ningaloo","Metro")
  COL.proP=rep("tansparent",length(COL.proP))
  
  #dummy to sort by maturity, date release and sex
  DATA=DATA[order(DATA$ReleaseDate,DATA$Sex,DATA$FL),]
  Mat.or.not=DATA[!duplicated(DATA$TagCode),match(c("TagCode","FL"),names(DATA))]
  
  #create table
  Table.hits.date1=with(DATA,table(TagCode,as.character(Date.local)))
  Table.hits.date=with(DATA,table(TagCode,Julian))
  
  First.Jul=colnames(Table.hits.date)[1] #this is the first julian date
  Hits.date=colnames(Table.hits.date1)
  
  Table.hits=table(DATA$TagCode)
  Range.hits=range(Table.hits)
  
  # add sex, size and release area
  these.ones=match(c("TagCode","Sex","Zone.rel","Julian.release"),names(DATA))
  noduplicates <- DATA[!duplicated(DATA$TagCode),these.ones]
  noduplicates$Sex=ifelse(as.character(noduplicates$Sex)=="","U",as.character(noduplicates$Sex))
  noduplicates=noduplicates[order(noduplicates$TagCode),]
  if(sum(row.names(Table.hits.date)==noduplicates$TagCode)==nrow(noduplicates))
  {
    Table.hits.date=cbind(noduplicates,as.data.frame.matrix(Table.hits.date))  
  }
  
  #sort detections
  Table.hits.date=merge(Table.hits.date,Mat.or.not,by="TagCode",all.x=T) 
  Table.hits.date=Table.hits.date[order(Table.hits.date$Sex,Table.hits.date$FL,
                                        Table.hits.date$Zone.rel,Table.hits.date$Julian.release),]
  
  Names=Table.hits.date[,match(c("TagCode","FL","Sex"),colnames(Table.hits.date))]
  Table.hits.date=Table.hits.date[,-match("FL",colnames(Table.hits.date))]
  
  
  Rel.date.Julian=Table.hits.date$Julian.release
  This.col=match(First.Jul,colnames(Table.hits.date)) 
  
  # create colors by array
  Table.tag.day.area=with(DATA,table(TagCode,Julian,as.character(Array)))
  DIM.area=dim(Table.tag.day.area)[3]
  
  #convert observations to levels (1=Ningaloo, 10=Perth, 5=Southern.lines)
  Replace.number=dimnames(Table.tag.day.area)[[3]]
  Replace.number=ifelse(Replace.number=="Ningaloo",1,ifelse(Replace.number=="Perth",10,
                                                            ifelse(Replace.number=="Southern.lines",5,NA)))        
  for(j in 1:DIM.area) Table.tag.day.area[,,j]=ifelse(Table.tag.day.area[,,j]>0,Replace.number[j],0)
  if(DIM.area==1) Table.tag.day.area=Table.tag.day.area[,,1]
  if(DIM.area==2)Table.tag.day.area=Table.tag.day.area[,,1]+Table.tag.day.area[,,2]
  if(DIM.area==3)Table.tag.day.area=Table.tag.day.area[,,1]+Table.tag.day.area[,,2]+Table.tag.day.area[,,3]
  
  #reorder matrix
  ID=match(Table.hits.date$TagCode,row.names(Table.tag.day.area))
  Table.tag.day.area=Table.tag.day.area[ID,]
  
  #specify colors
  Mat.colors=as.data.frame.matrix(Table.tag.day.area)
  Mat.colors=ifelse(Mat.colors==1,"grey55",ifelse(Mat.colors==10,CLS[2],ifelse(Mat.colors==5,CLS[3],NA)))
  Mat.colors.legends=c("Ningaloo","Perth","Southern Lines")
  Mat.colors.legends.colors=CLS
  
  # create pch by array
  Mat.pch=as.data.frame.matrix(Table.tag.day.area)
  Mat.pch=ifelse(Mat.pch==1,19,ifelse(Mat.pch==10,19,ifelse(Mat.pch==5,19,NA)))
  Mat.colors.legends.pch=c(19,19,19)
  
  colores=DATA[!duplicated(DATA$TagCode),match(c("TagCode","rilis.col"),names(DATA))]
  colores=merge(subset(Table.hits.date,select=TagCode),colores,by="TagCode")  
  ID=match(Table.hits.date$TagCode,colores$TagCode)
  colores=colores[ID,]
  colores$rilis.col=rep("black",length(colores$rilis.col))
  Y.leg=paste("Detected individuals (n=",n.detected[2] ,"males &",n.detected[1],"females)")
  bubble.plot.detections=function(x,y,z,scale,xlab,ylab,Rel.date,Mat.colors,Mat.pch,Mat.colors.legends,
                                  Mat.colors.legends.colors,Mat.colors.legends.pch,PROP,where.legend,
                                  Rilis.col,...) 
  {
    n=length(x); ny=length(y)
    xo=outer(x,rep(1,length=length(y)))
    yo=t(outer(y,rep(1,length=length(x))))
    
    if(PROP=="proportion") zo=z/max(z,na.rm=T)*scale
    if(PROP=="presence") zo=ifelse(z>0,1,0)
    
    matplot(xo,yo,type="n",xlab="",ylab="",xaxt='n',yaxt='n',xlim=c(0,max(x)))    
    mtext(xlab,1,line=1.25,cex=1.5)
    mtext(ylab,2,las=3,line=1.8,cex=1.5)
    
    for(i in 1:n)
    {
      points(xo[i,],yo[i,],cex=zo[,i],pch=Mat.pch[,i],col=Mat.colors[,i])
    }
    points(Rel.date,1:ny,col=Rilis.col,cex=1.85,pch="*")
    
    
    if(PROP=="proportion")
    {
      RANGES=unique(c(as.matrix(z)))
      RANGES=pretty(RANGES)
      RANGES=c(1,RANGES[RANGES>0])
      RANGES=RANGES[RANGES<=1500]
      RANGES.prop=RANGES/max(z,na.rm=T)*scale
      legend(where.legend,legend=c(RANGES,Mat.colors.legends),pch=c(Mat.colors.legends.pch),
             pt.cex=c(RANGES.prop,rep(1,length(Mat.colors.legends)),1),bty='n',
             col=c(rep(1,length(RANGES)),Mat.colors.legends.colors))
    }
    
    if(PROP=="presence")
    {
      legend(where.legend,legend=c(Mat.colors.legends),pch=c(Mat.colors.legends.pch),
             pt.cex=1.5,cex=1.25,bty='n',col=c(Mat.colors.legends.colors))
    }
  }
  bubble.plot.detections(as.numeric(colnames(Table.hits.date[,This.col:ncol(Table.hits.date)])),
                         1:nrow(Table.hits.date),
                         Table.hits.date[,This.col:ncol(Table.hits.date)],4,"Date",Y.leg,Rel.date.Julian,
                         Mat.colors,Mat.pch,rep("",3),rep("transparent",3),
                         Mat.colors.legends.pch,"presence","top",as.character(colores$rilis.col))
  FEM=subset(Names,Sex=="F")
  Legnd=FEM$FL
  axis(side = 2, at = 1:nrow(FEM), labels = Legnd, tcl = 0.5,
       cex.axis=CEX.AX,hadj=0.5,col.axis="grey40",col.ticks="grey40",las=1) # major ticks & labels
  MAL= subset(Names,Sex=="M") 
  Legnd=MAL$FL
  axis(side = 2, at = (nrow(FEM)+1):(nrow(FEM)+nrow(MAL)), labels = Legnd, tcl = 0.5,
       cex.axis=CEX.AX,hadj=0.5,col.axis="grey70",col.ticks="grey70",las=1) 
  
  U= subset(Names,Sex=="U") 
  if(nrow(U)>0)
  {
    Legnd=U$FL
    axis(side = 2, at = (nrow(FEM)+nrow(MAL)+1):(nrow(FEM)+nrow(MAL)+nrow(U)), labels = Legnd, tcl = 0.5,
         cex.axis=CEX.AX,hadj=0.5,col.axis=1,col.ticks=1,las=1) 
  }
  
  #Polygons
  for(pp in 1:length(Example_BP))
  {
    Polyd=subset(Table.hits.date,TagCode==Example_BP[pp])
    a=Polyd[1,-match(c("TagCode","Sex","Zone.rel" ,"Julian.release"),names(Polyd))]
    a[a>0]=1
    b=a*as.numeric(names(a))
    b=b[b>0]
    idd=match(Example_BP[pp],Table.hits.date$TagCode)
    if(names(Example_BP[pp])=="Resident")CLP=list(Red=1,Green=0,Blue=0)
    if(names(Example_BP[pp])=="Transient")CLP=list(Red=0,Green=1,Blue=0)
    if(names(Example_BP[pp])=="Wide ranging")CLP=list(Red=0,Green=0,Blue=1)
    PoLy(Polyd$Julian.release,1.025*tail(b,1),floor(idd*.985),ceiling(idd*1.015),CLP)
  }
  
  LEG.col=c(rgb(red=1,green=0,blue=0,alpha=0.2),rgb(red=0,green=1,blue=0,alpha=0.2),rgb(red=0,green=0,blue=1,alpha=0.2))
  legend(-5,nrow(Table.hits.date)*1.1,names(Example_BP),fill=LEG.col,bty='n',horiz=T,cex=1.25)
}
PoLy=function(x,x1,y,y1,CL) polygon(c(x,x1,x1,x),c(rep(y,2),rep(y1,2)),col=rgb(CL$Red,CL$Green,CL$Blue,alpha=0.2),border="transparent")
E.g._Dusky=c("DS.74","DS.13","DS.30")
names(E.g._Dusky)=c("Resident","Transient","Wide ranging")
E.g._Thickskin=c("SS.22","SS.61","SS.25")
names(E.g._Thickskin)=c("Resident","Transient","Wide ranging")
E.g._BP=list(Dusky=E.g._Dusky,Thickskin=E.g._Thickskin)

for (s in 1:length(Detections.species_BP))             
{ 
  # plot presence absence time line
  tiff(file=paste(hndl.BP,"Time.lines_",SPECIES[s],".tiff",sep=""),width = 2400, height = 2400,
       units = "px", res = 300,compression = "lzw")
  par(mai=c(.05,.5,.4,.01),oma=c(2,1,1,1),las=1,xpd=T,mgp=c(2,.8,0))
  store=function.timeline_BP(DATA=Detections.species[[s]],CEX.AX=CX[s],CEX.TOP=CEX.TOp[s],Example_BP=E.g._BP[[s]])
  axis(side = 1, at = Months.for.plot[[s]],labels = Months.for.plot.label[[s]], 
       tcl = 0.5,cex.axis=.85,padj=-1.5) # minor ticks
  hand=Months.for.plot[[s]]
  COL.t=rgb(.1,.1,.1,alpha=0.1)
  ADD.polys=function(H1,H2)polygon(x=c(H1,H2,H2,H1),y=POLYS[[s]],lwd=2,col=COL.t,border=COL.t)
  is.even <- function(x) x %% 2 == 0 #get even numbers
  Even=1:length(hand)
  Even=Even[is.even(1:length(hand))]
  dev.off()
}                                                      


#Bubble plot map    #incomplete: new data show small bubbles paralel to big bubbles in line 3
if(do.map=="YES")
{
  What.to.show="Individuals"
  
  Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
  xbat=sort(unique(Bathymetry$V1))
  ybat=sort(unique(Bathymetry$V2)) 
  reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
  
  
  fn.bubble.rec=function(ESPECIES1,ESPECIES2,ARRAY,PlotlonG,PlotlatT,
                         scale,LETR,Plus.Lon,Plus.Lat,Plus.Lon1,Plus.Lat1,COLOR,Loc1,Loc2,
                         Leg1,Leg2,Leg3,scale2)
  {
    Det=subset(Detections,!(TagCode.original%in%only.recaptured) & Recapture.hit=="NO")
    
    if(What.to.show=="Individuals")
    {
      Det$Dummy=with(Det,paste(TagCode,Species,Station))
      Det=Det[!duplicated(Det$Dummy),]
    }
    
    DATA=subset(Det,Species==ESPECIES1 & Array==ARRAY  &!(Area.release=="SA"))
    DATA2=subset(Det,Species==ESPECIES2 & Array==ARRAY &!(Area.release=="SA"))
    
    TABLE=table(DATA$Station)
    TABLE2=table(DATA2$Station)
    
    if(!What.to.show=="Individuals")
    {
      N=sum(TABLE)
      N2=sum(TABLE2)
      
    }
    if(What.to.show=="Individuals")
    {
      N=length(unique(DATA$TagCode))
      N2=length(unique(DATA2$TagCode))
    }
    
    MATRIX=matrix(names(TABLE))
    MATRIX=matrix(as.numeric(unlist(strsplit(MATRIX, split=" "))),ncol=2,byrow=T)
    LAT=MATRIX[,1]
    LONG=MATRIX[,2]
    zo=TABLE/sum(TABLE,na.rm=T)
    if(What.to.show=="Individuals") zo=TABLE/N
    
    if(nrow(TABLE2)>0)
    {
      MATRIX2=matrix(names(TABLE2))
      MATRIX2=matrix(as.numeric(unlist(strsplit(MATRIX2, split=" "))),ncol=2,byrow=T)
      LAT2=MATRIX2[,1]
      LONG2=MATRIX2[,2]
      zo2=TABLE2/sum(TABLE2,na.rm=T)
      if(What.to.show=="Individuals") zo2=TABLE2/N2   
    }
    
    
    plotMap(worldLLhigh, xlim=PlotlonG,ylim=PlotlatT,plt = c(.001, 1, 0.0775, 1),
            col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=PlotlatT,xlim=PlotlonG, zlim=c(-1,-100),
            nlevels = 2,labcex=1.25,lty = 1,col=c(COLOR,COLOR,"transparent"),add=T)  
    
    #add receivers
    points(STATIONS$longitude,STATIONS$latitude,col=1,pch=19)
    box(lwd=2)
    
    #add hits
    if(nrow(TABLE)>0)points(LONG-Plus.Lon,LAT*(1-Plus.Lat),pch=21,cex=zo*scale,col=COL.Sp[1],bg=COL.Sp.bg[1],lwd=2)
    if(nrow(TABLE2)>0)points(LONG2+Plus.Lon, LAT2*(1+Plus.Lat),pch=21,cex=zo2*scale,col=COL.Sp[2],bg=COL.Sp.bg[2],lwd=2)
    if(!(ARRAY=="Perth"))legend("topright",LETR,cex=3,bty='n')
    
    
    if(LETR=="2")
    {   
      if(!What.to.show=="Individuals") 
      {
        legend("bottomleft",legend=c(paste("Dusky (n= ",N," detections)",sep=""),
                                     paste("Sandbar (n= ",N2," detections)",sep="")),
               pch=21,col=COL.Sp,pt.bg=COL.Sp.bg ,pt.lwd=2,bty='n',cex=1.75,pt.cex=2)
      }
      if(What.to.show=="Individuals") 
      {
        legend("bottomleft",legend=c(paste("Dusky (n= ",N," sharks)",sep=""),
                                     paste("Sandbar (n= ",N2," sharks)",sep="")),
               pch=21,col=COL.Sp,pt.bg=COL.Sp.bg ,pt.lwd=2,bty='n',cex=1.75,pt.cex=2)
      }
      
    }
    
    if(LETR=="1" & ARRAY=="Ningaloo")
    {
      points(Loc1,Loc2,pch=21,cex=Leg1*scale2)
      points(Loc1,Loc2*1.000183,pch=21,cex=Leg2*scale2)
      points(Loc1,Loc2*1.000366,pch=21,cex=Leg3*scale2)
      boxed.labels(Loc1,Loc2*0.9992558,paste(Leg1*100,"%",sep=""), cex = .75,bg="white",border=NA)
      boxed.labels(Loc1,Loc2*0.9996425,paste(Leg2*100,"%",sep=""), cex = .75,bg="white",border=NA)
      boxed.labels(Loc1,Loc2*1.00001,paste(Leg3*100,"%",sep=""), cex = .75,bg="white",border=NA)
    }
    
  }
  
  
  fn.main.map=function(XLIM,YLIM,LONGSEQ,LATSEQ,Xtext,Ytext,TEXT,TEXT.cx)
  {
    plotMap(worldLLhigh, xlim=XLIM,ylim=YLIM,plt = c(.001, 1, 0.075, 1),
            col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    points(STATIONS$longitude,STATIONS$latitude,col=1,pch=19)
    axis(side = 1, at =LONGSEQ, labels = LONGSEQ, tcl = .5,las=1,cex.axis=1.5,padj=-.5)
    axis(side = 2, at = LATSEQ, labels = -LATSEQ,tcl = .5,las=2,cex.axis=1.5,hadj=.75)
    box(lwd=2)
    mtext("Latitude (ºS)",side=4,line=1,las=3,cex=1.5)
    mtext("Longitude (ºE)",side=1,line=2.5,cex=1.5)
    contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=YLIM,xlim=XLIM, zlim=c(-1,-200),
            nlevels = 4,labcex=1.,lty = 1,col=c(COLOR,COLOR,COLOR,COLOR,"transparent"),add=T)
    text(Xtext,Ytext,TEXT,cex=TEXT.cx)
  }
  
  def.par <- par(no.readonly = TRUE)  #to reset par to default
  COLOR="dark grey"
  COL.Sp=c(1,1)
  COL.Sp.bg=c("grey85","grey60")

  
  #_ Ningaloo only
  plotlong.lines=list(c(113.8,114.04),c(113.54,113.68),c(113.57,113.865))
  plotlat.lines=list(c(-22.00,-21.845),c(-22.64,-22.55),c(-23.16,-22.975))
  tiff(file=paste(hndl.BP,"Figure1.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mar = c(3, 3, 0, 0))
  m <- cbind(c(2, 3,4),c(1, 1,1))
  layout(m)
  #layout.show(4)
  fn.main.map(XLIM=plotlong[[1]],YLIM=c(-23.25,-21.75),LONGSEQ=Long.seq[[1]],LATSEQ=Lat.seq[[1]],
              Xtext=rbind(114,113.7,113.85),Ytext=rbind(-21.9,-22.6,-23.12),TEXT=c(1:3),3)
  
  #Add detections by species
  LETRA=list(as.character(1),as.character(2),as.character(3))
  for (j in 1:length(LETRA))
  {
    fn.bubble.rec(SPECIES[1],SPECIES[2],"Ningaloo",plotlong.lines[[j]],
                  plotlat.lines[[j]],10,LETRA[[j]],0,6e-04,0,8e-04,COLOR,Loc1=113.82,Loc2=-21.875,
                  Leg1=1,Leg2=.75,Leg3=.5,10)
  } 
  
  #Add West Oz
  WestOz=c(112,129,-36,-12)
  FIG=c(.70,.92,.01,.6)
  par(fig=FIG, new = T,xpd=T)
  plotMap(worldLLhigh, xlim=WestOz[1:2],ylim=WestOz[3:4],plt = c(.1, 1, 0.075, 1),
          bg="transparent",
          col="black",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  text(121,-25,("Western"),col="white", cex=2.25)
  text(121,-29,("Australia"),col="white", cex=2.25)
  box(col=COLOR)
  polygon(x=c(113,115,115,113),y=c(rep(-24,2),rep(-21,2)),lwd=3,
          border="white")
  
  dev.off()
  
}



#18.1.3 Kernel density       #incomplete.....
Do.Kernel="NO"
if(Do.Kernel=="YES")
{
  
  #kernel density function    
  fn.krnl=function(SPec)
  {
    Dat=subset(Detections,Species==SPec & Latitude>=Y.range[1] & Latitude<=Y.range[2])
    Tgs=unique(Dat$TagCode)
    n.Tgs=length(Tgs)
    
    Dat$Month1=factor(Dat$Month,levels=1:12)
    Dat$Day1=factor(Dat$Day,levels=1:31)
    a=vector('list',n.Tgs)
    names(a)=Tgs
    
    for(l in 1:n.Tgs)
    {
      q=subset(Dat,TagCode==Tgs[l])
      
      q=table(q$Month1,q$Day1)
      q[q>0]=1    
      q=as.data.frame.matrix(q)
      q$SUM=rowSums(q)
      q$TagCode=Tgs[l]
      q$Month=1:12
      a[[l]]=subset(q,select=c(SUM,TagCode,Month))  
    }
    
    a=do.call(rbind,a)
    a=subset(a,SUM>=Criteria) 
    Dat$Dummy=paste(Dat$TagCode,Dat$Month,sep="")
    Dat.Kernel=subset(Dat,Dummy%in%paste(a$TagCode,a$Month,sep=""),select=c(TagCode,Longitude,Latitude,Month1))
    This.Mn=sort(unique(Dat.Kernel$Month1))
    
    for(x in 1:length(This.Mn))
    {
      ss=subset(Dat.Kernel,Month1==This.Mn[x])
      loc <- ss[, c("Longitude", "Latitude")]
      id <- ss[, "TagCode"]
      ud <- kernelUD(loc, id, h="href",same4all = T)
      image(ud) # contour corresponds to the home ranges estimated at different probability
      # levels (i.e. the contour 90 corresponds to the 90 percent kernel home-range)
      ver <- getverticeshr(ud, 95) #95% kernel
      plot(ver,ylab="",xlab="",xaxt="n",yaxt="n")
      legend("bottomright",Mns[as.numeric(This.Mn[x])],bty='n',cex=2)
      
      ud.mean <- lapply(ud, mean)
      
    }
    return(list(ud.mean=ud.mean))
  }
  Store.Mean.kernel=vector("list",length(Tgt.spec))
  names(Store.Mean.kernel)=Tgt.spec
  for(q in 1:length(Tgt.spec))Store.Mean.kernel[[q]]=fn.krnl(Tgt.spec[q])
  
}


rm(Detections.species_BP)




#19. -- Dusky natal migration analysis ---

setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement/Natal_migration/Paper"))

#export tagged duskies 
write.csv(subset(TAGS,Species2=="Dusky" & Project.rel=="SMN",select=c(Sex2,ReleaseLength,DATE,
                              ReleaseLongitude2,ReleaseLatitude2)),"Dusky.releases.csv")


START.dusky.mig=as.POSIXlt("2011-06-01")

#19.1 time line
GREYSCALE="YES"
if(GREYSCALE=="YES")
{
  North.col="grey25"
  South.col="grey65"
  COL.kld=rgb(.9,.9,.9,alpha=0.1)
  COL.wrm=rgb(.1,.1,.1,alpha=0.1)
}
if(GREYSCALE=="NO")
{
  North.col="brown4"
  South.col="darkslategray4"
  COL.kld=rgb(.1,.1,.9,alpha=0.1)
  COL.wrm=rgb(.9,.1,.1,alpha=0.1)
}

fn.dusky.natal.mig=function(DATA,ARRAY,FL_50,HADJ,PADJ,Show.what,Y1)
{
  DATA=subset(DATA,Array%in%ARRAY & !TagCode.original%in%only.recaptured)
  DATA=subset(DATA,!is.na(FL))
  DATA$Sex=as.character(DATA$Sex)
  n.detected=DATA[!(duplicated(DATA$TagCode)),]
  n.detected=table(n.detected$Sex)
  DATA$Rel.area=with(DATA,ifelse(ReleaseLatitude>= (-26),"North","South"))
  DATA=subset(DATA,select=c(TagCode,FL,Sex,Zone.rel,ReleaseDate,Year.rel,Month.rel,
                            Date.local,Year, Month,Array,Julian.release,Julian,Rel.area))
  DATA$Mov.area=with(DATA,ifelse(Array=="Ningaloo","North","South"))
  DATA=DATA[order(DATA$Sex,DATA$FL,DATA$Date.local),]
  DATA$area=with(DATA,ifelse(Mov.area=="North",1,10))
  
  First.label=START.dusky.mig
  Last.label=max(DATA$Date.local)
  Labels=seq(First.label, Last.label, by = "6 month")
  Ticks=(Labels-First.label)/(24*3600)
  
  DATA$Julian=as.numeric(DATA$Date.local-First.label)
  DATA$Julian.release=as.numeric(DATA$ReleaseDate-First.label)
  
  
  # table of hits
  DATA$Julian=as.factor(DATA$Julian)
  DATA$TagCode=as.factor(DATA$TagCode)
  Table.hits.date=with(DATA,table(TagCode,Julian))
  First.Jul=colnames(Table.hits.date)[1] #this is the first julian date
  
  ORDR=as.data.frame.matrix(Table.hits.date)
  Lenzs=DATA[!duplicated(DATA$TagCode),match(c("TagCode","FL"),names(DATA))]
  ORDR$TagCode=rownames(ORDR)
  ORDR=merge(ORDR,Lenzs,by="TagCode")
  ORDR=ORDR[order(ORDR$FL),]
  ID=match(ORDR$TagCode,rownames(Table.hits.date))
  Table.hits.date=Table.hits.date[ID,]
  
  # create colors by array
  Table.tag.day.area.north=with(subset(DATA,Mov.area=='North'),table(TagCode,Julian))
  Table.tag.day.area.south=with(subset(DATA,Mov.area=='South'),table(TagCode,Julian))
  Table.tag.day.area.north[Table.tag.day.area.north>0]=1
  Table.tag.day.area.south[Table.tag.day.area.south>0]=10
  Table.tag.day.area=Table.tag.day.area.north+Table.tag.day.area.south
  
  ID=match(rownames(Table.hits.date),rownames(Table.tag.day.area))
  Table.tag.day.area=Table.tag.day.area[ID,]
  
  Mat.colors=as.data.frame.matrix(Table.tag.day.area)
  Mat.colors=ifelse(Mat.colors==1,North.col,ifelse(Mat.colors==10,South.col,NA))
  Mat.colors.legends=c("North","South")
  Mat.colors.legends.colors=c(North.col,South.col)
  
  # create pch by array
  Mat.pch=as.data.frame.matrix(Table.tag.day.area)
  Mat.pch=ifelse(Mat.pch==1,15,ifelse(Mat.pch==10,15,NA))
  Mat.colors.legends.pch=c(15,15)
  
  This.col=match(First.Jul,colnames(Table.hits.date))
  x=as.numeric(colnames(Table.hits.date[,This.col:ncol(Table.hits.date)]))
  y=1:nrow(Table.hits.date) 
  n=ncol(Table.hits.date)
  n.shks=nrow(Table.hits.date)
  Rels=DATA[!duplicated(DATA$TagCode),match(c("TagCode","FL","Julian.release","Rel.area"),names(DATA))]
  ID=match(rownames(Table.hits.date),as.character(Rels$TagCode))
  Rels=Rels[ID,]
  Rels$col=with(Rels,ifelse(Rel.area=="North",North.col,South.col))
  plot(x,x,ylim=c(y[1],y[length(y)]),ylab="",xlab="",xlim=c(0,x[length(x)]),yaxt='n',xaxt='n')
  legend("top",Mat.colors.legends,bty='n',pt.cex=2,pch=Mat.colors.legends.pch,
         col=Mat.colors.legends.colors,horiz=T,inset=c(0,-.08),cex=1.75)
  legend("topleft",paste("(FL @ 50% maturity= ",FL_50*100," cm)",sep=""),
         bty='n',inset=c(-0.1,-.06),cex=1.05)
  
  axis(side = 1, at = Ticks,labels = Labels,tcl = 0.5,cex.axis=1,padj=PADJ) 
  mtext("Date",1,line=2,cex=1.75)
  if(Show.what=="FL") mtext("FL (cm)",2,las=3,line=2.05,cex=1.75)
  if(Show.what=="FL")axis(side = 2, at = 1:n.shks,labels =ORDR$FL*100 ,tcl = 0.5,cex.axis=1,hadj=HADJ)
  if(Show.what=="TagCode") axis(side = 2, at = 1:n.shks,labels =ORDR$TagCode ,tcl = 0.5,cex.axis=1,hadj=HADJ)
  
  #Win-Spring polygons  (coolest water months)   
  Solsti=c(as.POSIXlt("2011-06-21"),seq(as.POSIXlt("2011-12-21"),max(DATA$Date.local),by="6 month"))
  First.pol=c(First.label,Solsti) 
  names(Solsti)=month(Solsti)
  names(Solsti)=ifelse(names(Solsti)==6,"Cold","Warm")
  ADD.polys=function(X1,X2,CL)polygon(x=c(X1,X2,X2,X1),y=c(Y1,Y1,max(y)*1.01,max(y)*1.01),lwd=2,col=CL,border="transparent")
  for(p in 1:length(Solsti))
  {
    if(names(Solsti[p])=="Cold") clr=COL.kld
    if(names(Solsti[p])=="Warm") clr=COL.wrm
    x1=as.numeric(Solsti[p]-First.label)
    if(p<length(Solsti))x2=as.numeric(Solsti[p+1]-First.label)
    if(p==length(Solsti))x2=x1+183
    ADD.polys(x1,x2,clr)
  }
  for(i in 1:n) points(rep(x[i],n.shks),1:n.shks,cex=1.5,pch=Mat.pch[,i],col=Mat.colors[,i])
  points(Rels$Julian.release,1:n.shks,cex=1.7,col=Rels$col,pch=21,bg="white")
  
}

tiff(file="Dusky.fem.tiff",width = 2400, height = 2400,
     units = "px", res = 300,compression = "lzw")
par(mai=c(.05,.5,.4,.01),oma=c(3,1,1,1),las=1,xpd=T,mgp=c(2,.1,0))
fn.dusky.natal.mig(DATA=subset(Detections.species[[1]],Sex=="F" & !Project.rel=="South.Australia"),
                   ARRAY=c("Ningaloo","Perth","Southern.lines"),
                   FL_50=FL_0.5[[1]][1],HADJ=1.15,PADJ=0.5,Show.what="FL",Y1=-.5)
dev.off()

tiff(file="Dusky.male.tiff",width = 2400, height = 2400,
     units = "px", res = 300,compression = "lzw")
par(mai=c(.05,.5,.4,.01),oma=c(3,1,1,1),las=1,xpd=T,mgp=c(2,1,0))
fn.dusky.natal.mig(subset(Detections.species[[1]],Sex=="M" & !Project.rel=="South.Australia"),
                   c("Ningaloo","Perth","Southern.lines"),
                   FL_0.5[[1]][2],HADJ=0.5,PADJ=-1,Show.what="FL",Y1=0.2)
dev.off()


tiff(file="Dusky.fem.TagCode.tiff",width = 2400, height = 2400,
     units = "px", res = 300,compression = "lzw")
par(mai=c(.05,.5,.4,.01),oma=c(3,1,1,1),las=1,xpd=T,mgp=c(2,.1,0))
fn.dusky.natal.mig(DATA=subset(Detections.species[[1]],Sex=="F" & !Project.rel=="South.Australia"),
                   ARRAY=c("Ningaloo","Perth","Southern.lines"),
                   FL_50=FL_0.5[[1]][1],HADJ=1.15,PADJ=0.5,Show.what="TagCode",Y1=-.5)
dev.off()

tiff(file="Dusky.male.TagCode.tiff",width = 2400, height = 2400,
     units = "px", res = 300,compression = "lzw")
par(mai=c(.05,.5,.4,.01),oma=c(3,1,1,1),las=1,xpd=T,mgp=c(2,1,0))
fn.dusky.natal.mig(subset(Detections.species[[1]],Sex=="M" & !Project.rel=="South.Australia"),
                   c("Ningaloo","Perth","Southern.lines"),
                   FL_0.5[[1]][2],HADJ=0.5,PADJ=-1,Show.what="TagCode",Y1=0.2)
dev.off()

#Depth by line 
#note: some lines seem to have the depth wrong... (e.g. Southern 3)
fun.corridor=function(D)
{
  D$Line=with(D,
              ifelse(Latitude< (-21.8) & Latitude>(-22.17),"Ningaloo 1",
                     ifelse(Latitude< (-22.52589) & Latitude>-(22.67555),"Ningaloo 2",
                            ifelse(Latitude<(-22.96274) & Latitude>  (-23.22162),"Ningaloo 3",
                                   ifelse(Latitude<(-31) & Latitude>  (-33),"Perth",
                                          ifelse(Latitude< (-34.06557) & Latitude>(-34.51456) & Longitude < 115.29,"Southern 1",
                                                 ifelse( Latitude< (-34.67056) & Latitude> (-35.43074) & Longitude > 116.14 & Longitude <116.8,"Southern 2",
                                                         ifelse(Latitude< (-34.81793) & Latitude>(-35.27299) & Longitude > 118.2450,"Southern 3",NA))))))))
  
  Ln=sort(unique(D$Line))
  
  for(l in 1:length(Ln))
  {
    x=subset(D, Line==Ln[l])
    x$Depth.bin=floor(x$Depth/10)*10
    x$Depth.bin.shk=with(x,paste(Depth.bin,TagCode))
    
    #Detections
    Det=aggregate(hit~Depth.bin,x,sum)
    plot(Det$Depth.bin,Det$hit/sum(Det$hit),ylab="",xlab="")
    if(l ==1) mtext("Number of detections",3)
    
    #Sharks
    Shk=aggregate(hit~Depth.bin,x[!duplicated(x$Depth.bin.shk),],sum)
    plot(Shk$Depth.bin,Shk$hit/sum(Shk$hit),ylab="",xlab="")
    legend("bottomright",Ln[l],bty='n')
    if(l ==1) mtext("Number of sharks",3)
  }
  
}
tiff(file="Depths.tiff",width = 2400, height = 2400,
     units = "px", res = 300,compression = "lzw")
par(mfrow=c(7,2),mai=c(.05,.5,.4,.01),oma=c(3,1,1,1),las=1,xpd=T,mgp=c(2,1,0))
fun.corridor(D=subset(Detections.species$Dusky,Longitude<119 &!is.na(Depth)))
dev.off()

  #19.2 Movement model
  #data
#removed DS.7 which was only detected once
dat.glm=subset(Detections.species$Dusky,Sex%in%c("F","M") &!is.na(FL) & !TagCode=="DS.7" &!Project.rel=="South.Australia",
     select=c(TagCode,Sex,FL,ReleaseDate,Year.rel,Month.rel,ReleaseLatitude,Date.local,Year,
              Month,Array,Depth,Dist.moved.conseq.det,Longitude,Latitude))

dat.glm$dummy=paste(dat.glm$TagCode,dat.glm$Date.local)
dat.glm=dat.glm[!duplicated(dat.glm$dummy),]  #keep only one detection per day per shark

dat.glm$yday.rec=dat.glm$Date.local$yday    #note: different years are combined (only days within a year used)
dat.glm$yday.rel=dat.glm$ReleaseDate$yday
SOLSTICE=as.POSIXlt(c("2011-03-21","2011-06-21","2011-09-21","2011-12-21"))$yday
names(SOLSTICE)=c("Autumn","Winter","Spring","Summer")


#season defined as astronomical solstice (summer, winter) and equinox (autumn, spring)
dat.glm$Sesn.rel=with(dat.glm,ifelse(yday.rel>=SOLSTICE[2] & yday.rel<SOLSTICE[4],"Wi.Sp",           
            ifelse(yday.rel<SOLSTICE[2] | yday.rel>=SOLSTICE[4],"Su.Au.",NA)))

dat.glm$Sesn.rec=with(dat.glm,ifelse(yday.rec>=SOLSTICE[2] & yday.rec<SOLSTICE[4],"Wi.Sp",           
            ifelse(yday.rec<SOLSTICE[2] | yday.rec>=SOLSTICE[4],"Su.Au.",NA)))

dat.glm$Zn.rel=with(dat.glm,ifelse(ReleaseLatitude>=(-26),"North","South"))
dat.glm$Zn.rec=with(dat.glm,ifelse(Array=="Ningaloo","North","South"))
Tags.d=unique(dat.glm$TagCode)

dat.glm$N=with(dat.glm,ifelse(Zn.rec=="North",1,0))
dat.glm$week=week(dat.glm$Date.local)

  #19.2.1 Plot migrating animals to determine seasons
dat=subset(dat.glm,select=c(TagCode,Sesn.rec,Zn.rec,FL,Sex,Date.local,Month,ReleaseDate,
                            Zn.rel,Dist.moved.conseq.det,N,week))


#Figure 2. Natal migration paper
what=dat
Seasonal=what
use.release="NO"
if(use.release=="YES")
{
  Dummy=with(what,table(TagCode,Zn.rel))
  This.north=names(Dummy[Dummy[,1]>0,1])
  Dummy=with(what,table(TagCode,Zn.rec))
  This.north=unique(c(names(Dummy[Dummy[,1]>0,1]),This.north))
}

#remove release zone to become independent of sampling effort
if(use.release=="NO")
{
  Dummy=with(what,table(TagCode,Zn.rec))
  This.north=names(Dummy[Dummy[,1]>0,1])
}


# Step2. determine seasonal parameters
Seasonal=subset(Seasonal,TagCode%in%This.north)  
Seasonal$TagCode=as.character(Seasonal$TagCode)
TG=unique(Seasonal$TagCode)
kEP=rep(NA,length(TG))
kep="round.trip"

#round-trip migrations (N-S-N or S-N-S)    
if(kep=="round.trip")
{
  for(s in 1:length(TG))
  {
    a=subset(Seasonal,TagCode==TG[s])
    a=a[order(a$Date.local,a$Zn.rec),]
    #a=a[order(a$Month,a$Zn.rec),]
    a$cum=0
    if(nrow(a)>1)
    {
      if(!a$Zn.rel[1]==a$Zn.rec[1]) a$cum[1]=1
      for(w in 2:nrow(a)) if(!a$Zn.rec[w]==a$Zn.rec[w-1]) a$cum[w]=1
    }
      
    x=NA
    if(sum(a$cum)>1) x=TG[s]
    kEP[s]=x
  }
  kEP=subset(kEP,!is.na(kEP))
}

#North and South occurrence
if(kep=="n.s")
{
  for(s in 1:length(TG))
  {
    a=subset(Seasonal,TagCode==TG[s])
    UNIK=unique(a$Zn.rec)
    x=NA
    if(length(UNIK)>1) x=TG[s]
    kEP[s]=x
  }
  kEP=subset(kEP,!is.na(kEP))
}

Seasonal=subset(Seasonal,TagCode%in%kEP)
Seasonal$S=with(Seasonal,ifelse(Zn.rec=="South",1,0))
Seasonal$day=Seasonal$Date.local$yday   #days from beginning of year


N.day=aggregate(cbind(N,S)~day,Seasonal,sum)
PT.mx=max(N.day[,2:3])
N.day.one=N.day
N.day.one[N.day.one>0]=1
N.day.one$day=N.day$day
N.day$P.north=N.day$N/(N.day$N+N.day$S)
plot(N.day$P.north)

By.week=Seasonal
By.week$dummy=with(By.week,paste(TagCode,week,Zn.rec))
By.week=By.week[!duplicated(By.week$dummy),]
a=Seasonal
a$year=year(a$Date.local)
N.week=aggregate(cbind(N,S)~week,a,sum) 
N.week.one=N.week
N.week.one[N.week.one>0]=1
N.week.one$week=N.week$week
N.week$P.north=N.week$N/(N.week$N+N.week$S)
plot(N.week$P.north)


By.month=Seasonal
By.month$dummy=with(By.month,paste(TagCode,Month,Zn.rec))
By.month=By.month[!duplicated(By.month$dummy),]
a=Seasonal
a$year=year(a$Date.local)
N.month=aggregate(cbind(N,S)~Month,a,sum) 
N.month.one=N.month
N.month.one[N.month.one>0]=1
N.month.one$Month=N.month$Month
N.month$P.north=N.month$N/(N.month$N+N.month$S)
plot(N.month$P.north)


#Temperature
Temp$date=paste(Temp$Year,"-",Temp$Month,"-","15",sep="")
Temp$date=as.POSIXlt(as.character(Temp$date))
Temp$week=week(Temp$date)
#Temp$Temperature=Temp$value
#Temp1=subset(Temp,Year%in%2011:2014)
Temp1=subset(Temp,Year%in%2011:2012)
Temp.N=subset(Temp1,Lat <=(-21) & Lat>= (-23) & Long>= 112.5 & Long <114 )
Temp.S=subset(Temp1,Lat <=(-32) & Lat>= (-35) & Long>= 113 & Long <119 )
Mean.temp.N=aggregate(Temperature~Month,Temp.N,mean)
Mean.temp.S=aggregate(Temperature~Month,Temp.S,mean)


#Day length
Margs=c(33.69,115.07)
Perth=c(31.95,115.86)
S.pos=c(mean(c(Margs[1],Perth[1])),mean(c(Margs[2],Perth[2])))
N.pos=c(21.93,114.13)

t0 <- as.POSIXct("2013-01-01 12:00:00", tz="UTC")
t <- seq.POSIXt(t0, by="1 week", length.out=53)


daylength <- function(lon, lat)
{
  t <- as.numeric(t)
  alt <- function(t) sunAngle(t, longitude=lon, latitude=lat)$altitude
  rise=set=rep(NA,length(t))
  for(i in 1:length(t))
  {
    rise[i]=uniroot(alt, lower=t[i]-86400/2, upper=t[i])$root
    set[i] <- uniroot(alt, lower=t[i], upper=t[i]+86400/2)$root
  }
  daylen=set - rise
  return(daylen/3600)
}


South.day.len=daylength(S.pos[2],S.pos[1])
North.day.len=daylength(N.pos[2],N.pos[1])

#solstice                 
AT=week(as.POSIXlt(c("2011-03-21","2011-06-21","2011-09-21","2011-12-21")))
SEaSon=c("Autum","Winter","Spring","Summer")

WEE=1:53
YLIM=c(.05,.3)

  #no map approach
tiff(file="Figure5.tiff",
     width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mai=c(.725,1,.01,1.5),oma=c(1,1,.1,.1),mgp=c(1,.6,0))
plot(WEE,WEE,col="transparent",cex=N.week$N, ylab="",ylim=YLIM,xaxt='n',yaxt='n',xlab="")
at=c(.25,.1)
symbols(N.week.one$week,rep(at[1],nrow(N.week.one)),
        rectangles=matrix(c(rep(1.25,nrow(N.week.one)),N.week$N),byrow=F,ncol=2),add=T,fg="red4",bg="red")
symbols(N.week.one$week,rep(at[2],nrow(N.week.one)),
        rectangles=matrix(c(rep(.75,nrow(N.week.one)),N.week$S),byrow=F,ncol=2),add=T,fg="blue",bg="deepskyblue3")
axis(2,at,c("North","South"),las=2,cex.axis=1.75)
axis(2,at-.0125,c("(22-23° S)","(32-35° S)"),las=2,cex.axis=1,col='transparent')
axis(1,WEE,F,las=2,tcl = -0.25)
axis(1,AT,AT,tcl = -0.5,cex.axis=1.75)
axis(1,AT,paste("(",SEaSon,")",sep=""),tcl = -0.5,cex.axis=1.5,line=1,col="transparent")
mtext("Location",2,las=3,outer=T,cex=2.5,line=-1)
mtext("Week of the year",1,outer=T,cex=2.5,line=-0.25)

X1=AT[2];X2=AT[4];Y1=YLIM[1];Y2=YLIM[2]
polygon(x=c(X1,X2,X2,X1),y=c(Y1,Y1,Y2,Y2),lwd=2,col=rgb(.1,.9,.1,alpha=0.1),border=1)


#add day length
par(new=T)
plot(WEE,South.day.len,axes=F,type='l',lwd=2,col=4,ann=F)
lines(WEE,North.day.len,lwd=3,col=2)
Rng=range(South.day.len)
att=seq(floor(Rng[1]),ceiling(Rng[2]),by=1)
axis(4,att,att,las=1,cex.axis=1.25)
n=2
mtext("Day length (hours)",4,line=n,cex=1.5)

#add temperature
par(new=T)
Rng=range(c(Mean.temp.S$Temperature,Mean.temp.N$Temperature))
YLIM=c(Rng[1],Rng[2])
plot(Mean.temp.S$Temperature,axes=F,type='l',lwd=2,col=4,ann=F,lty=2,ylim=YLIM)
lines(Mean.temp.N$Temperature,lwd=2,col=2,lty=2)

att=seq(floor(Rng[1]),ceiling(Rng[2]),by=1)
axis(4,att,att,las=1,line=n+2,cex.axis=1.25)
mtext("Average temperature (°C)",4,line=n+4,cex=1.5)

dev.off()

#Mapping approach
Do.map.natal.mig="YES"
COLOR="dark grey"

if(!exists("xbat"))
{
  Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
  xbat=sort(unique(Bathymetry$V1))
  ybat=sort(unique(Bathymetry$V2)) 
}
if(!exists("reshaped"))reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))


if(do.map=="NO" & Do.map.natal.mig=="YES") 
{
  OZ.lat=c(-44.5,-11)
  OZ.long=c(112,155)
  XX=c(112,119)
  YY=c(-35.5,-21.5)
  
  fn.axs=function()
  {
    axis(1,WEE,F,las=2,tcl = -0.25)
    axis(1,AT,F,tcl = -0.75)
  }
  fn.week.axs=function() 
  {
    axis(1,AT,AT,tcl = -0.75,cex.axis=1.25)
    #axis(1,AT,paste("(",SEaSon,")",sep=""),cex.axis=1.1,line=.75,col="transparent")
    axis(1,AT,SEaSon,cex.axis=1.1,line=.75,col="transparent")
  }
  fn.scale.bar=function(p1,p2,ADJ1,ADJ2,ADJ3,ADK4)
  {
    KM=distCosine(p1,p2)/1000
    segments(p1[1],p1[2],p2[1],p2[2],lwd=3)
    text(p1[1]*ADJ1,p1[2]*ADJ2,floor(KM),cex=1.5)
    text(p1[1]*ADJ3,p2[2]*ADK4,"kilometres",cex=1.5)
  }
  
  fn.det=function(WW,ZZ,CL1,TxT,LN,BY)
  {  
    WW$N=ifelse(WW$N==0,NA,WW$N)
    WW$S=ifelse(WW$S==0,NA,WW$S)
    z=match(ZZ,names(WW))
    YLIM=c(0,max(c(WW$N),c(WW$S),na.rm=T))
    plot(WW$week,WW[,z],type='h',ylab="",ylim=YLIM,xaxt='n',yaxt='n',xlab="",axes=F,col=CL1,lwd=5)
    fn.axs()
    mtext("Number of sharks",2,las=3,cex=1.3,line=LN)
    legend('topright',TxT,cex=2,bty='n',text.col=CL1)  
    att=seq(floor(YLIM[1]),ceiling(YLIM[2]),by=BY)
    axis(2,att,att,las=1,line=0,cex.axis=1.25)
    # X1=AT[2];X2=AT[4];Y1=YLIM[1];Y2=YLIM[2]*1.02
    # polygon(x=c(X1,X2,X2,X1),y=c(Y1,Y1,Y2,Y2),lwd=2,col=rgb(.1,.9,.1,alpha=0.1),border=1)
    box()
  }
  
  fn.T.day.len=function(Y1,Y2,TxT,BY,LN)
  {
    Rng=range(c(Y1,Y2))
    plot(Y1,type='l',lwd=3,col=4,ann=F,ylim=c(Rng[1],Rng[2]*1.03),axes=F)
    lines(Y2,lwd=3,col=2)
    att=seq(floor(Rng[1]),ceiling(Rng[2]),by=BY)
    axis(2,att,att,las=1,line=0,cex.axis=1.25)
    mtext(TxT,2,line=LN,cex=1.3)
  }
  
  fn.map.natal.mig=function(XLIM,YLIM,LONGSEQ,LATSEQ,EXM,PER,FIG,txt.Oz,cx.Oz)
  {
    plotMap(worldLLhigh, xlim=XLIM,ylim=YLIM,plt = c(.001, 1, 0.05, 1),
            col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    points(STATIONS$longitude,STATIONS$latitude,col=3,pch=19,cex=1.5)
    contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=YLIM,xlim=XLIM, zlim=c(-1,-200),
            nlevels = 4,labcex=1.,lty = 1,col=c(COLOR,COLOR,COLOR,COLOR,"transparent"),add=T)
    
    axis(side = 1, at =LONGSEQ, labels = LONGSEQ, tcl = .5,las=1,cex.axis=1.5,padj=-.5)
    axis(side = 4, at = LATSEQ, labels = -LATSEQ,tcl = .5,las=2,cex.axis=1.5,hadj=.25)
    box(lwd=1.5)
    mtext("Latitude (ºS)",side=4,line=2.25,las=3,cex=1.75)
    mtext("Longitude (ºE)",side=1,line=2,cex=1.75)
    
    #North and south polygons
    polygon(x=c(113.25,114.25,114.25,113.25),y=c(-23.25,-23.25,-21.75,-21.75),
            lwd=2,col=rgb(.9,.01,.01,alpha=.3),border=2)
    polygon(x=c(114.5,118.75,118.75,114.5),y=c(-35.35,-35.35,-31.5,-31.5),
            lwd=2,col=rgb(.01,.01,.09,alpha=.3),border=4)
    text(EXM[1],EXM[2],"Ningaloo",cex=2)
    text(PER[1],PER[2],"Perth",cex=2)
    
    fn.scale.bar(c(116,-25),c(117,-25),1.004,.9875,1.005,1.01)
    
    #Add Oz inset
    par(fig=FIG, new = T,mgp=c(.1,.4,0))
    plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
            col="grey60",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    box()
    polygon(x=c(XLIM,rev(XLIM)),y=c(YLIM[1],YLIM[1],YLIM[2],YLIM[2]),lwd=2,
            col="transparent")
            # col=rgb(.01,.9,.01,alpha=.3))
    text(txt.Oz[1],txt.Oz[2],("Australia"),col="black", cex=cx.Oz)
  }
  
  tiff(file="Figure_5.Natal_migration_paper.tiff",
       width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mar = c(3.2, 4, 0.95, 0),oma=c(1,.1,.1,2),mgp=c(1,.75,0))
  m <- cbind(1:2,rep(3,2))
  #m <- cbind(1:4,rep(5,4))
  layout(m)
  
  #Add detections north
  fn.det(WW=N.week,ZZ="N",CL1=2,"North",2.5,5)
  fn.week.axs()
  
  # #add day length
  # fn.T.day.len(South.day.len,North.day.len,"Day length (hours)",1,2.5)
  # fn.axs() 
  # fn.week.axs()
  # 
  # #Add temperature
  # fn.T.day.len(Mean.temp.S$Temperature,Mean.temp.N$Temperature,"Temperature (°C)",2,2.5)
  # axis(1,Mean.temp.S$Month,F,las=2,tcl = -0.25)
  # axis(1,c(1,3,6,9,12),c("Jan","Mar","Jun","Sep","Dec"),tcl = -0.75,cex.axis=1.5)
  # mtext("Month",1,outer=F,cex=1.5,line=2)
  # 
  
  #Add detections south
  fn.det(WW=N.week,ZZ="S",CL1=4,"South",2.5,5)
  fn.week.axs()
  mtext("Week of the year",1,outer=F,cex=1.3,line=3)
  par(xpd=NA)
  X1=AT[2];X2=AT[4];Y1=-.6;Y2=64
  #polygon(x=c(X1,X2,X2,X1),y=c(Y1,Y1,Y2,Y2),col=rgb(.2,.7,.2,alpha=0.1),border="transparent")
  
  
  #Add WA
  fn.map.natal.mig(XLIM=XX,YLIM=YY,LONGSEQ=XX[1]:XX[2],LATSEQ=seq(-35,-21,by=1),
                   EXM=c(114.75,-22.5),PER=c(116.5,-32),FIG=c(.75,.92,.35,.75),
                   txt.Oz=c(135,-25),cx.Oz=1.25)
  dev.off()

  
  
  #Week barplot  
  tiff(file="Figure_5.barplot.tiff",
       width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mar = c(3, 3, 0.95, 0),oma=c(1,.1,.1,2),mgp=c(1,.75,0),las=1)
  BP=barplot(rbind(N.week$N,N.week$S),col=c(North.col,South.col))
  axis(1,BP,F,tcl = -0.5,cex.axis=1.25)
  axis(1,BP[AT],AT,tcl = -0.75,cex.axis=1.25)
  axis(1,BP[AT],SEaSon,cex.axis=1.1,line=.85,col="transparent")
  mtext("Week of the year",1,outer=F,cex=1.5,line=2.75)
  box()
  mtext("Number of sharks",2,las=3,line=1.75,cex=1.75)
  legend("topright",c("North","South"),bty='n',fill=c(North.col,South.col),cex=1.75)
  dev.off()
  
  #barplot(rbind(N.month$N,N.month$S))  
}



#Show day length and T by year, north and south
Show.dy.ln="NO"
if(Show.dy.ln=="YES")
{
  #day length
  South.day.len.years=list(o11=NA,o12=NA,o13=NA,o14=NA)
  North.day.len.years=South.day.len.years
  YrSs=2011:2014
  
  for( i in 1:length(YrSs))
  {
    t0 <- as.POSIXct(paste(YrSs[i],"-01-01 12:00:00",sep=""), tz="UTC")
    t <- seq.POSIXt(t0, by="1 week", length.out=53)
    
    South.day.len.years[[i]]=daylength(S.pos[2],S.pos[1])
    North.day.len.years[[i]]=daylength(N.pos[2],N.pos[1])
  }
  Temp1.N=subset(Temp,Lat <=(-21) & Lat>= (-23) & Long>= 112.5 & Long <114 )   
  Temp1.S=subset(Temp,Lat <=(-32) & Lat>= (-35) & Long>= 113 & Long <119 )
  Mean.temp1.N=aggregate(Temperature~Month,Temp1.N,mean)
  Mean.temp1.S=aggregate(Temperature~Month,Temp1.S,mean)
  LTYPE=1:4
  tiff(file="Day length and Temp by year.tiff",
       width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(2,1),mar=c(2,1.75,2,0.1),oma=c(1,2.75,.1,0.1),mgp=c(2, 0.65, 0),las=1)
  plot(WEE,ylim=c(min(c(unlist(South.day.len.years),unlist(North.day.len.years))),
                  max(c(unlist(South.day.len.years),unlist(North.day.len.years)))),
       col='transparent',ylab='',xlab='',xaxt='n',cex.axis=1.25)
  for( i in 1:length(YrSs))
  {
    lines(WEE,South.day.len.years[[i]],lwd=1,col=4,lty=LTYPE[i])
    lines(WEE,North.day.len.years[[i]],lwd=1,col=2,lty=LTYPE[i])
  }
  legend("bottomright",paste(YrSs),lty=LTYPE,cex=1.25,bty='n')
  mtext("Day length (hours)",2,line=2.75,cex=1.75,las=3)
  fn.axs() 
  mtext("Week of the year",1,outer=F,cex=1.75,line=2.65)
  axis(1,AT,AT,tcl = -0.75,cex.axis=1.5)
  axis(1,AT,paste("(",SEaSon,")",sep=""),cex.axis=1.25,line=.75,col="transparent")
  
  #add temperature
  plot(1:12,ylim=c(min(c(Mean.temp1.N$Temperature,Mean.temp1.S$Temperature)),
                   max(c(Mean.temp1.N$Temperature,Mean.temp1.S$Temperature))),
       col='transparent',ylab='',xlab='',xaxt='n',cex.axis=1.25)
  for( i in 1:length(YrSs))
  {
    aaa=subset(Mean.temp1.S,Year==YrSs[i])
    lines(aaa$Month,aaa$Temperature,lwd=1,col=4,lty=LTYPE[i])
    
    aaa=subset(Mean.temp1.N,Year==YrSs[i])
    lines(aaa$Month,aaa$Temperature,lwd=1,col=2,lty=LTYPE[i])
    
  }
  mtext("Average temperature (°C)",2,line=2.75,cex=1.75,las=3)
  axis(1,c(1,3,6,9,12),c("Jan","Mar","Jun","Sep","Dec"),tcl = -0.75,cex.axis=1.5)
  mtext("Month",1,outer=F,cex=1.75,line=2)
  dev.off()
  
}


#Symbols approach
# fn.det=function(WW,ZZ,CL,CL1,YLIM,N,TxT,LN)
# {  
#   plot(WEE,WEE,col="transparent",ylab="",ylim=YLIM,xaxt='n',yaxt='n',xlab="",axes=F)
#   symbols(WW$week,rep(.5,nrow(WW)),
#           rectangles=matrix(c(rep(N,nrow(WW)),ZZ),byrow=F,ncol=2),add=T,fg=CL,bg=CL1)
#   fn.axs()
#   mtext(TxT,2,las=3,cex=1.35,line=LN)  
#   RNG=range(ZZ)
#   att=seq(floor(RNG[1]),ceiling(RNG[2]),by=2)
#   axis(2,att,att,las=1,line=0,cex.axis=1.25)
#   X1=AT[2];X2=AT[4];Y1=YLIM[1];Y2=YLIM[2]
#   polygon(x=c(X1,X2,X2,X1),y=c(Y1,Y1,Y2,Y2),lwd=2,col=rgb(.1,.9,.1,alpha=0.1),border=1)
# }
#fn.det(WW=N.week.one,ZZ=N.week$N,CL="red4",CL1="red",YLIM=c(0.44,0.56),N=1.25,"North",2.5)
#fn.det(WW=N.week.one,ZZ=N.week$S,CL="blue",CL1="deepskyblue3",YLIM=c(0.44,0.56),N=.75,"South",2.5)


#Paper data summary 
fn.summary.paper.natal=function(TaG,DaT)
{
  N.tagged=nrow(TaG)
  Sx=table(TaG$Sex2)
  Male.Female.ratio.tagged=paste(1,":",round(Sx[1]/Sx[2],1))
  Size.range.tagged.FL=range(TaG$Size,na.rm=T)
  N.fem.t=Sx[1]
  N.mal.t=Sx[2]
  N.hits=nrow(DaT)
  First.rel=as.character(unique(min(DaT$ReleaseDate)))
  Last.rel=as.character(unique(max(DaT$ReleaseDate)))
  First.det=as.character(unique(min(DaT$Date.local)))
  Last.det=as.character(unique(max(DaT$Date.local)))
    
  N.detected=length(unique(DaT$TagCode))
  
  DaT$dumy=with(DaT,paste(TagCode,Date.local))
  dum=DaT[!duplicated(DaT$dumy),]
  dum1=DaT[!duplicated(DaT$TagCode),]
  Sx=table(as.character(dum1$Sex))
  N.fem.d=Sx[1]
  N.mal.d=Sx[2]
  
  
  Male.Female.ratio.detected=paste(1,":",round(Sx[1]/Sx[2],1))
  
  Size.range.detected=range(dum1$FL*100,na.rm=T)
  
  dum$Zn.rec=with(dum,ifelse(Array=="Ningaloo","North","South"))
  TAb=table(dum$Zn.rec,dum$TagCode)
  TAb[TAb>0]=1
  AA=colSums(TAb)
  
  #Days detected
  UNI=unique(DaT$TagCode)
  Dusky.days.data=subset(Dist.days.monitored,Species=="Dusky" & TagCode%in%UNI)
  Dusky.days.detected=subset(Dusky.days.data,select=c(TagCode,days.det))
  MIN.days.det=min(Dusky.days.detected$days.det)
  MAX.days.det=max(Dusky.days.detected$days.det)
  MEDIAN.days.det=median(Dusky.days.detected$days.det)
  SD.days.det=sd(Dusky.days.detected$days.det)
  COMBINED.days.det=sum(Dusky.days.detected$days.det)
  
  #migration
  Number.detected.N.and.S=length(AA[AA>1])
  Number.migrating=length(unique(kEP)) 
  
  
  RET=c(N.tagged,Male.Female.ratio.tagged,Size.range.tagged.FL,N.detected,
        Male.Female.ratio.detected,Size.range.detected,
        MAX.days.det,MEDIAN.days.det,SD.days.det,COMBINED.days.det,
        Number.detected.N.and.S,Number.migrating,
        N.fem.t,N.mal.t,N.hits,N.fem.d,N.mal.d,
        First.rel,Last.rel,First.det,Last.det)
  names(RET)=c("n.tag","M:F.tag","FL.min.tag","FL.max.tag","n.det",
               "M:F.det","FL.min.det","FL.max.det",
               "max days det","median days det","SD.days.det","total days det all shks",
               "n.det.N&S","n.mig",
               "N.fem.t","N.mal.t","N.hits","N.fem.d","N.mal.d",
               "First.rel","Last.rel","First.det","Last.det")
  
  return(RET)
}

a=fn.summary.paper.natal(subset(TAGS,Project.rel=="SMN" & Species2=="Dusky"),
                         subset(Detections,Species=="Dusky" & !Project.rel=="South.Australia"))
#a=Table1[,match(c("Variable","Dusky"),colnames(Table1))]
write.csv(a,"Summary.data.csv",row.names=T)


#Select migrating duskies (round-trip migrants, used in Seasonal)
Migrating.Dusky=subset(Detections.species$Dusky,TagCode%in%unique(Seasonal$TagCode),
                       select=c(TagCode,Sex,ReleaseDate,ReleaseLatitude,ReleaseLongitude,
                                Array,Date.local,Time.local,Latitude,Longitude,Dist.moved.conseq.det,Dup))


Migrating.Dusky$Zn.rel=with(Migrating.Dusky,ifelse(ReleaseLatitude>=(-26),"North","South"))
Migrating.Dusky$Zn.rec=with(Migrating.Dusky,ifelse(Array=="Ningaloo","North","South"))

#Define receiver lines
Migrating.Dusky$Line=with(Migrating.Dusky,ifelse(Latitude>(-22.12),"Ning.1",
             ifelse(Latitude<=(-22.12)& Latitude>(-22.79),"Ning.2", 
            ifelse(Latitude<=(-22.79)& Latitude>(-23.298),"Ning.3", 
             ifelse(Latitude<=(-31.45)& Latitude>(-32.67),"Perth", 
              ifelse(Latitude<=(-32.67)& Latitude>(-34.765) & Longitude <116,"South.1", 
             ifelse(Latitude<=(-34.765)& Latitude>(-35.691) & Longitude >116 & Longitude <117,"South.2",
          ifelse(Latitude<=(-34.69)& Latitude>(-35.47) & Longitude >117,"South.3",NA))))))))

#Select hits at different lines plus the release
Migrating.Dusky=Migrating.Dusky[order(Migrating.Dusky$TagCode,Migrating.Dusky$Date.local),]
N.mig=nrow(Migrating.Dusky)

Migrating.Dusky$TagCode.prev=c(NA,Migrating.Dusky$TagCode[1:(N.mig-1)])

Migrating.Dusky$Zn.rec.prev=c(NA,Migrating.Dusky$Zn.rec[1:(N.mig-1)])
Migrating.Dusky$Zn.rec.prev=with(Migrating.Dusky,
                                 ifelse(TagCode==TagCode.prev,Zn.rec.prev,NA))
Migrating.Dusky$Latitude.prev=c(NA,Migrating.Dusky$Latitude[1:(N.mig-1)])
Migrating.Dusky$Latitude.prev=with(Migrating.Dusky,
                                   ifelse(TagCode==TagCode.prev,Latitude.prev,NA))
Migrating.Dusky$Longitude.prev=c(NA,Migrating.Dusky$Longitude[1:(N.mig-1)])
Migrating.Dusky$Longitude.prev=with(Migrating.Dusky,
                                    ifelse(TagCode==TagCode.prev,Longitude.prev,NA))

Migrating.Dusky$Array.prev=c(NA,Migrating.Dusky$Array[1:(N.mig-1)])
Migrating.Dusky$Array.prev=with(Migrating.Dusky,
                                ifelse(TagCode==TagCode.prev,Array.prev,NA))


Migrating.Dusky$Line.prev=c(NA,Migrating.Dusky$Line[1:(N.mig-1)])
Migrating.Dusky$Line.prev=with(Migrating.Dusky,
                               ifelse(TagCode==TagCode.prev,Line.prev,NA))

Migrating.Dusky$Latitude.prev=with(Migrating.Dusky,ifelse(Dup==FALSE,ReleaseLatitude,Latitude.prev))
Migrating.Dusky$Longitude.prev=with(Migrating.Dusky,ifelse(Dup==FALSE,ReleaseLongitude,Longitude.prev))


Migrating.Dusky$KEEP=with(Migrating.Dusky,ifelse(Dup==FALSE,"Keep",
                     ifelse(!is.na(Line.prev)& !Line.prev==Line,"Keep",
                    "Do-not")))
Migrating.Dusky=subset(Migrating.Dusky,KEEP=="Keep")


TG=unique(Migrating.Dusky$TagCode)
Round.trip=rep(NA,length(TG))
names(Round.trip)=TG

for(s in 1:length(TG))
{
  a=subset(Migrating.Dusky,TagCode==TG[s])
  a=a[order(a$Date.local,a$Zn.rec),]
  a$cum=0
  if(nrow(a)>1)
  {
    if(!a$Zn.rel[1]==a$Zn.rec[1]) a$cum[1]=1
    for(w in 2:nrow(a)) if(!a$Zn.rec[w]==a$Zn.rec[w-1]) a$cum[w]=1
  }
  x=sum(a$cum)
  Round.trip[s]=x/2
}

#Distance travelled between arrays
Migrating.Dusky.dist.trav=subset(Migrating.Dusky,select=c(TagCode,TagCode.prev,Array,Date.local,
            Dist.moved.conseq.det,Dup,Zn.rel,Zn.rec,Zn.rec.prev,ReleaseDate,ReleaseLongitude,
            ReleaseLatitude,Latitude,Longitude))

Migrating.Dusky.dist.trav$keep=with(Migrating.Dusky.dist.trav,
      ifelse(Dup==FALSE,"Keep",ifelse(TagCode.prev==TagCode & !Zn.rec==Zn.rec.prev,"Keep","do_not")))
Migrating.Dusky.dist.trav=subset(Migrating.Dusky.dist.trav,keep=="Keep")
Migrating.Dusky.dist.trav=Migrating.Dusky.dist.trav[order(Migrating.Dusky.dist.trav$TagCode,
              Migrating.Dusky.dist.trav$Date.local),]
Migrating.Dusky.dist.trav$TagCode.prev[1]=Migrating.Dusky.dist.trav$TagCode[1]

dummy=vector('list',length(TG))

for(s in 1:length(TG))
{
  a=subset(Migrating.Dusky.dist.trav,TagCode==TG[s])
  nn=1:nrow(a)-1
  a$Stage=with(a,ifelse(Dup==FALSE,NA,nn))
  dummy[[s]]=a
}
Migrating.Dusky.dist.trav=do.call(rbind,dummy)

#add release to account for recapture in different zone
dum=vector('list',length(TG))
for(i in 1:length(TG))
{
  a=subset(Migrating.Dusky.dist.trav,TagCode==TG[i])
  a$Added.Rel='NO'
  if(a$Zn.rel[1]==a$Zn.rec[1]) dum[[i]]=a
  
  if(!a$Zn.rel[1]==a$Zn.rec[1])
  {
    a$Zn.rec.prev[1]=a$Zn.rel[1]
    nn=1:sum(!is.na(a$Zn.rec.prev))
    a$Stage=nn
    a$Added.Rel="YES"
    dum[[i]]=a
  }
}

Migrating.Dusky.dist.trav=do.call(rbind,dum)


Cum.dist.Tra.mig.dus=aggregate(Dist.moved.conseq.det~TagCode,Migrating.Dusky.dist.trav,sum)
Round.trip=data.frame(Round.trip)
Round.trip$TagCode=rownames(Round.trip)
Cum.dist.Tra.mig.dus=merge(Cum.dist.Tra.mig.dus,Round.trip,by="TagCode")
Cum.dist.Tra.mig.dus$Mean.dist.travl=with(Cum.dist.Tra.mig.dus,Dist.moved.conseq.det/Round.trip)

qq=table(Cum.dist.Tra.mig.dus$Round.trip)
Migration.summary=data.frame(
  median.km.travl.per.round.mig=median(Cum.dist.Tra.mig.dus$Mean.dist.travl),
  sd.km.travl.per.round.mig=sd(Cum.dist.Tra.mig.dus$Mean.dist.travl),
  max.km.travl.in.total=max(Cum.dist.Tra.mig.dus$Dist.moved.conseq.det))
qq.t=t(as.matrix(qq))
colnames(qq.t)=paste(colnames(qq.t),"round trip (# sharks)")
Migration.summary=cbind(Migration.summary,qq.t)
write.csv(Migration.summary,"Migration.summary.csv",row.names=F)

Migrate.Tags=unique(Migrating.Dusky.dist.trav$TagCode)
Migrate.Tags=subset(Detections,Species=="Dusky" & TagCode%in%Migrate.Tags,select=c(TagCode,Sex))
Migrate.Tags=Migrate.Tags[!duplicated(Migrate.Tags$TagCode),]
Migrate.Tags=merge(Migrate.Tags,Cum.dist.Tra.mig.dus,by="TagCode")
Migrate.Tags=with(Migrate.Tags,table(Round.trip,as.character(Sex)))
write.csv(as.data.frame(Migrate.Tags),"Migration.summary.Tags.sex.csv",row.names=F)


#Plot distance and time between N-S movements for migrants
Nx=nrow(Migrating.Dusky.dist.trav)
dummy=as.character(Migrating.Dusky.dist.trav$Date.local)
Migrating.Dusky.dist.trav$Date.local.prev=c(NA,dummy[1:(Nx-1)])
Migrating.Dusky.dist.trav$Date.local.prev=with(Migrating.Dusky.dist.trav,
                                               ifelse(!is.na(Stage),Date.local.prev,NA))
Migrating.Dusky.dist.trav$Date.local.prev=with(Migrating.Dusky.dist.trav,
                                               ifelse(Added.Rel=="YES" & Stage==1,as.character(ReleaseDate),Date.local.prev))
Migrating.Dusky.dist.trav$Date.local.prev=as.POSIXlt(Migrating.Dusky.dist.trav$Date.local.prev)
Migrating.Dusky.dist.trav$Days.travl=with(Migrating.Dusky.dist.trav,Date.local-Date.local.prev)  

Migrating.Dusky.dist.trav=subset(Migrating.Dusky.dist.trav,!is.na(Days.travl))

Sum.dist=aggregate(Dist.moved.conseq.det~TagCode+Stage+Zn.rec.prev+Zn.rec,Migrating.Dusky.dist.trav,sum)
Sum.days=aggregate(Days.travl~TagCode+Stage+Zn.rec.prev+Zn.rec,Migrating.Dusky.dist.trav,sum)
Sum.dist=Sum.dist[order(Sum.dist$TagCode,Sum.dist$Stage),]
Sum.days=Sum.days[order(Sum.days$TagCode,Sum.days$Stage),]

Sum.dist$Dum=factor(with(Sum.dist,paste(Zn.rec.prev,Zn.rec)))
Sum.days$Dum=factor(with(Sum.days,paste(Zn.rec.prev,Zn.rec)))
Sum.days$Days.travl=as.numeric(Sum.days$Days.travl)


fn.poly=function(XLIM,YLIM,CL)polygon(x=c(XLIM,rev(XLIM)),y=c(YLIM[1],YLIM[1],YLIM[2],YLIM[2]),col=CL) 
fn.set.scence=function(YLM) plot(1:length(TG),1:length(TG),ylim=YLM,col='transparent',ann=F,xaxt='n',yaxt='n')
fn.add.polys=function(dat,var,YLIM,BY,BY1,MTEXT)
{
  id=match(var,names(dat))
  dat$CL=with(dat,ifelse(Dum=="North South",North.col,South.col))
  fn.set.scence(YLIM)
  for(i in 1:length(TG))
  {
    ss=subset(dat,TagCode==TG[i])
    ss$cum=cumsum(ss[,id])
    fn.poly(c(i,i+.5),c(0,ss[,id][1]),ss$CL[1])
    for(n in 2:nrow(ss)) fn.poly(c(i,i+.5),c(ss$cum[n-1],ss$cum[n]),ss$CL[n]) 
  }
  axis(1,(1:length(TG)+.25),F,tcl =0)
  axis(2,seq(0,YLIM[2],BY),F,tcl = .25)
  axis(2,seq(0,YLIM[2],BY1),las=1,tcl = .5,cex.axis=1.25)
  mtext(MTEXT,2,line=3,cex=1.4,las=3)
}

aa=Cum.dist.Tra.mig.dus$Round.trip
names(aa)=Cum.dist.Tra.mig.dus$TagCode
aa=rev(sort(aa))
TG=TG[match(names(aa),TG)]  #sort TG by number of displacements
 
tiff(file="Figure.Dist.Days.Mig.Dus.tiff",width=2000,height=2400,units="px",res=300,compression="lzw")

#Distance
par(mfcol=c(3,1),mar=c(2,1.95,1,0.1),oma=c(1.2,2.75,1,0.1),mgp=c(2, 0.25, 0),xpd=T,las=1)
fn.add.polys(Sum.dist,"Dist.moved.conseq.det",c(0,Migration.summary$max.km.travl.in.total),500,1000,"Distance travelled (km)")
legend('top',c("N-S movement","S-N movement"),bty='n',fill=c(North.col,South.col),cex=1.5)
legend('topright',"(a)",bty='n',cex=1.5)

#Time   
fn.add.polys(Sum.days,"Days.travl",
             c(0,1.03*max(as.numeric(aggregate(Days.travl~TagCode,Migrating.Dusky.dist.trav,sum)$Days.travl))),
             100,200,"Number of days")
#axis(1,(1:length(TG)+.25),1:length(TG),tcl =0,las=1,cex.axis=1.35,padj=.25)
legend('topright',"(b)",bty='n',cex=1.5)
  #add sex symbol
TG.sex=subset(Detections,TagCode%in%TG,select=c(TagCode,Sex,FL))
TG.sex=TG.sex[!duplicated(TG.sex$TagCode),]
TG.sex=TG.sex[match(TG,TG.sex$TagCode),]
Sym.vec=ifelse(TG.sex$Sex=="F","\\VE","\\MA")
#COL.SX=ifelse(TG.sex$Sex=="F","hotpink1","blue")
COL.SX=ifelse(TG.sex$Sex=="F","grey35","grey35")
text((1:length(TG)+.25),rep(-120,length(TG)),Sym.vec,vfont=c("sans serif","bold"),cex=1.95,col=COL.SX) 

Where.txt=aggregate(Days.travl~TagCode,Sum.days,sum)
Where.txt=Where.txt[match(TG,Where.txt$TagCode),]
text((1:length(TG)+.22),Where.txt$Days.travl*.9925,100*TG.sex$FL,cex=.9,pos=3)

mtext("Tagged individual",1,line=1.8,cex=1.4)


#ROM

Sum.dist$Dum=as.character(Sum.dist$Dum)
Sum.days$Dum=as.character(Sum.days$Dum)
Sum.ROM=merge(Sum.dist,Sum.days,by=c("TagCode","Stage","Zn.rec.prev","Zn.rec","Dum"),all=T)
Sum.ROM$ROM=Sum.ROM$Dist.moved.conseq.det/Sum.ROM$Days.travl

fn.add.ROM=function(dat,var,MTEXT)
{
  id=match(var,names(dat))
  dat$CL=with(dat,ifelse(Dum=="North South",North.col,South.col))
  TAB=with(dat,table(CL,floor(ROM/5)*5))
  LEG=colnames(TAB)
  for(q in 1:(length(LEG)-1))LEG[q]=paste(LEG[q],LEG[q+1],sep="-")
  LEG[length(LEG)]=paste(">",LEG[length(LEG)],sep="")
  barplot(TAB,col=rownames(TAB),names.arg=LEG,cex.axis=1.25,yaxt="n",cex.names=1.25)
  box()
  mtext(MTEXT,2,line=3,cex=1.4,las=3)
  YLIM=max(colSums(TAB))
  axis(2,seq(0,YLIM,5),F,tcl = .25)
  axis(2,seq(0,YLIM,10),las=1,tcl = .5,cex.axis=1.25)
}
fn.add.ROM(dat=Sum.ROM,var="ROM",MTEXT="Frequency")
mtext("Rate of movement (km/day)",1,line=2,cex=1.4)

legend('topright',"(c)",bty='n',cex=1.5)
dev.off()


#Figure timeline female and male together
# tiff(file="Figure Dusky.fem.male.tiff",width = 2400, height = 2400,
#      units = "px", res = 300,compression = "lzw")
# par(mai=c(.05,.5,.4,.01),oma=c(3,1,1,1),las=1,xpd=T,mgp=c(2,.1,0))
# m <- cbind(1:2,1:2)
# layout(m)
# fn.dusky.natal.mig(subset(Detections.species[[1]],Sex=="F"),c("Ningaloo","Perth","Southern.lines"),
#                    FL_0.5[[1]][1],HADJ=1.15,PADJ=0.5,c(-.5,-.5,39.3,39.3),Show.what="FL")
# fn.dusky.natal.mig(subset(Detections.species[[1]],Sex=="M"),c("Ningaloo","Perth","Southern.lines"),
#                    FL_0.5[[1]][2],HADJ=0.5,PADJ=-1,c(0.25,0.25,20.75,20.75),Show.what="FL")
# dev.off()


#Table of number of days detected per array   
dd=dat.glm  
dd$Day=1
dd=aggregate(Day~Array+Year+TagCode,dd,sum)
MIN.d=aggregate(Day~Array,dd,min)
MEAN.d=aggregate(Day~Array,dd,mean)
SD.d=aggregate(Day~Array,dd,sd)
MAX.d=aggregate(Day~Array,dd,max)
Table1=cbind(MIN.d,MEAN.d,SD.d,MAX.d)
Table1=Table1[,-c(3,5,7,9)]
colnames(Table1)=c("Array","Minimum","Mean","SD","Maximum")
write.csv(Table1,"Table1.csv",row.names=F)


#No detections   
Tag.dus=subset(TAGS,Species2=="Dusky" & Project.rel=="SMN"& !Code2%in%unique(Detections.species$Dusky$TagCode.original))
plot(Tag.dus$ReleaseLongitude2,Tag.dus$ReleaseLatitude2)
table(Tag.dus$Sex2)
hist(Tag.dus$ReleaseLength)

#Figure 1. Bubble of detections    
Tab.sks.arry=table(dat.glm$TagCode,dat.glm$Array)
Tab.sks.arry[Tab.sks.arry>0]=1
Tab.sks.arry=colSums( Tab.sks.arry)

fn.plt.mp=function(PlotlonG,PlotlatT,Add.depth,add.closure)
{
  plotMap(worldLLhigh, xlim=PlotlonG,ylim=PlotlatT,plt = c(.001, 1, 0.075, 1),
          col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  if(add.closure=="YES") plot(WA_Northern_Shark_2,add=T,col="grey85",border=1)
  if(Add.depth=="YES")
  {
    contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=YLIM,xlim=XLIM, zlim=c(-1,-200),
            nlevels = 3,labcex=.8,lty = 1,col=c("darkgrey","darkgrey","darkgrey","darkgrey","transparent"),add=T)
  }
  
}

#fn.poli=function(XLIM,YLIM,CL)polygon(x=c(XLIM,rev(XLIM)),y=c(YLIM[1],YLIM[1],YLIM[2],YLIM[2]),lwd=2,col=CL)
fn.poli=function(XLIM,YLIM,CL,LW)
{
  polygon(x=c(XLIM,rev(XLIM)),y=c(YLIM[1],YLIM[1],YLIM[2],YLIM[2]),lwd=LW,border=CL) 
}

L2=-35.3
L1=115
bubble.Dus.mig=function(scale,Pt.col,Pt.bg,Add.numb,nn,add.scale)
{
  Det=subset(Detections.species$Dusky,!Project.rel=="South.Australia")
  
  #Plot individuals
  Det$Dummy=with(Det,paste(TagCode,Species,Station))
  Det=Det[!duplicated(Det$Dummy),]
  DATA=Det
  TABLE=table(DATA$Station)
  N=length(unique(DATA$TagCode))
  MATRIX=matrix(names(TABLE))
  MATRIX=matrix(as.numeric(unlist(strsplit(MATRIX, split=" "))),ncol=2,byrow=T)
  LAT=MATRIX[,1]
  LONG=MATRIX[,2]
  zo=TABLE/N
  
  #add hits
  points(LONG,LAT,pch=21,cex=zo*scale,col=Pt.col,bg=Pt.bg)
  
  #add receivers
  points(STATIONS$longitude,STATIONS$latitude,col=1,pch=20,cex=1.25)
  
  #add number of sharks detected in array
  if(Add.numb=="YES")legend("bottomright",legend=paste("n= ",nn," individuals",sep=""),bty='n',cex=1.55)  

  if(add.scale=="YES")
  {
    points(L1,L2,pch=21,cex=(20/N)*scale,col=Pt.col,bg='transparent',lwd=2)
    points(L1,(L2+.11),pch=21,cex=(40/N)*scale,col=Pt.col,bg=Pt.bg,lwd=2.15)
    boxed.labels(L1,L2,'20',cex =1.15,col=Pt.col,bg='transparent',border=NA)
    boxed.labels(L1,L2+.25,'40',cex =1.15,col=Pt.col,bg='transparent',border=NA)
  }
}

colfunc <- colorRampPalette(c("yellow", "red"))
#if(GREYSCALE=="YES") colfunc <- colorRampPalette(c("grey70", "grey20"))

Heat.Dus.mig=function(Add.numb,nn,add.scale,where.num)
{
  Det=subset(Detections.species$Dusky,!Project.rel=="South.Australia")
  
  #Plot number of sharks
  Det$Dummy=with(Det,paste(TagCode,Species,Station))
  Det=Det[!duplicated(Det$Dummy),]
  DATA=Det
  TABLE=table(DATA$Station)
  N=length(unique(DATA$TagCode))
  MATRIX=matrix(names(TABLE))
  MATRIX=matrix(as.numeric(unlist(strsplit(MATRIX, split=" "))),ncol=2,byrow=T)
  LAT=MATRIX[,1]
  LONG=MATRIX[,2]
  Max.n.shark.station=max(TABLE)
  couleurs=colfunc(Max.n.shark.station)
  zo=couleurs[TABLE]
  
  #add receivers
  points(STATIONS$longitude,STATIONS$latitude,bg="white",pch=21,cex=1)
  
  #add number of sharks per receiver
  points(LONG,LAT,pch=21,cex=1.5,bg=zo,col=zo)
  
  #add number of sharks detected in array
  if(Add.numb=="YES")
  {
    legend(where.num,legend=paste(nn," sharks",sep=""),bty='n',text.font=2,cex=1.1,adj=c(0,-0.95))
    legend(where.num,legend=" detected by the array",bty='n',cex=1.1,text.font=2)
  }
      
  
  if(add.scale=="YES")
  {
    Leg.seq=seq(1,Max.n.shark.station,3)
    legend('topright',paste(Leg.seq),pch=19,col=couleurs[Leg.seq],horiz=T,bty='n',
           title="Number of sharks detected",cex=1.25)
    
  }
}

#c("blue","red","green","brown")
#  c("#0000ff22","#8B000022","#0000ff39","#8B000032")
Pt.bg="#8B000032"
Pt.col="brown"
COLOR="grey70"
#if(GREYSCALE=="YES")COLOR="grey95"


tiff(file="Figure 1.Dusky.natal.tiff",width = 2000, height = 2400,    
     units = "px", res = 300,compression = "lzw")
par(mar = c(3.2, 5, 0.95, 0),oma=c(1,.1,.1,1),mgp=c(1,.75,0))
m <- cbind(c(1:4,6,6),c(rep(5,4),6,6))
layout(m)

#Ningaloo
Point.list=list(first=list(c(113.88,-21.98)),
                second=list(c(113.5560,-22.61770),c(113.5824,-22.61595)),
                third=list(c(113.7138,-23.0309))
                )
txt.list=list(first=list(c(113.8824,-21.97483,50,45)),
              second=list(c(113.5560,-22.61770,200,55),
                          c(113.5824,-22.61595,100,55)),
              third=list(c(113.7138,-23.0309,50,55))
              )
Line.lab=c("Line 1","Line 2","Line 3")

plotlong.lines=list(c(113.8,114.04),c(113.54,113.68),c(113.57,113.865))
plotlat.lines=list(c(-22.00,-21.845),c(-22.64,-22.55),c(-23.16,-22.975))
Ad=c("NO","YES","NO")
for(j in 1:3) 
{
  fn.plt.mp(plotlong.lines[[j]],plotlat.lines[[j]],"YES","NO")
  
  #fn.poli(c(113,114.5),c(-23.5,-21.5),rgb(.01,.8,.1,alpha=.2))
  #bubble.Dus.mig(scale=20,Pt.col=Pt.col,Pt.bg=Pt.bg,Add.numb=Ad[j],nn=Tab.sks.arry[1],"NO")
  Heat.Dus.mig(Add.numb=Ad[j],nn=Tab.sks.arry[1],add.scale="NO",where.num="bottomright")
  mtext(Line.lab[j],2,las=3,0)
  #add depth lables
  if(j==3)
  {
     TxT=txt.list[[j]]
    for(t in 1:length(TxT))
    {
      points(Point.list[[j]][[t]][1],Point.list[[j]][[t]][2],cex=4,pch=22,col="white",bg="white")
      text(TxT[[t]][1],TxT[[t]][2],TxT[[t]][3],col="grey70",srt=TxT[[t]][4],cex=1.1)
    }
  }

  fn.poli(plotlong.lines[[j]],plotlat.lines[[j]],"grey30",5)
}

#Perth
fn.plt.mp(c(115.1,115.95),c(-32.3,-31.9),"YES","NO")
#fn.plt.mp(c(115,116.75),c(-32.5,-31.5),"YES","NO")

fn.poli(c(115.1,115.95),c(-32.3,-31.9),"grey30",5)
#fn.poli(c(115,116.75),c(-32.5,-31.5),rgb(.01,.01,.9,alpha=.2))
#bubble.Dus.mig(scale=20,Pt.col=Pt.col,Pt.bg=Pt.bg,Add.numb="YES",nn=Tab.sks.arry[2],"NO")
Heat.Dus.mig(Add.numb="YES",nn=Tab.sks.arry[2],add.scale="NO",where.num="bottomleft")
polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")

# points(115.2744,-32.40167,cex=3,pch=22,col="transparent",bg="white")
# text(115.2744,-32.40167,50,col="grey70",srt=55,cex=1.1)
# 
# points(115.1705,-32.39915,cex=3,pch=22,col="transparent",bg="white")
# text(115.1705,-32.39915,100,col="grey70",srt=55,cex=1.1)
# 
# points(115.0873,-32.30095,cex=3,pch=22,col="transparent",bg="white")
# text(115.0873,-32.30095,200,col="grey70",srt=55,cex=1.1)


fn.poli(c(115,116.75),c(-32.5,-31.5),'grey30',5)

#WA
XX=c(112,119)
YY=c(-35.5,-21.5)
LONGSEQ=XX[1]:XX[2]
LATSEQ=seq(ceiling(YY[1]),ceiling(YY[2]),by=1)
LATSEQ2=seq(ceiling(YY[1]),ceiling(YY[2]),by=2)
fn.plt.mp(XX,YY,"NO",'YES')
points(STATIONS$longitude,STATIONS$latitude,col=1,pch=20,cex=1.25)
axis(side = 1, at =LONGSEQ, labels = paste(LONGSEQ,"ºE",sep=""), tcl = .5,las=1,cex.axis=1.5,padj=-.5)
axis(side = 4, at = LATSEQ, labels =F,tcl = .5,las=2,cex.axis=1.5,hadj=.15)
axis(side = 4, at = LATSEQ2, labels =paste(-LATSEQ2,"ºS",sep=""),tcl = .5,las=2,cex.axis=1.5,hadj=.15)
box(lwd=1.5)
fn.scale.bar(c(116,-25.5),c(117,-25.5),1.004,.9875,1.005,1.01)
legend(115.5,-28,c("Receiver","Releases"),pch=c(19,4),col=c("black","grey15"),cex=1.25,bty='n')

#releases 
RELes=subset(TAGS,Species2=="Dusky" & Project.rel=="SMN")
points(RELes$ReleaseLongitude2,RELes$ReleaseLatitude2,pch="X",cex=1,col="grey15")

#turning points
points(Cape.Leuwin[1],Cape.Leuwin[2],pch=21,col=1,bg="white",cex=1.75)
points(Shark.bay[1],Shark.bay[2],pch=21,col=1,bg="white",cex=1.75)
points(Exmouth[1]+.2,Exmouth[2],pch=21,col=1,bg="white",cex=1.75)     


#polygons
fn.poli(c(113,114.5),c(-23.5,-21.5),'grey30',2)
fn.poli(c(115,116.75),c(-32.5,-31.5),'grey30',2)
fn.poli(c(114.5,118.7),c(-35.5,-34),'grey30',2)
# fn.poli(c(113,114.5),c(-23.5,-21.5),rgb(.01,.8,.1,alpha=.2))
# fn.poli(c(115,116.75),c(-32.5,-31.5),rgb(.01,.01,.9,alpha=.2))
# fn.poli(c(114.5,118.7),c(-35.5,-34),rgb(.9,.01,.01,alpha=.2))
text(117,-31,"Perth array",cex=1.65,col='grey10',font=2)
text(117.5,-33.25,"Southern",cex=1.65,col='grey10',font=2)
text(117.5,-33.75,"array",cex=1.65,col='grey10',font=2)
text(116.5,-22,"Ningaloo",cex=1.65,col='grey10',font=2)
text(116.5,-22.5,"array",cex=1.65,col='grey10',font=2)

#Southern lines
fn.plt.mp(c(114.5,118.7),c(-35.5,-34),"YES","NO")
fn.poli(c(114.5,118.7),c(-35.5,-34),'grey30',5)
#fn.poli(c(114.5,118.7),c(-35.5,-34),rgb(.9,.01,.01,alpha=.2))
#bubble.Dus.mig(scale=20,Pt.col=Pt.col,Pt.bg=Pt.bg,Add.numb="YES",nn=Tab.sks.arry[3],"YES")
Heat.Dus.mig(Add.numb="YES",nn=Tab.sks.arry[3],add.scale="YES",where.num="bottomleft")

par(xpd=NA)
arrows(117.42,-30.85,116.1,-31.18, length = 0.25, angle = 30,lwd=2,col='grey30')
arrows(117.35,-31,116.1,-31.96, length = 0.25, angle = 30,lwd=2,col='grey30')
arrows(117.36,-31.12,116.1,-32.73, length = 0.25, angle = 30,lwd=2,col='grey30')
arrows(117.7,-33.05,116.1,-33.73, length = 0.25, angle = 30,lwd=2,col='grey30')

#Inset Oz
par(fig=c(.75,.85,.85,.95), new = T,mgp=c(.1,.4,0))
fn.plt.mp(c(111,156),c(-45,-10),"NO","NO")
box(lwd=1)
polygon(x=c(XX,rev(XX)),y=c(YY[1],YY[1],YY[2],YY[2]),lwd=1.5,col="transparent",border=1)
dev.off()



RELes$Year=year(RELes$DATE)
RELes$Zone=with(RELes,ifelse(ReleaseLatitude2>=(-26),"North","South"))
A=table(as.character(RELes$Zone),as.character(RELes$Year))
#table(RELes$Sex2)

# QUANTIFY AND PREDICT MIGRATION

#19.2.2 Logistic glm for predicting size and timing of migration  

#note: Export data and run this analysis externally for sharing with Simon in GitHub
write.csv(dat.glm,"dat.glm.csv",row.names=F)

#check each shark used  
x1=c(0,max(dat.glm$yday.rec))
x2=c(0,max(dat.glm$week))
x3=as.numeric(c(min(dat.glm$Date.local),max(dat.glm$Date.local)))
pdf(handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement/Natal_migration/Each_shark/plot.pdf"))
for(n in 1:length(Tags.d))
{
  a=subset(dat.glm,TagCode==Tags.d[n],select=c(TagCode,Date.local,Year,Month,week,yday.rec,Sesn.rec,Array,N,Sex,FL))
  a=a[order(a$Date.local),]
  a$col=with(a,ifelse(N==1,2,4))
#  tiff(file=paste("C:/Matias/Analyses/Acoustic_tagging/FRDC/Outputs_movement/Natal_migration/Each_shark/",Tags.d[n],".tiff",sep=""),
#       width = 2000, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(3,1),mar=c(4,4,.1,.1),mgp=c(2,.6,0))
  plot(a$yday.rec,a$N,pch=19,xlab="ydays",ylab="North",col=a$col,cex.axis=1.25,cex.lab=1.5,cex=2,
       ylim=c(-0.1,1.1),xlim=x1)
  plot(a$week,a$N,pch=19,xlab="week",ylab="North",col=a$col,cex.axis=1.25,cex.lab=1.5,cex=2,
       ylim=c(-0.1,1.1),xlim=x2)
  legend("center",paste(Tags.d[n]," (",unique(a$Sex),", ",unique(a$FL)," cm)",sep=""),bty='n',cex=3)
  plot(a$Date.local,a$N,pch=19,xlab="date",ylab="North",col=a$col,cex.axis=1.25,cex.lab=1.5,cex=2,
       ylim=c(-0.1,1.1),xlim=x3)  
#  dev.off()
}
dev.off()


#Other modelling approaches trielled 
Do.other.approach="NO"
if(Do.other.approach=="YES")
{
  source(handl_OneDrive("Analyses/Acoustic_tagging/Dusky_migration/Git_dusky_migration/Dusky_migration.R"))  
  
  #19.2.3 GLM approach not used
  
  #19.2.3.1. GLM. presence North by FL, day length, temperature and Sex (incomplete) 
  do.glm.sex.day.len.FL='NO'
  if (do.glm.sex.day.len.FL=="YES")
  {
    
    daylength1 <- function(t,lon, lat)
    {
      t <- as.numeric(t)
      alt <- function(t) sunAngle(t, longitude=loN, latitude=laT)$altitude
      rise=set=rep(NA,length(t))
      for(i in 1:length(t))
      {
        loN=lon[i]
        laT=lat[i]
        rise[i]=uniroot(alt, lower=t[i]-86400/2, upper=t[i])$root
        set[i] <- uniroot(alt, lower=t[i], upper=t[i]+86400/2)$root
      }
      daylen=set - rise
      return(daylen/3600)
    }
    dat.glm$day.length=daylength1(dat.glm$Date.local,dat.glm$Longitude,dat.glm$Latitude) 
    
    #add Temp
    Mean.temp1.N$Zn.rec="North"
    Mean.temp1.S$Zn.rec="South"
    Mean.TemP=rbind(Mean.temp1.N,Mean.temp1.S)
    dat.glm=merge(dat.glm,Mean.TemP,by=c("Zn.rec","Year","Month"),all.x=T)
    
    dat.glm$Pos=with(dat.glm,ifelse(Zn.rec=="North",1,0))
    
    #It should be a GAM because day.length and Temperature are cyclical!!
    #Also, issues with Temperature!! see  plot(dat.glm$Temperature,dat.glm$Pos)
    MOD=glm(Pos~FL*as.factor(Sex)+day.length+Temperature, data=dat.glm, family="binomial", maxit=500)
    MOD2=glmer(Pos~FL*as.factor(Sex)+day.length +(1 | TagCode), data=dat.glm, family="binomial")
    
    
    
    with(dat.glm,table(N,week))
    
    mod1<- glm(N~FL*day.length*as.factor(Sex), data=dat.glm, family="binomial", maxit=500)
    mod2<- glm(N~FL*as.factor(week)+as.factor(Sex), data=dat.glm, family="binomial", maxit=500)
    
    new.day.len=c(rev(seq(10,14,1)),seq(11,14,1))
    NEWdat=expand.grid(FL=seq(1.5,3,.5),day.length=new.day.len,
                       Sex=factor("M",levels=c("F","M")))
    
    NEWdat$pred=predict(mod1,newdata=NEWdat,type='response')
    
    Len.Cls=unique(NEWdat$FL)
    plot(as.character(new.day.len),new.day.len,ylim=c(0,1))
    for(x in 1:length(Len.Cls))
    {
      a=subset(NEWdat,FL==Len.Cls[i])
      lines()
    }
    
  }
  
  
  #19.2.3.2. GLM. Number of detections per day by Size and Depth
  dat.glm$hit=1
  agg=aggregate(hit~FL+Depth+TagCode,dat.glm,sum)
  plot(agg$FL,agg$hit)
  plot(agg$Depth,agg$hit)
  if (do.glm.sex.day.len.FL=="YES")
  {
    mod1=glm(hit~ FL+Depth+TagCode,data=agg,family=Gamma(link=log))
    #mod2=glmer(hit~ FL+Depth+(1|TagCode),data=agg,family=Gamma(link=log)) 
  }
  cor.test(agg$FL,agg$hit,alternative="two.sided")
  
  
  #19.2.4. Population dynamics approach 
  Do.pop.dyn.apprach="NO"
  if(Do.pop.dyn.apprach=="YES")
  {
    #growth parameters
    Linf.f=374.4
    K.f=.0367
    Linf.m=336.5
    K.m=.045
    Lo=75.3
    age=0:35
    mid.FL.fem=Lo+(Linf.f-Lo)*(1-exp(-K.f*age))
    mid.FL.male=Lo+(Linf.m-Lo)*(1-exp(-K.m*age))
    plot(age,mid.FL.fem)
    lines(age,mid.FL.male)
    
    #data
    DAT=vector('list',length(Tags.d))
    for(i in 1:length(Tags.d))
    {
      a=subset(dat.glm,TagCode==Tags.d[i])
      a=a[order(a$Date.local),]
      a=subset(a,select=c(TagCode,Sex,FL,Year.rel,Sesn.rel,Zn.rel,Year,Sesn.rec,Zn.rec,Date.local))
      names(a)[match("Year",names(a))]="Year.rec"
      a.rel=a[1,]
      b=a[-1,]
      if(nrow(b)==1)
      {
        b$Year.rel=a.rel$Year.rec
        b$Sesn.rel=a.rel$Sesn.rec
        b$Zn.rel=a.rel$Zn.rec
      }
      if(nrow(b)>1)
      {
        n=nrow(b)-1
        b$Year.rel=c(a.rel$Year.rec,b$Year.rec[1:n])
        b$Sesn.rel=c(a.rel$Sesn.rec,b$Sesn.rec[1:n])
        b$Zn.rel=c(a.rel$Zn.rec,b$Zn.rec[1:n])
      }
      if(nrow(b)==0)  a=a.rel
      if(nrow(b)>0)
      {
        a=rbind(a.rel,b)
        a=a[order(a$Date.local),]
      } 
      DAT[[i]]=a
    }
    
    DAT=do.call(rbind,DAT)
    
    #aggregate observations by terms
    DAT$N=1
    DAT$bin.resp=with(DAT,ifelse(Zn.rec=='North',1,0))
    
    Agg.dat.glm=aggregate(N~TagCode+Sex+FL +Year.rel+Sesn.rel+Zn.rel+Year.rec+Sesn.rec+Zn.rec,DAT,sum)
    Agg.dat.glm=reshape(Agg.dat.glm,v.names = "N", idvar = c("TagCode","Sex","FL","Year.rel",
                                                             "Sesn.rel","Zn.rel","Year.rec","Sesn.rec"),timevar = "Zn.rec", direction = "wide")
    Agg.dat.glm[is.na(Agg.dat.glm)]=0
    Agg.dat.glm1=Agg.dat.glm[order(Agg.dat.glm$TagCode,Agg.dat.glm$Year.rel,Agg.dat.glm$Year.rec),]
    
    Agg.dat.glm1$Prop.rec.N=Agg.dat.glm1$N.North/(Agg.dat.glm1$N.North+Agg.dat.glm1$N.South)
    Agg.dat.glm1$Prop.rec.S=Agg.dat.glm1$N.South/(Agg.dat.glm1$N.North+Agg.dat.glm1$N.South)
    Agg.dat.glm1=Agg.dat.glm1[,-match(c("N.North","N.South"),names(Agg.dat.glm))]
    
    write.csv(Agg.dat.glm1,"Data.Simon.csv",row.names=F)
    write.csv(dat.glm,"dat.glm.csv",row.names=F)       
    
    dat.glm$Sex=factor(as.character(dat.glm$Sex),levels=c("F","M"))
    dat.glm$TagCode=as.factor(dat.glm$TagCode)
    
  }
  
  
  #19.2.5 run and source Simon's movement model 
  #source("C:/Matias/Analyses/Acoustic_tagging/FRDC/Outputs_movement/Natal_migration/Movement model.R")
  
  
  #19.2.6  Apply double logistic models to determine seasons
  Do.double.logistic="NO"
  if(Do.double.logistic=="YES")
  {
    Mal=subset(dat.glm,Sex=="M",select=c(TagCode,FL,Zn.rec,Date.local,Month,N,yday.rec,week,day.length,Temperature,Temp.anomaly))
    Fem=subset(dat.glm,Sex=="F",select=c(TagCode,FL,Zn.rec,Date.local,Month,N,yday.rec,week,day.length,Temperature,Temp.anomaly))   
    
    Calculate.Combined.process="NO"
    
    
    # a) Two tier approach: first calculate probability of moving, then prob of timing probability   
    if(Calculate.Combined.process=="NO")
    {
      #functions
      Logis.type='Inflection'
      #Logis.type='Percent.50'
      
      if(Logis.type=='Inflection')
      {
        fn.logis=function(dat,pmax,inflx,slop) pmax/(1+exp((dat-inflx)/slop))
        plot(fn.logis(100:300,1,180,-10))
        
        #initial par value
        #size at migration
        inflx.m=2.1
        slop.m=-.1
        Pmax.m=1
        
        inflx.f=2.4
        slop.f=-.1
        Pmax.f=1
        
        theta.m=c(inflx=inflx.m,slop=slop.m)
        theta.f=c(inflx=inflx.f,slop=slop.f)
        
        
        #time of migration
        inflx.up.m=160
        slop.m=5
        inflx.dwn.m=280
        P.max.m=.7
        plot(fn.logis(0:465,1,180,-5))
        
        theta1.m=c(inflx.up=inflx.up.m,slop=slop.m,inflx.dwn=inflx.dwn.m,P.max=P.max.m)
        theta1.f=theta1.m
        
      }
      
      if(Logis.type=='Percent.50')
      {
        fn.logis=function(dat,pmax,p50,p95) 1/(1+exp(-log(19)*(dat-p50)/(p95-p50)))
        plot(fn.logis(100:300,1,150,180)) 
        fn.logis.down=function(dat,pmax,p50,p95) 1/(1+exp(log(19)*(dat-p50)/(p95-p50)))
        
        #initial par value
        #size at migration
        L50.m=2.1
        L95.m=2.3
        Pmax.m=1
        
        L50.f=2.4
        L95.f=2.8
        Pmax.f=1
        
        theta.m=c(L50=L50.m,L95=L95.m)
        theta.f=c(L50=L50.f,L95=L95.f)
        
        
        #time of migration
        L50.up.m=rnorm(1,160,10)
        L95.up.m=rnorm(1,180,10)
        L50.dwn.m=rnorm(1,280,10)
        L95.dwn.m=rnorm(1,310,10)
        P.max.m=rnorm(1,.7,.1)
        
        theta1.m=log(c(L50.up=L50.up.m,L95.up=L95.up.m,L50.dwn=L50.dwn.m,L95.dwn=L95.dwn.m,P.max=P.max.m))
        theta1.f=theta1.m
      }
      
      Months=as.POSIXlt(c("2011-01-01","2011-02-01","2011-03-01","2011-04-01","2011-05-01","2011-06-01",
                          "2011-07-01","2011-08-01","2011-09-01","2011-10-01","2011-11-01","2011-12-01"))
      Mn=Months$yday
      names(Mn)=Months
      MonthS=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
      
      
      fn.natal.pars=function(what,theta,SEQ,theta1,use.release)
      {
        Seasonal=what
        if(use.release=="YES")
        {
          Dummy=with(what,table(TagCode,Zn.rel))
          This.north=names(Dummy[Dummy[,1]>0,1])
          Dummy=with(what,table(TagCode,Zn.rec))
          This.north=unique(c(names(Dummy[Dummy[,1]>0,1]),This.north))
        }
        
        #remove release zone to become independent of sampling effort
        if(use.release=="NO")
        {
          Dummy=with(what,table(TagCode,Zn.rec))
          This.north=names(Dummy[Dummy[,1]>0,1])
        }
        
        
        # Step1. Determine size at migration  
        what=what[!duplicated(what$TagCode),]
        what$N=with(what,ifelse(TagCode%in%This.north,1,0))  
        #fit model
        objfun=function(theta)
        {
          what$Prob.n = fn.logis(what$FL,1,theta[1],theta[2])
          what$Like.n <-  with(what,ifelse(N==1,Prob.n,(1-Prob.n)))
          what$log.like=log(what$Like.n+1e-100)
          ll=sum(-what$log.like)
          return(ll)
        }
        fit=optim(theta,objfun, hessian=T, control=c(trace=1, maxit=10000))   
        #mod1 <- nls(N~1/(1+exp((FL-a)/b)), start=c(a=2.1, b=-2), data=what)
        #summary(mod1)
        PAR=fit$par
        v_ob=solve(fit$hessian)  #variance covariance matrix
        std_ob=sqrt(diag(v_ob))
        OUT.size=data.frame(PAR,std_ob)
        
        #plot
        what=what[order(what$FL),]
        plot(what$FL,what$N,main=paste(unique(what$Sex),"(use.release=",use.release,")"),
             xlab="FL (m)",ylab="Proportion North",pch=19,cex=1.25)
        lines(SEQ,fn.logis(SEQ,1,PAR[1],PAR[2]),col=2,lwd=2)
        
        
        # Step2. determine seasonal parameters
        Seasonal=subset(Seasonal,TagCode%in%This.north)  
        Seasonal$S=with(Seasonal,ifelse(Zn.rec=="South",1,0))
        TG=unique(Seasonal$TagCode)
        kEP=rep(NA,length(TG))
        kep="n.s"
        
        #round-trip migrations (N-S-N or S-N-S)    
        if(kep=="round.trip")
        {
          for(s in 1:length(TG))
          {
            a=subset(Seasonal,TagCode==TG[s])
            a=a[order(a$Month,a$Zn.rec),]
            a$cum=0
            if(nrow(a)>1)for(w in 2:nrow(a)) if(!a$Zn.rec[w]==a$Zn.rec[w-1]) a$cum[w]=1
            x=NA
            if(sum(a$cum)>1) x=TG[s]
            kEP[s]=x
          }
          kEP=subset(kEP,!is.na(kEP))
        }
        
        #North and South occurrence
        if(kep=="n.s")
        {
          for(s in 1:length(TG))
          {
            a=subset(Seasonal,TagCode==TG[s])
            UNIK=unique(a$Zn.rec)
            x=NA
            if(length(UNIK)>1) x=TG[s]
            kEP[s]=x
          }
          kEP=subset(kEP,!is.na(kEP))
        }
        
        Seasonal=subset(Seasonal,TagCode%in%kEP)
        This.north=unique(Seasonal$TagCode)
        Seasonal$day=Seasonal$Date.local$yday   #days from beginning of year
        Seasonal$week=week(Seasonal$Date.local)
        Agg.wk=aggregate(cbind(N,S)~week+TagCode,Seasonal,sum)
        Agg.wk$Prop.N=with(Agg.wk,N/(N+S))
        Agg.wk=Agg.wk[order(Agg.wk$week),]
        plot(Agg.wk$week,Agg.wk$Prop.N,type='l')
        
        Dummy=with(Seasonal,table(day,Zn.rec,TagCode)) 
        LST=vector('list',length(This.north))
        names(LST)=This.north
        for(t in 1:length(This.north))
        {
          A=Dummy[,,t]
          indx=which(rowSums(A)>0)
          Dys=as.numeric(as.character(rownames(A)[indx]))
          LST[[t]]=data.frame(TagCode=This.north[t],day=Dys,North=A[indx,1],South=A[indx,2])
        }
        Dummy=do.call(rbind,LST)
        Dummy$North=with(Dummy,ifelse(North>1,1,North))
        Dummy$South=with(Dummy,ifelse(South>1,1,South))
        
        #Prob north
        #fit model
        objfun=function(theta1)
        {
          if(Logis.type=='Inflection') theta2=theta1
          if(Logis.type=='Percent.50') theta2=exp(theta1)
          
          Dummy$Prob.n.up = fn.logis(Dummy$day,theta2[4],theta2[1],-theta2[2])
          if(Logis.type=='Inflection') Dummy$Prob.n.dwn = fn.logis(Dummy$day,theta2[4],theta2[3],theta2[2])
          if(Logis.type=='Percent.50') Dummy$Prob.n.dwn = fn.logis.down(Dummy$day,theta2[4],theta2[3],theta2[2])
          Dummy$Prob.n=Dummy$Prob.n.up*Dummy$Prob.n.dwn
          Dummy$Like.n=with(Dummy,ifelse(North==1,Prob.n,(1-Prob.n)))
          Dummy$log.like=log(Dummy$Like.n+1e-100)
          ll=sum(-Dummy$log.like)
          return(ll)
        }  
        fit.month=optim(theta1,objfun, hessian=T, control=c(trace=1, maxit=10000))   
        v_ob=solve(fit.month$hessian)  #variance covariance matrix
        
        if(Logis.type=='Inflection')
        {
          PAR=fit.month$par
          std_ob=sqrt(diag(v_ob))
        }
        
        if(Logis.type=='Percent.50')
        {
          PAR=exp(fit.month$par)
          std_ob=exp(sqrt(diag(v_ob)))
        }
        
        OUT.time.north=data.frame(PAR,std_ob)
        
        #Plot
        Dummy=Dummy[order(Dummy$day),]
        plot(Dummy$day,Dummy$North,xlim=c(0,365),pch=19,cex=1.25,ylab="Proportion north",xlab="Day",
             xaxt='n')
        if(Logis.type=='Inflection')Preds=fn.logis(0:365,PAR[4],PAR[1],-PAR[2])*fn.logis(0:365,PAR[4],PAR[3],PAR[2])
        if(Logis.type=='Percent.50')Preds=fn.logis(0:365,PAR[5],PAR[1],PAR[2])*fn.logis.down(0:365,PAR[5],PAR[3],PAR[4])
        lines(0:365,Preds,col=2,lwd=2)
        axis(1,Mn,MonthS,tcl = .5)
        
        
        return(list(OUT.size=OUT.size,OUT.time.north=OUT.time.north))
      }
      
      
      SeQ=seq(1.6,2.7,by=.1)
      
      par(mfcol=c(2,1))
      Male.out=fn.natal.pars(what=Mal,theta=theta.m,SeQ,theta1=theta1.m,use.release="NO")
      
      
      Female.out=fn.natal.pars(what=Fem,theta=theta.f,SeQ,theta1=theta1.f,use.release="NO")
      
    }
    
    
    # b) Combined probability of moving and timing probability
    #note: the process is modelled properly but the data do not support the estimation
    if(Calculate.Combined.process=="YES")
    {
      what=Mal
      knife.edge="YES"
      if(knife.edge=="NO") theta.m=log(c(inflx=210,slop=.1,inflx.up.t=160,slop.t=1,inflx.dwn.t=280,P.max.t=.7))
      if(knife.edge=="YES") theta.m=log(c(inflx=210,inflx.up.t=160,inflx.dwn.t=280,P.max.t=.7))
      theta=theta.m
      use.release="NO"
      
      fn.natal.pars=function(what,theta,SEQ,use.release)
      {
        what$FL=100*what$FL
        
        Seasonal=what
        if(use.release=="YES")
        {
          Dummy=with(what,table(TagCode,Zn.rel))
          This.north=names(Dummy[Dummy[,1]>0,1])
          Dummy=with(what,table(TagCode,Zn.rec))
          This.north=unique(c(names(Dummy[Dummy[,1]>0,1]),This.north))
        }
        
        #remove release zone to become independent of sampling effort
        if(use.release=="NO")
        {
          Dummy=with(what,table(TagCode,Zn.rec))
          This.north=names(Dummy[Dummy[,1]>0,1])
        }
        
        
        # size at migration data  
        what=what[!duplicated(what$TagCode),]
        what$N=with(what,ifelse(TagCode%in%This.north,1,0))  
        
        # time of migration data
        Seasonal$day=Seasonal$Date.local$yday   #days from beginning of year
        Dummy=with(Seasonal,table(day,Zn.rec,TagCode)) 
        This.shk=unique(Seasonal$TagCode)
        LST=vector('list',length(This.shk))
        for(t in 1:length(This.shk))
        {
          A=Dummy[,,t]
          indx=which(rowSums(A)>0)
          Dys=as.numeric(as.character(rownames(A)[indx]))
          LST[[t]]=data.frame(TagCode=This.shk[t],day=Dys,North=A[indx,1],South=A[indx,2])
        }
        Dummy=do.call(rbind,LST)
        Dummy$North=with(Dummy,ifelse(North>1,1,North))
        Dummy$South=with(Dummy,ifelse(South>1,1,South))
        
        #fit model
        objfun=function(theta)
        {
          theta1=exp(theta)
          
          #size at migration process
          if(knife.edge=="YES") what$Prob.size = fn.logis(what$FL,1,theta1[1],-0.1)
          if(knife.edge=="NO") what$Prob.size = fn.logis(what$FL,1,theta1[1],-theta1[2])
          #     what$Like.n <-  with(what,ifelse(N==1,Prob.size,(1-Prob.size)))
          #     what$log.like=log(what$Like.n+1e-100)
          #     ll=sum(-what$log.like)
          
          #movement process
          #improve: grow sharks...
          Dummy=merge(Dummy,what[,match(c("TagCode","Prob.size"),names(what))],by="TagCode",all.x=T)
          if(knife.edge=="NO") 
          {
            Dummy$Prob.n.up = fn.logis(Dummy$day,theta1[6],theta1[3],-theta1[4])
            Dummy$Prob.n.dwn = fn.logis(Dummy$day,theta1[6],theta1[5],theta1[4])
          }
          
          if(knife.edge=="YES") 
          {
            Dummy$Prob.n.up = fn.logis(Dummy$day,theta1[4],theta1[2],-.1)
            Dummy$Prob.n.dwn = fn.logis(Dummy$day,theta1[4],theta1[3],.1)
          }
          
          Dummy$Prob.n=Dummy$Prob.n.up*Dummy$Prob.n.dwn
          Dummy$PROB=Dummy$Prob.size*Dummy$Prob.n
          Dummy$Like.n=with(Dummy,ifelse(North==1,PROB,(1-PROB)))
          Dummy$log.like=log(Dummy$Like.n+1e-100)
          ll=sum(-Dummy$log.like)
          return(ll)
        }
        fit=optim(theta,objfun, hessian=T, control=c(trace=1, maxit=10000))
        v_ob=solve(fit$hessian)  #variance covariance matrix
        std_ob=exp(sqrt(diag(v_ob)))
        PAR=exp(fit$par)
        OUT.size=data.frame(PAR,std_ob)
        
        #plot
        what=what[order(what$FL),]
        plot(what$FL,what$N,main=paste(unique(what$Sex),"(use.release=",use.release,")"),
             xlab="FL (cm)",ylab="Proportion North",pch=19,cex=1.25)
        #lines(what$FL,what$Prob.size)
        lines(SEQ,fn.logis(SEQ,1,PAR[1],-PAR[2]))
        
        Dummy=Dummy[order(Dummy$day),]
        plot(Dummy$day,Dummy$North,xlab="day",ylab="Proportion North",pch=19,cex=1.25)
        #lines(Dummy$day,Dummy$Prob.n.up,col=2)
        #lines(Dummy$day,Dummy$Prob.n.dwn,col=2)
        #lines(Dummy$day,Dummy$Prob.n,col=3)
        lines(0:365,fn.logis(0:365,PAR[6],PAR[3],-PAR[4])*fn.logis(0:365,PAR[6],PAR[5],PAR[4]),col=2)
        
        
        return(OUT)
      }
    }
    
  }
  
  
}

setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC"))

#19.3 Proportion of detections by month per array for dusky sharks

dummy1=subset(Detections,Species=="Dusky" & !is.na(Array))
dummy1$Array=with(dummy1,ifelse(Array=="Ningaloo","Ningaloo","Perth.or.Southern"))
Month.prop.time=aggregate(hit~Array+Month+TagCode+Sex,dummy1,sum)
Month.prop.time1=aggregate(hit~Month+TagCode+Sex,dummy1,sum)
names(Month.prop.time1)[4]="Max"
Month.prop.time=merge(Month.prop.time,Month.prop.time1,by=c("Month","TagCode","Sex"),all.x=T)
Month.prop.time$Prop=with(Month.prop.time,hit/Max)
Month.prop.time=subset(Month.prop.time,!Sex=="U")


# Month.prop.time.keep=aggregate(hit~Array+TagCode,dummy1,sum)
# Month.prop.time.keep$hit=1
# Month.prop.time.keep=aggregate(hit~TagCode,Month.prop.time.keep,sum)
# Month.prop.time.keep=subset(Month.prop.time.keep, hit>1)
# Month.prop.time=subset(Month.prop.time,TagCode%in%Month.prop.time.keep$TagCode)
fn.prop.mn=function(dat)
{
  sk=unique(dat$TagCode)
  sk1=subset(sk,sk%in%names(Arys.det))
  sk2=subset(sk,sk%in%Arys.det.Ning.only)
  sk3=subset(sk,sk%in%Arys.det.Perth.South.only)
  sk=c(sk1,sk2,sk3)
  Arys=unique(dat$Array)
  for (k in 1:length(sk))
  {
    aa=subset(dat,TagCode==sk[k])
    aa=aa[order(aa$Month),]
    a=reshape(subset(aa,select=c(Month,Array,Prop)),v.names="Prop",timevar="Month",
                idvar="Array",direction="wide")
    if(nrow(a)==1)
    {
      s=Arys[which(!Arys%in%a$Array)]
      b=a
      b$Array=s
      b[2:ncol(b)]=NA
      a=rbind(a,b)
    }
    colnames(a)[2:ncol(a)]=unique(aa$Month)
    a=a[order(a$Array),]
    Mn=1:12
    ID=Mn[which(!Mn%in%aa$Month)]
    add=matrix(nrow=nrow(a),ncol=length(ID))
    colnames(add)=ID
    this=as.matrix(a[,2:ncol(a)])
    if(ncol(this)==1)colnames(this)=colnames(a)[2]
    this=cbind(this,add)
    this=this[,match(Mn,colnames(this))]
    this[is.na(this)]=0
    row.names(this)=a$Array
    barplot(this,main=unique(aa$TagCode),col=c(3,4))
    box()
  }
}

Arys.det=with(Month.prop.time,table(TagCode,Array))
Arys.det[Arys.det>0]=1
a=as.data.frame.matrix(Arys.det)
a$TagCode=row.names(a)
Arys.det.Ning.only=subset(a, Ningaloo==1 & Perth.or.Southern==0,select=TagCode)[,1]
Arys.det.Perth.South.only=subset(a, Ningaloo==0 & Perth.or.Southern==1,select=TagCode)[,1]
Arys.det=rowSums(Arys.det)
Arys.det=Arys.det[Arys.det>1]

nn.F=subset(Month.prop.time,Sex=="F")
nn.F=unique(nn.F$TagCode)
#nn.F=subset(nn.F,!nn.F==29600) #drop single detection one for displaying purposes
nn.M=subset(Month.prop.time,Sex=="M")
nn.M=unique(nn.M$TagCode)
tiff(file="Outputs_movement/Migration.sync.females.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(7,7),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.65, 0),las=1)
fn.prop.mn(subset(Month.prop.time,Sex=="F" & TagCode%in%nn.F))
mtext("Proportion of detected time",2,las=3,line=0.5,outer=T,cex=1.5)
mtext("Month",1,outer=T,line=0.5,cex=1.5)
dev.off()


tiff(file="Outputs_movement/Migration.sync.males.tiff",width = 2000, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(6,5),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.65, 0),las=1)
fn.prop.mn(subset(Month.prop.time,Sex=="M"))
#mtext("Males",3,line=0.5,cex=1.25)
mtext("Proportion of detected time",2,las=3,line=0.5,outer=T,cex=1.5)
mtext("Month",1,outer=T,line=0.5,cex=1.5)
dev.off()



#20. --  Compare sex ratio, condition and size of individuals tagged in WA vs detected  ---
Condition$Species=as.character(Condition$SPECIES)
Condition$Species=with(Condition,ifelse(Species=="BW","Dusky",ifelse(Species=="TK","Thickskin",
                      ifelse(Species=="GM","Gummy",ifelse(Species=="WH","Whiskery",Species)))))
 
fn.compr=function(SP)
{
  det=subset(Detections,Species==SP & !(TagCode.original%in%only.recaptured))
  tgs=subset(TAGS, Project.rel=="SMN" & Species==SP )
  Cond=subset(Condition,Species==SP,select=c(ATAG,CONDITION))
  
  det=det[!duplicated(det$TagCode.original),]
  tgs=tgs[!duplicated(tgs$Code2),]
  names(Cond)[match("ATAG",names(Cond))]="Code2"
  tgs=merge(tgs,Cond,by="Code2",all.x=T)
  det=merge(det,Cond,by.x="TagCode.original",by.y="Code2",all.x=T)
  det$FL=100*det$FL
  
  #Sex 
  det$Sex=as.character(det$Sex)
  det=subset(det,!(Sex=="U"))
  SX <- as.table(cbind(table(as.character(tgs$Sex)),table(det$Sex)))
  dimnames(SX) <- list(gender = c("F","M"),group = c("Released","Detected"))
  Xsq <- chisq.test(SX)  
  M.F.Sx.ratio=list(rel=c(1,SX[1,1]/SX[2,1]),det=c(1,SX[1,2]/SX[2,2]))
  
  #Condition  
  Xsq.CN=Cond.ratio.1.2.3=CN=NULL
  z=table(as.character(tgs$CONDITION))
  if(length(z)>1)
  {
    a=table(as.character(det$CONDITION))
    S=c("1","2","3")
    if(!length(names(a))==length(S))
    {
      aa=S[which(!S%in%names(a))]
      Add=rep(0,length(aa))
      names(Add)=aa
      a=c(a,Add)
      a=a[sort(names(a))]
    }
    CN <- as.table(cbind(z,a))
    dimnames(CN) <- list(condition = c("1","2","3"),group = c("Released","Detected"))
    Xsq.CN <- chisq.test(CN)  
    Cond.ratio.1.2.3=list(rel=c(CN[1,1]/CN[3,1],CN[2,1]/CN[3,1],1),
                          det=c(CN[1,2]/CN[3,2],CN[2,2]/CN[3,2],1))
    
  }
  

  
  #Size distribution
  Kolmo=ks.test(tgs$Size,det$FL)
  
  return(list(Size.Dif=Kolmo,N.rel.N.det=SX,M.F.Sx.ratio=M.F.Sx.ratio,Xsq.sex=Xsq,
              Cond.rel.det=CN,Cond.ratio.1.2.3=Cond.ratio.1.2.3,Xsq.CN=Xsq.CN))
  
}

Comparo=vector('list',N.sp)
names(Comparo)=SPECIES
for(i in 1:length(SPECIES))Comparo[[i]]=fn.compr(SPECIES[i])

Comparison.table=data.frame(Species=SPECIES,Size=NA,Sex=NA,COND=NA)
for(i in 1:length(SPECIES))
{
  if(is.null(Comparo[[i]]$Xsq.CN$p.value)) Comparo[[i]]$Xsq.CN$p.value=NA
  Comparison.table[i,2:4]=c(Comparo[[i]]$Size.Dif$p.value,Comparo[[i]]$Xsq.sex$p.value,Comparo[[i]]$Xsq.CN$p.value)
}
write.csv(Comparison.table,"Outputs_movement/Size.sex.ratio.condition.csv",row.names=F)



#21. -- Proportion per array ---
Prop.Table=Comparo
Prop.per.array=function(SP)
{
  Det=subset(Detections,!(TagCode.original%in%only.recaptured) & Recapture.hit=="NO" & Species==SP)
  Det$Array=factor(Det$Array,levels=c("Ningaloo","Perth", "Southern.lines"))
  Tabl.hits=table(Det$Array)/sum(table(Det$Array))
  Det$Dummy=with(Det,paste(TagCode,Array))
  Det=Det[!duplicated(Det$Dummy),]
  Tabl.shks=table(Det$Array)/length(unique(Det$TagCode))
  
  return(list(Tabl.hits=Tabl.hits,Tabl.shks=Tabl.shks))
}
for(i in 1:length(SPECIES))Prop.Table[[i]]=Prop.per.array(SPECIES[i])
Prop.Table[[2]]$Tabl.hits=c(.98,.01,.01)  #round for displaying purposes
Tabl.hits=Tabl.shks=NULL
for(i in 1:length(SPECIES))
{
  Tabl.hits=rbind(Tabl.hits,Prop.Table[[i]]$Tabl.hits)
  Tabl.shks=rbind(Tabl.shks,Prop.Table[[i]]$Tabl.shks)
}

colnames(Tabl.hits)=c("Ningaloo","Perth","Southern Lines")
  
tiff(file="Outputs_movement/Proportion.det..Array.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mai=c(.8,1.25,.75,0.2))
barplot(t(Tabl.hits),horiz=T,xlim=c(0,1.025),col=CLS,cex.axis=1.5,
    legend.text=T,args.legend = list(x =1.05, y=5.5, bty = "n",horiz=T,cex=1.75,col=CLS))
box()
mtext("Proportion of detections",1,outer=T,cex=2,line=-1.5)
axis(2,c(.65,1.89,3.03,4.25),c("Dusky","Sandbar","Gummy","Whiskery"),cex.axis=1.5,las=1)
dev.off()



#22. -- Bubble plots of proportion of hits by station for each array ---

#Define location of three lines
  #Ningaloo
plotlong.lines=list(c(113.8,114.04),c(113.54,113.68),c(113.57,113.865))
plotlat.lines=list(c(-22.00,-21.845),c(-22.64,-22.55),c(-23.16,-22.975))

plotlong.lines.seq=list(c(113.90),c(113.35),c(113.70))
plotlat.lines.seq=list(c(-22.0,-21.9),c(-22.5,-22.4),c(-23.1,-23))

TEXTO=rbind(cbind(113.8245,-21.89),cbind(113.5429,-22.56379),cbind(113.59,-23.0))
LETRA=list(as.character(1),as.character(2),as.character(3))



  #Plotting functions               
#note: this doesn't plot SA detections
fn.bubble.rec=function(ESPECIES1,ESPECIES2,ESPECIES3,ESPECIES4,ARRAY,PlotlonG,PlotlatT,
                       scale,LETR,Plus.Lon,Plus.Lat,Plus.Lon1,Plus.Lat1,COLOR,Loc1,Loc2,Leg1,Leg2,Leg3,scale2)
{
  Det=subset(Detections,!(TagCode.original%in%only.recaptured) & Recapture.hit=="NO")
  
  if(What.to.show=="Individuals")
  {
    Det$Dummy=with(Det,paste(TagCode,Species,Station))
    Det=Det[!duplicated(Det$Dummy),]
  }
  
  DATA=subset(Det,Species==ESPECIES1 & Array==ARRAY  &!(Area.release=="SA"))
  DATA2=subset(Det,Species==ESPECIES2 & Array==ARRAY &!(Area.release=="SA"))
  DATA3=subset(Det,Species==ESPECIES3 & Array==ARRAY &!(Area.release=="SA"))
  DATA4=subset(Det,Species==ESPECIES4 & Array==ARRAY &!(Area.release=="SA"))
  
  TABLE=table(DATA$Station)
  TABLE2=table(DATA2$Station)
  TABLE3=table(DATA3$Station)
  TABLE4=table(DATA4$Station)
  
  if(!What.to.show=="Individuals")
  {
    N=sum(TABLE)
    N2=sum(TABLE2)
    N3=sum(TABLE3)
    N4=sum(TABLE4)
  }
  if(What.to.show=="Individuals")
  {
    N=length(unique(DATA$TagCode))
    N2=length(unique(DATA2$TagCode))
    N3=length(unique(DATA3$TagCode))
    N4=length(unique(DATA4$TagCode))
  }
  
  MATRIX=matrix(names(TABLE))
  MATRIX=matrix(as.numeric(unlist(strsplit(MATRIX, split=" "))),ncol=2,byrow=T)
  LAT=MATRIX[,1]
  LONG=MATRIX[,2]
  zo=TABLE/sum(TABLE,na.rm=T)
  if(What.to.show=="Individuals") zo=TABLE/N
  
  if(nrow(TABLE2)>0)
  {
    MATRIX2=matrix(names(TABLE2))
    MATRIX2=matrix(as.numeric(unlist(strsplit(MATRIX2, split=" "))),ncol=2,byrow=T)
    LAT2=MATRIX2[,1]
    LONG2=MATRIX2[,2]
    zo2=TABLE2/sum(TABLE2,na.rm=T)
    if(What.to.show=="Individuals") zo2=TABLE2/N2   
  }
  
  if(nrow(TABLE3)>0)
  {
    MATRIX3=matrix(names(TABLE3))
    MATRIX3=matrix(as.numeric(unlist(strsplit(MATRIX3, split=" "))),ncol=2,byrow=T)
    LAT3=MATRIX3[,1]
    LONG3=MATRIX3[,2]
    zo3=TABLE3/sum(TABLE3,na.rm=T)
    if(What.to.show=="Individuals") zo3=TABLE3/N3
  }
  
  if(nrow(TABLE4)>0)
  {
    MATRIX4=matrix(names(TABLE4))
    MATRIX4=matrix(as.numeric(unlist(strsplit(MATRIX4, split=" "))),ncol=2,byrow=T)
    LAT4=MATRIX4[,1]
    LONG4=MATRIX4[,2]
    zo4=TABLE4/sum(TABLE4,na.rm=T)
    if(What.to.show=="Individuals") zo4=TABLE4/N4
  }
  
  plotMap(worldLLhigh, xlim=PlotlonG,ylim=PlotlatT,plt = c(.001, 1, 0.0775, 1),
          col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=PlotlatT,xlim=PlotlonG, zlim=c(-1,-100),
          nlevels = 2,labcex=1.25,lty = 1,col=c(COLOR,COLOR,"transparent"),add=T)  
  
  #add receivers
  points(STATIONS$longitude,STATIONS$latitude,col=1,pch=19)
  box(lwd=2)
  
  #add hits
  if(nrow(TABLE)>0)points(LONG-Plus.Lon,LAT*(1-Plus.Lat),pch=21,cex=zo*scale,col=COL.Sp[1],bg=COL.Sp.bg[1],lwd=2)
  if(nrow(TABLE2)>0)points(LONG2+Plus.Lon, LAT2*(1+Plus.Lat),pch=21,cex=zo2*scale,col=COL.Sp[2],bg=COL.Sp.bg[2],lwd=2)
  if(nrow(TABLE3)>0)
  {
    if(ARRAY=="Perth")
    {
      LO=ifelse(LONG3<115.687,LONG3-Plus.Lon1,LONG3)
      LA=ifelse(LONG3<115.687,LAT3*(1-Plus.Lat1),LAT3)
      points(LO, LA,pch=21,cex=zo3*scale,col=COL.Sp[3],bg=COL.Sp.bg[3],lwd=2)      
    }else
    {
      points(LONG3-Plus.Lon1, LAT3*(1-Plus.Lat1),pch=21,cex=zo3*scale,col=COL.Sp[3],bg=COL.Sp.bg[3],lwd=2)
    }

    
  }
  if(nrow(TABLE4)>0)points(LONG4+Plus.Lon1, LAT4*(1+Plus.Lat1),pch=21,cex=zo4*scale,col=COL.Sp[4],bg=COL.Sp.bg[4],lwd=2)
  if(!(ARRAY=="Perth"))legend("topright",LETR,cex=3,bty='n')
  
  
  if(LETR=="2")
  {   
    if(!What.to.show=="Individuals") 
    {
      legend("bottomleft",legend=c(paste("dusky (n= ",N," detections)",sep=""),paste("sandbar (n= ",N2," detections)",sep=""),
                                   paste("gummy (n= ",N3," detections)",sep=""),paste("whiskery (n= ",N4," detections)",sep="")),
             pch=21,col=COL.Sp,pt.bg=COL.Sp.bg ,pt.lwd=2,bty='n',cex=1.25,pt.cex=2)
    }
    if(What.to.show=="Individuals") 
    {
      legend("bottomleft",legend=c(paste("dusky (n= ",N," sharks)",sep=""),paste("sandbar (n= ",N2," sharks)",sep=""),
                                   paste("gummy (n= ",N3," sharks)",sep=""),paste("whiskery (n= ",N4," sharks)",sep="")),
             pch=21,col=COL.Sp,pt.bg=COL.Sp.bg ,pt.lwd=2,bty='n',cex=1.25,pt.cex=2)
    }
    
  }
  
  if(LETR=="1" & ARRAY=="Ningaloo")
  {
    points(Loc1,Loc2,pch=21,cex=Leg1*scale2)
    points(Loc1,Loc2*1.000183,pch=21,cex=Leg2*scale2)
    points(Loc1,Loc2*1.000366,pch=21,cex=Leg3*scale2)
    boxed.labels(Loc1,Loc2*0.9992558,paste(Leg1*100,"%",sep=""), cex = .75,bg="white",border=NA)
    boxed.labels(Loc1,Loc2*0.9996425,paste(Leg2*100,"%",sep=""), cex = .75,bg="white",border=NA)
    boxed.labels(Loc1,Loc2*1.00001,paste(Leg3*100,"%",sep=""), cex = .75,bg="white",border=NA)
    
  }
  
  
  if(LETR=="3" & ARRAY=="Southern.lines")
  {
    points(Loc1,Loc2,pch=21,cex=Leg1*scale2)
    points(Loc1,Loc2*1.0003,pch=21,cex=Leg2*scale2)
    if(!What.to.show=="Individuals")
    {
      points(Loc1,Loc2*1.000575,pch=21,cex=Leg3*scale2)
      boxed.labels(Loc1,Loc2*0.9992,paste(Leg1*100,"%",sep=""), cex = .625,bg="white",border=NA)
      boxed.labels(Loc1,Loc2*0.99975,paste(Leg2*100,"%",sep=""), cex = .625,bg="white",border=NA)
      boxed.labels(Loc1,Loc2*1.0003,paste(Leg3*100,"%",sep=""), cex = .625,bg="white",border=NA)
      
    }
    
    if(What.to.show=="Individuals")
    {
      points(Loc1,Loc2*1.00068,pch=21,cex=Leg3*scale2)
      boxed.labels(Loc1,Loc2*0.9984,paste(Leg1*100,"%",sep=""), cex = .625,bg="white",border=NA)
      boxed.labels(Loc1,Loc2*0.9993,paste(Leg2*100,"%",sep=""), cex = .625,bg="white",border=NA)
      boxed.labels(Loc1,Loc2*1.0000,paste(Leg3*100,"%",sep=""), cex = .625,bg="white",border=NA)
      
    }
    
  }
  
  
  if(ARRAY=="Perth")
  {
    points(Loc1,Loc2,pch=21,cex=Leg1*scale)
    points(Loc1,Loc2*1.00043,pch=21,cex=Leg2*scale)
    points(Loc1,Loc2*1.000866,pch=21,cex=Leg3*scale)
    
    boxed.labels(Loc1,Loc2*0.99875,paste(Leg1*100,"%",sep=""), cex = .7,bg="white",border=NA)
    boxed.labels(Loc1,Loc2*0.9996,paste(Leg2*100,"%",sep=""), cex = .7,bg="white",border=NA)
    boxed.labels(Loc1,Loc2*1.00035,paste(Leg3*100,"%",sep=""), cex = .7,bg="white",border=NA)
    
  }
  
  
}

fn.main.map=function(XLIM,YLIM,LONGSEQ,LATSEQ,Xtext,Ytext,TEXT,TEXT.cx)
{
  plotMap(worldLLhigh, xlim=XLIM,ylim=YLIM,plt = c(.001, 1, 0.075, 1),
          col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  points(STATIONS$longitude,STATIONS$latitude,col=1,pch=19)
  axis(side = 1, at =LONGSEQ, labels = LONGSEQ, tcl = .5,las=1,cex.axis=1.5,padj=-.5)
  axis(side = 2, at = LATSEQ, labels = -LATSEQ,tcl = .5,las=2,cex.axis=1.5,hadj=.75)
  box(lwd=2)
  mtext("Latitude (ºS)",side=4,line=1,las=3,cex=1.5)
  mtext("Longitude (ºE)",side=1,line=1.85,cex=1.5)
  contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=YLIM,xlim=XLIM, zlim=c(-1,-200),
          nlevels = 4,labcex=1.,lty = 1,col=c(COLOR,COLOR,COLOR,COLOR,"transparent"),add=T)
  text(Xtext,Ytext,TEXT,cex=TEXT.cx)
}

def.par <- par(no.readonly = TRUE)  #to reset par to default
COLOR="dark grey"
COL.Sp=c("mediumpurple4","orange","green","brown")
COL.Sp.bg=c("#0000ff22","#8B000022","#0000ff39","#8B000032")

if(do.map=="YES")
{
  
    #Plot hits
  What.to.show="Detections"
  
  #_ Ningaloo
  tiff(file=paste("Outputs_movement/Figure4.Ningaloo.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mar = c(3, 3, 0, 0))
  m <- cbind(c(2, 3,4),c(1, 1,1))
  layout(m)
  #layout.show(4)
  fn.main.map(XLIM=plotlong[[1]],YLIM=c(-23.25,-21.75),LONGSEQ=Long.seq[[1]],LATSEQ=Lat.seq[[1]],
              Xtext=rbind(114,113.7,113.85),Ytext=rbind(-21.9,-22.6,-23.1),TEXT=c(1:3),2.5)
  text(114.128,-21.93,"Exmouth",cex=1.75)
  for (j in 1:3)
  {
    scale=50
    fn.bubble.rec(SPECIES[1],SPECIES[2],SPECIES[3],SPECIES[4],"Ningaloo",plotlong.lines[[j]],
      plotlat.lines[[j]],scale,LETRA[[j]],0,6e-04,0,8e-04,COLOR,Loc1=113.82,Loc2=-21.875,
      Leg1=.2,Leg2=.15,Leg3=.1,scale)
  } 
  dev.off()
  par(def.par)
  
  #_ Perth
  scale=10
  par(mar = c(2, 5, 3, 0.5),oma=c(2,2,2,1))
  tiff(file=paste("Outputs_movement/Figure4.Perth.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  fn.bubble.rec(SPECIES[1],SPECIES[2],SPECIES[3],SPECIES[4],"Perth",c(115,116),
    c(-32.5,-31.5),scale,"2",0,6e-04,0,8e-04,COLOR,Loc1=115.0911,Loc2=-31.65178,Leg1=.75,Leg2=.5,Leg3=.25,scale)
  polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
  polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
  text(115.82,-31.92,"Perth",cex=1.75)
   dev.off()
  
  
  #_ Southern Lines
  par(def.par)
  tiff(file=paste("Outputs_movement/Figure4.South.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mar = c(2, 5, 1, 0.5),oma=c(2,2,1,1))
  m <- rbind(c(1,1,1),c(2,3,4))
  layout(m)
  #layout.show(4)
  fn.main.map(XLIM=c(114.5,119.5),YLIM=c(-35.5,-33.5),LONGSEQ=seq(115,119,by=.5),
              LATSEQ=seq(-35.5,-33.5,by=.5),Xtext=rbind(115.1,116.5,118.32),
              Ytext=rbind(-34.2,-34.75,-34.78),TEXT=c("1","2","3"),2.5)
 # text(115.15,-34.36,"Augusta",cex=1.75)
 #  text(117.88,-35.01,"Albany",cex=1.75)
 text(115.125,-34.15,"(Cape Leeuwin)",cex=2.5,srt=35,pos=4)
 text(116.5,-34.65,"(Chatham island)",cex=2.5,srt=45,pos=4)
 text(118.35,-34.66,"(Bald island)",cex=2.5,srt=45,pos=4)
  long.line=rbind(c(114.58,115.11),c(116.25,116.75),c(118.25,118.75))
  lat.line=rbind(c(-34.41,-34.18),c(-35.5,-35),c(-35.35,-34.85))
  Plus.LONG=c(0,2.5e-02,1e-02)
  Plus.LAT=c(6e-04,0,0)
  Plus.LONG1=c(0,4e-02,5e-02)
  Plus.LAT1=c(1.5e-03,0,0)
  for (j in 1:3)
  {
    scale=30
    fn.bubble.rec(SPECIES[1],SPECIES[2],SPECIES[3],SPECIES[4],"Southern.lines",long.line[j,],
                  lat.line[j,],scale,LETRA[[j]],Plus.LONG[j],Plus.LAT[j],Plus.LONG1[j],
                  Plus.LAT1[j],COLOR,Loc1=118.6889,Loc2=-35.185,Leg1=.3,Leg2=.2,Leg3=.1,scale)
  } 
  dev.off()
  par(def.par)
  
  
  
  #Plot Individuals
  What.to.show="Individuals"
  
  #_ Ningaloo
  tiff(file=paste("Outputs_movement/Figure4.Ningaloo.individuals.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mar = c(3, 3, 0, 0))
  m <- cbind(c(2, 3,4),c(1, 1,1))
  layout(m)
  #layout.show(4)
  fn.main.map(XLIM=plotlong[[1]],YLIM=c(-23.25,-21.75),LONGSEQ=Long.seq[[1]],LATSEQ=Lat.seq[[1]],
              Xtext=rbind(114,113.7,113.85),Ytext=rbind(-21.9,-22.6,-23.1),TEXT=c(1:3))
  text(114.128,-21.93,"Exmouth",cex=1.75)
  for (j in 1:3)
  {
    fn.bubble.rec(SPECIES[1],SPECIES[2],SPECIES[3],SPECIES[4],"Ningaloo",plotlong.lines[[j]],
                  plotlat.lines[[j]],10,LETRA[[j]],0,6e-04,0,8e-04,COLOR,Loc1=113.82,Loc2=-21.875,
                  Leg1=1,Leg2=.75,Leg3=.5,10)
  } 
  dev.off()
  
  
  #_ Perth
  par(def.par)
  tiff(file=paste("Outputs_movement/Figure4.Perth.individuals.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  fn.bubble.rec(SPECIES[1],SPECIES[2],SPECIES[3],SPECIES[4],"Perth",c(115,116),
                c(-32.5,-31.5),10,"2",0,6e-04,0,8e-04,COLOR,Loc1=115.0911,Loc2=-31.65178,Leg1=.75,Leg2=.5,Leg3=.25,10)
  polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
  polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
  text(115.86,-31.95,"Perth",cex=1.75)
  
  dev.off()
  
  
  #_ Southern Lines
  par(def.par)
  tiff(file=paste("Outputs_movement/Figure4.South.individuals.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mar = c(2, 3, 1, 0.5))
  m <- rbind(c(1, 1,1),c(2, 3,4))
  layout(m)
  #layout.show(4)
  fn.main.map(XLIM=c(114.5,119.5),YLIM=c(-35.5,-33.5),LONGSEQ=seq(115,119,by=.5),
              LATSEQ=seq(-35.5,-33.5,by=.5),Xtext=rbind(115.1,116.5,118.32),
              Ytext=rbind(-34.2,-34.75,-34.78),TEXT=c(1:3))
  text(115.15,-34.36,"Augusta",cex=1.75)
  text(117.88,-35.01,"Albany",cex=1.75)
  long.line=rbind(c(114.58,115.11),c(116.25,116.75),c(118.25,118.75))
  lat.line=rbind(c(-34.41,-34.18),c(-35.5,-35),c(-35.35,-34.85))
  Plus.LONG=c(0,2.5e-02,1e-02)
  Plus.LAT=c(6e-04,0,0)
  Plus.LONG1=c(0,4e-02,5e-02)
  Plus.LAT1=c(1.5e-03,0,0)
  for (j in 1:3)
  {
    fn.bubble.rec(SPECIES[1],SPECIES[2],SPECIES[3],SPECIES[4],"Southern.lines",long.line[j,],
                  lat.line[j,],15,LETRA[[j]],Plus.LONG[j],Plus.LAT[j],Plus.LONG1[j],
                  Plus.LAT1[j],COLOR,Loc1=118.6889,Loc2=-35.185,Leg1=1,Leg2=.75,Leg3=.5,15)
  } 
  dev.off()
  par(def.par)
  
}


#23. -- Daily patterns ---
#note: drop Perth and south for sandbar and Perth for gummy for displaying purposes
Daily.fn=function(SP)
{
  Det=subset(Detections,!(TagCode.original%in%only.recaptured) & Recapture.hit=="NO" & Species==SP)
  if(SP=="Thickskin") Det=subset(Det,Array%in%"Ningaloo")
  if(SP=="Gummy") Det=subset(Det,Array%in%"Southern.lines")
  Det$Array=factor(Det$Array,levels=c("Ningaloo","Perth", "Southern.lines"))
  if(class(Det$Time.local)=="times") Det$hour=substr(Det$Time.local,1,2)
  if(!class(Det$Time.local)=="times")Det$hour=hours(Det$Time.local)
  Tabl.hits=table(Det$hour,Det$Array)
  SUM=colSums(Tabl.hits)
  for(j in 1:ncol(Tabl.hits)) if(SUM[j]>0) Tabl.hits[,j]=Tabl.hits[,j]/SUM[j]
  return(Tabl.hits)
}
Prop.daily=Prop.Table
for(i in 1:length(SPECIES))Prop.daily[[i]]=Daily.fn(SPECIES[i])

fn.day=function(dat)
{
  plot(dat[,3],type='l',lwd=2,col='transparent',ylim=c(0,max(dat)*1.01),ylab="",
       xlab="",xaxt='n',cex.axis=1.25)
  if(sum(dat[,2])>0)lines(dat[,2],lwd=2,col=CLS[2])
  if(sum(dat[,1])>0) lines(dat[,1],lwd=2,col=CLS[1])
  if(sum(dat[,3])>0)lines(dat[,3],lwd=2,col=CLS[3])
}

tiff(file="Outputs_movement/Daily.Array.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mai=c(.4,.1,.2,.5),oma=c(2,4,2,.1),mgp=c(1,.8,0),las=1)
NMS=c("Dusky","Sandbar" ,"Gummy","Whiskery")
for(i in 1:length(SPECIES))
{
  fn.day(Prop.daily[[i]])
  mtext(NMS[i],3,cex=1.25)
  if(i==1) legend('bottomleft',c("Ningaloo","Perth", "Southern Lines"),col=CLS,
                  lty=1,lwd=3,bty='n',cex=1.5)
  
  axis(1,0:23,F,tcl = -.25)
  axis(1,seq(0,23,by=3),F,tcl = -.5)
  if(i%in%c(2,4)) axis(1,seq(0,23,by=3),seq(0,23,by=3),cex.axis=1.25,tcl = -.5)
}
mtext("Proportion of detections",2,outer=T,line=2.5,las=3,cex=1.5)
mtext("Hour of day",1,outer=T,line=0,las=1,cex=1.5)
dev.off()



#24. -- Co-detection of individuals within same hour---
Co.detect=function(SP,intraspecies,where)
{
  if(!intraspecies=="YES")
  {
    Det=subset(Detections,!(TagCode.original%in%only.recaptured) & Recapture.hit=="NO" & Species%in%SP)
    if(class(Det$Time.local)=="times") Det$hour=substr(Det$Time.local,1,2)
    if(!class(Det$Time.local)=="times")Det$hour=hours(Det$Time.local)
    Det$Same=with(Det,paste(Station,Date.local,hour))    
    tabl=with(Det,table(Same,TagCode))
    tabl[tabl>0]=1
    a=rowSums(tabl)
  }
  
  if(intraspecies=="YES")
  {
    Det=subset(Detections,!(TagCode.original%in%only.recaptured) & 
                 Recapture.hit=="NO" & Species%in%SP & Array==where) 
    if(class(Det$Time.local)=="times") Det$hour=substr(Det$Time.local,1,2)
    if(!class(Det$Time.local)=="times")Det$hour=hours(Det$Time.local)
    Det$Same=with(Det,paste(Station,Date.local,hour))    
    Det$Same=as.factor(Det$Same)
    Det$TagCode=as.factor(Det$TagCode)
    
    tabl=with(Det,table(Same,TagCode,Species))
    b=tabl[,,1]
    b[b>0]=1
    b1=rowSums(b)
 
    b=tabl[,,2]
    b[b>0]=1
    b2=rowSums(b)
     
    a=b1+b2
    
  }
  
  return(a)
}

  #within species
Co.detec.same=Prop.daily
for(i in 1:length(SPECIES))Co.detec.same[[i]]=Co.detect(SPECIES[i],'NO','all')

  #between pairs of species
Pair.comp=list(c("Dusky","Thickskin"),c("Dusky","Gummy"),c("Dusky","Whiskery"),c("Gummy","Whiskery"))
Pair.where=c('Ningaloo','Southern.lines','Southern.lines','Southern.lines')
Co.detec.diff=vector('list',length(Pair.comp))
names(Co.detec.diff)=Pair.comp
for(i in 1:length(Pair.comp))Co.detec.diff[[i]]=Co.detect(Pair.comp[[i]],'YES',Pair.where[i])

fn.see.co.det=function(a)
{
  barplot(table(a),col=Singl.col,cex.axis=1.5,cex.lab=1.5,cex.names=1.5)
  box()
}

tiff(file="Outputs_movement/Co.detections_same.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mai=c(.4,.3,.2,.5),oma=c(2,4,2,.1),mgp=c(1,.6,0),las=1)
for(i in 1:length(SPECIES))
{
  fn.see.co.det(Co.detec.same[[i]])  
  mtext(NMS[i],3,cex=1.3)
}
mtext("Number of individuals detected together",1,outer=T,line=0,las=1,cex=1.75)
mtext("Frequency",2,outer=T,line=2.5,las=3,cex=1.75)
dev.off()

NMS.co=list(c("Dusky and Sandbar (Ningaloo)"),c("Dusky and Gummy (Southern lines)"),
            c("Dusky and Whiskery (Southern lines)"),c("Gummy and Whiskery (Southern lines)"))
tiff(file="Outputs_movement/Co.detections_different.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mai=c(.4,.3,.2,.5),oma=c(2,4,2,.1),mgp=c(1,.6,0),las=1)
for(i in 1:length(SPECIES))
{
  fn.see.co.det(Co.detec.diff[[i]])  
  mtext(NMS.co[[i]],3,cex=1.3)
}
mtext("Number of individuals detected together",1,outer=T,line=0,las=1,cex=1.75)
mtext("Frequency",2,outer=T,line=2.5,las=3,cex=1.75)
dev.off()




#25. --  Speed and distance travelled ---

    # Create data sets

  #check release and recapture area by species
Table.diff.Area=with(Detections,table(Area,Area.release,Species))

#create subset of data for individuals moving more than the migration threshold (Migration.threshold)
MIGRATE.list=MIGRATE.tags=vector('list',length(SPECIES))
names(MIGRATE.list)=names(MIGRATE.tags)=SPECIES
Fun.mig.dat=function(Thresh,SPEC)
  {
    dat=subset(Detections,Dist.moved.conseq.det>Thresh & Species ==SPEC)    
    unikTag=unique(dat$TagCode)
    dat=dat[order(dat$TagCode, dat$Date.local),]
    dat1=NULL
    for (p in 1:length(unikTag))
    {
      dd=subset(dat,TagCode==unikTag[p])
      conseq.det=subset(dd,Dist.moved.conseq.det>Thresh)
      conseq.det$TYPE="Conseq.det"
      thiss=conseq.det[!duplicated(conseq.det$DateTime.local.prev),]
      dat1=rbind(dat1,thiss)
    }
     return(dat1)
  }

for(i in 1:length(MIGRATE.list))
  {
    a=Fun.mig.dat(Migration.threshold[i],SPECIES[i])
    MIGRATE.list[[i]]=a
    MIGRATE.tags[[i]]=unique(a$TagCode)
  }


    #Long displacements summary tables 
Summary.mig=function(Dat)
{
  #Consecutive detections is unique so keep all
  Dat$Km.per.day=Dat$Dist.moved.conseq.det/Dat$days.conseq.det
  ID=c("TagCode","Dist.moved.conseq.det","days.conseq.det","Km.per.day")
  MaxDist=Dat[order(Dat$TagCode,Dat$Dist.moved.conseq.det),ID]

  Mean.speed=round(mean(MaxDist$Km.per.day,na.rm=T),2)
  Max.speed=max(MaxDist$Km.per.day,na.rm=T)
  
  Mean.dist=round(mean(MaxDist$Dist.moved.conseq.det,na.rm=T))
  Max.dist=max(MaxDist$Dist.moved.conseq.det,na.rm=T)
  
  MaxDist[,2:3]=round(MaxDist[,2:3])
  MaxDist[,4]=round(MaxDist[,4],2)
  
  N.detections=nrow(MaxDist)
  N.sharks=length(unique(MaxDist$TagCode))
   
  return(list(MaxDist=MaxDist,Stats=rbind(N.detections,N.sharks,Mean.dist,Max.dist,Mean.speed,Max.speed)))
}
# Summary.mig=function(Dat)
# {
#   #Consecutive detections is unique so keep all
#   Dat$Km.per.day=Dat$Dist.moved.conseq.det/Dat$days.conseq.det
#   
#   Max.distance=aggregate(Dist.moved.conseq.det~TagCode,Dat,max) 
#   ID.max=match(Max.distance$Dist.moved.conseq.det,Dat$Dist.moved.conseq.det)
#   MaxDist=Dat[ID.max,c("TagCode","Dist.moved.conseq.det","days.conseq.det","Km.per.day")]
#   
#   Max.speed=aggregate(Km.per.day~TagCode,Dat,max)
#   ID.max=match(Max.speed$Km.per.day,Dat$Km.per.day)
#   MaxSpeed=Dat[ID.max,c("TagCode","Dist.moved.conseq.det","days.conseq.det","Km.per.day")]
#   
#   return(list(MaxDist=MaxDist,MaxSpeed=MaxSpeed))
# }

store.MAXDIST=MIGRATE.list
for (i in 1:length(MIGRATE.list))
{
  MAXDIST=Summary.mig(MIGRATE.list[[i]])$MaxDist
  #MAXSPID=Summary.mig(MIGRATE.list[[i]])$MaxSpeed
  Thress=Migration.threshold[i]
  
  tiff(file=paste("Outputs_movement/Long.scale.displacement/MaxDispSpeed.",SPECIES[i],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(1,1),mar=c(2,2,1,1), oma=c(1,1.5,1,1),las=1)
  plot(MAXDIST$Dist.moved.conseq.det,MAXDIST$Km.per.day,pch=19,col=Singl.col,cex=1.5,ylab="",xlab="",
       ylim=c(0,max(MAXDIST$Km.per.day)))
  legend("topright",c("Maximum displacement",paste("(",Thress, "km threshold)")),bty="n",cex=1.15)
#   plot(MAXSPID$Dist.moved.conseq.det,MAXSPID$Km.per.day,pch=19,col=2,cex=1.5,ylab="",xlab="",
#        ylim=c(0,max(MAXSPID$Km.per.day)))
#   legend("topright","Maximum speed",bty="n",cex=1.15)
   mtext("Displacement (km)",1,line=0,outer=T,cex=1.5)
 # mtext("Speed (km/day)",2,line=0,outer=T,cex=1.5,las=3)
  dev.off()
  
  #max distance
  write.table(Summary.mig(MIGRATE.list[[i]])$MaxDist,file=
  paste("Outputs_movement/Long.scale.displacement/","Cons.Det.Max.dist(km)",SPECIES[i],".csv",sep=""),
              sep = ",", col.names = T,row.names = F,qmethod = "double")
   store.MAXDIST[[i]]=MAXDIST
}

#export table of stats
n=nrow(Summary.mig(MIGRATE.list[[i]])$Stats)
Sum=matrix(nrow=n,ncol=4)
for(i in 1:ncol(Sum)) Sum[,i]=Summary.mig(MIGRATE.list[[i]])$Stats
rownames(Sum)=rownames(Summary.mig(MIGRATE.list[[i]])$Stats)
colnames(Sum)=names(MIGRATE.list)

#stats
write.table(Sum,file="Outputs_movement/Long.scale.displacement/Summary.of.long-.distance.displacements.csv",
            sep = ",", col.names = T,row.names = T,qmethod = "double")


#compare the 4 species maximum consecutive detections
COLss=1:4
Xmax=max(c(store.MAXDIST[[1]]$days.conseq.det,store.MAXDIST[[2]]$days.conseq.det,
           store.MAXDIST[[3]]$days.conseq.det,store.MAXDIST[[4]]$days.conseq.det))
Ymax=max(c(store.MAXDIST[[1]]$Dist.moved.conseq.det,store.MAXDIST[[2]]$Dist.moved.conseq.det,
           store.MAXDIST[[3]]$Dist.moved.conseq.det,store.MAXDIST[[4]]$Dist.moved.conseq.det))

tiff(file="Outputs_movement/Long.scale.displacement/Max.dist.moved.Consec.detec.tiff",width = 2400, height = 2200,units = "px", res = 300,compression = "lzw")
par(mfcol=c(1,1),mar=c(2,3,1,1), oma=c(1,1.5,1,1),las=1,mgp=c(1,.6,0))
with(store.MAXDIST[[1]],plot(days.conseq.det,Dist.moved.conseq.det,pch=19,col=COLss[1]
                             ,ylim=c(0,Ymax),xlim=c(0,Xmax),cex=1.5,
                             ylab="",xlab=""))
for(i in 2:4) with(store.MAXDIST[[i]],points(days.conseq.det,Dist.moved.conseq.det,pch=19,col=COLss[i],cex=1.5))
legend("right",SPECIES,pch=19,col=COLss,bty='n',pt.cex=1.5,cex=1.5)
mtext("Max. distance moved between consec. detections (km)",2,line=0,outer=T,cex=1.5,las=3)
mtext("Days between consec. detections",1,line=0,outer=T,cex=1.5)
dev.off()


#Compare the 4 species maximum speed
Xmax=max(c(store.MAXDIST[[1]]$days.conseq.det,store.MAXDIST[[2]]$days.conseq.det,
           store.MAXDIST[[3]]$days.conseq.det,store.MAXDIST[[4]]$days.conseq.det))
Ymax=max(c(store.MAXDIST[[1]]$Km.per.day,store.MAXDIST[[2]]$Km.per.day,
           store.MAXDIST[[3]]$Km.per.day,store.MAXDIST[[4]]$Km.per.day))

# tiff(file="Outputs_movement/Long.scale.displacement/Max.Km.per.day.Consec.detec.tiff",width = 2400, height = 2200,units = "px", res = 300,compression = "lzw")
# par(mfcol=c(1,1),mar=c(2,3,1,1), oma=c(1,1.5,1,1),las=1,mgp=c(1,.6,0))
# with(store.MAXDIST[[1]],plot(days.conseq.det,Km.per.day,pch=19,col=COLss[1]
#                              ,ylim=c(0,Ymax),xlim=c(0,Xmax),cex=1.5,
#                              ylab="",xlab=""))
# for(i in 2:4) with(store.MAXDIST[[i]],points(days.conseq.det,Km.per.day,pch=19,col=COLss[i],cex=1.5))
# legend("right",SPECIES,pch=19,col=COLss,bty='n',pt.cex=1.5,cex=1.5)
# mtext("Max. speed between consec. detections (km/day)",2,line=0,outer=T,cex=1.5,las=3)
# mtext("Days between consec. detections",1,line=0,outer=T,cex=1.5)
# dev.off()


#Speed boxplot for different observation windows
Dist.ranges=c(1,10,50,100,10000)

Speed.boxplot=function(DAT,range1,LAB)
{
  DAT=DAT[order(DAT$TagCode,DAT$DateTime.local),]
  DAT=subset(DAT,Rec.hours.conseq.det=="diff" & hours.conseq.det>minHours)    
  a=DAT
  a$KM.hour=ifelse(a$hours.conseq.det>0,a$Dist.moved.conseq.det/a$hours.conseq.det,NA)                 
  a=subset(a,KM.hour>0,select=c(TagCode,KM.hour,Dist.moved.conseq.det))    
  
  
  names(a)[match("Dist.moved.conseq.det",names(a))]=c("Distance.km")
  a=subset(a,Distance.km>range1[1])
  a$Dis.range=cut(a$Distance.km,range1,include.lowest=T,right=F)
  
  boxplot(KM.hour~Dis.range,a,xaxt='n',col="grey70")
  legend("top",LAB,bty="n",cex=1.75)
}

tiff(file="Outputs_movement/Speed.boxplot.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,1),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.75, 0))
for (s in 1:length(MIGRATE.list))
{
  Speed.boxplot(subset(Detections,Species==SPECIES[s]),Dist.ranges,SPEC.nms[s])
}
mtext("Displacement category (km)",1,outer=T,cex=1.25)
mtext("Minimum speed (km/hour)",2,outer=T,las=3,cex=1.25,line=0.5)
dev.off()




    # Plot migration
#note: this is incomplete, too hard to see for duskies, too many points. Just present table
fun.plot.mig=function(DATA,tags,LONG,LAT)
{
  plotMap(worldLLhigh, xlim=LONG,ylim=LAT,plt = c(.001, 1, 0.075, 1),
          col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  axis(side = 1, at =seq(LONG[1],LONG[2]), labels = seq(LONG[1],LONG[2]), tcl = .5,las=1,cex.axis=1.5,padj=-.5)
  axis(side = 2, at = seq(round(LAT[1]),LAT[2]-1), labels = -(seq(round(LAT[1]),LAT[2]-1)),tcl = .5,las=2,cex.axis=1.5,hadj=.75)
  box(lwd=2)
  mtext("Latitude (ºS)",side=2,line=2.5,las=3,cex=1.5)
  mtext("Longitude (ºE)",side=1,line=1.55,cex=1.5)
  contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=LONG,xlim=LAT, zlim=c(-1,-200),
          nlevels = 4,labcex=1.,lty = 1,col=c(COLOR,COLOR,COLOR,COLOR,"transparent"),add=T)
  
  for(t in 1:length(tags)) 
  {
    Data=subset(DATA,TagCode==tags[t])
    
    Col=col.mig[t]
    
    ###
    LTY=1:nrow(Data)
    LWD=seq(1,3,length.out=nrow(Data))
    
    jit=seq(0,.005,length.out=nrow(Data))
    for (a in 1:nrow(Data))
    {
      if(a==1)
      {
        LonP=jitter(Data$ReleaseLongitude[1],jit[a])
        LatP=jitter(Data$ReleaseLatitude[1],jit[a])
        Lon=jitter(Data$Longitude[1],jit[a])
        Lat=jitter(Data$Latitude[1],jit[a])
      }
      
      if(a>1)
      {
        LonP=jitter(Data$Longitude.prev[a],jit[a])
        LatP=jitter(Data$Latitude.prev[a],jit[a])
        Lon=jitter(Data$Longitude[a],jit[a])
        Lat=jitter(Data$Latitude[a],jit[a])
      }
      
      
      if(LonP<=115 & LatP>(-34.3))
      {
        if(LatP>Lat)
        {
          arrows(LonP,LatP,LonP,LatP-.5,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lon<=115)arrows(Lon,Lat+.5,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lon>115&Lat>(-34.3))arrows(Lon,Lat+.5,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lon>115&Lat<(-34.3))arrows(Lon-.5,Lat,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
        }
        
        if(LatP<Lat)
        {
          if(Lon<=115)
          {
            arrows(LonP,LatP,LonP,LatP+.5,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
            arrows(Lon,Lat-.5,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          }
          if(Lon>115)
          {
            arrows(LonP,LatP,LonP+.5,LatP,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
            arrows(Lon-.5,Lat,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          }
        }
      }
      
      if(LonP>115& LatP<=(-34.3))
      {
        if(LonP>Lon)
        {
          arrows(LonP,LatP,LonP-.5,LatP,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lat<=(-34.3))arrows(Lon,Lat,Lon-.5,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lat>-(34.3))arrows(Lon,Lat-.5,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
        }
        
        if(LonP<Lon)
        {
          arrows(LonP,LatP,LonP+.5,LatP,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          arrows(Lon-.5,Lat,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
        }
      }
      
      if(LonP>115& LatP>(-34.3))
      {
        if(LonP>Lon)
        {
          if(Lat>(-22))arrows(LonP,LatP,LonP-.5,LatP-.5,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lat<(-22))arrows(LonP,LatP,LonP,LatP-.5,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lat<=(-34.3))arrows(Lon-.5,Lat,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lat>(-34.3))arrows(Lon,Lat+.5,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          
        }
        if(LonP<Lon)
        {
          if(Lat>LatP)arrows(LonP,LatP,LonP+.5,LatP+.5,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lat<LatP)arrows(LonP+.5,LatP+.5,LonP,LatP,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
          if(Lat<(-34.3))arrows(Lon-.5,Lat,Lon,Lat,col=Col,length=.15,lwd=LWD[a],lty=LTY[a])
        }
      }
      
    }
    ###
    
    points(Data$ReleaseLongitude,Data$ReleaseLatitude,col=Col,pch=19,cex=1.35)
  }
  
}


# 
# for (i in 1:length(MIGRATE.list))
# {
#   par(mar = c(3.5, 3, 0, 0))
#   tiff(file=paste("Outputs_movement/Migration.",SPECIES[i],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
#   col.mig=rainbow(length(MIGRATE.tags[[i]]))
#   fun.plot.mig(MIGRATE.list[[i]],MIGRATE.tags[[i]],LONG=c(112,125),LAT=c(-35.5,-21))
#   dev.off()
# }



Speed.dist.hist=function(DAT,range1,plotWhat,MAX1,MAX2,MAX3,r)
{
  DAT=DAT[order(DAT$TagCode,DAT$DateTime.local),]
  DAT=subset(DAT,Rec.hours.conseq.det=="diff" & hours.conseq.det>minHours)    
  a=DAT
  
  if(plotWhat=="speed")
  {
    a$KM.hour=ifelse(a$hours.conseq.det>0,a$Dist.moved.conseq.det/a$hours.conseq.det,NA)                 
    a=subset(a,KM.hour>0,select=c(TagCode,KM.hour,Dist.moved.conseq.det))    
  }
  
  if(plotWhat=="distance")
  {
    a=subset(a,Dist.moved.conseq.det>0,select=c(TagCode,Dist.moved.conseq.det))  
  }
  
  
  names(a)[match("Dist.moved.conseq.det",names(a))]=c("Distance.km")
  
  first=subset(a,Distance.km<=range1[1] )
  N.first=nrow(first)
  mean.first=round(mean(first[,2],na.rm=T),r)
  median.first=round(median(first[,2],na.rm=T),r)
  SD.first=round(sd(first[,2],na.rm=T),r)
  
  second=subset(a,Distance.km<=range1[2] & Distance.km>range1[1])
  N.second=nrow(second)
  mean.second=round(mean(second[,2],na.rm=T),r)
  median.second=round(median(second[,2],na.rm=T),r)
  SD.second=round(sd(second[,2],na.rm=T),r)
  
  third=subset(a,Distance.km>range1[2])
  N.third=nrow(third)
  mean.third=round(mean(third[,2],na.rm=T),r)
  median.third=round(median(third[,2],na.rm=T),r)
  SD.third=round(sd(third[,2],na.rm=T),r)
  
  leg=function(N,MEAN,MED,SD)
  {
    legend("topright",c(paste("N=",N),paste("mean=",MEAN),paste("median=",MED),paste("SD=",SD)),bty='n',cex=1.25)
  }
  
  #plot
  hist(first[,2],main="",breaks=100,col=Singl.col,xlab="",ylab="",cex.axis=AX,cex.lab=AX,xlim=c(0,MAX1))
  leg(N.first,mean.first,median.first,SD.first)
  box()
  if(s==1)mtext(paste("Displacement <=",range1[1],"km"),3)
  
  hist(second[,2],main="",breaks=100,col=Singl.col,xlab="",ylab="",cex.axis=AX,cex.lab=AX,xlim=c(0,MAX2))
  leg(N.second,mean.second,median.second,SD.second)
  box()
  if(s==1)mtext(paste(range1[1],"km < Displacement <=",range1[2],"km"),3)
  
  hist(third[,2],main="",breaks=100,col=Singl.col,xlab="",ylab="",cex.axis=AX,cex.lab=AX,xlim=c(0,MAX3))
  leg(N.third,mean.third,median.third,SD.third)
  box()
  if(s==1)mtext(paste("Displacement >",range1[2],"km"),3)
  mtext(SPEC.nms[s],4,las=3,line=0.5)
  
}

Distance.ranges=list(c(10,50),c(10,50),c(10,50),c(10,50))
AX=1

#Speed distributions for different observation windows
tiff(file="Outputs_movement/Speed.dist.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfrow=c(4,3),mai=c(.25,.4,.2,.2),oma=c(1.2,1,.2,.2),las=1,mgp=c(1,.5,0))
for (s in 1:length(MIGRATE.list))
{
  Speed.dist.hist(subset(Detections,Species==SPECIES[s]),
                  Distance.ranges[[s]],"speed",4,4,4,3)
}
mtext("km/hour",1,outer=T,cex=1.5)
mtext("Frequency",2,outer=T,las=3,cex=1.5,line=-0.65)
dev.off()


#Displacement distributions for different observation windows
tiff(file="Outputs_movement/Displacement.dist.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfrow=c(4,3),mai=c(.2,.3,.2,.2),oma=c(1.2,1,.2,.2),las=1,mgp=c(1,.5,0))
for (s in 1:length(MIGRATE.list))
{
  Speed.dist.hist(subset(Detections,Species==SPECIES[s]),
                  Distance.ranges[[s]],"distance",10,70,2200,0)
}
mtext("Displacement (km)",1,outer=T,cex=1.25)
mtext("Frequency",2,outer=T,las=3,cex=1.25,line=-0.5)
dev.off()

  #by species speed distribution    
#note:  Analysis done for consequtive detections in different receivers more than 1 hour apart
#       This analysis is minimum possible speed (it underestimates speed as it assumes straight-line movement);
#       if detections are not really consecutive, eg shark is detected, wanders off outside range and the detected again
Speed.dist=function(DAT,range1)
{
  DAT=DAT[order(DAT$TagCode,DAT$DateTime.local),]
  DAT=subset(DAT,Rec.hours.conseq.det=="diff" & hours.conseq.det>minHours)    
  a=DAT
  a$KM.hour=ifelse(a$hours.conseq.det>0,a$Dist.moved.conseq.det/a$hours.conseq.det,NA)                 
  a=subset(a,KM.hour>0,select=c(TagCode,KM.hour,Dist.moved.conseq.det,hours.conseq.det))     
  names(a)[3:4]=c("Distance.km","Hours")
  speed=a
  MAX=round(max(speed$KM.hour))
  first=subset(speed,Distance.km<=range1[1] )
  max.first=round(max(first$KM.hour,na.rm=T),2)
  second=subset(speed,Distance.km<=range1[2] & Distance.km>range1[1])
  max.snd=round(max(second$KM.hour,na.rm=T),2)
  third=subset(speed,Distance.km>range1[2])
  max.thrd=round(max(third$KM.hour,na.rm=T),2)
  par(mfcol=c(3,1))
  hist(first$KM.hour,main=paste("Displacement <=",range1[1],"km"),breaks=100,col='grey',xlab="km/hour",
       cex.axis=1.5,cex.lab=1.5,xlim=c(0,MAX))
  legend("topright",paste("max speed=",max.first,"km/hour"),bty='n')
  box()
  if(nrow(second)>0)
{
    hist(second$KM.hour,main=paste(range1[1],"km > Displacement <=",range1[2],"km"),breaks=100,col='grey',
         xlab="km/hour",cex.axis=1.5,cex.lab=1.5,xlim=c(0,MAX))
    legend("topright",paste("max speed=",max.snd,"km/hour"),bty='n')
}
  box()
  hist(third$KM.hour,main=paste("Displacement >",range1[2],"km"),breaks=100,col='grey',xlab="km/hour",
       cex.axis=1.5,cex.lab=1.5,xlim=c(0,MAX))
  legend("topright",paste("max speed=",max.thrd,"km/hour"),bty='n')
  box()
}

Distance.ranges=list(c(10,100),c(10,100),c(10,50),c(10,50))

for (i in 1:length(MIGRATE.list))
{
  tiff(file=paste("Outputs_movement/Speed.dist.",SPECIES[i],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  Speed.dist(subset(Detections,Species==SPECIES[i]),Distance.ranges[[i]])
  dev.off()
}


#Plot speed, distance travelled and time between consecutive detections 
Store.speed.dist.day=store.MAXDIST
fn.speed.dist.days=function(dat)
{
  dat=subset(dat,hours.conseq.det>minHours)
  dat=subset(dat,Dist.moved.conseq.det>0 & days.conseq.det>0)
  dat$SPEED=dat$Dist.moved.conseq.det/dat$days.conseq.det  
  return(dat)
}
for(i in 1:N.sp) Store.speed.dist.day[[i]]=fn.speed.dist.days(Detections.species[[i]])

fn.plot=function(what,BKS,COLSS,MAIN)
{
  b=hist(what,breaks=BKS,col=COLSS,ylab="",xlab="",main=MAIN)
  n=max(what)*0.95
  y=max(b$counts)/2
  
  text(n,y,round(max(what)),col=COLSS,cex=1.25)
  arrows(n,y*0.8,n,0,col=COLSS)
  box()
}
SPECS=c("Dusky shark","Sandbar shark","Gummy shark","Whiskery shark")

#Distance
tiff(file="Outputs_movement/Distance.dist.all.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.75, 0),las=1)
for(i in 1:length(SPECS))fn.plot(Store.speed.dist.day[[i]]$Dist.moved.conseq.det,50,COLss[i],SPECS[i])
mtext("Minimum dist. moved between consec. det. (km)",1,outer=T,cex=1.5)
mtext("Frequency",2,outer=T,las=3,line=0.7,cex=1.5)
dev.off()

#Days between consecutive detections
tiff(file="Outputs_movement/Days.all.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.75, 0),las=1)
for(i in 1:length(SPECS))fn.plot(Store.speed.dist.day[[i]]$days.conseq.det,50,COLss[i],SPECS[i])
mtext("Days between consec. det.",1,outer=T,cex=1.5)
mtext("Frequency",2,outer=T,las=3,line=0.7,cex=1.5)
dev.off()

#speed (km/day)
tiff(file="Outputs_movement/Speed.all.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.75, 0),las=1)
for(i in 1:length(SPECS))fn.plot(Store.speed.dist.day[[i]]$SPEED,50,COLss[i],SPECS[i])
mtext("Minimum speed (km/day)",1,outer=T,cex=1.5)
mtext("Frequency",2,outer=T,las=3,line=0.7,cex=1.5)
dev.off()



#26. -- Cumulative Distance ---
all.tags=unique(Detections$TagCode)
Store.Cum.Dis=vector('list',length(all.tags))
for(i in 1:length(all.tags))
{
  a=subset(Detections,TagCode==all.tags[i])
  a$Cum.Dist=sum(a$Dist.moved.conseq.det)
  Store.Cum.Dis[[i]]=a[1,match(c("TagCode","Species","Sex","FL","Cum.Dist"),names(a))]
}
Store.Cum.Dis=do.call(rbind,Store.Cum.Dis)


#27. -- Patterns in residency, area use and speed ---


#Random Forest function

  #Get roaming data
fn.roaming=function(SPec,SPEC.nm,Tot.num.lines)
{
  Dat=subset(Detections,Species==SPec & !Recapture.hit=="YES")
  Dat$Line=with(Dat,
       ifelse(Latitude>=(-22.2),"Ning.L1",
       ifelse(Latitude<(-22.2) & Latitude>(-22.8),"Ning.L2",
       ifelse(Latitude<(-22.8) & Latitude>(-23.5),"Ning.L3",
       ifelse(Latitude<(-31.9) & Latitude>(-32.111) & Longitude>115.1255 & Longitude<115.5438,"Perth.outer",
       ifelse(Latitude<(-31.75) & Latitude>(-32.435) & Longitude>=115.5438 & Longitude<115.8,"Perth.inner",
       ifelse(Latitude<(-33.17) & Latitude>(-34.616) & Longitude>=114.4 & Longitude<115.424,"South.L1",
       ifelse(Latitude<(-32) & Longitude>116.3 & Longitude <=116.78,"South.L2",
       ifelse(Latitude<(-32) & Longitude>116.8 & Longitude <=120,"South.L3",NA
       )))))))))
  Dat$dummy=1
  
  
  #Data for random forest
  Dat$LineTagCode=paste(Dat$TagCode,Dat$Year,Dat$Line)
  Dat=Dat[!duplicated(Dat$LineTagCode),]
  Roaming=aggregate(dummy~TagCode+Sex+FL+Year+Zone.rel+Species,Dat,sum)
  Roaming$Roaming=Roaming$dummy/Tot.num.lines
  Roaming.RF=subset(Roaming,Sex%in%c("F","M"))
  Roaming.RF$Sex=as.character(Roaming.RF$Sex)
  Roaming.RF=Roaming.RF[,-match("dummy",names(Roaming.RF))]
  Roaming.RF$Total=Tot.num.lines
  
  #Data for simple figure
  Dat$LineTagCode=paste(Dat$TagCode,Dat$Line)
  Dat=Dat[!duplicated(Dat$LineTagCode),]
  Roaming=aggregate(dummy~TagCode+Sex,Dat,sum)
  if(SPec%in%c("Gummy","Whiskery"))Tot.num.lines=5
  if(SPec%in%c("Dusky","Thickskin"))Tot.num.lines=8
  Roaming$Roaming=Roaming$dummy/Tot.num.lines
  Roaming=subset(Roaming,Sex%in%c("F","M"))  
  Roam=aggregate(Roaming~Sex,Roaming, mean)
  Roam.SD=aggregate(Roaming~Sex,Roaming, sd)
  return(list(Roam=Roam,Roam.SD=Roam.SD,Data.RF=Roaming.RF))
}
Store.Roaming=vector('list',length(SPEC.nms))
names(Store.Roaming)=SPEC.nms
NNSp=1:length(All.spec)
TOT.LINES=c(7,7,7,7)
names(TOT.LINES)=All.spec
for(p in NNSp) Store.Roaming[[p]]=fn.roaming(All.spec[p],SPEC.nms[p],TOT.LINES[p])
Data.RF.Roaming=rbind(Store.Roaming$Dusky$Data.RF,
                      Store.Roaming$Sandbar$Data.RF,
                      Store.Roaming$Gummy$Data.RF,
                      Store.Roaming$Whiskery$Data.RF)

  #Put residency data in right format

a=subset(Time.mon.zone,!(Species%in%c("Gummy","Whiskery")&Latitude>(-26)))
a$Year=as.character(substr(a$Date.local,1,4))
b=aggregate(Days~zone+TagCode+Species+Sex+Year,a,sum,na.rm=T)
Time.monit.yr=aggregate(Days~TagCode+Year,b,sum)
names(Time.monit.yr)[3]="days.mon"
b=merge(b,Time.monit.yr,by=c("TagCode","Year"),all.x=T)
b$prop.time.in.zn=b$Days/b$days.mon
b=b[,-match("Days",names(b))]  
Add.FL=subset(Detections,TagCode%in%unique(b$TagCode),select=c(TagCode,FL,Zone.rel))
Add.FL=Add.FL[!duplicated(Add.FL$TagCode),]
b=merge(b,Add.FL,by=c("TagCode"),all.x=T)
names(b)[match("prop.time.in.zn",names(b))]="Residency"
Data.RF.Residency=b
rm(a,b)
Data.RF.Residency=subset(Data.RF.Residency,!Zone.rel=="SA")

Anova.and.Dev.exp=function(MODEL)
{
  #Anovas
  Anova.tab=anova(MODEL, test = "Chisq")
  
  #Deviance explained
  
  #By each term
  n=2:length(Anova.tab$Deviance)
  Term.dev.exp=100*(Anova.tab$Deviance[n]/MODEL$null.deviance)
  names(Term.dev.exp)=rownames(Anova.tab)[n]
  
  #By full model
  Dev.exp=sum(Term.dev.exp)
  
  return(list(Anova=Anova.tab,Dev.exp=Dev.exp,Term.Dev.exp=Term.dev.exp))
}
fn.RF.importance=function(fit,Predictors,DATA)
{
  #Extract predictor importance
  Importance=importance(fit, type=1)
  Importance=Importance[rev(order(Importance[,1])),]
  
  #Extract relative importance of each predictor level
  Lvl.Rel.Imp=vector('list',length(Predictors))
  names(Lvl.Rel.Imp)=Predictors
  for(p in 1:length(Predictors)) Lvl.Rel.Imp[[p]]=partialPlot(fit,DATA,Predictors[p],plot=F)
  
}
Fn.models=function(DATA,Response,Predictors,Ntree,WEIGHT,add.inter)
{
  set.seed(4543)
  
  DATA=subset(DATA,!Sex=="U")
  #Set factors
  for(f in 1:length(Factors)) DATA[,match(Factors[f],names(DATA))]=as.factor(DATA[,match(Factors[f],names(DATA))])
  
  #remove NAs in response variable
  id=match(Response,colnames(DATA))
  DATA=DATA[!is.na(DATA[,id]),]
  
  # write formula 
  Preds=paste(Predictors, collapse="+")
  if(add.inter=="YES")
  {
    Ids=match(c("Species","Zone.rel"),Predictors)
    Inter=paste(Predictors[Ids], collapse="*")
    Preds=paste(c(Inter,Predictors[-Ids]), collapse="+")
    
  }
  Formula <- as.formula(paste(Response, "~", Preds,collapse=NULL))
  
  
  #1. fit binomial GLM with
  DATA$weights=DATA[,match(WEIGHT,names(DATA))]
  GLM=glm(Formula, weights = weights,family = binomial, data = DATA)
  #GLM=glmer(Residency ~ Year + Species + Sex + Zone.rel +zone + FL+ (1 | TagCode), weights = weights,family = binomial, data = DATA)
  
  # #2.fit random forest
  # fit= randomForest(Formula,data=DATA,importance=T,na.action=na.roughfix, ntree=Ntree,proximity=T)###can add more trees
  # 
  # # SHOWS ERROR ASSOCIATED WITH VARIABLES (does error stabilize?)
  # plot(fit, log="y",main=Response)
  # 
  # #Display predictor importance
  # varImpPlot(fit,main=Response)
  # 
  # #Extract variance explained (as %)
  # Var.exp=fit$rsq 
  # Var.exp=100*Var.exp[length(Var.exp)]
  
  # return(list(GLM=GLM,RF=fit,DATA=DATA,Perc.Var.exp=Var.exp))
  return(list(GLM=GLM,DATA=DATA))
}

#note: "TagCode" has too many levels for random forests


  #1. Roaming
PREDS=c("Species","Zone.rel","FL","Sex","Year")
Factors=c("Species","Year","Sex","Zone.rel")
Roaming.model=Fn.models(DATA=Data.RF.Roaming,Response="Roaming",
                        Predictors=PREDS,Ntree=500,WEIGHT="Total",add.inter="NO")

    #1.1 deviance explained
#Roaming.RF.var.exp=Roaming.model$Perc.Var.exp
Roaming.GLM.var.exp=Anova.and.Dev.exp(Roaming.model$GLM)


  #2. Residency
PREDS.res=c(PREDS)
Factors=c(Factors)
Res.zones=c("Closed.ningaloo","West","Closed.metro","Zone1","Zone2")
names(Res.zones)=c("Ningaloo","WCDGDLF","Metro","Zone1","Zone2")
Residency.model=vector('list',length(Res.zones))
names(Residency.model)=Res.zones
for(z in 1:length(Res.zones)) Residency.model[[z]]=Fn.models(DATA=subset(Data.RF.Residency,zone==Res.zones[z]),
                       Response="Residency",Predictors=PREDS.res,Ntree=500,WEIGHT="days.mon",add.inter="NO")


    #2.1 deviance explained
#Residency.RF.var.exp=Residency.model$Perc.Var.exp
Residency.GLM.var.exp=Anova.table.Res=Residency.model
for(z in 1:length(Res.zones)) Residency.GLM.var.exp[[z]]=Anova.and.Dev.exp(Residency.model[[z]]$GLM)


ANOV.table=function(D)
{
  TAB=as.data.frame(D$Anova)
  TAB=TAB[-1,]
  TAB$Term=rownames(TAB)
  Per.dev.exp=data.frame(Term=names(D$Term.Dev.exp),Per.dev.exp=D$Term.Dev.exp)
  TAB=merge(TAB,Per.dev.exp,by="Term")
  TAB=TAB[,match(c("Term","Per.dev.exp","Pr(>Chi)"),names(TAB))]
  TAB=TAB[order(-TAB$Per.dev.exp),]
  TAB[,2:3]=round(TAB[,2:3],2)
  TAB[TAB[,2]<0.01,2]="<0.001"
  TAB[TAB[,3]<0.01,3]="<0.001"
  Total=TAB[1,]
  Total$Term="Total"
  Total$Per.dev.exp=round(D$Dev.exp,2)
  Total$"Pr(>Chi)"=""
  return(rbind(TAB,Total))
}
Anova.table.Roam=ANOV.table(D=Roaming.GLM.var.exp)
colnames(Anova.table.Roam)[2:3]=paste("Roam_",colnames(Anova.table.Roam)[2:3],sep="")
Anova.table.Roam=Anova.table.Roam[match(c("Species","Zone.rel","Year","FL","Sex","Total"),Anova.table.Roam$Term),]

for(z in 1:length(Res.zones)) 
{
  Anova.table.Res[[z]]=ANOV.table(D=Residency.GLM.var.exp[[z]])
  colnames(Anova.table.Res[[z]])[2:3]=paste("Res_",colnames(Anova.table.Res[[z]])[2:3],sep="")
  Anova.table.Res[[z]]=Anova.table.Res[[z]][match(c("Species","Zone.rel","Year","FL","Sex","Total"),Anova.table.Res[[z]]$Term),]
}


#Plot RF term level importance
#plot(Roaming.model$Lvl.Rel.Imp[[p]]$x,Roaming.model$Lvl.Rel.Imp[[p]]$y)



#28. -- MEPS figures ---
#All figures and tables for MEPS paper
MEPS.paper=function(do.MEPS.outptus)
{
  PaZ=getwd()
  if(do.MEPS.outptus=="YES")
  {
    setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement/MEPS_paper"))
    
    #Table 1
    Tab1=Table1[,-ncol(Table1)]
    

    #add distance and speed table
    Summary.mig=function(Dat)
    {
      #Consecutive detections is unique so keep all
      Dat$Km.per.day=Dat$Dist.moved.conseq.det/Dat$days.conseq.det
      ID=c("TagCode","Dist.moved.conseq.det","days.conseq.det","Km.per.day")
      MaxDist=Dat[order(Dat$TagCode,Dat$Dist.moved.conseq.det),ID]
      
      Mean.speed=round(mean(MaxDist$Km.per.day,na.rm=T),2)
      sd.speed=round(sd(MaxDist$Km.per.day,na.rm=T),2)
      Max.speed=round(max(MaxDist$Km.per.day,na.rm=T),2)
      
      Mean.dist=round(mean(MaxDist$Dist.moved.conseq.det,na.rm=T))
      sd.dist=round(sd(MaxDist$Dist.moved.conseq.det,na.rm=T))
      Max.dist=round(max(MaxDist$Dist.moved.conseq.det,na.rm=T),2)
      
      MaxDist[,2:3]=round(MaxDist[,2:3])
      MaxDist[,4]=round(MaxDist[,4],2)
      
      N.detections=nrow(MaxDist)
      N.sharks=length(unique(MaxDist$TagCode))
      
      return(list(MaxDist=MaxDist,Stats=rbind(N.detections,N.sharks,Mean.dist,sd.dist,
                                              Max.dist,Mean.speed,sd.speed,Max.speed)))
    }
    
    store.MAXDIST=MIGRATE.list
    for (i in 1:length(MIGRATE.list))
    {
      MAXDIST=Summary.mig(MIGRATE.list[[i]])$MaxDist
      Thress=Migration.threshold[i]
      #max distance
      store.MAXDIST[[i]]=MAXDIST
    }
    
    #export table of stats
    n=nrow(Summary.mig(MIGRATE.list[[i]])$Stats)
    Sum=matrix(nrow=n,ncol=4)
    for(i in 1:ncol(Sum)) Sum[,i]=Summary.mig(MIGRATE.list[[i]])$Stats
    rownames(Sum)=rownames(Summary.mig(MIGRATE.list[[i]])$Stats)
    colnames(Sum)=names(MIGRATE.list)
    Sum1=cbind(paste(rownames(Sum),"_long_distance",sep=""),Sum)
    Tab1=rbind(Tab1,Sum1)
    
    Scenarios.tbl(WD=getwd(),Tbl=Tab1,Doc.nm="Table1",caption=NA,paragph=NA,
                  HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                  Zebra='NO',Zebra.col='grey60',Grid.col='black',
                  Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
                  HEDR=c('Variable','Dusky','Thickskin','Gummy','Whiskery'),
                  HEDR.cols=c(1,1,1,1,1),HEDR2=NA,HEDR3=NA)
    
    #Figure 1. Map
    #define coordinates of polygons
    N.WA.long=c(North.WA.long[2], North.WA.long[2], North.WA.long[1], North.WA.long[1])
    N.WA.lat=c(North.WA.lat[2], North.WA.lat[1], North.WA.lat[1], North.WA.lat[2])
    S.WA.long=c(South.WA.long[2], South.WA.long[2], South.WA.long[1], South.WA.long[1])
    S.WA.lat=c(South.WA.lat[2], South.WA.lat[1], South.WA.lat[1], South.WA.lat[2])
    
    #define Perth and Rotnest
    Perth=c(115.866,-31.95)
    Rotnest=c(115.50,-32.02)
    
    #bathymetry
    Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
    xbat=sort(unique(Bathymetry$V1))
    ybat=sort(unique(Bathymetry$V2)) 
    if(!exists("reshaped"))
    {
      reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
      rm(Bathymetry)
    }
    
    #legends
    Letter.leg=c("A",NA,"B")
    Letter.leg.coor=cbind(c(113.61,NA,114.45),c(-21.61,NA,-31.85))
    m <- rbind(c(0.0385, 0.3, 0, 1),
               c(0.35, 1, 0.5, 1),
               c(0.3, 1, 0, 0.5))
    
    tiff(file="Figure1.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    split.screen(m)
    for(i in 1:3)
    {
      screen(i)
      par(mar = c(0, 0, 0, 0),mgp=c(.1, 0.15, 0))
      if(i==2) plotlong[[i]][1]=109
      plotMap(worldLLhigh, xlim=plotlong[[i]],ylim=plotlat[[i]],plt = c(.001, 1, 0.075, 1),
              col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
      if(i==1)
      {
        text(114,-22.6,"Ningaloo",cex=1.4)
        text(114,-22.7,"Reef",cex=1.4)
        text(114,-22.8,"array",cex=1.4)        
       # text(113.94,-23.134,"Coral Bay",cex=1)
        text(114.04,-22.03,"North West Cape",srt=75,cex=1)
     
      }
        
      if(i==2)
      {
        #add zones
        plot(WA_Northern_Shark,ylim=c(-39,-11),xlim=c(108,130),add=T,col="grey50")
        plot(JA_Northern_Shark,add=T,col="grey70")
        plot(WA_Northern_Shark_2,add=T,col="grey90")
        plot(WCDGDLL,add=T,col="white")
        plot(SDGDLL_zone1,add=T,col="grey80")
        plot(SDGDLL_zone2,add=T,col="grey60")
        text(117.65,-25.5,"Shark Bay",cex=.9)
        
        segments(113,-40,124.75,-40,lwd=4)
        text(130,-40,"1000 km",cex=1.5)
        #library(SDMTools)#add scale bar
        #Scalebar(x=113,y=123.6,distance=1000)
        
        text(126,-13.4,"JANSF",srt=10,cex=1)
        text(119,-17,"WANCSF",srt=45,cex=0.85)
        text(112.225,-22,"Ningaloo",cex=1,srt=85)
        text(112.2,-29.95,"WCDGDLF",cex=1,srt=85)
        text(114,-35.6,"Zone 1",cex=1,srt=330)
        text(122.5,-36.5,"Zone 2",cex=1)
        
        
        polygon(x=N.WA.long,y=N.WA.lat,lwd=2)
        polygon(x=S.WA.long,y=S.WA.lat,lwd=2)
        text(133,-25,("Australia"),col="black", cex=2)
        mtext("Latitude (ºS)",side=2,line=0,las=3,cex=1.5)
        mtext("Longitude (ºE)",side=1,line=0,cex=1.5)
        text(115.98,-22.6,("A"),col="black", cex=1.75)
        text(119,-33,("B"),col="black", cex=1.75)
        
        text(120.25,-29.5,"Metro")
        arrows(117.75,-30,115.68,-31.2,lwd=1.5,length=0.1)
      }
      if(!i==2)
      {
        points(Receiverlong[[i]],Receiverlat[[i]],col='black',pch=19,cex=0.6)  #receiver location
        axis(side = 1, at =Long.seq[[i]], labels = Long.seq[[i]], tcl = .5,las=1,cex.axis=0.9)
        axis(side = 2, at = Lat.seq[[i]], labels = -Lat.seq[[i]],tcl = .5,las=2,cex.axis=0.9)
        box(lwd=2)
        contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=plotlat[[i]],xlim=plotlong[[i]], zlim=c(-1,-300),
                nlevels = 3,labcex=1.25,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
        text(Letter.leg.coor[i,1],Letter.leg.coor[i,2],Letter.leg[i],cex=2.5)
      }
      
      if(i==3)
      { 
        #polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="red")
        points(115.5,-32.0,pch=19,col='dark grey',cex=1.125)
        text(115,-32.75,"Rottnest",cex=1,srt=60)
        text(115.2,-32.9,"Island",cex=1,srt=60)
        arrows(115.4078,-32.3714,115.4932,-32.129,lwd=1.5,length=0.1)
        #polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
        text(117.75,Perth[2],"Perth array",cex=1.4)
        text(117.1,-34.2,"Southern",cex=1.4)
        text(117.15,-34.5,"Lines array",cex=1.4)
        
        text(115.85,-33.8,"Cape Leeuwin",cex=1,srt=30)
        text(116.17,-31.95,"Perth",cex=1)
        text(117.88,-34.85,"Albany",cex=1)
      }
      
      
    }
    close.screen(all = TRUE)
    dev.off()
    
    #Figure 2. prop monitored time detected
    fn.n.days.monitor=function(SPec,SPEC.nm)
    {
      Dat=subset(Detections,Species==SPec & !(TagCode.original%in%only.recaptured))
      Tgs=unique(Dat$TagCode)
      n.Tgs=length(Tgs)
      b=subset(Total.time.monitored,TagCode%in%Tgs,select=c(TagCode,days.mon))
      a=matrix(nrow=length(Tgs),ncol=2)
      for(e in 1:n.Tgs)
      {
        qq=subset(Dat,TagCode==Tgs[e])
        a[e,]=cbind(Tgs[e],length(unique(qq$Date.local)))
      }
      a=as.data.frame(a) 
      names(a)=c("TagCode","Days.detected")
      
      a=merge(a,b,by="TagCode")
      a$Days.detected=as.numeric(as.character(a$Days.detected))
      a$Prop=a$Days.detected/a$days.mon
      a=a[order(a$days.mon),]
      plot(1:n.Tgs,a$Prop,type='h',ylab="",xlab="",xaxt="n",lwd=3,col=2,ylim=c(0,1.15))
      axis(1,1:n.Tgs,a$TagCode,las=3,cex.axis=.65)
      mtext(SPEC.nm,3,cex=0.9)
      text(1:n.Tgs,a$Prop,a$days.mon,pos = 3,cex=.7,srt=45)
    }
    SPEC.nMs=paste(SPEC.nms,"shark")
    
    tiff(file="Figure2.tiff",width = 2400, height = 1600,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(4,1),mar=c(2,1.75,2,0.1),oma=c(2,2,0.001,0.1),mgp=c(2, 0.75, 0))
    for(p in 1:length(All.spec)) fn.n.days.monitor(All.spec[p],SPEC.nMs[p])
    mtext("Porportion of days detected",2,outer=T,cex=1.15,line=0.5)
    mtext("Tag ID",1,outer=T,cex=1.15,line=0.825)
    dev.off()
    
    #Figure 3. time per zone 
    in.colors="NO"
    
    COL.prop=c("lightseagreen","seagreen4","lightgreen","olivedrab4","olivedrab3","mediumseagreen")
    if(in.colors=="NO")
    {
      colfunc <- colorRampPalette(c("grey20","grey95"))
      #COL.prop=rev(colfunc(length(COL.prop)))
      COL.prop=c("black","grey90","grey30","grey80","grey50","grey70")
      
    }
    names(COL.prop)=Zns
    
      #all individuals
    fn.plot.prop.MEPS=function(what)
    {
      a=subset(Species.time.zone,Species==what)
      id=match(Zns,names(a))
      a[,id]=a[,id]/rowSums(a[,id])
      crap=as.data.frame(as.matrix(COL.prop))
      names(crap)="Col"
      crap$Zone.rel=rownames(crap)
      a=merge(a,crap,by="Zone.rel")
      a$Sort=with(a,ifelse(Zone.rel=="North",1,
                           ifelse(Zone.rel=="Closed.ningaloo",2,
                                  ifelse(Zone.rel=="West",3,
                                         ifelse(Zone.rel=="Closed.metro",4,
                                                ifelse(Zone.rel=="Zone1",5,6))))))
      a=a[order(a$Sort),]
      
      r=barplot(t(a[,id]), col = COL.prop,horiz=T,beside=F,yaxt='n',cex.axis=1.25)
      box()
      points(rep(-0.045,length(r)),r,pch=15,cex=1.5,col=as.character(a$Col))
      
    }
    
    fn.par=function()par(mar = c(1, .9, 1.5, 0.2),oma=c(1.8,1.2,0.1,0.5),xpd=TRUE,mgp=c(1,.5,0))
    fn.LEG=function(wat) legend("top",wat,inset=c(0,-.15),bty='n',cex=1.25)
    
    m <- rbind(c(0, .475, .55, 1),
               c(.5, 1, .55, 1),
               c(0, .475, .1, .53),
               c(.5, 1, .1, .53),
               c(0, 1, 0, .075))
    tiff(file="Figure 3_all_indiv.tiff",width = 2400, height = 2200,units = "px", res = 300,compression = "lzw")
    split.screen(m)
    
    screen(1)
    fn.par()
    fn.plot.prop.MEPS("Dusky")
    
    fn.LEG('Dusky shark')
    screen(3)
    fn.par()
    fn.plot.prop.MEPS("Gummy")
    fn.LEG("Gummy shark")
    screen(2)
    fn.par()
    fn.plot.prop.MEPS("Thickskin")
    fn.LEG("Sandbar shark")
    screen(4)
    fn.par()
    fn.plot.prop.MEPS("Whiskery")
    fn.LEG("Whiskery shark")
    
    screen(5)
    fn.par()
    par(mar = c(0, 0, 0, 0),oma=c(0,0,0,0),xpd=TRUE,mgp=c(1,.5,0))
    legend("center",Zns.leg,bty='n',pt.cex=2,pch=15,col=COL.prop,text.width = max(sapply(Zns.leg, strwidth)),
           horiz=T,cex=1.1) 
    mtext("             Release zone",2,line=-1.25,cex=1.5,outer=T)
    mtext("Proportion of time",1,line=-3.75,cex=1.5)
    close.screen(all = TRUE)
    dev.off()
    
    
    #Individuals combined 
    colfunc <- colorRampPalette(c("grey15","white"))
    couleurs=rev(colfunc(N.int))
    
    tiff(file="Figure3.tiff",width = 2000, height = 2000,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(2,2),mar = c(1, 1, 1, 1),oma=c(3,4,0,0),xpd=T,mgp=c(1,.7,0))
    write.csv(Image.mig.fn("Dusky","Dusky shark"),"Prop.time.image.Dusky.csv")
    write.csv(Image.mig.fn("Thickskin","Sandbar shark"),"Prop.time.image.Sandbar.csv")    
    write.csv(Image.mig.fn("Gummy","Gummy shark"),"Prop.time.image.Gummy.csv")
    write.csv(Image.mig.fn("Whiskery","Whiskery shark"),"Prop.time.image.Whiskery.csv")
    SQ=seq(BREAKS[1],BREAKS[length(BREAKS)],.1)
    color.legend(6.5,0.5,7,6.5,SQ,rect.col=couleurs,gradient="y",col=1,cex=1)
    mtext("Zone released",2,outer=T,line=2.5,cex=1.5)
    mtext("Zone detected",1,outer=T,line=1.5,cex=1.5)
    dev.off()
    
    
    #Figure 4. Movement among zones
    fn.plt=function(what,xLAB)
    {
      a=barplot(what,horiz=T,xlab="",main=xLAB,ylab="",col=1,cex.main=1.5,cex.axis=1.25)
      axis(2,a,c("Dusky","Sandbar","Gummy","Whiskery"),las=1,cex.axis=1.35)
      box()
    }
    tiff(file="Figure4.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(4,1),mai=c(.5,.75,.2,.1),oma=c(.2,1,.2,.1),mgp=c(1,.7,0))
    
    #Movement among adjance zones
    for(s in 1:length(sps))Time.adj.zn[[s]]=fn.time.zone(subset(Tagging.pop.dyn,Species==sps[s]),ADJ="YES")
    fn.plt(cbind(Time.adj.zn[[1]][2],Time.adj.zn[[2]][2],Time.adj.zn[[3]][2],Time.adj.zn[[4]][2]),
           "Number of individuals that moved to an adjacent zone in one month")
    fn.plt(cbind(Time.adj.zn[[1]][3],Time.adj.zn[[2]][3],Time.adj.zn[[3]][3],Time.adj.zn[[4]][3]),
           "Number of individuals that moved to an adjacent zone in one year")
    
    #Movement among non-adjance zones    
    for(s in 1:length(sps))Time.adj.zn[[s]]=fn.time.zone(subset(Tagging.pop.dyn,Species==sps[s]),ADJ="NO")
    fn.plt(cbind(Time.adj.zn[[1]][2],Time.adj.zn[[2]][2],Time.adj.zn[[3]][2],Time.adj.zn[[4]][2]),
           "Number of individuals that moved to a non-adjacent zone in one month")
    fn.plt(cbind(Time.adj.zn[[1]][3],Time.adj.zn[[2]][3],Time.adj.zn[[3]][3],Time.adj.zn[[4]][3]),
           "Number of individuals that moved to a non-adjacent zone in one year")
    
    dev.off()


    #Figure Suppl. 3. Proportion staying in zone
    Prop.stay.fn=function(what)
    {
      a=subset(Species.time.zone,Species==what)  
      a$North.stay=with(a,ifelse(North==1 & Zone.rel=="North",1,0))
      a$Closed.ningaloo.stay=with(a,ifelse(Closed.ningaloo==1 & Zone.rel=="Closed.ningaloo",1,0))
      a$West.stay=with(a,ifelse(West==1 & Zone.rel=="West",1,0))
      a$Closed.metro.stay=with(a,ifelse(Closed.metro==1 & Zone.rel=="Closed.metro",1,0))
      a$Zone1.stay=with(a,ifelse(Zone1==1 & Zone.rel=="Zone1",1,0))
      a$Zone2.stay=with(a,ifelse(Zone2==1 & Zone.rel=="Zone2",1,0))
      
      add.Sex=subset(Detections,Species==what & TagCode%in%unique(a$TagCode))
      add.Sex=add.Sex[!duplicated(add.Sex$TagCode),match(c("TagCode","Sex"),names(add.Sex))]
      a=merge(a,add.Sex,by="TagCode",all.x=T)
      a=subset(a,Sex%in%c("F","M"))
      a$Dummy=1  
      Prop.stay=aggregate(cbind(North.stay,Closed.ningaloo.stay,
                                West.stay, Closed.metro.stay, Zone1.stay, Zone2.stay)~Sex,a,sum)
      Tot=aggregate(Dummy~Sex,a,sum)  
      names(Tot)[2]="Total"
      Prop.stay$tot.stay=rowSums(Prop.stay[,-match("Sex",colnames(Prop.stay))])
      Prop.stay=merge(Prop.stay,Tot,by="Sex",all.x=T)
      Prop.stay$prop=Prop.stay$tot.stay/Prop.stay$Total
      Prop.stay$Species=what
      return(Prop.stay[,match(c("Species","Sex","Total","prop"),names(Prop.stay))])
    }
    Prop.stay.zone=vector('list',length(SPEC.nms))
    for(i in 1:length(SPEC.nms))Prop.stay.zone[[i]]=Prop.stay.fn(All.spec[i])
    Prop.stay.zone=do.call(rbind,Prop.stay.zone)
    Prop.stay.zone.r=reshape(Prop.stay.zone,v.names = c("prop","Total"), idvar = c("Sex"),
                             timevar = c("Species"), direction = "wide")
    Prop.stay.zone.r.prop=Prop.stay.zone.r[,c(1:2,4,6,8)]
    Prop.stay.zone.r.tot=Prop.stay.zone.r[,c(1,3,5,7,9)]
    
    tiff(file="Prop staying zone.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(1,1),mai=c(.5,.65,.2,.1),oma=c(.2,1,.2,.1),mgp=c(1,.7,0),las=1)
    
    mids=barplot(as.matrix(Prop.stay.zone.r.prop[,-1]),beside=T,col=c("pink","blue"),
                 xaxt="n",ylim=c(0,1.075),cex.axis=1.25)
    axis(1,colMeans(mids),paste(SPEC.nms,"shark"),cex.axis=1.25)
    box()
    text(c(mids),c(unlist(Prop.stay.zone.r.prop[,-1])),c(unlist(Prop.stay.zone.r.tot[,-1])),pos=3,cex=1.25)
    mtext("Proportion staying within release zone",2,outer=F,cex=1.5,line=2.5,las=3)
    legend("topright",c("Female","Male"),fill=c("pink","blue"),bty="n",cex=1.5)
    dev.off()

    #Figure Residency and Roaming index
    Col.fem="grey20"
    Col.male="grey75"
    #Col.fem="pink"
    #Col.male="blue"
    fn.SD.bar=function(Y,Y.SD,CLs,X)
    {
      segments(X,Y,X,Y+Y.SD,col=CLs,lwd=2.5)
      segments(X,Y,X,Y-Y.SD,col=CLs,lwd=2.5)
    }
    
    fn.residency=function(SPec,SPEC.nm)
    {
      Dat=subset(Species.time.zone,Species==SPec)
      add.Sex=subset(Detections,Species==SPec & TagCode%in%unique(Dat$TagCode))
      add.Sex=add.Sex[!duplicated(add.Sex$TagCode),match(c("TagCode","Sex"),names(add.Sex))]
      Dat=merge(Dat,add.Sex,by="TagCode",all.x=T)
      Dat=subset(Dat,Sex%in%c("F","M"))
      
      Dat$Residency.to.zone=NA
      for (qq in 1:nrow(Dat))
      {
        M=match(Dat$Zone.rel[qq],colnames(Dat))
        Dat$Residency.to.zone[qq]=Dat[qq,M]
      }
      Dat$Sex=as.character(Dat$Sex)
      Dat$Zone.rel=as.character(Dat$Zone.rel)
      Dat$Zone.rel = factor(Dat$Zone.rel,
                            levels=c("North","Closed.ningaloo","West","Closed.metro",  
                                     "Zone1","Zone2"))
      
      Res=tapply(Dat$Residency.to.zone, list(Sex=Dat$Sex, Zone.rel=Dat$Zone.rel), mean) 
      Res.SD=tapply(Dat$Residency.to.zone, list(Sex=Dat$Sex, Zone.rel=Dat$Zone.rel), sd)
      
      
      plot(1:ncol(Res),Res[1,],pch=19,cex=2.5,col=Col.fem,cex.axis=1.5,
           ylim=c(0,1),xaxt="n",las=1,ylab="",xlab="")
      fn.SD.bar(Res[1,],Res.SD[1,],Col.fem,1:ncol(Res))
      xx=jitter(1:ncol(Res))
      points(xx,Res[2,],pch=19,cex=2.5,col=Col.male)
      fn.SD.bar(Res[2,],Res.SD[2,],Col.male,xx)
      axis(1,1:ncol(Res),F)
      mtext(SPEC.nm,3,line=0.25,cex=1.15)
      
      #boxplot
      #     Dat$Sex.Zone=paste(Dat$Zone.rel,Dat$Sex)
      #     
      #     Dat$Sex.Zone = factor(Dat$Sex.Zone,
      #           levels=c(  "North F","North M",
      #                      "Closed.ningaloo F", "Closed.ningaloo M",                        
      #                      "West F","West M", "Closed.metro F","Closed.metro M" ,  
      #                      "Zone1 F", "Zone1 M","Zone2 F","Zone2 M"))
      #     
      #     boxplot(Residency.to.zone~Sex.Zone,Dat,col=1:10)
      
      #violin plot
      #     library(vioplot)
      #     Uniks=unique(Dat$Sex.Zone)
      #     yalist=vector('list',length(Uniks))
      #     for(x in 1:length(Uniks))yalist[[x]]=subset(Dat,Sex.Zone==Uniks[x])
      #     CL=rainbow(length(yalist))
      #     plot(0,0,type="n",xlim=c(0,length(yalist)),ylim=c(0,1),xaxt='n',xlab="",ylab="Pc [%]")
      #     for (i in 1:length(yalist))
      #     { vioplot(na.omit(yalist[[i]]$Residency.to.zone), at = i, add = T, col = CL[i]) }
      #     axis(side=1,at=1:3,labels=3:1)
      
      
      
      
    }
    
    set.seed(999)
    tiff(file="Figure.Residency and roaming.tiff",width = 1800, height = 2400,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(5,1),mar=c(2,3,1,2),oma=c(2,2,1,0.1),mgp=c(2, 0.75, 0))
    
    #Residency index
    for(p in 1:length(All.spec)) fn.residency(All.spec[p],SPEC.nMs[p])
    axis(1,1:6,c("WANCSF","Ningaloo","WCDGDLF","Metro","Zone1","Zone2 "),cex.axis=1.5)
    #mtext("Release zone",1,line=2.5,cex=1.75)
    mtext("                          Residency index",2,line=0,cex=1.25,outer=T)
    
    #Roaming index
    plot(1:length(All.spec),cex.axis=1.5,ylim=c(0,1),col="transparent",xaxt="n",las=1,ylab="",xlab="")
    for(p in NNSp)
    {
      points(NNSp[p],Store.Roaming[[p]]$Roam$Roaming[1],pch=19,col=Col.fem,cex=2.5)
      fn.SD.bar(Store.Roaming[[p]]$Roam$Roaming[1],Store.Roaming[[p]]$Roam.SD$Roaming[1],Col.fem,NNSp[p])
      
      xx=jitter(NNSp[p])
      points(xx,Store.Roaming[[p]]$Roam$Roaming[2],pch=19,col=Col.male,cex=2.5)
      fn.SD.bar(Store.Roaming[[p]]$Roam$Roaming[2],Store.Roaming[[p]]$Roam.SD$Roaming[2],Col.male,xx)
      
    }
    axis(1,NNSp,c("Dusky","Sandbar","Gummy","Whiskery"),cex.axis=1.5)
    mtext("Roaming index",2,line=3,cex=1.25)
    legend("topright",c("female","male"),fill=c(Col.fem,Col.male),bty="n",cex=1.5)
    
    dev.off()
    
    
    #Anova table Residency and Roaming
    write.csv(Anova.table.Roam,"Anova.table.Roam.csv",row.names=F)
    for(z in 1:length(Res.zones)) write.csv(Anova.table.Res[[z]],paste("Anova.table.Res_",Res.zones[[z]],".csv",sep=""),row.names=F)
    # Scenarios.tbl(WD=getwd(),Tbl=Anova.table.Roam,Doc.nm="Anova.table",
    #               caption=NA,paragph=NA,HdR.col='black',HdR.bg='white',
    #               Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
    #               Zebra='NO',Zebra.col='grey60',Grid.col='black',
    #               Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
    #               HEDR=c('Term','Residency','Roaming'),
    #               HEDR.cols=c(1,2,2),
    #               HEDR2=c("","Per.dev.exp","P","Per.dev.exp","P"),HEDR3=NA)
    

    # fn.plt.glm2=function(MOD,DAT,What,Predictors,Axis.Labl,PLOT)
    # {
    #   Predictors=Predictors[-match(What,Predictors)]
    #   NEW1=factor(as.character(unique(DAT[,match(What[1],names(DAT))])),levels=levels(DAT[,match(What[1],names(DAT))]))
    #   NEW2=factor(as.character(unique(DAT[,match(What[2],names(DAT))])),levels=levels(DAT[,match(What[2],names(DAT))]))
    #   NEW=expand.grid(NEW1,NEW2)
    #   names(NEW)=What
    #   Add=DAT[,match(Predictors,names(DAT))]
    #   stOrE=rep(NA,ncol(Add))
    #   for(n in 1:ncol(Add))
    #   {
    #     aa=sort(table(Add[,n])) 
    #     stOrE[n]=names(aa[length(aa)])
    #   }
    #   
    #   Add2=as.data.frame(t(stOrE))
    #   colnames(Add2)=colnames(Add)
    #   for(n in 1:length(stOrE))
    #   {
    #     Cls=class(Add[,n])
    #     if(Cls=='numeric') Add2[,n]=as.numeric(as.character(Add2[,n]))
    #     if(Cls=='factor') Add2[,n]=factor(Add2[,n],levels=levels(DAT[,match(Predictors[n],names(DAT))]))
    #   }
    #   NEW=cbind(NEW,Add2)
    #   
    #   for(n in 1:length(stOrE))
    #   {
    #     Cls=class(Add[,n])
    #     if(Cls=='numeric') NEW[,2+n]=as.numeric(as.character(NEW[,2+n]))
    #   }
    #   
    #   PREDs=predict(MOD,newdata=NEW,type='response',se.fit=T)
    #   
    #   NEW$PRed=PREDs$fit
    #   NEW$PRed.SE=PREDs$se.fit
    #   
    #   NEW$Dummy=paste(NEW$Species,NEW$zone,sep="_")
    #   ID=match(c("Thickskin_West","Thickskin_Closed.metro","Thickskin_Zone1","Thickskin_Zone2",
    #              "Gummy_North","Gummy_Closed.ningaloo","Gummy_West",
    #              "Whiskery_North","Whiskery_Closed.ningaloo","Whiskery_West","Whiskery_Closed.metro"),NEW$Dummy)
    #   NEW$PRed[ID]=NEW$PRed.SE[ID]=NA
    #   
    #   
    #   
    #   if(PLOT=="points")
    #   {
    #     X=1:nrow(NEW)
    #     Y=PREDs$fit 
    #     Y.SD=PREDs$se.fit
    #     
    #     plot(1:nrow(NEW),PREDs$fit,ylab="",xlab="",xaxt="n",pch=19,cex=2,cex.axis=1.5,
    #          ylim=c(min(Y-Y.SD),max(Y+Y.SD)))
    #     segments(X,Y,X,Y+Y.SD,col=1,lwd=2)
    #     segments(X,Y,X,Y-Y.SD,col=1,lwd=2)
    #     axis(1,1:nrow(NEW),Axis.Labl,cex.axis=1.75)
    #   }
    #   
    #   if(PLOT=="bars")
    #   {
    #     SP=c("Dusky","Thickskin","Gummy","Whiskery")
    #     names(SP)=c("Dusky","Sandbar","Gummy","Whiskery")
    #     for(X in 1:length(SP))
    #     {
    #       aa=subset(NEW,Species==SP[X])
    #       aa=aa[match(c("North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2"),aa$zone),]
    #       B=barplot(aa$PRed,ylab="",xlab="",xaxt="n",cex.axis=1.5,main="",ylim=c(0,max(NEW$PRed,na.rm=T)))
    #       box()
    #       axis(1,B,F,cex.axis=1.75)
    #       legend("topleft",names(SP)[X],cex=1.6,bty='n')
    #     }
    #     
    #     axis(1,B,Axis.Labl,cex.axis=1.4)
    #     
    #   }
    #   
    # }
    
      #Residency
    fn.plt.glm2=function(MOD,DAT,What,Predictors,PLOT)
    {
      Predictors=Predictors[-match(What,Predictors)]
      NEW1=factor(as.character(unique(DAT[,match(What[1],names(DAT))])),levels=levels(DAT[,match(What[1],names(DAT))]))
      NEW=data.frame(NEW1)
      names(NEW)=What
      
      Add=DAT[,match(Predictors,names(DAT))]
      Add=Add[nrow(NEW),]
      for(nn in 1:nrow(NEW))
      {
        xx=subset(DAT,Species==NEW$Species[nn])
        dummy=xx[,match(Predictors,names(xx))]
        qq=dummy[1,]
        for(n in 1:ncol(dummy))
        {
          ss=dummy[,n]
          if(is.factor(ss))
          {
            TaBla=sort(table(ss))
            qq[,n]=names(TaBla[length(TaBla)])
          }
          if(is.numeric(ss)) qq[,n]=mean(ss)
        }
        Add[nn,]=qq
      }
      NEW=cbind(NEW,Add)
      PREDs=predict(MOD,newdata=NEW,type='response',se.fit=T)
      NEW$PRed=PREDs$fit
      NEW$PRed.SE=PREDs$se.fit
      
      NEW$Species=as.character(NEW$Species)
      NEW$Species=with(NEW,ifelse(Species=="Thickskin","Sandbar",Species))
      
      Allsp=c("Dusky","Sandbar","Gummy","Whiskery")
      Msn=Allsp[which(!Allsp%in%NEW$Species)]
      AdD=as.data.frame(matrix(rep(NA,ncol(NEW)*length(Msn)),ncol=ncol(NEW)))
      names(AdD)=names(NEW)
      AdD$Species=Msn
      NEW=rbind(NEW,AdD)
      NEW=NEW[order(NEW$Species),]
      
      
      X=1:nrow(NEW)
      Y=NEW$PRed 
      Y.SD=2*NEW$PRed.SE
      
      if(PLOT=="points")
      {
        plot(1:nrow(NEW),PREDs$fit,ylab="",xlab="",xaxt="n",pch=19,cex=2,cex.axis=1.5,
             ylim=c(0,1))
        segments(X,Y,X,Y+Y.SD,col=1,lwd=2)
        segments(X,Y,X,Y-Y.SD,col=1,lwd=2)
        axis(1,1:nrow(NEW),NEW$Species,cex.axis=1.75)
      }
      
      if(PLOT=="bars")
      {
        a=barplot(Y,ylab="",xlab="",xaxt="n",pch=19,cex=2,cex.axis=1.5,ylim=c(0,1))
        segments(a,Y,a,Y+Y.SD,col=1,lwd=2)
        segments(a,Y,a,Y-Y.SD,col=1,lwd=2)
        axis(1,a,NEW$Species,cex.axis=1.5)
        legend("topleft",names(Res.zones)[z],cex=1.5,bty='n')
        
      }
      box()
    }
    
    tiff(file="GLM.Residency_bars.tiff",width = 1800, height = 2400,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(5,1),mar=c(2,3.75,.1,.1),oma=c(.7,2,1,0.1),mgp=c(2, 0.75, 0),las=1)
    for(z in 1:length(Res.zones))  fn.plt.glm2(MOD=Residency.model[[z]]$GLM,DAT=Residency.model[[z]]$DATA,
                                               What=c("Species"),Predictors=PREDS.res,PLOT="bars")
    mtext("Residency (±SE)",2,0,cex=1.5,las=3,outer=T)
    dev.off()
    

     #Roaming
    fn.plt.glm=function(MOD,DAT,What,Predictors,Axis.Labl,PLOT)
    {
      Predictors=Predictors[-match(What,Predictors)]
      NEW1=factor(as.character(unique(DAT[,match(What[1],names(DAT))])),levels=levels(DAT[,match(What[1],names(DAT))]))
      NEW=data.frame(NEW1)
      names(NEW)=What
      
      Add=DAT[,match(Predictors,names(DAT))]
      Add=Add[nrow(NEW),]
      for(nn in 1:nrow(NEW))
      {
        xx=subset(DAT,Species==NEW$Species[nn])
        dummy=xx[,match(Predictors,names(xx))]
        qq=dummy[1,]
        for(n in 1:ncol(dummy))
        {
          ss=dummy[,n]
          if(is.factor(ss))
          {
            TaBla=sort(table(ss))
            qq[,n]=names(TaBla[length(TaBla)])
          }
          if(is.numeric(ss)) qq[,n]=mean(ss)
        }
        Add[nn,]=qq
      }
      NEW=cbind(NEW,Add)
      PREDs=predict(MOD,newdata=NEW,type='response',se.fit=T)
      NEW$PRed=PREDs$fit
      NEW$PRed.SE=PREDs$se.fit
      
      NEW$Species=as.character(NEW$Species)
      NEW$Species=with(NEW,ifelse(Species=="Thickskin","Sandbar",Species))
      
      Allsp=c("Dusky","Sandbar","Gummy","Whiskery")
      Msn=Allsp[which(!Allsp%in%NEW$Species)]
      AdD=as.data.frame(matrix(rep(NA,ncol(NEW)*length(Msn)),ncol=ncol(NEW)))
      names(AdD)=names(NEW)
      AdD$Species=Msn
      NEW=rbind(NEW,AdD)
      NEW=NEW[order(NEW$Species),]
      
      
      X=1:nrow(NEW)
      Y=NEW$PRed 
      Y.SD=2*NEW$PRed.SE
      
      if(PLOT=="points")
      {
        plot(1:nrow(NEW),PREDs$fit,ylab="",xlab="",xaxt="n",pch=19,cex=2,cex.axis=1.5,
             ylim=c(0,1.2*max(Y+Y.SD)))
        segments(X,Y,X,Y+Y.SD,col=1,lwd=2)
        segments(X,Y,X,Y-Y.SD,col=1,lwd=2)
        axis(1,1:nrow(NEW),NEW$Species,cex.axis=1.75)
      }
      
      if(PLOT=="bars")
      {
        a=barplot(Y,ylab="",xlab="",xaxt="n",pch=19,cex=2,cex.axis=1.5,ylim=c(0,1.2*max(Y+Y.SD)))
        segments(a,Y,a,Y+Y.SD,col=1,lwd=2)
        segments(a,Y,a,Y-Y.SD,col=1,lwd=2)
        axis(1,a,NEW$Species,cex.axis=1.5)
      }
      box()
    }
    
    tiff(file="GLM.roaming_bars.tiff",width = 1800, height = 2400,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(1,1),mar=c(2,0,.1,.1),oma=c(.7,4,1,0.1),mgp=c(2, 0.75, 0),las=1)
    fn.plt.glm(MOD=Roaming.model$GLM,DAT=Roaming.model$DATA,What="Species",Predictors=PREDS,
               Axis.Labl=c("Dusky","Sandbar","Gummy","Whiskery"),PLOT="bars")
    mtext("Roaming (±SE)",2,2.7,cex=1.75,las=3,outer=T)
    dev.off()
    
  }

}

MEPS.paper(do.MEPS.outptus="NO")



#29. -- Fisheries Oceanography figures ---
FO.paper=function(do.FO.outputs)
{
  if(do.FO.outputs=="YES")
  {
    setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement/FO_paper"))
    
    #Figure 1. Map
    #define coordinates of polygons
    N.WA.long=c(North.WA.long[2], North.WA.long[2], North.WA.long[1], North.WA.long[1])
    N.WA.lat=c(North.WA.lat[2], North.WA.lat[1], North.WA.lat[1], North.WA.lat[2])
    S.WA.long=c(South.WA.long[2], South.WA.long[2], South.WA.long[1], South.WA.long[1])
    S.WA.lat=c(South.WA.lat[2], South.WA.lat[1], South.WA.lat[1], South.WA.lat[2])
    
    #define Perth and Rotnest
    Perth=c(115.866,-31.95)
    Rotnest=c(115.50,-32.02)
    
    #bathymetry
    Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
    xbat=sort(unique(Bathymetry$V1))
    ybat=sort(unique(Bathymetry$V2)) 
    if(!exists("reshaped"))
    {
      reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
      rm(Bathymetry)
    }
    
    #legends
    Letter.leg=c("A",NA,"B")
    Letter.leg.coor=cbind(c(113.61,NA,114.45),c(-21.61,NA,-31.85))
    m <- rbind(c(0.0385, 0.3, 0, 1),
               c(0.35, 1, 0.5, 1),
               c(0.3, 1, 0, 0.5))
    
    tiff(file="Figure1.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    split.screen(m)
    for(i in 1:3)
    {
      screen(i)
      par(mar = c(0, 0, 0, 0),mgp=c(.1, 0.15, 0))
      if(i==2) plotlong[[i]][1]=109
      plotMap(worldLLhigh, xlim=plotlong[[i]],ylim=plotlat[[i]],plt = c(.001, 1, 0.075, 1),
              col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
      if(i==1)
      {
        text(114,-22.6,"Ningaloo",cex=1.4)
        text(114,-22.7,"Reef",cex=1.4)
        text(114,-22.8,"array",cex=1.4)        
        # text(113.94,-23.134,"Coral Bay",cex=1)
        text(114.04,-22.03,"North West Cape",srt=75,cex=1)
        
      }
      
      if(i==2)
      {
        #add zones
        # plot(WA_Northern_Shark,ylim=c(-39,-11),xlim=c(108,130),add=T,col="grey90",border="grey90")
        # plot(JA_Northern_Shark,add=T,,col="grey90",border="grey90")
        # plot(WA_Northern_Shark_2,add=T,,col="grey90",border="grey90")
        # plot(WCDGDLL,add=T,,col="grey90",border="grey90")
        # plot(SDGDLL_zone1,add=T,,col="grey90",border="grey90")
        # plot(SDGDLL_zone2,add=T,,col="grey90",border="grey90")
        text(117.65,-25.5,"Shark Bay",cex=.9)
        
        segments(113,-40,124.75,-40,lwd=4)
        text(130,-40,"1000 km",cex=1.5)
        #library(SDMTools)#add scale bar
        #Scalebar(x=113,y=123.6,distance=1000)
        
        
        
        polygon(x=N.WA.long,y=N.WA.lat,lwd=2)
        polygon(x=S.WA.long,y=S.WA.lat,lwd=2)
        text(133,-25,("Australia"),col="black", cex=2)
        mtext("Latitude (ºS)",side=2,line=0,las=3,cex=1.5)
        mtext("Longitude (ºE)",side=1,line=0,cex=1.5)
        text(115.98,-22.6,("A"),col="black", cex=1.75)
        text(119,-33,("B"),col="black", cex=1.75)
        
        text(120.25,-29.5,"Metro")
        arrows(117.75,-30,115.68,-31.2,lwd=1.5,length=0.1)
      }
      if(!i==2)
      {
        points(Receiverlong[[i]],Receiverlat[[i]],col='black',pch=19,cex=0.6)  #receiver location
        axis(side = 1, at =Long.seq[[i]], labels = Long.seq[[i]], tcl = .5,las=1,cex.axis=0.9)
        axis(side = 2, at = Lat.seq[[i]], labels = -Lat.seq[[i]],tcl = .5,las=2,cex.axis=0.9)
        box(lwd=2)
        contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=plotlat[[i]],xlim=plotlong[[i]], zlim=c(-1,-300),
                nlevels = 3,labcex=1.25,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
        text(Letter.leg.coor[i,1],Letter.leg.coor[i,2],Letter.leg[i],cex=2.5)
      }
      
      if(i==3)
      { 
        #polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="red")
        points(115.5,-32.0,pch=19,col='dark grey',cex=1.125)
        text(115,-32.75,"Rottnest",cex=1,srt=60)
        text(115.2,-32.9,"Island",cex=1,srt=60)
        arrows(115.4078,-32.3714,115.4932,-32.129,lwd=1.5,length=0.1)
        #polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
        text(117.75,Perth[2],"Perth array",cex=1.4)
        text(117.1,-34.2,"Southern",cex=1.4)
        text(117.15,-34.5,"Lines array",cex=1.4)
        
        text(115.85,-33.8,"Cape Leeuwin",cex=1,srt=30)
        text(116.17,-31.95,"Perth",cex=1)
        text(117.88,-34.85,"Albany",cex=1)
      }
      
      
    }
    close.screen(all = TRUE)
    dev.off()
    
    
    
    
  }
}
FO.paper(do.FO.outputs="NO")









  
 
# #MOVE TO END
# ##################################################
#dunno what this for....
# fn.data.mov.rat=function(DATA)
# {
#   DATA=DATA[order(DATA$TagCode,DATA$DateTime.local),]
#   unicas=unique(DATA$TagCode)
#   LISTA=vector("list",length=length(unicas))
#   for(j in 1:length(unicas))
#   {
#     Data=subset(DATA,TagCode==unicas[j])
#     Data$Same=c(1,abs(diff(Data$SerialNumber)))    
#     Data=subset(Data,Same>0)
#     Data=Data[,-match("Same",names(Data))]
#     
#     if(nrow(Data)>1)
#     {
#       STOP=nrow(Data)-1
#       Data$Latitude2=c(Data$ReleaseLatitude[1],Data$Latitude[1:STOP])
#       Data$Longitude2=c(Data$ReleaseLongitude[1],Data$Longitude[1:STOP])
#       Data$Date2=c(Data$ReleaseDate[1],Data$DateTime.local[1:STOP])
#       
#       #Hours at liberty
#       Data$Hours.at.large=as.numeric((Data$DateTime.local-Data$Date2)/(60*60))
#       #minimum distance moved following the shortest path (in km). 
#       Data$dist.trav_km=with(Data,
#                               ifelse(Longitude>=115.15 & Longitude<116.425 & Longitude2<115.15,
#                                    (distCosine(cbind(Longitude,Latitude),Cape.Leuwin)+distCosine(Cape.Leuwin,Shark.bay)+
#                                    distCosine(Shark.bay,cbind(Longitude2,Latitude2)))/1000,
#                               ifelse(Longitude>=116.425 & Longitude2<115.15,
#                                   (distCosine(cbind(Longitude,Latitude),Mid.point)+distCosine(Mid.point,Cape.Leuwin)+
#                                   distCosine(Cape.Leuwin,Shark.bay)+distCosine(Shark.bay,cbind(Longitude2,Latitude2)))/1000,
#                               distCosine(cbind(Longitude,Latitude),cbind(Longitude2,Latitude2))/1000)))
#       Data$dist.trav_km=Data$dist.trav_km-(2*detec.range/1000)
#       #speed
#       Data$speed_km_h=with(Data,(dist.trav_km/Hours.at.large))
#     }
#     
#     if(nrow(Data)<=1)
#     {
#       Data=as.data.frame(matrix(nrow=1,ncol=ncol(LISTA[[j-1]])))
#       colnames(Data)=colnames(LISTA[[j-1]])
#     }
#     LISTA[[j]]=Data
#   }
#   LISTA=do.call(rbind,LISTA)
#   id=which(is.na(LISTA$speed_km_h))
#   id=c(id,which(LISTA$Hours.at.large<=(Delta.t/60)))
#   id=c(id,which(LISTA$dist.trav_km<=0))
#   LISTA=LISTA[-id,]
#   return(LISTA)  
# }
# for(i in 1:N.sp)
#   {
#     Data=fn.data.mov.rat(Detections.species[[i]])
# 
#     SPECIES[i]
#   }


# Figure 7. Daily Patterns.       THIS IS WRONG, NOT RIGHT WAY OF DISPLAYING, SEE TIME
#                                 SERIES ANALYSIS!!!!!!


# Detections$hit=ifelse(is.na(Detections$Station),0,1)    #create dummy hit column
# hours(test[2])
# minutes(test[2])
# Detections$HrsMins=hours(Detections$Time.local)+(minutes(Detections$Time.local)/60)  #create hours an minutes
# #Detections$HrsMins=ifelse(Detections$HrsMins==0,NA,Detections$HrsMins)
# 
# #create range of hours
# horas=1:23
# horas.range=vector()
# for (i in 0:length(horas-2)) horas.range=c(horas.range,paste(i,":","00","-",i+1,":","00",sep=""))
# 
# #daily hits
# plot.daily.hits=function(Detections,ESPECIES)
# {
#   DATA=Detections
#   Sharks=unique(DATA$TagCode)
#   table.ID.by.time=table(floor(DATA$HrsMins))
#   
#   a=0:23                                             #add dummy for missing hours
#   if(length(table.ID.by.time)<length(a))
#   {
#     non.detec.rec=match(colnames(table.ID.by.time),a)
#     non.detec.rec=a[-non.detec.rec[!is.na(non.detec.rec)]]
#     non.detec.rec=non.detec.rec[!is.na(non.detec.rec)]
#     table.non.detec=as.table(matrix(0,nrow(table.ID.by.time),length(non.detec.rec)))
#     colnames(table.non.detec)=non.detec.rec;rownames(table.non.detec)=rownames(table.ID.by.time)
#     table.ID.by.time=cbind(table.ID.by.time,table.non.detec)
#     sortnames=a
#     if(nrow(table.ID.by.time)>1)
#     {
#       table.ID.by.time=table.ID.by.time[,match(sortnames,colnames(table.ID.by.time))]
#     }
#     #table.ID.by.time=table.ID.by.time[!(rownames(table.ID.by.time)=="       "),]
#   }
#   Time.names=horas.range
#   
#   #plot
#   maxY= round(max(table.ID.by.time)*1.1)
#   tiff(file=paste("Outputs_movement/",ESPECIES,".Figure4.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
#   bp=barplot(table.ID.by.time,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="",
#              main="",font.main=1,ylim=c(0,maxY),las=2,cex.names=1.5,cex.lab=1.5,cex.axis=1.5,axis.lty=4)
#   axis(side = 1, at = bp, labels = Time.names, tcl = -0.7,las=2,cex.axis=1.3)
#   mtext("Detection numbers",side=2,outer=T,line=-1,font=1,cex=1.5)
#   mtext("Time",side=1,outer=T,line=5.5,font=1,cex=1.5)
#   legend("topright",paste("n=",length(Sharks)," sharks",sep=""),bty="n",cex=1.5)
#   box()
#   dev.off()
# }
# 
# for (i in 1:N.sp) plot.daily.hits(Detections.species[[i]],SPECIES[i])
# 
# par(mfrow=c(2,5),mai=c(0.2,0.2,0,0))
# for(i in 1:length(this))
# {
#   DATA=subset(Detections,TagCode==this[i])
#   #Sharks=unique(DATA$TagCode)
#   table.ID.by.time=table(floor(DATA$HrsMins))
#   
#   a=0:23                                             #add dummy for missing hours
#   if(length(table.ID.by.time)<length(a))
#   {
#     non.detec.rec=match(as.numeric(names(table.ID.by.time)),a)
#     non.detec.rec=a[-non.detec.rec[!is.na(non.detec.rec)]]
#     non.detec.rec.names=non.detec.rec[!is.na(non.detec.rec)]
#     non.detec.rec=rep(0,length(non.detec.rec))
#     names(non.detec.rec)=non.detec.rec.names
#     table.ID.by.time=c(table.ID.by.time,non.detec.rec)
#     table.ID.by.time=table.ID.by.time[match(a,names(table.ID.by.time))]
#   }
#   Time.names=horas.range
#   
#   #plot
#   maxY= round(max(table.ID.by.time)*1.1)
#   #tiff(file=paste("Outputs_movement/",ESPECIES,".Figure4.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
#   bp=barplot(table.ID.by.time,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="",
#              main="",font.main=1,ylim=c(0,maxY),las=2,cex.names=1.5,cex.lab=1.5,cex.axis=1.5,axis.lty=4)
#   axis(side = 1, at = bp, labels = Time.names, tcl = -0.7,las=2,cex.axis=1.3)
#   mtext("Detection numbers",side=2,outer=T,line=-1,font=1,cex=1.5)
#   mtext("Time",side=1,outer=T,line=5.5,font=1,cex=1.5)
#   legend("topright",paste("n=",length(Sharks)," sharks",sep=""),bty="n",cex=1.5)
#   box()
#   #dev.off()
#   
# }

# fn.bubble.rec=function(Detections,ESPECIES1,ESPECIES2,ESPECIES3,ESPECIES4,ARRAY,PlotlonG,PlotlatT,
#                        scale,LETR,Plus.Lon,Plus.Lat,Plus.Lon1,Plus.Lat1,COLOR,Adjst.sym,Adjst.leg,
#                        heightAdj)
# {
#   DATA=subset(Detections,Species==ESPECIES1 & Array==ARRAY)
#   DATA2=subset(Detections,Species==ESPECIES2 & Array==ARRAY)
#   DATA3=subset(Detections,Species==ESPECIES3 & Array==ARRAY)
#   DATA4=subset(Detections,Species==ESPECIES4 & Array==ARRAY)
#   
#   N=N2=N3=N4=0
#   
#   TABLE=table(DATA$Station)
#   TABLE2=table(DATA2$Station)
#   TABLE3=table(DATA3$Station)
#   TABLE4=table(DATA4$Station)
#   
#   MATRIX=matrix(names(TABLE))
#   MATRIX=matrix(as.numeric(unlist(strsplit(MATRIX, split=" "))),ncol=2,byrow=T)
#   LAT=MATRIX[,1]
#   LONG=MATRIX[,2]
#   zo=TABLE/sum(TABLE,na.rm=T)*scale
#   N=sum(TABLE)
#   radius <- sqrt( zo/ pi )
#   
#   if(nrow(TABLE2)>0)
#   {
#     MATRIX2=matrix(names(TABLE2))
#     MATRIX2=matrix(as.numeric(unlist(strsplit(MATRIX2, split=" "))),ncol=2,byrow=T)
#     LAT2=MATRIX2[,1]
#     LONG2=MATRIX2[,2]
#     zo2=TABLE2/sum(TABLE2,na.rm=T)*scale
#     N2=sum(TABLE2)
#     radius2 <- sqrt( zo2/ pi )
#   }
# 
#   if(nrow(TABLE3)>0)
#   {
#     MATRIX3=matrix(names(TABLE3))
#     MATRIX3=matrix(as.numeric(unlist(strsplit(MATRIX3, split=" "))),ncol=2,byrow=T)
#     LAT3=MATRIX3[,1]
#     LONG3=MATRIX3[,2]
#     zo3=TABLE3/sum(TABLE3,na.rm=T)*scale
#     N3=sum(TABLE3)
#     radius3 <- sqrt( zo3/ pi )
#   }
#   if(nrow(TABLE4)>0)
#   {
#     MATRIX4=matrix(names(TABLE4))
#     MATRIX4=matrix(as.numeric(unlist(strsplit(MATRIX4, split=" "))),ncol=2,byrow=T)
#     LAT4=MATRIX4[,1]
#     LONG4=MATRIX4[,2]
#     zo4=TABLE4/sum(TABLE4,na.rm=T)*scale
#     N4=sum(TABLE4)
#     radius4 <- sqrt( zo4/ pi )
#   }
#   
#   
#   plotMap(worldLLhigh, xlim=PlotlonG,ylim=PlotlatT,plt = c(.001, 1, 0.075, 1),
#           col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
#   contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=PlotlatT,xlim=PlotlonG, zlim=c(-1,-100),
#           nlevels = 2,labcex=1.25,lty = 1,col=c(COLOR,COLOR,"transparent"),add=T)  
# 
#   if(!(LETR=="2"))
#   {
#     symbols(LONG-Plus.Lon, LAT*(1-Plus.Lat), circles=radius, inches=0.25, fg="blue", bg="#0000ff22",add=T,
#           xlim=PlotlonG,ylim=PlotlatT,lwd=2)
#     if(nrow(TABLE2)>0)
#     {
#       symbols(LONG2+Plus.Lon, LAT2*(1+Plus.Lat), circles=radius2, inches=0.25, fg="red", bg="#8B000022",add=T,
#               xlim=PlotlonG,ylim=PlotlatT,lwd=2) 
#     }
#     if(nrow(TABLE3)>0)
#     {
#       symbols(LONG3-Plus.Lon1, LAT3*(1-Plus.Lat1), circles=radius3, inches=0.25, fg="green", bg="#0000ff39",add=T,
#               xlim=PlotlonG,ylim=PlotlatT,lwd=2) 
#     }    
#     if(nrow(TABLE4)>0)
#       {
#         symbols(LONG4+Plus.Lon1, LAT4*(1+Plus.Lat1), circles=radius4, inches=0.25, fg="brown", bg="#8B000032",add=T,
#               xlim=PlotlonG,ylim=PlotlatT,lwd=2)
#       }
#   }
# 
#   if(LETR=="2")
#   {
#     symbols(LONG-Plus.Lon, LAT*(1-Plus.Lat), circles=radius, inches=0.25, fg="blue", bg="#0000ff22",add=T,
#             xlim=PlotlonG,ylim=PlotlatT,lwd=2)
#     if(nrow(TABLE2)>0)
#     {
#       symbols(LONG2+Plus.Lon, LAT2*(1+Plus.Lat), circles=radius2, inches=0.25, fg="red", bg="#8B000022",add=T,
#               xlim=PlotlonG,ylim=PlotlatT,lwd=2) 
#     }
#     if(nrow(TABLE3)>0)
#     {
#       symbols(LONG3-Plus.Lon1, LAT3*(1-Plus.Lat1), circles=radius3, inches=0.25, fg="green", bg="#0000ff39",add=T,
#               xlim=PlotlonG,ylim=PlotlatT,lwd=2) 
#     }
#     if(nrow(TABLE4)>0)
#     {
#       symbols(LONG4+Plus.Lon1, LAT4*(1+Plus.Lat1), circles=radius4, inches=0.25, fg="brown", bg="#8B000032",add=T,
#               xlim=PlotlonG,ylim=PlotlatT,lwd=2)
#     }
#     
#     #create legend
#     legPop <- c(5,10,15)
#     legRad <- sqrt(legPop/pi)
#     hin <- par('pin')[2]
#     ylim <- PlotlatT
#     burgPerInch <- ( ylim[2] - ylim[1] ) / hin
#     radPerInch <- max(radius)/0.25
#     #heightAdj <- legRad/radPerInch*burgPerInch
#     #heightAdj= c(0.006114971,0.008647875,0.010591441)
#     Here.X=rep(PlotlonG[1]*Adjst.sym,3)
#     Here.Y=rep(PlotlatT[2]*Adjst.sym,3)
#     
#     symbols(Here.X , Here.Y - heightAdj,circles = legRad,inches = 0.25, add = TRUE,lwd=2)
#     text(Here.X*Adjst.leg, Here.Y- heightAdj, c('5%', '10%', '15%'), cex = 1.1)
#         
#     legend("bottomleft",legend=c(paste("dusky (n= ",N," detections)",sep=""),paste("sandbar (n= ",N2," detections)",sep="")
#                                  ,paste("gummy (n= ",N3," detections)",sep=""),paste("whiskery (n= ",N4," detections)",sep="")),
#            pch=21,col=c("blue","red","green","brown"),pt.bg=c("#0000ff22","#8B000022","#0000ff39","#8B000032") ,pt.lwd=2,bty='n',cex=1.25,
#            pt.cex=2)
#   }
#   points(STATIONS$longitude,STATIONS$latitude,col=1,pch=19)
#   box(lwd=2)
# 
#   legend("topright",LETR,cex=3,bty='n')
# }
# #3.4.1 Ningaloo
# tiff(file=paste("Outputs_movement/Figure4.Ningaloo.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
# par(mar = c(3, 3, 0, 0))
# m <- cbind(c(2, 3,4),c(1, 1,1))
# layout(m)
# #layout.show(4)
# 
# fn.main.map(XLIM=plotlong[[1]],YLIM=c(-23.25,-21.75),LONGSEQ=Long.seq[[1]],LATSEQ=Lat.seq[[1]],
#             Xtext=rbind(114,113.7,113.85),Ytext=rbind(-21.9,-22.6,-23.1),TEXT=c(1:3))
# text(114.128,-21.93,"Exmouth",cex=1.75)
# for (j in 1:3)
# {
#   fn.bubble.rec(Detections,SPECIES[1],SPECIES[2],SPECIES[3],SPECIES[4],"Ningaloo",plotlong.lines[[j]],
#                 plotlat.lines[[j]],5,LETRA[[j]],0,6e-04,0,8e-04,COLOR,1.0001,1.000009,c(0.0015,0.0131,0.0295))
# } 
# dev.off()
# 
# 
# 
# #3.4.2 Perth
# par(def.par)
# tiff(file=paste("Outputs_movement/Figure4.Perth.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
# fn.bubble.rec(Detections,SPECIES[1],SPECIES[2],SPECIES[3],SPECIES[4],"Perth",c(115,116),
#               c(-32.5,-31.5),5,"2",0,6e-04,0,8e-04,COLOR,1.0005,1.00065,c(0.01,0.073,0.16))
# polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
# polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
# text(115.86,-31.95,"Perth",cex=1.75)
# 
# dev.off()
# 
# 
# #3.4.3 Southern Lines
# par(def.par)
# tiff(file=paste("Outputs_movement/Figure4.South.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
# par(mar = c(3, 3, 0, 0))
# m <- rbind(c(1, 1,1),c(2, 3,4))
# layout(m)
# #layout.show(4)
# fn.main.map(XLIM=c(114.5,119.2),YLIM=c(-35.5,-33.5),LONGSEQ=seq(115,119,by=.5),
#             LATSEQ=seq(-35.5,-33.5,by=.5),Xtext=rbind(115.1,116.5,118.32),
#             Ytext=rbind(-34.2,-34.75,-34.78),TEXT=c(1:3))
# text(115.15,-34.36,"Augusta",cex=1.75)
# text(117.88,-35.01,"Albany",cex=1.75)
# 
# 
# long.line=rbind(c(114.58,115.11),c(116.25,116.75),c(118.25,118.75))
# lat.line=rbind(c(-34.41,-34.18),c(-35.5,-35),c(-35.35,-34.85))
# Plus.LONG=c(0,2.5e-02,1e-02)
# Plus.LAT=c(6e-04,0,0)
# Plus.LONG1=c(0,4e-02,5e-02)
# Plus.LAT1=c(1.5e-03,0,0)
# for (j in 1:3)
# {
#   fn.bubble.rec(Detections,SPECIES[1],SPECIES[2],SPECIES[3],SPECIES[4],"Southern.lines",long.line[j,],
#                 lat.line[j,],5,LETRA[[j]],Plus.LONG[j],Plus.LAT[j],Plus.LONG1[j],
#                 Plus.LAT1[j],COLOR,1.0039,0.99999,c(0,0.063,0.13))
# } 
# dev.off()
# par(def.par)




# 
# #1. Manipulate Perth receivers
# SMN_VR2_2010$Lat=ifelse(SMN_VR2_2010$Lat>0,-SMN_VR2_2010$Lat,SMN_VR2_2010$Lat)
# thesecolumns=match(c("Station No","Lat","Long","Serial No"),names(SMN_VR2_2010))
# SMN_VR2_2010=SMN_VR2_2010[,thesecolumns]
# thesecolumnsVR4=match(c("Station No","Lat","Long","VR4G Serial No"),names(SMN_VR4G_2009))
# thesecolumnsOTN09a=match(c("Station No","Lat","Long","Serial No"),names(OTN_2009a))
# SMN_VR4G_2009=SMN_VR4G_2009[,thesecolumnsVR4]
# SMN_VR4G_2009=SMN_VR4G_2009[!(is.na(SMN_VR4G_2009$Lat)),] #remove NA rows
# OTN_2009a=OTN_2009a[,thesecolumnsOTN09a]
# names(SMN_VR4G_2009)=names(OTN_2009a)=names(SMN_VR2_2010)
# SMN_VR2_2010$RecType="VR2.SMN"
# SMN_VR4G_2009$RecType="VR4.SMN"
# OTN_2009a$RecType="VR2.OTN"
# Perth=rbind(SMN_VR2_2010,OTN_2009a,SMN_VR4G_2009)
# Perth=Perth[!(is.na(Perth$Lat)),]
# 
# #2. Manipulate tagging data and convert UTC detection time to Perth local time     (+8hours)
# Tagging.Data.Ningaloo$DS="Ningaloo"
# colnames(Tagging.Data.Ningaloo)[1]="Date.and.Time..UTC."
# Tagging.Data.Ningaloo$Date.time=as.POSIXlt(as.POSIXct(as.character(Tagging.Data.Ningaloo$"Date.and.Time..UTC."),
# tz="GMT"), "Australia/Perth")
# Tagging.Data.Ningaloo$Date=as.Date(Tagging.Data.Ningaloo$Date.time)
# 
# Tagging.Data.Perth$DS="Perth"
# Tagging.Data.Perth$Date.time=strptime(as.character(Tagging.Data.Perth$"Date.and.Time..UTC."),format='%d/%m/%Y %H:%M')
# Tagging.Data.Perth$Date.time=as.POSIXlt(as.POSIXct(Tagging.Data.Perth$Date.time,tz="GMT"), "Australia/Perth")
# Tagging.Data.Perth$Date=as.Date(Tagging.Data.Perth$Date.time)
# 
# Tagging.Data.SouthWest$DS="SouthWest"
# Tagging.Data.SouthWest$Date.time=strptime(as.character(Tagging.Data.SouthWest$"Date.and.Time..UTC."),format='%d/%m/%Y %H:%M')
# Tagging.Data.SouthWest$Date.time=as.POSIXlt(as.POSIXct(Tagging.Data.SouthWest$Date.time,tz="GMT"), "Australia/Perth")
# Tagging.Data.SouthWest$Date=as.Date(Tagging.Data.SouthWest$Date.time)
# 
# Tagging.Data=rbind(Tagging.Data.Ningaloo,Tagging.Data.Perth,Tagging.Data.SouthWest)
# 
#       #extract tagged individuals by species
# tagged.species=table(Tag.Deployed$Species)
# Species=unique(as.character(Tag.Deployed$Species))
# 
# 
# #3. Extract receiver and tag numbers
# Tagging.Data$Receiver.Number=unlist(strsplit(as.character(Tagging.Data$Receiver), split="-"))[seq(2,2*nrow(Tagging.Data),2)]
# Tagging.Data$Tag=unlist(strsplit(as.character(Tagging.Data$Transmitter), split="-"))[seq(3,3*nrow(Tagging.Data),3)]
# 
# 
# Tagging.Data$hit=ifelse(is.na(Tagging.Data$Receiver.Number)==T,0,1)    #create dummy hit column
# Tagging.Data$HrsMins=as.POSIXlt(Tagging.Data$Date.time)$hour+(as.POSIXlt(Tagging.Data$Date.time)$min/60)  #create hours an minutes
# Tagging.Data$HrsMins=ifelse(Tagging.Data$HrsMins==0,NA,Tagging.Data$HrsMins)
# 
# 
# 
# #4. Merge receiver position and hits                            (FIX FOR REAL DATA FROM PERTH AND SOUTH WEST)
# these.cols.Ning=match(c("Receiver (S/N)","Recovery_latitude","Recovery_longitude","Comments"),names(Ningaloo))
# Ning.Rec=Ningaloo[,these.cols.Ning]
# these.cols.Perth=match(c("Serial No","Lat","Long","RecType"),names(Perth))
# Perth.Rec=Perth[,these.cols.Perth]
# these.cols.SouthWest=match(c("Mooring.Number","Lat","Long","Line"),names(SouthWest))
# SouthWest.Rec=SouthWest[,these.cols.SouthWest]
# names(Perth.Rec)=names(SouthWest.Rec)=names(Ning.Rec)
# 
# All.Rec=rbind(Ning.Rec,Perth.Rec,SouthWest.Rec)
# 
# 
# these.cols.Tagging=match(c("Date.time","Date","Receiver.Number","Tag","DS","hit","HrsMins"),names(Tagging.Data))
# 
# Tagging.Data=merge(Tagging.Data[,these.cols.Tagging],All.Rec,by.y="Receiver (S/N)",
# by.x="Receiver.Number",all.x=T)
# 
# 
# 
# #5. Select relevant tag numbers
# these.tags=unique(Tag.Deployed$ATAG.NO)
# Tagging.Data=subset(Tagging.Data,Tag%in%these.tags)
# 
# #6. Manipulate tagging data
#              Tagging.Data$"Recovery_latitude"=ifelse(Tagging.Data$"Recovery_latitude">0,-Tagging.Data$"Recovery_latitude",Tagging.Data$"Recovery_latitude")
# Tagging.Data=Tagging.Data[!(duplicated(paste(Tagging.Data$Receiver.Number,Tagging.Data$Date.time,Tagging.Data$Tag))),]  #remove duplicates
# 
# maxdate=strptime(max(Tagging.Data$Date.time), "%Y-%m-%d")    #min and max dates of data
# mindate=strptime(min(Tagging.Data$Date.time), "%Y-%m-%d")
# 
# Tagging.Data=Tagging.Data[order(Tagging.Data$Date.time),]     #order by date time
# 
# 
# 
#   #extract tagging data for each species
# Tag.Deployed.SandBar=subset(Tag.Deployed,Species=="TK")
# Tag.Deployed.Dusky=subset(Tag.Deployed,Species=="BW")
# #Tag.Deployed.Gummy=subset(Tag.Deployed,Species=="GM")                                          #UPDATE
# #Tag.Deployed.Whiskery=subset(Tag.Deployed,Species=="WH")
# 
# 
# 
# #7. Expand hits to have all days within studied period
# YEARS=unique(1900 + as.POSIXlt(Tagging.Data$Date.time)$year)
# year=NULL
# getDays <- function(year)
# {
#   seq(as.Date(paste(year, "-01-01", sep="")), as.Date(paste(year, "-12-31", sep="")), by="+1 day")
# }
# alldays=lapply(YEARS,getDays)
# 
# 
# #8. Create useful vectors  and data frames
# Receivers.Ning=unique(Ningaloo$"Receiver (S/N)")
# Receivers.Per=unique(Perth$"Serial No")
# Receivers.SW=unique(SouthWest$Mooring.Number)
# 
# 
#   #add species name
# these.Tag.Deployed=match(c("ATAG.NO","Species"),names(Tag.Deployed))
# Tagging.Data=merge(Tagging.Data,Tag.Deployed[,these.Tag.Deployed],by.x="Tag", by.y="ATAG.NO",all.x=T)
# 
# 
#   #compare tagged vs detected individuals
# Detected.hits.sp=table(Tagging.Data$Species)
# Deployed.sp=table(Tag.Deployed$Species)
# 
# Detected.hits.sp.and.tag=table(Tagging.Data$Species,Tagging.Data$Tag)
# Deployed.sp.and.tag=table(Tag.Deployed$Species,Tag.Deployed$ATAG.NO)
# 
# Detected.hits.numbers.sp.=ifelse(Detected.hits.sp.and.tag>0,1,0)
# Detected.hits.numbers.sp.=rowSums(Detected.hits.numbers.sp.)
# Deployed.numbers.sp.=rowSums(Deployed.sp.and.tag)
# 
# 
#   #create data frames by species
# Tagging.Data.SandBar=subset(Tagging.Data,Species=="TK")
# Tagging.Data.Dusky=subset(Tagging.Data,Species=="BW")
# #Tagging.Data.Gummy=subset(Tagging.Data,Species=="xxxx")                            #UPDATE with gummy and whiskery code
# #Tagging.Data.Whiskery=subset(Tagging.Data,Species=="xxxx")
# 
# 
# Useful.Vecs=function(whichshark)
# {
#   Ning=subset(whichshark,DS=="Ningaloo")
#   Per=subset(whichshark,DS=="Perth")
#   SoWe=subset(whichshark,DS=="SouthWest")
# 
#   by.area=function(AREA)
#   {
#     Sharks=unique(AREA$Tag)     #unique sharks
#     n.shark=length(Sharks)
#     Shark.Times=unique(AREA$Date.time)     #unique date and time
#     Shark.colors=rainbow(n.shark)       #colors for these sharks if detected simultaneously
#     Shark.colors=sample(Shark.colors, n.shark, replace = FALSE)       #make sequential colors different
#     return(list(Sharks=Sharks,n.shark=n.shark,Shark.Times=Shark.Times,Shark.colors=Shark.colors))
#    }
#    Ning=by.area(Ning)
#    Per=by.area(Per)
#    SoWe=by.area(SoWe)
#    return(list(Ningaloo=Ning,Perth=Per,SouthWest=SoWe))
# }
# 
# Useful.Vecs.SandBar=Useful.Vecs(Tagging.Data.SandBar)
# Useful.Vecs.Dusky=Useful.Vecs(Tagging.Data.Dusky)
# #Useful.Vecs.Gummy=Useful.Vecs(Tagging.Data.Gummy)                                                     #UPDATE
# #Useful.Vecs.Whiskery=Useful.Vecs(Tagging.Data.Whiskery)
# 
# 
#   #add dummy shark order for plotting
# Tagging.Data.SandBar$ID.number=match(Tagging.Data.SandBar$Tag,Useful.Vecs.SandBar$Sharks)
# Tagging.Data.Dusky$ID.number=match(Tagging.Data.Dusky$Tag,Useful.Vecs.Dusky$Sharks)
# #Tagging.Data.Gummy$ID.number=match(Tagging.Data.Gummy$Tag,Useful.Vecs.Gummy$Sharks)
# #Tagging.Data.Whiskery$ID.number=match(Tagging.Data.Whiskery$Tag,Useful.Vecs.Whiskery$Sharks)
# 
#     #add missing dates
# add.missing.dates=function(whichspecies)
# {
#   Tagging.Data.byday=whichspecies
#   maxdate=strptime(max(Tagging.Data.byday$Date), "%Y-%m-%d")    #min and max dates of data
#   mindate=strptime(min(Tagging.Data.byday$Date), "%Y-%m-%d")
# 
#   #add continuous dates
#   DATES=alldays[[1]][1]
#   for (i in 1:length(alldays))
#   {
#     DATES=c(DATES,alldays[[i]])
#   }
#   DATES=DATES[-1]
#   matchDates=unique(match(Tagging.Data.byday$Date,DATES))
# 
#   DATES=DATES[-matchDates]
#   DATES=subset(DATES,DATES>=as.Date(mindate))
#   DATES=subset(DATES,DATES<=as.Date(maxdate))
# 
#   dummyDATES=data.frame(Tag="       ",Receiver.Number=NA,Date.time=NA,Date=DATES,DS=NA,hit=NA,HrsMins=NA,
#   "Recovery_latitude"=NA,"Recovery_longitude"=NA,Comments=NA,Species=NA,ID.number=NA)
#   names(dummyDATES)=names(Tagging.Data.byday)
# 
#   Tagging.Data.byday=rbind(Tagging.Data.byday,dummyDATES)
# 
#   Tagging.Data.byday=Tagging.Data.byday[order(Tagging.Data.byday$Date),]
# 
#   return(Tagging.Data.byday)
# }
# 
# Tagging.Data.SandBar=add.missing.dates(Tagging.Data.SandBar)
# Tagging.Data.Dusky=add.missing.dates(Tagging.Data.Dusky)
# ##Tagging.Data.byday.Gummy=add.missing.dates(Tagging.Data.Gummy)                                                     #UPDATE
# ##Tagging.Data.byday.Whiskery=add.missing.dates(Tagging.Data.Whiskery)
# 

# 
# 
# #    #extract hits by day, shark, receiver and data set
# #hits.day=function(whichspecies)
# #{
# #  whichspecies$Hits=1
# #  Tagging.Data.byday=aggregate(Hits ~ Date + DS + Receiver.Number + Tag + Recovery_latitude + Recovery_longitude, data = whichspecies, sum)
# #  maxdate=strptime(max(Tagging.Data.byday$Date), "%Y-%m-%d")    #min and max dates of data
# #  mindate=strptime(min(Tagging.Data.byday$Date), "%Y-%m-%d")
# #
# #  #add continuous dates
# #  DATES=alldays[[1]][1]
# #  for (i in 1:length(alldays))
# #  {
# #    DATES=c(DATES,alldays[[i]])
# #  }
# #  DATES=DATES[-1]
# #  matchDates=unique(match(Tagging.Data.byday$Date,DATES))
# #
# #  DATES=DATES[-matchDates]
# #  DATES=subset(DATES,DATES>=as.Date(mindate))
# #  DATES=subset(DATES,DATES<=as.Date(maxdate))
# #
# #  dummyDATES=data.frame(Date=DATES,DS=NA,Receiver.Number=NA,Tag="       ","Recovery_latitude"=NA,
# #  "Recovery_longitude"=NA,"Hits"=NA)
# #  names(dummyDATES)=names(Tagging.Data.byday)
# #
# #  Tagging.Data.byday=rbind(Tagging.Data.byday,dummyDATES)
# #
# #  Tagging.Data.byday=Tagging.Data.byday[order(Tagging.Data.byday$Date),]
# #
# #  return(Tagging.Data.byday)
# #
# #}
# #
# #Tagging.Data.byday.SandBar=hits.day(Tagging.Data.SandBar)
# #Tagging.Data.byday.Dusky=hits.day(Tagging.Data.Dusky)
# ##Tagging.Data.byday.Gummy=hits.day(Tagging.Data.Gummy)                                                     #UPDATE
# ##Tagging.Data.byday.Whiskery=hits.day(Tagging.Data.Whiskery)
# #
# 
# 
# 
# #9. General bar charts
# 
#     #--hits by day, receivers and time--
# # get major and minor multiples for choosing labels
# ntick.Days=20  #control the number of days between labels
# ntick.Receiver=ntick.Time=2
# MinorMajorTick=function(ntick)
# { if (ntick < 16) mult = c(2, 2)
#     else if(ntick < 41) mult = c(5, 5)
#     else if (ntick < 101) mult = c(10, 5)     else mult = c(20, 5)
#     return(mult=mult)
# }
# 
# #create range of hours
# horas=1:23
# horas.range=vector()
# for (i in 0:length(horas-2)) horas.range=c(horas.range,paste(i,":","00","-",i+1,":","00",sep=""))
# 
# #plotting function
# I.plot.You=function(database,Useful,ActiveFile)
# {
#   Especie=as.character(unique(database$Species))
#   Especie=Especie[!(is.na(Especie))]
#   
#   database.Ningaloo=subset(database,!(Tag=="       ")& DS=="Ningaloo")
#   database.Perth=subset(database,!(Tag=="       ")& DS=="Perth")
#   database.SouthWest=subset(database,!(Tag=="       ")& DS=="SouthWest")
# 
#   Useful.Ningaloo=Useful[[1]]
#   Useful.Perth=Useful[[2]]
#   Useful.SouthWest=Useful[[3]]
# 
# 
# 
#   #hits by day
#   plot.hits.by.day=function(Sub.database,Sub.Useful)
#   {
#       #create table
#     table.ID.by.Date=table(Sub.database$Tag,Sub.database$Date)
#     date.names=colnames(table.ID.by.Date)   #drop some dates for labelling
#     Shark.colors=Sub.Useful$Shark.colors
#     Sharks=Sub.Useful$Sharks
#     thesecolors=Shark.colors[match(rownames(table.ID.by.Date),Sharks)]  #choose colors to keep it consistent among datasets
#        #select ticks
#     mult.Days=MinorMajorTick(ntick.Days)
#     datalabels <- seq(from = 1, along.with = date.names)
#     label.index <- which(datalabels %% mult.Days[1] == 0)
#     minor.index = which(datalabels %% mult.Days[2] == 0) # Draw all ticks:
#       #plot
#     maxY= round(max(colSums(table.ID.by.Date))*1.05)
#     #bp=barplot(table.ID.by.Date,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="",args.legend = list(bty="n",cex=1.1),
#     bp=barplot(table.ID.by.Date,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="",
#     main="",font.main=1,legend=F,
#     #main=paste(unique(Sub.database$DS),"(",ActiveFile,")",sep=""),cex.main=1.0,font.main=1,legend=T,
#     ylim=c(0,maxY),las=2,col=thesecolors,cex.names=1.5,cex.lab=1.5,cex.axis=1.5,axis.lty=4)
#     mtext(paste(unique(Sub.database$DS),"(",ActiveFile,")",sep=""),3,line=-1.5,cex=1.3)
#     axis(side = 1, at = bp, labels = FALSE, tcl = -0.2) # Draw minor ticks:
#     axis(side = 1, at = bp[minor.index], labels = FALSE, tcl = -0.5) # Draw major ticks & labels:
#     axis(side = 1, at = bp[label.index], labels = date.names[label.index], tcl = -0.7,las=2,cex.axis=1.3)
#     if(unique(Sub.database$DS)=="Perth")
#     {mtext("Detection numbers",side=2,outer=T,line=-1,font=1,cex=1.5)}
#     mtext("Date",side=1,outer=T,line=-0.2,font=1,cex=1.5)
#     legend("topright",paste("n=",length(Sharks)," sharks",sep=""),bty="n",cex=1.5)
#     box()
#   }
# 
#   png(file=paste("hits by day.",ActiveFile,".png",sep=""),width=1200,height=900)
#   if(Especie%in%c("TK","BW"))
#   {
#     par(mfcol=c(3,1),mar=c(6,4,1,1), oma=c(1,1,1,1))
#     plot.hits.by.day(database.Ningaloo,Useful.Ningaloo)
#     plot.hits.by.day(database.Perth,Useful.Perth)
#     plot.hits.by.day(database.SouthWest,Useful.SouthWest)
#   }
#   if(Especie%in%c("GM","WH"))
#   {
#     par(mfcol=c(2,1),mar=c(6,4,1,1), oma=c(1,1,1,1))
#     plot.hits.by.day(database.Perth,Useful.Perth)
#     plot.hits.by.day(database.SouthWest,Useful.SouthWest)
# 
#   }
#   dev.off()
# 
# 
#   #hits by receiver
#    plot.hits.by.receiver=function(Sub.database,Sub.Useful,Receiver.Area)
#   {
#      #create table
#     table.ID.by.Receiver=table(Sub.database$Tag,Sub.database$Receiver.Number)
#     a=Receiver.Area
# 
#     if(ncol(table.ID.by.Receiver)<length(a))
#     {
#       non.detec.rec=match(colnames(table.ID.by.Receiver),a)
#       non.detec.rec=a[-non.detec.rec[!is.na(non.detec.rec)]]
#       non.detec.rec=non.detec.rec[!is.na(non.detec.rec)]
#       table.non.detec=as.table(matrix(0,nrow(table.ID.by.Receiver),length(non.detec.rec)))
#       colnames(table.non.detec)=non.detec.rec;rownames(table.non.detec)=rownames(table.ID.by.Receiver)
#       table.ID.by.Receiver=cbind(table.ID.by.Receiver,table.non.detec)
#       sortnames=sort(colnames(table.ID.by.Receiver))
#       if(nrow(table.ID.by.Receiver)>1)
#       {
#         table.ID.by.Receiver=table.ID.by.Receiver[,match(sortnames,colnames(table.ID.by.Receiver))]
#       }
#     }
#     Receiver.names=colnames(table.ID.by.Receiver)   #drop some dates for labelling
#     Shark.colors=Sub.Useful$Shark.colors
#     Sharks=Sub.Useful$Sharks
#     thesecolors=Shark.colors[match(rownames(table.ID.by.Receiver),Sharks)]  #choose colors to keep it consistent among datasets
#        #select ticks
#     mult.Receiver=MinorMajorTick(ntick.Receiver)
#     datalabels <- seq(from = 1, along.with = Receiver.names)
#     label.index <- which(datalabels %% mult.Receiver[1] == 0)
#     minor.index = which(datalabels %% mult.Receiver[2] == 0) # Draw all ticks:
#       #plot
#     maxY= round(max(colSums(table.ID.by.Receiver))*1.05)
#     #bp=barplot(table.ID.by.Receiver,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="",args.legend = list(bty="n",cex=1.1),
#     bp=barplot(table.ID.by.Receiver,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="",
#     #main=paste(unique(Sub.database$DS),"(",ActiveFile,")",sep=""),cex.main=1.0,font.main=1,legend=T,
#     main="",font.main=1.5,legend=F,
#     ylim=c(0,maxY),las=2,col=thesecolors,cex.names=1.5,cex.lab=1.5,cex.axis=1.5,axis.lty=4)
#     mtext(paste(unique(Sub.database$DS),"(",ActiveFile,")",sep=""),3,line=-1.5,cex=1.3)
#     axis(side = 1, at = bp, labels = FALSE, tcl = -0.2) # Draw minor ticks:
#     axis(side = 1, at = bp[minor.index], labels = FALSE, tcl = -0.5) # Draw major ticks & labels:
#     axis(side = 1, at = bp[label.index], labels = Receiver.names[label.index], tcl = -0.7,las=2,cex.axis=1.5)
#     if(unique(Sub.database$DS)=="Perth")
#     {mtext("Detection numbers",side=2,outer=T,line=-1,font=1,cex=1.5)}
#     mtext("Receiver",side=1,outer=T,line=-0.2,font=1,cex=1.5)
#     legend("topright",paste("n=",length(Sharks)," sharks",sep=""),bty="n",cex=1.5)
#     box()
#   }
#   
#   png(file=paste("hits by receiver.",ActiveFile,".png",sep=""),width=1200,height=900)
#   if(Especie%in%c("TK","BW"))
#   {
#     par(mfcol=c(3,1),mar=c(6,4,1,1), oma=c(1,1,1,1))
#     plot.hits.by.receiver(database.Ningaloo,Useful.Ningaloo,Receivers.Ning)
#     plot.hits.by.receiver(database.Perth,Useful.Perth,Receivers.Per)
#     plot.hits.by.receiver(database.SouthWest,Useful.SouthWest,Receivers.SW)
#   }
#   if(Especie%in%c("GM","WH"))
#   {
#     par(mfcol=c(2,1),mar=c(6,4,1,1), oma=c(1,1,1,1))
#     plot.hits.by.receiver(database.Perth,Useful.Perth,Receivers.Per)
#     plot.hits.by.receiver(database.SouthWest,Useful.SouthWest,Receivers.SW)
#    }
#   dev.off()
# 
# 
#  #daily hits
#   plot.daily.hits=function(Sub.database,Sub.Useful)
#   {
#     table.ID.by.time=table(Sub.database$Tag,floor(Sub.database$HrsMins))
# 
#        a=0:23                                             #add dummy for missing hours
#     if(ncol(table.ID.by.time)<length(a))
#     {
#       non.detec.rec=match(colnames(table.ID.by.time),a)
#       non.detec.rec=a[-non.detec.rec[!is.na(non.detec.rec)]]
#       non.detec.rec=non.detec.rec[!is.na(non.detec.rec)]
#       table.non.detec=as.table(matrix(0,nrow(table.ID.by.time),length(non.detec.rec)))
#       colnames(table.non.detec)=non.detec.rec;rownames(table.non.detec)=rownames(table.ID.by.time)
#       table.ID.by.time=cbind(table.ID.by.time,table.non.detec)
#       sortnames=a
#       if(nrow(table.ID.by.time)>1)
#       {
#         table.ID.by.time=table.ID.by.time[,match(sortnames,colnames(table.ID.by.time))]
#       }
#       #table.ID.by.time=table.ID.by.time[!(rownames(table.ID.by.time)=="       "),]
#     }
#     Time.names=horas.range
#     Shark.colors=Sub.Useful$Shark.colors
#     Sharks=Sub.Useful$Sharks
#     thesecolors=Shark.colors[match(rownames(table.ID.by.time),Sharks)]  #choose colors to keep it consistent among datasets
# 
#       #plot
#     maxY= round(max(colSums(table.ID.by.time))*1.05)
#     bp=barplot(table.ID.by.time,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="",
#     #bp=barplot(table.ID.by.time,beside=FALSE,xaxt='n', space = 0,xlab="",ylab="",args.legend = list(bty="n",cex=1.1),
#     main="",font.main=1,
# #    main=paste(unique(Sub.database$DS),"(",ActiveFile,")",sep=""),cex.main=1.0,font.main=1,legend=T,
#     ylim=c(0,maxY),las=2,col=thesecolors,cex.names=1.5,cex.lab=1.5,cex.axis=1.5,axis.lty=4)
#     mtext(paste(unique(Sub.database$DS),"(",ActiveFile,")",sep=""),3,line=-1.5,cex=1.3)
#     if(unique(Sub.database$DS)=="SouthWest")
#     {axis(side = 1, at = bp, labels = Time.names, tcl = -0.7,las=2,cex.axis=1.3)}
#     if(unique(Sub.database$DS)=="Perth")
#     {mtext("Detection numbers",side=2,outer=T,line=-1,font=1,cex=1.5)}
#     mtext("Time",side=1,outer=T,line=5.5,font=1,cex=1.5)
#     legend("topright",paste("n=",length(Sharks)," sharks",sep=""),bty="n",cex=1.5)
#     box()
#   }
# 
#   png(file=paste("hits by time.",ActiveFile,".png",sep=""),width=1200,height=900)
#   if(Especie%in%c("TK","BW"))
#   {
#     par(mfcol=c(3,1),mar=c(2,4,1,1), oma=c(7,1,1,1))
#     plot.daily.hits(database.Ningaloo,Useful.Ningaloo)
#     plot.daily.hits(database.Perth,Useful.Perth)
#     plot.daily.hits(database.SouthWest,Useful.SouthWest)
#   }
#   if(Especie%in%c("GM","WH"))
#   {
#     par(mfcol=c(2,1),mar=c(2,4,1,1), oma=c(7,1,1,1))
#     plot.daily.hits(database.Perth,Useful.Perth)
#     plot.daily.hits(database.SouthWest,Useful.SouthWest)
#    }
#   dev.off()
# }
# 
# # submit functions to export graphs
# I.plot.You(Tagging.Data.SandBar,Useful.Vecs.SandBar,"Sandbar")
# I.plot.You(Tagging.Data.Dusky,Useful.Vecs.Dusky,"Dusky")
# #I.plot.You(Tagging.Data.Gummy,Useful.Vecs.Gummy,"Gummy")                                                      #UPDATE!!!!!!!
# #I.plot.You(Tagging.Data.Whiskery,Useful.Vecs.Whiskery,"Whiskery")
# 
# 
# ##10. Create general map of study site                                          #MISSING MAP!!!!!!!
# #data(worldLLhigh)
# #plotlat=c(-32.4,-31.7)
# #plotlong=c(115.2,115.9)
# #
# #Perth=c(115.866,-31.95)
# #Rotnest=c(115.50,-32.02)
# #
# ##small inset map
# #insetOz <- function()
# #{
# #  opar <- par(mai = c(1.5,2,1,1))
# #  on.exit(par(opar))
# #  par(mar=rep(.1,4),xaxt="n",yaxt="n",plt=par("plt"))
# #  plotMap(worldLLhigh, xlim=c(110,155), ylim=c(-44.5,-11),col="light grey", axes=F, xlab="", ylab="")
# #  text(133,-25,("Australia"),col="black", cex=1)
# #  points(153.2199,-11.53294,col="white",pch=21,bg="white")
# #  edgeX=c(113,117,117,113)
# #  edgeY=c(-29,-29,-32,-32)
# #  polygon(x=edgeX,y=edgeY,lwd=2)
# #  polygon(x=c(112,155,155,112),y=c(-44.5,-44.5,-11,-11),lwd=2)
# #}
# #
# ##large map with stations
# #png(file="map of receivers.png",width=600,height=600)
# ##X11(width=14,height=14)		#control width and height of graphic
# #plotMap(worldLLhigh, xlim=plotlong,ylim=plotlat, plt = c(0.025, .99, 0.1, 0.99),
# #col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",tckLab=F)
# #axis(side = 1, at = c(115.4,115.6,115.8), labels = c("115.4","115.6","115.8"), tcl = 0,las=1,cex.axis=0.9)
# #axis(side = 2, at = c(-31.8,-31.9,-32,-32.1,-32.2,-32.3), labels = c("31.8","31.9","32","32.1","32.2","32.3"),
# #tcl = 0,las=2,cex.axis=0.9)
# #polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
# #polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
# #points(c(SMN_VR2[['2009']]$Long,SMN_VR2[['2010']]$Long),c(SMN_VR2[['2009']]$Lat,SMN_VR2[['2010']]$Lat),cex=1.11,
# #pch=19,col="black")
# #points(SMN_VR4G[['2009']]$Long,SMN_VR4G[['2009']]$Lat,cex=1.11,pch=17,col="gray42")     #plot VR4s
# #points(OTN_VR2[['2009a']]$Long,OTN_VR2[['2009a']]$Lat,cex=1.075,pch=19,col="black")      #plot all VR2s (SMN and OTN)
# #text(115.79,-31.97,("Perth"),col="black", cex=1.15)
# #text(115.51,-31.97, ("Rottnest Is."),col="black", cex=1.15)
# #legend(115.3,-32.3,c("VR2","VR4"),pch=c(19,17),bty="n",col=c("black","gray42"),text.col="black",cex=1.1)
# #mtext("Latitude (ºS)",side=2,outer=T,line=-2.5,font=1,las=0,cex=1)
# #mtext("Longitude (ºE)",side=1,outer=T,line=-2.5,font=1,las=1,cex=1)
# #vp <- baseViewports()
# #pushViewport(vp$inner,vp$figure,vp$plot)
# #pushViewport(viewport(x=0.1,y=0.99,width=.30,height=.30,just=c("left","top")))
# #par(fig=gridFIG(),new=T)
# #insetOz()
# #dev.off()
# ##X11()	#reset with and height to defaut
# #
# #
# 
# 
# 
# #11 Create presence of tagged sharks thru time
# 
# PresenceTagTime=function(database,tagshark,Useful,ActiveFile)
# {
#   database=subset(database,!(Tag=="       "))
# 
#   Especie=as.character(unique(database$Species))
#   Especie=Especie[!(is.na(Especie))]
# 
#   database.Ningaloo=subset(database,DS=="Ningaloo")
#   database.Perth=subset(database,DS=="Perth")
#   database.SouthWest=subset(database,DS=="SouthWest")
# 
#   Useful.Ningaloo=Useful[[1]]
#   Useful.Perth=Useful[[2]]
#   Useful.SouthWest=Useful[[3]]
#   
#   tagshark=tagshark[order(tagshark$ATAG.NO),]
# 
#        #date range
#   daterange=c(min(as.numeric(tagshark$DateDeployed)),max(as.numeric(database$Date)))
#   datadailyticks <- seq(from = min(daterange), max(daterange))
# 
#       #get major and minor multiples for choosing labels
#   #date.names=unique(sort(database$Date))
#   #date.names=c(seq(min(tagshark$Date),date.names[1]-1, "days"),date.names)   #add tagging date
#   date.names=as.Date(datadailyticks, origin="1970-01-01")
# 
#   mult.Days=7       #days in between labels
#   datalabels <- seq(from = min(as.numeric(date.names)), along.with = date.names)
#   label.index = which(datalabels %% mult.Days == 0)
# 
# 
# 
#    plot.detec.thru.time=function(Sub.database,Sub.Useful)
#   {
#     # add non-detected but tagged sharks
#     detectedsharks=match(as.numeric(unique(Sub.database$Tag)),tagshark$ATAG.NO)
#     nondetectedsharks=tagshark$ATAG.NO[-detectedsharks]
#     Dummy.noDetec=data.frame(Tag=nondetectedsharks,Receiver.Number=NA,Date.time=NA,Date=NA,DS=unique(Sub.database$DS),hit=NA,HrsMins=NA,Recovery_latitude=NA,
#     Recovery_longitude=NA,Comments=NA,Species=unique(Sub.database$Species),ID.number=NA)
#     Sub.database=rbind(Sub.database,Dummy.noDetec)
#     Sub.database=Sub.database[order(as.numeric(Sub.database$Tag),Sub.database$Date),]
# 
#     Atshark=sort(unique(Sub.database$Tag[!is.na(Sub.database$Tag)]))
#     Sharks=Sub.Useful$Sharks
# 
#     with(Sub.database,{
#       plot(Date,as.factor(Tag),type="n", ylab="",xlab="",yaxt='n',xaxt='n',cex.main=1,
#       main="",xlim=daterange)
#       points(Date,as.factor(Tag),pch=19,col="black",cex=1.25)
#       axis(side = 2, at =as.factor(Atshark) , labels = Atshark, tcl = -0.5,las=2,cex.axis=0.75)
#       axis(side = 1, at = datadailyticks, labels = FALSE, tcl = -0.2) # Draw minor ticks
#       axis(side = 1, at = datadailyticks[label.index], labels = F,tcl = -0.5) # Draw major ticks & labels
#       mtext(paste(unique(DS),"(",ActiveFile,")",sep=""),3,line=-1.5,cex=1.1)
#       if(unique(DS)=="SouthWest")
#       {
#         axis(side = 1, at = datadailyticks[label.index], labels = date.names[label.index],
#         tcl = -0.5,las=2,cex.axis=0.9) # Draw major ticks & labels
#       }
#       if(unique(Sub.database$DS)=="Perth")
#       {mtext("Tag ID code",side=2,outer=T,line=-1,font=1,cex=1.25)}
#       mtext("Date",side=1,outer=T,line=5,font=1,cex=1.25)
#       legend("topright",paste("n=",length(Sharks)," sharks",sep=""),bty="n",cex=1.5)
#     })
#     points(tagshark$DateDeployed,as.factor(tagshark$ATAG.NO),pch=17,col="gray42",cex=1.25)
# #    if(length(Atshark)==length(Sharks))points(tagshark$DateDeployed,as.factor(tagshark$ATAG.NO),pch=17,col="gray42",cex=1.25)
# #    if(length(Atshark)<length(Sharks))points(tagshark$DateDeployed,as.factor(as.character(tagshark$ATAG.NO)),pch=17,col="gray42",cex=1.25)
#   }
# 
#   png(file=paste("PresenceTagThruTime.",ActiveFile,".png",sep=""),width=600,height=600)
#   if(Especie%in%c("TK","BW"))
#   {
#     par(mfcol=c(3,1),mar=c(1,4,0.1,1), oma=c(6,1,0.1,1))
#     #par(mfcol=c(3,1),mar=c(6,4,1,1), oma=c(1,1,1,1))
#     plot.detec.thru.time(database.Ningaloo,Useful.Ningaloo)
#     plot.detec.thru.time(database.Perth,Useful.Perth)
#     plot.detec.thru.time(database.SouthWest,Useful.SouthWest)
#   }
#   if(Especie%in%c("GM","WH"))
#   {
#     par(mfcol=c(2,1),mar=c(6,4,1,1), oma=c(1,1,1,1))
#     plot.detec.thru.time(database.Perth,Useful.Perth)
#     plot.detec.thru.time(database.SouthWest,Useful.SouthWest)
#    }
#   dev.off()
# }
# 
# # submit functions to export graphs
# PresenceTagTime(Tagging.Data.SandBar,Tag.Deployed.SandBar,Useful.Vecs.SandBar,"Sandbar")
# PresenceTagTime(Tagging.Data.Dusky,Tag.Deployed.Dusky,Useful.Vecs.Dusky,"Dusky")
# #PresenceTagTime(Tagging.Data.Gummy,Tag.Deployed.Gummy,Useful.Vecs.Gummy,"Gummy")                                                      #UPDATE!!!!!!!
# #PresenceTagTime(Tagging.Data.Whiskery,Tag.Deployed.Whiskery,Useful.Vecs.Whiskery,"Whiskery")
# 
# 
# 
# #--RANDOM WALK SIMULATION--
#   #1. Create list by shark and species, of consecutive hits and for Ningaloo
# list.by.shark=function(database)
# {
#   database=subset(database,DS=="Ningaloo")
# 
#   n.shark=unique(as.numeric(database$Tag));n.shark=n.shark[!is.na(n.shark)]
#   singlehit=table(database$Tag)
# #  singlehit=which(singlehit==1)
#   singlehit=which(singlehit>1)
#   n.shark=n.shark[match(names(singlehit),n.shark)]
# #  n.shark=n.shark[-match(names(singlehit),n.shark)]
#   byShark=list()
#   for (i in 1:length(n.shark))
#   {
#       datos=subset(database,Tag==n.shark[i])
#       datos=datos[order(datos$Date.time),]
# 
#       #calculate speed, distance, bearing and delta t
#       datos$bearing=datos$speed=datos$distance.m=datos$delta.t.sec=NA
#       for(j in 2:nrow(datos))           #MISSING:  REMOVE FOR LOOP, SEE #4. OF SIMULATION EVALUATION!!
#       {
#         if(!(datos$Receiver.Number[j]==datos$Receiver.Number[j-1]) & !(is.na(datos$Recovery_latitude[j]) | is.na(datos$Recovery_latitude[j-1])))
#         {
#           #random sample of receiver location considering detection range     #NOTE: change to triangular SEE #4. OF SIMULATION EVALUATION!
#           datos$Recovery_latitude[j]=rnorm(1,mean=datos$Recovery_latitude[j],sd=detec.range)
#           datos$Recovery_longitude[j]=rnorm(1,mean=datos$Recovery_longitude[j],sd=detec.range)
# 
#           datos$delta.t.sec[j]=as.numeric(datos$Date.time[j]-datos$Date.time[j-1])*60     #delta t
#           Start=c(datos$Recovery_longitude[j-1],datos$Recovery_latitude[j-1])
#           End=c(datos$Recovery_longitude[j],datos$Recovery_latitude[j])
#           datos$distance.m[j]=distCosine(Start,End)                     #distance in metres (using Great Circle distance, see package "geosphere")
#           datos$bearing[j]=bearing(Start,End)                           #bearing in degrees (N = 0 and 360, E = 90, S = 180, and W = 270 degrees)
#           datos$speed[j]=  datos$distance.m[j]/datos$delta.t.sec[j]                          #speed, in m/s
# 
#           if(datos$delta.t.sec[j]>delta.t.threshold)
#           {
#             datos$delta.t.sec[j]=NA
#             datos$distance.m[j]=NA
#             datos$speed[j]=NA
#             datos$bearing[j]=NA
#           }
# 
#         }
#       }
# 
#       byShark[[i]]=datos
#   }
#   return(byShark)
# }
# 
# Tag.SandBar.byshark=list.by.shark(Tagging.Data.SandBar)
# Tag.Dusky.byshark=list.by.shark(Tagging.Data.Dusky)
# #Tag.Gummy.byshark=list.by.shark(Tagging.Data.Gummy)
# #Tag.Whiskery.byshark=list.by.shark(Tagging.Data.Whiskery)
# 
# 
#   #2. Draw random samples of speed and bearing          (SEE SIMS 2010 PAGE 428) TEST IF LOG FREQUNCY VS LOG MOVE STEP SLOPE = 2 (LEVY)
# #2.1. estimate distributions
# Distributions=function(database)
# {
#   newdata=do.call(rbind,database)
#   newdata$distance.m=ifelse(newdata$speed>speed.threshold,NA,newdata$distance.m)
#   newdata$speed=ifelse(newdata$speed>speed.threshold,NA,newdata$speed)
#   distance.m=newdata$distance.m[!is.na(newdata$distance.m)]
#   speed=newdata$speed[!is.na(newdata$speed)]
#   bearing=newdata$bearing[!is.na(newdata$bearing)]
#   
#   par(mfcol=c(3,1),mai=c(.8,.75,.2,.2),omi=c(.2,.5,.2,.5))
#   #displacement
#   hist(distance.m,freq =F,breaks=20,xlab="displacement (m)",main="")
# 
#   #speed
#   hist(speed,freq =F,breaks=20,xlab="speed (m/s)",main="")
#   #lines(density(speed),col=2)
#   fit.speed=fitdistr(speed, "exponential")    #exponential distribution gave best fit
#   lines(density(rexp(100,fit.speed$estimate[1])),col=2)
#   legend("topright",c("exponential"),col=c(2),lty=1,bty='n',cex=1.3)
#   #  delta=0;mygamma= 1      #pars of Levy distribution
#   #  fit = vglm(speed ~ 1, levy(idelta=delta, igamma=mygamma),as.data.frame(speed), trace=TRUE)
#   #  a=rlevy(100, coef(fit)[1], coef(fit)[2])
# 
#   #Bearing
#   hist(bearing,freq =F,breaks=360/2,xlab="bearing (0=N, 90=E, 180=S, 270=W)",main="")
#   #lines(density(bearing),col=2)
#   fit.bearing <- normalmixEM(bearing, lambda = 0.5, mu = c(130, 310), sigma = 10)
#   lines(density(rnormmix(1000, lambda=fit.bearing$lambda, mu=fit.bearing$mu, sigma=fit.bearing$sigma)),col=3)
#   legend("topright",c("mixed normal"),col=c(3),lty=1,bty='n',cex=1.3)
#   
#   return(list(speed.dist=fit.speed,bearing.dist=fit.bearing))
# }
# 
# Speed.dist.Sandbar=Distributions(Tag.SandBar.byshark)$speed.dist
# Bearing.dist.Sandbar=Distributions(Tag.SandBar.byshark)$bearing.dist
# 
# Speed.dist.Dusky=Distributions(Tag.Dusky.byshark)$speed.dist
# Bearing.dist.Dusky=Distributions(Tag.Dusky.byshark)$bearing.dist
# 
# #Speed.dist.Gummy=Distributions(Tag.Gummy.byshark)$speed.dist
# #Bearing.dist.Gummy=Distributions(Tag.Gummy.byshark)$bearing.dist
# 
# #Speed.dist.Whiskery=Distributions(Tag.Whiskery.byshark)$speed.dist
# #Bearing.dist.Whiskery=Distributions(Tag.Whiskery.byshark)$bearing.dist
# 
# 
#    #ACA!!!!
# 
# #2.2. draw samples and create walk
#  Samples=function(n,fit.speed,fit.bearing,Time,position)
#  {
#     #draw random samples         
#   speed.sample=rexp(n,fit.speed$estimate[1])
#   bearing.sample=rnormmix(n, lambda=fit.bearing$lambda, mu=fit.bearing$mu, sigma=fit.bearing$sigma)
#   bearing.sample=bearing.sample[bearing.sample<=360]
# 
#  distance=speed.sample*Time         #(in metres)
# 
#  Step.size=quantile(distance,0.95)/(60*1853)    #max length of step (in degrees)
#  names(Step.size)=NULL
#  
#  position=as.numeric(position)
# 
#  # compute path
# n <- 5
# rw <- matrix(0, ncol = 2, nrow = n)
# # generate the indices to set the deltas
# indx <- cbind(seq(n), sample(c(1, 2), n, TRUE))
# 
# # now set the values
# rw[indx] <- sample(c(-1, 1), n, TRUE)
# # cumsum the columns
# rw[,1] <- cumsum(rw[, 1])
# rw[, 2] <- cumsum(rw[, 2]) 
#   plot(0,type="n",xlab="x",ylab="y",main="Random Walk Simulation In Two Dimensions",col=1:10,xlim=range(rw[,1]),ylim=range(rw[,2]))
#    # use 'segments' to color each path
#   segments(head(rw[, 1], -1)
#          , head(rw[, 2], -1)
#          , tail(rw[, 1], -1)
#          , tail(rw[, 2], -1)
#         ,  col = rainbow(nrow(rw) -1)  # a range of colors
#         )
#   end<-cbind(rw[n,1],rw[n,2])
#   start<-cbind(0,0)
#   points(start,pch=16,col="green", cex = 3)
#   points(end,pch=16,col="red", cex = 3)
# 
# 
#       
#     #Create biased random walk
# 
#  next.point=destPoint(position, b=bearing.sample[1], d=distance[1])
# 
# 
#  plot(simm.brown(date = 1:10, x0 = c(119, -10), h = 1, id = "A1", burst = 1))
#  
#  
#  }
# 
# 
#  N=1000 #number of random samples
#  Seed.location=data.frame(long=Tag.Deployed.SandBar$MID.LONG,lat=Tag.Deployed.SandBar$MID.LAT) #initial position of shark
#  TimeStep=60     #time step of random walk, in seconds
# 
#  Random.Walk.Sandbar=Samples(N,Speed.dist.Sandbar,Bearing.dist.Sandbar,TimeStep,Seed.location[1,])
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #MISSING!!!!!!
# ##10.  Create table and figures of movement summaries                                                             #MISSING, CHANGE ACCORDINGLY!!!!!!
# #
# #table.function=function(tagshark,database,ActiveFile)
# #{
# #  database=subset(database,!(Tag=="       "))
# #  Atshark=sort(unique(database$ID[!is.na(database$ID)]))
# #
# #  Tabla=data.frame(Transmitter.number=tagshark$ID,Total.length.cm=tagshark$Length,Sex=tagshark$Sex,Release.date=tagshark$Date,
# #  Release.location=tagshark$Location)
# #
# #  hits=table(database$ID)
# #  receivers=rowSums(ifelse(table(database$ID,database$"Receiver S/N")>0,1,0))
# #  days.monitored=as.numeric(max(max(database$Date),as.Date("2011-11-01"))-tagshark$Date)
# #  days.detected=rowSums(ifelse(table(database$ID,database$Date)>0,1,0))
# #  Min.cons.det.days=Max.cons.det.days=days.start.end.det=Total.cons.det.days=NULL
# #  for(i in 1:length(Atshark))
# #  {
# #    datos=subset(database,ID==Atshark[i])
# #    datos=datos[!(duplicated(datos$Date)),]
# #    Breaks <- c(0, which(diff(datos$Date) != 1), length(datos$Date))
# #    consdays=sapply(seq(length(Breaks) - 1),function(i) length(datos$Date[(Breaks[i] + 1):Breaks[i+1]]))
# #    Min.cons.det.days=rbind(Min.cons.det.days,min(consdays))
# #    Max.cons.det.days=rbind(Max.cons.det.days,max(consdays))
# #    Total.cons.det.days=c(Total.cons.det.days,consdays)
# #    if(length(consdays)==1) days.start.end.det=rbind(days.start.end.det,NA)
# #    if(length(consdays)>1) days.start.end.det=rbind(days.start.end.det,max(datos$Date)-min(datos$Date))
# #  }
# #
# #  Residency=days.detected/days.start.end.det
# #
# #  Tabla$Number.of.hits=hits
# #  Tabla$Number.of.receivers=receivers
# #  Tabla$Number.of.days.monitored=days.monitored
# #  Tabla$Number.of.days.detected=days.detected
# #  #Tabla$Min.Consecutive.days.present=Min.cons.det.days
# #  #Tabla$Max.Consecutive.days.present=Max.cons.det.days
# #  Tabla$Temp.Residence.time.percent=round(100*Residency,0)
# #
# #  write.table(Tabla,paste("Table_1.",ActiveFile,".csv",sep=""), col.names=T,row.names=F, sep=",")
# #
# #
# #      #frequency of continuous detections
# #  missingdays=match(1:max(Total.cons.det.days),Total.cons.det.days)
# #  names(missingdays)=1:max(Total.cons.det.days)
# #  missingdays=missingdays[is.na(missingdays)]
# #  TableConsDy=table(Total.cons.det.days)
# #  TableConsDy=c(TableConsDy,missingdays)
# #  TableConsDy=TableConsDy[match(1:max(Total.cons.det.days),names(TableConsDy))]
# #  #maxY= round(max(table.ID.by.Date)*1.2)
# #  png(file=paste("frequency.cont.det.",ActiveFile,".png",sep=""),width=600,height=600)
# #  bp=barplot(TableConsDy, space = 0,xlab="Days continuously detected",ylab="Frequency",args.legend = list(bty="n",cex=1.1),
# #  main=paste("Continuous detections (",ActiveFile,")",sep=""),cex.main=1.0,font.main=1,ylim=c(0,max(TableConsDy,na.rm=T)+10),
# #  las=1,col="gray42",cex.names=1.0,cex.lab=1.1,cex.axis=1.1,axis.lty=1)
# #  box()
# #  dev.off()
# #}
# #
# #table.function(Tag.Deployed,WhitePointers,"OTN & SMN")       #all receivers
# #table.function(Tag.Deployed.SMN,WhitePointers.SMN,"SMN")         #SMN only
# #table.function(Tag.Deployed.VR4s,WhitePointers.VR4s,"VR4s")       #VR4s only
# #
# #
# #

# fn.timeline.hits=function(Detections,ESPECIES)
# { 
#   DATA=Detections
#   #create table
#   Table.hits.date1=with(DATA,table(TagCode,as.character(Date.local)))
#   Table.hits.date=with(DATA,table(TagCode,Julian))
#   First.Jul=colnames(Table.hits.date)[1] #this is first julian date
#   Hits.date=colnames(Table.hits.date1)
#   
#   Table.hits=table(DATA$TagCode)
#   Range.hits=range(Table.hits)
#   
#   # add sex, size and release area
#   these.ones=match(c("TagCode","Sex","Area.release","ReleaseDate"),names(DATA))
#   noduplicates <- DATA[!duplicated(DATA$TagCode),these.ones]
#   noduplicates$Sex=ifelse(as.character(noduplicates$Sex)=="","U",as.character(noduplicates$Sex))
#   noduplicates=noduplicates[order(noduplicates$TagCode),]
#   Table.hits.date=cbind(noduplicates,as.data.frame.matrix(Table.hits.date))
#   
#   if(!ESPECIES=="Thickskin") Rel.date.Julian=as.numeric(round((Table.hits.date$ReleaseDate-First.release)))
#   if(ESPECIES=="Thickskin") Rel.date.Julian=as.numeric(round((Table.hits.date$ReleaseDate-First.release)/(3600*24)))
#   This.col=match(First.Jul,colnames(Table.hits.date)) 
#   
#   # create colors by area
#   Table.tag.day.area=with(DATA,table(TagCode,Julian,as.character(Area)))
#   DIM.area=dim(Table.tag.day.area)[3]
#   
#   Replace.number=dimnames(Table.tag.day.area)[[3]]
#   Replace.number=ifelse(Replace.number=="North.WA",1,ifelse(Replace.number=="South.WA",5,
#           ifelse(Replace.number=="SA",10,NA)))        #area conversion (1:NW.WA, 5:SW.WA, 10:SA)
#   
#   #convert observations to levels of areas 
#   for(j in 1:DIM.area) Table.tag.day.area[,,j]=ifelse(Table.tag.day.area[,,j]>0,Replace.number[j],0)
#   
#   if(DIM.area==1) Table.tag.day.area=Table.tag.day.area[,,1]
#   if(DIM.area==2)Table.tag.day.area=Table.tag.day.area[,,1]+Table.tag.day.area[,,2]
#   if(DIM.area==3)Table.tag.day.area=Table.tag.day.area[,,1]+Table.tag.day.area[,,2]+Table.tag.day.area[,,3]
#   
#   Mat.colors=as.data.frame.matrix(Table.tag.day.area)
#   if(sum(rownames(Mat.colors)==Table.hits.date$TagCode)==nrow(Mat.colors))
#   {
#     Mat.colors$Area.release=as.character(Table.hits.date$Area.release)
#     Mat.colors$Area.release=ifelse(Mat.colors$Area.release=="North.WA",2,
#                                    ifelse(Mat.colors$Area.release=="South.WA",8,NA))
#   }
#   Add.This=Mat.colors$Area.release  
#   Mat.colors=Mat.colors[,-match("Area.release",colnames(Mat.colors))]
#   Mat.colors=Add.This*Mat.colors
#   
#   Mat.colors=ifelse(Mat.colors==2,"black",
#                     ifelse(Mat.colors==20,"brown3",
#                            ifelse(Mat.colors==10,"forestgreen",
#                                   ifelse(Mat.colors==8,"blue",
#                                          ifelse(Mat.colors==80,"salmon",
#                                                 ifelse(Mat.colors==40,"red",NA))))))
#   
#   Mat.colors.legends=c("Rel NW.WA-Det NW.WA","Rel NW.WA-Det SA","Rel NW.WA-Det SW.WA",
#                        "Rel SW.WA-Det NW.WA","Rel SW.WA-Det SA","Rel SW.WA-Det SW.WA")
#   Mat.colors.legends.colors=c("black","brown3","forestgreen","blue","salmon","red")
#   
#   #Remove as more areas get hits
#   if(ESPECIES=="Dusky")
#   {
#     Mat.colors.legends=Mat.colors.legends[c(1,3,4)]
#     Mat.colors.legends.colors=Mat.colors.legends.colors[c(1,3,4)]
#   }
#   if(ESPECIES=="Thickskin")
#   {
#     Mat.colors.legends=Mat.colors.legends[c(1,2)]
#     Mat.colors.legends.colors=Mat.colors.legends.colors[c(1,2)]
#   }
#   if(ESPECIES%in%c("Gummy","Whiskery"))
#   {
#     Mat.colors.legends=Mat.colors.legends[4:6]
#     Mat.colors.legends.colors=Mat.colors.legends.colors[4:6]
#   }  
#   
#   # plot proportions
#   tiff(file=paste("Outputs_movement/Figure3.",ESPECIES,".prop.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
#   par(mai=c(.8,.8,.01,.01),las=1,mgp=c(2, 1, 0))
#   
#   bubble.plot.detections(as.numeric(colnames(Table.hits.date[,This.col:ncol(Table.hits.date)])),
#                          1:nrow(Table.hits.date),
#                          Table.hits.date[,This.col:ncol(Table.hits.date)],3,"Date","Tag",Rel.date.Julian,
#                          Mat.colors,Mat.colors.legends,Mat.colors.legends.colors,"proportion","bottomleft")
#   axis(side = 1, at = Months.for.plot,labels = Months.for.plot.label, tcl = -0.3,cex.axis=.75,padj=-2) # minor ticks
#   axis(side = 1, at = Start.of.Year, labels = Start.of.Year.label, tcl = -0.75) # major ticks & labels
#   
#   axis(side = 2, at = 1:nrow(Table.hits.date), labels = Table.hits.date$TagCode, tcl = -0.5,cex.axis=.5,hadj=0.75) # major ticks & labels
#   dev.off()
#   
#   
#   # plot presence absence
#   tiff(file=paste("Outputs_movement/Figure3.",ESPECIES,".pres_abs.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
#   par(mai=c(.8,.8,.01,.01),las=1,mgp=c(2, 1, 0))
#   
#   bubble.plot.detections(as.numeric(colnames(Table.hits.date[,This.col:ncol(Table.hits.date)])),
#                          1:nrow(Table.hits.date),
#                          Table.hits.date[,This.col:ncol(Table.hits.date)],3,"Date","Tag",Rel.date.Julian,
#                          Mat.colors,Mat.colors.legends,Mat.colors.legends.colors,"presence","bottomleft")
#   axis(side = 1, at = Months.for.plot,labels = Months.for.plot.label, tcl = -0.3,cex.axis=.75,padj=-2) # minor ticks
#   axis(side = 1, at = Start.of.Year, labels = Start.of.Year.label, tcl = -0.75) # major ticks & labels
#   axis(side = 2, at = 1:nrow(Table.hits.date), labels = Table.hits.date$TagCode, tcl = -0.5,cex.axis=.5,hadj=0.75) # major ticks & labels
#   dev.off()
# }
# 
# for (i in 1:N.sp) fn.timeline.hits(Detections.species[[i]],SPECIES[i])


######################################################################
#DELETE ALL THIS ONCE SM NETWORK DATABASE IS RUNNING

#   #1. load data
#     #1.1. listening stations                                                              #UPDATE
# #note: these files need to be updated as receiver deployment progresses
# 
# path="M:/Fisheries Research/FinFish/Shark/Braccini/Ningaloo October 2011 Data/"
# 
#       #Ningaloo
# channel <- odbcConnectExcel(paste(path,"Field Sheets/Field Sheets October 2011",sep=""))
# Ningaloo<- sqlFetch(channel,"Recoveries_Field_Sheets")
# close(channel)
# if (is.na(Ningaloo$"Receiver (S/N)"[1])) Ningaloo=Ningaloo[-1,]
# Ningaloo=Ningaloo[!(is.na(Ningaloo$"Receiver (S/N)")),]
# Ningaloo$"Recovery latitude"=-Ningaloo$"Recovery latitude"
# colnames(Ningaloo)[match(c("Recovery latitude","Recovery longitude"),names(Ningaloo))]=c("Recovery_latitude" , "Recovery_longitude")
# 
#       #SouthWest                                                                          #MISSING: UPDATE WITH REAL LOCATIONS 
#                                                                                                 # (see "Analysis conv tagging data.r")
# SouthWest<- read.csv("H:/Matias WA Fisheries/Data/Mapping/AcousticLinesSWA.csv")
# 
#       #Perth
# channel1 <- odbcConnectExcel("M:/Fisheries Research/FinFish/Shark/Braccini/Shark MonitoringOTN/SMN_OTN Running Sheet 2011 copy")
# SMN_VR2_2010<- sqlFetch(channel1,"SMN VR2W 2010")
# SMN_VR4G_2009<- sqlFetch(channel1,"VR4G Deployed 2009")
# SMN_VR4G_2009$Long[2]=115.73338       #fill in R import f... up
# SMN_VR4G_2009$Long[17]=115.47654
# OTN_2009a<- sqlFetch(channel1,"OTN 2009a")
# close(channel1)
# 
# 
#     #1.2. shark detections                                                              #UPDATE
# #note: these files need to be updated as tagging progresses
# 
# all.csv=list.files(path,pattern="csv")    #get all .csv files
# Tagging.Data.Ningaloo=NULL
# for(i in 1:length(all.csv))
# {
#   Tagging.Data.Ningaloo<- rbind(Tagging.Data.Ningaloo,read.csv(paste(path,all.csv[i],sep="")))
# }
# 
# path.Perth="M:/Fisheries Research/FinFish/Shark/Braccini/PerthDummy/"
# all.csv.Perth=list.files(path.Perth,pattern="csv")    #get all .csv files
# Tagging.Data.Perth=NULL
# for(i in 1:length(all.csv.Perth))
# {
#   Tagging.Data.Perth<- rbind(Tagging.Data.Perth,read.csv(paste(path.Perth,all.csv.Perth[i],sep="")))
# }
# 
# path.SW="M:/Fisheries Research/FinFish/Shark/Braccini/SouthWestDummy/"
# all.csv.SW=list.files(path.SW,pattern="csv")    #get all .csv files
# Tagging.Data.SouthWest=NULL
# for(i in 1:length(all.csv.SW))
# {
#   Tagging.Data.SouthWest<- rbind(Tagging.Data.SouthWest,read.csv(paste(path.SW,all.csv.SW[i],sep="")))
# }
# 
# 
#     #1.3. tag deployment data                                                              #UPDATE
# #note: these files need to be updated as tagging progresses
# Tag.Deployed=read.csv("H:/Matias WA Fisheries/Data/Tags_ID/AcousticsDeployed WA JunAug2011.csv")
# Tag.Deployed$DateDeployed=as.Date(Tag.Deployed$DateDeployed,"%d/%m/%Y")      #convert to date
# Tag.Deployed$MID.LAT=-Tag.Deployed$MID.LAT

#------ PROCEDURE SECTION ------

#----1. Manipulate data and create data frames----

#1.1. Create Detections file
#note: detections data frame with detection data from SMN and AATAMS

# #tidy up datasets
# these.AATAMS.vars=c("ID","Species2","Sex2","receiver.ID","Release.Date","ReleaseLatitude2",
#                     "ReleaseLongitude2","latitude","longitude","DateTime.local","Date.local","Time.local")
# matched.AATAMS.vars=match(these.AATAMS.vars,names(AATAMS))
# AATAMS=AATAMS[,matched.AATAMS.vars]
# 
# not.this.SMN=match(c("X"),names(SMN))
# not.this.SMN=not.this.SMN[!is.na(not.this.SMN)]
# if(length(not.this.SMN)>0)SMN=SMN[,-not.this.SMN]
# 
# 
# SMN$Project="SMN"
# AATAMS$Project="AATAMS"
# AATAMS$Depth=NA
# 
# re.arranged.cols=match(c("ID","Species2","Sex2","receiver.ID","Release.Date","ReleaseLatitude2","ReleaseLongitude2",
#                          "latitude","longitude","Depth","DateTime.local","Date.local","Time.local","Project"),
#                        names(AATAMS))
# AATAMS=AATAMS[,re.arranged.cols]
# 
# colnames(AATAMS)=colnames(SMN)
# 
# if(class(SMN$ReleaseDate)=="factor") SMN$ReleaseDate=as.POSIXlt(as.character(SMN$ReleaseDate),format="%d-%b-%y")
# 
# if(class(AATAMS$ReleaseDate)=="factor") AATAMS$ReleaseDate=as.POSIXlt(as.character(AATAMS$ReleaseDate))
# 
# AATAMS$SerialNumber=as.numeric(substring(as.character(AATAMS$SerialNumber), first=6, last = 30))
# 
# #merge data files
# Detections=rbind(SMN,AATAMS)
# rm(list=c("AATAMS","SMN"))
################################################################################################



# fun.plot.mig=function(DATA,tags,LONG,LAT)
# {
#   plotMap(worldLLhigh, xlim=LONG,ylim=LAT,plt = c(.001, 1, 0.075, 1),
#           col=COLOR,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
#   axis(side = 1, at =seq(LONG[1],LONG[2]), labels = seq(LONG[1],LONG[2]), tcl = .5,las=1,cex.axis=1.5,padj=-.5)
#   axis(side = 2, at = seq(round(LAT[1]),LAT[2]-1), labels = -(seq(round(LAT[1]),LAT[2]-1)),tcl = .5,las=2,cex.axis=1.5,hadj=.75)
#   box(lwd=2)
#   mtext("Latitude (ºS)",side=2,line=2.5,las=3,cex=1.5)
#   mtext("Longitude (ºE)",side=1,line=1.55,cex=1.5)
#   contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=LONG,xlim=LAT, zlim=c(-1,-200),
#           nlevels = 4,labcex=1.,lty = 1,col=c(COLOR,COLOR,COLOR,COLOR,"transparent"),add=T)
#   
#   for(t in 1:length(tags)) 
#   {
#     Data=subset(DATA,TagCode==tags[t])
#     
#     Col=col.mig[t]
#   
#     if(unique(Data$ReleaseLongitude<115.15 & Data$ReleaseLatitude>(-25)))
#     {
#       Data$Longitude=jitter(Data$Longitude,0.2)
#       Data$Latitude=jitter(Data$Latitude,.2)
#       segments(Data$ReleaseLongitude,Data$ReleaseLatitude,Shark.bay[1],Shark.bay[2],col=Col,lwd=2)
#       arrows(Shark.bay[1],Shark.bay[2],Data$Longitude,Data$Latitude,col=Col,lwd=2,length=.15)
#     }
#     if(Data$ReleaseLongitude<115.15 & Data$ReleaseLatitude<(-25))
#     {
#       arrows(Data$ReleaseLongitude,Data$ReleaseLatitude,Data$Longitude,Data$Latitude,col=Col,lwd=2,length=.15)
#     }
#     if(Data$ReleaseLongitude>=115.15 & Data$ReleaseLongitude<116.425)
#     {
#       segments(Data$ReleaseLongitude,Data$ReleaseLatitude,Cape.Leuwin[1],Cape.Leuwin[2],col=Col,lwd=2)
#       segments(Cape.Leuwin[1],Cape.Leuwin[2],Shark.bay[1],Shark.bay[2],col=Col,lwd=2)
#       arrows(Shark.bay[1],Shark.bay[2],Data$Longitude,Data$Latitude,col=Col,lwd=2,length=.15)
#     }  
#     if(Data$ReleaseLongitude>=116.425)
#     {
#       segments(Data$ReleaseLongitude,Data$ReleaseLatitude,Mid.point[1],Mid.point[2],col=Col,lwd=2)
#       segments(Mid.point[1],Mid.point[2],Cape.Leuwin[1],Cape.Leuwin[2],col=Col,lwd=2)
#       segments(Cape.Leuwin[1],Cape.Leuwin[2],Shark.bay[1],Shark.bay[2],col=Col,lwd=2)
#       arrows(Shark.bay[1],Shark.bay[2],Data$Longitude,Data$Latitude,col=Col,lwd=2,length=.15)
#     }
#     points(Data$ReleaseLongitude,Data$ReleaseLatitude,col=Col,pch=19,cex=1.35)
#   }
#   
# }

#19.2.4 Dusky shark natal migration Old GLM approach 
#note:incomplete
Do.glm="NO"
if(Do.glm=="YES")
{
  
  dat.glm$N=with(dat.glm,ifelse(Zn.rec=="North",1,0))
  
  daylength1 <- function(t,lon, lat)
  {
    t <- as.numeric(t)
    alt <- function(t) sunAngle(t, longitude=loN, latitude=laT)$altitude
    rise=set=rep(NA,length(t))
    for(i in 1:length(t))
    {
      loN=lon[i]
      laT=lat[i]
      rise[i]=uniroot(alt, lower=t[i]-86400/2, upper=t[i])$root
      set[i] <- uniroot(alt, lower=t[i], upper=t[i]+86400/2)$root
    }
    daylen=set - rise
    return(daylen/3600)
  }
  dat.glm$day.length=daylength1(dat.glm$Date.local,dat.glm$Longitude,dat.glm$Latitude) 
  
  
  #add Temp
  Mean.temp1.N$Zn.rec="North"
  Mean.temp1.S$Zn.rec="South"
  Mean.TemP=rbind(Mean.temp1.N,Mean.temp1.S)
  dat.glm=merge(dat.glm,Mean.TemP,by=c("Zn.rec","Month"),all.x=T)
  
  #Temp anomaly
  Mean.T.S=median(subset(dat.glm,Zn.rec=="South",select=c(Temperature))[,1])
  Mean.T.N=median(subset(dat.glm,Zn.rec=="North",select=c(Temperature))[,1])
  dat.glm$Temp.anomaly=with(dat.glm,ifelse(Zn.rec=="North",Temperature-Mean.T.N,Temperature-Mean.T.S))
  
  par(mfcol=c(2,1))
  with(dat.glm,plot(yday.rec,day.length))
  with(dat.glm,boxplot(Temperature~Zn.rec))
  
  par(mfcol=c(3,1))
  with(dat.glm,plot(day.length,N))
  with(dat.glm,plot(Temperature,N))
  with(dat.glm,plot(Temp.anomaly,N))
  with(dat.glm,plot(Month,Temp.anomaly))
  
  pairs(dat.glm[,match(c("Temperature","day.length"),names(dat.glm))])
  pairs(dat.glm[,match(c("Temp.anomaly","day.length"),names(dat.glm))])
  
  dat.glm$Sex=factor(as.character(dat.glm$Sex),levels=c("F","M"))
  dat.glm$TagCode=as.factor(dat.glm$TagCode)
  
  MOD=glm(N~FL*Sex+day.length+Temp.anomaly, data=dat.glm, family="binomial", maxit=500)
  MOD.random=glmer(N~FL*Sex+day.length+Temp.anomaly +(1 | TagCode), data=dat.glm, family="binomial")
  
  #term significance
  Anova=anova(MOD, test = "Chisq")
  Anova.random=anova(MOD.random, test = "Chisq")
  
  
  dat.glm$pred=predict(MOD.random,type='response')
  plot(dat.glm$N)
  points(dat.glm$pred,col=2)
  
  #predict new data
  fn.do.new.dat=function(varTerm,varTermValues,FACTRS,SEQ_len)
  {
    Rang=range(varTermValues)
    Var=seq(Rang[1],Rang[2],length.out=SEQ_len)
    
    other.vars=c("FL","Sex","day.length","Temp.anomaly","TagCode")
    other.vars=other.vars[-match(varTerm,other.vars)]
    fixedTerms=dat.glm[,match(other.vars,names(dat.glm))]
    
    #mean of continuous
    cont.vrs=other.vars[-match(FACTRS,other.vars)]
    This.vars=rep(NA,length(cont.vrs))
    for(s in 1:length(cont.vrs))This.vars[s]=mean(subset(fixedTerms,select=cont.vrs[s])[,1])
    
    #most common of factors
    Most.common=function(what)
    {
      ID=match(what,names(fixedTerms))
      Tabla=sort(table(fixedTerms[,ID]))
      return(names(Tabla[length(Tabla)]))
    }
    This.Fctrs=rep(NA,length(FACTRS))
    for(s in 1:length(FACTRS))This.Fctrs[s]=Most.common(FACTRS[s])
    
    Fixed.vars=data.frame(t(This.vars))
    Fixed.fctrs=data.frame(t(This.Fctrs))
    Fixed=cbind(Fixed.vars,Fixed.fctrs)
    names(Fixed)=c(cont.vrs,FACTRS)
    
    ID=match(FACTRS,names(Fixed))
    for(q in 1:length(ID))
    {
      Fixed[,ID[q]]=factor(Fixed[,ID[q]],levels=levels(dat.glm[,match(FACTRS[q],names(dat.glm))]))
    }
    
    D=cbind(Var,Fixed)
    names(D)[1]=varTerm
    return(D)
  }
  
  NEWDATA_Temp.anom=fn.do.new.dat("Temp.anomaly",dat.glm$Temp.anomaly,FACTRS=c("Sex","TagCode"),SEQ_len=100)
  NEWDATA_FL=fn.do.new.dat("FL",dat.glm$FL,FACTRS=c("Sex","TagCode"),SEQ_len=100)
  NEWDATA_day.length=fn.do.new.dat("day.length",dat.glm$day.length,FACTRS=c("Sex","TagCode"),SEQ_len=100)
  
  Pred_Temp.anom=predict(MOD.random,newdata=NEWDATA_Temp.anom,type='response')
  Pred_FL=predict(MOD.random,newdata=NEWDATA_FL,type='response')
  Pred_day.length=predict(MOD.random,newdata=NEWDATA_day.length,type='response')
  
  plot(NEWDATA_FL$FL,Pred_FL,type='l')
  plot(NEWDATA_day.length$day.length,Pred_day.length,type='l')
  plot(NEWDATA_Temp.anom$Temp.anomaly+(Mean.T.N+Mean.T.S)/2,Pred_Temp.anom,type='l')
  
  #sinusoidal model NOT REALLY, just see plot of N and the predictors!!!!
  
  
  MOD.sinusoidal <- glm(N ~ sin(2*pi*day.length)+
                          sin(2*pi*Temp.anomaly)+Sex, data=dat.glm, family="binomial")
  #   MOD.sinusoidal <- glm(N ~ sin(2*pi*day.length)+cos(2*pi*day.length)+
  #                           sin(2*pi*Temp.anomaly)+cos(2*pi*Temp.anomaly)+Sex, data=dat.glm, family="binomial")
  
  Pred.sinu_Temp.anom=predict(MOD.sinusoidal,newdata=NEWDATA_Temp.anom,type='response')
  Pred.sinu_day.length=predict(MOD.sinusoidal,newdata=NEWDATA_day.length,type='response')
  
  plot(NEWDATA_Temp.anom$Temp.anomaly+(Mean.T.N+Mean.T.S)/2,Pred.sinu_Temp.anom,type='l')
  plot(NEWDATA_day.length$day.length,Pred.sinu_day.length,type='l')
}
