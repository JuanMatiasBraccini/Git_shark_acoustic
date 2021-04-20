#Manually get AATMAS depth
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')


Bathymetry_120=read.table(handl_OneDrive("Data/Mapping/get_data112_120.cgi"))
Bathymetry_138=read.table(handl_OneDrive("Data/Mapping/get_data120.05_138.cgi"))
Bathymetry=rbind(Bathymetry_120,Bathymetry_138)
Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
xbat=sort(unique(Bathymetry$V1))
ybat=sort(unique(Bathymetry$V2))

#reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))



setwd(handl_OneDrive("Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/AATAMS_receiver_location"))

AATAMS=read.csv("AATAMS_receivers.manipulated_25_09_2012.csv")

a=subset(Detections,Project=="AATAMS")   #note: "Detections" is loaded in in "Analyses_FRDC.R"

CLOS=paste("gray",round(seq(10,90,length.out=15)),sep="")
#CLOS=c("yellow","green","coral","chartreuse4","chocolate4","red","azure4","orange","blue","black")
fn.check=function(Ylim,Xlim,MAIN)
{
  plot(AATAMS$longitude,AATAMS$latitude,ylim=Ylim,xlim=Xlim,pch=19,col=2,cex=1.5,main=MAIN,ylab="Latitude",xlab="longitude")
  contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=Ylim,xlim=Xlim, zlim=c(-1,-300),
          #nlevels = length(CLOS),labcex=1,lty = 1,col=CLOS,add=T)
  levels=seq(-150,-10,10),col=CLOS,add=T)
}

fn.dec=function(LAT,out)
{
  b=subset(a,Latitude<=(LAT[1]) & Latitude>=(LAT[2]))
  b$dupli=paste(b$Longitude,b$Latitude)
  b=b[!duplicated(b$dupli),]
  points(b$Longitude,b$Latitude,col="chartreuse4",cex=3)
  text(b$Longitude,b$Latitude*out,b$Depth,cex=1,col="chartreuse4")
  print(range(b$Longitude))
}

#line 1
fn.check(Ylim=c(-21.905,-21.85),Xlim=c(113.89, 113.94),"Line 1")
fn.dec(c(-21.85,-21.905),0.9997)

#grid line 1
fn.check(Ylim=c(-22.01,-21.94),Xlim=c(113.89, 113.95),"grid line 1")
fn.dec(c(-21.94,-22.01),0.9997)

#line 2
fn.check(Ylim=c(-22.615, -22.575),Xlim=c(113.58, 113.67),"Line 2")
fn.dec(c(-22.575,-22.615),0.99975)


#line 3
fn.check(Ylim=c(-23.13, -23.10),Xlim=c(113.60, 113.765),"Line 3")
fn.dec(c(-23.10,-23.13),0.99975)


#line 3 Bays 
fn.check(Ylim=c(-23.1, -22.95),Xlim=c(113.76, 113.81),"Line 3 bay")
fn.dec(c(-22.95,-23.1),0.99975)




#Get AATAMS receiver serial number by line to get depth from website
fn.dec=function(LAT)
{
  b=subset(a,Latitude<=(LAT[1]) & Latitude>=(LAT[2]))
  b$dupli=paste(b$Longitude,b$Latitude)
  b=b[!duplicated(b$dupli),]
  b=b[order(b$SerialNumber,b$Date.local),]
  print(b[,match(c("SerialNumber","Date.local"),names(b))])
  
  plot(b$Longitude,b$Latitude,col="chartreuse4",cex=3)
  text(b$Longitude,b$Latitude,b$SerialNumber,cex=1,col=2)
  
}
fn.dec(c(-21.94,-22.01))
