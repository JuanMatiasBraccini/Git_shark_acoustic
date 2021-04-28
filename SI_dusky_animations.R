rm(list=ls(all=TRUE))   #clear log

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

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

setwd(handl_OneDrive("Analyses\\Acoustic_tagging\\Dusky_migration"))

#---INPUT PARAMETERS----

n.trajectories=10 # number of simulated trajectories (i.e. number of tagged sharks)  
Months=12
Days=30
Samp.per.age=Months*Days  
Max.Age=3
Apply.mortality="YES"
if(Apply.mortality=="NO")Z=seq(0,0,length.out=Max.Age)else  #total mortality at age
Z=seq(.3,.01,length.out=Max.Age)
Age.move=2
n.samples=Max.Age*Samp.per.age
min.speed=0.1 # (in m/s)
max.speed.CRW=0.7 #(~60 km per day) 
Inc.speed=3.5  #increase in cruising speed when migration
speed.CRW=c(min.speed,max.speed.CRW)

#how to display trajectories
Show.what="arrows"
#Show.what="points"
Time=60*60*24  #calculate distance moved over 1 day (in seconds)
show.n=10   #show 10 points at a time to reduce number of frames

#Growth pars
Linf=374.4
K=.0367
Lo=75
TL=function(age) Lo+(Linf-Lo)*(1-exp(-K*age))

# Boundary that defines possible movement (avoid land)
boundary=matrix(c(c(
  114.0673, 113.5158, 113.6430, 113.2188, 113.6430, 112.8369, 113.6430, 114.0249, 114.8310, 115.5946,
  115.6795, 114.9158,
  115.0007, 115.9341 ,117.3766, 118.3099, 118.3099, 117.2493, 115.8916, 115.1704, 114.6613, 114.4491, 
  114.8734, 115.1280,
  114.8734, 114.2370, 113.4733, 112.5, 112.5, 112.7279,
  114.0673, 115.0007, 116.4007, 116.8250, 116.1886, 115.5522, 115.1280, 114.6613, 114.2794, 114.0673
),
  c(
    -21.63177, -22.45589, -23.42986, -24.25398, -25.37779, -25.41525, -26.72635, -27.88762,
    -29.61079, -31.78348, -33.24443,
    -33.43173, -34.25585, -34.78030, -35.07998, -34.89268, -35.49204, -35.49204, -35.15490,
    -34.70538, -34.36823, -33.43173,
    -33.16951, -32.42030, -31.33396, -29.12381, -27.66286, -25.34033, -23.39240, -21.44447,
    -20.84510, -20.17082, -19.98352, -20.5079, -20.84510 ,-21.25716, -21.63177, -21.78161, -21.63177, -21.63177
  )),ncol=2)
colnames(boundary)=c("lon",'lat')


#Seeding location of each trajectory within boundary
SW.corner=boundary[boundary[,1]<117 &boundary[,2]<(-34),]
a=SpatialPolygons(list(Polygons(list(Polygon(SW.corner)), "x")))  
Seed.location=as.data.frame(spsample(a, n = n.trajectories, "random"))
names(Seed.location)=c("Long","Lat")



#----PROCEDURE---- 

#4. Create trajectories
fn.rndmx=function(a,b) rmixedvm(1,a*pi/180,b*pi/180,1, 1,0.5)*180/pi
fn.berng=function(x)with(x,
         ifelse(year>=Age.move & month%in%c(5,6) 
              & next.point$lat[s-1]<=(-26)& next.point$lat[s-1]>(-34),
              fn.rndmx(310,360),
          ifelse(year>=Age.move & month%in%c(5,6) &
               next.point$lat[s-1]<=(-34),
               fn.rndmx(260,280),
        ifelse(year>=Age.move & month%in%c(11,12) 
               & next.point$lat[s-1]>(-25),
               fn.rndmx(210,230),
         ifelse(year>=Age.move & month%in%c(11,12) 
                & next.point$lat[s-1]>(-32)& next.point$lat[s-1]<=(-25),
                fn.rndmx(135,180),bearing)))))


fn.spid=function(x)with(x,
         ifelse(year>=Age.move & month%in%c(5,6) 
                & next.point$lat[s-1]<=(-24)& next.point$lat[s-1]>(-34),
                runif(1,speed.CRW[1]*Inc.speed,speed.CRW[2]*Inc.speed),
         ifelse(year>=Age.move & month%in%c(5,6) &
                next.point$lat[s-1]<=(-34),
                runif(1,speed.CRW[1]*Inc.speed,speed.CRW[2]*Inc.speed),
          ifelse(year>=Age.move & month%in%c(11,12) 
                & next.point$lat[s-1]>(-25),
                runif(1,speed.CRW[1]*Inc.speed,speed.CRW[2]*Inc.speed),
          ifelse(year>=Age.move & month%in%c(11,12) 
                & next.point$lat[s-1]>(-32)& next.point$lat[s-1]<=(-25),
                runif(1,speed.CRW[1]*Inc.speed,speed.CRW[2]*Inc.speed),speed)))))

#Create trajectory based on correlated random walk
Population=vector('list',n.trajectories)
Yrs=1:Max.Age
yrs=length(Yrs)
  #select surviving trajectories
N=rep(NA,yrs)
N[1]=n.trajectories
for(y in 2:yrs) N[y]=round(N[y-1]*exp(-Z[y-1]))
id=vector('list',length(N)-1)
id[[1]]=1:n.trajectories
for(ii in 2:length(N)) id[[ii]]=sample(id[[ii-1]],N[ii])
nCol=2
  #Loop over each time step
system.time(for(P in 1:length(Population))
{
  MAT=as.data.frame(matrix(rep(NA,n.samples*nCol),ncol=nCol))
  colnames(MAT)=c("year","month")
  MAT$year=rep(1:Max.Age,each=Samp.per.age)
  MAT$month=rep(rep(1:Months,each=Days),Max.Age)
  Age.indx=as.list(seq(Samp.per.age,Max.Age*Samp.per.age,Samp.per.age))
  idx=which(unlist(id)==P)
  
  next.point=matrix(rep(NA,2*nrow(MAT)),ncol=2)
  next.point=as.data.frame(next.point)
  colnames(next.point)=c("lon","lat")
  next.point[1,]=as.numeric(Seed.location[P,])
  
  #loop over time step until death
  for(s in 2:Age.indx[[length(idx)]])
  {
    #Control speed and bearing depending on whether migrating or cruising
    speed=runif(1,speed.CRW[1],speed.CRW[2])
    speed=fn.spid(MAT[s,])
    distance=min(70000,speed*Time)
    
    bearing=rvm(1,0,k=0.010)*180/pi
    bearing=fn.berng(MAT[s,])
    next.point[s,]=destPoint(next.point[s-1,], b=bearing, d=distance)
    
    #recalculate if position outside boundaries 
    Boundary=boundary
    if(MAT[s,]$year<Age.move) Boundary=boundary[boundary[,2]<(-30),]
    if(MAT[s,]$year==Age.move & MAT[s,]$month<6) Boundary=boundary[boundary[,2]<(-30),]
    print(s)
    test=inout(matrix(next.point[s,],ncol=2),Boundary,bound=T)
    if(test==F)
    { repeat
    {
      speed=runif(1,speed.CRW[1],speed.CRW[2])
      speed=fn.spid(MAT[s,])
      distance=min(70000,speed*Time)
      bearing=rvm(1,0,k=0.010)*180/pi
      bearing=fn.berng(MAT[s,])
      next.point[s,]=destPoint(next.point[s-1,], b=bearing, d=distance)
      test=inout(matrix(next.point[s,],ncol=2),Boundary,bound=T)
      if(test==T)break
    }
    }
  }
  next.point$year=MAT$year
  next.point$month=MAT$month
  next.point$ID=P
  next.point$Age=rep(1:Max.Age,each=Samp.per.age)
  Population[[P]]=next.point
})
CRW.trajectory=Population



#5. Animate trajectories
p <- ggmap(get_map(c(116,-30),maptype="satellite", zoom = 5))
p<-p+labs(x="Longitude",y="Latitude")+
  scale_x_continuous(limits = c(min(boundary[,1]*.995), max(boundary[,1])*1.02), expand = c(0, 0))+
  scale_y_continuous(limits = c(min(boundary[,2]*1.02), max(boundary[,2])*.995), expand = c(0, 0))+
  theme(axis.text=element_text(size=12),
      axis.title=element_text(size=14,face="bold"))

#Add circle (life span)
circ.cl="white"
fn.circ=function(pts)
{
  circulo=data.frame(sin = sin(pts), cos = cos(pts))
  circulo$x=mean(117:119)+circulo$sin
  circulo$y=mean(-26:-22)+circulo$cos
  return(circulo)
}
p=p+geom_path(data = fn.circ(pts=seq(0, 2*pi, l = 100)), aes(x, y),size=5,col="grey",alpha=0.75)
p=p+annotate("text", x = mean(117:119), y=-22.25, label = "Life span",col=circ.cl,size=6)


colfunc <- colorRampPalette(c("red", "yellow"))

COL.scale=data.frame(ID=1:length(CRW.trajectory),CLR=colfunc(length(CRW.trajectory)))
COL.scale$CLR=as.character(COL.scale$CLR)
D=do.call(rbind,CRW.trajectory)
D=merge(D,COL.scale,by="ID")
time=sort(unique(D$Age)) 

D$Class=with(D,
    ifelse(Age<Age.move,"Small individuals",
    ifelse(Age>=Age.move,"Large individuals",NA)))
D$Size=TL(D$Age+D$month/13)
D$Size=1.2*D$Size/max(D$Size)

system.time({saveGIF({
  suppressWarnings(for(x in time)
  {
    DF <- subset(D,Age==x)
    N=unique(DF$ID)
    indx=match(N,DF$ID)
    indx=indx[!duplicated(indx)]-1
    Agevec=rep(indx,each=show.n)
    Inc=seq(1,Samp.per.age,by=show.n)
    for(n in 1:length(Inc))
    {
      ii=Agevec+(Inc[n]:(show.n*n))
      df=DF[ii,]
      if(Show.what=="points")df$alpha=rep(seq(1,0.5,length.out=show.n),n.trajectories)
      df=subset(df,!is.na(lat))
      TTL=paste("Population size=",length(unique(df$ID)))
      stop=Inc[n]+(x-1)*Samp.per.age
      if(Show.what=="points")
      {
        print(p+ggtitle(TTL)+
                geom_point(data = df,mapping=aes(x = lon, y = lat,size=df$Size,alpha= df$alpha),color=df$CLR,show.legend=F)+
            geom_path(data = fn.circ(pts=seq(0, stop*2*pi/n.samples, l = 100)), aes(x, y),size=5,col=circ.cl))
      }
      if(Show.what=="arrows")
      {
        print(p+ggtitle(TTL)+
          geom_path(data = df,mapping=aes(x = lon, y = lat,color=df$CLR),size=.7,
                          arrow=arrow(angle = 15,length=unit(0.15, "inches")),alpha=0.9,show.legend=F)+
          geom_path(data = fn.circ(pts=seq(0, stop*2*pi/n.samples, l = 100)), aes(x, y),size=5,col=circ.cl))
      }
    }
  })
},movie.name="Animation_Population.gif",interval=0.3,loop =1)})   #takes 0.4 sec per nrow of D



#Zoom in into Western Australia
Aus.lng=c(130,130,130,120)
Aus.lat=c(-20,-20,-25,-25)
Lngs=seq(-33,Aus.lng[1],by=5)
N.Shks.Int=8
ZOOMS=c(rep(5,N.Shks.Int),4,3,2,rep(2,length(Lngs)),
        2,3,4,5)
Lngs=c(rep(-35,11),Lngs)
COOR=data.frame(lon=c(Lngs,Aus.lng),
                lat=c(rep(-7,length(Lngs)),Aus.lat))
NN=1
saveGIF({
  
  for(x in 1:nrow(COOR))
  {
    p <- ggmap(get_map(COOR[x,],maptype="satellite", zoom = ZOOMS[x])) 
    if(x%in%1:N.Shks.Int) p=p+ annotate("text", x = -35, y=-5, label = "Sharks",col="white",size=8)+
        annotate("text", x = -35, y=-7, label = "International",col="white",size=8)
    if(x==nrow(COOR)) p=p+ annotate("text", x = 125, y=-26, label = "Western Australia",col="white",size=8)
    print(p+theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14,face="bold")))
  }
},movie.name="Animation_Zoom.gif",interval=.1,loop =NN) 



#Plotting
PlT="NO"
if(PlT=="YES")
{
  XLIM=c(112,120)
  YLIM=c(-36,-19)
  bb=Samp.per.age/2
  ID.list=seq(1,n.samples,by=bb)
  ColL=rainbow(length(ID.list))
  i=1
  plotMap(worldLLhigh, xlim=XLIM, ylim=YLIM,col="grey85")
  for(l in 1:(length(ID.list)-1))
  {
    Id=ID.list[l]:ID.list[l+1]
    with(CRW.trajectory[[i]],text(lon[Id],lat[Id],Id,col=ColL[l],cex=.5))
  }
}