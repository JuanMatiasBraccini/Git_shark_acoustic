##### SIMULATION EVALUATION OF SHARK MOVEMENT AND DETECTION BY CURTAINS OF ACOUSTIC ARRAYS #########
#notes: script for simulation movement of sharks according to different individual-based movement
#      models and testing if a curtain array can pick each movement model.
#      Steps:
#           1. Generate random distributions of key movement parameters: speed and bearing
#           2. Draw samples from distributions and simulate movement trajectories
#           3. Overlap  receiver curtains and detect hits 
#           4. Calculate speeds and bearings from detections 
#           5. Compute detected speed and bearing distributions 
#           6. Compare simulated and observed distributions


#sources: Nams 2006; Reynolds & Rhodes 2009; Tremblay et al 2009 

#notes:   consider incorporating depth as another parameter to sample from.  
      
#sensitivity: test the effects of tagging different number of sharks, monitoring for different number
              # of days, having 2 or 3 curtains, having sharks swimming a different speeds

#assumptions: all individuals move independently, accordingly to a specific movement model
#             no influence of habitat on movement (i.e. habitat is homogeneous, other than directed
#             movement that implies and implicit heterogeneous habitat)
#             minimum overlap between receivers with no gap
#             detection range is fixed across receivers
#             no receiver failure



#### Library section ####

rm(list=ls(all=TRUE))
setwd("H:/Matias WA Fisheries/Analyses/Acoustic_tagging/Simulation evaluation")



library(MASS)           #for fitting distributions
library(geosphere)      #for spatial statistics and trigonometry in a sphere
library(mixtools)       #for fitting bimodal normal distributions
library(prada)          #for fitting bivariate normal distribution
# To install "prada"
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("prada")
library(adehabitatLT)   #for different types of random walks
library(VGAM)           #for fitting a Levy distribution
library(rmutil)         #for drawing samples of Levy distribution

library(CircStats)      #for circular distributions
library(splancs)        #for testing if location is within boundaries
library(triangle)       #for triangular distribution


memory.limit(3900)   #set memory limit to maximum size

#### Data section ####


#### Parameter section ####

# sampling interval for simulated movement (in minutes)
#note: 10 mins sampling assumes that at most, a shark will move at most 600 in a straight line,
#      considering the speed distribution used
samp.int=10



#number of days monitored
ndays=365

# number of steps in each trajectory (in samp.int units)
n.steps=ndays*24*60/samp.int

# number of random samples from each pdfs
n.samples=n.steps 

# time step of trajectory (in seconds)
Time=60*samp.int

# number of simulated trajectories (i.e. number of tagged sharks)
n.trajectories=10

# speed range (in m/s)
min.speed=0.1
max.speed=1
max.speed.CRW=0.5   #lower max speed for CRW to constrain movement to same scale of other walks
max.speed.Levy=5    #larger max speed for Levy to allow movement within same scale of othe walks

# seeding location of each trajectory (assume single tagging event so release loc are similar)
Lat.range=c(10,10.1)
Long.range=c(110,110.1)
Seed.location.Lat=runif (n.trajectories,Lat.range[1],Lat.range[2])
Seed.location.Long=runif (n.trajectories,Long.range[1],Long.range[2])
Seed.location=data.frame(Long=Seed.location.Long,Lat=Seed.location.Lat)

# tolerance of movement (to constrain movement within polygon) 
max.mov.long=1/3  #in degrees. Average width of continental shelves is about 80 km (50 mi). Wikipedia.
max.mov.lat=1
Long.range.extended=c(Long.range[1]-max.mov.long,Long.range[2]+max.mov.long)
Lat.range.extended=c(Lat.range[1]-max.mov.lat,Lat.range[2]+max.mov.lat)
boundary=matrix(c(Long.range.extended,rev(Long.range.extended),rep(Lat.range.extended[1],2),
                  rep(Lat.range.extended[2],2)),ncol=2)


# speed (in m/s) and bearing (in degrees, with 0º being true N) ranges
  #Brownian random walk
speed.Brown=c(min.speed,max.speed)
bearing.Brown=c(0, 360)

  #Correlated random walk
speed.CRW=c(min.speed,max.speed.CRW)
crit.angle=45     #critical angle (see Reynolds and Rhodes 2009). Maximum changing angle from previous step
bearing.CRW=c(-crit.angle,crit.angle)

  #Levy flights
mu=2    #corresponding to optimal Levy flight search patterns (Bradshaw et al 2007)
speed.Levy=c(min.speed,max.speed.Levy)
length.out.speed=1000
speed.range=seq(speed.Levy[1],speed.Levy[2],length.out=length.out.speed)

  #Directed turning direction
speed.DTD=5   
bearing.DTD=c(60*pi/180,240*pi/180)
kappa.DTD=10
prop=0.5


# receiver pars
n.curtain=2         #number of receiver curtains          
detect.range=500    #receiver detection range (in metres)
ping.freq=60        #pinging frequency (in seconds)

n.intermedios=(60*samp.int/ping.freq)-1 #number of steps in between 10 mins steps

maxtime=20       #maximum time (in minutes) for consecutive hits (from trigonometry?)

pmax=1            #parameters of detection range sigmodeal pdf (see How & De Lestang in press)
D95=detect.range
D50=0.74*D95

#dist.btwn.curt=2*detect.range/(60*1853)    #(in degrees). twice detection range, giving marginal 
                                                          # overlap, it creates huge dataset in #3
dist.btwn.curt=2*detect.range/(60*1853)    #10 times detection range




#### Procedure section ####


#set up a loop for each walk as cannot do them all at once due to memory limitations

Start.syst=date()
Modelled.Walks=c("Brownian","CRW","Levy","DTD")

for(w in 1:length(Modelled.Walks))
{
# 1. Generate random distributions

if(Modelled.Walks[w]=="Brownian")
{
   #1.1.Brownian random walk
  pdf.speed.Brown=matrix(c(runif(n.samples*n.trajectories,speed.Brown[1],speed.Brown[2])),ncol=n.trajectories)
  pdf.bearing.Brown=matrix(c(runif(n.samples*n.trajectories,bearing.Brown[1],bearing.Brown[2])),ncol=n.trajectories)
  
  first.pdf.bearing.Brown=pdf.bearing.Brown[1,]
  nrow.pdf.bearing.Brown=nrow(pdf.bearing.Brown)
  ncol.pdf.bearing.Brown=ncol(pdf.bearing.Brown)
}
 

if(Modelled.Walks[w]=="CRW")
{
   #1.2.Correlated random walk
  #note: bearing in each step depends on the previous step plus a random term
  pdf.speed.CRW=matrix(c(runif(n.samples*n.trajectories,speed.CRW[1],speed.CRW[2])),ncol=n.trajectories)
  Rand.term.bearing=matrix(c(runif(n.samples*n.trajectories,bearing.CRW[1],bearing.CRW[2])),ncol=n.trajectories)
  pdf.bearing.CRW=matrix(NA,nrow=nrow.pdf.bearing.Brown,ncol=ncol.pdf.bearing.Brown)
  pdf.bearing.CRW[1,]=first.pdf.bearing.Brown

  for(i in 2:nrow(pdf.bearing.CRW))
  {
    pdf.bearing.CRW[i,]=pdf.bearing.CRW[i-1,]+Rand.term.bearing[i,]
    pdf.bearing.CRW[i,]=ifelse(pdf.bearing.CRW[i,]>360,pdf.bearing.CRW[i,]-360,pdf.bearing.CRW[i,])
    pdf.bearing.CRW[i,]=ifelse(pdf.bearing.CRW[i,]<0,360+pdf.bearing.CRW[i,],pdf.bearing.CRW[i,])
  }
  rm(Rand.term.bearing)
}

if(Modelled.Walks[w]=="Levy")
{
  #1.3.Truncated Levy walk (walk, not flight, see Reynolds and Rhodes 2009)
  #note: this is a truncated power-law distribution, as in Plank and Codling 2009 Eq. 2
  dist.speed.Levy=((mu-1)/(speed.Levy[1]^(1-mu)-speed.Levy[2]^(1-mu)))*(speed.range^-mu)
  dist.speed.Levy=dist.speed.Levy/sum(dist.speed.Levy)
  pdf.speed.Levy=matrix(c(sample(speed.range, n.samples*n.trajectories,replace=T,prob =dist.speed.Levy)),
                       ncol=n.trajectories)
  pdf.bearing.Levy=matrix(c(runif(n.samples*n.trajectories,bearing.Brown[1],bearing.Brown[2])),
                        ncol=n.trajectories)  # bearing is drawn from uniform distribution (Plank and
                                              # Codling 2009, Bradshaw et al 2007)
}

if(Modelled.Walks[w]=="DTD")
{
  #1.4.Directed turning direction
  #note: modelled through a biased random walk with very steep bimodal von Mises distribution for bearing
  #       add 0.1 to have this as minimum speed
  pdf.speed.DTD=matrix(rexp(n.samples*n.trajectories,speed.DTD)+0.1,ncol=n.trajectories)
  pdf.bearing.DTD=matrix(c(rmixedvm(n.samples*n.trajectories, bearing.DTD[1], bearing.DTD[2], kappa.DTD,
                        kappa.DTD, prop)*180/pi),ncol=n.trajectories)
}



# 2. Draw samples from pdfs and generate trajectories

  #2.1. Function for drawing samples and constructing trajectory
Samples=function(speed.sample,bearing.sample,Time)
{
    #calculate distance moved per time step (in metres)
  distance=speed.sample*Time

    #create dummy to store number of iterations outside boundaries
  prop.out<<-NULL  

    #compute paths using recursive calculation of position
  simulate=function(x)
  {
    position=as.numeric(Seed.location[x,])
    bb=bearing.sample[,x]
    dd=distance[,x]

    next.point=NULL
    next.point=rbind(next.point,position)
    for (j in 2:n.steps)
      {
        #calculate position for each step based on bearing and speed
        datos=destPoint(next.point[j-1,], b=bb[j], d=dd[j])
        
        #recalculate if position outside boundaries 
        traj=rbind(next.point[j-1,],datos)
        test=inout(traj,boundary,bound=T)
        prop.out<<-c(prop.out,sum(test[2]==F))
        
        if(test[2]==F)
        { repeat 
            {
              datos=destPoint(next.point[j-1,], b=rvm(1,0,k=0.010)*180/pi, d=dd[j])
              traj=rbind(next.point[j-1,],datos)
              test=inout(traj,boundary,bound=T)
              if(test[2]==T)break
            }
        }
        #combine
        next.point=rbind(next.point,datos)
      }
    rownames(next.point)=1:n.steps
    return(next.point)
  }

  paths <- lapply(1:n.trajectories, simulate)
  return(paths)
}  

  #2.2. Function for calculating distance travelled per trajectory
dist.travel=function(walk)
{
  Dist.Travel_km=Area.Travel_km2=NULL 
  for(i in 1:length(walk))
  {
    Next.point=as.data.frame(walk[[i]])
    next.point1=Next.point[-1,]
    colnames(next.point1)=c("Long1","Lat1")
    next.point=Next.point[-n.steps,]
    Next.point=cbind(next.point,next.point1)
    Next.point$dist.trav_m=distCosine(Next.point[,1:2],Next.point[,3:4])
    Dist.Travel_km=rbind(Dist.Travel_km,sum(Next.point$dist.trav_m)/1000)
    Area.Travel_km2=rbind(Area.Travel_km2,areaPolygon(Next.point[,1:2])/1000^2)
  }
  return(list(Dist.Travel_km=Dist.Travel_km,Area.Travel_km2=Area.Travel_km2))
}


  #2.3. Function for Plotting trajectories
traj.col=rainbow(n.trajectories)
traj.col=traj.col[sample(n.trajectories)]#resample to avoid having too similar colors together
plot.trajectory=function(walk,walkType)
{  
  dummy=do.call(rbind,walk)
  Long=range(dummy[,1])
  Lat=range(dummy[,2])
  plot(0,type="n",xlab="Long (º)",ylab="Lat (º)",main=walkType,xlim=Long,ylim=Lat)
  for(i in 1:length(walk))
  {
    Next.point=walk[[i]]
    next.point1=Next.point[-1,]
    next.point=Next.point[-n.steps,]
    arrows(next.point[,1],next.point[,2],next.point1[,1],next.point1[,2],col = traj.col[i],length=0.05)
    points(Seed.location[i,1],Seed.location[i,2],pch=15,col=traj.col[i], cex = 1.5) #start
    points(Next.point[n.steps,1],Next.point[n.steps,2],pch=17,col=traj.col[i], cex = 1.5)#end
  }
} 


  #2.4. Function for subsampling walks to generate different sampling periods
 Subsample=function(samp.period,walk)
{
  period=samp.period*24*60/samp.int
  datos=walk[1:period,]
  return(list(datos))  
}


  #2.5. Run function that creates trajectories
  #-- Brownian random walk --
if(Modelled.Walks[w]=="Brownian")
{
  Brownian.trajectory=Samples(pdf.speed.Brown,pdf.bearing.Brown,Time)
  Brownian.out=prop.out
  Brownian.out=100*sum(Brownian.out)/length(Brownian.out)    #percentage of iterations outside boundaries
    #export distributions and remove to free up memory
  write.csv(pdf.speed.Brown, file="pdf.speed.Brown.csv",row.names=F)
  write.csv(pdf.bearing.Brown, file="pdf.bearing.Brown.csv",row.names=F)
  rm(list=c("pdf.speed.Brown","pdf.bearing.Brown"))  
}

if(Modelled.Walks[w]=="CRW")
{
  #-- Correlated random walk --
  CRW.trajectory=Samples(pdf.speed.CRW,pdf.bearing.CRW,Time)
  CRW.out=prop.out
  CRW.out=100*sum(CRW.out)/length(CRW.out)
    #export distributions and remove to free up memory
  write.csv(pdf.speed.CRW, file="pdf.speed.CRW.csv",row.names=F)
  write.csv(pdf.bearing.CRW, file="pdf.bearing.CRW.csv",row.names=F)
  rm(list=c("pdf.speed.CRW","pdf.bearing.CRW"))  

}

if(Modelled.Walks[w]=="Levy")
{
  #-- Truncated Levy walk --
  Levy.trajectory=Samples(pdf.speed.Levy,pdf.bearing.Levy,Time)
  Levy.out=prop.out
  Levy.out=100*sum(Levy.out)/length(Levy.out)
    #export distributions and remove to free up memory
  write.csv(pdf.speed.Levy, file="pdf.speed.Levy.csv",row.names=F)
  write.csv(pdf.bearing.Levy, file="pdf.bearing.Levy.csv",row.names=F)
  rm(list=c("pdf.speed.Levy","pdf.bearing.Levy"))   
}

if(Modelled.Walks[w]=="DTD")
{
  #-- Directed turning direction --
  DTD.trajectory=Samples(pdf.speed.DTD,pdf.bearing.DTD,Time)
  DTD.out=prop.out
  DTD.out=100*sum(DTD.out)/length(DTD.out) 
    #export distributions and remove to free up memory
  write.csv(pdf.speed.DTD, file="pdf.speed.DTD.csv",row.names=F)
  write.csv(pdf.bearing.DTD, file="pdf.bearing.DTD.csv",row.names=F)
  rm(list=c("pdf.speed.DTD","pdf.bearing.DTD"))  
}



# Comment out when needing to subsample walks (SENSITIVITY TESTS)
#   #2.6. Run function to subsample walks 
# #-- Brownian random walk --
# Brownian.trajectory.30days.10min=mapply(Subsample,30,Brownian.trajectory)
# Brownian.trajectory.60days.10min=mapply(Subsample,60,Brownian.trajectory)
# Brownian.trajectory.6months.10min=mapply(Subsample,183,Brownian.trajectory)
# 
# #-- Correlated random walk --
# CRW.trajectory.30days.10min=mapply(Subsample,30,CRW.trajectory)
# CRW.trajectory.60days.10min=mapply(Subsample,60,CRW.trajectory)
# CRW.trajectory.6months.10min=mapply(Subsample,183,CRW.trajectory)
# 
# #-- Truncated Levy walk --
# Levy.trajectory.30days.10min=mapply(Subsample,30,Levy.trajectory)
# Levy.trajectory.60days.10min=mapply(Subsample,60,Levy.trajectory)
# Levy.trajectory.6months.10min=mapply(Subsample,183,Levy.trajectory)
# 
# #-- Directed turning direction --
# DTD.trajectory.30days.10min=mapply(Subsample,30,DTD.trajectory)
# DTD.trajectory.60days.10min=mapply(Subsample,60,DTD.trajectory)
# DTD.trajectory.6months.10min=mapply(Subsample,183,DTD.trajectory)
# 

# Comment out when required
#   #2.7. Run function that calculates distance and area travelled
#   #-- Brownian random walk --
# Brownian.dist.travel=dist.travel(Brownian.trajectory)
# 
#   #-- Correlated random walk --
# CRW.dist.travel=dist.travel(CRW.trajectory)
# 
#   #-- Truncated Levy walk --
# Levy.dist.travel=dist.travel(Levy.trajectory)
# 
#   #-- Directed turning direction --
# DTD.dist.travel=dist.travel(DTD.trajectory)
# 
# 
#   #2.8. Run function that plots trajectories
#   #-- Brownian random walk --
# plot.trajectory(Brownian.trajectory,"Brownian walk")
# 
#   #-- Correlated random walk --
# plot.trajectory(CRW.trajectory,"CRW")
# 
#   #-- Truncated Levy walk --
# plot.trajectory(Levy.trajectory,"LRW")
# 
#   #-- Directed turning direction --
# plot.trajectory(DTD.trajectory,"DTM")
# 


# 3. Overlap two receiver curtains over trajectories and detect hits

  #3.1. Overlap receiver curtains
    #calculate number of receivers per curtain
#Range.long=range(do.call(rbind,c(Brownian.trajectory,CRW.trajectory,Levy.trajectory,DTD.trajectory))[,1])
Range.long=Long.range.extended
Dist.rec.long=(Range.long[2]-Range.long[1])*111.2*1000  #in metres  
n.rec.per.curtain=ceiling(Dist.rec.long/(detect.range*2))

    #Determine position of receivers
LAT.Curtain1=mean(Seed.location.Lat)+dist.btwn.curt
LAT.Curtain2=mean(Seed.location.Lat)-dist.btwn.curt
LONG.Curt=seq(Range.long[1],Range.long[2],by=(detect.range*2)/(111.2*1000))
Curtain1=data.frame(Curt=1,Rec.num=1:length(LONG.Curt),Long=LONG.Curt,Lat=rep(LAT.Curtain1,length(LONG.Curt)))
Curtain2=data.frame(Curt=2,Rec.num=(length(LONG.Curt)+1):(2*length(LONG.Curt)),Long=LONG.Curt,Lat=rep(LAT.Curtain2,length(LONG.Curt)))
Curtains=rbind(Curtain1,Curtain2)
rm(list=c("Curtain1"))

write.csv(Curtain2, file="Curtain2.csv", row.names=F)


  #3.2 Calculate detections per receiver
#note: this calculates the position for every ping.freq (i.e. interpolate within the 10 min intervals)
#     for positions within a maximum distance from receivers and then calculates number of hits with 
#     a given probability based on the distance from the receiver

#detection probability (in m)
detect.prob=function(x)pmax/(1+exp(log(19)*((x-D50)/(D95-D50))))  #sigmoideal prob
# detect.prob=function(x)1-ptriangle(x, a=0, b=detect.range, c=0.000001)  #linearly decaying prob


threshold.BRW.CRW=detect.range*0.1+detect.range    #threshold calculated +/- 10% detection range
threshold.LRW=2000+detect.range       #these thresholds are larger due to larger displacement (due to 
                                      #larger max speed)
threshold.DTM=2000+detect.range

THRESHOLD=threshold.BRW.CRW

hits=function(walk)
{
  simulate=function(x)
  {
    #create matrices for vectorising
    START=as.data.frame(x)
    END=START[-1,]
    START=START[-nrow(START),]
    START.END=cbind(START,END)
    colnames(START.END)[3:4]=c("long1","lat1")
    #keep only those positions within threshold
    START.END=subset(START.END,lat>=(LAT.Curtain2-(THRESHOLD/(60*1853))) &
      lat<=(LAT.Curtain2+(THRESHOLD/(60*1853))) |
      lat<=(LAT.Curtain1+(THRESHOLD/(60*1853))) & lat>=(LAT.Curtain1-(THRESHOLD/(60*1853))))
    colnames(START.END)[3:4]=c("lon","lat")

    if(nrow(START.END)>0)
    {
      #extract time
      time.i=as.numeric(rownames(START.END))

      #interpolate intermediate points  (i.e. location by minute)
      gci <- gcIntermediate(START.END[,1:2], START.END[,3:4],n=n.intermedios)
      
      #combine interpolated points and times
      Expand=function(i)
      {
        #create matrices for vectorising
        Datos=rbind(START.END[i,1:2],as.data.frame(gci[[i]]),START.END[i,3:4])
        Datos$time=seq(time.i[i]*samp.int,length.out=nrow(Datos))
        Datos=as.data.frame(cbind(Datos[-nrow(Datos),],Datos[-1,]))
        Datos$delta.t=Datos[,6]-Datos[,3]
        Datos=Datos[,-3]
        return(Datos)  
      }

      combined=lapply(1:length(gci), Expand)
      rm(gci)
      START.END=do.call(rbind,combined)
      rm(combined)

      rownames(START.END)=START.END$time      #crucial for allocating time step time

      #re-calculate distance from all receivers for each trajectory point
      Final.Distance=apply(START.END, MARGIN=1, function(x) distCosine(c(x[1],x[2]),Curtains[,3:4]))
      Final.Distance=t(Final.Distance)
      rm(list=c("START.END","START","END"))
      
      #select receivers detecting sharks    
      Final.Distance.in.out=ifelse(Final.Distance<=detect.range,1,0)
      xxx=which(rowSums(Final.Distance.in.out)>0)
      yyy=which(colSums(Final.Distance.in.out)>0)                               
      these.receivers=as.data.frame(Final.Distance[xxx,yyy])
      Final.Distance.in.out=as.data.frame(Final.Distance.in.out[xxx,yyy])
      colnames(these.receivers)=colnames(Final.Distance.in.out)=yyy
      
      #remove double hits
      doube.hit=NULL
      double.hit=which(rowSums(Final.Distance.in.out)>1)
      if(length(double.hit)>0)
      {
        mincols=apply(these.receivers[double.hit,],1,function(x) match(min(x),x))        
        Final.Distance.in.out[double.hit,]=0
        for(i in 1:length(mincols))
        {
          Final.Distance.in.out[double.hit[i],mincols[i]]=1
        }
      }

      #calculate hits by sampling with given probability based on distance from receiver
        #first calculate probability given triangular distribution
      if(ncol(these.receivers)>0)
      {
        probs=matrix(detect.prob(as.matrix(these.receivers)),nrow=nrow(these.receivers),
              ncol=ncol(these.receivers))
       
          #then sample with the given probability calculated in previous step
        Probabilities=matrix(mapply(probs,FUN=function(x)rbinom(1, 1,prob=x)),nrow=nrow(probs),
               ncol=ncol(probs))
    
          #then multiply (element-wise) matrices
        Detections=Final.Distance.in.out*Probabilities
          
        #remove again 0 rows and columns
        Detections=Detections[which(rowSums(Detections)>0),which(colSums(Detections)>0)]
  
        #assign position of the receivers that detected sharks
         Dummy=matrix(rep(as.numeric(colnames(Detections)),each=nrow(Detections)),nrow=nrow(Detections),
                ncol=ncol(Detections))
        Detections.pos=Detections*Dummy
        Detections.pos=rowSums(Detections.pos)
        longs.lats=Curtains[Detections.pos,]
        rownames(longs.lats)=rownames(Detections)
        

         #add corresponding time  (in minutes)
         longs.lats$time=as.numeric(rownames(longs.lats))
        return(longs.lats)
        }
      }

   }  
  hits.Data=0
  hits.Data <- lapply(walk, simulate)
  return(hits.Data)
}

if(Modelled.Walks[w]=="Brownian")
{
  #-- Brownian random walk --
  Brownian.hits=hits(Brownian.trajectory)
  #export trajectories for plotting and freeing up memory
  if(n.steps==52560)
  {
    write.csv(Brownian.trajectory[[1]][1:20000,], file="Brownian.trajectory.csv", row.names=F)
  }
  rm(Brownian.trajectory)  
}

if(Modelled.Walks[w]=="CRW")
{
  #-- Correlated random walk --
  CRW.hits=hits(CRW.trajectory)
  if(n.steps==52560)
  {
  write.csv(CRW.trajectory[[1]][1:10000,], file="CRW.trajectory.csv", row.names=F)    
  }
  rm(CRW.trajectory)  
}

if(Modelled.Walks[w]=="Levy")
{
  #-- Truncated Levy walk --
  Levy.hits=hits(Levy.trajectory)
  if(n.steps==52560)
  {
    write.csv(Levy.trajectory[[1]][1:10000,], file="Levy.trajectory.csv", row.names=F)    
  }
  rm(Levy.trajectory)  
}

if(Modelled.Walks[w]=="DTD")
{
  #-- Directed turning direction --
  DTD.hits=hits(DTD.trajectory)
  if(n.steps==52560)
  {
    write.csv(DTD.trajectory[[1]][1:10000,], file="DTD.trajectory.csv", row.names=F)    
  }
  rm(DTD.trajectory)
}


# 4. Calculate speeds and bearings from detections

#calculate likely long and lat given detection range and max movement probabilities
Dist.sequence=seq(0,detect.range,by=1)
these.probs=pmax/(1+exp(log(19)*((Dist.sequence-D50)/(D95-D50))))
these.probsx2=c(rev(these.probs),these.probs[-1])
ntimes=length(these.probs)

speed.bear=function(walk)
{
  these=c(bearing=NA,speed=NA)
  datos=walk
  if(nrow(as.data.frame(datos))>0)
  {
    #Vectorise data to avoid looping 
    datos=data.frame(datos[-nrow(datos),],datos[-1,])

    #calculate if hits occur at different receiver
    datos$Different=with(datos,ifelse(Rec.num==Rec.num.1,"NO","YES"))
  
    #calculate consecutive time                      
    datos$delta.t.1=datos$time.1-datos$time
    
    #extract consecutive hits at different receivers
    datos=subset(datos,datos$delta.t.1<=maxtime & datos$Different=="YES")
    
    
    if(nrow(datos)>0)
    {
      #allocate appropriate speed distribution
       if(Modelled.Walks[w]=="Brownian")
       {
         speed=runif(nrow(datos),speed.Brown[1],speed.Brown[2])
       }
       if(Modelled.Walks[w]=="CRW")
       {
         speed=runif(nrow(datos),speed.CRW[1],speed.CRW[2])
       }
       if(Modelled.Walks[w]=="Levy")
       {
         speed=sample(speed.range, nrow(datos),replace=T,prob =dist.speed.Levy)
       }
       if(Modelled.Walks[w]=="DTD")
       {
         speed=rexp(nrow(datos),speed.DTD)+0.1  
       }

      #remove the "different" column
      datos=datos[,-(match("Different",colnames(datos)))]
      
      #calculate likely maximum distance between receivers given speed distribution
      datos$max.distance=datos$delta.t.1*60*speed
     
      resampled=function(x)     #sigmoidal distribution, conditioned to max.distance
        { 
          x=as.matrix(x)
          dist=distCosine(c(x[3],x[4]),cbind(x[8],x[9]))
          extreme=x[12]/2
          start.long=which.min(abs(Dist.sequence-(abs(detect.range-extreme))))
          if(x[3]>x[8])
          {
            LONG.seq=seq(x[3]-(start.long/(60*1853)),x[3]-detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
            LONGS=sample(LONG.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
            LONG1.seq=seq(x[8]+(start.long/(60*1853)),x[8]+detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
            LONGS.1=sample(LONG1.seq,10000,prob=these.probs[start.long:detect.range],replace=T)

            if(x[4]==x[9])
            {
              LAT.seq=seq(x[4]-(detect.range/(60*1853)),x[4]+detect.range/(60*1853),
                         length.out=length(these.probsx2))
              LATS=sample(LAT.seq,10000,prob=these.probsx2,replace=T)
              LAT1.seq=seq(x[9]-(detect.range/(60*1853)),x[9]+detect.range/(60*1853),
                         length.out=length(these.probsx2))
              LATS.1=sample(LAT1.seq,10000,prob=these.probsx2,replace=T)
            }

            if(x[4]>x[9])
            {
              LAT.seq=seq(x[4]-(start.long/(60*1853)),x[4]-detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LATS=sample(LAT.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
              LAT1.seq=seq(x[9]+(start.long/(60*1853)),x[9]+detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LATS.1=sample(LAT1.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
            }

            if(x[4]<x[9])
            {
              LAT.seq=seq(x[4]+(start.long/(60*1853)),x[4]+detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LATS=sample(LAT.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
              LAT1.seq=seq(x[9]-(start.long/(60*1853)),x[9]-detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LAT.1=sample(LAT1.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
            }
          }
          
          if(x[3]<x[8])
          {
            LONG.seq=seq(x[3]+(start.long/(60*1853)),x[3]+detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
            LONGS=sample(LONG.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
            LONG1.seq=seq(x[8]-(start.long/(60*1853)),x[8]-detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
            LONGS.1=sample(LONG1.seq,10000,prob=these.probs[start.long:detect.range],replace=T)

            if(x[4]==x[9])
            {
              LAT.seq=seq(x[4]-(detect.range/(60*1853)),x[4]+detect.range/(60*1853),
                         length.out=length(these.probsx2))
              LATS=sample(LAT.seq,10000,prob=these.probsx2,replace=T)
              LAT1.seq=seq(x[9]-(detect.range/(60*1853)),x[9]+detect.range/(60*1853),
                         length.out=length(these.probsx2))
              LATS.1=sample(LAT1.seq,10000,prob=these.probsx2,replace=T)
            }

            if(x[4]>x[9])
            {
              LAT.seq=seq(x[4]-(start.long/(60*1853)),x[4]-detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LATS=sample(LAT.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
              LAT1.seq=seq(x[9]+(start.long/(60*1853)),x[9]+detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LATS.1=sample(LAT1.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
            }

            if(x[4]<x[9])
            {
              LAT.seq=seq(x[4]+(start.long/(60*1853)),x[4]+detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LATS=sample(LAT.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
              LAT1.seq=seq(x[9]-(start.long/(60*1853)),x[9]-detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LAT.1=sample(LAT1.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
            }
          }

          if(x[3]==x[8])
          {
            LONG.seq=seq(x[3]-(detect.range/(60*1853)),x[3]+detect.range/(60*1853),
                         length.out=length(these.probsx2))
            LONGS=sample(LONG.seq,10000,prob=these.probsx2,replace=T)
            LONG1.seq=seq(x[8]-(detect.range/(60*1853)),x[8]+detect.range/(60*1853),
                         length.out=length(these.probsx2))
            LONGS.1=sample(LONG1.seq,10000,prob=these.probsx2,replace=T)
      
            if(x[4]>x[9])
            {
              LAT.seq=seq(x[4]-(start.long/(60*1853)),x[4]-detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LATS=sample(LAT.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
              LAT1.seq=seq(x[9]+(start.long/(60*1853)),x[9]+detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LATS.1=sample(LAT1.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
            }
            
            if(x[4]<x[9])
            {
              LAT.seq=seq(x[4]+(start.long/(60*1853)),x[4]+detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LATS=sample(LAT.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
              LAT1.seq=seq(x[9]-(start.long/(60*1853)),x[9]-detect.range/(60*1853),
                         length.out=length(Dist.sequence[start.long:detect.range]))
              LAT.1=sample(LAT1.seq,10000,prob=these.probs[start.long:detect.range],replace=T)
            }
          }
          
          dist=distCosine(cbind(LONGS,LATS),cbind(LONGS.1,LATS.1))
          nn=length(LONGS)
          DATA=data.frame(Curt=rep(x[1],nn),Rec.num=rep(x[2],nn),Long=LONGS,Lat=LATS,time=rep(x[5],nn),
                          Curt.1=rep(x[6],nn),Rec.num.1=rep(x[7],nn),Long.1=LONGS.1,Lat.1=LATS.1,
                          time.1=rep(x[10],nn),delta.t.1=rep(x[11],nn),max.distance=rep(x[12],nn))
          DATA=DATA[dist<=x[12],]
          return(DATA)
        }

      #resampled=function(x) rtriangle(1, a=x-detect.range/(60*1853), b=x+detect.range/(60*1853), c=x)
      DATOS=apply(datos,1,resampled)
      datos=do.call(rbind,DATOS)
      rm(DATOS)
         
      #calculate distance, speed and bearing between consecutive hits   
      datos$distance.m=with(datos,distCosine(cbind(Long,Lat),cbind(Long.1,Lat.1))) #distance in metres
      datos$bearing=with(datos,bearing(cbind(Long,Lat),cbind(Long.1,Lat.1)))    #bearing in degrees (N = 0, E = 90)
      datos$speed=datos$distance.m/(datos$delta.t.1*ping.freq)            #speed, in m/s
  
      these=data.frame(bearing=datos$bearing,speed=datos$speed)
    }
  }
  return(these)
}


if(Modelled.Walks[w]=="Brownian")
{
  #-- Brownian random walk --
  Brownian.obs.speed.bear=lapply(Brownian.hits,speed.bear)
  rm(Brownian.hits)
  
}

if(Modelled.Walks[w]=="CRW")
{
  #-- Correlated random walk --
  CRW.obs.speed.bear=lapply(CRW.hits,speed.bear)
  rm(CRW.hits)
}

if(Modelled.Walks[w]=="Levy")
{
  #-- Truncated Levy walk --
  Levy.obs.speed.bear=lapply(Levy.hits,speed.bear)
  rm(Levy.hits)
}

if(Modelled.Walks[w]=="DTD")
{
  #-- Directed turning direction --
  DTD.obs.speed.bear=lapply(DTD.hits,speed.bear)
  rm(DTD.hits)
}



# 5. Compute overall observed distributions of speed and bearing
Distributions=function(database)
{
  newdata=do.call(rbind,database)

  #remove nonsense speeds
#  newdata=subset(newdata,newdata$speed<=max.speed)

  #remove NAs
  newdata=newdata[!(is.na(newdata[,1])),]

  return(newdata)  
}

if(Modelled.Walks[w]=="Brownian")
{
  #-- Brownian random walk --
  Brownian.obs.distribution=Distributions(Brownian.obs.speed.bear)
  #export and free up space
  write.csv(Brownian.obs.distribution, file="Brownian.obs.distribution.csv", row.names=F)
  rm(list=c("Brownian.obs.speed.bear","Brownian.obs.distribution"))
}

if(Modelled.Walks[w]=="CRW")
{
  #-- Correlated random walk --
  CRW.obs.distribution=Distributions(CRW.obs.speed.bear)
  write.csv(CRW.obs.distribution, file="CRW.obs.distribution.csv", row.names=F)
  rm(list=c("CRW.obs.speed.bear","CRW.obs.distribution"))
}

if(Modelled.Walks[w]=="Levy")
{
  #-- Truncated Levy walk --
  Levy.obs.distribution=Distributions(Levy.obs.speed.bear)
  write.csv(Levy.obs.distribution, file="Levy.obs.distribution.csv", row.names=F)
  rm(list=c("Levy.obs.speed.bear","Levy.obs.distribution"))
}


if(Modelled.Walks[w]=="DTD")
{
  #-- Directed turning direction --
  DTD.obs.distribution=Distributions(DTD.obs.speed.bear)
  write.csv(DTD.obs.distribution, file="DTD.obs.distribution.csv", row.names=F)
  rm(list=c("DTD.obs.speed.bear","DTD.obs.distribution"))
}



# 6. Compare simulated and observed distributions


}

End.syst=date()
print(Start.syst)
print(End.syst)


