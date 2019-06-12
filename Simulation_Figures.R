#Script for producing Figures from Simulation Evaluation of acoustic data

#MISSING:Fig 1.C. show full lines but only hits within detection range


#DATA SECTION
setwd("C:/Matias/Analyses/Acoustic_tagging/Simulation evaluation")

#Fig1.A
pdf.speed.Brown=as.matrix(read.csv(file="pdf.speed.Brown.csv"))
pdf.speed.CRW=as.matrix(read.csv(file="pdf.speed.CRW.csv"))
pdf.speed.Levy=as.matrix(read.csv(file="pdf.speed.Levy.csv"))
pdf.speed.DTD=as.matrix(read.csv(file="pdf.speed.DTD.csv"))
pdf.bearing.Brown=as.matrix(read.csv(file="pdf.bearing.Brown.csv"))
pdf.bearing.CRW=as.matrix(read.csv(file="pdf.bearing.CRW.csv"))
pdf.bearing.Levy=as.matrix(read.csv(file="pdf.bearing.Levy.csv"))
pdf.bearing.DTD=as.matrix(read.csv(file="pdf.bearing.DTD.csv"))

#Fig1.B
Brownian.trajectory=as.matrix(read.csv(file="Brownian.trajectory.csv"))
CRW.trajectory=as.matrix(read.csv(file="CRW.trajectory.csv"))
Levy.trajectory=as.matrix(read.csv(file="Levy.trajectory.csv"))
DTD.trajectory=as.matrix(read.csv(file="DTD.trajectory.csv"))

#Fig1.C
Walk.Fig1C=read.csv(file='e.g.walk.Fig1C.csv',sep = ",")
Curtain2=read.csv(file='Curtain2.csv',sep = ",")
Curtain2$Lat=9.9
  
#Fig1.D
Brownian.obs.distribution=read.csv(file="Brownian.obs.distribution.csv")
CRW.obs.distribution=read.csv(file="CRW.obs.distribution.csv")
Levy.obs.distribution=read.csv(file="Levy.obs.distribution.csv")
DTD.obs.distribution=read.csv(file="DTD.obs.distribution.csv")

#PROCEDURE SECTION


def.par <- par(no.readonly = TRUE) # save default, for resetting...

#Figure 1.
Shark.Traj=c("Shark 1")
nf <- layout(matrix(c(1,1,2,2,3,4,5,6,0,7,7,0,8,8,9,9), 4, 4, byrow = TRUE), 
             widths=rep(19,4), heights=rep(20,5), TRUE)
layout.show(nf)
par(oma=c(0.1,0.1,0.1,0.1))
LEGENDS=c("Brownian","CRW","LRW","DTM")

# Figure 1.A)
par(mar=c(0.2,0.3,0.2,0.3))  
plot(density(c(pdf.speed.Brown)),xlab="",ylab="",main="",xlim=c(0.1,1),xaxt='n',
     yaxt='n',ylim=c(0,7),lwd=2)
lines(density(c(pdf.speed.CRW)),lty=2,lwd=2,col="gray46")
lines(density(c(pdf.speed.Levy)),lty=1,col="gray46",lwd=3)
lines(density(c(pdf.speed.DTD)),lty=4, lwd=3)
legend("topright",LEGENDS,col=c(1,"gray46","gray46",1),lty=c(1,2,1,4),bty='n',cex=0.9,lwd=c(2,2,2,3))
title("Speed (m/s)",line = -0.8,cex.main =0.9)
text(0.01,6.072, "A", xpd=NA, pos = 3,cex=1.75)



plot(density(c(pdf.bearing.Brown)),xlab="",ylab="",main="",xlim=c(0,360),xaxt='n',
     yaxt='n',ylim=c(0,0.0105),lwd=2)
lines(density(c(pdf.bearing.CRW)),col="transparent")
lines(density(c(pdf.bearing.Levy)),col="transparent")
lines(density(c(pdf.bearing.DTD)),lty=4, lwd=3)
legend("topright",LEGENDS,col=c(1,"white","white",1),lty=c(1,2,1,4),bty='n',cex=0.9,lwd=c(2,2,2,3))
title("Bearing (degrees)",line = -0.8,cex.main =0.9)



# Figure1.B)
Long.Fig1.walk=c(109.85,110.6); Lat.Fig1.walk=c(9.37,10.3)
plot.Fig1.walk.=function(walk,Title)
{  
    Next.point=walk
    LENGTH=nrow(walk)
    next.point1=Next.point[-1,]
    next.point=Next.point[-nrow(walk),]
    Lim.x=c(min(Next.point[,1]),max(Next.point[,1]))
    Lim.y=c(min(Next.point[,2]),max(Next.point[,2]))

    if(!(Title==Shark.Traj[1])) 
    {
      plot(0,type="n",xlab="",ylab="",main="",xlim=Lim.x,ylim=Lim.y,xaxt='n',yaxt='n')
      #title(Title,line = -0.8,cex.main =0.9)
      legend("topleft",Title,bty="n")
      segments(next.point[,1],next.point[,2],next.point1[,1],next.point1[,2],lty=1,col=1,lwd=0.8)
      points(walk[1,1],walk[1,2],pch=15,col=2, cex = 1.5) #start
      points(walk[LENGTH,1],walk[LENGTH,2],pch=17,col=3, cex = 1.5)#end
    
      if(Title==LEGENDS[1])text(109.8875,10.1, "B", xpd=NA,pos = 3,cex=1.75)
    }
    
    if(Title==Shark.Traj[1]) 
    {
      Lim.x=c(109.8774, 109.9013)
      Lim.y=c( 9.899622,9.900333)
      plot(0,type="n",xlab="",ylab="",main="",xlim=Lim.x,ylim=Lim.y,xaxt='n',yaxt='n')
      title("",line = -0.8,cex.main =0.9)
      points(Curtain2$Long,Curtain2$Lat,pch=21,col=4, cex = 15, bg="yellow")
#      points(Curtain2$Long,Curtain2$Lat,pch=21,col=4, cex = 23.5, bg="yellow")
      points(Curtain2$Long,Curtain2$Lat,pch=21,col=4, cex = 2,bg=4)
      arrows(next.point[,1],next.point[,2],next.point1[,1],next.point1[,2],lty=1,col=1,lwd=1.75,
             length=0.25)
      box()
      text(min(Lim.x)*0.99986,max(Lim.y)*0.99999, "C",xpd=NA,pos = 3,cex=1.75)
    }
} 

plot.Fig1.walk.(Brownian.trajectory,LEGENDS[1])
plot.Fig1.walk.(CRW.trajectory,LEGENDS[2])
plot.Fig1.walk.(Levy.trajectory,LEGENDS[3])
plot.Fig1.walk.(DTD.trajectory,LEGENDS[4])




#Figure 1.C)
plot.Fig1.walk.(Walk.Fig1C,Shark.Traj[1])




#Figure 1.D) replaces with observations
plot(density(Brownian.obs.distribution$speed),xlab="",ylab="",main="",xlim=c(0.1,1),xaxt='n',
     yaxt='n',ylim=c(0,3.5),lwd=2)
lines(density(CRW.obs.distribution$speed),lty=2,lwd=2,col="gray46")
lines(density(Levy.obs.distribution$speed),lty=1,col="gray46",lwd=3)
lines(density(DTD.obs.distribution$speed),lty=4, lwd=3)
legend("topright",LEGENDS,col=c(1,"gray46","gray46",1),lty=c(1,2,1,4),bty='n',cex=0.9,lwd=c(2,2,2,3))
title("Speed (m/s)",line = -0.8,cex.main =0.9)
text(0.0066,3.0, "D", xpd=NA, pos = 3,cex=1.75)


plot(density(Brownian.obs.distribution$bearing),xlab="",ylab="",main="",xlim=c(0,360),xaxt='n',
     yaxt='n',ylim=c(0,.005),lwd=2)
lines(density(CRW.obs.distribution$bearing),lty=2,lwd=2,col="gray46")
lines(density(Levy.obs.distribution$bearing),lty=1,col="gray46",lwd=3)
lines(density(DTD.obs.distribution$bearing),lty=4, lwd=3)
legend("topright",LEGENDS,col=c(1,"gray46","gray46",1),lty=c(1,2,1,4),bty='n',cex=0.9,lwd=c(2,2,2,3))
title("Bearing (degrees)",line = -0.8,cex.main =0.9)




par(def.par)#- reset to default
