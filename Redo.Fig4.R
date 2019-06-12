#Re do Simons
setwd("C:/Matias/Analyses/Acoustic_tagging/Dusky_migration/Re.do.Fig.4")
library(RODBC)
channel <- odbcConnectExcel2007("FL.xlsx") 
FL<- sqlFetch(channel,"dat")
close(channel)

channel <- odbcConnectExcel2007("Month.xlsx") 
Month<- sqlFetch(channel,"dat")
close(channel)


CI.pol=function(X,Ylow,Yhigh,COL,BORDER)
{
  XX=c(X,tail(X, 1),rev(X),X[1])
  YY <- c(Ylow, tail(Yhigh, 1), rev(Yhigh), Ylow[1])
  polygon(XX,YY,col=COL,border=BORDER)
}

male.col="grey40"
fem.col="grey40"
col.bg.F=rgb(.1,.1,.1,alpha=0.1)
col.bg.M=rgb(.1,.1,.1,alpha=0.1)

#fem.col=rgb(.6,0,0.1,alpha=0.6)
#col.bg.F=rgb(.75,0,0,alpha=0.1)
#male.col="blue"
#col.bg.M=rgb(0,0,0.75,alpha=0.1)

setwd("C:/Matias/Analyses/Acoustic_tagging/FRDC/Outputs_movement/Natal_migration/Paper")
tiff(file="Figure_prob_movement N_S_Simon.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfrow=c(2,2),mar=c(3,3,1.5,.1), oma=c(1,1.75,.1,1),las=1,mgp=c(1,0.8,0))

#Size effect

  #females
#January
a=subset(FL,sex=="F" & month==1)
plot(a$length,a$median,type='l',lwd=2,col=fem.col,xlab="",ylab="",cex.axis=1.5,ylim=c(0,1))
CI.pol(X=a$length,Ylow=a$lower95,Yhigh=a$upper95,col.bg.F,"transparent")
legend("topleft",c("January","August"),lty=c(2,1),lwd=3,col=c(fem.col,fem.col),bty="n",cex=1.5)

#August
a=subset(FL,sex=="F" & month==8)
lines(a$length,a$median,lty=2,lwd=2,col=fem.col)
CI.pol(X=a$length,Ylow=a$lower95,Yhigh=a$upper95,col.bg.F,"transparent")
#mtext("Migration probability",2,las=3,cex=1.75,line=3.1)
mtext("Females",3,cex=1.5,line=0)

#males
#January
a=subset(FL,sex=="M" & month==1)
plot(a$length,a$median,type='l',lwd=2,col=male.col,xlab="",ylab="",cex.axis=1.5,ylim=c(0,1))
CI.pol(X=a$length,Ylow=a$lower95,Yhigh=a$upper95,col.bg.M,"transparent")
legend("topleft",c("January","August"),lty=c(2,1),lwd=3,col=c(male.col,male.col),bty="n",cex=1.5)

#August
a=subset(FL,sex=="M" & month==8)
lines(a$length,a$median,lty=2,lwd=2,col=male.col)
CI.pol(X=a$length,Ylow=a$lower95,Yhigh=a$upper95,col.bg.M,"transparent")
mtext("Fork length (cm)                                                    ",1,cex=1.5,line=2.35)
mtext("Males",3,cex=1.5,line=0)



#Month effect
  
#females
a=subset(Month,Sex=="F")
plot(a$Month,a$Median,type='l',lwd=2,col=fem.col,xlab="",ylab="",cex.axis=1.5,ylim=c(0,1))
CI.pol(X=a$Month,Ylow=a$PER_Low_95,Yhigh=a$PER_Upp_95,col.bg.F,"transparent")
#mtext("Prob. of occurring north",2,las=3,cex=1.75,line=3.1)

#males
a=subset(Month,Sex=="M")
plot(a$Month,a$Median,type='l',lwd=2,col=male.col,xlab="",ylab="",cex.axis=1.5,ylim=c(0,1))
CI.pol(X=a$Month,Ylow=a$PER_Low_95,Yhigh=a$PER_Upp_95,col.bg.M,"transparent")

mtext("Month                                                      ",1,cex=1.5,line=2.35)
mtext("Prob. of migrating north",2,las=3,cex=1.5,line=0,outer=T)
dev.off()
