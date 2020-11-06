#Bubble Plot
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
  mtext(ylab,2,las=3,line=2.25,cex=1.5)

  for(i in 1:n)
    {
      points(xo[i,],yo[i,],cex=zo[,i],pch=Mat.pch[,i],col=Mat.colors[,i])
    }
  points(Rel.date,1:ny,col=Rilis.col,cex=1.25,pch="R")
  

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

