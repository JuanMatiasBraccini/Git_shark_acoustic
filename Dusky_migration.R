#Modelling dusky shark migration thru logistic glms

#note: use glms to quantify
#               1. the size at which males and females start migration
#               2. the timing of migration and the proportion of males and females undertaking migrations


#Read in data
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

setwd(handl_OneDrive("Analyses/Acoustic_tagging/FRDC/Outputs_movement/Natal_migration/Paper"))
dat.glm=read.csv("dat.glm.csv")
library(lme4)

#some functions
Most.common=function(what)
{
  ID=match(what,names(fixedTerms))
  Tabla=sort(table(fixedTerms[,ID]))
  return(names(Tabla[length(Tabla)]))
}
CI.pol=function(X,Ylow,Yhigh,COL,BORDER)
{
  XX=c(X,tail(X, 1),rev(X),X[1])
  YY <- c(Ylow, tail(Yhigh, 1), rev(Yhigh), Ylow[1])
  polygon(XX,YY,col=COL,border=BORDER)
}
fn.logis=function(dat,pmax,inflx,slop) pmax/(1+exp((dat-inflx)/slop))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Deviance.explained.R"))
fn.plt.pred.M=function(Preds,CL.lin,CL.back)
{
  Fit=Preds$fit
  SEs=Preds$se.fit
  Fit.plusSE=Fit+1.96*SEs
  Fit.minSE=Fit-1.96*SEs
  
  smooth.fem=smooth.spline(1:12,Fit,df=4)
  smooth.fem.plusSE=smooth.spline(1:12,Fit.plusSE,df=4)
  smooth.fem.minSE=smooth.spline(1:12,Fit.minSE,df=4)
   
  lines(smooth.fem$x,smooth.fem$y,lwd=3,col=CL.lin)
  CI.pol(X=smooth.fem.plusSE$x,Ylow=smooth.fem.minSE$y,
         Yhigh=smooth.fem.plusSE$y,CL.back,CL.lin)
}
fn.plt.pred.M.spline.only=function(Preds,CL.lin)
{
  Fit=Preds$fit  
  smooth.fem=smooth.spline(1:12,Fit,df=4)  
  lines(smooth.fem$x,smooth.fem$y,lwd=3,col=CL.lin,lty=2)
}

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/MS.Office.outputs.R"))
exp.tabl=function(Tbl,Doc.nm)
{
  fn.word.table(WD=getwd(),TBL=Tbl,Doc.nm=Doc.nm,caption=NA,paragph=NA,
          HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
          Zebra='NO',Zebra.col='grey60',Grid.col='black',
          Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")  
}
 
fn.see.month.glm.dat=function(DATA)
{
  DATA$CL=with(DATA,ifelse(S==1,"blue","red"))
  TGG=unique(DATA$TagCode)
  plot(1:12,1:1:12,ylim=c(0,length(TGG)),col="transparent",yaxt="n",ylab="",xlab="Month")
  for(p in 1:length(TGG))    
  {
    qq=subset(DATA,TagCode==TGG[p])
    points(qq$Month,rep(p,nrow(qq)),col=qq$CL,pch=19)
  }
  axis(2,1:length(TGG),TGG,las=1,cex=0.85)  
}

fn.which.50=function(Range,Target) which(abs(Range-Target)==min(abs(Range-Target)))
fn.CI=function(DAT1,DAT2,DAT3)
{
  ID=fn.which.50(DAT1,0.5)
  FL_50=round(100*DAT3[ID])
  ID=fn.which.50(DAT1-1.96*DAT2,0.5)
  FL_50_low_SE=round(100*DAT3[ID])
  ID=fn.which.50(DAT1+1.96*DAT2,0.5)
  FL_50_up_SE=round(100*DAT3[ID])
  
  return(data.frame(FL_50_low=FL_50_up_SE,FL_50=FL_50,FL_50_up=FL_50_low_SE))
}

#1. the size at which males and females start migration

dat.FL=subset(dat.glm,select=c(TagCode,Sex,FL,N))
TAb=table(dat.FL$TagCode,dat.FL$N)
TAb=ifelse(TAb>0,1,0)
TAb=as.data.frame(TAb)
TAb$TagCode=rownames(TAb)
names(TAb)[1:2]=c("S","N")
never.north=subset(TAb,S==1 & N==0)$TagCode
Never.north=subset(dat.FL,TagCode%in%never.north)
Never.north=Never.north[!duplicated(Never.north$TagCode),]
Moved.north=subset(dat.FL,!TagCode%in%never.north)
Moved.north=Moved.north[!duplicated(Moved.north$TagCode),]
dat.FL.glm=rbind(Never.north,Moved.north)
  # adjust two large males that were release in the North but only detected south
dat.FL.glm$N=with(dat.FL.glm,ifelse(FL%in%c(2.46,2.58) & Sex=="M",1,N)) 

  #see data used in model
with(subset(dat.FL.glm,Sex=="F"),plot(FL,N,pch=19,col="pink"))
with(subset(dat.FL.glm,Sex=="M"),points(FL,jitter(N,0.5),pch=19,col="blue"))

  #run models
Rand.eff="NO"  #note: for random effect, cannot predict SE!
MOD.FL=glm(N~FL*Sex, data=dat.FL.glm, family="binomial")
MOD.FL.no.sex=glm(N~FL, data=dat.FL.glm, family="binomial")

if(Rand.eff=="YES")
{
  MOD.FL=glmer(N~FL*Sex+(1 | TagCode), data=dat.FL.glm, family="binomial")
  MOD.FL.no.sex=glmer(N~FL+(1 | TagCode), data=dat.FL.glm, family="binomial")
}

Anova=anova(MOD.FL, test = "Chisq")
Anova=as.data.frame(Anova[1:5])
Anova$Terms=row.names(Anova)
Anova=Anova[-1,match(c("Terms","Df","Deviance","Pr(>Chi)"),names(Anova))]
Anova$Dev.exp.by.model=Dsquared(MOD.FL,adjust=T)*100
exp.tabl(Anova,"Anova.FL")

  #predict new data
SEQ_len=100
Rang=range(dat.FL.glm$FL)
Var=seq(Rang[1],Rang[2],length.out=SEQ_len)     
This.Fctrs=names(sort(table(dat.glm$TagCode)))  #get most common of factors
This.Fctrs=This.Fctrs[length(This.Fctrs)]
New.dat=expand.grid(FL=Var,TagCode=This.Fctrs,Sex="F")
New.dat$TagCode=factor(New.dat$TagCode,levels(dat.glm$TagCode))
New.dat$Sex=factor(New.dat$Sex,levels(dat.glm$Sex))
New.dat.male=New.dat
New.dat.male$Sex=factor("M",levels(dat.glm$Sex))

Pred_FL.female=predict(MOD.FL,newdata=New.dat,type='response',se.fit=T)
Pred_FL.male=predict(MOD.FL,newdata=New.dat.male,type='response',se.fit=T)
Pred_FL=predict(MOD.FL.no.sex,newdata=New.dat,type='response',se.fit=T)
New.dat.FL=New.dat

fem.col=rgb(.6,0,0.1,alpha=0.6)

Get_50_mig=data.frame(FL=New.dat.FL$FL,Prob_fem=Pred_FL.female$fit,Prob_fem_SE=Pred_FL.female$se.fit,
                      Prob_male=Pred_FL.male$fit,Prob_male_SE=Pred_FL.male$se.fit)


FL_F_50=fn.CI(Get_50_mig$Prob_fem,Get_50_mig$Prob_fem_SE,Get_50_mig$FL)
FL_M_50=fn.CI(Get_50_mig$Prob_male,Get_50_mig$Prob_male_SE,Get_50_mig$FL)
FL_50=rbind(FL_F_50,FL_M_50)

FL_50=data.frame(Sex=c("F","M"),FL_50)

write.csv(FL_50,"Predicted_FL_50%_migration.csv",row.names=F)



#  2. the timing of migration and the proportion of males and females undertaking migrations

dat.glm$S=with(dat.glm,ifelse(Zn.rec=="South",1,0))
dat.FL.timming=subset(dat.glm,!TagCode%in%never.north,select=c(TagCode,Sex,FL,N,S,week,Month)) 
dat.FL.timming$Month=factor(dat.FL.timming$Month)

  #see model data
par(mai=c(.45,1,.1,.1),oma=c(1,.1,.1,.1),mgp=c(.9,.4,0))
fn.see.month.glm.dat(subset(dat.FL.timming,Sex=="F"))
fn.see.month.glm.dat(subset(dat.FL.timming,Sex=="M"))

  #run model
MOdel=glm(S~Sex+Month, data=dat.FL.timming, family="binomial")

if(Rand.eff=="YES")
{
  MOdel=glmer(S~Sex+Month+(1 | TagCode), data=dat.FL.timming, family="binomial") 
}


Anova=anova(MOdel, test = "Chisq")
Anova=as.data.frame(Anova[1:5])
Anova$Terms=row.names(Anova)
Anova=Anova[-1,match(c("Terms","Df","Deviance","Pr(>Chi)"),names(Anova))]
Anova$Dev.exp.by.model=Dsquared(MOdel,adjust=T)*100
exp.tabl(Anova,"Anova.Month")



Mns=as.factor(1:12)
New.dat=expand.grid(Month=Mns,TagCode=This.Fctrs,Sex="F")
Pred.Month.fem=predict(MOdel,newdata=New.dat,type='response',se.fit=T)
New.dat=expand.grid(Month=Mns,TagCode=This.Fctrs,Sex="M")
Pred.Month.male=predict(MOdel,newdata=New.dat,type='response',se.fit=T)

Show.prob.north="YES"
if(Show.prob.north=="YES")
{
  Pred.Month.male$fit=1-Pred.Month.male$fit
  Pred.Month.fem$fit=1-Pred.Month.fem$fit
}

Max.prob.North.fem=max(Pred.Month.fem$fit)
Max.prob.North.male=max(Pred.Month.male$fit)

add.spline="NO"
add.spline.ontop="YES"

if(add.spline=="NO")tiff(file="Figure_prob_movement N_S.tiff",width = 1400, height = 2400,units = "px", res = 300,compression = "lzw")
if(add.spline=="YES")tiff(file="Figure_prob_movement N_S_spline.tiff",width = 1400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,1),mar=c(3,3,1,.1), oma=c(1,1.75,.1,1),las=1,mgp=c(1,0.8,0))

#Size effect 
show.comb.SEX="NO"
if(show.comb.SEX=="NO")
{
  
  #males
  plot(New.dat.FL$FL*100,Pred_FL.male$fit,type='l',lwd=2,col="blue",xlab="",
       ylab="",cex.axis=1.5)
  CI.pol(X=New.dat.FL$FL*100,Ylow=(Pred_FL.male$fit-1.96*Pred_FL.male$se.fit),
         Yhigh=(Pred_FL.male$fit+1.96*Pred_FL.male$se.fit),rgb(0,0,0.75,alpha=0.1),rgb(0,0,0.75,alpha=0.3))
  
  #females
  lines(New.dat.FL$FL*100,Pred_FL.female$fit,lwd=2,col=fem.col)
  CI.pol(X=New.dat.FL$FL*100,Ylow=(Pred_FL.female$fit-1.96*Pred_FL.female$se.fit),
         Yhigh=(Pred_FL.female$fit+1.96*Pred_FL.female$se.fit),rgb(.75,0,0,alpha=0.1),rgb(.75,0,0,alpha=0.3))
  legend("topleft",c("male","female"),lty=1,lwd=3,col=c("blue",fem.col),bty="n",cex=1.5)

}
if(show.comb.SEX=="YES")
{
  plot(New.dat.FL$FL*100,Pred_FL$fit,type='l',lwd=3,col="black",xlab="",
       ylab="",cex.axis=1.75)
  CI.pol(X=New.dat.FL$FL*100,Ylow=(Pred_FL$fit-1.96*Pred_FL$se.fit),
         Yhigh=(Pred_FL$fit+1.96*Pred_FL$se.fit),rgb(.5,.5,.5,alpha=0.2),rgb(.5,.5,.5,alpha=0.5))
}
mtext("Migration probability",2,las=3,cex=2,line=3.1)
mtext("Fork length (cm)",1,cex=2,line=2.35)


#Month effect
  #straight prediction
if(add.spline=="NO")
{
  plot(1:12,Pred.Month.male$fit,type='l',lwd=2,col="blue",xlab="",
       ylab="",cex.axis=1.5,ylim=c(0,1))
  CI.pol(X=1:12,Ylow=(Pred.Month.male$fit-1.96*Pred.Month.male$se.fit),
         Yhigh=(Pred.Month.male$fit+1.96*Pred.Month.male$se.fit),rgb(0,0,0.75,alpha=0.1),rgb(0,0,0.75,alpha=0.3))
  if(add.spline.ontop=="YES") fn.plt.pred.M.spline.only(Pred.Month.male,"black")
  
  #females
  lines(1:12,Pred.Month.fem$fit,lwd=2,col=fem.col)
  CI.pol(X=1:12,Ylow=(Pred.Month.fem$fit-1.96*Pred.Month.fem$se.fit),
         Yhigh=(Pred.Month.fem$fit+1.96*Pred.Month.fem$se.fit),rgb(.75,0,0,alpha=0.1),rgb(.75,0,0,alpha=0.3))
  if(add.spline.ontop=="YES") fn.plt.pred.M.spline.only(Pred.Month.fem,"red")
}
  #with smoothing spline
if(add.spline=="YES")
{
  plot(1:12,1:12,ylim=c(0,1),pch=19,cex=1.5,ylab="",xlab="",cex.axis=1.75,col="transparent")
  fn.plt.pred.M(Pred.Month.fem,rgb(.75,0,0,alpha=0.3),rgb(.75,0,0,alpha=0.1))
  fn.plt.pred.M(Pred.Month.male,rgb(0,0,0.75,alpha=0.3),rgb(0,0,0.75,alpha=0.1))
}
if(Show.prob.north=="NO") mtext("Prob. of ocurring south",2,las=3,cex=2,line=3)
if(Show.prob.north=="YES") mtext("Prob. of ocurring north",2,las=3,cex=2,line=3)
mtext("Month",1,cex=2,line=2.25)
dev.off()



