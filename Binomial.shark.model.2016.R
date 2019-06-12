 setwd("C:/Minor fisheries/Sharks/Dusky_migration")
library(beepr)
library(dplyr)
library(magrittr)
library(runjags)
dir()
# open up data
dat <- read.csv("rawrawdata1.csv")
head(dat)

## Female sharks
dat %<>% mutate(N=ifelse(Zn.rec=='North',1,0))
## Data.frame of release locations to use later
rel <- dat[!(duplicated(dat$TagCode)) & dat$Sex=='F',c('TagCode','Month.rel', 'FL', 'Zn.rel')] %>%
  mutate(N=ifelse(Zn.rel=='North',1,0)) %>% mutate(occurance='releas')
head(rel)

## monthy whether a shark was south or not
tdat <- dat  %>% 
  mutate(len10=floor(FL*10)) %>% 
  group_by(Sex,Month,len10) %>% 
  summarise(loc=max(N)) %>% 
  as.data.frame()
head(tdat)

## Examine data
par(mfrow=c(4,3))
for (i in 1:12){
  with(tdat[tdat$Sex=='F' & tdat$Month==i,], plot(len10, loc, pch=16, col=2, xlim=c(15,30), ylim=c(0,1.1), main=i))
  with(tdat[tdat$Sex=='M' & tdat$Month==i,], points(len10, loc+0.05, pch=16, col=4))
}

## Summarise into a monthly format and limit to girls
tdat <- dat  %>% 
  mutate(fmon=as.factor(Month))%>% 
  subset(Sex=='F') %>% 
  group_by(TagCode,fmon,FL) %>% 
  summarise(loc=max(N)) %>% 
  mutate(occurance='recap') %>%
  as.data.frame()

head(tdat)

head(rel)
## Tweak rel in to just girls
rel %<>% rename(fmon=Month.rel,loc=N) %>% select(-Zn.rel)

## Add the recapture and release data together
tdat <- rbind(tdat,rel)

## monthy whether a shark was south or not
tdat2 <- tdat  %>% 
  mutate(len10=floor(FL*10)) %>% 
  group_by(fmon,len10) %>% 
  summarise(loc=max(loc)) %>% 
  as.data.frame()
head(tdat2)

## Examine data
par(mfrow=c(4,3))
for (i in 1:12){
  with(tdat2[tdat2$fmon==i,], plot(len10, loc, pch=16, col=2, xlim=c(15,30), ylim=c(0,1.1), main=i))
}

## make the model matrix with cat as a factor
X <- model.matrix(~ fmon, data = tdat)                      
K <- ncol(X)  #Number of columns
head(X)

#  make data for random effects - one effect for each animal:
Tag <- as.numeric(as.factor(as.character(tdat$TagCode)))
Tagids <- sort(unique(tdat$TagCode))
Nre <- length(unique(tdat$TagCode))

## make model data input
JAGS.data <- list(Migrate = tdat$loc,
                  len     = tdat$FL*100,
                  X       = X,
                  #Tag     = Tag,
                  #Nre     = Nre,
                  N       = nrow(X),
                  K       = ncol(X))

###################################################
# JAGS model code
sink("JAGSmodRE.txt")
cat("
    model{
    #1A. Priors alpha 
    mx[1] ~ dunif(0, 1) 
    for (i in 2:K) { mx[i] ~ dunif((-mx[1]), (1-mx[1]))}

    #1B. Priors beta 
    inflec ~ dunif(100, 400) 
    #1C. Priors gamma 
    slope ~ dunif(-100, -1) 

    # Priors for random effects and Sigma_tag Looks confusing becasue I am using a Half-Cauchy distribution
    #num           ~ dnorm(0, 0.0016)          
    #denom         ~ dnorm(0, 1)               
    #Sigma_tag   <- abs(num / denom)
    #tau_Plot <- 1 / (Sigma_tag * Sigma_tag)
    #for (i in 1:Nre) { a[i] ~ dnorm(0, tau_Plot)} 
    
    # Likelihood
    for (i in 1:N) {
    mxV[i] <-     inprod(mx[], X[i,])  
    Pi[i] <- mxV[i] / (1+exp((len[i]-inflec)/slope)) #+ a[Tag[i]]
    Migrate[i] ~ dbern(Pi[i])  ## Likelihood
    }          
    }
    ",fill = TRUE)
sink()
#####################################

# Initial values 
inits  <- function () {
  list(
    mx =      c(0.8, rep(0, (ncol(X)-1))),
    inflec  = runif(1,200,250), #c(220, rep(0, (ncol(X)-1))),
    slope  = runif(1,-50,-5))} #c(-50, rep(0, (ncol(X)-1))))}#,
    #num   = rnorm(1, 0, 25), 
    #denom = rnorm(1, 0, 1))}

# Parameters to monitor
params <- c("mx","inflec","slope")#, "Sigma_tag")

## Run the jags model and save in r1 
r1 <- run.jags("JAGSmodRE.txt", monitor=params, data=JAGS.data,inits=inits, 
                      n.chains=3, method = "rjparallel", thin=3,
                      sample=200000)


## Examine mixing
#plot(r1, plot.type = "trace", col=1:3)
#plot(r1, plot.type = "histogram")
print(r1)


beep(7)

## Extend the model with heaps more runs
#r1 <- extend.jags(r1, sample=5000)
#print(r1)

# if you want to save the model
# dput(r1, 'Female run')
# if you want to re-load the model
# r1 <- (dget('Female run'))


## Convert to include offsets and then make summaries
cns <- r1$mcmc
cns <- as.array(cns)
cnames <- colnames(cns)
DIM <- dim(cns)
cns <- aperm(cns,c(1,3,2))  ## Change the order of the array
dim(cns)<- c(DIM[1]*DIM[3], DIM[2])       ## Change the dimensions to a matrix
colnames(cns) <- cnames
cns <- as.data.frame(cns)
head(cns)

##new method for tweaking output
out2 <- data.frame(tmp=rep(NA, nrow(cns)))
tmp <- expand.grid(fmon=unique(tdat$fmon))
tmp2 <- model.matrix(~fmon , data = tmp)
nms <- tdat$fmon
i <- 1
for(i in 1:length(unique(nms))){
  out2$tmp <- as.vector(tmp2[i,] %*% t(cns[,grepl('mx', colnames(cns), ignore.case = T)]))
  Name <- as.character(unique(nms)[i])
  Name <- gsub(' ','_',Name)
  Name <- paste('Mon',Name,sep='')
  colnames(out2)[grepl('tmp',colnames(out2))]  <- Name
}


head(out2)

(paradj <- (apply(out2, c(2), quantile, probs=c(0.005, 0.025, 0.125,0.5,0.875,0.975, 0.995))))

(par1 <- t(paradj))

par(mfrow=c(1,1), mar=c(7,7,5,5),las=1)
plot(1:nrow(par1), par1[,'50%'], xlim=c(0.5,(nrow(par1)+0.5)), type='o', ylim=c(0,1.1), pch=16, xlab='Month', ylab='Proportion migrating',axes=F, col=2)
polygon(c(1:nrow(par1), nrow(par1):1), c(par1[,"0.5%"], rev(par1[,"99.5%"])), col=rgb(1,0,0,0.2), border=F)
polygon(c(1:nrow(par1), nrow(par1):1), c(par1[,"2.5%"], rev(par1[,"97.5%"])), col=rgb(1,0,0,0.2), border=F)
lines(1:nrow(par1), par1[,'50%'],col=2)
axis(1,1:nrow(par1), month.abb) 
axis(2)
legend('top', bty='n', fill=c(rgb(1,0,0,0.2),rgb(1,0,0,0.4)), legend=c('95%', '99%'), ncol=2)

# rand <- cns[,grepl('a[', colnames(cns), fixed=T)]
# (rand <- t(apply(as.matrix(rand), c(2), quantile, probs=c(0.005, 0.025,0.5,0.975, 0.995))))
# 
# plot(1:nrow(rand), rand[,'50%'], ylim=c(-.02,0.04))
# arrows(1:nrow(rand), rand[,'2.5%'],y1=rand[,'97.5%'],code=3,angle=90,length=0.05)
## Record female results

## Function to pull curves out of mcmc chains
estout <- function(x){   tmp <- out2[,loc]/(1+exp((x-cns[,'inflec'])/cns[,'slope']))
                         return((c(mean(tmp), quantile(tmp,probs = c(0.025, 0.975)))  ) )   }

i <- 8
loc <- paste('Mon',i,sep="")
lens <- seq(100,300,10)
Fcurves <- t(apply(as.matrix(lens),1,estout))
i <- 1
loc <- paste('Mon',i,sep="")
lens <- seq(100,300,10)
Fcurves1 <- t(apply(as.matrix(lens),1,estout))

par(mfrow=c(1,1), mar=c(7,7,5,5),las=1)
plot(lens, Fcurves[,1], type='l', axes=F,ylim=c(0,1))
polygon(c(lens,rev(lens)), c(Fcurves[,2],rev(Fcurves[,3])), border=F, col=grey(0.2,0.2))
axis(1); axis(2,seq(0,1,0.2))


Fpar1 <- par1
Fpar1a <- print(r1)
#############################################################################################
###### Males
## Remove female info to increase memory
## rm(r1)
## Data.frame of release locations to use later
rel <- dat[!(duplicated(dat$TagCode)) & dat$Sex=='M',c('TagCode','Month.rel', 'FL', 'Zn.rel')] %>%
  mutate(N=ifelse(Zn.rel=='North',1,0)) %>% mutate(occurance='releas')
head(rel)

## Summarise into a monthly format and limit to girls
tdat <- dat  %>% 
  mutate(fmon=as.factor(Month))%>% 
  subset(Sex=='M') %>% 
  group_by(TagCode,fmon,FL) %>% 
  summarise(loc=max(N)) %>% 
  mutate(occurance='recap') %>%
  as.data.frame()

head(tdat)
head(rel)
## Tweak rel in to just girls
rel %<>% rename(fmon=Month.rel,loc=N) %>% select(-Zn.rel)

tdat <- rbind(tdat,rel)

## monthy whether a shark was south or not
tdat2 <- tdat  %>% 
  mutate(len10=floor(FL*10)) %>% 
  group_by(fmon,len10) %>% 
  summarise(loc=max(loc)) %>% 
  as.data.frame()
head(tdat2)

## Examine data
par(mfrow=c(4,3))
for (i in 1:12){
  with(tdat2[tdat2$fmon==i,], plot(len10, loc, pch=16, col=4, xlim=c(15,30), ylim=c(0,1.1), main=i))
}

tdat$fmon <- factor(tdat$fmon,as.character(c(5:12,1:4)))

## make the model matrix with cat as a factor
X <- model.matrix(~ fmon, data = tdat)                      
K <- ncol(X)  #Number of columns
head(X)
###################################################
# JAGS model code
sink("JAGSmodRE.txt")
cat("
    model{
    #1A. Priors alpha 
    mx[1] ~ dunif(0, 1) 
    for (i in 2:K) { mx[i] ~ dunif((-1-mx[1]), (1-mx[1]))}

    #1B. Priors beta 
    inflec ~ dunif(100, 400) 
    #1C. Priors gamma 
    slope ~ dunif(-100, -1) 

    # Priors for random effects and Sigma_tag Looks confusing becasue I am using a Half-Cauchy distribution
    #num           ~ dnorm(0, 0.0016)          
    #denom         ~ dnorm(0, 1)               
    #Sigma_tag   <- abs(num / denom)
    #tau_Plot <- 1 / (Sigma_tag * Sigma_tag)
    #for (i in 1:Nre) { a[i] ~ dnorm(0, tau_Plot)} 
    
    # Likelihood
    for (i in 1:N) {
    mxV[i] <-     inprod(mx[], X[i,])  
    Pi[i] <- mxV[i] / (1+exp((len[i]-inflec)/slope)) #+ a[Tag[i]]
    Migrate[i] ~ dbern(Pi[i])  ## Likelihood
    }          
    }
    ",fill = TRUE)
sink()
#####################################

#  make data for random effects - one effect for each animal:
Tag <- as.numeric(as.factor(as.character(tdat$TagCode)))
Tagids <- sort(unique(tdat$TagCode))
Nre <- length(unique(tdat$TagCode))

## make model data input
JAGS.data <- list(Migrate = tdat$loc,
                  len     = tdat$FL*100,
                  X       = X,
                  #Tag     = Tag,
                  #Nre     = Nre,
                  N       = nrow(X),
                  K       = ncol(X))

# Initial values 
inits  <- function () {
  list(
    mx =      c(1, rep(0, (ncol(X)-1))),
    inflec  = runif(1,200,250), #c(220, rep(0, (ncol(X)-1))),
    slope  = runif(1,-50,-5))} #c(-50, rep(0, (ncol(X)-1))))}#,
    #num   = rnorm(1, 0, 25), 
    #denom = rnorm(1, 0, 1))}

# Parameters to monitor
params <- c("mx","inflec","slope")#, "Sigma_tag")

## Run the jags model and save in r1 
r2 <- run.jags("JAGSmodRE.txt", monitor=params, data=JAGS.data,inits=inits, 
               n.chains=3, method = "rjparallel", thin=3,
               sample=200000)

## Examine mixing
#plot(r2)
#plot(r2, plot.type = "trace", col=1:3)
#plot(r2, plot.type = "histogram")
print(r2)
beep(7)

## Extend the model with heaps more runs
#r2 <- extend.jags(r2, sample=100000, add.monitor=c('Sigma_tag','a')  )

# if you want to save the model
# dput(r2, 'Binomial MCCM run months')
# if you want to re-load the model
# r2 <- (dget('Binomial MCCM run months'))

## Convert to include offsets and then make summaries
cns <- r2$mcmc
cns <- as.array(cns)
cnames <- colnames(cns)
DIM <- dim(cns)
cns <- aperm(cns,c(1,3,2))  ## Change the order of the array
dim(cns)<- c(DIM[1]*DIM[3], DIM[2])       ## Change the dimensions to a matrix
colnames(cns) <- cnames
cns <- as.data.frame(cns)
head(cns)

##new method for tweaking output
out2 <- data.frame(tmp=rep(NA, nrow(cns)))
tmp <- expand.grid(fmon=unique(tdat$fmon))
tmp2 <- model.matrix(~fmon , data = tmp)
nms <- tdat$fmon
i <- 1
for(i in 1:length(unique(nms))){
  out2$tmp <- as.vector(tmp2[i,] %*% t(cns[,grepl('mx', colnames(cns), ignore.case = T)]))
  Name <- as.character((unique(nms))[i])
  Name <- gsub(' ','_',Name)
  Name <- paste('Mon',Name,sep='')
  colnames(out2)[grepl('tmp',colnames(out2))]  <- Name
}

head(out2)
## Resort into month order
out2 <- out2[,paste('Mon',1:12,sep='')]

(paradj <- (apply(out2, c(2), quantile, probs=c(0.005,0.025,0.125,0.5,0.875,0.975, 0.995))))

(par1 <- t(paradj))

# par(mfrow=c(1,2), mar=c(7,4,5,0),las=1)
# plot(1:nrow(Fpar1), Fpar1[,'50%'], xlim=c(0.5,(nrow(Fpar1)+0.5)), type='o', ylim=c(0,1.1), pch=16, xlab='Month', ylab='Proportion migrating',axes=F, col=2)
# polygon(c(1:nrow(Fpar1), nrow(Fpar1):1), c(Fpar1[,"0.5%"], rev(Fpar1[,"99.5%"])), col=rgb(1,0,0,0.2), border=F)
# polygon(c(1:nrow(Fpar1), nrow(Fpar1):1), c(Fpar1[,"2.5%"], rev(Fpar1[,"97.5%"])), col=rgb(1,0,0,0.2), border=F)
# lines(1:nrow(Fpar1), Fpar1[,'50%'],col=2)
# axis(1,1:nrow(par1), month.abb) 
# axis(2)
# plot(1:nrow(par1), par1[,'50%'], xlim=c(0.5,(nrow(par1)+0.5)), type='o', ylim=c(0,1.1), pch=16, xlab='Month', ylab='Proportion migrating',axes=F, col=4)
# polygon(c(1:nrow(par1), nrow(par1):1), c(par1[,"0.5%"], rev(par1[,"99.5%"])), col=rgb(0,0,1,0.2), border=F)
# polygon(c(1:nrow(par1), nrow(par1):1), c(par1[,"2.5%"], rev(par1[,"97.5%"])), col=rgb(0,0,1,0.2), border=F)
# lines(1:nrow(par1), par1[,'50%'],col=4,type='o',pch=16)
# axis(1,1:nrow(par1), month.abb) 
# axis(2)
# legend('top', bty='n', fill=c(rgb(0,0,1,0.2),rgb(0,0,1,0.4)), legend=c('95%', '99%'), ncol=2)


tiff('Fig 4b.tif', width=3000, height=2500, res=300,compression = 'lzw')
par(mfrow=c(1,2), mar=c(7,4,5,0),las=1)
plot(1:nrow(Fpar1), Fpar1[,'50%'], xlim=c(0.5,(nrow(Fpar1)+0.5)), type='o', ylim=c(0,1.1), pch=16, xlab='Month', ylab='Proportion migrating',axes=F, col=1)
polygon(c(1:nrow(Fpar1), nrow(Fpar1):1), c(Fpar1[,"12.5%"], rev(Fpar1[,"87.5%"])), col=grey(0.2,0.2), border=F)
polygon(c(1:nrow(Fpar1), nrow(Fpar1):1), c(Fpar1[,"2.5%"], rev(Fpar1[,"97.5%"])), col=grey(0.2,0.2), border=F)
lines(1:nrow(Fpar1), Fpar1[,'50%'],col=1)
axis(1,1:nrow(par1), month.abb) 
axis(2)
mtext(side=3,at=1,text='Females',line=-3,font=2)
plot(1:nrow(par1), par1[,'50%'], xlim=c(0.5,(nrow(par1)+0.5)), type='o', ylim=c(0,1.1), pch=16, xlab='Month', ylab='Proportion migrating',axes=F, col=1)
polygon(c(1:nrow(par1), nrow(par1):1), c(par1[,"12.5%"], rev(par1[,"87.5%"])), col=grey(0.2,0.2), border=F)
polygon(c(1:nrow(par1), nrow(par1):1), c(par1[,"2.5%"], rev(par1[,"97.5%"])), col=grey(0.2,0.2), border=F)
lines(1:nrow(par1), par1[,'50%'],col=1,type='o',pch=16)
axis(1,1:nrow(par1), month.abb) 
axis(2)
mtext(side=3,at=1,text='Males',line=-3,font=2)
par(mfrow=c(1,1), mar=c(1,1,5,1),las=1, new=T)
plot(1,1,type='n',xlab='',ylab='',axes=F)
legend('top', bty='n', fill=c(grey(0.2,0.4),grey(0.2,0.2)), legend=c('75%', '95%'), ncol=2)
dev.off()

# rand <- cns[,grepl('a[', colnames(cns), fixed=T)]
# (rand <- t(apply(as.matrix(rand), c(2), quantile, probs=c(0.005, 0.025,0.5,0.975, 0.995))))
# 
# plot(1:nrow(rand), rand[,'50%'], ylim=c(-.02,0.04))
# arrows(1:nrow(rand), rand[,'2.5%'],y1=rand[,'97.5%'],code=3,angle=90,length=0.05)


i <- 8
loc <- paste('Mon',i,sep="")
lens <- seq(100,300,10)
Mcurves <- t(apply(as.matrix(lens),1,estout))
i <- 1
loc <- paste('Mon',i,sep="")
Mcurves1 <- t(apply(as.matrix(lens),1,estout))

tiff('Fig 4a.tif', width=3000, height=2500, res=300,compression = 'lzw')
##Females
par(mfrow=c(1,2), mar=c(7,4,5,0),las=1)
plot(lens, Fcurves[,1], type='l', axes=F,ylim=c(0,1),ylab='Proportion Migrating',xlab="Total length (cm)")
polygon(c(lens,rev(lens)), c(Fcurves[,2],rev(Fcurves[,3])), border=F, col=grey(0.2,0.2))
lines(lens, Fcurves1[,1], type='l', lty=3)
polygon(c(lens,rev(lens)), c(Fcurves1[,2],rev(Fcurves1[,3])), border=F, col=grey(0.7,0.2))
axis(1); axis(2,seq(0,1,0.2))
mtext(side=3,at=110,text='Females',line=-1,font=2)
##Males
plot(lens, Mcurves[,1], type='l', axes=F,ylim=c(0,1),ylab='Proportion Migrating',xlab="Total length (cm)")
polygon(c(lens,rev(lens)), c(Mcurves[,2],rev(Mcurves[,3])), border=F, col=grey(0.2,0.2))
lines(lens, Mcurves1[,1], type='l', lty=3)
polygon(c(lens,rev(lens)), c(Mcurves1[,2],rev(Mcurves1[,3])), border=F, col=grey(0.7,0.2))
axis(1); axis(2,seq(0,1,0.2))
mtext(side=3,at=110,text='Males',line=-1,font=2)
par(mfrow=c(1,1), mar=c(1,1,3,1),las=1, new=T)
plot(1,1,type='n',xlab='',ylab='',axes=F)
legend('top', bty='n', lty=c(3,1), legend=c('January', 'August'), ncol=2)
dev.off()