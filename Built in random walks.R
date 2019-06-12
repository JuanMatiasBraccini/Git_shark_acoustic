
test=as.ltraj(Brownian.trajectory[[1]],typeII=F,id=as.factor("test"))
b=hbrown(test)

test1=as.ltraj(CRW.trajectory[[1]],typeII=F,id=as.factor("test"))
b1=hbrown(test1)

Nstep=52560
par(mfcol=c(2,1),mar=c(3,3,1,1))
plot(simm.brown(1:Nstep, x0 = c(109, 10),h=b),addpoints = FALSE,main="Brownian",xlim=c(108,111),
     ylim=c(9,11))
plot(simm.crw(1:Nstep, x0 = c(109, 10),h=b1),addpoints = FALSE,main="CRW",xlim=c(108,111),
     ylim=c(9,11))

Max=520
set.seed(411)
w <- simm.levy(1:Max, mu = 1.5, burst = "mu = 1.5")
u <- simm.levy(1:Max, mu = 2, burst = "mu = 2")
v <- simm.levy(1:Max, mu = 2.5, burst = "mu = 2.5")
x <- simm.levy(1:Max, mu = 3, burst = "mu = 3")
par(mfrow=c(2,2))
lapply(list(w,u,v,x), plot, perani=FALSE)