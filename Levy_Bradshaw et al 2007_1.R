###############################################
## Power-law tail probability density function
###############################################
#R code generates an idealized Lévy flight power-law tail probability density function based on µ = 2
# and = 1. The code also provides a sampling 'check' function to determine µ from the relationship 
# between the log of the step length frequencies and the log of the step lengths (corrected for the 
# appropriate binning procedure as discussed in Appendix A and the main text of the paper)


############################################################
## Author: Corey J. A. Bradshaw (corey.bradshaw@cdu.edu.au) 
############################################################

######################
## Date: 14 June 2006
######################

## Set mu and alpha values
mu <- 2
alpha <- 1

## Set maximum and minimum daily step lengths (km)
max.daily.len <- 5
min.daily.len <- 0.1

## Set number of iterations and create step length vector
iter <- 10000 ## number of iterations
st.len.vec <- seq(min.daily.len,max.daily.len,length.out=1000)

## Create probability vector
pr.lenfr.vec <- 10^(log10(alpha) - (mu*log10(st.len.vec)))
pr.lenfr.vec <- pr.lenfr.vec/sum(pr.lenfr.vec)

## Plot probability density function
plot(st.len.vec,pr.lenfr.vec,type="l",xlab="step length (km)",ylab="Pr")

## Check step length sampling function
## random step length from power-tail distribution (above)

len.check <- 0
  for (i in 1:10000) {
		len.check[i] <- sample(st.len.vec,1,replace=T,prob=pr.lenfr.vec)
	}

len.check <- len.check[len.check > min.daily.len]

## create breaks vector (set for 100 km maximum daily step length)
break.st <- 0
break.int <- 0.15

break.vec <- 2^(seq(break.st,2*log(max.daily.len),break.int))
break.max.sub <- which(break.vec >= (1.3*max.daily.len))
break.vec <- break.vec[-break.max.sub]
break.vec <- c(break.vec)

## Make plots
par(mfrow=c(1,2),pty="s")
hist.check <- (hist(len.check,br=break.vec,xlab="step length bin",main=""))
ch.len.fr <- hist.check$counts
ch.len.fr <- ifelse(ch.len.fr == 0,NA,ch.len.fr)
delta.k <- hist.check$breaks[2:length(hist.check$breaks)] - hist.check$breaks[1:(length(hist.check$breaks)-1)]
hist.check2 <- na.omit(data.frame(ch.len.fr,delta.k))
freq.delta.k <- hist.check2$ch.len.fr/hist.check2$delta.k
plot(log10(hist.check2$delta.k),log10(freq.delta.k),xlab="log10 X", ylab="log10 (N(X)/bin width)")
abline(lm(log10(freq.delta.k) ~ log10(hist.check2$delta.k)))
mu.pred <- as.numeric((lm(log10(freq.delta.k) ~ log10(hist.check2$delta.k)))$coefficients[2])
mu.pred.rnd <- round(mu.pred,3)
mtext(paste("m"," = ",mu.pred.rnd),side=3,line=-2,font=5,cex=1.2)
par(mfrow=c(1,1))

## End