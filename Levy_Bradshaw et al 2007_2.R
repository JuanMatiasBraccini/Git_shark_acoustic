###############################################
## Lévy flight track simulation based on the 
## power-law tail probability density function
###############################################

#uses the idealized Lévy flight power-law tail probability density function to simulate animal 
# movement tracks over 365 days assuming random turn angles after each step


############################################################
## Author: Corey J. A. Bradshaw (corey.bradshaw@cdu.edu.au) 
############################################################

######################
## Date: 14 June 2006
######################

## Set number of days tracked
days.tracked <- 365 ## set by user

mu.vec <- levy.check <- 0

## Set number of iterations
mu.iter <- 1000

## create breaks vector (for 100 km maximum step length)
break.st <- 0
break.int <- 0.001*days.tracked

for (i in 1:mu.iter) {

  x.prime <- y.prime <- x.vec <- y.vec <- len.vec <- 0

	for (d in 1:days.tracked) {
	
		x <- x.prime
		y <- y.prime

		## random step length from power-law tail distribution (above)
		len <- sample(st.len.vec,1,replace=T,prob=pr.lenfr.vec)
		len

		## Turn angle
		theta <- runif(1,0,359)
		theta.rad <- theta * (pi/180)

		## new coordinates
		x.prime <- x + (len * cos(theta.rad))
		y.prime <- y + (len * sin(theta.rad))

		x.vec[d] <- x.prime
		x.vec
		y.vec[d] <- y.prime
		y.vec

		if (len < min.daily.len) len <- min.daily.len
		if (len > max.daily.len) len <- max.daily.len
		len.vec[d] <- len
		len.vec
	}

	x.vec <- c(0,x.vec)
	y.vec <- c(0,y.vec)
	xy.dat <- data.frame(x.vec,y.vec)
	colnames(xy.dat) <- c("x","y")

	## create breaks vector
	break.vec <- 2^(seq(break.st,2*log(max(len.vec)),break.int))
	break.max.sub <- which(break.vec >= (1.3*(max(len.vec))))
	break.vec <- break.vec[-break.max.sub]
	if (break.vec[length(break.vec)] < max(len.vec)) break.vec <- c(break.vec,max(len.vec))
	
	par(mfrow=c(1,3),pty="s")
	plot(x.vec,y.vec,pch=4,cex=0.6,xlab="Easting (km)",ylab="Northing (km)",type="l",cex.axis=1.3,cex.lab=1.5)

	hist.fit <- (hist(len.vec,br=break.vec,xlab="step length bin",main="",cex.axis=1.3,cex.lab=1.5))
	len.fr <- hist.fit$counts
	len.fr <- ifelse(len.fr == 0,NA,len.fr)
	delta.k <- hist.fit$breaks[2:length(hist.fit$breaks)] - hist.fit$breaks[1:(length(hist.fit$breaks)-1)]
	hist.fit2 <- na.omit(data.frame(len.fr,delta.k))
	freq.delta.k <- hist.fit2$len.fr/hist.fit2$delta.k

	plot(log10(hist.fit2$delta.k),log10(freq.delta.k),xlab="log step length",ylab="log frequency",cex.axis=1.3,cex.lab=1.5)
	abline(lm(log10(freq.delta.k) ~ log10(hist.fit2$delta.k)))
	mu.pred <- as.numeric((lm(log10(freq.delta.k) ~ log10(hist.fit2$delta.k)))$coefficients[2])
	mu.pred
	mu.pred.disp <- abs(as.numeric(round(mu.pred,2)))
	mtext(paste("m"," = ",mu.pred.disp),side=3,line=-2,font=5,cex=1.2)

	par(mfrow=c(1,1))

	if (mu.pred <= -1 & mu.pred >= -3) levy.check[i] <- 1

	mu.vec[i] <- abs(mu.pred)

	print(i)

}

prop.levy <- sum(levy.check)/mu.iter
prop.levy

## Plot mu histogram over individuals simulated
hist(mu.vec,br=(mu.iter/10),xlab="mu",main=paste("mean mu = ",round(mean(mu.vec),2)))
mu.lo <- as.numeric(quantile(mu.vec,probs=0.025))
mu.hi <- as.numeric(quantile(mu.vec,probs=0.975))
abline(v=mu.lo,col="red")
abline(v=mu.hi,col="red")
abline(v=mean(mu.vec),type="solid",col="red")
 
## End