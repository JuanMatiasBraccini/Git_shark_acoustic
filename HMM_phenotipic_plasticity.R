#Hidden Markov Model for modelling shark movement
library(moveHMM)

#--Data set up
# data.frame(ID=, Longitude (or Easting)=,Latitude (or Northing)=,other covariates)


#Simple example
head(elk_data)
#The easting and northing values are expressed in meters in the data  for the step length
#The last variable is the distance of the animal to water, 
# which for illustration purposes we want to include in the model as a covariate


#transform the coordinates into km
elk_data$Easting <- elk_data$Easting/1000
elk_data$Northing <- elk_data$Northing/1000

#Compute step length and turning angles
data <- prepData(elk_data,type="UTM",coordNames=c("Easting","Northing"))
#type specifies whether the coordinates are easting/northing (type="UTM") or longitude/latitude (type="LL") values
head(data)

#missing values are imputed using the closest non-missing value, by default the previous one if it is available

summary(data)

#plot maps of animals' tracks
plot(data,compact=T)
# The resulting map, and the steps and angles graphs for each animal are displayed
#The time series of step lengths is one way to check the data for outliers.


#--Fitting the HMM model
# model arguments:
#. nbStates=2, i.e. we fit a 2-state HMM to the data
#               state 1 involving relatively short steps and many turnings (hence the choice of a small initial value for the mean of the gamma step length distribution
#                 and an initial value of pi for the mean turning angle) 
#               state 2 involving longer steps and fewer turnings (hence the choice of a larger initial value for the mean of the gamma step length distribution
#                 and an initial value of 0 for the mean turning angle)


# . beta0=NULL and delta0=NULL, i.e. we use the default values for the initial values beta0 and delta0;
# . formula=???dist water, i.e. the transition probabilities are functions of the covariate "dist water";
# . stepDist="gamma", to model the step lengths with the gamma distribution (note that it is the
#                   default, so we do not need to explicitely specify it);
# . angleDist="vm", to model the turning angles with the von Mises distribution (default);
# . angleMean=NULL, because we want to estimate the mean of the angle distribution (default);
# . stationary=FALSE, as due to the covariates the process is not stationary (default).

#argument knownStates  makes it possible to set some values of the state process to fixed values, prior to fitting the model

# also need to specify initial values for the parameters of the state-dependent distributions;
# the algorithm might not find the global optimum of the likelihood function if the initial parameters are poorly
# chosen. The initial parameters should be specified in two vectors, stepPar0 (for the step distribution)
# and anglePar0 (for the angle distribution).
# to test different sets of initial values, possibly
# chosen randomly. By comparing the resulting estimates for the different initial values used, one usually
# obtains a good feeling for any potential sensitivity of the numerical search to its chosen starting point.


#Zero-inflation (as described in Section 4.1.2) must be included in the step length distribution if
#some steps are of length exactly zero (which is the case for the elk data). To do so, another parameter
# is added to the step distribution: its mass on zero.


#standardize the covariate values before fitting the model for numerical stability
data$dist_water <-(data$dist_water-mean(data$dist_water))/sd(data$dist_water)


## initial parameters for gamma and von Mises distributions
mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
sigma0 <- c(0.1,1) # step SD
zeromass0 <- c(0.1,0.05) # step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)
## call to fitting function
system.time({
  m <- fitHMM(data=data,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,formula=~dist_water)
            })

#see max like estimates
m

#Get confidence intervals
CIs=CI(m)
CIs$stepPar$upper


#plot model
plot(m)
# This outputs:
# . an histogram of step lengths of all animals, with the fitted state-dependent densities,
# . an histogram of turning angles of all animals, with the fitted state-dependent densities,
# . plots of the transition probabilities as functions of the covariate considered,
# . a map of each animal's track, colored by states.


#--Deciphering the model

#a) Viterbi algorithm 
# To globally decode the state process, the Viterbi algorithm is implemented in
#  the function viterbi. This function outputs the most likely sequence of states to have generated the
#  observation, under the fitted model. 
# Below are the most probable states for the first 25 observations of the first individual
states <- viterbi(m)
states[1:25]

#To get more accurate information on the state process, it is possible to compute
# the state probabilities for each observation
# This returns a matrix with as many columns as there are states in the model, and as many rows as 
# there are observations

#The elements of the matrix are defined as stateProbs(m)[t,j] = Pr(St = j)
sp <- stateProbs(m)
head(sp)

#The state with highest probability according to stateProbs might not be the same as the state in
# the most probable sequence returned by the Viterbi algorithm. This is because the Viterbi algorithm
# performs "global decoding", whereas the state probabilities are "local decoding"

# visualize the results of viterbi and stateProbs
plotStates(m,animals="elk-115")



#--Model selection using AIC to determine number of states
#  now fit a 3-state HMM to the data, and want to compare the AICs of the 2-state and 3-state models.

# initial parameters
mu0 <- c(0.1,0.5,3)
sigma0 <- c(0.05,0.5,1)
zeromass0 <- c(0.05,0.0001,0.0001)
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(pi,pi,0)
kappa0 <- c(1,1,1)
anglePar0 <- c(angleMean0,kappa0)

# fit the 3-state model
m3 <- fitHMM(data=data,nbStates=3,stepPar0=stepPar0,
             anglePar0=anglePar0,formula=~dist_water)
#Cmpare them:
  AIC(m,m3)
#In terms of AIC, the 3-state model is favoured over the 2-state model in this example

  
#--Model diagnostics
  
# pseudo-residuals (a.k.a. quantile residuals)
  pr <- pseudoRes(m)
  plotPR(m)
  
  
#-- plot tracking data on a satellite image
library(rgdal)
utmcoord <- SpatialPoints(cbind(data$x*1000,data$y*1000),proj4string=CRS("+proj=utm +zone=17"))
llcoord <- spTransform(utmcoord,CRS("+proj=longlat"))
lldata <- data.frame(ID=data$ID,x=attr(llcoord,"coords")[,1],y=attr(llcoord,"coords")[,2])
# In the code above, we need to multiply the UTM coordinates by 1000, as we had divided them
# by 1000 earlier to work with distances in kilometres. In the function SpatialPoints, we indicate
# +zone=17, because the data come from the UTM zone 17.
plotSat(lldata,zoom=8)
  
  
  