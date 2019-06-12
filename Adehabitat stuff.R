library(adehabitat)

#EXPLORE ALL THE FUNCTIONS IN THIS PACKAGE!!!


#Minimum convex polygon
data(puechabon)
locs <- puechabon$locs
cp <- mcp(locs[,4:5], locs[,1])
## Plot the home ranges
opar <- par(mar = c(0,0,0,0))
area.plot(cp)
## ... And the relocations
points(locs[,4:5], pch = 16, col = as.numeric(locs[,1]))
par(opar)
## Computation of the home-range size:
if (require(gpclib)) {
  cuicui1 <- mcp.area(locs[,4:5], locs[,1])
  plot(cuicui1)
}


MCP=function(DAT,SPEC)
{
  dat=subset(DAT,SPECIES%in%SPEC)
  ID.pos=match(c("LAT","LONG"),names(dat))
  ID.spec=match("SPECIES",names(dat))
  cp <- mcp(dat[,ID.pos], dat[,ID.spec])
  opar <- par(mar = c(0,0,0,0))
  area.plot(cp)
  points(dat[,ID.pos], pch = 16, col = as.numeric(dat[,ID.spec]))
  par(opar)
  
}
#MCP(Data.monthly,TARGETS)



data(bear)
## compute the sequence of dates at which the UD is to be
## estimated
vv <- seq(min(bear[[1]]$date), max(bear[[1]]$date), length=50)
head(vv)



#SOME KERNEL DENSITIES

## estimates the UD at each time point
re <- lapply(1:length(vv), function(i) {
  ## estimate the UD. We choose a smoothing parameter of
  ## 1000 meters for X and Y coordinates, and of 72 hours
  ## for the time (after a visual exploration)
  uu <- kernelkc(bear, h = c(1000,1000,72*3600),
                 tcalc= vv[i])
  ## now, we show the result
  ## potentially, we could type
  ##
  ## jpeg(paste("UD", i, ".jpg", sep=""))
  ##
  ## to store the figures in a file, and then to build a
  ## movie with the resulting files:
  ##
  image(uu[[1]], main=vv[i])
  ## highlight the 95 percent home range
  hh <- getvolumeUDk(uu)
  contour(hh[[1]], levels=95, col="red",
          drawlabels=FALSE, add=TRUE)
})



## Or, just show the home range:
re <- lapply(1:length(vv), function(i) {
  uu <- kernelkc(bear, h = c(1000,1000,72*3600),
                 tcalc= vv[i])
  pc <- getverticeshrk(uu, lev=95)
  plot(pc, xlim=c(510000, 530000),
       ylim=c(6810000, 6825000), main=vv[i])
})




data(puechabon)
loc <- puechabon$locs[, c("X", "Y")]
id <- puechabon$locs[, "Name"]
## Estimation of UD for the four animals
(ud <- kernelUD(loc, id))
image(ud) ## Note that the contours
## corresponds to values of probability density
udvol <- getvolumeUD(ud)
image(udvol)
## Here, the contour corresponds to the
## home ranges estimated at different probability
## levels (i.e. the contour 90 corresponds to the 90 percent
## kernel home-range)
## udvol describes, for each cell of the grid,
## the smaller home-range to which it belongs