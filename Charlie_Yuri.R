#--SCRIPT FOR ANALYSING CROSS-JURISDICTIONAL MOVEMENTS OF DUSKY & BRONZE WHALERS--

library(ozmaps)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(circlize)
library(cmocean)

library(lubridate)
library(geosphere)
library(chron) 
library(mgcv)
library(mgcViz)
library(vioplot)
library(beanplot)
library(PBSmapping)
library(data.table)
options(stringsAsFactors = FALSE,dplyr.summarise.inform = FALSE)


# Data section ------------------------------------------------------------

# Load datasets and sort variables:
User='Matias'
#User='Yuri'

if(User=='Matias')
{
  handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')
  hndl.in=function(x) paste(handl_OneDrive("Analyses/Acoustic_tagging/For Charlie/Data"),x,sep='/')
  hndl.out=function(x) paste(handl_OneDrive("Analyses/Acoustic_tagging/For Charlie/Results"),x,sep='/')
}
  
if(User=='Yuri')
{
  hndl.in= function(x) paste("Final Data",x,sep='/')
  hndl.out= function(x) paste("Output",x,sep='/')
}
  
tag <- read.csv(hndl.in("Tagging_tot.csv"))
tag$ReleaseDate2 <- as.Date(tag$ReleaseDate2, format = "%Y-%m-%d")

df <- fread(hndl.in('Total_detections_2020_12_14.csv'),data.table=FALSE) 
df$Timestamp.UTC <- as.POSIXct(df$Timestamp.UTC, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")


# Yuri's code -------------------------------------------------------------


# Analyse ammounts of station/specific detections: IMOS x Non-IMOS
df$Station2 <- paste(df$Latitude, df$Longitude, sep = "_")
recs <- unique(df$Station2)
receiver <- NULL
n <- NULL
DATA <- NULL
Lat <- NULL
Lon <- NULL
for (i in 1:length(recs)) {
	aux <- subset(df, Station2 == recs[i])
	aux1 <- subset(aux, Species == "Carcharhinus brachyurus") # Bronze whalers
	aux2 <- subset(aux, Species == "Carcharhinus obscurus")   # Dusky whalers
	n <- c(n, nrow(aux1), nrow(aux2))

	if (nrow(aux1) == 0) {
		DATA <- c(DATA, aux2$Data[1], aux2$Data[1])
		Lat <- c(Lat, aux2$Latitude[1], aux2$Latitude[1])
		Lon <- c(Lon, aux2$Longitude[1], aux2$Longitude[1])
		receiver <- c(receiver, rep(paste0("Receiver.", i), 2))
	}

	if (nrow(aux2) == 0) {
		DATA <- c(DATA, aux1$Data[1], aux1$Data[1])
		Lat <- c(Lat, aux1$Latitude[1], aux1$Latitude[1])
		Lon <- c(Lon, aux1$Longitude[1], aux1$Longitude[1])
		receiver <- c(receiver, rep(paste0("Receiver.", i), 2))
	}

	if (nrow(aux1) > 0 & nrow(aux2) > 0) {
		DATA <- c(DATA, aux1$Data[1], aux2$Data[1])
		Lat <- c(Lat, aux1$Latitude[1], aux2$Latitude[1])
		Lon <- c(Lon, aux1$Longitude[1], aux2$Longitude[1])
		receiver <- c(receiver, rep(paste0("Receiver.", i), 2))
	}	
}
df.detections <- data.frame(Receiver = receiver, Species = c("Bronze", "Dusky"), Data = DATA, Latitude = Lat, Longitude = Lon, Detections = n)
rm(aux, aux1, aux2, i, DATA, n, Lat, Lon, receiver, recs)


# Plot species-specific maps:
oz_sf_states <- ozmap_data("states")  

plot1 <- ggplot() + theme_minimal() + 
	geom_sf(data = oz_sf_states, color = "darkgray", fill = "lightgray", alpha = 0.5, size = 0.2) +
	geom_point(data = subset(df.detections, Species == "Bronze" & Detections > 0 & Data == "IMOS"), aes(x = Longitude, y = Latitude, size = Detections), alpha = 0.5) +
	xlim(c(113, 155)) + labs(x = "Longitude", y = "Latitude", title = paste0("Bronze whalers (IMOS = ", nrow(subset(df, Species == "Carcharhinus brachyurus" & Data == "IMOS")), ")"))

plot2 <- ggplot() + theme_minimal() + 
	geom_sf(data = oz_sf_states, color = "darkgray", fill = "lightgray", alpha = 0.5, size = 0.2) +
	geom_point(data = subset(df.detections, Species == "Bronze" & Detections > 0 & Data == "Non-IMOS"), aes(x = Longitude, y = Latitude, size = Detections), alpha = 0.5) +
	xlim(c(113, 155)) + labs(x = "Longitude", y = "Latitude", title = paste0("Bronze whalers (Non-IMOS = ", nrow(subset(df, Species == "Carcharhinus brachyurus" & Data == "Non-IMOS")), ")"))

plot3 <- ggplot() + theme_minimal() + 
	geom_sf(data = oz_sf_states, color = "darkgray", fill = "lightgray", alpha = 0.5, size = 0.2) +
	geom_point(data = subset(df.detections, Species == "Dusky" & Detections > 0 & Data == "IMOS"), aes(x = Longitude, y = Latitude, size = Detections), alpha = 0.5) +
	xlim(c(113, 155)) + labs(x = "Longitude", y = "Latitude", title = paste0("Dusky whalers (IMOS = ", nrow(subset(df, Species == "Carcharhinus obscurus" & Data == "IMOS")), ")"))

plot4 <- ggplot() + theme_minimal() + 
	geom_sf(data = oz_sf_states, color = "darkgray", fill = "lightgray", alpha = 0.5, size = 0.2) +
	geom_point(data = subset(df.detections, Species == "Dusky" & Detections > 0 & Data == "Non-IMOS"), aes(x = Longitude, y = Latitude, size = Detections), alpha = 0.5) +
	xlim(c(113, 155)) + labs(x = "Longitude", y = "Latitude", title = paste0("Dusky whalers (Non-IMOS = ", nrow(subset(df, Species == "Carcharhinus obscurus" & Data == "Non-IMOS")), ")"))

ggarrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
ggsave(hndl.out("Map_detections_2020-12-14.png"), width = 20, height = 20, units = "cm")
rm(plot1, plot2, plot3, plot4, df.detections)


#====================================================#
# Randomized analyses of inter-state shark movements #
#====================================================#
#'
#' - Description:
#' Run a total of 1000 simulations per species using 50% of all tagged dusky/bronze whaler sharks 
#' to investigate number of inter-state movements using I) IMOS, II) Non-IMOS, and III) both receiver types!

# Check that state listing is correct:
#note: takes too long to plot all points
do.this=FALSE
if(do.this)
{
  ggplot() + theme_minimal() + 
    geom_sf(data = oz_sf_states, color = "darkgray", fill = "lightgray", alpha = 0.5, size = 0.2) +
    geom_point(data = df, aes(x = Longitude, y = Latitude, colour = State), alpha = 0.5)
  ggsave(hndl.out("Map1.png"), width = 20, height = 20, units = "cm")
}

# Run randomized network analysis: jump to 'Load saved simulation results' to load simulation results!
do.random=FALSE
# Empty datasets to store the data: 
bronze.imos <- NULL
bronze.nonimos <- NULL
bronze.total <- NULL
dusky.imos <- NULL
dusky.nonimos <- NULL
dusky.total <- NULL

if(do.random)
{
  for (simu in 1:1000) {
    cat(paste("Running random simulation number:", simu), fill = 1)
    
    ## SIMULATION STARTS
    
    # Select 50% of random sharks from each spp:
    aux.bronze <- sample(tag$Code2[tag$Species2 == "Carcharhinus brachyurus"], # only this species
                         length(tag$Code2[tag$Species2 == "Carcharhinus brachyurus"]) / 2) # use 50% of tagged sharks
    aux.dusky <- sample(tag$Code2[tag$Species2 == "Carcharhinus obscurus"], # only this species
                        length(tag$Code2[tag$Species2 == "Carcharhinus obscurus"]) / 2) # use 50% of tagged sharks
    
    # Subset total detection data and check which sharks were actually detected:
    df.bronze <- subset(df, TagCode %in% aux.bronze)
    detect.bronze <- unique(df.bronze$TagCode)
    df.dusky <- subset(df, TagCode %in% aux.dusky)
    detect.dusky <- unique(df.dusky$TagCode)
    
    # Subset datasets per each receiver type:
    bronze1 <- subset(df.bronze, Data == "IMOS")
    bronze2 <- subset(df.bronze, Data == "Non-IMOS")
    dusky1 <- subset(df.dusky, Data == "IMOS")
    dusky2 <- subset(df.dusky, Data == "Non-IMOS")
    
    # Analyze how many sharks moved inter-state:
    
    # Bronze whalers:
    bronze.save1 <- NULL # IMOS
    bronze.save2 <- NULL # Non-IMOS
    bronze.save3 <- NULL # TOTAL
    for (i in 1:length(detect.bronze)) {
      aux1 <- subset(bronze1, TagCode == detect.bronze[i])
      aux2 <- subset(bronze2, TagCode == detect.bronze[i])
      aux3 <- subset(df.bronze, TagCode == detect.bronze[i])
      
      if (nrow(aux1) == 0){
        bronze.save1 <- c(bronze.save1, "FALSE")
      } else {
        tagged <- tag$Tagging.State[tag$Code2 == detect.bronze[i]]
        detected <- unique(aux1$State)
        if (length(detected[-which(detected == tagged)]) > 0){
          bronze.save1 <- c(bronze.save1, "TRUE")
        } else {
          bronze.save1 <- c(bronze.save1, "FALSE")
        }
      }
      
      if (nrow(aux2) == 0){
        bronze.save2 <- c(bronze.save2, "FALSE")
      } else {
        tagged <- tag$Tagging.State[tag$Code2 == detect.bronze[i]]
        detected <- unique(aux2$State)
        if (length(detected[-which(detected == tagged)]) > 0){
          bronze.save2 <- c(bronze.save2, "TRUE")
        } else {
          bronze.save2 <- c(bronze.save2, "FALSE")
        }
      }
      
      if (nrow(aux3) == 0){
        bronze.save3 <- c(bronze.save3, "FALSE")
      } else {
        tagged <- tag$Tagging.State[tag$Code2 == detect.bronze[i]]
        detected <- unique(aux3$State)
        if (length(detected[-which(detected == tagged)]) > 0){
          bronze.save3 <- c(bronze.save3, "TRUE")
        } else {
          bronze.save3 <- c(bronze.save3, "FALSE")
        }
      }
    }
    bronze.imos <- c(bronze.imos, round((length(which(bronze.save1 == "TRUE")) / length(detect.bronze) * 100), digits = 1))
    bronze.nonimos <- c(bronze.nonimos, round((length(which(bronze.save2 == "TRUE")) / length(detect.bronze) * 100), digits = 1))
    bronze.total <- c(bronze.total, round((length(which(bronze.save3 == "TRUE")) / length(detect.bronze) * 100), digits = 1))
    
    
    # Dusky whalers:
    dusky.save1 <- NULL # IMOS
    dusky.save2 <- NULL # Non-IMOS
    dusky.save3 <- NULL # TOTAL
    for (i in 1:length(detect.dusky)) {
      aux1 <- subset(dusky1, TagCode == detect.dusky[i])
      aux2 <- subset(dusky2, TagCode == detect.dusky[i])
      aux3 <- subset(df.dusky, TagCode == detect.dusky[i])
      
      if (nrow(aux1) == 0){
        dusky.save1 <- c(dusky.save1, "FALSE")
      } else {
        tagged <- tag$Tagging.State[tag$Code2 == detect.dusky[i]]
        detected <- unique(aux1$State)
        if (length(detected[-which(detected == tagged)]) > 0){
          dusky.save1 <- c(dusky.save1, "TRUE")
        } else {
          dusky.save1 <- c(dusky.save1, "FALSE")
        }
      }
      
      if (nrow(aux2) == 0){
        dusky.save2 <- c(dusky.save2, "FALSE")
      } else {
        tagged <- tag$Tagging.State[tag$Code2 == detect.dusky[i]]
        detected <- unique(aux2$State)
        if (length(detected[-which(detected == tagged)]) > 0){
          dusky.save2 <- c(dusky.save2, "TRUE")
        } else {
          dusky.save2 <- c(dusky.save2, "FALSE")
        }
      }
      
      if (nrow(aux3) == 0){
        dusky.save3 <- c(dusky.save3, "FALSE")
      } else {
        tagged <- tag$Tagging.State[tag$Code2 == detect.dusky[i]]
        detected <- unique(aux3$State)
        if (length(detected[-which(detected == tagged)]) > 0){
          dusky.save3 <- c(dusky.save3, "TRUE")
        } else {
          dusky.save3 <- c(dusky.save3, "FALSE")
        }
      }
    }
    dusky.imos <- c(dusky.imos, round((length(which(dusky.save1 == "TRUE")) / length(detect.dusky) * 100), digits = 1))
    dusky.nonimos <- c(dusky.nonimos, round((length(which(dusky.save2 == "TRUE")) / length(detect.dusky) * 100), digits = 1))
    dusky.total <- c(dusky.total, round((length(which(dusky.save3 == "TRUE")) / length(detect.dusky) * 100), digits = 1))
  }
  
  df.simu.bronze <- data.frame(Simulation = c(1:1000), Freq.movements = c(bronze.imos, bronze.nonimos, bronze.total),
                               Type = c(rep("IMOS", 1000), rep("Non-IMOS", 1000), rep("Total", 1000)))
  df.simu.dusky <- data.frame(Simulation = c(1:1000), Freq.movements = c(dusky.imos, dusky.nonimos, dusky.total),
                              Type = c(rep("IMOS", 1000), rep("Non-IMOS", 1000), rep("Total", 1000)))
  
  write.csv(df.simu.bronze, hndl.out("Simu_bronze.csv"), row.names = FALSE)
  write.csv(df.simu.dusky, hndl.out("Simu_dusky.csv"), row.names = FALSE)
  
}

# Load saved simulation results: 2020-12-14
if(!do.random)
{
  df.simu.bronze <- read.csv("Final data/Simu_bronze.csv")
  df.simu.dusky <- read.csv("Final data/Simu_dusky.csv")
}

mod1 <- glm((Freq.movements/100) ~ Type, family = binomial, data = df.simu.bronze)
summary(mod1)
mod2 <- glm((Freq.movements/100) ~ Type, family = binomial, data = df.simu.dusky)
summary(mod2)

plot1 <- ggplot() + labs(title = "Bronze whalers", y = "Frequency of inter-state movements (%)", x = "Receiver type") +
	geom_boxplot(data = df.simu.bronze, aes(x = Type, y = Freq.movements, fill = Type)) + theme_bw() +
	scale_fill_manual(values = c('white', 'white', 'darkorange')) +
	theme(legend.position = "none")
plot2 <- ggplot() + labs(title = "Dusky whalers", y = "", x = "Receiver type") +
	geom_boxplot(data = df.simu.dusky, aes(x = Type, y = Freq.movements, fill = Type)) + theme_bw() +
	scale_fill_manual(values = c('white', 'darkorange', 'white')) +
	theme(legend.position = "none")

ggarrange(plot1, plot2, ncol = 2)
ggsave("Output/Randomized_simu.png", width = 25, height = 10, units = "cm")

rm(aux.bronze, aux.dusky, df.bronze, df.dusky, detect.bronze, detect.dusky,  # Clear working environment
	bronze1, bronze2, dusky1, dusky2, bronze.imos, bronze.nonimos, bronze.total,
	dusky.imos, dusky.nonimos, dusky.total, i, aux1, aux2, aux3,
	bronze.save1, bronze.save2, bronze.save3,
	tagged, detected, dusky.save1, dusky.save2, dusky.save3, simu,
	plot1, plot2, mod1, mod2)



#===========================================#
# Network analyses of inter-state movements #
#===========================================#

#================#
# Bronze Whalers #
#================#
df.bronze <- subset(df, Species == "Carcharhinus brachyurus")

# Get State connectivity data
trans <- unique(df.bronze$TagCode)
df.connect <- NULL
for (i in 1:length(trans)) {
	aux <- subset(df.bronze, TagCode == trans[i])
	tag.loc <- tag$Tagging.State[tag$Code2 == trans[i]]

	if (nrow(aux) > 1) {
		aux <- aux[order(aux$Timestamp.UTC), ]
		state1 <- NULL
		state2 <- NULL
		for(ii in 2:nrow(aux)) {
			if (ii == 1) {
				if (tag.loc != aux$State[ii - 1]) {
				state1 <- c(state1, tag.loc)
				state2 <- c(state2, aux$State[ii - 1])
				}
			}
			aux.first <- aux$State[ii - 1]
			aux.second <- aux$State[ii]
			if (aux.first != aux.second) {
				state1 <- c(state1, aux.first)
				state2 <- c(state2, aux.second)
			}
		}
		df.save <- data.frame(State1 = state1, State2 = state2, Transmitter = rep(trans[i], length(state1)))
		if (nrow(df.save) > 0)
			df.connect <- rbind(df.connect, df.save)
	}
}

# Identify number of connections = outcoming(+) incoming (-)
factors <- c("WA", "SA", "VIC", "NSW")
df.comb <- expand.grid(factors, factors)
df.comb$N <- NA
for (i in 1:nrow(df.comb)) {
	df.comb$N[i] <- length(which(df.connect$State1 == df.comb$Var1[i] & df.connect$State2 == df.comb$Var2[i]))
}
df.comb <- df.comb[-which(df.comb$N == 0), ]

# Identify sector widths: number of connections
widths.pos <- NULL
widths.neg <- NULL
for (i in 1:length(factors)) {
	aux1 <- sum(df.comb$N[df.comb$Var1 == factors[i]])
	aux2 <- sum(df.comb$N[df.comb$Var2 == factors[i]])

	widths.pos <- c(widths.pos, aux1)
	widths.neg <- c(widths.neg, aux2)
}
df.width <- data.frame(Factor = factors, Positive = widths.pos, Negative = widths.neg)

# Transform connections to proportions:
tot <- sum(df.comb$N)
df.width$Tot <- df.width$Positive + df.width$Negative
df.width$Width <- NA
for (i in 1:nrow(df.width)) {
	df.width$Width[i] <- df.width$Tot[i] / tot
}
df.width$Zero <- NA
for (i in 1:nrow(df.width)) {
	df.width$Zero[i] <- df.width$Positive[i] / tot
	if (df.width$Positive[i] > df.width$Negative[i]) 
		df.width$Zero[i] <- df.width$Zero[i] * -1
}

# Identify sector proportion for plotting the robbons:
df.width$Positive.length <- 1 - df.width$Zero
df.width$Negative.length <- 2 - df.width$Positive.length
df.width$Positive.prop <- df.width$Positive.length / df.width$Positive
df.width$Negative.prop <- (-1 - df.width$Zero) / df.width$Negative

# Create auxiliary point coordinates for robbon:
df.comb$Zero1 <- NA
df.comb$Point1 <- NA
df.comb$Zero2 <- NA
df.comb$Point2 <- NA
df.save <- NULL
out.fact <- unique(df.comb$Var1)
for (i in 1:length(out.fact)) {
	aux <- subset(df.comb, Var1 == out.fact[i])
	zero1 <- df.width$Zero[df.width$Factor == out.fact[i]]
	point1 <- zero1 + (df.width$Positive.prop[df.width$Factor == out.fact[i]] * aux$N[1])
	aux$Zero1[1] <- zero1
	aux$Point1[1] <- point1
	if (nrow(aux) > 1) {
		for (ii in 2:nrow(aux)) {
			zero1 <- point1
			point1 <- point1 + (df.width$Positive.prop[df.width$Factor == out.fact[i]] * aux$N[ii])
			aux$Zero1[ii] <- zero1
			aux$Point1[ii] <- point1
		}
	}
	df.save <- rbind(df.save, aux)
}
df.save2 <- NULL
inc.fact <- unique(df.save$Var2)
for (i in 1:length(inc.fact)) {
	aux <- subset(df.save, Var2 == inc.fact[i])
	zero2 <- df.width$Zero[df.width$Factor == inc.fact[i]]
	point2 <- zero2 + (df.width$Negative.prop[df.width$Factor == inc.fact[i]] * aux$N[1])
	aux$Zero2[1] <- zero2
	aux$Point2[1] <- point2
	if (nrow(aux) > 1) {
		for (ii in 2:nrow(aux)) {
			zero2 <- point2
			point2 <- point2 + (df.width$Negative.prop[df.width$Factor == inc.fact[i]] * aux$N[ii])
			aux$Zero2[ii] <- zero2
			aux$Point2[ii] <- point2
		}
	}
	df.save2 <- rbind(df.save2, aux)
}
df.comb <- df.save2


# Auxiliar object for colors
df.col <- data.frame(Factor = c("WA", "SA", "VIC", "NSW"), Colour = rainbow(4))
df.col$Colour <- paste0(df.col$Colour, 80)


### CIRCULAR PLOT
factors <- factor(c("WA", "SA", "VIC", "NSW"), levels = c("WA", "SA", "VIC", "NSW")) # Initialize circular plot layout

jpeg(hndl.out("Network_bronze.png"), width = 13, height = 13, units = "cm", quality = 500, bg = NA, res = 500)
par(mar = c(2, 2, 2, 2))
circos.par(canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.2, 1.2))
circos.initialize(factors, xlim = c(-1, 1), sector.width = df.width$Width)
# Plot tracks: Australian States
circos.track(ylim = c(0, 1), track.height = 0.05, bg.col = rainbow(4), bg.border = NA)
# Add axis
circos.axis(h = "top", sector.index = "WA", major.at = c(df.width$Zero[df.width$Factor == "WA"], 1), minor.ticks = 0, 
	labels = c(0, df.width$Positive[df.width$Factor == "WA"]), labels.cex = 0.8)
circos.text(x = 0, y = 10, labels = "WA", sector.index = "WA", niceFacing = TRUE)
circos.axis(h = "top", sector.index = "SA", major.at = c(df.width$Zero[df.width$Factor == "SA"], 1), minor.ticks = 0, 
	labels = c(0, df.width$Positive[df.width$Factor == "SA"]), labels.cex = 0.8)
circos.text(x = 0, y = 10, labels = "SA", sector.index = "SA", niceFacing = TRUE)
circos.axis(h = "top", sector.index = "VIC", major.at = c(df.width$Zero[df.width$Factor == "VIC"], 1), minor.ticks = 0, 
	labels = c(0, 1), labels.cex = 0.8)
circos.text(x = 0, y = 10, labels = "VIC", sector.index = "VIC", niceFacing = TRUE)
circos.axis(h = "top", sector.index = "NSW", major.at = c(df.width$Zero[df.width$Factor == "NSW"], 1), minor.ticks = 0, 
	labels = c(0, df.width$Positive[df.width$Factor == "NSW"]), labels.cex = 0.8)
circos.text(x = 0, y = 10, labels = "NSW", sector.index = "NSW", niceFacing = TRUE)

# Add robbons
out.fact <- unique(df.comb$Var1)
for (i in 1:length(out.fact)) {
	aux <- subset(df.comb, Var1 == out.fact[i])
	for (ii in 1:nrow(aux)) {
		circos.link(sector.index1 = aux$Var1[ii], point1 = c(aux$Zero1[ii], aux$Point1[ii]),
			sector.index2 = aux$Var2[ii], point2 = c(aux$Zero2[ii], aux$Point2[ii]),
			col = df.col$Colour[as.character(df.col$Factor) == as.character(aux$Var1[ii])])
	}
}
circos.clear()
dev.off()


#===============#
# Dusky Whalers #
#===============#
df.dusky <- subset(df, Species == "Carcharhinus obscurus")

# Get State connectivity data
trans <- unique(df.dusky$TagCode)
df.connect <- NULL
for (i in 1:length(trans)) {
	aux <- subset(df.dusky, TagCode == trans[i])
	tag.loc <- tag$Tagging.State[tag$Code2 == trans[i]]

	if (nrow(aux) > 1) {
		aux <- aux[order(aux$Timestamp.UTC), ]
		state1 <- NULL
		state2 <- NULL
		for(ii in 2:nrow(aux)) {
			if (ii == 1) {
				if (tag.loc != aux$State[ii - 1]) {
				state1 <- c(state1, tag.loc)
				state2 <- c(state2, aux$State[ii - 1])
				}
			}
			aux.first <- aux$State[ii - 1]
			aux.second <- aux$State[ii]
			if (aux.first != aux.second) {
				state1 <- c(state1, aux.first)
				state2 <- c(state2, aux.second)
			}
		}
		df.save <- data.frame(State1 = state1, State2 = state2, Transmitter = rep(trans[i], length(state1)))
		if (nrow(df.save) > 0)
			df.connect <- rbind(df.connect, df.save)
	}
}


# Identify number of connections = outcoming(+) incoming (-)
factors <- c("WA", "SA", "VIC", "NSW")
df.comb <- expand.grid(factors, factors)
df.comb$N <- NA
for (i in 1:nrow(df.comb)) {
	df.comb$N[i] <- length(which(df.connect$State1 == df.comb$Var1[i] & df.connect$State2 == df.comb$Var2[i]))
}
df.comb <- df.comb[-which(df.comb$N == 0), ]

# Identify sector widths: number of connections
widths.pos <- NULL
widths.neg <- NULL
for (i in 1:length(factors)) {
	aux1 <- sum(df.comb$N[df.comb$Var1 == factors[i]])
	aux2 <- sum(df.comb$N[df.comb$Var2 == factors[i]])

	widths.pos <- c(widths.pos, aux1)
	widths.neg <- c(widths.neg, aux2)
}
df.width <- data.frame(Factor = factors, Positive = widths.pos, Negative = widths.neg)

# Transform connections to proportions:
tot <- sum(df.comb$N)
df.width$Tot <- df.width$Positive + df.width$Negative
df.width$Width <- NA
for (i in 1:nrow(df.width)) {
	df.width$Width[i] <- df.width$Tot[i] / tot
}
df.width$Zero <- NA
for (i in 1:nrow(df.width)) {
	df.width$Zero[i] <- df.width$Positive[i] / tot
	if (df.width$Positive[i] > df.width$Negative[i]) 
		df.width$Zero[i] <- df.width$Zero[i] * -1
}

# Identify sector proportion for plotting the robbons:
df.width$Positive.length <- 1 - df.width$Zero
df.width$Negative.length <- 2 - df.width$Positive.length
df.width$Positive.prop <- df.width$Positive.length / df.width$Positive
df.width$Negative.prop <- (-1 - df.width$Zero) / df.width$Negative

# Create auxiliary point coordinates for robbon:
df.comb$Zero1 <- NA
df.comb$Point1 <- NA
df.comb$Zero2 <- NA
df.comb$Point2 <- NA
df.save <- NULL
out.fact <- unique(df.comb$Var1)
for (i in 1:length(out.fact)) {
	aux <- subset(df.comb, Var1 == out.fact[i])
	zero1 <- df.width$Zero[df.width$Factor == out.fact[i]]
	point1 <- zero1 + (df.width$Positive.prop[df.width$Factor == out.fact[i]] * aux$N[1])
	aux$Zero1[1] <- zero1
	aux$Point1[1] <- point1
	if (nrow(aux) > 1) {
		for (ii in 2:nrow(aux)) {
			zero1 <- point1
			point1 <- point1 + (df.width$Positive.prop[df.width$Factor == out.fact[i]] * aux$N[ii])
			aux$Zero1[ii] <- zero1
			aux$Point1[ii] <- point1
		}
	}
	df.save <- rbind(df.save, aux)
}
df.save2 <- NULL
inc.fact <- unique(df.save$Var2)
for (i in 1:length(inc.fact)) {
	aux <- subset(df.save, Var2 == inc.fact[i])
	zero2 <- df.width$Zero[df.width$Factor == inc.fact[i]]
	point2 <- zero2 + (df.width$Negative.prop[df.width$Factor == inc.fact[i]] * aux$N[1])
	aux$Zero2[1] <- zero2
	aux$Point2[1] <- point2
	if (nrow(aux) > 1) {
		for (ii in 2:nrow(aux)) {
			zero2 <- point2
			point2 <- point2 + (df.width$Negative.prop[df.width$Factor == inc.fact[i]] * aux$N[ii])
			aux$Zero2[ii] <- zero2
			aux$Point2[ii] <- point2
		}
	}
	df.save2 <- rbind(df.save2, aux)
}
df.comb <- df.save2


# Auxiliar object for colors
df.col <- data.frame(Factor = c("WA", "SA", "VIC", "NSW"), Colour = rainbow(4))
df.col$Colour <- paste0(df.col$Colour, 80)


### CIRCULAR PLOT
factors <- factor(c("WA", "SA", "VIC", "NSW"), levels = c("WA", "SA", "VIC", "NSW")) # Initialize layout
df.width$Width[3:4] <- 0.1

jpeg(hndl.out("Network_dusky.png"), width = 13, height = 13, units = "cm", quality = 500, bg = NA, res = 500)
par(mar = c(2, 2, 2, 2))
circos.par(canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.2, 1.2))
circos.initialize(factors, xlim = c(-1, 1), sector.width = df.width$Width)
# Plot tracks: Australian States
circos.track(ylim = c(0, 1), track.height = 0.05, bg.col = rainbow(4), bg.border = NA)
# Add axis
circos.axis(h = "top", sector.index = "WA", major.at = c(df.width$Zero[df.width$Factor == "WA"], 1), minor.ticks = 0, 
	labels = c(0, df.width$Positive[df.width$Factor == "WA"]), labels.cex = 0.8)
circos.text(x = 0, y = 10, labels = "WA", sector.index = "WA", niceFacing = TRUE)
circos.axis(h = "top", sector.index = "SA", major.at = c(df.width$Zero[df.width$Factor == "SA"], 1), minor.ticks = 0, 
	labels = c(0, df.width$Positive[df.width$Factor == "SA"]), labels.cex = 0.8)
circos.text(x = 0, y = 10, labels = "SA", sector.index = "SA", niceFacing = TRUE)
circos.axis(h = "top", sector.index = "VIC", major.at = c(df.width$Zero[df.width$Factor == "VIC"], 1), minor.ticks = 0, 
	labels = c(0, 1), labels.cex = 0.8)
circos.text(x = 0, y = 10, labels = "VIC", sector.index = "VIC", niceFacing = TRUE)
circos.axis(h = "top", sector.index = "NSW", major.at = c(df.width$Zero[df.width$Factor == "NSW"], 1), minor.ticks = 0, 
	labels = c(0, df.width$Positive[df.width$Factor == "NSW"]), labels.cex = 0.8)
circos.text(x = 0, y = 10, labels = "NSW", sector.index = "NSW", niceFacing = TRUE)

# Add robbons
out.fact <- unique(df.comb$Var1)
for (i in 1:length(out.fact)) {
	aux <- subset(df.comb, Var1 == out.fact[i])
	for (ii in 1:nrow(aux)) {
		circos.link(sector.index1 = aux$Var1[ii], point1 = c(aux$Zero1[ii], aux$Point1[ii]),
			sector.index2 = aux$Var2[ii], point2 = c(aux$Zero2[ii], aux$Point2[ii]),
			col = df.col$Colour[as.character(df.col$Factor) == as.character(aux$Var1[ii])])
	}
}
circos.clear()
dev.off()




# Matias' code ------------------------------------------------------------------

#1. Get conventional tagging data
Conv.Tagging=read.csv(hndl.in('Conv.Tagging.csv'))

#2. Parameters Section
  #define if doing exploratory analysis
do.expl="NO"

  #define arbitrary turning points
Cape.Leuwin=cbind(114.969,-34.459)
Mid.point=cbind(116.425,-35.043)
SA.Border=cbind(129.040964,-31.744653)
Eyre=c(135.694043,-34.883332)
Wilsons.prom=c(146.4, -39.1)
Eyre_SA.Border=distGeo(Eyre,SA.Border)
SA.Border_Mid.point=distGeo(SA.Border,Mid.point)
Mid.point_Cape.Leuwin=distGeo(Mid.point,Cape.Leuwin)

  #Spatial reference
Adelaide=c(138.6,-34.928)

UTC.WA=8

#3. Acoustic tags

#Data manipulations
Dat=df%>%
    left_join(tag%>%
             dplyr::select(Code2,ReleaseDate2,ReleaseSite2,
                           ReleaseLongitude2,ReleaseLatitude2),
             by=c("TagCode" = "Code2"))%>%
    rename(ReleaseDate=ReleaseDate2,
           ReleaseSite=ReleaseSite2,
           ReleaseLongitude=ReleaseLongitude2,
           ReleaseLatitude=ReleaseLatitude2,
           Station.type=Data)%>%
    mutate(ReleaseLongitude=ifelse(is.na(ReleaseLongitude) & ReleaseSite=='Spencer Gulf',
                                   136.98,ReleaseLongitude),
           ReleaseLatitude=ifelse(is.na(ReleaseLatitude) & ReleaseSite=='Spencer Gulf',
                                   -34.30,ReleaseLatitude),
           Datetime=Timestamp.UTC-(UTC.WA*3600),
           Mn=month(Datetime),
           Yr=year(Datetime),
           N=1,
           Rel.state=ifelse(ReleaseLongitude<=129,'WA',
                        ifelse(ReleaseLongitude>129 & ReleaseLongitude<=141,'SA',NA)),
           Organisation=ifelse(Rel.state=="SA",'Flinders/SARDI',
                        ifelse(Rel.state=="WA",'WA Fisheries',NA)),
           zone=ifelse(Latitude>(-26) & Longitude<115,"Ningaloo",
                ifelse(Latitude>(-33) & Latitude<=(-26) & Longitude<116.5,"WC",
                ifelse(Latitude<=(-33) & Longitude<116.5,"Zone1",
                ifelse(Longitude>=116.5 & Longitude<129 & Latitude<=(-26),"Zone2",
                ifelse(Longitude>=129 & Longitude<135.7,"SA.west",
                ifelse(Longitude>=135.6 & Longitude<=141,"SA.east",
                ifelse(Longitude>141 & Longitude<=149.7,"Vic",  
                ifelse(Longitude>149.7,"NSW",
                       NA)))))))))

  #get previous stuff
Dat=Dat%>%
    arrange(TagCode,Datetime) %>%
    mutate(TagCode.prev=lag(TagCode,1),
           Datetime.prev=lag(Datetime,1),
           Longitude.prev=lag(Longitude,1),
           Latitude.prev=lag(Latitude,1),
           State.prev=lag(State,1),
           Datetime.prev=as.POSIXct(ifelse(is.na(TagCode.prev)|!TagCode==TagCode.prev,
                                           as.character(ReleaseDate),as.character(Datetime.prev)),tz = "UTC"),
           Longitude.prev=ifelse(is.na(TagCode.prev)|!TagCode==TagCode.prev,
                                 ReleaseLongitude,Longitude.prev),
           Latitude.prev=ifelse(is.na(TagCode.prev)|!TagCode==TagCode.prev,ReleaseLatitude,Latitude.prev),
           State.prev=ifelse(is.na(TagCode.prev)|!TagCode==TagCode.prev,Rel.state,State.prev),
           TagCode.prev=ifelse(is.na(TagCode.prev),TagCode,TagCode.prev),
           zone.prev=ifelse(Latitude.prev>(-26) & Longitude.prev<115,"Ningaloo",
                     ifelse(Latitude.prev>(-33) & Latitude.prev<=(-26) & Longitude.prev<116.5,"WC",
                     ifelse(Latitude.prev<=(-33) & Longitude.prev<116.5,"Zone1",
                     ifelse(Longitude.prev>=116.5 & Longitude.prev<129 & Latitude.prev<=(-26),"Zone2",
                     ifelse(Longitude.prev>=129 & Longitude.prev<135.7,"SA.west",
                     ifelse(Longitude.prev>=135.6 & Longitude.prev<=141,"SA.east",
                     ifelse(Longitude.prev>141 & Longitude.prev<=149.7,"Vic",  
                     ifelse(Longitude.prev>149.7,"NSW",
                     NA)))))))),
           Same.station=ifelse(TagCode==TagCode.prev & 
                                 Longitude==Longitude.prev &
                                 Latitude==Latitude.prev  ,"YES","NO"),
           Time=ifelse(TagCode==TagCode.prev,difftime(Datetime,Datetime.prev,units="mins"),NA))

 #keep only those crossing jurisdictions
Dat=Dat%>%
  mutate(Crossing=ifelse(State.prev==State,"NO","YES"))

Cross.tags=Dat%>%filter(Crossing=="YES")%>%distinct(TagCode)%>%pull(TagCode)

Dat=Dat%>%filter(TagCode%in%Cross.tags)


# Report table of tagcodes by species and state
Tab1= group_by(Dat, TagCode, Species,Organisation,State) %>%
  summarise(sum = sum(N)) %>%
  as.data.frame()
write.csv(Tab1,hndl.out('Table.1_Number.detections.by.tag.and.organisation.csv'),row.names = F)


Rel.dat=Dat%>%
  distinct(TagCode,.keep_all = TRUE)%>%
  dplyr::select(TagCode,ReleaseDate,ReleaseLatitude,ReleaseLongitude)%>%
  mutate(ReleaseState=
              ifelse(ReleaseLongitude>129 & ReleaseLongitude<=141,"SA",             
              ifelse(ReleaseLongitude<=129,"WA",
              ifelse(ReleaseLongitude>141 & ReleaseLongitude<=149.7,"Vic",  
              ifelse(ReleaseLongitude>149.7,"NSW",
              NA)))))


# Data set scenarios
Dat.scen1=Dat   # full data set
Dat.scen2=subset(Dat, Station.type=='Non-IMOS') # only state

# Some functions
    #proportion of time per jurisdiction (straight line movement assumption)
fn.prop.time.jur=function(d)  
{
  d=rbind(data.frame(Datetime=d$Datetime.prev[1],
                     Longitude=d$Longitude.prev[1],
                     Latitude=d$Latitude.prev[1],
                     zone=d$zone.prev[1],
                     Juris=d$Juris.prev[1],
                     Rel.state=d$Rel.state[1]),
          d%>%
            dplyr::select(Datetime,Longitude,Latitude,zone,Juris,Rel.state))
  delta.t=max(d$Datetime)-min(d$Datetime)
  Rel.state=unique(d$Rel.state)
  d=d %>% mutate(date=date(Datetime)) %>%
    arrange(Datetime) %>%
    distinct(date,zone,.keep_all = TRUE) %>%
    mutate(Juris.prev=lag(Juris,1),
           Longitude.prev=lag(Longitude,1),
           Latitude.prev=lag(Latitude,1),
           date.prev=lag(date,1),
           delta.t=date-date.prev,
           Same.juris=ifelse(Juris==Juris.prev,"YES",
                             ifelse(!Juris==Juris.prev,"NO",NA)))
  d.diff.jur=subset(d,Same.juris=="NO")     
  d=subset(d,!Same.juris=="NO"|is.na(Same.juris)) 
  if(nrow(d.diff.jur)>0)
  {
    b1=vector('list',nrow(d.diff.jur))
    for(n in 1:length(b1))
    {
      dummy=with(d.diff.jur[n,],gcIntermediate(c(Longitude.prev,Latitude.prev),
                                               c(Longitude,Latitude),n=as.numeric(delta.t)-1, addStartEnd=T))
      dummy=data.frame(Longitude=dummy[,1],Latitude=dummy[,2])
      dummy$date=seq(d.diff.jur[n,'date.prev'],d.diff.jur[n,'date'],by='day')
      b1[[n]]=dummy
    }
    b1=do.call(rbind,b1)
    d.new=rbind(b1,subset(d,select=names(b1))) %>%
      arrange(date) %>%
      mutate(Juris=ifelse(Longitude>129 & Longitude<=141,"SA",             
                   ifelse(Longitude<=129,"WA",
                   ifelse(Longitude>141 & Longitude<=149.7,"Vic",  
                   ifelse(Longitude>149.7,"NSW",
                   NA)))),
             Juris.prev=lag(Juris,1),
             date.prev=lag(date,1),
             delta.t=date-date.prev)
    id=which(!d.new$Juris==d.new$Juris.prev)
    d.new$index=NA
    d.new$index[1:(id[1]-1)]=1
    if(length(id)>1)for(x in 2:length(id)) d.new$index[id[x-1]:(id[x]-1)]=x
    d.new$index[id[length(id)]:nrow(d.new)]=length(id)+1
    Tab1=group_by(subset(d.new,!is.na(delta.t)), index, Juris) %>% 
      summarise(Total.days = as.numeric(sum(delta.t))) %>%
      as.data.frame() 
    Tab1$Prop=with(Tab1,Total.days/sum(Total.days))
  }
  if(nrow(d.diff.jur)==0)
  {
    d=subset(d,!is.na(Juris))
    Tab1=data.frame(index=1,Juris=unique(d$Juris),Total.days=delta.t,Prop=1)
  }
  Tab1$Rel.state=Rel.state
  return(Tab1)
}

    #proportion of time per jurisdiction
fn.plt.prop.time=function(d,CL1,CL2,CL3,CL4)
{
  plot(1:length(d),xlim=c(0,1),ylim=c(0,length(d)),col="transparent",ylab="",
       xlab="",yaxt='n',xaxt='n',bty='n')
  for(n in 1:length(d))
  {
    d[[n]]$col=with(d[[n]],ifelse(Juris=="WA",CL1,ifelse(Juris=="SA",CL2,ifelse(Juris=='Vic',CL3,CL4))))
    d[[n]]$cum.prop=cumsum(d[[n]]$Prop)
    if(nrow(d[[n]])>1)Total.days=sum(d[[n]]$Total.days)
    if(nrow(d[[n]])==1)Total.days=d[[n]]$Total.days
    
    for(r in 1:nrow(d[[n]]))
    {
      if(r==1) x1=0 else
        x1=d[[n]]$cum.prop[r-1]
      x2=d[[n]]$cum.prop[r]
      polygon(x=c(x1,rep(x2,2),x1),y=c(n-0.5,n-0.5,n,n),col=d[[n]]$col[r])
      text(mean(c(x1,x2)),mean(c(n-0.5,n)),
           paste(round(d[[n]]$Prop[r],2)*100,'%',sep=''),cex=.75)   
    }
    Clr=ifelse(d[[n]]$Rel.state[1]=="WA",CL1,CL2)
    text(0,n-.25,d[[n]]$Rel.state[1],pos=2,col=Clr,font=2)
    text(1,n-.25,round(Total.days),pos=4,cex=.85)
    text(1,n-.25,names(d)[n],pos=3,col='firebrick',cex=.7)
    
  }
}

    #residency
fn.plt.residency=function(d,CL1,CL2,CL3,CL4)
{
  plot(1:length(d),xlim=c(0,1),ylim=c(0,length(d)),col="transparent",ylab="",
       xlab="",yaxt='n',xaxt='n',bty='n')
  for(n in 1:length(d))
  {
    d[[n]]$col=with(d[[n]],ifelse(Juris=="WA",CL1,ifelse(Juris=="SA",CL2,ifelse(Juris=='Vic',CL3,CL4))))
    d[[n]]$cum.prop=cumsum(d[[n]]$Prop)
    
    for(r in 1:nrow(d[[n]]))
    {
      if(r==1) x1=0 else
        x1=d[[n]]$cum.prop[r-1]
      x2=d[[n]]$cum.prop[r]
      polygon(x=c(x1,rep(x2,2),x1),y=c(n-0.5,n-0.5,n,n),col=d[[n]]$col[r])
      text(mean(c(x1,x2)),mean(c(n-0.5,n)),
           paste(round(d[[n]]$Prop[r],2)*100,'%',sep=''),cex=.75)   
    }
    Clr=ifelse(d[[n]]$Rel.state[1]=="WA",CL1,CL2)
    text(0,n-.25,d[[n]]$Rel.state[1],pos=2,col=Clr,font=2)
    text(1,n-.4,names(d)[n],pos=3,col='forestgreen',cex=.85)
    
  }
}

    #functions for plotting seasonality
fn.ggplot=function(dd,Y,X,Y.lab,X.lab)
{
  p=ggplot(dd,aes_string(y=Y,x=X,col='TagCode'))+
    geom_line()+geom_point()+scale_y_reverse()+
    labs(y = Y.lab,x = X.lab)
  if(X=='Mn')
  {
    p=p+scale_x_continuous(breaks = seq(1,12,2),
                           label = month.abb[seq(1,12,2)])
  }
  p
}

    #function for predicting GAM
pred.fun=function(MOD,DaT,Predictor)   
{
  Id.p=match(Predictor,names(DaT))
  Min.x=min(DaT[,Id.p])
  Max.x=max(DaT[,Id.p])
  Tab=sort(table(DaT$TagCode))
  newd=expand.grid(X=seq(Min.x,Max.x,length.out=100),
                   TagCode=factor(names(Tab)[length(Tab)],
                                  levels(DaT$TagCode)))
  colnames(newd)[match('X',names(newd))]=Predictor
  
  pred <- predict.gam(MOD,newd,type = "response", se.fit = TRUE)
  
  sig=sigma(MOD)  #log transformation bias correction
  newd=newd%>%
    mutate(fit=pred$fit,
           se.fit=pred$se.fit,
           rel.fit=exp(fit)*exp(sig^2/2)/mean(exp(fit)*exp(sig^2/2)),
           rel.upper=exp(fit+1.96*se.fit)*exp(sig^2/2)/mean(exp(fit)*exp(sig^2/2)),
           rel.lower=exp(fit-1.96*se.fit)*exp(sig^2/2)/mean(exp(fit)*exp(sig^2/2)))
  return(newd)
}

    #function for plotting GAM predictions
plt.pred.fn=function(newd)
{
  with(newd,
       {
         plot(Mn,rel.fit,ylab="",xlab="",col="transparent",
              ylim=c(min(rel.lower),max(rel.upper)))
         polygon(x=c(Mn,rev(Mn)),
                 y=c(rel.lower,rev(rel.upper)),
                 col='grey80')
         lines(Mn,rel.fit,lwd=2)
       })
}

    #function for wrapping scenarios     
fun.run.scen=function(Dat,SCEN)
{
  #1. create useful objects
  state=unique(Dat$State)
  TAG=unique(Dat$TagCode)
  Sp=unique(Dat$Species)
  TAG.species=vector('list',length(Sp))
  names(TAG.species)=Sp
  for(t in 1:length(TAG.species)) TAG.species[[t]]=unique(subset(Dat,Species==Sp[t])$TagCode)
  
  #2. preliminary stuff
  Lon.range=range(Dat$Longitude,na.rm=T)
  Lon.range[2]=ifelse(Lon.range[2]<129,140,Lon.range[2])
  Lat.range=range(Dat$Latitude,na.rm=T)
  Lat.range[1]=ifelse(Lat.range[1]>(-35),-38,Lat.range[1])
  XLIM=c(Lon.range[1],Lon.range[2])
  YLIM=c(Lat.range[1],Lat.range[2])
  if(do.expl=="YES")
  {
    fn.plt1=function(TG)
    {
      a=subset(Dat,TagCode==TG)
      plot(a$Longitude,a$Latitude,ylim=YLIM,xlim=XLIM,las=1,ylab="Lat",xlab="Long",
           main=paste(unique(a$Species),"_",unique(a$TagCode)," (n=",nrow(a),
                      " detections)",sep=""))
    }
    pdf(hndl.out("Exploratory.pdf")) 
    sapply(TAG,fn.plt1)
    dev.off() 
  }
  
  #3. straight line distances (in km) between consecutive detections
  #       applying algorithm to avoid going over land
  Dat$Distance=distGeo(Dat[,c("Longitude.prev","Latitude.prev")],
                       Dat[,c("Longitude","Latitude")])
  Dat$Distance.c=Dat$Distance
  Dat$Distance.c=ifelse(Dat$zone.prev=='SA.east' & Dat$zone=='Zone2',
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Eyre)+
                        Eyre_SA.Border+
                        distGeo(SA.Border,Dat[,c("Longitude","Latitude")]),
                  ifelse(Dat$zone.prev=='SA.east' & Dat$zone=='Zone1' &
                        Dat$Longitude>= Cape.Leuwin[1] & Dat$Latitude<= Cape.Leuwin[2],
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Eyre)+
                        Eyre_SA.Border+SA.Border_Mid.point+
                        distGeo(Mid.point,Dat[,c("Longitude","Latitude")]),   
                  ifelse(Dat$zone.prev=='SA.east' & Dat$zone%in%c('Zone1','WC') & 
                        Dat$Latitude>Cape.Leuwin[2],
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Eyre)+
                        Eyre_SA.Border+SA.Border_Mid.point+Mid.point_Cape.Leuwin+
                        distGeo(Cape.Leuwin,Dat[,c("Longitude","Latitude")]), 
                  ifelse(Dat$zone.prev=='Zone2' & Dat$zone=='SA.east',
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],SA.Border)+
                        Eyre_SA.Border+distGeo(Eyre,Dat[,c("Longitude","Latitude")]),   
                  ifelse(Dat$zone.prev=='Zone2' & Dat$zone=='Zone1' &
                        Dat$Longitude>= Cape.Leuwin[1] & Dat$Latitude<= Cape.Leuwin[2],
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Mid.point)+
                        distGeo(Mid.point,Dat[,c("Longitude","Latitude")]),
                  ifelse(Dat$zone.prev=='Zone2' & Dat$zone%in%c('Zone1','WC') & 
                        Dat$Latitude>Cape.Leuwin[2],
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Mid.point)+
                        Mid.point_Cape.Leuwin+distGeo(Cape.Leuwin,Dat[,c("Longitude","Latitude")]),
                  ifelse(Dat$zone.prev%in%c('Zone1','WC') & Dat$zone=='SA.east' & 
                        Dat$Latitude.prev> Cape.Leuwin[2],
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Cape.Leuwin)+
                        Mid.point_Cape.Leuwin+SA.Border_Mid.point+Eyre_SA.Border+
                        distGeo(Eyre,Dat[,c("Longitude","Latitude")]), 
                  ifelse(Dat$zone.prev%in%c('Zone1') & Dat$zone=='SA.east' & 
                        Dat$Latitude.prev<= Cape.Leuwin[2],
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Mid.point)+
                        SA.Border_Mid.point+Eyre_SA.Border+
                        distGeo(Eyre,Dat[,c("Longitude","Latitude")]), 
                   ifelse(Dat$zone.prev%in%c('Zone1','WC') & Dat$zone=='Zone2' &
                        Dat$Longitude< Cape.Leuwin[1] & Dat$Latitude> Cape.Leuwin[2],
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Cape.Leuwin)+
                        Mid.point_Cape.Leuwin+
                        distGeo(Mid.point,Dat[,c("Longitude","Latitude")]),
                   ifelse(Dat$zone.prev%in%c('Zone1') & Dat$zone=='Zone2' &
                        Dat$Longitude>= Cape.Leuwin[1] & Dat$Latitude<= Cape.Leuwin[2],
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Mid.point)+
                        distGeo(Mid.point,Dat[,c("Longitude","Latitude")]),
                   ifelse(Dat$zone.prev%in%c('Zone1','WC') & Dat$zone=='Zone1' &
                        Dat$Longitude< Cape.Leuwin[1] & Dat$Latitude> Cape.Leuwin[2],
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Cape.Leuwin)+
                        distGeo(Cape.Leuwin,Dat[,c("Longitude","Latitude")]),
                  ifelse(Dat$zone.prev=='NSW' & Dat$zone=='SA.east',
                        distGeo(Dat[,c("Longitude.prev","Latitude.prev")],Wilsons.prom)+
                        distGeo(Wilsons.prom,Dat[,c("Longitude","Latitude")]),
                    Dat$Distance.c))))))))))))
  
  Dat$Distance=with(Dat,ifelse(TagCode==TagCode.prev,Distance/1000,NA)) 
  Dat$Distance.c=with(Dat,ifelse(TagCode==TagCode.prev,Distance.c/1000,NA)) 
  
  
  #4. Distribution of displacements and ROM across jurisdictions
  Dat= Dat%>% 
    mutate(Juris=ifelse(zone%in%c("Ningaloo","WC","Zone1","Zone2"),"WA",
                 ifelse(zone%in%c("SA.east","SA.west"),"SA",
                 ifelse(zone=="NSW","NSW",
                 ifelse(zone=="Vic","Vic",NA)))),
           Juris.prev=ifelse(zone.prev%in%c("Ningaloo","WC","Zone1","Zone2"),"WA",
                      ifelse(zone.prev%in%c("SA.east","SA.west"),"SA",
                      ifelse(zone=="NSW","NSW",
                      ifelse(zone=="Vic","Vic",NA)))),
           Same.juris=ifelse(Juris==Juris.prev,"YES",
                      ifelse(!Juris==Juris.prev,"NO",NA)))
  
  Cros.jur=Dat%>%filter(Same.juris=="NO")%>%
    dplyr::select(Species,TagCode,TagCode.prev,Juris,Juris.prev,
                  Longitude,Latitude,Longitude.prev,Latitude.prev,
                  Datetime,Datetime.prev,Distance.c,Time)%>%
    mutate(Time=ifelse(is.na(Time),difftime(Datetime,Datetime.prev,units='mins'),Time))
  
  Cros.jur$Distance.c=ifelse(is.na(Cros.jur$Distance.c),distGeo(Cros.jur[,c("Longitude.prev","Latitude.prev")],
                                                                Cros.jur[,c("Longitude","Latitude")])/1000,Cros.jur$Distance.c)
  Cros.jur=Cros.jur%>%mutate(ROM=Distance.c/(Time/(24*60)))  # km/day
  
  tiff(file=hndl.out(paste('Figure.2_Cross.juris_hist_figure_',SCEN,'.tiff',sep='')),
       width=2400,height=2400,units="px",res=300,compression="lzw+p")
  par(mfrow=c(2,1),mar=c(1,2,1,.1),oma=c(.5,3,.1,.5),las=1,
      mgp=c(2,.5,0),cex.axis=1.25,cex.lab=1.5,xpd=T)
  
  #distance
  vioplot(Distance.c~Species,Cros.jur, horizontal=F, col="gray",ylab='')
  mtext("Distance (km)",2,line=3.25,las=3,cex=1.75)
  
  #ROM
  vioplot(ROM~Species,Cros.jur, horizontal=F, col="gray",ylab='')
  mtext("ROM (km/day)",2,line=3.25,las=3,cex=1.75)
  dev.off()
  
  ROM.sumery=tapply(Cros.jur$ROM, Cros.jur$Species, summary)
  Dist.sumery=tapply(Cros.jur$Distance.c, Cros.jur$Species, summary)
  Tab2<-rbind(do.call(rbind,Dist.sumery),do.call(rbind,ROM.sumery))
  rownames(Tab2)=paste(c(rep('Distance',2),rep('ROM',2)),rownames(Tab2),sep='_')
  write.csv(Tab2,paste('Table.2_Cross.juris_summary_',SCEN,'.csv',sep=''))
  
  
  #5. cross jurisdictional displacements  
  Prop.time=vector('list',length(TAG))
  names(Prop.time)=TAG
  for(i in 1:length(TAG))
  {
    Prop.time[[i]]=fn.prop.time.jur(d=subset(Dat,TagCode==TAG[i],
                                             select=c(Datetime,Longitude,Latitude,zone,Juris,Rel.state,
                                                      Datetime.prev,Longitude.prev,Latitude.prev,zone.prev,Juris.prev)))
  }
  
  tiff(file=hndl.out(paste('Figure.3_Prop.time.in.jurisdiction_',SCEN,'.tiff',sep='')),width=1800,height=2400,units="px",res=300,
       compression="lzw+p")
  layout(matrix(c(1,2),ncol=1), heights=c(2,1))
  par(mar=c(1,1,1,1.25),oma=c(.1,.1,.1,.5),las=1,
    mgp=c(1.25,.35,0),cex.axis=1.1,cex.lab=1.25,xpd=T)
  for(t in 1:length(TAG.species)) 
  {
    id=which(names(Prop.time)%in%TAG.species[[t]])
    fn.plt.prop.time(d=Prop.time[id],CL1=Cols[1],CL2=Cols[2],CL3=Cols[3],CL4=Cols[4])
    if(t==2)legend('bottom',names(Cols),fill=Cols,horiz=T,bty='n')
    mtext(names(TAG.species)[t],3,cex=1.5)
  }
  mtext("Percent of time per zone",1,line=0,cex=1.5)
  dev.off()
  
  
  #6. seasonal patterns plot distance from release thru time 
  Spics=Dat%>%distinct(TagCode,.keep_all = T)%>%dplyr::select(TagCode,Species)
  Seasonal=Dat%>%
  #  left_join(subset(Rel.dat,select=c(TagCode,ReleaseDate,
  #                                    ReleaseLatitude,ReleaseLongitude)),by='TagCode')%>%
    arrange(TagCode,Datetime)%>%
    mutate(Delta.t=as.numeric(difftime(Datetime,ReleaseDate,units='days')),
           AdelaideLon=Adelaide[1],
           AdelaideLat=Adelaide[2])
  Seasonal$Dist.frm.rel=ifelse(Seasonal$TagCode==Seasonal$TagCode.prev,
                               distGeo(Seasonal[,c("ReleaseLongitude","ReleaseLatitude")],
                                       Seasonal[,c("Longitude","Latitude")])/1000,NA)
  Seasonal$Dist.frm.Adld=ifelse(Seasonal$TagCode==Seasonal$TagCode.prev,
                                distGeo(Seasonal[,c("AdelaideLon","AdelaideLat")],
                                        Seasonal[,c("Longitude","Latitude")])/1000,NA)
  
  Ril.dt=subset(Rel.dat,TagCode%in%unique(Seasonal$TagCode))
  Add.relis=Seasonal[1:length(unique(Seasonal$TagCode)),]
  Add.relis[,]=NA
  Add.relis=Add.relis%>%
    mutate(TagCode=Ril.dt$TagCode,
           TagCode.prev=Ril.dt$TagCode,
           Datetime=Ril.dt$ReleaseDate,
           Dist.frm.rel=0,
           ReleaseLatitude=Ril.dt$ReleaseLatitude,
           ReleaseLongitude=Ril.dt$ReleaseLongitude,
           AdelaideLon=Adelaide[1],
           AdelaideLat=Adelaide[2],
           Delta.t=0)%>%
    dplyr::select(-Species)%>%
    left_join(Spics,by='TagCode')%>%
    dplyr::select(names(Seasonal))
  Add.relis$Dist.frm.Adld=ifelse(Add.relis$TagCode==Add.relis$TagCode.prev,
                                 distGeo(Add.relis[,c("AdelaideLon","AdelaideLat")],
                                         Add.relis[,c("ReleaseLongitude","ReleaseLatitude")])/1000,NA)
  
  Seasonal=rbind(Seasonal,Add.relis)%>%
    arrange(TagCode,Delta.t)%>%
    dplyr::select(TagCode,TagCode.prev,Species,Datetime,
                  Dist.frm.rel,Dist.frm.Adld,Delta.t)%>%
    filter(!is.na(Dist.frm.rel))%>%
    mutate(Mn=month(Datetime),
           Yr=year(Datetime),
           ln.Dist.frm.rel=log(Dist.frm.rel+1e-4),
           ln.Dist.frm.Adld=log(Dist.frm.Adld),
           scaled.Dist.frm.Adld=scale(Dist.frm.Adld))
  
  #6.1 Bronzie    
  Seasonal.bronzie=Seasonal%>%
    filter(Species=='Carcharhinus brachyurus')
  Sisonls.bronzie=c(27698, 29454, 29542, 29587, 30717, 30870, 30894, 30992, 31000, 31003, 33180)
  Seasonal.bronzie.Mod=Seasonal.bronzie%>%
    filter(TagCode%in%Sisonls.bronzie)%>%
    mutate(TagCode=as.factor(TagCode))%>%
    arrange(TagCode,Datetime)
  #check data first
  tiff(file=hndl.out(paste('Figure.4_Seasonality_Distance.fom.Adelaide_Bronzie_',SCEN,'.tiff',sep='')),
       width=2400,height=1800,units="px",res=300,compression="lzw+p")
  fn.ggplot(dd=Seasonal.bronzie.Mod,
            Y='Dist.frm.Adld',
            X='Mn',
            Y.lab="Distance from Adelaide (km)",
            X.lab="Month")
  dev.off()
  
  #fit model
  Mod.bronzie=gam(ln.Dist.frm.Adld~s(Mn,bs='cc', k = 12)+s(TagCode,bs='re'),data=Seasonal.bronzie.Mod)
  
  #predict normalized data
  Bronzie.preds=pred.fun(MOD=Mod.bronzie,DaT=Seasonal.bronzie.Mod,Predictor="Mn")
  
  #Plot predictions
  tiff(file=hndl.out(paste('Figure.5_Seasonality_Distance.fom.Adelaide_Bronzie_model.pred_',SCEN,'.tiff',sep='')),
       width=2400,height=1800,units="px",res=300,compression="lzw+p")
  plt.pred.fn(Bronzie.preds)
  dev.off()  
  
  #6.2 Dusky
  Seasonal.dusky=Seasonal%>%
    filter(Species=='Carcharhinus obscurus')
  Sisonls.dusky=c(49144,52644,49146)
  Seasonal.dusky.Mod=Seasonal.dusky%>%
    filter(TagCode%in%Sisonls.dusky)%>%
    mutate(TagCode=as.factor(TagCode))%>%
    arrange(TagCode,Datetime)
  
  #check data first
  tiff(file=hndl.out(paste('Figure.4_Seasonality_Distance.fom.Adelaide_Dusky_',SCEN,'.tiff',sep='')),
       width=2400,height=1800,units="px",res=300,compression="lzw+p")
  fn.ggplot(dd=Seasonal.dusky.Mod,
            Y='Dist.frm.Adld',
            X='Mn',
            Y.lab="Distance from Adelaide (km)",
            X.lab="Month")
  dev.off()
  
  #fit model
  n.knots=length(table(Seasonal.dusky.Mod$Mn))
  Mod.dusky=gam(ln.Dist.frm.Adld~s(Mn,bs='cc', k = n.knots)+s(TagCode,bs='re'),data=Seasonal.dusky.Mod)
  
  #predict normalized data
  Dusky.preds=pred.fun(MOD=Mod.dusky,DaT=Seasonal.dusky.Mod,Predictor="Mn")
  
  #Plot predictions
  tiff(file=hndl.out(paste('Figure.5_Seasonality_Distance.fom.Adelaide_Dusky_model.pred_',SCEN,'.tiff',sep='')),
       width=2400,height=1800,units="px",res=300,compression="lzw+p")
  plt.pred.fn(newd=Dusky.preds)
  dev.off()
  
  #7. Residency
  Residency=Dat%>%
    arrange(TagCode,Datetime)%>%
    mutate(event=ifelse(TagCode==TagCode.prev & Same.juris=="NO",1,0),
           event=ifelse(TagCode==TagCode.prev,cumsum(event),NA),
           TagCode.Juris=paste(TagCode,Juris,Juris.prev))%>%
    group_by(TagCode) %>%
    mutate(Tot.time = difftime(max(Datetime),min(Datetime.prev),units='mins'))%>%   
    group_by(TagCode,event,TagCode.Juris)%>%
    mutate(Time.in.zn=ifelse(TagCode==TagCode.prev & Same.juris=="YES",
                             difftime(max(Datetime),min(Datetime.prev),units='mins'),NA))%>%
    mutate(Residency=Time.in.zn/as.numeric(Tot.time))%>%
    distinct(TagCode,event,TagCode.Juris,Residency,.keep_all = T)%>%
    dplyr::select(TagCode,Species,Rel.state,Juris,Residency,event,TagCode.Juris)%>%
    filter(!is.na(Residency))%>%
    data.frame()%>%mutate(Event.Juris=paste(event,Juris))
  
  Prop.time.Res=vector('list',length(unique(Residency$Species)))
  names(Prop.time.Res)=unique(Residency$Species)
  for(j in 1:length(Prop.time.Res))
  {
    d=Residency%>%filter(Species==names(Prop.time.Res)[j])
    this=unique(d$TagCode)
    dummy=vector('list',length(this))
    names(dummy)=this
    for(xx in 1:length(dummy))
    {
      dummy[[xx]]=d%>%filter(TagCode==names(dummy)[xx])%>%
        mutate(index=1:length(event),
               Prop=Residency,
               Total.days='')%>%
        dplyr::select(index,Juris,Total.days,Prop,Rel.state)
    }
    Prop.time.Res[[j]]=dummy
  }
  
  tiff(file=hndl.out(paste('Figure.6_Residency_',SCEN,'.tiff',sep='')),
       width=1800,height=2400,units="px",res=300,compression="lzw+p")
  layout(matrix(c(1,2),ncol=1), heights=c(2,1))
  par(mar=c(1,1,1,1.25),oma=c(.1,.1,.1,.5),las=1,
      mgp=c(1.25,.35,0),cex.axis=1.1,cex.lab=1.25,xpd=T)
  #par(mfcol=c(2,1),mar=c(1,1,1,1.25),oma=c(.1,.1,.1,.5),las=1,
  #    mgp=c(1.25,.35,0),cex.axis=1.1,cex.lab=1.25,xpd=T)
  for(j in 1:length(Prop.time.Res))
  {
    fn.plt.residency(d=Prop.time.Res[[j]],CL1=Cols[1],CL2=Cols[2],CL3=Cols[3],CL4=Cols[4])
    if(j==2)legend('bottom',names(Cols),fill=Cols,horiz=T,bty='n')
    mtext(names(Prop.time.Res)[j],3,cex=1.5)
  }
  mtext("Residency",1,line=0,cex=1.5)
  dev.off()
  
}

# Run scenarios
Cols=c("steelblue","pink2","forestgreen","orange")
names(Cols)=c("WA","SA","Vic","NSW")

fun.run.scen(Dat=Dat.scen1,SCEN="Scen1")
fun.run.scen(Dat=Dat.scen2,SCEN="Scen2")



# 4. Conventional tags

  #summary of releases and recaptures by jurisdiction
TAB.conv.rel.rec=Conv.Tagging%>%
  group_by(Rel.juris,Rec.juris,COMMON_NAME,Recaptured)%>%
  summarise(n=n())%>%
  filter(!is.na(Rec.juris))%>%
  arrange(COMMON_NAME,Rel.juris,Recaptured)%>%
  data.frame%>%
  mutate(Rec.juris=ifelse(Recaptured=="No","N/A",Rec.juris))
write.csv(TAB.conv.rel.rec,hndl.out("Table.3_Conventional.tag_summary.csv"),row.names = F)

  #Plot recaptures
data(worldLLhigh)

library(shape)
SA.WA.lat=c(-40,-29)
SA.WA.long=c(112,140)
tiff(hndl.out("Figure.7.Conventional.tag.map.tiff"),width=2400,height=1200,
     units="px",res=300,compression="lzw")
par(mar=c(.1,.1,.1,1),oma=c(.1,.1,.1,1))
plotMap(worldLLhigh, xlim=SA.WA.long,ylim=SA.WA.lat,plt = c(.1, 1, 0.075, 1),
        col="grey85",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
box()
with(subset(Conv.Tagging,Recaptured=="Yes" & Same.Juris=="No"),{
  Arrows(Long.rels,Lat.rels,Long.rec,Lat.rec,col=Col.sp,lwd=2, 
         arr.type="curved", arr.width=.2)
  points(Long.rels,Lat.rels,pch=19,col=Col.sp,cex=1.5)
  
})
axis(side = 1, at =SA.WA.long[1]:SA.WA.long[2], labels = F, tcl = 0.15)
axis(side = 2, at = SA.WA.lat[2]:SA.WA.lat[1], labels = F,tcl =0.15)
n=seq(SA.WA.long[1],SA.WA.long[2],4)
axis(side = 1, at =n, labels = n, tcl = 0.3,padj=-1.75)
n=seq(SA.WA.lat[1]+1,SA.WA.lat[2],4)
axis(side = 2, at = n, labels = abs(n),las=2,tcl =0.3,hadj=.2)

OZ.lat=c(-44,-10)
OZ.long=c(112,154)
par(fig=c(.5,.85,.125,.50), new = T,mgp=c(.1,.4,0))
plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
        col="grey85",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
box()
polygon(x=c(SA.WA.long,rev(SA.WA.long)),
        y=c(rep(SA.WA.lat[1],2),rep(SA.WA.lat[2],2)),lwd=1.5,col=rgb(.1,.1,.1,alpha=.5))
text(134,-22.5,("Australia"),col="black", cex=1.3)
mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=-1,font=1,las=0,cex=1.2,outer=T)
mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=-2,font=1,las=0,cex=1.2,outer=T)
dev.off()

  #Size and time at liberty of those recaptured in SA
Siz.sex.cross.conv=Conv.Tagging%>%
  filter(Recaptured=="Yes" & Same.Juris=="No")%>%
  mutate(Date.rel=as.POSIXct(paste(Yr.rel,Mn.rel,Day.rel,sep="-")),
         Date.rec=as.POSIXct(paste(Yr.rec,Mn.rec,Day.rec,sep="-")),
         time.at.liberty=Date.rec-Date.rel)%>%
  dplyr::select(Species,COMMON_NAME,Lat.rels,Long.rels,Lat.rec,Long.rec,
                Rel.juris,Rec.juris,Rel_FL,Sex,Date.rel,
                Date.rec,time.at.liberty)%>%
  arrange(Species,Date.rel,Sex)
write.csv(Siz.sex.cross.conv,hndl.out("Table.4_Conventional.tag_size.time.liberty_recaptured in SA.csv"),row.names = F)
