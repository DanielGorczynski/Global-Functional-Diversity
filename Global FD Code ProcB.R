### Code for Gorczynski et al. "Tropical mammal functional diversity ###
### increases with productivity but decreases with anthropogenic disturbance ###


##Calculate FD (FRich and FDis) for all sites
##both realized community and species pools
library(tidyr)
library(dplyr)

##Import species traits
spp.pool.traits <- read.csv("Full trait list bi 1 kg.csv")
rownames(spp.pool.traits) <- spp.pool.traits$Species
spp.pool.traits$Species <- NULL

##Import Species pool list
Species1 <- read.csv("Species1.csv")

#Import Occupancy data
All.TEAM.Occupancy <- read.csv("All TEAM Occupancy.csv")
colnames(All.TEAM.Occupancy) <- c("Site","Species","2007","2008","2009","2010","2011","2012","2013","2014")
Occ <- gather(All.TEAM.Occupancy, "Year", "Occupancy", 3:10)
#Include next line to test for unweighted functional dispersion
#Occ$Occupancy <- 1
Species1$Occupancy <- 1
Species1$Year <- "Spp Pool"
Species1 <- rbind(Species1,Occ)
Species1$Community <- paste(Species1$Site, Species1$Year)
Species1 <- Species1[,-c(2,4)]
Species1 <- Species1 %>% spread(Community,Occupancy)
Species1[is.na(Species1)] <- 0
Species1 <- filter(Species1, Species1$Species %in% rownames(spp.pool.traits))
rownames(Species1) <- Species1$Species
Species1$Species <- NULL
Species1 <- t(Species1)
Species1 <- Species1[which(rowSums(Species1) > 0),] 


spp.pool.traits <- spp.pool.traits[order(rownames(spp.pool.traits)),]
spp.pool.traits$Body.Mass <- log(spp.pool.traits$Body.Mass)
spp.pool.traits$Litter.Size <- log(spp.pool.traits$Litter.Size)
##Include to calculate FRich with locally extripated species
##BCI## row 7
#White-lipped peccary- col 226
#Species1[7,153] <- 1
#Giant Anteater- col 146
#Species1[7,111] <- 1
#Jaguar- col 159
#Species1[7,120] <- 1

##NAK## row 40
#Tiger- col 161
#Species1[40,122] <- 1
#Leopard- col 160
#Species1[40,121] <- 1

##BIF## row 13
#Buffalo- col 222
#Species1[13,149] <- 1
#Leopard- col 160
#Species1[13,121] <- 1
#Giant forest hog- col 97
#Species1[13,77] <- 1

##KRP## row 36
#Leopard- col 160
#Species1[36,121] <- 1
#Golden Cat- col 19
#Species1[36,16] <- 1
#Giant Pangolin- col 214
#Species1[36,145] <- 1

library(FD)
x <- dbFD(spp.pool.traits, Species1, w = c(15,3,3,3,3,3,15,15,5,5,5,15))

##dbFD function without weighting to treat traits equally
#x <- dbFD(spp.pool.traits, Species1)


##Import table with predictor variables 
TEAM.Site.Funct.bi <- read.csv("TEAM Site Funct bi.csv")

##Correlation test for predictor variables
Cor.test <- TEAM.Site.Funct.bi[,c(2,3,4,5,7,13,14,15)]
cor(Cor.test)

#Select data and response variable
dat <- TEAM.Site.Funct.bi



##Scale Predictor Variables
dat$ndvi <- scale(dat$ndvi)
dat$Landcover.div <- scale(dat$Landcover.div) 
dat$hum.dens.ZOI <- scale(dat$hum.dens.ZOI) 
dat$FPoolRich <- scale(dat$FPoolRich) 
dat$FR.Lost <- scale(dat$FR.Lost)

#Bayesian Linear Modeling
library(rstan)
library(brms)
library(ggplot2)
library(ggdist)
library(tidybayes)

#######Functional richness model

model.FR1 <- brm(formula = Frich ~ FPoolRich + Landcover.div + ndvi + hum.dens.ZOI + FR.Lost + Africa + Asia + Madagascar, 
             data = dat, control = list(adapt_delta = 0.99, max_treedepth = 15), family = "weibull", iter = 4000)

summary(model.FR1)


##Plot Credible Intervals for Functional Richness
x <- model.FR1 %>%
  spread_draws(b_Intercept, b_FPoolRich, b_Landcover.div, b_ndvi, b_hum.dens.ZOI,b_FR.Lost, b_Africa, b_Asia, b_Madagascar)%>%
  gather(key = "b_", value = "Est" ) #%>%
x <- x[-c(1:24000),]
x$b_ <- factor(x$b_, levels = c("b_ndvi", "b_Landcover.div",
                                "b_hum.dens.ZOI","b_FR.Lost",
                                "b_FPoolRich","b_Africa",
                                "b_Asia","b_Madagascar","b_Intercept"))
  ggplot(data = x, aes(y = b_, x = Est, fill = b_)) +
  stat_halfeye()+
    geom_vline(xintercept = 0)+
    labs(title="Functional Richness model",
          x ="Beta effect size", y = "Predictor variables")+
    scale_y_discrete(labels = c("b_ndvi" = "NDVI", "b_Landcover.div" = "Landcover Diversity",
                                "b_hum.dens.ZOI" = "Human Density","b_FR.Lost"="Extinction",
                                "b_FPoolRich"="Species Pools","b_Africa"="Africa",
                                "b_Asia"="Asia","b_Madagascar"="Madagascar","b_Intercept"="Intercept"))+
    scale_fill_manual(values = c("grey85","grey85","grey85","grey35","grey85","grey85","grey85","grey35","grey85"))+
    theme_classic()+
    theme(legend.position = "none")

##Posterior predictive checks
  pp_check(model.FR1, type ="dens_overlay", nsamples = 50)
  pp_check(model.FR1, type ="error_hist", nsamples = 9)
  
  
#######Functional dispersion model

model.FDis <- brm(formula = Fdis ~ FPoolRich + Landcover.div + ndvi + FR.Lost + hum.dens.ZOI + Africa + Asia + Madagascar, 
             data = dat,  control = list(adapt_delta = 0.99, max_treedepth = 20), family = "weibull", iter = 4000)

summary(model.FDis)


##Plot Credible Intervals for Functional dispersion

x <- model.FDis %>%
  spread_draws(b_Intercept, b_FPoolRich, b_Landcover.div, b_ndvi, b_hum.dens.ZOI,b_FR.Lost, b_Africa, b_Asia, b_Madagascar)%>%
  gather(key = "b_", value = "Est" ) #%>%
x <- x[-c(1:24000),]
x$b_ <- factor(x$b_, levels = c("b_ndvi", "b_Landcover.div",
                                "b_hum.dens.ZOI","b_FR.Lost",
                                "b_FPoolRich","b_Africa",
                                "b_Asia","b_Madagascar","b_Intercept"))
ggplot(data = x, aes(y = b_, x = Est, fill = b_)) +
  stat_halfeye()+
  geom_vline(xintercept = 0)+
  labs(title="Functional Dispersion model",
       x ="Beta", y = "Predictor variables")+
  scale_y_discrete(labels = c("b_ndvi" = "NDVI", "b_Landcover.div" = "Landcover Diversity",
                              "b_hum.dens.ZOI" = "Human Density","b_FR.Lost"="Extinction",
                              "b_FPoolRich"="Species Pools","b_Africa"="Africa",
                              "b_Asia"="Asia","b_Madagascar"="Madagascar", "b_Intercept"="Intercept"))+
  scale_fill_manual(values = c("grey35","grey85","grey85","grey85","grey85","grey85","grey85","grey85","grey85"))+
  theme_classic()+
  theme(legend.position = "none")

##Posterior predictive checks
pp_check(model.FDis, type ="dens_overlay", nsamples = 50)
pp_check(model.FDis, type ="error_hist", nsamples = 9)


#####Comparison of NDVI effect for weighted and unweighted FDis metrics

model.FDis <- brm(formula = Fdis ~ FPoolRich + Landcover.div + ndvi + FR.Lost + hum.dens.ZOI + Africa + Asia + Madagascar, 
                  data = dat,  control = list(adapt_delta = 0.999, max_treedepth = 20), family = "weibull", iter = 4000)

model.FDisUW <- brm(formula = FDisUW ~ FPoolRich + Landcover.div + ndvi + hum.dens.ZOI + FR.Lost + Africa + Asia + Madagascar, 
             data    = dat,  control = list(adapt_delta = 0.999, max_treedepth = 20), iter = 4000)

model.FDisSPUW <- brm(formula = FDisSPUW ~ FPoolRich + Landcover.div + ndvi + hum.dens.ZOI + FR.Lost + Africa + Asia + Madagascar, 
                      data    = dat,  control = list(adapt_delta = 0.999, max_treedepth = 20), iter = 4000)

x.FDis <- model.FDis %>%
  spread_draws(b_ndvi)%>%
  gather(key = "b_", value = "Est" ) #%>%
x.FDis <- x.FDis[-c(1:24000),]
x.FDis$metric <- "FDis"

x.FDisUW <- model.FDisUW %>%
  spread_draws(b_ndvi)%>%
  gather(key = "b_", value = "Est" ) #%>%
x.FDisUW <- x.FDisUW[-c(1:24000),]
x.FDisUW$metric <- "FDisUW"

x.FDisSPUW <- model.FDisSPUW %>%
  spread_draws(b_ndvi)%>%
  gather(key = "b_", value = "Est" ) #%>%
x.FDisSPUW <- x.FDisSPUW[-c(1:24000),]
x.FDisSPUW$metric <- "FDisSPUW"

x.FDis.Comp <- rbind(x.FDis,x.FDisUW,x.FDisSPUW)

x.FDis.Comp$metric <- factor(x.FDis.Comp$metric, 
                             levels = c("FDisSPUW", "FDisUW",
                                        "FDis" ))

##Plot for Figure 3
ggplot(data = x.FDis.Comp, aes(y = metric, x = Est, fill = metric)) +
  stat_halfeye()+
  geom_vline(xintercept = 0)+
  labs(title="Effect of NDVI on different functional dispersion metrics",
       x ="Beta effect size", y = "Metric")+
  scale_y_discrete(labels = c("FDis" = "Occupancy weighted functional dispersion",
                              "FDisUW" = "Realized community unweighted functional dispersion",
                              "FDisSPUW" = "Species pool unweighted functional dispersion"))+
  scale_fill_manual(values = c( "grey85","grey85","grey35"))+
  theme_classic()+
  theme(legend.position = "none")






###############Figures######################
##Figure 1###
##Plots of first three dimmensions of FRich for all sites 
library(tidyverse)
library(FactoMineR)
library(geometry)
library(rgl)
library(FD)

spp.pool.traits <- read.csv("Full trait list bi 1 kg.csv", row.names = 1)
spp.pool.traits$Body.Mass <- log(spp.pool.traits$Body.Mass)
spp.pool.traits$Litter.Size <- log(spp.pool.traits$Litter.Size)
x.dist <- gowdis(spp.pool.traits, w = c(15,3,3,3,3,3,15,15,5,5,5,15))
x.dist2 <- sqrt(x.dist)
x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
#Spp.PCA <- x.pco$li[,1:2]
Spp.PCA <- x.pco$li[,1:3]
Spp.PCA <- as.data.frame(Spp.PCA)
Spp.PCA$Species <- rownames(Spp.PCA)

##All species plot##
site.conv <- convhulln(Spp.PCA[1:3], output.options = TRUE)
rgl.open()
rgl.bg(color = "white")
plot(site.conv, color = "grey", alpha = 0.2)
rgl.points(Spp.PCA[,1],Spp.PCA[,2],Spp.PCA[,3], color = "black", size = 2)
rgl.points(-0.004688967,	0.426497406,	-0.07075241, color = "green", size = 10) #elephant
rgl.points(-0.412017842,	0.185362652,	0.20420967, color = "orange", size = 10) #chimpanzee
rgl.points(0.248052810,	-0.037742627,	-0.15164242, color = "red", size = 10) #tiger
rgl.points(0.119801317,	0.098036581,	0.24726129, color = "blue", size = 10) #Cephalophus ogilbyi
rgl.points(-0.094972307,	-0.36533246,	-0.14583557, color = "pink", size = 10) #Didelphis marsupialis
rgl.lines(c(-.5, .5), c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(-.5,.5), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(0, 0), c(-.5,.5), color = "black")
#rgl.postscript("All_Species.pdf",fmt="pdf")
rgl.snapshot("All_Species.png", fmt ="png")



######
###Figure 2####

library(tidyverse)
library(FactoMineR)
library(geometry)
library(rgl)
library(FD)
library(tidyr)
library(dplyr)
library(vegan)
##Functional richness 3D site-level plots##
#Call communtiy lists
BBS.Community.traits <- read.csv("Site Trait and Occ data/BBS/BBS Community traits.csv")
BCI.Community.traits <- read.csv("Site Trait and Occ data/BCI/BCI Community traits.csv")
BIF.Community.traits <- read.csv("Site Trait and Occ data/BIF/BIF Community traits.csv")
CAX.Community.traits <- read.csv("Site Trait and Occ data/CAX/CAX Community traits.csv")
COU.Community.traits <- read.csv("Site Trait and Occ data/COU/COU Community traits.csv")
CSN.Community.traits <- read.csv("Site Trait and Occ data/CSN/CSN Community traits.csv")
KRP.Community.traits <- read.csv("Site Trait and Occ data/KRP/KRP Community traits.csv")
NAK.Community.traits <- read.csv("Site Trait and Occ data/NAK/NAK Community traits.csv")
NNN.Community.traits <- read.csv("Site Trait and Occ data/NNN/NNN Community traits.csv")
PSH.Community.traits <- read.csv("Site Trait and Occ data/PSH/PSH Community traits.csv")
RNF.Community.traits <- read.csv("Site Trait and Occ data/RNF/RNF Community traits.csv")
UDZ.Community.traits <- read.csv("Site Trait and Occ data/UDZ/UDZ Community traits.csv")
VB.Community.traits <- read.csv("Site Trait and Occ data/VB/VB Community traits.csv")
YAN.Community.traits <- read.csv("Site Trait and Occ data/YAN/YAN Community traits.csv")
YAS.Community.traits <- read.csv("Site Trait and Occ data/YAS/YAS Community traits.csv")


spp.pool.traits <- read.csv("Full trait list bi 1 kg.csv", row.names = 1)
spp.pool.traits$Body.Mass <- log(spp.pool.traits$Body.Mass)
spp.pool.traits$Litter.Size <- log(spp.pool.traits$Litter.Size)
x.dist <- gowdis(spp.pool.traits, w = c(15,3,3,3,3,3,15,15,5,5,5,15))
x.dist2 <- sqrt(x.dist)
x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
#Spp.PCA <- x.pco$li[,1:2]
Spp.PCA <- x.pco$li[,1:3]
Spp.PCA <- as.data.frame(Spp.PCA)
Spp.PCA$Species <- rownames(Spp.PCA)
Site.PCA <- filter(Spp.PCA, row.names(Spp.PCA) %in% YAS.Community.traits$Species)

site.conv <- convhulln(Site.PCA[1:3], output.options = TRUE)
rgl.open()
rgl.bg(color = "white")
plot(site.conv, color = "#FF6600", alpha = 0.35) ##Find color in Hull.Color.csv 
rgl.points(Site.PCA[,1],Site.PCA[,2],Site.PCA[,3], color = "black", size = 3)

##Extirpated species
##BCI
#rgl.points(-0.041552291,	0.35048043,	-0.10586191, color = "grey45", size = 10) #White-lipped peccary
#rgl.points(0.221552697,	0.04363491,	0.03526927, color = "grey45", size = 10) #Giant anteater
#rgl.points(-0.100327721,	-0.220308258,	0.02300786, color = "grey45", size = 10) #Jaguar
##BIF
#rgl.points(-0.054550847,	-0.261151009,	-0.05237290, color = "grey45", size = 10) #Leopard
#rgl.points(-0.015123127,	0.37834339,	-0.11112738, color = "grey45", size = 10) #Buffalo
#rgl.points(-0.015686861,	0.29048518,	-0.16486485, color = "grey45", size = 10) #Giant forest hog
##KRP
#rgl.points(-0.054550847,	-0.261151009,	-0.05237290, color = "grey45", size = 10) #Leopard
#rgl.points(-0.062212154,	-0.284508503,	-0.06660082, color = "grey45", size = 10) #African Golden Cat
#rgl.points(0.222536297,	-0.03705216,	0.07619492, color = "grey45", size = 10) #Giant pangolin
##NAK
#rgl.points(-0.054550847,	-0.261151009,	-0.05237290, color = "grey45", size = 10) #Leopard
#rgl.points(0.248052810,	-0.037742627,	-0.15164242, color = "grey45", size = 10) #Tiger
rgl.snapshot("YAS.ext.Rich.png",fmt="png")

##Functional dispersion site-level plots##
#Call communtiy lists
BBS.Community.traits <- read.csv("Site Trait and Occ data/BBS/BBS Community traits.csv")
BCI.Community.traits <- read.csv("Site Trait and Occ data/BCI/BCI Community traits.csv")
BIF.Community.traits <- read.csv("Site Trait and Occ data/BIF/BIF Community traits.csv")
CAX.Community.traits <- read.csv("Site Trait and Occ data/CAX/CAX Community traits.csv")
COU.Community.traits <- read.csv("Site Trait and Occ data/COU/COU Community traits.csv")
CSN.Community.traits <- read.csv("Site Trait and Occ data/CSN/CSN Community traits.csv")
KRP.Community.traits <- read.csv("Site Trait and Occ data/KRP/KRP Community traits.csv")
NAK.Community.traits <- read.csv("Site Trait and Occ data/NAK/NAK Community traits.csv")
NNN.Community.traits <- read.csv("Site Trait and Occ data/NNN/NNN Community traits.csv")
PSH.Community.traits <- read.csv("Site Trait and Occ data/PSH/PSH Community traits.csv")
RNF.Community.traits <- read.csv("Site Trait and Occ data/RNF/RNF Community traits.csv")
UDZ.Community.traits <- read.csv("Site Trait and Occ data/UDZ/UDZ Community traits.csv")
VB.Community.traits <- read.csv("Site Trait and Occ data/VB/VB Community traits.csv")
YAN.Community.traits <- read.csv("Site Trait and Occ data/YAN/YAN Community traits.csv")
YAS.Community.traits <- read.csv("Site Trait and Occ data/YAS/YAS Community traits.csv")

spp.pool.traits <- read.csv("Full trait list bi 1 kg.csv", row.names = 1)
spp.pool.traits$Body.Mass <- log(spp.pool.traits$Body.Mass)
spp.pool.traits$Litter.Size <- log(spp.pool.traits$Litter.Size)
x.dist <- gowdis(spp.pool.traits, w = c(15,3,3,3,3,3,15,15,5,5,5,15))
x.dist2 <- sqrt(x.dist)
x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
#Spp.PCA <- x.pco$li[,1:2]
Spp.PCA <- x.pco$li[,1:3]
Spp.PCA <- as.data.frame(Spp.PCA)
Spp.PCA$Species <- rownames(Spp.PCA)
Site.PCA <- filter(Spp.PCA, 
      row.names(Spp.PCA) %in% YAS.Community.traits$Species) ###

##List for folders of TEAM sites
file.list = list.files("./Site Trait and Occ data")
Species <- matrix(c("Species","Site"), ncol = 2)
for (i in 1:length(file.list)){
  fgdb = paste0("Site Trait and Occ data/", file.list[i],"/",file.list[i]," spp pool.csv")
  spp.pool <- read.csv(fgdb, row.names=1)
  spp <- as.matrix(rownames(spp.pool))
  spp <- cbind(spp, file.list[i])
  Species <- rbind(Species,spp)
}  
Species <- as.data.frame(Species)  
colnames(Species) <- c("Species","Site")
Species <- Species[-1,]
All.TEAM.Occupancy <- read.csv("All TEAM Occupancy.csv")
All.TEAM.Occupancy$means <- rowMeans(All.TEAM.Occupancy[,3:10], na.rm = TRUE)
for (i in 1:nrow(All.TEAM.Occupancy)){
ifelse(is.na(All.TEAM.Occupancy[i,10]),ifelse(
  is.na(All.TEAM.Occupancy[i,9]),
  All.TEAM.Occupancy[i,12] <- All.TEAM.Occupancy[i,8], 
  All.TEAM.Occupancy[i,12] <- All.TEAM.Occupancy[i,9]),
  All.TEAM.Occupancy[i,12] <- All.TEAM.Occupancy[i,10])
}
Occ <- All.TEAM.Occupancy[,c(1,2,12)]
#Occ$means <- 1
Occ$Year <- "Comm Pres"
colnames(Occ) <- c("Site","Species","Occupancy","Year")
Site.Occ <- filter(Occ, Occ$Site == "YAS") ###
Site.PCA <- merge(Site.PCA,Site.Occ, by ="Species")
Norm.Occ <- Site.PCA[,c(1,6)]
rownames(Norm.Occ) <- Norm.Occ$Species
Norm.Occ$Species <- NULL
Norm.Occ <- t(Norm.Occ)
deco<- decostand(Norm.Occ, method = "total")
deco <- t(deco)
Site.PCA$Occupancy <- deco
Site.PCA <- mutate(Site.PCA, Occupancy = 100*Occupancy)

##FDis 3D plots##

site.conv <- convhulln(Site.PCA[2:4], output.options = TRUE)
rgl.open()
rgl.bg(color = "white")
plot(site.conv, color = "#345ECE", alpha = 0.35) #####Find color in Hull.Color.csv 
for (i in 1:nrow(Site.PCA)){
rgl.points(Site.PCA[i,2],Site.PCA[i,3],Site.PCA[i,4], color = "black", size = Site.PCA[i,6] )
           }
#rgl.lines(c(-10, 10), c(0, 0), c(0, 0), color = "black")
#rgl.lines(c(0, 0), c(-10,10), c(0, 0), color = "red")
#rgl.lines(c(0, 0), c(0, 0), c(-10,10), color = "green")
rgl.snapshot("YAS_FDis.png",fmt="png") ###

###Other aspects of map figures
##Making map
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_coastline(scale = "medium", returnclass = "sf")


Site.Coords <- read.csv("Site Coords.csv")
Site.Coor_sf <- Site.Coords %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

ggplot(data = world) +
  geom_sf(size = 0.05) +
  theme_classic() +
  geom_sf(data = Site.Coor_sf, color = 'black', shape = 8, size = 2)

##Making color ramp for maps
colfunc <- colorRampPalette(c('yellow','orange','red'))
colfunc(165)
plot(rep(1,165),col=colfunc(165),pch=19,cex=3,yaxt='n', ann=FALSE)


colfunc <- colorRampPalette(c('yellow','green','dodgerblue','purple4'))
colfunc(73)
plot(rep(1,73),col=colfunc(73),pch=19,cex=3,yaxt='n', ann=FALSE)


####
##Figure 3
##See above

###
##Figure 1S
##Compare trait weighting vs non-weighting

TEAM.Site.Funct.bi <- read.csv("TEAM Site Funct bi.csv")
dat <- TEAM.Site.Funct.bi


##Standardize Predictor Variables
dat$ndvi <- scale(dat$ndvi)
dat$Landcover.div <- scale(dat$Landcover.div)
dat$hum.dens.ZOI <- scale(dat$hum.dens.ZOI)
dat$FPoolRich <- scale(dat$FPoolRich)
dat$Fpool.uw <- scale(dat$Fpool.uw)
dat$FR.Lost <- scale(dat$FR.Lost)
dat$FR.Lost.uw <- scale(dat$FR.Lost.uw)

library(rstan)
library(brms)
library(ggplot2)
library(tidyr)
library(ggdist)
library(tidybayes)

##Functional richness

model.FR1 <- brm(formula = Frich ~ FPoolRich + Landcover.div + ndvi + hum.dens.ZOI + FR.Lost + Africa + Asia +Madagascar, 
                 data = dat, control = list(adapt_delta = 0.99, max_treedepth = 15), family = "weibull", iter = 4000)

model.FR1.uw <- brm(formula = Frich.uw ~ Fpool.uw + Landcover.div + ndvi + hum.dens.ZOI + FR.Lost.uw + Africa + Asia + Madagascar, 
                    data = dat, control = list(adapt_delta = 0.99, max_treedepth = 15), family = "weibull", iter = 4000)

summary(model.FR1)
summary(model.FR1.uw)


##Plot Credible Intervals for Functional Richness
x <- model.FR1 %>%
  spread_draws(b_Intercept, b_FPoolRich, b_Landcover.div, b_ndvi, b_hum.dens.ZOI,b_FR.Lost, b_Africa, b_Asia, b_Madagascar)%>%
  gather(key = "b_", value = "Est" ) #%>%
x <- x[-c(1:32000),]
x$b_ <- factor(x$b_, levels = c("b_ndvi", "b_Landcover.div",
                                "b_hum.dens.ZOI","b_FR.Lost",
                                "b_FPoolRich","b_Africa",
                                "b_Asia","b_Madagascar"))
ggplot(data = x, aes(y = b_, x = Est, fill = b_)) +
  stat_halfeye()+
  geom_vline(xintercept = 0)+
  labs(title="Functional Richness model",
       x ="Beta effect size", y = "Predictor variables")+
  scale_y_discrete(labels = c("b_ndvi" = "NDVI", "b_Landcover.div" = "Landcover Diversity",
                              "b_hum.dens.ZOI" = "Human Density","b_FR.Lost"="Extinction",
                              "b_FPoolRich"="Species Pools","b_Africa"="Africa",
                              "b_Asia"="Asia","b_Madagascar"="Madagascar"))+
  scale_fill_manual(values = c("grey85","grey85","grey85","red","grey85","grey85","grey85","red"))+
  theme_classic()+
  theme(legend.position = "none")


x <- model.FR1.uw %>%
  spread_draws(b_Intercept, b_Fpool.uw, b_Landcover.div, b_ndvi, b_hum.dens.ZOI,b_FR.Lost.uw, b_Africa, b_Asia, b_Madagascar)%>%
  gather(key = "b_", value = "Est" ) #%>%
x <- x[-c(1:32000),]
x$b_ <- factor(x$b_, levels = c("b_ndvi", "b_Landcover.div",
                                "b_hum.dens.ZOI","b_FR.Lost.uw",
                                "b_Fpool.uw","b_Africa",
                                "b_Asia","b_Madagascar"))
ggplot(data = x, aes(y = b_, x = Est, fill = b_)) +
  stat_halfeye()+
  geom_vline(xintercept = 0)+
  labs(title="Functional Richness model (unweighted traits)",
       x ="Beta effect size", y = "Predictor variables")+
  scale_y_discrete(labels = c("b_ndvi" = "NDVI", "b_Landcover.div" = "Landcover Diversity",
                              "b_hum.dens.ZOI" = "Human Density","b_FR.Lost.uw"="Extinction",
                              "b_Fpool.uw"="Species Pools","b_Africa"="Africa",
                              "b_Asia"="Asia","b_Madagascar"="Madagascar")) +
  scale_fill_manual(values = c("grey85","grey85","grey85","red","grey85","grey85","grey85","red"))+
  theme_classic()+
  theme(legend.position = "none")


##Functional Dispersion

model.FDis <- brm(formula = Fdis ~ FPoolRich + Landcover.div + ndvi + FR.Lost + hum.dens.ZOI + Africa + Asia + Madagascar, 
                  data = dat,  control = list(adapt_delta = 0.99, max_treedepth = 20), family = "weibull", iter = 4000)
summary(model.FDis)

model.FDis.uw <- brm(formula = Fdis.uw ~ Fpool.uw + Landcover.div + ndvi + hum.dens.ZOI + FR.Lost.uw + Africa + Asia + Madagascar, 
                     data    = dat,  control = list(adapt_delta = 0.99, max_treedepth = 20), family = "weibull", iter = 4000)
summary(model.FDis.uw)

x <- model.FDis %>%
  spread_draws(b_Intercept, b_FPoolRich, b_Landcover.div, b_ndvi, b_hum.dens.ZOI,b_FR.Lost, b_Africa, b_Asia, b_Madagascar)%>%
  gather(key = "b_", value = "Est" ) #%>%
x <- x[-c(1:32000),]
x$b_ <- factor(x$b_, levels = c("b_ndvi", "b_Landcover.div",
                                "b_hum.dens.ZOI","b_FR.Lost",
                                "b_FPoolRich","b_Africa",
                                "b_Asia","b_Madagascar"))
ggplot(data = x, aes(y = b_, x = Est, fill = b_)) +
  stat_halfeye()+
  geom_vline(xintercept = 0)+
  labs(title="Functional Dispersion model",
       x ="Beta", y = "Predictor variables")+
  scale_y_discrete(labels = c("b_ndvi" = "NDVI", "b_Landcover.div" = "Landcover Diversity",
                              "b_hum.dens.ZOI" = "Human Density","b_FR.Lost"="Extinction",
                              "b_FPoolRich"="Species Pools","b_Africa"="Africa",
                              "b_Asia"="Asia","b_Madagascar"="Madagascar"))+
  scale_fill_manual(values = c("red", "grey85","grey85","grey85","grey85","grey85","grey85","grey85"))+
  theme_classic()+
  theme(legend.position = "none")

x <- model.FDis.uw %>%
  spread_draws(b_Intercept, b_Fpool.uw, b_Landcover.div, b_ndvi, b_hum.dens.ZOI,b_FR.Lost.uw, b_Africa, b_Asia, b_Madagascar)%>%
  gather(key = "b_", value = "Est" ) #%>%
x <- x[-c(1:32000),]
x$b_ <- factor(x$b_, levels = c("b_ndvi", "b_Landcover.div",
                                "b_hum.dens.ZOI","b_FR.Lost.uw",
                                "b_Fpool.uw","b_Africa",
                                "b_Asia","b_Madagascar"))
ggplot(data = x, aes(y = b_, x = Est, fill = b_)) +
  stat_halfeye()+
  geom_vline(xintercept = 0)+
  labs(title="Functional Dispersion model (unweighted traits)",
       x ="Beta", y = "Predictor variables")+
  scale_y_discrete(labels = c("b_ndvi" = "NDVI", "b_Landcover.div" = "Landcover Diversity",
                              "b_hum.dens.ZOI" = "Human Density","b_FR.Lost.uw"="Extinction",
                              "b_Fpool.uw"="Species Pools","b_Africa"="Africa",
                              "b_Asia"="Asia","b_Madagascar"="Madagascar"))+
  scale_fill_manual(values = c("red","grey85","grey85","red","grey85","grey85","grey85","red"))+
  theme_classic()+
  theme(legend.position = "none")
