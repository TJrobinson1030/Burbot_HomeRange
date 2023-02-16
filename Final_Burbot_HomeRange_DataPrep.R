# MONTHLY HOME RANGE CALCULATIONS
rm(list=ls(all=TRUE))

# NEED THESE PACKAGES TO RUN ALL OF THE CODE IN THIS SCRIPT
library(readr)
library(adehabitatHR)
library(tidyverse)
library(anytime)
library(PBSmapping)
library(lubridate)
library(RchivalTag)
library(akima)

#
citation("adehabitatHR")

#Read Data 
BurbLocsHR <- read_csv("E:/BSUGradBurbot/DATA/FishData/FinalFilteredData/FinalBurbLocsHR_CLIPPED.csv")

#refomat text date from depth data as a factor to POSIXct
BurbLocsHR$datetime = anytime(as.factor(BurbLocsHR$datetime))
str(BurbLocsHR)
tz(BurbLocsHR$datetime)

BurbLocsHR$sunrise = anytime(as.factor(BurbLocsHR$sunrise))
tz(BurbLocsHR$sunrise)

BurbLocsHR$sunset = anytime(as.factor(BurbLocsHR$sunset))
tz(BurbLocsHR$sunset)

BurbLocsHR$UTC = anytime(as.factor(BurbLocsHR$UTC))
tz(BurbLocsHR$UTC)

#force timezone to UTC
BurbLocsHR <- BurbLocsHR %>%
  mutate(UTC = force_tz(UTC, "UTC"))
tz(BurbLocsHR$UTC)

#rm(BurbLocsHR)
setwd("E:/BSUGradBurbot/DATA/FishData/FinalFilteredData/HomeRange_AllSeasons")

SUM_TAG =
  BurbLocsHR %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#Start of data used for hot spot analysis
BurbLocsHR1 <- subset(BurbLocsHR, DEPTH != "NA")

BurbLocsHR_NegDepth <- BurbLocsHR1[BurbLocsHR1$DEPTH >= 0,]

BurbLoc_Dep <- sample_n(BurbLocsHR_NegDepth, 20)

#write.csv(BurbLocsHR_NegDepth, "BurbLoc_Dep_all.csv")

# see how many locations per fish during time period
BurbLocsHR_depthcount <- BurbLocsHR1 %>% group_by(TRANSMITTER) %>% summarise(n = n())
BurbLocsHR_negdepthcount <- BurbLocsHR_NegDepth %>% group_by(TRANSMITTER) %>% summarise(n = n())

depthcount <- left_join(BurbLocsHR_depthcount,BurbLocsHR_negdepthcount, by = "TRANSMITTER")

################################################################################
######### BEGIN ESTIMATING HOME RANGE ############
################################################################################

################################################################################
# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- BurbLocsHR[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(BurbLocsHR[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb <- SpatialPointsDataFrame(coords      = coords,
                               data        = data, 
                               proj4string = crs)
#rm(burb)
#create grid to calculate volume and area of home range
#https://stackoverflow.com/questions/41683905/grid-too-small-for-kernelud-getverticeshr-adehabitathr-home-range-estimation
xmin <- 314450.03
xmax <- 321589.03
ymin <- 5217400.97
ymax <- 5225680.97

x <- seq(xmin, xmax, by = 10.)
y <- seq(ymin, ymax, by = 10.)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)

#calculate kernel density estimate for all fish
kud <- kernelUD(burb[,1], h="href", grid = xy)
#rm(kud)

#print all fish estimates
image(kud)

#vectorized home range estimate kud
homerange <- getverticeshr(kud, unin = "m",unout = "km2")
class(homerange)
plot(homerange,col = 1:58)
#rm(homerange)

#rasterized home range estimate 95 kud
vud <- getvolumeUD(kud)
vud
#rm(vud)

## Set up graphical parameters
par(mfrow=c(1,2))
par(mar=c(0,0,2,0))

## The output of kernelUD for the animal
image(kud[[5]])
title("Output of kernelUD")

## Convert into a suitable data structure for the use of contour
xyz <- as.image.SpatialGridDataFrame(kud[[5]])
contour(xyz, add=TRUE)

## and similarly for the output of getvolumeUD
par(mar=c(0,0,2,0))
image(vud[[5]])
title("Output of getvolumeUD")
xyzv <- as.image.SpatialGridDataFrame(vud[[5]])
contour(xyzv, add=TRUE)

## store the volume under the UD (as computed by getvolumeUD)
## of the animal in fud
fud <- vud[[5]]

## store the value of the volume under UD in a vector hr90
hr90 <- as.data.frame(fud)[,1]
## if hr90 is <= 90 then the pixel belongs to the home range
## (takes the value 1, 0 otherwise)
hr90 <- as.numeric(hr90 <= 90)

## Converts into a data frame
hr90 <- data.frame(hr90)

## Converts to a SpatialPixelsDataFrame
coordinates(hr90) <- coordinates(vud[[1]])
gridded(hr90) <- TRUE

## display the results
image(hr90)
#rm(fud)

#write 95 kud as data frame
HR <- as.data.frame(homerange)
#rm(HR)

#calculate kernel area estimates
N <- kernel.area(kud, percent = seq(50,95, by = 5), unout = "km2")
N <- as.data.frame(t(N))
colnames(N) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs <- cbind(fish,N)
HR_Burbs <- HR_Burbs %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN, HR_50,HR_75,HR_80,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs,"BurbHR_AllStudy.csv")

################################################################################
#get vertices for each area estimate
burb.all.95.kde <- getverticeshr(kud, percent = 95, unin = "m", unout = "km2")
burb.all.90.kde <- getverticeshr(kud, percent = 90, unin = "m", unout = "km2")
burb.all.85.kde <- getverticeshr(kud, percent = 85, unin = "m", unout = "km2")
burb.all.80.kde <- getverticeshr(kud, percent = 80, unin = "m", unout = "km2")
burb.all.75.kde <- getverticeshr(kud, percent = 75, unin = "m", unout = "km2")
burb.all.70.kde <- getverticeshr(kud, percent = 70, unin = "m", unout = "km2")
burb.all.65.kde <- getverticeshr(kud, percent = 65, unin = "m", unout = "km2")
burb.all.60.kde <- getverticeshr(kud, percent = 60, unin = "m", unout = "km2")
burb.all.55.kde <- getverticeshr(kud, percent = 55, unin = "m", unout = "km2")
burb.all.50.kde <- getverticeshr(kud, percent = 50, unin = "m", unout = "km2")
#rm(burb.all.95.kde)

#figure out what transmitters were alive
TRANSMITTER <- row.names(burb.all.95.kde)

#page that provided below function to clip away land from HR estimates
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
#To remove land from your 95% kde, first read in the shapefile of the U
#S (e.g"NOS80k.shp") using the readOGR function:
#If you have previously installed these packages, load the libraries:

library(adehabitatHR)
library(rgdal)
library(raster)
library(rgeos)  

######################################################################
  #https://gis.stackexchange.com/questions/64537/clip-polygon-and-retain-data
us <- readOGR("E:/BSUGradBurbot/DATA/FishData/FinalFilteredData/ShapeFiles/LakeAreaE.shp")

#Transform/reproject to WGS84, so that it has the same CRS as everything else:
us.pro <- spTransform(us,CRS("+init=epsg:32615"))

##########
# create function that is used to remove land from hr measurements
HR_CLIP <- function(data){
  #https://stat.ethz.ch/pipermail/r-sig-geo/2020-October/028458.html
  proj4string(data) <- proj4string(us.pro)
  data <- spTransform(data, proj4string(us.pro))
  data <- spTransform(x = data, CRSobj = "+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs"); us.pro<- spTransform(x = us.pro, CRSobj = "+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs")
  
  #Remove land using the gdifference function:
  burb.all.kde1 <- gDifference(data, us.pro, byid=TRUE,drop_lower_td=TRUE)
  
  #Turn the new kde into a spatial polygon data frame, so that it can be outputted as a shapefile:
  burb.all.kde1 <- as(burb.all.kde1, "SpatialPolygonsDataFrame")
  
  #get all fishes areas for HR
  area_Km2 = sapply(slot(burb.all.kde1, "polygons"), slot, "area")
  
  #combine trans and area estimats
  HR = cbind.data.frame(TRANSMITTER,area_Km2)
  
  #convert to km2
  HR$area_Km2 <- (HR$area_Km2) / 1000000
}

#run land clip for each kde  
HR95 <- HR_CLIP(burb.all.95.kde)
HR90 <- HR_CLIP(burb.all.90.kde)
HR85 <- HR_CLIP(burb.all.85.kde)
HR80 <- HR_CLIP(burb.all.85.kde)
HR75 <- HR_CLIP(burb.all.85.kde)
HR70 <- HR_CLIP(burb.all.85.kde)
HR65 <- HR_CLIP(burb.all.85.kde)
HR60 <- HR_CLIP(burb.all.85.kde)
HR55 <- HR_CLIP(burb.all.85.kde)
HR50 <- HR_CLIP(burb.all.50.kde)
#rm(HR50)

#combine each kde estimate without land
study_HR <- as.data.frame(cbind(TRANSMITTER,HR95,HR90,HR85,HR80,HR75,HR70,HR65,HR60,HR55,HR50))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR <- merge(fish,study_HR, by = "TRANSMITTER")
study_HR1 <- study_HR %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50,HR85,HR90,HR95)

#write.csv(study_HR1,"BurbHR_TRIMMED.csv")

##################
#All above code is how you can get hr estimate across the study period for each fish
#code takes a very long time to run and needs to be ran one line at a time
#below is code that is used for season and below that is for each 25 day period
##################

################################################################################
############# ICE ONE - ICE OFF SEASON #########################################
################################################################################

Ice_Season <- BurbLocsHR

Ice19 <- Ice_Season[Ice_Season$datetime >= "2019-04-10" & Ice_Season$datetime <= "2019-04-27",]
Ice20 <- Ice_Season[Ice_Season$datetime >= "2019-12-03" & Ice_Season$datetime <= "2020-04-29",]
Open19 <- Ice_Season[Ice_Season$datetime >= "2019-04-27" & Ice_Season$datetime <= "2019-12-03",]
Open20 <- Ice_Season[Ice_Season$datetime >= "2020-04-29" & Ice_Season$datetime <= "2020-07-29",]

Ice <-rbind(Ice19,Ice20)
Open <-rbind(Open19,Open20) 
#rm(Open)

################################################################################
#################### ICE ICE ICE ###############################################
################################################################################
# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Ice[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Ice[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_ice <- SpatialPointsDataFrame(coords      = coords,
                               data        = data, 
                               proj4string = crs)
#rm(burb_ice)

#calculate kernel density estimate for all fish during ice period
kud_ice <- kernelUD(burb_ice[,1], h="href", grid = xy)
#rm(kud_ice)

#print all fish estimates
image(kud_ice)

#calculate kernel area estimates
N_ice <- kernel.area(kud, percent = seq(50,95, by = 5), unout = "km2")
N_ice <- as.data.frame(t(N_ice))
colnames(N_ice) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_ice)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Ice <- cbind(fish,N_ice)
HR_Burbs_Ice <- HR_Burbs_Ice %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN, HR_50,HR_75,HR_80,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Ice,"BurbHR_Ice_Study.csv")
################################################################################
burb.ice.95.kde <- getverticeshr(kud_ice, percent = 95, unin = "m", unout = "km2")
burb.ice.90.kde <- getverticeshr(kud_ice, percent = 90, unin = "m", unout = "km2")
burb.ice.85.kde <- getverticeshr(kud_ice, percent = 85, unin = "m", unout = "km2")
burb.ice.50.kde <- getverticeshr(kud_ice, percent = 50, unin = "m", unout = "km2")
#rm(burb.ice.90.kde)

# hr clipped to lake area
HR95_ice <- HR_CLIP(burb.ice.95.kde)
HR90_ice <- HR_CLIP(burb.ice.90.kde)
HR85_ice <- HR_CLIP(burb.ice.85.kde)
HR50_ice <- HR_CLIP(burb.ice.50.kde)

study_HR_ice <- as.data.frame(cbind(TRANSMITTER,HR95_ice,HR90_ice,HR85_ice,HR50_ice))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_ice <- merge(fish,study_HR_ice, by = "TRANSMITTER")
study_HR_ice <- study_HR_ice %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_ice,HR85_ice,HR90_ice,HR95_ice)

#write.csv(study_HR_ice,"BurbHR_Ice_TRIMMED.csv")

#############################################
########## OPEN WATER #######################
#############################################

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Open[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Open[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_open <- SpatialPointsDataFrame(coords      = coords,
                                   data        = data, 
                                   proj4string = crs)
#rm(burb_open)

#calculate kernel density estimate for all fish during open water
kud_open <- kernelUD(burb_open[,1], h="href", grid = xy)
#rm(kud_open)

#print all fish estimates
image(kud_open)

#calculate kernel area estimates
N_open <- kernel.area(kud_open, percent = seq(50,95, by = 5), unout = "km2")
N_open <- as.data.frame(t(N_open))
colnames(N_open) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_open)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_open <- cbind(fish,N_open)
HR_Burbs_open <- HR_Burbs_open %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN, HR_50,HR_75,HR_80,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_open,"BurbHR_Open_Study.csv")

################################################################################
burb.open.95.kde <- getverticeshr(kud_open, percent = 95, unin = "m", unout = "km2")
burb.open.90.kde <- getverticeshr(kud_open, percent = 90, unin = "m", unout = "km2")
burb.open.85.kde <- getverticeshr(kud_open, percent = 85, unin = "m", unout = "km2")
burb.open.50.kde <- getverticeshr(kud_open, percent = 50, unin = "m", unout = "km2")
#rm(burb.open.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.open.95.kde)

HR95_open <- HR_CLIP(burb.open.95.kde)
HR90_open <- HR_CLIP(burb.open.90.kde)
HR85_open <- HR_CLIP(burb.open.85.kde)
HR50_open <- HR_CLIP(burb.open.50.kde)

study_HR_open <- as.data.frame(cbind(TRANSMITTER,HR95_open,HR90_open,HR85_open,HR50_open))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_open <- merge(fish,study_HR_open, by = "TRANSMITTER")
study_HR_open <- study_HR_open %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_open,HR85_open,HR90_open,HR95_open)

#write.csv(study_HR_open,"BurbHR_Open_TRIMMED.csv")


################################################################################
############### 25 day span home range across study period #####################
################################################################################
#make new file for date range hr estimates
HR_25day <- BurbLocsHR
HR_25day$datetime = anytime(as.factor(HR_25day$datetime))
str(HR_25day)

# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE1 <- as.Date("2019-04-25")

DATE1 - 12
DATE1 + 12

DATE2 <- as.Date("2019-04-13")
DATE3 <- as.Date("2019-05-07")

Apr25_19 <- myfunc(DATE2,DATE3)    

################################
# April 2019
################################
Apr2019_sum_tag =
  Apr25_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Apr25_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Apr25_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Apr25_19 <- SpatialPointsDataFrame(coords      = coords,
                                   data        = data, 
                                   proj4string = crs)
#rm(burb_Apr25_19)

#calculate kernel density estimate for all fish
kud_Apr25_19 <- kernelUD(burb_Apr25_19[,1], h="href", grid = xy)
#rm(kud_Apr25_19)

#print all fish estimates
#image(kud_Apr25_19)

#calculate kernel area estimates
N_Apr25_19 <- kernel.area(kud_Apr25_19, percent = seq(50,95, by = 5), unout = "km2")
N_Apr25_19 <- as.data.frame(t(N_Apr25_19))
colnames(N_Apr25_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Apr25_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Apr25_19 <- cbind(fish,N_Apr25_19)
HR_Burbs_Apr25_19 <- HR_Burbs_Apr25_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Apr25_19,"BurbHR_Apr25_2019.csv")

################################################################################
burb.Apr25_19.95.kde <- getverticeshr(kud_Apr25_19, percent = 95, unin = "m", unout = "km2")
burb.Apr25_19.90.kde <- getverticeshr(kud_Apr25_19, percent = 90, unin = "m", unout = "km2")
burb.Apr25_19.85.kde <- getverticeshr(kud_Apr25_19, percent = 85, unin = "m", unout = "km2")
burb.Apr25_19.50.kde <- getverticeshr(kud_Apr25_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Apr25_19.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Apr25_19.95.kde)

HR95_Apr25_19 <- HR_CLIP(burb.Apr25_19.95.kde)
HR90_Apr25_19 <- HR_CLIP(burb.Apr25_19.90.kde)
HR85_Apr25_19 <- HR_CLIP(burb.Apr25_19.85.kde)
HR50_Apr25_19 <- HR_CLIP(burb.Apr25_19.50.kde)

#rm(HR95_Apr25_19)
#merge to data frame
study_HR_Apr25_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Apr25_19,HR90_Apr25_19,HR85_Apr25_19,HR50_Apr25_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Apr25_19 <- merge(fish,study_HR_Apr25_19, by = "TRANSMITTER")
study_HR_Apr25_19 <- study_HR_Apr25_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Apr25_19,HR85_Apr25_19,HR90_Apr25_19,HR95_Apr25_19)

#write.csv(study_HR_Apr25_19,"BurbHR_Apr25_2019_TRIMMED.csv")

################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_May <- as.Date("2019-05-25")

DATE_May - 12
DATE_May + 12

DATE_May1 <- as.Date("2019-05-13")
DATE_May2 <- as.Date("2019-06-06")

May25_19 <- myfunc(DATE_May1,DATE_May2)    

################################
# May 2019
################################
May2019_sum_tag =
  May25_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- May25_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(May25_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_May25_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_May25_19)

#calculate kernel density estimate for all fish
kud_May25_19 <- kernelUD(burb_May25_19[,1], h="href", grid = xy)
#rm(kud_May25_19)

#print all fish estimates
image(kud_May25_19)

#calculate kernel area estimates
N_May25_19 <- kernel.area(kud_May25_19, percent = seq(50,95, by = 5), unout = "km2")
N_May25_19 <- as.data.frame(t(N_May25_19))
colnames(N_May25_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_May25_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_May25_19 <- cbind(fish,N_May25_19)
HR_Burbs_May25_19 <- HR_Burbs_May25_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_May25_19,"BurbHR_May25_2019.csv")

################################################################################

#get area HR for each percentage
burb.May25_19.95.kde <- getverticeshr(kud_May25_19, percent = 95, unin = "m", unout = "km2")
burb.May25_19.90.kde <- getverticeshr(kud_May25_19, percent = 90, unin = "m", unout = "km2")
burb.May25_19.85.kde <- getverticeshr(kud_May25_19, percent = 85, unin = "m", unout = "km2")
burb.May25_19.50.kde <- getverticeshr(kud_May25_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.May25_19.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.May25_19.95.kde)

HR95_May25_19 <- HR_CLIP(burb.May25_19.95.kde)
HR90_May25_19 <- HR_CLIP(burb.May25_19.90.kde)
HR85_May25_19 <- HR_CLIP(burb.May25_19.85.kde)
HR50_May25_19 <- HR_CLIP(burb.May25_19.50.kde)

study_HR_May25_19 <- as.data.frame(cbind(TRANSMITTER,HR95_May25_19,HR90_May25_19,HR85_May25_19,HR50_May25_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_May25_19 <- merge(fish,study_HR_May25_19, by = "TRANSMITTER")
study_HR_May25_19 <- study_HR_May25_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_May25_19,HR85_May25_19,HR90_May25_19,HR95_May25_19)

#write.csv(study_HR_May25_19,"BurbHR_May25_2019_TRIMMED.csv")

################################################################################
################################################################################

################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Jun <- as.Date("2019-06-25")

DATE_Jun - 12
DATE_Jun + 12

DATE_Jun1 <- as.Date("2019-06-13")
DATE_Jun2 <- as.Date("2019-07-07")

Jun25_19 <- myfunc(DATE_Jun1,DATE_Jun2)    

################################
# Jun 2019
################################
Jun2019_sum_tag =
  Jun25_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish with less than 30 locations in the time period
Jun25_19 <- Jun25_19[!Jun25_19$TRANSMITTER == "492mm_F_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Jun25_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Jun25_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Jun25_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Jun25_19)
#calculate kernel density estimate for all fish
kud_Jun25_19 <- kernelUD(burb_Jun25_19[,1], h="href", grid = xy)
#rm(kud_Jun25_19)
#print all fish estimates
image(kud_Jun25_19)

#calculate kernel area estimates
N_Jun25_19 <- kernel.area(kud_Jun25_19, percent = seq(50,95, by = 5), unout = "km2")
N_Jun25_19 <- as.data.frame(t(N_Jun25_19))
colnames(N_Jun25_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Jun25_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fishJun <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Jun25_19 <- cbind(fishJun,N_Jun25_19)
HR_Burbs_Jun25_19 <- HR_Burbs_Jun25_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Jun25_19,"BurbHR_Jun25_2019.csv")
################################################################################
burb.Jun25_19.95.kde <- getverticeshr(kud_Jun25_19, percent = 95, unin = "m", unout = "km2")
burb.Jun25_19.90.kde <- getverticeshr(kud_Jun25_19, percent = 90, unin = "m", unout = "km2")
burb.Jun25_19.85.kde <- getverticeshr(kud_Jun25_19, percent = 85, unin = "m", unout = "km2")
burb.Jun25_19.50.kde <- getverticeshr(kud_Jun25_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Jun25_19.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Jun25_19.95.kde)

HR95_Jun25_19 <- HR_CLIP(burb.Jun25_19.95.kde)
HR90_Jun25_19 <- HR_CLIP(burb.Jun25_19.90.kde)
HR85_Jun25_19 <- HR_CLIP(burb.Jun25_19.85.kde)
HR50_Jun25_19 <- HR_CLIP(burb.Jun25_19.50.kde)

study_HR_Jun25_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Jun25_19,HR90_Jun25_19,HR85_Jun25_19,HR50_Jun25_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

fishJun <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]

#condense all data in to file with HR estimates and fish info
study_HR_Jun25_19 <- merge(fishJun,study_HR_Jun25_19, by = "TRANSMITTER")
study_HR_Jun25_19 <- study_HR_Jun25_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Jun25_19,HR85_Jun25_19,HR90_Jun25_19,HR95_Jun25_19)

#write.csv(study_HR_Jun25_19,"BurbHR_Jun25_2019_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Jul <- as.Date("2019-07-25")

DATE_Jul - 12
DATE_Jul + 12

DATE_Jul1 <- as.Date("2019-07-13")
DATE_Jul2 <- as.Date("2019-08-06")

Jul25_19 <- myfunc(DATE_Jul1,DATE_Jul2)    

################################
# Jul 2019
################################
Jul2019_sum_tag =
  Jul25_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish with less than 30 locations in the time period or causing failure in estimation
Jul25_19 <- Jul25_19[!Jul25_19$TRANSMITTER == "492mm_F_a_t",]
Jul25_19 <- Jul25_19[!Jul25_19$TRANSMITTER == "566mm_U_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Jul25_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Jul25_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Jul25_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Jul25_19)
#calculate kernel density estimate for all fish
kud_Jul25_19 <- kernelUD(burb_Jul25_19[,1], h="href", grid = xy)
#rm(kud_Jul25_19)
#print all fish estimates
image(kud_Jul25_19)

#calculate kernel area estimates
N_Jul25_19 <- kernel.area(kud_Jul25_19, percent = seq(50,95, by = 5), unout = "km2")
N_Jul25_19 <- as.data.frame(t(N_Jul25_19))
colnames(N_Jul25_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Jul25_19)
#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

fishJul <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]
fishJul <- fishJul[!fishJul$TRANSMITTER == "566mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Jul25_19 <- cbind(fishJul,N_Jul25_19)
HR_Burbs_Jul25_19 <- HR_Burbs_Jul25_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

HR_Burbs_Jul25_19_566 <- HR_Burbs_Jul25_19[HR_Burbs_Jul25_19$TRANSMITTER == "566mm_U_a_t",]

#write.csv(HR_Burbs_Jul25_19,"BurbHR_Jul25_2019.csv")
################################################################################
burb.Jul25_19.95.kde <- getverticeshr(kud_Jul25_19, percent = 95, unin = "m", unout = "km2")
burb.Jul25_19.90.kde <- getverticeshr(kud_Jul25_19, percent = 90, unin = "m", unout = "km2")
burb.Jul25_19.85.kde <- getverticeshr(kud_Jul25_19, percent = 85, unin = "m", unout = "km2")
burb.Jul25_19.50.kde <- getverticeshr(kud_Jul25_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Jul25_19.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Jul25_19.95.kde)

HR95_Jul25_19 <- HR_CLIP(burb.Jul25_19.95.kde)
HR90_Jul25_19 <- HR_CLIP(burb.Jul25_19.90.kde)
HR85_Jul25_19 <- HR_CLIP(burb.Jul25_19.85.kde)
HR50_Jul25_19 <- HR_CLIP(burb.Jul25_19.50.kde)
#rm(HR95_Jul25_19)

study_HR_Jul25_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Jul25_19,HR90_Jul25_19,HR85_Jul25_19,HR50_Jul25_19))

#condense all data in to file with HR estimates and fish info
study_HR_Jul25_19 <- merge(fishJul,study_HR_Jul25_19, by = "TRANSMITTER")
study_HR_Jul25_19 <- study_HR_Jul25_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Jul25_19,HR85_Jul25_19,HR90_Jul25_19,HR95_Jul25_19)

names(HR_Burbs_Jul25_19_566) <- c("TRANSMITTER", "SEX", "LENGTH", "WEIGHT", "BASIN", "HR50_Jul25_19", "HR85_Jul25_19", "HR90_Jul25_19","HR95_Jul25_19")

study_HR_Jul25_19_566 <- rbind(study_HR_Jul25_19,HR_Burbs_Jul25_19_566)

#write.csv(study_HR_Jul25_19_566,"BurbHR_Jul25_2019_TRIMMED.csv")

################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Aug <- as.Date("2019-08-25")

DATE_Aug - 12
DATE_Aug + 12

DATE_Aug1 <- as.Date("2019-08-13")
DATE_Aug2 <- as.Date("2019-09-06")

Aug25_19 <- myfunc(DATE_Aug1,DATE_Aug2)    

################################
# Aug 2019
################################
Aug2019_sum_tag =
  Aug25_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish with less than 30 locations in the time period or causing failure in estimation
Aug25_19 <- Aug25_19[!Aug25_19$TRANSMITTER == "566mm_U_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Aug25_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Aug25_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Aug25_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Aug25_19)
#calculate kernel density estimate for all fish
kud_Aug25_19 <- kernelUD(burb_Aug25_19[,1], h="href", grid = xy)
#rm(kud_Aug25_19)

#print all fish estimates
image(kud_Aug25_19)

#calculate kernel area estimates
N_Aug25_19 <- kernel.area(kud_Aug25_19, percent = seq(50,95, by = 5), unout = "km2")
N_Aug25_19 <- as.data.frame(t(N_Aug25_19))
colnames(N_Aug25_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Aug25_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#remove fish with less than 30 locations in the time period or causing failure in estimation
fish <- fish[!fish$TRANSMITTER == "566mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Aug25_19 <- cbind(fish,N_Aug25_19)
HR_Burbs_Aug25_19 <- HR_Burbs_Aug25_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Aug25_19,"BurbHR_Aug25_2019.csv")
################################################################################
burb.Aug25_19.95.kde <- getverticeshr(kud_Aug25_19, percent = 95, unin = "m", unout = "km2")
burb.Aug25_19.90.kde <- getverticeshr(kud_Aug25_19, percent = 90, unin = "m", unout = "km2")
burb.Aug25_19.85.kde <- getverticeshr(kud_Aug25_19, percent = 85, unin = "m", unout = "km2")
burb.Aug25_19.50.kde <- getverticeshr(kud_Aug25_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Aug25_19.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Aug25_19.95.kde)

HR95_Aug25_19 <- HR_CLIP(burb.Aug25_19.95.kde)
HR90_Aug25_19 <- HR_CLIP(burb.Aug25_19.90.kde)
HR85_Aug25_19 <- HR_CLIP(burb.Aug25_19.85.kde)
HR50_Aug25_19 <- HR_CLIP(burb.Aug25_19.50.kde)
#rm(HR95_Aug25_19)

study_HR_Aug25_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Aug25_19,HR90_Aug25_19,HR85_Aug25_19,HR50_Aug25_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

fish <- fish[!fish$TRANSMITTER == "566mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
study_HR_Aug25_19 <- merge(fish,study_HR_Aug25_19, by = "TRANSMITTER")
study_HR_Aug25_19 <- study_HR_Aug25_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Aug25_19,HR85_Aug25_19,HR90_Aug25_19,HR95_Aug25_19)

#write.csv(study_HR_Aug25_19,"BurbHR_Aug25_2019_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Sep <- as.Date("2019-09-25")

DATE_Sep - 12
DATE_Sep + 12

DATE_Sep1 <- as.Date("2019-09-13")
DATE_Sep2 <- as.Date("2019-10-07")

Sep25_19 <- myfunc(DATE_Sep1,DATE_Sep2)    

################################
# Sep 2019
################################
Sep2019_sum_tag =
  Sep25_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish that have less than 30 decections
Sep25_19 <- Sep25_19[!Sep25_19$TRANSMITTER == "367mm_U_a_t",]
Sep25_19 <- Sep25_19[!Sep25_19$TRANSMITTER == "566mm_U_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Sep25_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Sep25_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Sep25_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Sep25_19)
#calculate kernel density estimate for all fish
kud_Sep25_19 <- kernelUD(burb_Sep25_19[,1], h="href", grid = xy)
#rm(kud_Sep25_19)
#print all fish estimates
image(kud_Sep25_19)

#calculate kernel area estimates
N_Sep25_19 <- kernel.area(kud_Sep25_19, percent = seq(50,95, by = 5), unout = "km2")
N_Sep25_19 <- as.data.frame(t(N_Sep25_19))
colnames(N_Sep25_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")

#rm(N_Sep25_19)
#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

fish <- fish[!fish$TRANSMITTER == "367mm_U_a_t",]
fish <- fish[!fish$TRANSMITTER == "566mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Sep25_19 <- cbind(fish,N_Sep25_19)
HR_Burbs_Sep25_19 <- HR_Burbs_Sep25_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Sep25_19,"BurbHR_Sep25_19.csv")
################################################################################
burb.Sep25_19.95.kde <- getverticeshr(kud_Sep25_19, percent = 95, unin = "m", unout = "km2")
burb.Sep25_19.90.kde <- getverticeshr(kud_Sep25_19, percent = 90, unin = "m", unout = "km2")
burb.Sep25_19.85.kde <- getverticeshr(kud_Sep25_19, percent = 85, unin = "m", unout = "km2")
burb.Sep25_19.50.kde <- getverticeshr(kud_Sep25_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Sep25_19.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Sep25_19.95.kde)

HR95_Sep25_19 <- HR_CLIP(burb.Sep25_19.95.kde)
HR90_Sep25_19 <- HR_CLIP(burb.Sep25_19.90.kde)
HR85_Sep25_19 <- HR_CLIP(burb.Sep25_19.85.kde)
HR50_Sep25_19 <- HR_CLIP(burb.Sep25_19.50.kde)
#rm(HR95_Sep25_19)
study_HR_Sep25_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Sep25_19,HR90_Sep25_19,HR85_Sep25_19,HR50_Sep25_19))

#condense all data in to file with HR estimates and fish info
study_HR_Sep25_19 <- merge(fish,study_HR_Sep25_19, by = "TRANSMITTER")
study_HR_Sep25_19 <- study_HR_Sep25_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Sep25_19,HR85_Sep25_19,HR90_Sep25_19,HR95_Sep25_19)

#write.csv(study_HR_Sep25_19,"BurbHR_Sep25_19_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Oct <- as.Date("2019-10-25")

DATE_Oct - 12
DATE_Oct + 12

DATE_Oct1 <- as.Date("2019-10-13")
DATE_Oct2 <- as.Date("2019-11-06")

Oct25_19 <- myfunc(DATE_Oct1,DATE_Oct2)    

################################
# Oct 2019
################################
Oct2019_sum_tag =
  Oct25_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Oct25_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Oct25_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Oct25_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Oct25_19)
#calculate kernel density estimate for all fish
kud_Oct25_19 <- kernelUD(burb_Oct25_19[,1], h="href", grid = xy)
#rm(kud_Oct25_19)
#print all fish estimates
image(kud_Oct25_19)

#calculate kernel area estimates
N_Oct25_19 <- kernel.area(kud_Oct25_19, percent = seq(50,95, by = 5), unout = "km2")
N_Oct25_19 <- as.data.frame(t(N_Oct25_19))
colnames(N_Oct25_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Oct25_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Oct25_19 <- cbind(fish,N_Oct25_19)
HR_Burbs_Oct25_19 <- HR_Burbs_Oct25_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Oct25_19,"BurbHR_Oct25_19.csv")
################################################################################
burb.Oct25_19.95.kde <- getverticeshr(kud_Oct25_19, percent = 95, unin = "m", unout = "km2")
burb.Oct25_19.90.kde <- getverticeshr(kud_Oct25_19, percent = 90, unin = "m", unout = "km2")
burb.Oct25_19.85.kde <- getverticeshr(kud_Oct25_19, percent = 85, unin = "m", unout = "km2")
burb.Oct25_19.50.kde <- getverticeshr(kud_Oct25_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Oct25_19.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Oct25_19.95.kde)

HR95_Oct25_19 <- HR_CLIP(burb.Oct25_19.95.kde)
HR90_Oct25_19 <- HR_CLIP(burb.Oct25_19.90.kde)
HR85_Oct25_19 <- HR_CLIP(burb.Oct25_19.85.kde)
HR50_Oct25_19 <- HR_CLIP(burb.Oct25_19.50.kde)
#rm(HR95_Oct25_19)

study_HR_Oct25_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Oct25_19,HR90_Oct25_19,HR85_Oct25_19,HR50_Oct25_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Oct25_19 <- merge(fish,study_HR_Oct25_19, by = "TRANSMITTER")
study_HR_Oct25_19 <- study_HR_Oct25_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Oct25_19,HR85_Oct25_19,HR90_Oct25_19,HR95_Oct25_19)

#write.csv(study_HR_Oct25_19,"BurbHR_Oct25_19_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Nov <- as.Date("2019-11-25")

DATE_Nov - 12
DATE_Nov + 12

DATE_Nov1 <- as.Date("2019-11-13")
DATE_Nov2 <- as.Date("2019-12-07")

Nov25_19 <- myfunc(DATE_Nov1,DATE_Nov2)    

################################
# Nov 2019
################################
Nov2019_sum_tag =
  Nov25_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Nov25_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Nov25_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Nov25_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Nov25_19)
#calculate kernel density estimate for all fish
kud_Nov25_19 <- kernelUD(burb_Nov25_19[,1], h="href", grid = xy)
#rm(kud_Nov25_19)
#print all fish estimates
image(kud_Nov25_19)

#calculate kernel area estimates
N_Nov25_19 <- kernel.area(kud_Nov25_19, percent = seq(50,95, by = 5), unout = "km2")
N_Nov25_19 <- as.data.frame(t(N_Nov25_19))
colnames(N_Nov25_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Nov25_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Nov25_19 <- cbind(fish,N_Nov25_19)
HR_Burbs_Nov25_19 <- HR_Burbs_Nov25_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Nov25_19,"BurbHR_Nov25_19.csv")
################################################################################

burb.Nov25_19.95.kde <- getverticeshr(kud_Nov25_19, percent = 95, unin = "m", unout = "km2")
burb.Nov25_19.90.kde <- getverticeshr(kud_Nov25_19, percent = 90, unin = "m", unout = "km2")
burb.Nov25_19.85.kde <- getverticeshr(kud_Nov25_19, percent = 85, unin = "m", unout = "km2")
burb.Nov25_19.50.kde <- getverticeshr(kud_Nov25_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Nov25_19.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf
TRANSMITTER <- row.names(burb.Nov25_19.95.kde)

#The Kernel Density Estimates overlap with land!
################################################################################

HR95_Nov25_19 <- HR_CLIP(burb.Nov25_19.95.kde)
HR90_Nov25_19 <- HR_CLIP(burb.Nov25_19.90.kde)
HR85_Nov25_19 <- HR_CLIP(burb.Nov25_19.85.kde)
HR50_Nov25_19 <- HR_CLIP(burb.Nov25_19.50.kde)
#rm(HR95_Nov25_19)

study_HR_Nov25_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Nov25_19,HR90_Nov25_19,HR85_Nov25_19,HR50_Nov25_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Nov25_19 <- merge(fish,study_HR_Nov25_19, by = "TRANSMITTER")
study_HR_Nov25_19 <- study_HR_Nov25_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Nov25_19,HR85_Nov25_19,HR90_Nov25_19,HR95_Nov25_19)

#write.csv(study_HR_Nov25_19,"BurbHR_Nov25_19_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Dec <- as.Date("2019-12-25")

DATE_Dec - 12
DATE_Dec + 12

DATE_Dec1 <- as.Date("2019-12-13")
DATE_Dec2 <- as.Date("2020-01-06")

Dec25_19 <- myfunc(DATE_Dec1,DATE_Dec2)    

################################
# Dec 2019
################################
Dec2019_sum_tag =
  Dec25_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Dec25_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Dec25_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Dec25_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Dec25_19)
#calculate kernel density estimate for all fish
kud_Dec25_19 <- kernelUD(burb_Dec25_19[,1], h="href", grid = xy)
#rm(kud_Dec25_19)

#print all fish estimates
image(kud_Dec25_19)

#calculate kernel area estimates
N_Dec25_19 <- kernel.area(kud_Dec25_19, percent = seq(50,95, by = 5), unout = "km2")
N_Dec25_19 <- as.data.frame(t(N_Dec25_19))
colnames(N_Dec25_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Dec25_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Dec25_19 <- cbind(fish,N_Dec25_19)
HR_Burbs_Dec25_19 <- HR_Burbs_Dec25_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Dec25_19,"BurbHR_Dec25_19.csv")
################################################################################
burb.Dec25_19.95.kde <- getverticeshr(kud_Dec25_19, percent = 95, unin = "m", unout = "km2")
burb.Dec25_19.90.kde <- getverticeshr(kud_Dec25_19, percent = 90, unin = "m", unout = "km2")
burb.Dec25_19.85.kde <- getverticeshr(kud_Dec25_19, percent = 85, unin = "m", unout = "km2")
burb.Dec25_19.50.kde <- getverticeshr(kud_Dec25_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Dec25_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Dec25_19.95.kde)

HR95_Dec25_19 <- HR_CLIP(burb.Dec25_19.95.kde)
HR90_Dec25_19 <- HR_CLIP(burb.Dec25_19.90.kde)
HR85_Dec25_19 <- HR_CLIP(burb.Dec25_19.85.kde)
HR50_Dec25_19 <- HR_CLIP(burb.Dec25_19.50.kde)
#rm(HR95_Dec25_19)

study_HR_Dec25_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Dec25_19,HR90_Dec25_19,HR85_Dec25_19,HR50_Dec25_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Dec25_19 <- merge(fish,study_HR_Dec25_19, by = "TRANSMITTER")
study_HR_Dec25_19 <- study_HR_Dec25_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Dec25_19,HR85_Dec25_19,HR90_Dec25_19,HR95_Dec25_19)

#write.csv(study_HR_Dec25_19,"BurbHR_Dec25_19_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Jan20 <- as.Date("2020-01-25")

DATE_Jan20 - 12
DATE_Jan20 + 12

DATE_Jan201 <- as.Date("2020-01-13")
DATE_Jan202 <- as.Date("2020-02-06")

Jan25_20 <- myfunc(DATE_Jan201,DATE_Jan202)    

################################
# Jan 2020
################################
Jan2020_sum_tag =
  Jan25_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Jan25_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Jan25_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Jan25_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Jan25_20)
#calculate kernel density estimate for all fish
kud_Jan25_20 <- kernelUD(burb_Jan25_20[,1], h="href", grid = xy)
#rm(kud_Jan25_20)

#print all fish estimates
image(kud_Jan25_20)

#calculate kernel area estimates
N_Jan25_20 <- kernel.area(kud_Jan25_20, percent = seq(50,95, by = 5), unout = "km2")
N_Jan25_20 <- as.data.frame(t(N_Jan25_20))
colnames(N_Jan25_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Jan25_20)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "524mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Jan25_20 <- cbind(fish,N_Jan25_20)
HR_Burbs_Jan25_20 <- HR_Burbs_Jan25_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Jan25_20,"BurbHR_Jan25_20.csv")

################################################################################
burb.Jan25_20.95.kde <- getverticeshr(kud_Jan25_20, percent = 95, unin = "m", unout = "km2")
burb.Jan25_20.90.kde <- getverticeshr(kud_Jan25_20, percent = 90, unin = "m", unout = "km2")
burb.Jan25_20.85.kde <- getverticeshr(kud_Jan25_20, percent = 85, unin = "m", unout = "km2")
burb.Jan25_20.50.kde <- getverticeshr(kud_Jan25_20, percent = 50, unin = "m", unout = "km2")
#rm(burb.Jan25_20.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
#needed for step in HR_clip function
TRANSMITTER <- row.names(burb.Jan25_20.95.kde)

HR95_Jan25_20 <- HR_CLIP(burb.Jan25_20.95.kde)
HR90_Jan25_20 <- HR_CLIP(burb.Jan25_20.90.kde)
HR85_Jan25_20 <- HR_CLIP(burb.Jan25_20.85.kde)
HR50_Jan25_20 <- HR_CLIP(burb.Jan25_20.50.kde)
#rm(HR95_Jan25_20)

study_HR_Jan25_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Jan25_20,HR90_Jan25_20,HR85_Jan25_20,HR50_Jan25_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "524mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
study_HR_Jan25_20 <- merge(fish,study_HR_Jan25_20, by = "TRANSMITTER")
study_HR_Jan25_20 <- study_HR_Jan25_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Jan25_20,HR85_Jan25_20,HR90_Jan25_20,HR95_Jan25_20)

#write.csv(study_HR_Jan25_20,"BurbHR_Jan25_20_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Feb20 <- as.Date("2020-02-25")

DATE_Feb20 - 12
DATE_Feb20 + 12

DATE_Feb201 <- as.Date("2020-02-13")
DATE_Feb202 <- as.Date("2020-03-08")

Feb25_20 <- myfunc(DATE_Feb201,DATE_Feb202)    

################################
# Feb 2020
################################
Feb2020_sum_tag =
  Feb25_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

Feb25_20 <- Feb25_20[!Feb25_20$TRANSMITTER == "524mm_U_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Feb25_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Feb25_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Feb25_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Feb25_20)
#calculate kernel density estimate for all fish
kud_Feb25_20 <- kernelUD(burb_Feb25_20[,1], h="href", grid = xy)
#rm(kud_Feb25_20)
#print all fish estimates
image(kud_Feb25_20)

#calculate kernel area estimates
N_Feb25_20 <- kernel.area(kud_Feb25_20, percent = seq(50,95, by = 5), unout = "km2")
N_Feb25_20 <- as.data.frame(t(N_Feb25_20))
colnames(N_Feb25_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Feb25_20)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "524mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Feb25_20 <- cbind(fish,N_Feb25_20)
HR_Burbs_Feb25_20 <- HR_Burbs_Feb25_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Feb25_20,"BurbHR_Feb25_20.csv")
################################################################################
burb.Feb25_20.95.kde <- getverticeshr(kud_Feb25_20, percent = 95, unin = "m", unout = "km2")
burb.Feb25_20.90.kde <- getverticeshr(kud_Feb25_20, percent = 90, unin = "m", unout = "km2")
burb.Feb25_20.85.kde <- getverticeshr(kud_Feb25_20, percent = 85, unin = "m", unout = "km2")
burb.Feb25_20.50.kde <- getverticeshr(kud_Feb25_20, percent = 50, unin = "m", unout = "km2")
#rm(burb.Feb25_20.95.kde)
#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
#needed for step in HR_clip function
TRANSMITTER <- row.names(burb.Feb25_20.95.kde)

HR95_Feb25_20 <- HR_CLIP(burb.Feb25_20.95.kde)
HR90_Feb25_20 <- HR_CLIP(burb.Feb25_20.90.kde)
HR85_Feb25_20 <- HR_CLIP(burb.Feb25_20.85.kde)
HR50_Feb25_20 <- HR_CLIP(burb.Feb25_20.50.kde)
#rm(HR95_Feb25_20)

study_HR_Feb25_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Feb25_20,HR90_Feb25_20,HR85_Feb25_20,HR50_Feb25_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "524mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
study_HR_Feb25_20 <- merge(fish,study_HR_Feb25_20, by = "TRANSMITTER")
study_HR_Feb25_20 <- study_HR_Feb25_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Feb25_20,HR85_Feb25_20,HR90_Feb25_20,HR95_Feb25_20)

#write.csv(study_HR_Feb25_20,"BurbHR_Feb25_20_TRIMMED.csv")

################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Mar20 <- as.Date("2020-03-25")

DATE_Mar20 - 12
DATE_Mar20 + 12

DATE_Mar201 <- as.Date("2020-03-13")
DATE_Mar202 <- as.Date("2020-04-06")

Mar25_20 <- myfunc(DATE_Mar201,DATE_Mar202)    

################################
# Mar 2020
################################
Mar2020_sum_tag =
  Mar25_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Mar25_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Mar25_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Mar25_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Mar25_20)

#calculate kernel density estimate for all fish
kud_Mar25_20 <- kernelUD(burb_Mar25_20[,1], h="href", grid = xy)
#rm(kud_Mar25_20)
#print all fish estimates
image(kud_Mar25_20)

#rasterized home range estimate 95 kud
vud <- getvolumeUD(kud_Mar25_20)
vud
#rm(vud)
## Set up graphical parameters
par(mfrow=c(1,2))
par(mar=c(0,0,2,0))
## The output of kernelUD for the animal
image(kud_Mar25_20[[22]])
title("Output of kernelUD")
## Convert into a suitable data structure for the use of contour
xyz <- as.image.SpatialGridDataFrame(kud_Mar25_20[[22]])
contour(xyz, add=TRUE)
## and similarly for the output of getvolumeUD
par(mar=c(0,0,2,0))
image(vud[[22]])
title("Output of getvolumeUD")
xyzv <- as.image.SpatialGridDataFrame(vud[[22]])
contour(xyzv, add=TRUE)
## store the volume under the UD (as computed by getvolumeUD)
## of the animal in fud
fud <- vud[[22]]
## store the value of the volume under UD in a vector hr90
hr90 <- as.data.frame(fud)[,1]
## if hr90 is <= 90 then the pixel belongs to the home range
## (takes the value 1, 0 otherwise)
hr90 <- as.numeric(hr90 <= 90)
## Converts into a data frame
hr90 <- data.frame(hr90)
## Converts to a SpatialPixelsDataFrame
coordinates(hr90) <- coordinates(vud[[1]])
gridded(hr90) <- TRUE
## display the results
image(hr90)

#calculate kernel area estimates
N_Mar25_20 <- kernel.area(kud_Mar25_20, percent = seq(50,95, by = 5), unout = "km2")
N_Mar25_20 <- as.data.frame(t(N_Mar25_20))
colnames(N_Mar25_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Mar25_20)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Mar25_20 <- cbind(fish,N_Mar25_20)
HR_Burbs_Mar25_20 <- HR_Burbs_Mar25_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Mar25_20,"BurbHR_Mar25_20.csv")
################################################################################
burb.Mar25_20.95.kde <- getverticeshr(kud_Mar25_20, percent = 95, unin = "m", unout = "km2")
burb.Mar25_20.90.kde <- getverticeshr(kud_Mar25_20, percent = 90, unin = "m", unout = "km2")
burb.Mar25_20.85.kde <- getverticeshr(kud_Mar25_20, percent = 85, unin = "m", unout = "km2")
burb.Mar25_20.50.kde <- getverticeshr(kud_Mar25_20, percent = 50, unin = "m", unout = "km2")
#rm(burb.Mar25_20.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
#needed for step in HR_clip function
TRANSMITTER <- row.names(burb.Mar25_20.95.kde)

HR95_Mar25_20 <- HR_CLIP(burb.Mar25_20.95.kde)
HR90_Mar25_20 <- HR_CLIP(burb.Mar25_20.90.kde)
HR85_Mar25_20 <- HR_CLIP(burb.Mar25_20.85.kde)
HR50_Mar25_20 <- HR_CLIP(burb.Mar25_20.50.kde)
#rm(HR95_Mar25_20)

study_HR_Mar25_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Mar25_20,HR90_Mar25_20,HR85_Mar25_20,HR50_Mar25_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Mar25_20 <- merge(fish,study_HR_Mar25_20, by = "TRANSMITTER")
study_HR_Mar25_20 <- study_HR_Mar25_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Mar25_20,HR85_Mar25_20,HR90_Mar25_20,HR95_Mar25_20)

setwd("E:/BSUGradBurbot/DATA/FishData/FinalFilteredData/HomeRange_AllSeasons/LakeTrim")

#write.csv(study_HR_Mar25_20,"BurbHR_Mar25_20_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Apr20 <- as.Date("2020-04-25")

DATE_Apr20 - 12
DATE_Apr20 + 12

DATE_Apr201 <- as.Date("2020-04-13")
DATE_Apr202 <- as.Date("2020-05-07")

Apr25_20 <- myfunc(DATE_Apr201,DATE_Apr202)    

################################
# Apr 2020
################################
Apr2020_sum_tag =
  Apr25_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Apr25_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Apr25_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Apr25_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
rm(burb_Apr25_20)

#calculate kernel density estimate for all fish
kud_Apr25_20 <- kernelUD(burb_Apr25_20[,1], h="href", grid = xy)
#rm(kud_Apr25_20)

#print all fish estimates
image(kud_Apr25_20)

#calculate kernel area estimates
N_Apr25_20 <- kernel.area(kud_Apr25_20, percent = seq(50,95, by = 5), unout = "km2")
N_Apr25_20 <- as.data.frame(t(N_Apr25_20))
colnames(N_Apr25_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Apr25_20)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Apr25_20 <- cbind(fish,N_Apr25_20)
HR_Burbs_Apr25_20 <- HR_Burbs_Apr25_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Apr25_20,"BurbHR_Apr25_20.csv")
################################################################################
burb.Apr25_20.95.kde <- getverticeshr(kud_Apr25_20, percent = 95, unin = "m", unout = "km2")
burb.Apr25_20.90.kde <- getverticeshr(kud_Apr25_20, percent = 90, unin = "m", unout = "km2")
burb.Apr25_20.85.kde <- getverticeshr(kud_Apr25_20, percent = 85, unin = "m", unout = "km2")
burb.Apr25_20.50.kde <- getverticeshr(kud_Apr25_20, percent = 50, unin = "m", unout = "km2")
#rm(burb.Apr25_20.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
#needed for step in HR_clip function
TRANSMITTER <- row.names(burb.Apr25_20.95.kde)

HR95_Apr25_20 <- HR_CLIP(burb.Apr25_20.95.kde)
HR90_Apr25_20 <- HR_CLIP(burb.Apr25_20.90.kde)
HR85_Apr25_20 <- HR_CLIP(burb.Apr25_20.85.kde)
HR50_Apr25_20 <- HR_CLIP(burb.Apr25_20.50.kde)
#rm(HR95_Apr25_20)

study_HR_Apr25_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Apr25_20,HR90_Apr25_20,HR85_Apr25_20,HR50_Apr25_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Apr25_20 <- merge(fish,study_HR_Apr25_20, by = "TRANSMITTER")
study_HR_Apr25_20 <- study_HR_Apr25_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Apr25_20,HR85_Apr25_20,HR90_Apr25_20,HR95_Apr25_20)

#write.csv(study_HR_Apr25_20,"BurbHR_Apr25_20_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_May20 <- as.Date("2020-05-25")

DATE_May20 - 12
DATE_May20 + 12

DATE_May201 <- as.Date("2020-05-13")
DATE_May202 <- as.Date("2020-06-06")

May25_20 <- myfunc(DATE_May201,DATE_May202)    

################################
# May 2020
################################
May2020_sum_tag =
  May25_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish with less than 30 detections
May25_20 <- May25_20[!May25_20$TRANSMITTER == "492mm_F_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- May25_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(May25_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_May25_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_May25_20)

#calculate kernel density estimate for all fish
kud_May25_20 <- kernelUD(burb_May25_20[,1], h="href", grid = xy)
#rm(kud_May25_20)

#print all fish estimates
image(kud_May25_20)

#calculate kernel area estimates
N_May25_20 <- kernel.area(kud_May25_20, percent = seq(50,95, by = 5), unout = "km2")
N_May25_20 <- as.data.frame(t(N_May25_20))
colnames(N_May25_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_May25_20)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

fish <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_May25_20 <- cbind(fish,N_May25_20)
HR_Burbs_May25_20 <- HR_Burbs_May25_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_May25_20,"BurbHR_May25_20.csv")
################################################################################
burb.May25_20.95.kde <- getverticeshr(kud_May25_20, percent = 95, unin = "m", unout = "km2")
burb.May25_20.90.kde <- getverticeshr(kud_May25_20, percent = 90, unin = "m", unout = "km2")
burb.May25_20.85.kde <- getverticeshr(kud_May25_20, percent = 85, unin = "m", unout = "km2")
burb.May25_20.50.kde <- getverticeshr(kud_May25_20, percent = 50, unin = "m", unout = "km2")
#rm(burb.May25_20.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
#needed for step in HR_clip function
TRANSMITTER <- row.names(burb.May25_20.95.kde)

HR95_May25_20 <- HR_CLIP(burb.May25_20.95.kde)
HR90_May25_20 <- HR_CLIP(burb.May25_20.90.kde)
HR85_May25_20 <- HR_CLIP(burb.May25_20.85.kde)
HR50_May25_20 <- HR_CLIP(burb.May25_20.50.kde)
#rm(HR95_May25_20)

study_HR_May25_20 <- as.data.frame(cbind(TRANSMITTER,HR95_May25_20,HR90_May25_20,HR85_May25_20,HR50_May25_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]

#condense all data in to file with HR estimates and fish info
study_HR_May25_20 <- merge(fish,study_HR_May25_20, by = "TRANSMITTER")
study_HR_May25_20 <- study_HR_May25_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_May25_20,HR85_May25_20,HR90_May25_20,HR95_May25_20)

#write.csv(study_HR_May25_20,"BurbHR_May25_20_TRIMMED.csv")


################################################################################
################################################################################
############### 25 day span home range across study period #####################
################################################################################
################################################################################

HR_25day <- BurbLocsHR

# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Apr19_1 <- as.Date("2019-04-13")

DATE_Apr19_1 - 12
DATE_Apr19_1 + 12

DATE_Apr191 <- as.Date("2019-04-01")
DATE_Apr192 <- as.Date("2019-04-25")

Apr13_19 <- myfunc(DATE_Apr191,DATE_Apr192)    

################################
# April 2019
################################
Apr2019_sum_tag_1 =
  Apr13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Apr13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Apr13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Apr13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Apr13_19)
#calculate kernel density estimate for all fish
kud_Apr13_19 <- kernelUD(burb_Apr13_19[,1], h="href", grid = xy)
#rm(kud_Apr13_19)

#print all fish estimates
image(kud_Apr13_19)

#calculate kernel area estimates
N_Apr13_19 <- kernel.area(kud_Apr13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Apr13_19 <- as.data.frame(t(N_Apr13_19))
colnames(N_Apr13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Apr13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "504mm_F_a_t",]
fish <- fish[!fish$TRANSMITTER == "625mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Apr13_19 <- cbind(fish,N_Apr13_19)
HR_Burbs_Apr13_19 <- HR_Burbs_Apr13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Apr13_19,"BurbHR_Apr13_19.csv")

################################################################################
burb.Apr13_19.95.kde <- getverticeshr(kud_Apr13_19, percent = 95, unin = "m", unout = "km2")
burb.Apr13_19.90.kde <- getverticeshr(kud_Apr13_19, percent = 90, unin = "m", unout = "km2")
burb.Apr13_19.85.kde <- getverticeshr(kud_Apr13_19, percent = 85, unin = "m", unout = "km2")
burb.Apr13_19.50.kde <- getverticeshr(kud_Apr13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Apr13_19.50.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Apr13_19.95.kde)

HR95_Apr13_19 <- HR_CLIP(burb.Apr13_19.95.kde)
HR90_Apr13_19 <- HR_CLIP(burb.Apr13_19.90.kde)
HR85_Apr13_19 <- HR_CLIP(burb.Apr13_19.85.kde)
HR50_Apr13_19 <- HR_CLIP(burb.Apr13_19.50.kde)
#rm(HR95_Apr13_19)

study_HR_Apr13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Apr13_19,HR90_Apr13_19,HR85_Apr13_19,HR50_Apr13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Apr13_19 <- merge(fish,study_HR_Apr13_19, by = "TRANSMITTER")
study_HR_Apr13_19 <- study_HR_Apr13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Apr13_19,HR85_Apr13_19,HR90_Apr13_19,HR95_Apr13_19)

#write.csv(study_HR_Apr13_19,"BurbHR_Apr13_19_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_May19_1 <- as.Date("2019-05-13")

DATE_May19_1 - 12
DATE_May19_1 + 12

DATE_May191 <- as.Date("2019-05-01")
DATE_May192 <- as.Date("2019-05-25")

May13_19 <- myfunc(DATE_May191,DATE_May192)

################################
# May 2019
################################
May2019_sum_tag_1 =
  May13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- May13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(May13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_May13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_May13_19)
#calculate kernel density estimate for all fish
kud_May13_19 <- kernelUD(burb_May13_19[,1], h="href", grid = xy)
#rm(kud_May13_19)
#print all fish estimates
image(kud_May13_19)

#calculate kernel area estimates
N_May13_19 <- kernel.area(kud_May13_19, percent = seq(50,95, by = 5), unout = "km2")
N_May13_19 <- as.data.frame(t(N_May13_19))
colnames(N_May13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_May13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_May13_19 <- cbind(fish,N_May13_19)
HR_Burbs_May13_19 <- HR_Burbs_May13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_May13_19,"BurbHR_May13_19.csv")
################################################################################
#get area HR for each percentage
burb.May13_19.95.kde <- getverticeshr(kud_May13_19, percent = 95, unin = "m", unout = "km2")
burb.May13_19.90.kde <- getverticeshr(kud_May13_19, percent = 90, unin = "m", unout = "km2")
burb.May13_19.85.kde <- getverticeshr(kud_May13_19, percent = 85, unin = "m", unout = "km2")
burb.May13_19.50.kde <- getverticeshr(kud_May13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.May13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.May13_19.95.kde)

HR95_May13_19 <- HR_CLIP(burb.May13_19.95.kde)
HR90_May13_19 <- HR_CLIP(burb.May13_19.90.kde)
HR85_May13_19 <- HR_CLIP(burb.May13_19.85.kde)
HR50_May13_19 <- HR_CLIP(burb.May13_19.50.kde)
#rm(HR95_May13_19)

study_HR_May13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_May13_19,HR90_May13_19,HR85_May13_19,HR50_May13_19))

#condense all data in to file with HR estimates and fish info
study_HR_May13_19 <- merge(fish,study_HR_May13_19, by = "TRANSMITTER")
study_HR_May13_19 <- study_HR_May13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_May13_19,HR85_May13_19,HR90_May13_19,HR95_May13_19)

#write.csv(study_HR_May13_19,"BurbHR_May13_19_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Jun19_1 <- as.Date("2019-06-13")

DATE_Jun19_1 - 12
DATE_Jun19_1 + 12

DATE_Jun191 <- as.Date("2019-06-01")
DATE_Jun192 <- as.Date("2019-06-25")

Jun13_19 <- myfunc(DATE_Jun191,DATE_Jun192)

################################
# June 2019
################################
Jun2019_sum_tag_1 =
  Jun13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

Jun13_19 <- Jun13_19[!Jun13_19$TRANSMITTER == "492mm_F_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Jun13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Jun13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Jun13_19 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)
#rm(burb_Jun13_19)

#calculate kernel density estimate for all fish
kud_Jun13_19 <- kernelUD(burb_Jun13_19[,1], h="href", grid = xy)
#rm(kud_Jun13_19)

#print all fish estimates
image(kud_Jun13_19)

#calculate kernel area estimates
N_Jun13_19 <- kernel.area(kud_Jun13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Jun13_19 <- as.data.frame(t(N_Jun13_19))
colnames(N_Jun13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Jun13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Jun13_19 <- cbind(fish,N_Jun13_19)
HR_Burbs_Jun13_19 <- HR_Burbs_Jun13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Jun13_19,"BurbHR_Jun13_19.csv")
################################################################################
#get area HR for each percentage
burb.Jun13_19.95.kde <- getverticeshr(kud_Jun13_19, percent = 95, unin = "m", unout = "km2")
burb.Jun13_19.90.kde <- getverticeshr(kud_Jun13_19, percent = 90, unin = "m", unout = "km2")
burb.Jun13_19.85.kde <- getverticeshr(kud_Jun13_19, percent = 85, unin = "m", unout = "km2")
burb.Jun13_19.50.kde <- getverticeshr(kud_Jun13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Jun13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Jun13_19.95.kde)

HR95_Jun13_19 <- HR_CLIP(burb.Jun13_19.95.kde)
HR90_Jun13_19 <- HR_CLIP(burb.Jun13_19.90.kde)
HR85_Jun13_19 <- HR_CLIP(burb.Jun13_19.85.kde)
HR50_Jun13_19 <- HR_CLIP(burb.Jun13_19.50.kde)
#rm(HR95_Jun13_19)

study_HR_Jun13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Jun13_19,HR90_Jun13_19,HR85_Jun13_19,HR50_Jun13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Jun13_19 <- merge(fish,study_HR_Jun13_19, by = "TRANSMITTER")
study_HR_Jun13_19 <- study_HR_Jun13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Jun13_19,HR85_Jun13_19,HR90_Jun13_19,HR95_Jun13_19)

#write.csv(study_HR_Jun13_19,"BurbHR_Jun13_19_TRIMMED.csv")

################################################################################

# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Jul19_1 <- as.Date("2019-07-13")

DATE_Jul19_1 - 12
DATE_Jul19_1 + 12

DATE_Jul191 <- as.Date("2019-07-01")
DATE_Jul192 <- as.Date("2019-07-25")

Jul13_19 <- myfunc(DATE_Jul191,DATE_Jul192)

################################
# July 2019
################################
Jul2019_sum_tag_1 =
  Jul13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish with less than 30 detections
Jul13_19 <- Jul13_19[!Jul13_19$TRANSMITTER == "492mm_F_a_t",]
Jul13_19 <- Jul13_19[!Jul13_19$TRANSMITTER == "527mm_M_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Jul13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Jul13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Jul13_19 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)
#rm(burb_Jul13_19)

#calculate kernel density estimate for all fish
kud_Jul13_19 <- kernelUD(burb_Jul13_19[,1], h="href", grid = xy)
#rm(kud_Jul13_19)

#print all fish estimates
image(kud_Jul13_19)

#calculate kernel area estimates
N_Jul13_19 <- kernel.area(kud_Jul13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Jul13_19 <- as.data.frame(t(N_Jul13_19))
colnames(N_Jul13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Jul13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]
fish <- fish[!fish$TRANSMITTER == "527mm_M_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Jul13_19 <- cbind(fish,N_Jul13_19)
HR_Burbs_Jul13_19 <- HR_Burbs_Jul13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Jul13_19,"BurbHR_Jul13_19.csv")
################################################################################
#get area HR for each percentage
burb.Jul13_19.95.kde <- getverticeshr(kud_Jul13_19, percent = 95, unin = "m", unout = "km2")
burb.Jul13_19.90.kde <- getverticeshr(kud_Jul13_19, percent = 90, unin = "m", unout = "km2")
burb.Jul13_19.85.kde <- getverticeshr(kud_Jul13_19, percent = 85, unin = "m", unout = "km2")
burb.Jul13_19.50.kde <- getverticeshr(kud_Jul13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Jul13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Jul13_19.95.kde)

HR95_Jul13_19 <- HR_CLIP(burb.Jul13_19.95.kde)
HR90_Jul13_19 <- HR_CLIP(burb.Jul13_19.90.kde)
HR85_Jul13_19 <- HR_CLIP(burb.Jul13_19.85.kde)
HR50_Jul13_19 <- HR_CLIP(burb.Jul13_19.50.kde)
#rm(HR95_Jul13_19)

study_HR_Jul13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Jul13_19,HR90_Jul13_19,HR85_Jul13_19,HR50_Jul13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Jul13_19 <- merge(fish,study_HR_Jul13_19, by = "TRANSMITTER")
study_HR_Jul13_19 <- study_HR_Jul13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Jul13_19,HR85_Jul13_19,HR90_Jul13_19,HR95_Jul13_19)

#write.csv(study_HR_Jul13_19,"BurbHR_Jul13_19_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Aug19_1 <- as.Date("2019-08-13")

DATE_Aug19_1 - 12
DATE_Aug19_1 + 12

DATE_Aug191 <- as.Date("2019-08-01")
DATE_Aug192 <- as.Date("2019-08-25")

Aug13_19 <- myfunc(DATE_Aug191,DATE_Aug192)

################################
# August 2019
################################
Aug2019_sum_tag_1 =
  Aug13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish with less than 30 detections
Aug13_19 <- Aug13_19[!Aug13_19$TRANSMITTER == "492mm_F_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Aug13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Aug13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Aug13_19 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)
#rm(burb_Aug13_19)

#calculate kernel density estimate for all fish
kud_Aug13_19 <- kernelUD(burb_Aug13_19[,1], h="href", grid = xy)
#rm(kud_Aug13_19)

#print all fish estimates
image(kud_Aug13_19)

#calculate kernel area estimates
N_Aug13_19 <- kernel.area(kud_Aug13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Aug13_19 <- as.data.frame(t(N_Aug13_19))
colnames(N_Aug13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Aug13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Aug13_19 <- cbind(fish,N_Aug13_19)
HR_Burbs_Aug13_19 <- HR_Burbs_Aug13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Aug13_19,"BurbHR_Aug13_19.csv")
################################################################################
#get area HR for each percentage
burb.Aug13_19.95.kde <- getverticeshr(kud_Aug13_19, percent = 93, unin = "m", unout = "km2")
#lack of new locations resulting in failure of 95 kud
burb.Aug13_19.90.kde <- getverticeshr(kud_Aug13_19, percent = 90, unin = "m", unout = "km2")
burb.Aug13_19.85.kde <- getverticeshr(kud_Aug13_19, percent = 85, unin = "m", unout = "km2")
burb.Aug13_19.50.kde <- getverticeshr(kud_Aug13_19, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Aug13_19.90.kde)

HR95_Aug13_19 <- HR_CLIP(burb.Aug13_19.95.kde)
HR90_Aug13_19 <- HR_CLIP(burb.Aug13_19.90.kde)
HR85_Aug13_19 <- HR_CLIP(burb.Aug13_19.85.kde)
HR50_Aug13_19 <- HR_CLIP(burb.Aug13_19.50.kde)

study_HR_Aug13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Aug13_19,HR90_Aug13_19,HR85_Aug13_19,HR50_Aug13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Aug13_19 <- merge(fish,study_HR_Aug13_19, by = "TRANSMITTER")
study_HR_Aug13_19 <- study_HR_Aug13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Aug13_19,HR85_Aug13_19,HR90_Aug13_19,HR95_Aug13_19)

#write.csv(study_HR_Aug13_19,"BurbHR_Aug13_19_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Sep19_1 <- as.Date("2019-09-13")

DATE_Sep19_1 - 12
DATE_Sep19_1 + 12

DATE_Sep191 <- as.Date("2019-09-01")
DATE_Sep192 <- as.Date("2019-09-25")

Sep13_19 <- myfunc(DATE_Sep191,DATE_Sep192)

################################
# September 2019
################################
Sep2019_sum_tag_1 =
  Sep13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish that results in KUD error
Sep13_19 <- Sep13_19[!Sep13_19$TRANSMITTER == "566mm_U_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Sep13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Sep13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Sep13_19 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)
#rm(burb_Sep13_19)

#calculate kernel density estimate for all fish
kud_Sep13_19 <- kernelUD(burb_Sep13_19[,1], h="href", grid = xy)
#rm(kud_Sep13_19)

#print all fish estimates
image(kud_Sep13_19)

#calculate kernel area estimates
N_Sep13_19 <- kernel.area(kud_Sep13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Sep13_19 <- as.data.frame(t(N_Sep13_19))
colnames(N_Sep13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Sep13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "566mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Sep13_19 <- cbind(fish,N_Sep13_19)
HR_Burbs_Sep13_19 <- HR_Burbs_Sep13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Sep13_19,"BurbHR_Sep13_19.csv")

################################################################################
#get area HR for each percentage
burb.Sep13_19.95.kde <- getverticeshr(kud_Sep13_19, percent = 95, unin = "m", unout = "km2")
burb.Sep13_19.90.kde <- getverticeshr(kud_Sep13_19, percent = 90, unin = "m", unout = "km2")
burb.Sep13_19.85.kde <- getverticeshr(kud_Sep13_19, percent = 85, unin = "m", unout = "km2")
burb.Sep13_19.50.kde <- getverticeshr(kud_Sep13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Sep13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Sep13_19.95.kde)

HR95_Sep13_19 <- HR_CLIP(burb.Sep13_19.95.kde)
HR90_Sep13_19 <- HR_CLIP(burb.Sep13_19.90.kde)
HR85_Sep13_19 <- HR_CLIP(burb.Sep13_19.85.kde)
HR50_Sep13_19 <- HR_CLIP(burb.Sep13_19.50.kde)
#rm(HR95_Sep13_19)

study_HR_Sep13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Sep13_19,HR90_Sep13_19,HR85_Sep13_19,HR50_Sep13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Sep13_19 <- merge(fish,study_HR_Sep13_19, by = "TRANSMITTER")
study_HR_Sep13_19 <- study_HR_Sep13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Sep13_19,HR85_Sep13_19,HR90_Sep13_19,HR95_Sep13_19)

#write.csv(study_HR_Sep13_19,"BurbHR_Sep13_19_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Oct19_1 <- as.Date("2019-10-13")

DATE_Oct19_1 - 12
DATE_Oct19_1 + 12

DATE_Oct191 <- as.Date("2019-10-01")
DATE_Oct192 <- as.Date("2019-10-25")

Oct13_19 <- myfunc(DATE_Oct191,DATE_Oct192)

################################
# October 2019
################################
Oct2019_sum_tag_1 =
  Oct13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Oct13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Oct13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Oct13_19 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)
#rm(burb_Oct13_19)

#calculate kernel density estimate for all fish
kud_Oct13_19 <- kernelUD(burb_Oct13_19[,1], h="href", grid = xy)
#rm(kud_Oct13_19)

#print all fish estimates
image(kud_Oct13_19)

#calculate kernel area estimates
N_Oct13_19 <- kernel.area(kud_Oct13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Oct13_19 <- as.data.frame(t(N_Oct13_19))
colnames(N_Oct13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Oct13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Oct13_19 <- cbind(fish,N_Oct13_19)
HR_Burbs_Oct13_19 <- HR_Burbs_Oct13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Oct13_19,"BurbHR_Oct13_19.csv")
################################################################################
#get area HR for each percentage
burb.Oct13_19.95.kde <- getverticeshr(kud_Oct13_19, percent = 95, unin = "m", unout = "km2")
burb.Oct13_19.90.kde <- getverticeshr(kud_Oct13_19, percent = 90, unin = "m", unout = "km2")
burb.Oct13_19.85.kde <- getverticeshr(kud_Oct13_19, percent = 85, unin = "m", unout = "km2")
burb.Oct13_19.50.kde <- getverticeshr(kud_Oct13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Oct13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Oct13_19.95.kde)

HR95_Oct13_19 <- HR_CLIP(burb.Oct13_19.95.kde)
HR90_Oct13_19 <- HR_CLIP(burb.Oct13_19.90.kde)
HR85_Oct13_19 <- HR_CLIP(burb.Oct13_19.85.kde)
HR50_Oct13_19 <- HR_CLIP(burb.Oct13_19.50.kde)
#rm(HR95_Oct13_19)

study_HR_Oct13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Oct13_19,HR90_Oct13_19,HR85_Oct13_19,HR50_Oct13_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Oct13_19 <- merge(fish,study_HR_Oct13_19, by = "TRANSMITTER")
study_HR_Oct13_19 <- study_HR_Oct13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Oct13_19,HR85_Oct13_19,HR90_Oct13_19,HR95_Oct13_19)

#write.csv(study_HR_Oct13_19,"BurbHR_Oct13_19_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Nov19_1 <- as.Date("2019-11-13")

DATE_Nov19_1 - 12
DATE_Nov19_1 + 12

DATE_Nov191 <- as.Date("2019-11-01")
DATE_Nov192 <- as.Date("2019-11-25")

Nov13_19 <- myfunc(DATE_Nov191,DATE_Nov192)

################################
# November 2019
################################
Nov2019_sum_tag_1 =
  Nov13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Nov13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Nov13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Nov13_19 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)
#rm(burb_Nov13_19)

#calculate kernel density estimate for all fish
kud_Nov13_19 <- kernelUD(burb_Nov13_19[,1], h="href", grid = xy)
#rm(kud_Nov13_19)

#print all fish estimates
image(kud_Nov13_19)

#calculate kernel area estimates
N_Nov13_19 <- kernel.area(kud_Nov13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Nov13_19 <- as.data.frame(t(N_Nov13_19))
colnames(N_Nov13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Nov13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Nov13_19 <- cbind(fish,N_Nov13_19)
HR_Burbs_Nov13_19 <- HR_Burbs_Nov13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Nov13_19,"BurbHR_Nov13_19.csv")
################################################################################
#get area HR for each percentage
burb.Nov13_19.95.kde <- getverticeshr(kud_Nov13_19, percent = 95, unin = "m", unout = "km2")
burb.Nov13_19.90.kde <- getverticeshr(kud_Nov13_19, percent = 90, unin = "m", unout = "km2")
burb.Nov13_19.85.kde <- getverticeshr(kud_Nov13_19, percent = 85, unin = "m", unout = "km2")
burb.Nov13_19.50.kde <- getverticeshr(kud_Nov13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Nov13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Nov13_19.95.kde)

HR95_Nov13_19 <- HR_CLIP(burb.Nov13_19.95.kde)
HR90_Nov13_19 <- HR_CLIP(burb.Nov13_19.90.kde)
HR85_Nov13_19 <- HR_CLIP(burb.Nov13_19.85.kde)
HR50_Nov13_19 <- HR_CLIP(burb.Nov13_19.50.kde)
#rm(HR95_Nov13_19)

study_HR_Nov13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Nov13_19,HR90_Nov13_19,HR85_Nov13_19,HR50_Nov13_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Nov13_19 <- merge(fish,study_HR_Nov13_19, by = "TRANSMITTER")
study_HR_Nov13_19 <- study_HR_Nov13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Nov13_19,HR85_Nov13_19,HR90_Nov13_19,HR95_Nov13_19)

#write.csv(study_HR_Nov13_19,"BurbHR_Nov13_19_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Dec19_1 <- as.Date("2019-12-13")

DATE_Dec19_1 - 12
DATE_Dec19_1 + 12

DATE_Dec191 <- as.Date("2019-12-01")
DATE_Dec192 <- as.Date("2019-12-25")

Dec13_19 <- myfunc(DATE_Dec191,DATE_Dec192)

################################
# December 2019
################################
Dec2019_sum_tag_1 =
  Dec13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Dec13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Dec13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Dec13_19 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)
#rm(burb_Dec13_19)

#calculate kernel density estimate for all fish
kud_Dec13_19 <- kernelUD(burb_Dec13_19[,1], h="href", grid = xy)
#rm(kud_Dec13_19)

#print all fish estimates
image(kud_Dec13_19)

#calculate kernel area estimates
N_Dec13_19 <- kernel.area(kud_Dec13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Dec13_19 <- as.data.frame(t(N_Dec13_19))
colnames(N_Dec13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Dec13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Dec13_19 <- cbind(fish,N_Dec13_19)
HR_Burbs_Dec13_19 <- HR_Burbs_Dec13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Dec13_19,"BurbHR_Dec13_19.csv")
################################################################################
#get area HR for each percentage
burb.Dec13_19.95.kde <- getverticeshr(kud_Dec13_19, percent = 95, unin = "m", unout = "km2")
burb.Dec13_19.90.kde <- getverticeshr(kud_Dec13_19, percent = 90, unin = "m", unout = "km2")
burb.Dec13_19.85.kde <- getverticeshr(kud_Dec13_19, percent = 85, unin = "m", unout = "km2")
burb.Dec13_19.50.kde <- getverticeshr(kud_Dec13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Dec13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf
#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Dec13_19.95.kde)

HR95_Dec13_19 <- HR_CLIP(burb.Dec13_19.95.kde)
HR90_Dec13_19 <- HR_CLIP(burb.Dec13_19.90.kde)
HR85_Dec13_19 <- HR_CLIP(burb.Dec13_19.85.kde)
HR50_Dec13_19 <- HR_CLIP(burb.Dec13_19.50.kde)
#rm(HR95_Dec13_19)

study_HR_Dec13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Dec13_19,HR90_Dec13_19,HR85_Dec13_19,HR50_Dec13_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Dec13_19 <- merge(fish,study_HR_Dec13_19, by = "TRANSMITTER")
study_HR_Dec13_19 <- study_HR_Dec13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Dec13_19,HR85_Dec13_19,HR90_Dec13_19,HR95_Dec13_19)

#write.csv(study_HR_Dec13_19,"BurbHR_Dec13_19_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Jan20_1 <- as.Date("2020-01-13")

DATE_Jan20_1 - 12
DATE_Jan20_1 + 12

DATE_Jan201 <- as.Date("2020-01-01")
DATE_Jan202 <- as.Date("2020-01-25")

Jan13_20 <- myfunc(DATE_Jan201,DATE_Jan202)

################################
# January 2020
################################
Jan2020_sum_tag_1 =
  Jan13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Jan13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Jan13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Jan13_20 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)

#calculate kernel density estimate for all fish
kud_Jan13_20 <- kernelUD(burb_Jan13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_Jan13_20)

#calculate kernel area estimates
N_Jan13_20 <- kernel.area(kud_Jan13_20, percent = seq(50,95, by = 5), unout = "km2")
N_Jan13_20 <- as.data.frame(t(N_Jan13_20))
colnames(N_Jan13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Jan13_20 <- cbind(fish,N_Jan13_20)
HR_Burbs_Jan13_20 <- HR_Burbs_Jan13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Jan13_20,"BurbHR_Jan13_20.csv")
################################################################################
#get area HR for each percentage
burb.Jan13_20.95.kde <- getverticeshr(kud_Jan13_20, percent = 95, unin = "m", unout = "km2")
burb.Jan13_20.90.kde <- getverticeshr(kud_Jan13_20, percent = 90, unin = "m", unout = "km2")
burb.Jan13_20.85.kde <- getverticeshr(kud_Jan13_20, percent = 85, unin = "m", unout = "km2")
burb.Jan13_20.50.kde <- getverticeshr(kud_Jan13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Jan13_20.95.kde)

HR95_Jan13_20 <- HR_CLIP(burb.Jan13_20.95.kde)
HR90_Jan13_20 <- HR_CLIP(burb.Jan13_20.90.kde)
HR85_Jan13_20 <- HR_CLIP(burb.Jan13_20.85.kde)
HR50_Jan13_20 <- HR_CLIP(burb.Jan13_20.50.kde)

study_HR_Jan13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Jan13_20,HR90_Jan13_20,HR85_Jan13_20,HR50_Jan13_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Jan13_20 <- merge(fish,study_HR_Jan13_20, by = "TRANSMITTER")
study_HR_Jan13_20 <- study_HR_Jan13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Jan13_20,HR85_Jan13_20,HR90_Jan13_20,HR95_Jan13_20)

#write.csv(study_HR_Jan13_20,"BurbHR_Jan13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Feb20_1 <- as.Date("2020-02-13")

DATE_Feb20_1 - 12
DATE_Feb20_1 + 12

DATE_Feb201 <- as.Date("2020-02-01")
DATE_Feb202 <- as.Date("2020-02-25")

Feb13_20 <- myfunc(DATE_Feb201,DATE_Feb202)

################################
# February 2020
################################
Feb2019_sum_tag_1 =
  Feb13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# remove fish with less than 30 locations
Feb13_20 <- Feb13_20[!Feb13_20$TRANSMITTER == "524mm_U_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Feb13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Feb13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Feb13_20 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)

#calculate kernel density estimate for all fish
kud_Feb13_20 <- kernelUD(burb_Feb13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_Feb13_20)

#calculate kernel area estimates
N_Feb13_20 <- kernel.area(kud_Feb13_20, percent = seq(50,95, by = 5), unout = "km2")
N_Feb13_20 <- as.data.frame(t(N_Feb13_20))
colnames(N_Feb13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "524mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Feb13_20 <- cbind(fish,N_Feb13_20)
HR_Burbs_Feb13_20 <- HR_Burbs_Feb13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Feb13_20,"BurbHR_Feb13_20.csv")
################################################################################
#get area HR for each percentage
burb.Feb13_20.95.kde <- getverticeshr(kud_Feb13_20, percent = 95, unin = "m", unout = "km2")
burb.Feb13_20.90.kde <- getverticeshr(kud_Feb13_20, percent = 90, unin = "m", unout = "km2")
burb.Feb13_20.85.kde <- getverticeshr(kud_Feb13_20, percent = 85, unin = "m", unout = "km2")
burb.Feb13_20.50.kde <- getverticeshr(kud_Feb13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Feb13_20.95.kde)

HR95_Feb13_20 <- HR_CLIP(burb.Feb13_20.95.kde)
HR90_Feb13_20 <- HR_CLIP(burb.Feb13_20.90.kde)
HR85_Feb13_20 <- HR_CLIP(burb.Feb13_20.85.kde)
HR50_Feb13_20 <- HR_CLIP(burb.Feb13_20.50.kde)

study_HR_Feb13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Feb13_20,HR90_Feb13_20,HR85_Feb13_20,HR50_Feb13_20))


#condense all data in to file with HR estimates and fish info
study_HR_Feb13_20 <- merge(fish,study_HR_Feb13_20, by = "TRANSMITTER")
study_HR_Feb13_20 <- study_HR_Feb13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Feb13_20,HR85_Feb13_20,HR90_Feb13_20,HR95_Feb13_20)

#write.csv(study_HR_Feb13_20,"BurbHR_Feb13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Mar20_1 <- as.Date("2020-03-13")

DATE_Mar20_1 - 12
DATE_Mar20_1 + 12

DATE_Mar201 <- as.Date("2020-03-01")
DATE_Mar202 <- as.Date("2020-03-25")

Mar13_20 <- myfunc(DATE_Mar201,DATE_Mar202)

################################
# March 2020
################################
Mar2020_sum_tag_1 =
  Mar13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Mar13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Mar13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Mar13_20 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)

#calculate kernel density estimate for all fish
kud_Mar13_20 <- kernelUD(burb_Mar13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_Mar13_20)

#calculate kernel area estimates
N_Mar13_20 <- kernel.area(kud_Mar13_20, percent = seq(50,95, by = 5), unout = "km2")
N_Mar13_20 <- as.data.frame(t(N_Mar13_20))
colnames(N_Mar13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Mar13_20 <- cbind(fish,N_Mar13_20)
HR_Burbs_Mar13_20 <- HR_Burbs_Mar13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Mar13_20,"BurbHR_Mar13_20.csv")
################################################################################
#get area HR for each percentage
burb.Mar13_20.95.kde <- getverticeshr(kud_Mar13_20, percent = 95, unin = "m", unout = "km2")
burb.Mar13_20.90.kde <- getverticeshr(kud_Mar13_20, percent = 90, unin = "m", unout = "km2")
burb.Mar13_20.85.kde <- getverticeshr(kud_Mar13_20, percent = 85, unin = "m", unout = "km2")
burb.Mar13_20.50.kde <- getverticeshr(kud_Mar13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Mar13_20.95.kde)

HR95_Mar13_20 <- HR_CLIP(burb.Mar13_20.95.kde)
HR90_Mar13_20 <- HR_CLIP(burb.Mar13_20.90.kde)
HR85_Mar13_20 <- HR_CLIP(burb.Mar13_20.85.kde)
HR50_Mar13_20 <- HR_CLIP(burb.Mar13_20.50.kde)

study_HR_Mar13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Mar13_20,HR90_Mar13_20,HR85_Mar13_20,HR50_Mar13_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Mar13_20 <- merge(fish,study_HR_Mar13_20, by = "TRANSMITTER")
study_HR_Mar13_20 <- study_HR_Mar13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Mar13_20,HR85_Mar13_20,HR90_Mar13_20,HR95_Mar13_20)

#write.csv(study_HR_Mar13_20,"BurbHR_Mar13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Apr20_1 <- as.Date("2020-04-13")

DATE_Apr20_1 - 12
DATE_Apr20_1 + 12

DATE_Apr201 <- as.Date("2020-04-01")
DATE_Apr202 <- as.Date("2020-04-25")

Apr13_20 <- myfunc(DATE_Apr201,DATE_Apr202)

################################
# April 2020
################################
Apr2020_sum_tag_1 =
  Apr13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Apr13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Apr13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Apr13_20 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)

#calculate kernel density estimate for all fish
kud_Apr13_20 <- kernelUD(burb_Apr13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_Apr13_20)

#calculate kernel area estimates
N_Apr13_20 <- kernel.area(kud_Apr13_20, percent = seq(50,95, by = 5), unout = "km2")
N_Apr13_20 <- as.data.frame(t(N_Apr13_20))
colnames(N_Apr13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Apr13_20 <- cbind(fish,N_Apr13_20)
HR_Burbs_Apr13_20 <- HR_Burbs_Apr13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Apr13_20,"BurbHR_Apr13_20.csv")
################################################################################
#get area HR for each percentage
burb.Apr13_20.95.kde <- getverticeshr(kud_Apr13_20, percent = 95, unin = "m", unout = "km2")
burb.Apr13_20.90.kde <- getverticeshr(kud_Apr13_20, percent = 90, unin = "m", unout = "km2")
burb.Apr13_20.85.kde <- getverticeshr(kud_Apr13_20, percent = 85, unin = "m", unout = "km2")
burb.Apr13_20.50.kde <- getverticeshr(kud_Apr13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Apr13_20.95.kde)

HR95_Apr13_20 <- HR_CLIP(burb.Apr13_20.95.kde)
HR90_Apr13_20 <- HR_CLIP(burb.Apr13_20.90.kde)
HR85_Apr13_20 <- HR_CLIP(burb.Apr13_20.85.kde)
HR50_Apr13_20 <- HR_CLIP(burb.Apr13_20.50.kde)

study_HR_Apr13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Apr13_20,HR90_Apr13_20,HR85_Apr13_20,HR50_Apr13_20))

#condense all data in to file with HR estimates and fish info
study_HR_Apr13_20 <- merge(fish,study_HR_Apr13_20, by = "TRANSMITTER")
study_HR_Apr13_20 <- study_HR_Apr13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Apr13_20,HR85_Apr13_20,HR90_Apr13_20,HR95_Apr13_20)

#write.csv(study_HR_Apr13_20,"BurbHR_Apr13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_May20_1 <- as.Date("2020-05-13")

DATE_May20_1 - 12
DATE_May20_1 + 12

DATE_May201 <- as.Date("2020-05-01")
DATE_May202 <- as.Date("2020-05-25")

May13_20 <- myfunc(DATE_May201,DATE_May202)

################################
# May 2020
################################
May2020_sum_tag_1 =
  May13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- May13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(May13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_May13_20 <- SpatialPointsDataFrame(coords      = coords,
                                          data        = data, 
                                          proj4string = crs)

#calculate kernel density estimate for all fish
kud_May13_20 <- kernelUD(burb_May13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_May13_20)

#calculate kernel area estimates
N_May13_20 <- kernel.area(kud_May13_20, percent = seq(50,95, by = 5), unout = "km2")
N_May13_20 <- as.data.frame(t(N_May13_20))
colnames(N_May13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_May13_20 <- cbind(fish,N_May13_20)
HR_Burbs_May13_20 <- HR_Burbs_May13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_May13_20,"BurbHR_May13_20.csv")
################################################################################
#get area HR for each percentage
burb.May13_20.95.kde <- getverticeshr(kud_May13_20, percent = 95, unin = "m", unout = "km2")
burb.May13_20.90.kde <- getverticeshr(kud_May13_20, percent = 90, unin = "m", unout = "km2")
burb.May13_20.85.kde <- getverticeshr(kud_May13_20, percent = 85, unin = "m", unout = "km2")
burb.May13_20.50.kde <- getverticeshr(kud_May13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.May13_20.95.kde)

HR95_May13_20 <- HR_CLIP(burb.May13_20.95.kde)
HR90_May13_20 <- HR_CLIP(burb.May13_20.90.kde)
HR85_May13_20 <- HR_CLIP(burb.May13_20.85.kde)
HR50_May13_20 <- HR_CLIP(burb.May13_20.50.kde)

study_HR_May13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_May13_20,HR90_May13_20,HR85_May13_20,HR50_May13_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_May13_20 <- merge(fish,study_HR_May13_20, by = "TRANSMITTER")
study_HR_May13_20 <- study_HR_May13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_May13_20,HR85_May13_20,HR90_May13_20,HR95_May13_20)

write.csv(study_HR_May13_20,"BurbHR_May13_20_TRIMMED.csv")
################################################################################

################################################################################
################################################################################
############### 25 day span home range across study period #####################
################################################################################
################################################################################

HR_25day <- BurbLocsHR

# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Apr19_2 <- as.Date("2019-04-25")

DATE_Apr19_2 - 12
DATE_Apr19_2 + 12

DATE_Apr193 <- as.Date("2019-04-01")
DATE_Apr194 <- as.Date("2019-04-25")

Apr13_19 <- myfunc(DATE_Apr191,DATE_Apr192)    

################################
# April 2019
################################
Apr2019_sum_tag_1 =
  Apr13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Apr13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Apr13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Apr13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Apr13_19)
#calculate kernel density estimate for all fish
kud_Apr13_19 <- kernelUD(burb_Apr13_19[,1], h="href", grid = xy)
#rm(kud_Apr13_19)

#print all fish estimates
image(kud_Apr13_19)

#calculate kernel area estimates
N_Apr13_19 <- kernel.area(kud_Apr13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Apr13_19 <- as.data.frame(t(N_Apr13_19))
colnames(N_Apr13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Apr13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "504mm_F_a_t",]
fish <- fish[!fish$TRANSMITTER == "625mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Apr13_19 <- cbind(fish,N_Apr13_19)
HR_Burbs_Apr13_19 <- HR_Burbs_Apr13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Apr13_19,"BurbHR_Apr13_19.csv")

################################################################################
burb.Apr13_19.95.kde <- getverticeshr(kud_Apr13_19, percent = 95, unin = "m", unout = "km2")
burb.Apr13_19.90.kde <- getverticeshr(kud_Apr13_19, percent = 90, unin = "m", unout = "km2")
burb.Apr13_19.85.kde <- getverticeshr(kud_Apr13_19, percent = 85, unin = "m", unout = "km2")
burb.Apr13_19.50.kde <- getverticeshr(kud_Apr13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Apr13_19.50.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Apr13_19.95.kde)

HR95_Apr13_19 <- HR_CLIP(burb.Apr13_19.95.kde)
HR90_Apr13_19 <- HR_CLIP(burb.Apr13_19.90.kde)
HR85_Apr13_19 <- HR_CLIP(burb.Apr13_19.85.kde)
HR50_Apr13_19 <- HR_CLIP(burb.Apr13_19.50.kde)
#rm(HR95_Apr13_19)

study_HR_Apr13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Apr13_19,HR90_Apr13_19,HR85_Apr13_19,HR50_Apr13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Apr13_19 <- merge(fish,study_HR_Apr13_19, by = "TRANSMITTER")
study_HR_Apr13_19 <- study_HR_Apr13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Apr13_19,HR85_Apr13_19,HR90_Apr13_19,HR95_Apr13_19)

#write.csv(study_HR_Apr13_19,"BurbHR_Apr13_19_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_May19_1 <- as.Date("2019-05-13")

DATE_May19_1 - 12
DATE_May19_1 + 12

DATE_May191 <- as.Date("2019-05-01")
DATE_May192 <- as.Date("2019-05-25")

May13_19 <- myfunc(DATE_May191,DATE_May192)

################################
# May 2019
################################
May2019_sum_tag_1 =
  May13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- May13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(May13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_May13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_May13_19)
#calculate kernel density estimate for all fish
kud_May13_19 <- kernelUD(burb_May13_19[,1], h="href", grid = xy)
#rm(kud_May13_19)
#print all fish estimates
image(kud_May13_19)

#calculate kernel area estimates
N_May13_19 <- kernel.area(kud_May13_19, percent = seq(50,95, by = 5), unout = "km2")
N_May13_19 <- as.data.frame(t(N_May13_19))
colnames(N_May13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_May13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_May13_19 <- cbind(fish,N_May13_19)
HR_Burbs_May13_19 <- HR_Burbs_May13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_May13_19,"BurbHR_May13_19.csv")
################################################################################
#get area HR for each percentage
burb.May13_19.95.kde <- getverticeshr(kud_May13_19, percent = 95, unin = "m", unout = "km2")
burb.May13_19.90.kde <- getverticeshr(kud_May13_19, percent = 90, unin = "m", unout = "km2")
burb.May13_19.85.kde <- getverticeshr(kud_May13_19, percent = 85, unin = "m", unout = "km2")
burb.May13_19.50.kde <- getverticeshr(kud_May13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.May13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.May13_19.95.kde)

HR95_May13_19 <- HR_CLIP(burb.May13_19.95.kde)
HR90_May13_19 <- HR_CLIP(burb.May13_19.90.kde)
HR85_May13_19 <- HR_CLIP(burb.May13_19.85.kde)
HR50_May13_19 <- HR_CLIP(burb.May13_19.50.kde)
#rm(HR95_May13_19)

study_HR_May13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_May13_19,HR90_May13_19,HR85_May13_19,HR50_May13_19))

#condense all data in to file with HR estimates and fish info
study_HR_May13_19 <- merge(fish,study_HR_May13_19, by = "TRANSMITTER")
study_HR_May13_19 <- study_HR_May13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_May13_19,HR85_May13_19,HR90_May13_19,HR95_May13_19)

#write.csv(study_HR_May13_19,"BurbHR_May13_19_TRIMMED.csv")
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Jun19_1 <- as.Date("2019-06-13")

DATE_Jun19_1 - 12
DATE_Jun19_1 + 12

DATE_Jun191 <- as.Date("2019-06-01")
DATE_Jun192 <- as.Date("2019-06-25")

Jun13_19 <- myfunc(DATE_Jun191,DATE_Jun192)

################################
# June 2019
################################
Jun2019_sum_tag_1 =
  Jun13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

Jun13_19 <- Jun13_19[!Jun13_19$TRANSMITTER == "492mm_F_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Jun13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Jun13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Jun13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Jun13_19)

#calculate kernel density estimate for all fish
kud_Jun13_19 <- kernelUD(burb_Jun13_19[,1], h="href", grid = xy)
#rm(kud_Jun13_19)

#print all fish estimates
image(kud_Jun13_19)

#calculate kernel area estimates
N_Jun13_19 <- kernel.area(kud_Jun13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Jun13_19 <- as.data.frame(t(N_Jun13_19))
colnames(N_Jun13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Jun13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Jun13_19 <- cbind(fish,N_Jun13_19)
HR_Burbs_Jun13_19 <- HR_Burbs_Jun13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Jun13_19,"BurbHR_Jun13_19.csv")
################################################################################
#get area HR for each percentage
burb.Jun13_19.95.kde <- getverticeshr(kud_Jun13_19, percent = 95, unin = "m", unout = "km2")
burb.Jun13_19.90.kde <- getverticeshr(kud_Jun13_19, percent = 90, unin = "m", unout = "km2")
burb.Jun13_19.85.kde <- getverticeshr(kud_Jun13_19, percent = 85, unin = "m", unout = "km2")
burb.Jun13_19.50.kde <- getverticeshr(kud_Jun13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Jun13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Jun13_19.95.kde)

HR95_Jun13_19 <- HR_CLIP(burb.Jun13_19.95.kde)
HR90_Jun13_19 <- HR_CLIP(burb.Jun13_19.90.kde)
HR85_Jun13_19 <- HR_CLIP(burb.Jun13_19.85.kde)
HR50_Jun13_19 <- HR_CLIP(burb.Jun13_19.50.kde)
#rm(HR95_Jun13_19)

study_HR_Jun13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Jun13_19,HR90_Jun13_19,HR85_Jun13_19,HR50_Jun13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Jun13_19 <- merge(fish,study_HR_Jun13_19, by = "TRANSMITTER")
study_HR_Jun13_19 <- study_HR_Jun13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Jun13_19,HR85_Jun13_19,HR90_Jun13_19,HR95_Jun13_19)

#write.csv(study_HR_Jun13_19,"BurbHR_Jun13_19_TRIMMED.csv")

################################################################################

# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Jul19_1 <- as.Date("2019-07-13")

DATE_Jul19_1 - 12
DATE_Jul19_1 + 12

DATE_Jul191 <- as.Date("2019-07-01")
DATE_Jul192 <- as.Date("2019-07-25")

Jul13_19 <- myfunc(DATE_Jul191,DATE_Jul192)

################################
# July 2019
################################
Jul2019_sum_tag_1 =
  Jul13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish with less than 30 detections
Jul13_19 <- Jul13_19[!Jul13_19$TRANSMITTER == "492mm_F_a_t",]
Jul13_19 <- Jul13_19[!Jul13_19$TRANSMITTER == "527mm_M_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Jul13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Jul13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Jul13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Jul13_19)

#calculate kernel density estimate for all fish
kud_Jul13_19 <- kernelUD(burb_Jul13_19[,1], h="href", grid = xy)
#rm(kud_Jul13_19)

#print all fish estimates
image(kud_Jul13_19)

#calculate kernel area estimates
N_Jul13_19 <- kernel.area(kud_Jul13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Jul13_19 <- as.data.frame(t(N_Jul13_19))
colnames(N_Jul13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Jul13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]
fish <- fish[!fish$TRANSMITTER == "527mm_M_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Jul13_19 <- cbind(fish,N_Jul13_19)
HR_Burbs_Jul13_19 <- HR_Burbs_Jul13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Jul13_19,"BurbHR_Jul13_19.csv")
################################################################################
#get area HR for each percentage
burb.Jul13_19.95.kde <- getverticeshr(kud_Jul13_19, percent = 95, unin = "m", unout = "km2")
burb.Jul13_19.90.kde <- getverticeshr(kud_Jul13_19, percent = 90, unin = "m", unout = "km2")
burb.Jul13_19.85.kde <- getverticeshr(kud_Jul13_19, percent = 85, unin = "m", unout = "km2")
burb.Jul13_19.50.kde <- getverticeshr(kud_Jul13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Jul13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Jul13_19.95.kde)

HR95_Jul13_19 <- HR_CLIP(burb.Jul13_19.95.kde)
HR90_Jul13_19 <- HR_CLIP(burb.Jul13_19.90.kde)
HR85_Jul13_19 <- HR_CLIP(burb.Jul13_19.85.kde)
HR50_Jul13_19 <- HR_CLIP(burb.Jul13_19.50.kde)
#rm(HR95_Jul13_19)

study_HR_Jul13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Jul13_19,HR90_Jul13_19,HR85_Jul13_19,HR50_Jul13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Jul13_19 <- merge(fish,study_HR_Jul13_19, by = "TRANSMITTER")
study_HR_Jul13_19 <- study_HR_Jul13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Jul13_19,HR85_Jul13_19,HR90_Jul13_19,HR95_Jul13_19)

#write.csv(study_HR_Jul13_19,"BurbHR_Jul13_19_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Aug19_1 <- as.Date("2019-08-13")

DATE_Aug19_1 - 12
DATE_Aug19_1 + 12

DATE_Aug191 <- as.Date("2019-08-01")
DATE_Aug192 <- as.Date("2019-08-25")

Aug13_19 <- myfunc(DATE_Aug191,DATE_Aug192)

################################
# August 2019
################################
Aug2019_sum_tag_1 =
  Aug13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish with less than 30 detections
Aug13_19 <- Aug13_19[!Aug13_19$TRANSMITTER == "492mm_F_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Aug13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Aug13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Aug13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Aug13_19)

#calculate kernel density estimate for all fish
kud_Aug13_19 <- kernelUD(burb_Aug13_19[,1], h="href", grid = xy)
#rm(kud_Aug13_19)

#print all fish estimates
image(kud_Aug13_19)

#calculate kernel area estimates
N_Aug13_19 <- kernel.area(kud_Aug13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Aug13_19 <- as.data.frame(t(N_Aug13_19))
colnames(N_Aug13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Aug13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "492mm_F_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Aug13_19 <- cbind(fish,N_Aug13_19)
HR_Burbs_Aug13_19 <- HR_Burbs_Aug13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Aug13_19,"BurbHR_Aug13_19.csv")
################################################################################
#get area HR for each percentage
burb.Aug13_19.95.kde <- getverticeshr(kud_Aug13_19, percent = 93, unin = "m", unout = "km2")
#lack of new locations resulting in failure of 95 kud
burb.Aug13_19.90.kde <- getverticeshr(kud_Aug13_19, percent = 90, unin = "m", unout = "km2")
burb.Aug13_19.85.kde <- getverticeshr(kud_Aug13_19, percent = 85, unin = "m", unout = "km2")
burb.Aug13_19.50.kde <- getverticeshr(kud_Aug13_19, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Aug13_19.90.kde)

HR95_Aug13_19 <- HR_CLIP(burb.Aug13_19.95.kde)
HR90_Aug13_19 <- HR_CLIP(burb.Aug13_19.90.kde)
HR85_Aug13_19 <- HR_CLIP(burb.Aug13_19.85.kde)
HR50_Aug13_19 <- HR_CLIP(burb.Aug13_19.50.kde)

study_HR_Aug13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Aug13_19,HR90_Aug13_19,HR85_Aug13_19,HR50_Aug13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Aug13_19 <- merge(fish,study_HR_Aug13_19, by = "TRANSMITTER")
study_HR_Aug13_19 <- study_HR_Aug13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Aug13_19,HR85_Aug13_19,HR90_Aug13_19,HR95_Aug13_19)

#write.csv(study_HR_Aug13_19,"BurbHR_Aug13_19_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Sep19_1 <- as.Date("2019-09-13")

DATE_Sep19_1 - 12
DATE_Sep19_1 + 12

DATE_Sep191 <- as.Date("2019-09-01")
DATE_Sep192 <- as.Date("2019-09-25")

Sep13_19 <- myfunc(DATE_Sep191,DATE_Sep192)

################################
# September 2019
################################
Sep2019_sum_tag_1 =
  Sep13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

#remove fish that results in KUD error
Sep13_19 <- Sep13_19[!Sep13_19$TRANSMITTER == "566mm_U_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Sep13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Sep13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Sep13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Sep13_19)

#calculate kernel density estimate for all fish
kud_Sep13_19 <- kernelUD(burb_Sep13_19[,1], h="href", grid = xy)
#rm(kud_Sep13_19)

#print all fish estimates
image(kud_Sep13_19)

#calculate kernel area estimates
N_Sep13_19 <- kernel.area(kud_Sep13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Sep13_19 <- as.data.frame(t(N_Sep13_19))
colnames(N_Sep13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Sep13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "566mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Sep13_19 <- cbind(fish,N_Sep13_19)
HR_Burbs_Sep13_19 <- HR_Burbs_Sep13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Sep13_19,"BurbHR_Sep13_20.csv")

################################################################################
#get area HR for each percentage
burb.Sep13_19.95.kde <- getverticeshr(kud_Sep13_19, percent = 95, unin = "m", unout = "km2")
burb.Sep13_19.90.kde <- getverticeshr(kud_Sep13_19, percent = 90, unin = "m", unout = "km2")
burb.Sep13_19.85.kde <- getverticeshr(kud_Sep13_19, percent = 85, unin = "m", unout = "km2")
burb.Sep13_19.50.kde <- getverticeshr(kud_Sep13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Sep13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Sep13_19.95.kde)

HR95_Sep13_19 <- HR_CLIP(burb.Sep13_19.95.kde)
HR90_Sep13_19 <- HR_CLIP(burb.Sep13_19.90.kde)
HR85_Sep13_19 <- HR_CLIP(burb.Sep13_19.85.kde)
HR50_Sep13_19 <- HR_CLIP(burb.Sep13_19.50.kde)
#rm(HR95_Sep13_19)

study_HR_Sep13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Sep13_19,HR90_Sep13_19,HR85_Sep13_19,HR50_Sep13_19))

#condense all data in to file with HR estimates and fish info
study_HR_Sep13_19 <- merge(fish,study_HR_Sep13_19, by = "TRANSMITTER")
study_HR_Sep13_19 <- study_HR_Sep13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Sep13_19,HR85_Sep13_19,HR90_Sep13_19,HR95_Sep13_19)

#write.csv(study_HR_Sep13_19,"BurbHR_Sep13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Oct19_1 <- as.Date("2019-10-13")

DATE_Oct19_1 - 12
DATE_Oct19_1 + 12

DATE_Oct191 <- as.Date("2019-10-01")
DATE_Oct192 <- as.Date("2019-10-25")

Oct13_19 <- myfunc(DATE_Oct191,DATE_Oct192)

################################
# October 2019
################################
Oct2019_sum_tag_1 =
  Oct13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Oct13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Oct13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Oct13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Oct13_19)

#calculate kernel density estimate for all fish
kud_Oct13_19 <- kernelUD(burb_Oct13_19[,1], h="href", grid = xy)
#rm(kud_Oct13_19)

#print all fish estimates
image(kud_Oct13_19)

#calculate kernel area estimates
N_Oct13_19 <- kernel.area(kud_Oct13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Oct13_19 <- as.data.frame(t(N_Oct13_19))
colnames(N_Oct13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Oct13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Oct13_19 <- cbind(fish,N_Oct13_19)
HR_Burbs_Oct13_19 <- HR_Burbs_Oct13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Oct13_19,"BurbHR_Oct13_20.csv")
################################################################################
#get area HR for each percentage
burb.Oct13_19.95.kde <- getverticeshr(kud_Oct13_19, percent = 95, unin = "m", unout = "km2")
burb.Oct13_19.90.kde <- getverticeshr(kud_Oct13_19, percent = 90, unin = "m", unout = "km2")
burb.Oct13_19.85.kde <- getverticeshr(kud_Oct13_19, percent = 85, unin = "m", unout = "km2")
burb.Oct13_19.50.kde <- getverticeshr(kud_Oct13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Oct13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Oct13_19.95.kde)

HR95_Oct13_19 <- HR_CLIP(burb.Oct13_19.95.kde)
HR90_Oct13_19 <- HR_CLIP(burb.Oct13_19.90.kde)
HR85_Oct13_19 <- HR_CLIP(burb.Oct13_19.85.kde)
HR50_Oct13_19 <- HR_CLIP(burb.Oct13_19.50.kde)
#rm(HR95_Oct13_19)

study_HR_Oct13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Oct13_19,HR90_Oct13_19,HR85_Oct13_19,HR50_Oct13_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Oct13_19 <- merge(fish,study_HR_Oct13_19, by = "TRANSMITTER")
study_HR_Oct13_19 <- study_HR_Oct13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Oct13_19,HR85_Oct13_19,HR90_Oct13_19,HR95_Oct13_19)

#write.csv(study_HR_Oct13_19,"BurbHR_Oct13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Nov19_1 <- as.Date("2019-11-13")

DATE_Nov19_1 - 12
DATE_Nov19_1 + 12

DATE_Nov191 <- as.Date("2019-11-01")
DATE_Nov192 <- as.Date("2019-11-25")

Nov13_19 <- myfunc(DATE_Nov191,DATE_Nov192)

################################
# November 2019
################################
Nov2019_sum_tag_1 =
  Nov13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Nov13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Nov13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Nov13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Nov13_19)

#calculate kernel density estimate for all fish
kud_Nov13_19 <- kernelUD(burb_Nov13_19[,1], h="href", grid = xy)
#rm(kud_Nov13_19)

#print all fish estimates
image(kud_Nov13_19)

#calculate kernel area estimates
N_Nov13_19 <- kernel.area(kud_Nov13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Nov13_19 <- as.data.frame(t(N_Nov13_19))
colnames(N_Nov13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Nov13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Nov13_19 <- cbind(fish,N_Nov13_19)
HR_Burbs_Nov13_19 <- HR_Burbs_Nov13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Nov13_19,"BurbHR_Nov13_20.csv")
################################################################################
#get area HR for each percentage
burb.Nov13_19.95.kde <- getverticeshr(kud_Nov13_19, percent = 95, unin = "m", unout = "km2")
burb.Nov13_19.90.kde <- getverticeshr(kud_Nov13_19, percent = 90, unin = "m", unout = "km2")
burb.Nov13_19.85.kde <- getverticeshr(kud_Nov13_19, percent = 85, unin = "m", unout = "km2")
burb.Nov13_19.50.kde <- getverticeshr(kud_Nov13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Nov13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Nov13_19.95.kde)

HR95_Nov13_19 <- HR_CLIP(burb.Nov13_19.95.kde)
HR90_Nov13_19 <- HR_CLIP(burb.Nov13_19.90.kde)
HR85_Nov13_19 <- HR_CLIP(burb.Nov13_19.85.kde)
HR50_Nov13_19 <- HR_CLIP(burb.Nov13_19.50.kde)
#rm(HR95_Nov13_19)

study_HR_Nov13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Nov13_19,HR90_Nov13_19,HR85_Nov13_19,HR50_Nov13_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Nov13_19 <- merge(fish,study_HR_Nov13_19, by = "TRANSMITTER")
study_HR_Nov13_19 <- study_HR_Nov13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Nov13_19,HR85_Nov13_19,HR90_Nov13_19,HR95_Nov13_19)

#write.csv(study_HR_Nov13_19,"BurbHR_Nov13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Dec19_1 <- as.Date("2019-12-13")

DATE_Dec19_1 - 12
DATE_Dec19_1 + 12

DATE_Dec191 <- as.Date("2019-12-01")
DATE_Dec192 <- as.Date("2019-12-25")

Dec13_19 <- myfunc(DATE_Dec191,DATE_Dec192)

################################
# December 2019
################################
Dec2019_sum_tag_1 =
  Dec13_19 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Dec13_19[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Dec13_19[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Dec13_19 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)
#rm(burb_Dec13_19)

#calculate kernel density estimate for all fish
kud_Dec13_19 <- kernelUD(burb_Dec13_19[,1], h="href", grid = xy)
#rm(kud_Dec13_19)

#print all fish estimates
image(kud_Dec13_19)

#calculate kernel area estimates
N_Dec13_19 <- kernel.area(kud_Dec13_19, percent = seq(50,95, by = 5), unout = "km2")
N_Dec13_19 <- as.data.frame(t(N_Dec13_19))
colnames(N_Dec13_19) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")
#rm(N_Dec13_19)

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Dec13_19 <- cbind(fish,N_Dec13_19)
HR_Burbs_Dec13_19 <- HR_Burbs_Dec13_19 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Dec13_19,"BurbHR_Dec13_20.csv")
################################################################################
#get area HR for each percentage
burb.Dec13_19.95.kde <- getverticeshr(kud_Dec13_19, percent = 95, unin = "m", unout = "km2")
burb.Dec13_19.90.kde <- getverticeshr(kud_Dec13_19, percent = 90, unin = "m", unout = "km2")
burb.Dec13_19.85.kde <- getverticeshr(kud_Dec13_19, percent = 85, unin = "m", unout = "km2")
burb.Dec13_19.50.kde <- getverticeshr(kud_Dec13_19, percent = 50, unin = "m", unout = "km2")
#rm(burb.Dec13_19.95.kde)

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf
#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Dec13_19.95.kde)

HR95_Dec13_19 <- HR_CLIP(burb.Dec13_19.95.kde)
HR90_Dec13_19 <- HR_CLIP(burb.Dec13_19.90.kde)
HR85_Dec13_19 <- HR_CLIP(burb.Dec13_19.85.kde)
HR50_Dec13_19 <- HR_CLIP(burb.Dec13_19.50.kde)
#rm(HR95_Dec13_19)

study_HR_Dec13_19 <- as.data.frame(cbind(TRANSMITTER,HR95_Dec13_19,HR90_Dec13_19,HR85_Dec13_19,HR50_Dec13_19))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Dec13_19 <- merge(fish,study_HR_Dec13_19, by = "TRANSMITTER")
study_HR_Dec13_19 <- study_HR_Dec13_19 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Dec13_19,HR85_Dec13_19,HR90_Dec13_19,HR95_Dec13_19)

#write.csv(study_HR_Dec13_19,"BurbHR_Dec13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Jan20_1 <- as.Date("2020-01-13")

DATE_Jan20_1 - 12
DATE_Jan20_1 + 12

DATE_Jan201 <- as.Date("2020-01-01")
DATE_Jan202 <- as.Date("2020-01-25")

Jan13_20 <- myfunc(DATE_Jan201,DATE_Jan202)

################################
# January 2020
################################
Jan2020_sum_tag_1 =
  Jan13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Jan13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Jan13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Jan13_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)

#calculate kernel density estimate for all fish
kud_Jan13_20 <- kernelUD(burb_Jan13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_Jan13_20)

#calculate kernel area estimates
N_Jan13_20 <- kernel.area(kud_Jan13_20, percent = seq(50,95, by = 5), unout = "km2")
N_Jan13_20 <- as.data.frame(t(N_Jan13_20))
colnames(N_Jan13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Jan13_20 <- cbind(fish,N_Jan13_20)
HR_Burbs_Jan13_20 <- HR_Burbs_Jan13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Jan13_20,"BurbHR_Jan13_20.csv")
################################################################################
#get area HR for each percentage
burb.Jan13_20.95.kde <- getverticeshr(kud_Jan13_20, percent = 95, unin = "m", unout = "km2")
burb.Jan13_20.90.kde <- getverticeshr(kud_Jan13_20, percent = 90, unin = "m", unout = "km2")
burb.Jan13_20.85.kde <- getverticeshr(kud_Jan13_20, percent = 85, unin = "m", unout = "km2")
burb.Jan13_20.50.kde <- getverticeshr(kud_Jan13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Jan13_20.95.kde)

HR95_Jan13_20 <- HR_CLIP(burb.Jan13_20.95.kde)
HR90_Jan13_20 <- HR_CLIP(burb.Jan13_20.90.kde)
HR85_Jan13_20 <- HR_CLIP(burb.Jan13_20.85.kde)
HR50_Jan13_20 <- HR_CLIP(burb.Jan13_20.50.kde)

study_HR_Jan13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Jan13_20,HR90_Jan13_20,HR85_Jan13_20,HR50_Jan13_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Jan13_20 <- merge(fish,study_HR_Jan13_20, by = "TRANSMITTER")
study_HR_Jan13_20 <- study_HR_Jan13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Jan13_20,HR85_Jan13_20,HR90_Jan13_20,HR95_Jan13_20)

#write.csv(study_HR_Jan13_20,"BurbHR_Jan13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Feb20_1 <- as.Date("2020-02-13")

DATE_Feb20_1 - 12
DATE_Feb20_1 + 12

DATE_Feb201 <- as.Date("2020-02-01")
DATE_Feb202 <- as.Date("2020-02-25")

Feb13_20 <- myfunc(DATE_Feb201,DATE_Feb202)

################################
# February 2020
################################
Feb2019_sum_tag_1 =
  Feb13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# remove fish with less than 30 locations
Feb13_20 <- Feb13_20[!Feb13_20$TRANSMITTER == "524mm_U_a_t",]

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Feb13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Feb13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Feb13_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)

#calculate kernel density estimate for all fish
kud_Feb13_20 <- kernelUD(burb_Feb13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_Feb13_20)

#calculate kernel area estimates
N_Feb13_20 <- kernel.area(kud_Feb13_20, percent = seq(50,95, by = 5), unout = "km2")
N_Feb13_20 <- as.data.frame(t(N_Feb13_20))
colnames(N_Feb13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")
fish <- fish[!fish$TRANSMITTER == "524mm_U_a_t",]

#condense all data in to file with HR estimates and fish info
HR_Burbs_Feb13_20 <- cbind(fish,N_Feb13_20)
HR_Burbs_Feb13_20 <- HR_Burbs_Feb13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Feb13_20,"BurbHR_Feb13_20.csv")
################################################################################
#get area HR for each percentage
burb.Feb13_20.95.kde <- getverticeshr(kud_Feb13_20, percent = 95, unin = "m", unout = "km2")
burb.Feb13_20.90.kde <- getverticeshr(kud_Feb13_20, percent = 90, unin = "m", unout = "km2")
burb.Feb13_20.85.kde <- getverticeshr(kud_Feb13_20, percent = 85, unin = "m", unout = "km2")
burb.Feb13_20.50.kde <- getverticeshr(kud_Feb13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Feb13_20.95.kde)

HR95_Feb13_20 <- HR_CLIP(burb.Feb13_20.95.kde)
HR90_Feb13_20 <- HR_CLIP(burb.Feb13_20.90.kde)
HR85_Feb13_20 <- HR_CLIP(burb.Feb13_20.85.kde)
HR50_Feb13_20 <- HR_CLIP(burb.Feb13_20.50.kde)

study_HR_Feb13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Feb13_20,HR90_Feb13_20,HR85_Feb13_20,HR50_Feb13_20))


#condense all data in to file with HR estimates and fish info
study_HR_Feb13_20 <- merge(fish,study_HR_Feb13_20, by = "TRANSMITTER")
study_HR_Feb13_20 <- study_HR_Feb13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Feb13_20,HR85_Feb13_20,HR90_Feb13_20,HR95_Feb13_20)

#write.csv(study_HR_Feb13_20,"BurbHR_Feb13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Mar20_1 <- as.Date("2020-03-13")

DATE_Mar20_1 - 12
DATE_Mar20_1 + 12

DATE_Mar201 <- as.Date("2020-03-01")
DATE_Mar202 <- as.Date("2020-03-25")

Mar13_20 <- myfunc(DATE_Mar201,DATE_Mar202)

################################
# March 2020
################################
May2019_sum_tag_1 =
  Mar13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Mar13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Mar13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Mar13_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)

#calculate kernel density estimate for all fish
kud_Mar13_20 <- kernelUD(burb_Mar13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_Mar13_20)

#calculate kernel area estimates
N_Mar13_20 <- kernel.area(kud_Mar13_20, percent = seq(50,95, by = 5), unout = "km2")
N_Mar13_20 <- as.data.frame(t(N_Mar13_20))
colnames(N_Mar13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Mar13_20 <- cbind(fish,N_Mar13_20)
HR_Burbs_Mar13_20 <- HR_Burbs_Mar13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Mar13_20,"BurbHR_Mar13_20.csv")
################################################################################
#get area HR for each percentage
burb.Mar13_20.95.kde <- getverticeshr(kud_Mar13_20, percent = 95, unin = "m", unout = "km2")
burb.Mar13_20.90.kde <- getverticeshr(kud_Mar13_20, percent = 90, unin = "m", unout = "km2")
burb.Mar13_20.85.kde <- getverticeshr(kud_Mar13_20, percent = 85, unin = "m", unout = "km2")
burb.Mar13_20.50.kde <- getverticeshr(kud_Mar13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Mar13_20.95.kde)

HR95_Mar13_20 <- HR_CLIP(burb.Mar13_20.95.kde)
HR90_Mar13_20 <- HR_CLIP(burb.Mar13_20.90.kde)
HR85_Mar13_20 <- HR_CLIP(burb.Mar13_20.85.kde)
HR50_Mar13_20 <- HR_CLIP(burb.Mar13_20.50.kde)

study_HR_Mar13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Mar13_20,HR90_Mar13_20,HR85_Mar13_20,HR50_Mar13_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_Mar13_20 <- merge(fish,study_HR_Mar13_20, by = "TRANSMITTER")
study_HR_Mar13_20 <- study_HR_Mar13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Mar13_20,HR85_Mar13_20,HR90_Mar13_20,HR95_Mar13_20)

#write.csv(study_HR_Mar13_20,"BurbHR_Mar13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_Apr20_1 <- as.Date("2020-04-13")

DATE_Apr20_1 - 12
DATE_Apr20_1 + 12

DATE_Apr201 <- as.Date("2020-04-01")
DATE_Apr202 <- as.Date("2020-04-25")

Apr13_20 <- myfunc(DATE_Apr201,DATE_Apr202)

################################
# April 2020
################################
Apr2020_sum_tag_1 =
  Apr13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- Apr13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(Apr13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_Apr13_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)

#calculate kernel density estimate for all fish
kud_Apr13_20 <- kernelUD(burb_Apr13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_Apr13_20)

#calculate kernel area estimates
N_Apr13_20 <- kernel.area(kud_Apr13_20, percent = seq(50,95, by = 5), unout = "km2")
N_Apr13_20 <- as.data.frame(t(N_Apr13_20))
colnames(N_Apr13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_Apr13_20 <- cbind(fish,N_Apr13_20)
HR_Burbs_Apr13_20 <- HR_Burbs_Apr13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_Apr13_20,"BurbHR_Apr13_20.csv")
################################################################################
#get area HR for each percentage
burb.Apr13_20.95.kde <- getverticeshr(kud_Apr13_20, percent = 95, unin = "m", unout = "km2")
burb.Apr13_20.90.kde <- getverticeshr(kud_Apr13_20, percent = 90, unin = "m", unout = "km2")
burb.Apr13_20.85.kde <- getverticeshr(kud_Apr13_20, percent = 85, unin = "m", unout = "km2")
burb.Apr13_20.50.kde <- getverticeshr(kud_Apr13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.Apr13_20.95.kde)

HR95_Apr13_20 <- HR_CLIP(burb.Apr13_20.95.kde)
HR90_Apr13_20 <- HR_CLIP(burb.Apr13_20.90.kde)
HR85_Apr13_20 <- HR_CLIP(burb.Apr13_20.85.kde)
HR50_Apr13_20 <- HR_CLIP(burb.Apr13_20.50.kde)

study_HR_Apr13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_Apr13_20,HR90_Apr13_20,HR85_Apr13_20,HR50_Apr13_20))

#condense all data in to file with HR estimates and fish info
study_HR_Apr13_20 <- merge(fish,study_HR_Apr13_20, by = "TRANSMITTER")
study_HR_Apr13_20 <- study_HR_Apr13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_Apr13_20,HR85_Apr13_20,HR90_Apr13_20,HR95_Apr13_20)

#write.csv(study_HR_Apr13_20,"BurbHR_Apr13_20_TRIMMED.csv")
################################################################################
################################################################################
# function for subseting data by 25 days spans on the 25th of the month
myfunc <- function(x,y){HR_25day[HR_25day$datetime >= x & HR_25day$datetime <= y,]}

DATE_May20_1 <- as.Date("2020-05-13")

DATE_May20_1 - 12
DATE_May20_1 + 12

DATE_May201 <- as.Date("2020-05-01")
DATE_May202 <- as.Date("2020-05-25")

May13_20 <- myfunc(DATE_May201,DATE_May202)

################################
# May 2020
################################
May2020_sum_tag_1 =
  May13_20 %>%
  group_by(TRANSMITTER) %>%
  summarise(first = first(datetime), last=last(datetime),n = n())

# prepare coordinates, data, and proj4string
#https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
coords <- May13_20[ , c("EASTING", "NORTHING")]   # coordinates
data   <- as.data.frame(May13_20[ , c("TRANSMITTER","datetime","HPE","TEMP","DEPTH","daytime","SEX","LENGTH","WEIGHT","BASIN")])          # data
crs    <- CRS("+init=epsg:32615") # proj4string of coords

# make the SpatialPointsDataFrame object
burb_May13_20 <- SpatialPointsDataFrame(coords      = coords,
                                        data        = data, 
                                        proj4string = crs)

#calculate kernel density estimate for all fish
kud_May13_20 <- kernelUD(burb_May13_20[,1], h="href", grid = xy)

#print all fish estimates
image(kud_May13_20)

#calculate kernel area estimates
N_May13_20 <- kernel.area(kud_May13_20, percent = seq(50,95, by = 5), unout = "km2")
N_May13_20 <- as.data.frame(t(N_May13_20))
colnames(N_May13_20) <- c("HR_50","HR_55","HR_60","HR_65","HR_70","HR_75","HR_80","HR_85","HR_90","HR_95")


#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
HR_Burbs_May13_20 <- cbind(fish,N_May13_20)
HR_Burbs_May13_20 <- HR_Burbs_May13_20 %>% 
  dplyr::select(TRANSMITTER, SEX, LENGTH, WEIGHT, BASIN,HR_50,HR_85,HR_90,HR_95)

#write.csv(HR_Burbs_May13_20,"BurbHR_May13_20.csv")
################################################################################
#get area HR for each percentage
burb.May13_20.95.kde <- getverticeshr(kud_May13_20, percent = 95, unin = "m", unout = "km2")
burb.May13_20.90.kde <- getverticeshr(kud_May13_20, percent = 90, unin = "m", unout = "km2")
burb.May13_20.85.kde <- getverticeshr(kud_May13_20, percent = 85, unin = "m", unout = "km2")
burb.May13_20.50.kde <- getverticeshr(kud_May13_20, percent = 50, unin = "m", unout = "km2")

#https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/data-expeditions-instructions-Lohmann-Jacoby.pdf

#The Kernel Density Estimates overlap with land!
################################################################################
TRANSMITTER <- row.names(burb.May13_20.95.kde)

HR95_May13_20 <- HR_CLIP(burb.May13_20.95.kde)
HR90_May13_20 <- HR_CLIP(burb.May13_20.90.kde)
HR85_May13_20 <- HR_CLIP(burb.May13_20.85.kde)
HR50_May13_20 <- HR_CLIP(burb.May13_20.50.kde)

study_HR_May13_20 <- as.data.frame(cbind(TRANSMITTER,HR95_May13_20,HR90_May13_20,HR85_May13_20,HR50_May13_20))

#Import fish information
fish <- read_csv("E:/BSUGradBurbot/DATA/FishData/HomeRange_32FishInfo.csv")

#condense all data in to file with HR estimates and fish info
study_HR_May13_20 <- merge(fish,study_HR_May13_20, by = "TRANSMITTER")
study_HR_May13_20 <- study_HR_May13_20 %>% 
  dplyr::select(TRANSMITTER,SEX,LENGTH,WEIGHT,BASIN,HR50_May13_20,HR85_May13_20,HR90_May13_20,HR95_May13_20)

#write.csv(study_HR_May13_20,"BurbHR_May13_20_TRIMMED.csv")

################################################################################
################################################################################
############### END HOME RANGE ESITIMATE BY 25 DAY PERIODS #####################
################################################################################
################################################################################

#####################################################################
#################### OLD CLIPPING ATTEMPT ###########################
########### MIGHT BE ABLE TO GET TO WORK WITH OTHER DATA ############
#####################################################################

################################
library(rgeos)          # For readOGR(); gDifference(); gUnion()
library(rgdal)          # For readOGR()
library(sp)             # For extent(); Polygon(); Polygons(); SpatialPolygons(); coordinates(); over()
library(ggplot2)        # For fortify(); ggplot()
library(gridExtra)      # For grid.arrange()
library(maptools)       # For unionSpatialPolygons(); proj4string()
library(deldir)         # For deldir()
library(rmapshaper)     # For ms_simplify()
library(magrittr)       # 

str(pw.all.95.fkde)
polys <- pw.all.95.fkde@polygons
length(polys)

# Attributes of the second polygon
this.poly <- polys[[29]]
str(this.poly)
# Convert to SpatialPolygons
SP.this.poly <- SpatialPolygons(list(this.poly))
ggplot(SP.this.poly, aes(x = long, y = lat, group = group)) +
  geom_polygon(colour = "black", size = 0.2, fill = "blue") +
  coord_map()

#############
#Clip out anything on land in California
clip <- gDifference(burb.Feb13_20.95.kde, us.pro, byid=TRUE)
row.names(clip) <- gsub(" .*", "", row.names(clip))


all_spdfs_together_cropped <- SpatialPolygonsDataFrame(clip,
                                                       data=burb.Feb13_20.95.kde@data)
#check to make sure it worked
plot(all_spdfs_together_cropped %>% filter(threshold==1), border="red")
plot(ca_coastline, add=TRUE)
########
#https://gis.stackexchange.com/questions/64537/clip-polygon-and-retain-data
row.names(pw.all.95.fkde1) <- row.names(burb.Feb13_20.95.kde)

spdf3 <- SpatialPolygonsDataFrame(pw.all.95.fkde1, data = pw.all.95.fkde1@polygons)

pw.all.95.fkde1@data<-pw.all.95.fkde@data

plot(spdf3, col=1:58)

###############
#http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS2_MergingSpatialData_part2_GeometricManipulations.html
# Convert from list to DataFrame and name result
result1 <- as.data.frame(pw.all.95.fkde1@polygons)
colnames(result) <- c("area")

# ids are stored as row-names, so make into a column.
ids <- as.character(rownames(result))

# Split the id into two columns based on the space in the middle.  Note the
# one place this could break is if the original names had spaces...
new.id.columns <- t(as.data.frame(strsplit(ids, " ")))
colnames(new.id.columns) <- c("id", "area")

# Now your results are a data.frame with separate columns for each id
# component!
result <- cbind(result, new.id.columns)
result
