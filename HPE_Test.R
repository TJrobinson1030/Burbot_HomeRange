################################################################################
################################################################################
######## HPE EVALUATION FOR BAD MEDICINE LAKE ARRAY ############################
######## CODE FROM MECKLEY ET AL 2014               ############################
################################################################################
################################################################################

rm(list=ls(all=TRUE))

require(PBSmapping)
require(ggplot2)
require(grid)
require(gridExtra)
require(sp)
require(readr)

#1 Fixed Tag Data 
# DATA FOR HERE IS FOUND IN EACH SECTION

#4 Fish VPS Data
Trial <- read_csv("D:/BSUGradBurbot/DATA/HPE_Data/TRANSMITTER-410mm_F_a_t-CALC-POSITIONS.csv", 
                  col_types = cols(DATETIME = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
str(Trial)

#5 Receiver File (Receiver locations during 3 20 day periods)
# CHANGE WHEN RUNNING ALL FOR DIFFERENT DATE RANGES
# RECEIVER LOCATION FOR APRIL
rec <- read_csv("D:/BSUGradBurbot/DATA/HPE_Data/Apr10_30RecLocation.csv")
str(rec)
# RECEIVER LOCATION FOR JUNE
rec <- read_csv("D:/BSUGradBurbot/DATA/HPE_Data/June20_Jul10RecLocation.csv")
str(rec)
#RECEIVER LOCATION FOR DECEMBER
rec <- read_csv("D:/BSUGradBurbot/DATA/HPE_Data/Dec21_Jan11RecLocation.csv")
str(rec)

#Pull in barrier
bm_outline <- read_csv("D:/BSUGradBurbot/DATA/HPE_Data/bm_outline.csv")
str(bm_outline)

################################################################################
################## REF TAG 1 APRIL DATA ########################################
################################################################################
# Fixed Tag Data
ref1 <- read_csv("D:/BSUGradBurbot/DATA/HPE_Data/TRANSMITTER-BMRef01-CALC-POSITIONS.csv")
str(ref1)

# RUN THE OBJECT YOU WANT TO EXTRACT DATA BY DATE WHEN YOU WANT TO CALCULATE HPE TO HPEm FOR EACH DATE RANGE
ref1 <- ref1[ref1$DATETIME >= "2019-04-10" & ref1$DATETIME <= "2019-04-30",]
ref1 <- ref1[ref1$DATETIME >= "2019-06-20" & ref1$DATETIME <= "2019-07-10",]
ref1 <- ref1[ref1$DATETIME >= "2019-12-21" & ref1$DATETIME <= "2020-01-11",]

#Known Easting and Northing (UTM) for fixed tag location if known...
 #...otherwise leave as NA
ref1East <-NA
ref1North <-NA
#tag transmission rate
transRateStat <- 600 #stationary tag average transmission rate
#transRateMobile <- 980 #tag drag tag test average transmission rate
#Select HPE filter, can be left at default until evaluating the array
HPEfilter <- 15
#Map Ploting Limits , a consistent cord-cartesian cutoff for each grid, look...
#...at study site and identify minimum and maximum easting and Northing
xAxisMin <-316153
xAxisMax <-320248  
yAxisMin <-5219142
yAxisMax <-5224226
#Error Range of interest, just leave this for most applications
ErrorRange <- c(5:30)

#Accuracy Goals: These are deterimined by your reserach questions and
#...analysis needs.
#Oultier Maximum
OutlierGoal <- 15

#Accuracy Goal for all Points
accGoal <- 3

#Goal for average accuracy if studying path lengths (1 order of magnitude...
#...less than max path length)
avgAccGoal<- 2.77

#Setup for Evaluations
HPELow <- 5 #Lowest HPE to consider in plot.
HPEHigh<- 25 #Higest HPE to consider in plot.

#Length of Tag Drag in Seconds
#tagDragLength <- 4860 #recorded during the test

#UTC time zone difference to correct for variable times
#timeZoneDiff=4*60*60

# Convert Latitude longitude to UTM
tagll <- data.frame(X = ref1$LON, Y = ref1$LAT)
attr(tagll, "projection") <- "LL"
xy <- convUL(tagll) * 1000

#Rename Postions
tag1Easting<- xy[,1]
tag1Northing<- xy[,2]
ref1$Easting<- tag1Easting
ref1$Northing<- tag1Northing

# Median Point if a measured point is unavailable
ref1East <- ifelse(is.na(ref1East), median(ref1$Easting),
                          + ref1East)
ref1North <- ifelse(is.na(ref1North), median(ref1$Northing),
                         + ref1North)
#calculate distance between each point and the average
ref1$error <- sqrt( (tag1Easting-(ref1East))^2 +
                           + (tag1Northing-(ref1North))^2)

#What is the fix rate based on transmit rate and test length
#convert date time
tagDateTime<- as.POSIXct((ref1$DATETIME), tz='UTC')
testLength<- as.numeric(max(tagDateTime))-as.numeric(min(tagDateTime))

#number of transmissions that should occur
possibleFixes <- (as.numeric(max(tagDateTime))-
                      + as.numeric(min(tagDateTime)))/transRateStat
#total number of transmisisons heard
actualFixes <- length(ref1$DATETIME)

#Fix Rate
fixRateStat <- (actualFixes/possibleFixes)*100

#Process tag drag data
#tagDragLat<-tagDrag$LAT
#tagDragLon<-tagDrag$LON
#tagDragDateTime<-tagDrag$DATETIME
#tagDragDateTime<-as.POSIXct((tagDragDateTime),tz='UTC')-timeZoneDiff
#Note that I subtract sections to standardize time (4hrs)
#tagDragHPE<-tagDrag$HPE
#tagDragdepth<-tagDrag$DEPTH
#tagDragDRX<-tagDrag$DRX
#convert latitude and longitude to UTM for tag drag
#tagll <- data.frame(X = tagDragLon, Y = tagDragLat)
#attr(tagll, "projection") <- "LL"
#xy <- convUL(tagll) * 1000
#Easting and Northing for tags
#tagDragEasting<- xy[,1]
#tagDragNorthing<- xy[,2]
#Format receiver locations
#preparing receiver data for use in plots
#tagll <- data.frame(X = rec2010$Longitude, Y=rec2010$Latitude)
#attr(tagll, "projection") <- "LL"
#xyrec <- convUL(tagll) * 1000
#rec2010$Easting<- xyrec[,1]
#rec2010$Northing<- xyrec[,2]
#Create a new tag dataframe to work from these variables
#tagDrag <-data.frame("Easting"= tagDragEasting, "Northing" =
#                           + tagDragNorthing,"dateTime" = tagDragDateTime, "HPE" =tagDragHPE,
#                         + "depth" = tagDragdepth , "DRX" = tagDragDRX)
#Pull out necessary columns from the GPS file
#gpsLat<- gpsDrag$Latitude
#gpsLon<- gpsDrag$Longitude
#gpsTime<-gpsDrag$Time
#gpsDate<-gpsDrag$Date
#gpsDateTime <- strptime(paste(gpsDate, gpsTime),
#                      + format="%m/%d/%Y %H:%M:%S", tz= 'UTC')
#Convert Gps positions from Lat./Lon. to UTM
#gpsll <- data.frame(X = gpsLon, Y = gpsLat)
#attr(gpsll, "projection") <- "LL"
#xy <- convUL(gpsll) * 1000
#gpsEasting<-xy[,1]
#gpsNorthing<-xy[,2]
#FIX RATE
#number of transmissions that should occur (drag times 13:12:59-13:57:24,
#...14:22:19-15:10:08, 15:35:08-15:59:42) which is 2665 s, 709 s, 1474 s,
#...corresponding to 162 expected transmissions.
#possibleFixesMobile <- tagDragLength/transRateMobile
#total number of transmisisons heard
#actualFixesMobile <- length(tagDrag$dateTime)
#Fix Rate
#fixRateMobile <- (actualFixesMobile/possibleFixesMobile)
#Estimating tag drag accuracy
# Estimate your true position from your GPS data
#tagDrag$gpsEasting <- approx(x=gpsDateTime,
 #                          + y=gpsEasting, xout=tagDrag$dateTime)$y
#tagDrag$gpsNorthing <- approx(x=gpsDateTime,
#                            + y=gpsNorthing, xout=tagDrag$dateTime)$y
#5 Calculate the Euclidian distance between esimated True and actual measured...
#...tag positions
#tagDrag$eucDist <- sqrt(((tagDrag$Easting-tagDrag$gpsEasting)^2)+
#                             + ((tagDrag$Northing-tagDrag$gpsNorthing)^2))
#gps <- data.frame(DateTime= gpsDateTime, Easting=gpsEasting,
#                   + Northing=gpsNorthing)
#estimateTimeOffset <- function(tOffset, tagDrag, gps){
  + # Estimate your true position from your GPS data
#    + tagDrag$gpsEasting <- approx(x=gps$DateTime, y=gps$Easting,
#                                   + xout=tagDrag$dateTime+tOffset)$y
#    + tagDrag$gpsNorthing <- approx(x=gps$DateTime,y=gps$Northing,
#                                    + xout=tagDrag$dateTime+tOffset)$y
#    +
#      + tagDrag$eucDist <- sqrt((tagDrag$Easting-tagDrag$Easting)^2 +
#                                  + ((tagDrag$Northing-tagDrag$gpsNorthing)^2))
#      +
#        + return(sum(tagDrag$eucDist^2,na.rm=T))
#      + }
#use optim to find tOffset that minimizes the sum of squared error
#offset.optim <- optim(0, fn=estimateTimeOffset, tag=tagDrag, gps=gps,
#                          + method="BFGS")
#tOffset <- offset.optim$par
#6 now calculate gps positions based on adusted tag position timestamps:
# this is because the Vemco times are not a perfect match.
# Estimate your true position from your GPS data
#tagDrag$gpsEasting2 <- approx(x=gpsDateTime, y=gps$Easting,
#                                  + xout=tagDrag$dateTime+tOffset)$y
#tagDrag$gpsNorthing2 <- approx(x=gpsDateTime,y=gps$Northing,
#                                 + xout=tagDrag$dateTime+tOffset)$y
#tagDrag$eucDist2 <- sqrt(((tagDrag$Easting-tagDrag$gpsEasting2)^2)+
#                             + ((tagDrag$Northing-tagDrag$gpsNorthing2)^2))
#Output a list of variables for accessing Array using stationary tags

#For Stationary Test
#Fix Rate %
fixRateStat

#Average Accuracy (m)
avgAccStat <- mean(ref1$error)
avgAccStat

#Median Accuracy (m)
medianAccStat <- median(ref1$error)
medianAccStat

#Proportion with goal accuracy
propAccStat<-length(which( ref1$error < accGoal ))/
+ length(ref1$error)
propAccStat

#Number of points with greater than outlier size
outlierCountStat <- length(which( ref1$error > OutlierGoal ))
outlierCountStat

#For tag drag Test
#Fix Rate
#fixRateMobile

#Average Accuracy
#meanAccMobile<- mean(na.omit(tagDrag$eucDist2))
#meanAccMobile

#Median Accuracy
#medianAccMobile <- median(na.omit(tagDrag$eucDist2))
#medianAccMobile

#Proportion with Goal accuracy
#propAccMobile <-length(which( tagDrag$eucDist2 < accGoal ))/
#  + length(tagDrag$eucDist2)
#propAccMobile

#Number of points with greater than outlier size (in our study 15 m)
#outlierCountMobile <- length(which( tagDrag$eucDist2 > OutlierGoal ))
#outlierCountMobile

#For Fixed Tag first
#What percentage of positions had and error below this value
percentGood1<- data.frame(ErrorRange)
percentGood1$goodperc<-NA

for( i in 1: length(percentGood1$ErrorRange)){percentGood1$goodperc [i] <- ((length(ref1$error[ref1$error < percentGood1$ErrorRange[i]])) / (length(ref1$error)))*100}

#tag drag Tag

#percentGood2<- data.frame(ErrorRange)
#percentGood2$goodperc<-NA
#for( i in 1: length(percentGood2$ErrorRange)){
#  + #What percentage of positions had error < value
#    + percentGood2$goodperc[i] <- ((length(tagDrag$eucDist2[tagDrag$eucDist2 <
#                                                              + percentGood2$ErrorRange[i]])) / (length(tagDrag$eucDist2)))*100
#    + }

#percent good
loc1 <- percentGood1
loc1$type<-c("Stationary Test 1")
#mobile1<- percentGood2
#mobile1$type<-c("Tag Drag Test")
AllperGood <- rbind(loc1)

plot0a <- ggplot(data=AllperGood, aes(x=ErrorRange,y=round(goodperc, digits=0)))+
  geom_point(size = 4)+
  geom_line()+
  theme_bw()

######## DO NOT RUN THIS SECTION ########
#plot0a <- ggplot(data=AllperGood, aes(x=ErrorRange,y=round(goodperc, digits=0))) +
#   geom_point(size= 4) +
#   geom_line() + 
#   theme(element_text(size=25),
#            axis.text.y = element_text(angle=90, hjust= 0.5),
#            axis.title.y = element_text(vjust= -0.02 ),
#            axis.title.x = element_text(vjust= -0.02 ),
#            plot.margin = unit(c(1,1.5,1,1.5), "cm"),
#            panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
#            legend.background = element_rect(fill = "transparent", colour = NA),
#            legend.position = c(0.6, 0.5),
#            legend.title=element_blank()) +
#  scale_y_continuous("Percent of all Positions",limits=c(0,100)) +
#  scale_x_continuous("Measured Error (m)",breaks= round(seq(min(AllperGood$ErrorRange), max(AllperGood$ErrorRange), by = 1),1)) +
#  coord_cartesian( xlim=c(0,29.5))
###############################################################  

#1. calculate the average error for each 1 unit HPE Bin.
#Bin by HPE
#create increments to bin by
breaks= seq(0, max(ref1$HPE),by=1)
# create a column with bins
ref1$bin <- cut(ref1$HPE,breaks)
#2. For each 1 m bin, the average HPE is calculated("binmean")
#use taplly to bin the error and calculate a min
binMean <- tapply(ref1$error,ref1$bin,mean)
#count the number of HPEs that fit into each bin value
binNum <- tapply(ref1$error,ref1$bin,length)
binNum[is.na(binNum)] <- 0
#New dataframe for bin data,cacluate xe and ye
bin<- data.frame(binMean, binNum)
#3. Every HPE value in a given 1 m bin has a corresponding HPEm. Each of...
#...these HPEm distances is composed of 2 elements: the error (difference...
#...between the calculated and measured position) in the X direction, and...
#...the error in the Y direction.

#I will refer to these error values as Xe and Ye, respectively.
#error in the x direction
ref1$xe <- sqrt((tag1Easting-(ref1East))^2)
#error in the y direction
ref1$ye <- sqrt((tag1Northing-(ref1North))^2)
#4. For each bin, the standard deviations of Xe and Ye are calculated.
bin$xeSd <- tapply(ref1$xe,ref1$bin,sd)
bin$yeSd <- tapply(ref1$ye,ref1$bin,sd)
#5. To convert the 2-dimensional standard deviations calculated in #4 into a...
#...single measure, the 2DRMS error is calculated from the standard...
#...deviations of Xe and Ye
bin$RMS2d <- 2*sqrt((bin$xeSd)^2 + (bin$yeSd)^2)
#6. Now create a line plot and a dataframe just for the...
#...numbers we need, (ie when we have at least 10 tag...
#...transmissions and an HPE less than 21)
bin$avgHPE <- tapply(ref1$HPE,ref1$bin,mean)
smallBin <- bin[ which(bin$binNum > 10),]
smallBin <- smallBin[which(smallBin$avgHPE < 25 ),]
res3 <- lm(smallBin$RMS2d ~ smallBin$avgHPE)

#PLOT THE POINT SPREAD
#Changing a bunch of theme elements to make it pretty.
  
xyplot<- data.frame(tag1Easting, tag1Northing)
plot1a<- qplot(tag1Easting, tag1Northing , alpha= 1/100,data = xyplot, xlab='Easting (m)', ylab= 'Northing (m) ') + theme_bw()
plot1a<- plot1a + geom_point(data=xyplot,aes(ref1East,ref1North),col="white", shape = 18, size= 4)
plot1a<- plot1a +theme(text = element_text(size=15),plot.title = element_text(vjust= 2), 
                       axis.text.y = element_text(angle=90, hjust= 0.5), 
                       axis.title.y = element_text(vjust= -0.02), 
                       axis.title.x = element_text(vjust= -0.02), 
                       plot.margin = unit(c(1,1.5,1,1.5), "cm"))
plot1a<- plot1a+theme(legend.position="none")#remove legend
plot1a<- plot1a + coord_cartesian(xlim=c(319525 ,319625),
                                  ylim=c(5223550,5223675))

#Graph of HPE to error relationship
plot2a<- qplot(ref1$HPE, ref1$error , alpha= 1/100,
               data = xyplot, xlab= 'HPE', ylab= 'Measured Error (m)')+
  scale_y_continuous(limits = c(0, 30)) +
  scale_x_continuous(limits = c(5, 30)) + theme_bw() +
  theme(plot.margin=unit(c(1,1.5,1,1.5), "cm"),panel.grid.major=element_blank(),
         panel.grid.minor=element_blank())+ theme(legend.position="none")
plot2a<- plot2a +theme(text = element_text(size=20),
                       plot.title = element_text(vjust= 2), axis.title.y = element_text(vjust= -0.02),
                       axis.title.x = element_text(vjust= -0.02))
plot2a <- plot2a + geom_point(data=smallBin, shape= 19,
                              size= 5,col='black', bg= 'black',alpha= 1, aes(x = avgHPE, y = RMS2d)) +
                              theme(legend.position = "none")
plot2a<- plot2a + geom_point(data=smallBin, shape= 19, size= 4.2,col='white',
                             alpha= 1, aes(x = avgHPE, y = RMS2d)) + theme(legend.position = "none")
plot2a<- plot2a + geom_point(data=smallBin, shape= 4, size = 3,color='red',
                             aes(x = avgHPE, y = RMS2d)) +
  theme(legend.position = "none")

plot2a<- plot2a +geom_abline(data=res3, col='white',
                             aes(intercept=res3$coefficients[1], slope=res3$coefficients[2] ), size=1)
plot2a<- plot2a +geom_abline(data=res3, col='black', linetype='dashed',
                             aes(intercept=res3$coefficients[1], slope=res3$coefficients[2]), size= 1)
#Add model equation and r2 to graph
#function to add linear model to graph
#Setup for a linear model to be added to a ggplot
lm_eqn = function(m) {
     l <- list(a = format(coef(m)[1], digits = 2),
                 b = format(abs(coef(m)[2]), digits = 2),
                 r2 = format(summary(m)$r.squared, digits = 3));
     if (coef(m)[2] >= 0) {
       eq <- substitute(italic(y) == a + b %.% italic(x)*","
                          ~~italic(r)^2~"="~r2,l)
       } else {
         eq <- substitute(italic(y) == a - b %.% italic(x)*","
                            ~~italic(r)^2~"="~r2,l)
         }
     as.character(as.expression(eq));
     }
#calls in linear model function for linear model res3
plot2a=plot2a + annotate("text",x = 15, y = 20,
                         label = lm_eqn(res3), size= 4, parse=TRUE)

#What type of HPE cutoff should we use if we know what level of error we want
#!#! YOU MUST CHANGE accGoal to your desired error
#now evaluate how HPE cutoffs work if we want a specific error value
HPE<- c(5:HPEHigh)# we want a list of all possible HPE cutoffs
AllHPE <- data.frame(HPE)
AllHPE$ptsReject <- NA# number of points removed
AllHPE$ptsRejectP <- NA# percetage of all ponits
AllHPE$ptsRetain <- NA# number of poitns retained
AllHPE$ptsRetainP <- NA# percetage of all ponits
AllHPE$incorrectReject <- NA# good= error < accGoal and unacceptably...
#...erroneous= error > accGoal
AllHPE$incorrectRejectP <- NA
AllHPE$incorrectRetain<- NA
AllHPE$incorrectRetainP <- NA
AllHPE$incorrectRetainPvsRetain <-NA # incorretly retained/all retained *100
AllHPE$correctReject <- NA
AllHPE$correctRejectP <- NA

AllHPE$correctRetain <- NA
AllHPE$correctRetainP <- NA
AllHPE$goodDataLossP <- NA# incorrectReject/goodpts *100
AllHPE$badDataRetainP <- NA# incorrectRetain/AllHPE$ptsRetain * 100
AllHPE$sdDEV<-NA
AllHPE$avgErr<-NA
AllHPE$maxErr <- NA

# Number of all points that were unacceptably erroneous
badpts <-length(ref1$error[ref1$error > accGoal ])
# % of all points were unacceptably erroneous
badptsP <-(length(ref1$error[ref1$error > accGoal ])/(length(
    + ref1$error)))*100
#Number of all points that are good
goodpts <-length(ref1$error[ref1$error <= accGoal ])
#% of all points that are good
goodptsP <- (length(ref1$error[ref1$error <= accGoal ])/(length(ref1$error)))*100

for( i in 1: length(AllHPE$HPE)){AllHPE$ptsReject[i] <- length( which(ref1$HPE > AllHPE$HPE[i]))
    AllHPE$ptsRejectP[i] <- length( which(ref1$HPE > AllHPE$HPE[i]))/(length(ref1$error))*100 
    #% of all points dropped
    AllHPE$ptsRetain[i] <- length( which(ref1$HPE < AllHPE$HPE[i]))
    AllHPE$ptsRetainP[i] <- length( which(ref1$HPE < AllHPE$HPE[i]))/(length(ref1$error))*100 
    #% of all points dropped
    #incorrect reject
    AllHPE$incorrectReject[i] <- length( which( (ref1$HPE[ ref1$error <= accGoal ] > AllHPE$HPE[i]) == TRUE))
        #divided by total with acceptable error
        AllHPE$incorrectRejectP[i] <-(AllHPE$incorrectReject[i]/length(ref1$error[ref1$error <= accGoal ]) )*100
        #incorrect retain
          AllHPE$incorrectRetain[i] <-length( which((ref1$HPE[ ref1$error > accGoal ] < AllHPE$HPE[i]) == TRUE))
            #divided by accurate points
            AllHPE$incorrectRetainP[i] <-(AllHPE$incorrectRetain[i]/length(ref1$error[ref1$error > accGoal ]))*100
            #divided by total retained by filter
            AllHPE$incorrectRetainPvsRetain[i] <-(AllHPE$incorrectRetain[i]/AllHPE$ptsRetain[i] )*100
            #Correctly Rejected
            AllHPE$correctReject[i] <-(length(which((ref1$HPE[ref1$error > accGoal] > AllHPE$HPE[i]) == TRUE)))
            #divided by sum of all positions
            AllHPE$correctRejectP[i]<- (AllHPE$correctReject[i]/AllHPE$ptsRetain[i])*100
            #Correctly Retained
            AllHPE$correctRetain[i] <- length(which((ref1$HPE[ ref1$error < accGoal] < AllHPE$HPE[i]) == TRUE))
            #divided by sum of all positions
            AllHPE$correctRetainP[i]<-(AllHPE$correctRetain[i]/length(ref1$error))*100
            AllHPE$goodDataLossP[i] <- (AllHPE$incorrectReject[i] / goodpts)*100 
            #incorrectReject/goodpts *100
            AllHPE$badDataRetainP[i] <- (AllHPE$incorrectRetain[i]/AllHPE$ptsRetain[i])*100 
            # incorrectRetain/AllHPE$ptsRetain * 100
            AllHPE$avgErr[i] <- mean(ref1$error[ ref1$HPE < AllHPE$HPE[i]]) 
            #mean of all error for each hpe
            AllHPE$sdDEV [i]<- sd(ref1$error[ ref1$HPE < AllHPE$HPE[i]]) 
            #Standard deviation of all error for each HPE
            AllHPE$maxErr[i] <- max(ref1$error[ ref1$HPE < AllHPE$HPE[i]]) 
            #mean of all error for each hpe
            AllHPE$count15[i] <- length(which( (ref1$HPE[ ref1$error > 15] < AllHPE$HPE[i]) == TRUE ))
            }
start <- HPELow
end <- HPEHigh
AllHPEperf1 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$correctRejectP[start:end])
AllHPEperf1$var <-"Correct Rejection"
AllHPEperf2 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$incorrectRetainP[start:end])
AllHPEperf2$var <-"Incorrect Retainment"
AllHPEperf3 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$incorrectRejectP[start:end])
AllHPEperf3$var <-"Incorrect Rejection"
#Combine all seconts for graph
AllHPEperf <- rbind(AllHPEperf1,AllHPEperf3,AllHPEperf2)

# PLOT FOR INCORRECTLY REJECTED AND RETAINED POINTS
#note legend is alphabetical but graph is by dataframe order
plot3a <- qplot( AllHPE$incorrectRetainPvsRetain, AllHPE$incorrectRejectP,size=3, xlab= '% Incorrecly Retained (of all retained positions)',ylab= ' % Incorrectly Rejected (of all acceptable positions)' ) + 
  theme_bw() + theme( plot.margin = unit(c(1,1.5,1,1.5), "cm"),panel.grid.major=element_blank(), panel.grid.minor=element_blank())

plot3a<- plot3a +geom_line(data=AllHPE, aes(x=AllHPE$incorrectRetainPvsRetain, y=AllHPE$incorrectRejectP),color= 'black', size=0.25)
plot3a<- plot3a + theme_bw() + theme( text = element_text(size=15), #plot.title = element_text(vjust= 2),
                                        axis.text.y = element_text(angle=90, hjust= 0.4),
                                        axis.title.y = element_text(vjust= .3),
                                        axis.title.x = element_text(vjust= -0.02 ),
                                        plot.margin = unit(c(1,1.5,1,1.5), "cm"),
                                        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                        legend.position = c(1.5, 0.8),
                                        legend.background = element_rect(fill = "transparent",colour = NA),
                                        legend.title=element_blank() )
plot3a <-plot3a+ geom_text(cex=6,aes(x=(AllHPE$incorrectRetainPvsRetain[3:15]),
                                       y=(AllHPE$incorrectRejectP[3:15])+4 , label=AllHPE$HPE[3:15], group=NULL),)                      

# PLOT FOR AVERAGE ERROR IN METERS AND THE ASSOCIATED HPE VALUE
plot4a <- ggplot( data= AllHPE[start:end,],
                   aes(x = AllHPE$HPE[start:end])) +
   geom_point(aes(y = AllHPE$avgErr[start:end],
                    shape= "Mean Error" ), size= 4, color= 'black') +
   geom_point(aes(y = AllHPE$maxErr[start:end],
                    shape="Maximum Error" ), size= 4, color='black') +
   geom_line(aes(y = OutlierGoal, ), size= 1, lty= 'dashed') +
   theme_bw() + theme( text = element_text(size=25),
                         plot.title = element_text(vjust= 2),
                         axis.text.y = element_text(angle=90, hjust= 0.5),
                         axis.title.y = element_text(vjust= -0.02 ),
                         axis.title.x = element_text(vjust= -0.02 ),
                         plot.margin = unit(c(1,1.5,1,1.5), "cm"),
                         panel.grid.major=element_blank(),
                         panel.grid.minor=element_blank(),
                         legend.position = c(0.15, 0.95),
                         legend.background = element_rect(
                           fill = "transparent",colour = NA),
                         legend.title=element_blank() ) +
   xlab('HPE' ) +
 coord_cartesian( ylim=c(0,2.5)) +
   scale_y_continuous("Error (m)",breaks= round(seq(min(0),
                                                      max(2.5), by = 0.5),1)) +
   scale_x_continuous("HPE",breaks=round(seq(min(AllHPE$HPE),
                                               max(AllHPE$HPE), by = 1),1)) +
   geom_text(cex=6,aes(x=AllHPE$HPE[1:21],y= (29.25),
                         label=AllHPE$count15[1:21], group=NULL),)
plot4a


########################################
####### ALL CODE FOR MOVING TAG DATA########
###### WE DO NOT HAVE THIS ##############

#with(tagDrag, eucDist-eucDist2)
#badHPE<- tagDrag$HPE[tagDrag$HPE >= HPEfilter ]
#badeucDist3 <-tagDrag$eucDist2[tagDrag$HPE >= HPEfilter]
#badHPENorthing<- tagDrag$Northing[tagDrag$HPE >= HPEfilter ]
#badHPEEasting<- tagDrag$Easting[tagDrag$HPE >= HPEfilter ]
#badErrorNorthing<- tagDrag$Northing[tagDrag$eucDist2 >= accGoal ]
#badErrorEasting<- tagDrag$Easting[tagDrag$eucDist2 >= accGoal ]
#eucDist3<- tagDrag$eucDist2[tagDrag$eucDist2 >= accGoal ]
#incorrectlyRejected<- (tagDrag[tagDrag$eucDist2 <6 & tagDrag$HPE>8,])
#correctlyRejected<- na.omit(tagDrag[tagDrag$eucDist2 >6 & tagDrag$HPE>8,],)
#incorrectlyRetained<- na.omit(tagDrag[tagDrag$eucDist2 >6 & tagDrag$HPE<8,],)

#correctlyRetained <- na.omit(tagDrag[tagDrag$eucDist2 < 6 & tagDrag$HPE<8,],)
#badError<- data.frame("Easting"=badErrorEasting,"Northing"=badErrorNorthing,
#                      + "eucDist3" = eucDist3, "HPE" = tagDrag$HPE[tagDrag$eucDist2 >= accGoal])
#badError <- na.omit(badError)
#New frame here
#tag<-tagDrag[ is.na(tagDrag$gpsEasting)==FALSE,]

# identify unacceptably erroneous points in the plot(must change HPE)
#plot5a<- qplot(correctlyRetained$gpsEasting2, correctlyRetained$gpsNorthing2,
#                 + xlab= ' Easting (m)', ylab= 'Northing (m) ') + theme_bw() + geom_point(data=
#                                                                                              + coast,aes(coast$Easting,coast$Northing),col="black", shape = 18, size= 5)
#plot5a<- plot5a +theme(text = element_text(size=20),
#                         + plot.title = element_text(vjust= 2),
#                         + axis.text.y = element_text(angle=90, hjust= 0.5),
#                         + axis.title.y = element_text(vjust=0.4 ),
#                         + axis.title.x = element_text(vjust= -0.01 ),
#                         + legend.position="none" ) #make the axis
#receivers
#plot5a<- plot5a + geom_point(data=rec2010, aes(x=Easting, y=Northing), shape= 3)
#all locations with a high HPE
#plot5a <- plot5a + geom_point(data= incorrectlyRejected,aes(x=Easting,
#                                                            + y=Northing), shape = 19, size= 5, col='Orange') # circled orange
#plot5a <- plot5a + geom_point(data= correctlyRejected,aes(x=Easting,
#                                                          + y=Northing), shape = 19, size= 5, col='Blue') # circled blue
#plot5a <- plot5a + geom_point(data= incorrectlyRetained,aes(x=Easting,
#                                                            + y=Northing), shape = 19, size= 5, col='red') # circled red
#plot5a <- plot5a + geom_point(data=correctlyRetained,aes(x=Easting,y=Northing),
#                              + shape = 19, size= 5, col='black') # circled black
#plot5a<- plot5a +coord_cartesian(xlim=c(xAxisMin,xAxisMax ),
#                                 + ylim=c(yAxisMin,yAxisMax))
#plot5b<- qplot(correctlyRetained$HPE, correctlyRetained$eucDist2 , alpha= 1/100,
#               + xlab= 'HPE', ylab= 'Euclidian distance (m)' )+ scale_y_continuous(limits =
#                                                                                       + c(0,21)) + scale_x_continuous(limits = c(0, 21)) + theme_bw() + theme(
#                                                                                         + plot.margin = unit(c(1,1.5,1,1.5), "cm"), panel.grid.major=element_blank(),
#                                                                                         + panel.grid.minor=element_blank())+ theme(legend.position="none") +
#  + scale_colour_identity()
#plot5b<- plot5b +theme(text = element_text(size=20),
#                       + plot.title = element_text(vjust= 2),
#                       + axis.title.y = element_text(vjust= -0.02 ),
#                       + axis.title.x = element_text(vjust= -0.02 ) ) #make the axis
#plot5b<- plot5b + geom_point(data= incorrectlyRejected,aes(x=HPE,y=eucDist2,
#                                                           + col="orange"), alpha=1,shape = 19, size= 4)
#plot5b<- plot5b + geom_point(data= correctlyRejected,aes(x=HPE, y=eucDist2,
#                                                         + col="blue"), shape = 19, size= 4)
#plot5b<- plot5b + geom_point(data= incorrectlyRetained,aes(x=HPE, y= eucDist2,
#                                                           + col="red"), shape = 19, size= 4)
########################################################################################################


#Area of Receiver Coverage
tagll <- data.frame(X = rec$LONG, Y = rec$LAT)
attr(tagll, "projection") <- "LL"
xyrec <- convUL(tagll) * 1000
rec$EASTING<- xyrec[,1]
rec$NORTHING<- xyrec[,2]

#Converting fish location to UTM
tagll <- data.frame(X = Trial$LON, Y = Trial$LAT)
attr(tagll, "projection") <- "LL"
xyrec <- convUL(tagll) * 1000
Trial$EASTING<- xyrec[,1]
Trial$NORTHING<- xyrec[,2]

#Evaluate Actual fish trackin data (HPE FILTER)
plot6a <- qplot(Trial$EASTING, Trial$NORTHING, xlab = ' Easting (m)', ylab = 'Northing (m)', size=I(1)) + theme_bw()
plot6a <- plot6a + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                   ylim=c(yAxisMin,yAxisMax))
plot6a <- plot6a +theme(text = element_text(size=15),
                        plot.title = element_text(vjust= 2),
                        axis.text.y = element_text(angle=90, hjust= 0.5),
                        axis.title.y = element_text(vjust= -0.02 ),
                        axis.title.x = element_text(vjust= -0.02 ),
                        plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot6a <- plot6a + geom_point(data=Trial[Trial$HPE < 3,], aes(x=EASTING,
                                                         y=NORTHING),col="blue", shape = 20, size= I(1))
plot6a <- plot6a + geom_point(data=Trial[Trial$HPE >= 3 &
                                          Trial$HPE < 5 ,],aes(x=EASTING, y=NORTHING),col="#56B4E9",
                                shape = 20, size= I(1))
plot6a <- plot6a + geom_point(data=Trial[Trial$HPE >= 5 &
                                          Trial$HPE < 10 ,],aes(x=EASTING, y=NORTHING),col="green",
                                shape = 20, size= I(1))
plot6a <- plot6a + geom_point(data=Trial[Trial$HPE >= 10 &
                                          Trial$HPE < 15 ,],aes(x=EASTING, y=NORTHING),col="yellow",
                                shape = 20,size= I(1))
plot6a <- plot6a + geom_point(data=Trial[Trial$HPE >= 15 &
                                          Trial$HPE < 20 ,],aes(x=EASTING, y=NORTHING),col="#FFCC33",
                                shape = 20,size= I(1))
plot6a <- plot6a + geom_point(data=Trial[Trial$HPE >= 20 ,],aes(x=EASTING,
                                                            y=NORTHING),col="red", shape = 20, size= I(1))
plot6a <- plot6a + geom_point(data=rec,aes(EASTING,NORTHING),fill="black",
                              col="white", shape = 24, size= 3, alpha=1)
plot6a <- plot6a + geom_point(data=bm_outline,aes(e,n),col="black",
                              shape = 18, size= 1)

xydata<-data.frame(X=Trial$EASTING, Y=Trial$NORTHING)

#create a polygon that covers the area of the data
fishPoly<- calcConvexHull(xydata, keepExtra=FALSE)

#Polygon(fishPoly)
#calcArea, calcCentroid are likley both useful

#Create a column in pos data.frame whethere in or outside array
Trial$array<- point.in.polygon(Trial$EASTING,Trial$NORTHING, rec$EASTING,
                             + rec$NORTHING, mode.checked= TRUE)
#plot of array polygon, with fish positions remvoed, 0 means not in array

#2 Filter data FOR FIRST SET OF PLOTS
cutoff <- HPEfilter
posDrop <- nrow(Trial[(Trial$HPE > cutoff),])
posDropArrayT <- ((Trial$HPE[Trial$array != 0]) > cutoff)
posDropArray <- length(subset(posDropArrayT, posDropArrayT == TRUE))
DropOutT <- ((Trial$HPE[Trial$array == 0]) > cutoff)
posDropOut <- length(subset( DropOutT, DropOutT == TRUE))
posTotal1<- nrow(Trial)
posTotal2<- length(Trial$array[Trial$array == 0])
posTotal3<- length(Trial$array[Trial$array != 0])
percDropFil1 <- (posDrop/posTotal1)*100 # % DroppedTota;
percADropFil1 <- (posDropArray/posTotal2)*100 #% dropped in the array
percODropFil1 <- (posDropOut/posTotal3)*100 #% dropped outside of the array
#Now identify the silhoutte of the point cluster using calcConvexHull
xydata1<-data.frame(X=Trial$EASTING[Trial$HPE < cutoff ], Y=Trial$NORTHING[Trial$HPE <
                                                                     + cutoff ])
#create a polygon that covers the area of the data
fishPoly1<- calcConvexHull(xydata1, keepExtra=FALSE)
head(xydata1)
plot(Trial$EASTING,Trial$NORTHING)
points(fishPoly1$X, fishPoly1$Y, col='red')
# #calcArea, calcCentroid are likley both useful
areaFil1 <- calcArea(fishPoly1) #calculates area of polygon

#Plot 1- all fish positions
plot7a <- qplot(Trial$EASTING, Trial$NORTHING, xlab= ' Easting (m)', ylab=
                 'Northing (m)',alpha=1/100, size=I(1)) + theme_bw()
plot7a <- plot7a + coord_cartesian(xlim=c(xAxisMin,xAxisMax), ylim=
                                    c(yAxisMin,yAxisMax))
plot7a <- plot7a + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                               shape = 17, size= 2, alpha=1)
plot7a <- plot7a +theme(text = element_text(size=13),
                         plot.title = element_text(vjust= 2),
                         axis.text.y = element_text(angle=90, hjust= 0.5),
                         axis.title.y = element_text(vjust= -0.02 ),
                         axis.title.x = element_text(vjust= -0.02 ),
                         plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7a <- plot7a + theme(legend.position="none") #remove legend
#plot 2- all positions with HPE < cutoff
plot7b <- qplot( Trial$EASTING[Trial$HPE<= cutoff], Trial$NORTHING[Trial$HPE<= cutoff
                                                              ], xlab= ' Easting (m)', ylab= 'Northing (m) ', alpha=1/100, size=I(1)) +
   theme_bw()
plot7b<- plot7b + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                  ylim=c(yAxisMin,yAxisMax))
plot7b<- plot7b + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                             shape = 17, size= 2, alpha=1)
plot7b<- plot7b +theme(text = element_text(size=13),
                        plot.title = element_text(vjust= 2),
                        axis.text.y = element_text(angle=90, hjust= 0.5),
                        axis.title.y = element_text(vjust= -0.02 ),
                        axis.title.x = element_text(vjust= -0.02 ),
                        plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7b<- plot7b + theme(legend.position="none") #remove legend
#plot 3- all positions with HPE > cutoff)
plot7c <- qplot( Trial$EASTING[Trial$HPE> cutoff], Trial$NORTHING[Trial$HPE> cutoff],
                  xlab= ' Easting (m)', ylab= 'Northing (m) ',alpha=1/100, size=I(1)) +
   theme_bw()
plot7c<- plot7c + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                   ylim=c(yAxisMin,yAxisMax))
plot7c<- plot7c + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                              shape = 17, size= 2, alpha=1)
plot7c<- plot7c +theme(text = element_text(size=13),
                        plot.title = element_text(vjust= 2),
                        axis.text.y = element_text(angle=90, hjust= 0.5),
                        axis.title.y = element_text(vjust= -0.02 ),
                        axis.title.x = element_text(vjust= -0.02 ),
                        plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7c<- plot7c + theme(legend.position="none") #remove legend


################################################################################
################################################################################
###########  REF TAG 2 DATA  #############################################
################################################################################
################################################################################

# RUN THE OBJECT YOU WANT TO EXTRACT DATA BY DATE WHEN YOU WANT TO CALCULATE HPE TO HPEm FOR EACH DATE RANGE
ref2 <- read_csv("D:/BSUGradBurbot/DATA/HPE_Data/TRANSMITTER-BMRef02-CALC-POSITIONS.csv")
str(ref2)
ref2 <- ref2[ref2$DATETIME >= "2019-04-10" & ref2$DATETIME <= "2019-04-30",]
ref2 <- ref2[ref2$DATETIME >= "2019-06-20" & ref2$DATETIME <= "2019-07-10",]
ref2 <- ref2[ref2$DATETIME >= "2019-12-21" & ref2$DATETIME <= "2020-01-11",]

#Known Easting and Northing (UTM) for fixed tag location if known...
#...otherwise leave as NA
ref2East <-NA
ref2North <-NA

#tag transmission rate
transRateStat <- 600 #stationary tag average transmission rate

#transRateMobile <- 980 #tag drag tag test average transmission rate
#Select HPE filter, can be left at default until evaluating the array
HPEfilter <- 17
#Map Ploting Limits , a consistent cord-cartesian cutoff for each grid, look...
#...at study site and identify minimum and maximum easting and Northing
xAxisMin <-316153
xAxisMax <-320248  
yAxisMin <-5219142
yAxisMax <-5224226
#Error Range of interest, just leave this for most applications
ErrorRange <- c(5:30)

#Accuracy Goals: These are deterimined by your reserach questions and
#...analysis needs.
#Oultier Maximum
OutlierGoal <- 12

#Accuracy Goal for all Points
accGoal <- 6

#Goal for average accuracy if studying path lengths (1 order of magnitude...
#...less than max path length)
avgAccGoal<- 2.77

#Setup for Evaluations
HPELow <- 5 #Lowest HPE to consider in plot.
HPEHigh<- 30 #Higest HPE to consider in plot.

#Length of Tag Drag in Seconds
#tagDragLength <- 4860 #recorded during the test

#UTC time zone difference to correct for variable times
#timeZoneDiff=4*60*60

# Convert Latitude longitude to UTM
tagll <- data.frame(X = ref2$LON, Y = ref2$LAT)
attr(tagll, "projection") <- "LL"
xy <- convUL(tagll) * 1000

#Rename Postions
tag2Easting<- xy[,1]
tag2Northing<- xy[,2]
ref2$Easting<- tag2Easting
ref2$Northing<- tag2Northing

# Median Point if a measured point is unavailable
ref2East <- ifelse(is.na(ref2East), median(ref2$Easting),
                   + ref2East)
ref2North <- ifelse(is.na(ref2North), median(ref2$Northing),
                    + ref2North)
#calculate distance between each point and the average
ref2$error <- sqrt( (tag2Easting-(ref2East))^2 +
                         + (tag2Northing-(ref2North))^2)

#What is the fix rate based on transmit rate and test length
#convert date time
tagDateTime<- as.POSIXct((ref2$DATETIME), tz='UTC')
testLength<- as.numeric(max(tagDateTime))-as.numeric(min(tagDateTime))

#number of transmissions that should occur
possibleFixes <- (as.numeric(max(tagDateTime))-
                    + as.numeric(min(tagDateTime)))/transRateStat
#total number of transmisisons heard
actualFixes <- length(ref2$DATETIME)

#Fix Rate
fixRateStat <- (actualFixes/possibleFixes)*100

#####################################
#### DO NOT RUN THIS SECTION ########
#####################################
#Process tag drag data
#tagDragLat<-tagDrag$LAT
#tagDragLon<-tagDrag$LON
#tagDragDateTime<-tagDrag$DATETIME
#tagDragDateTime<-as.POSIXct((tagDragDateTime),tz='UTC')-timeZoneDiff
#Note that I subtract sections to standardize time (4hrs)
#tagDragHPE<-tagDrag$HPE
#tagDragdepth<-tagDrag$DEPTH
#tagDragDRX<-tagDrag$DRX
#convert latitude and longitude to UTM for tag drag
#tagll <- data.frame(X = tagDragLon, Y = tagDragLat)
#attr(tagll, "projection") <- "LL"
#xy <- convUL(tagll) * 1000
#Easting and Northing for tags
#tagDragEasting<- xy[,1]
#tagDragNorthing<- xy[,2]
#Format receiver locations
#preparing receiver data for use in plots
#tagll <- data.frame(X = rec2010$Longitude, Y=rec2010$Latitude)
#attr(tagll, "projection") <- "LL"
#xyrec <- convUL(tagll) * 1000
#rec2010$Easting<- xyrec[,1]
#rec2010$Northing<- xyrec[,2]
#Create a new tag dataframe to work from these variables
#tagDrag <-data.frame("Easting"= tagDragEasting, "Northing" =
#                           + tagDragNorthing,"dateTime" = tagDragDateTime, "HPE" =tagDragHPE,
#                         + "depth" = tagDragdepth , "DRX" = tagDragDRX)
#Pull out necessary columns from the GPS file
#gpsLat<- gpsDrag$Latitude
#gpsLon<- gpsDrag$Longitude
#gpsTime<-gpsDrag$Time
#gpsDate<-gpsDrag$Date
#gpsDateTime <- strptime(paste(gpsDate, gpsTime),
#                      + format="%m/%d/%Y %H:%M:%S", tz= 'UTC')
#Convert Gps positions from Lat./Lon. to UTM
#gpsll <- data.frame(X = gpsLon, Y = gpsLat)
#attr(gpsll, "projection") <- "LL"
#xy <- convUL(gpsll) * 1000
#gpsEasting<-xy[,1]
#gpsNorthing<-xy[,2]
#FIX RATE
#number of transmissions that should occur (drag times 13:12:59-13:57:24,
#...14:22:19-15:10:08, 15:35:08-15:59:42) which is 2665 s, 709 s, 1474 s,
#...corresponding to 162 expected transmissions.
#possibleFixesMobile <- tagDragLength/transRateMobile
#total number of transmisisons heard
#actualFixesMobile <- length(tagDrag$dateTime)
#Fix Rate
#fixRateMobile <- (actualFixesMobile/possibleFixesMobile)
#Estimating tag drag accuracy
# Estimate your true position from your GPS data
#tagDrag$gpsEasting <- approx(x=gpsDateTime,
#                          + y=gpsEasting, xout=tagDrag$dateTime)$y
#tagDrag$gpsNorthing <- approx(x=gpsDateTime,
#                            + y=gpsNorthing, xout=tagDrag$dateTime)$y
#5 Calculate the Euclidian distance between esimated True and actual measured...
#...tag positions
#tagDrag$eucDist <- sqrt(((tagDrag$Easting-tagDrag$gpsEasting)^2)+
#                             + ((tagDrag$Northing-tagDrag$gpsNorthing)^2))
#gps <- data.frame(DateTime= gpsDateTime, Easting=gpsEasting,
#                   + Northing=gpsNorthing)
#estimateTimeOffset <- function(tOffset, tagDrag, gps){
+ # Estimate your true position from your GPS data
  #    + tagDrag$gpsEasting <- approx(x=gps$DateTime, y=gps$Easting,
  #                                   + xout=tagDrag$dateTime+tOffset)$y
  #    + tagDrag$gpsNorthing <- approx(x=gps$DateTime,y=gps$Northing,
  #                                    + xout=tagDrag$dateTime+tOffset)$y
  #    +
  #      + tagDrag$eucDist <- sqrt((tagDrag$Easting-tagDrag$Easting)^2 +
  #                                  + ((tagDrag$Northing-tagDrag$gpsNorthing)^2))
  #      +
  #        + return(sum(tagDrag$eucDist^2,na.rm=T))
  #      + }
  #use optim to find tOffset that minimizes the sum of squared error
#offset.optim <- optim(0, fn=estimateTimeOffset, tag=tagDrag, gps=gps,
#                          + method="BFGS")
#tOffset <- offset.optim$par
#6 now calculate gps positions based on adusted tag position timestamps:
# this is because the Vemco times are not a perfect match.
# Estimate your true position from your GPS data
#tagDrag$gpsEasting2 <- approx(x=gpsDateTime, y=gps$Easting,
#                                  + xout=tagDrag$dateTime+tOffset)$y
#tagDrag$gpsNorthing2 <- approx(x=gpsDateTime,y=gps$Northing,
#                                 + xout=tagDrag$dateTime+tOffset)$y
#tagDrag$eucDist2 <- sqrt(((tagDrag$Easting-tagDrag$gpsEasting2)^2)+
#                             + ((tagDrag$Northing-tagDrag$gpsNorthing2)^2))
#Output a list of variables for accessing Array using stationary tags
################################################################################


#For Stationary Test
#Fix Rate %
fixRateStat

#Average Accuracy (m)
avgAccStat <- mean(ref2$error)
avgAccStat

#Median Accuracy (m)
medianAccStat <- median(ref2$error)
medianAccStat

#Proportion with goal accuracy
propAccStat<-length(which( ref2$error < accGoal ))/
  + length(ref2$error)
propAccStat

#Number of points with greater than outlier size
outlierCountStat <- length(which( ref2$error > OutlierGoal ))
outlierCountStat

#For tag drag Test
#Fix Rate
#fixRateMobile

#Average Accuracy
#meanAccMobile<- mean(na.omit(tagDrag$eucDist2))
#meanAccMobile

#Median Accuracy
#medianAccMobile <- median(na.omit(tagDrag$eucDist2))
#medianAccMobile

#Proportion with Goal accuracy
#propAccMobile <-length(which( tagDrag$eucDist2 < accGoal ))/
#  + length(tagDrag$eucDist2)
#propAccMobile

#Number of points with greater than outlier size (in our study 15 m)
#outlierCountMobile <- length(which( tagDrag$eucDist2 > OutlierGoal ))
#outlierCountMobile

#For Fixed Tag first
#What percentage of positions had and error below this value
percentGood1<- data.frame(ErrorRange)
percentGood1$goodperc<-NA

for( i in 1: length(percentGood1$ErrorRange)){percentGood1$goodperc [i] <- ((length(ref2$error[ref2$error < percentGood1$ErrorRange[i]])) / (length(ref2$error)))*100}

#tag drag Tag
#percentGood2<- data.frame(ErrorRange)
#percentGood2$goodperc<-NA
#for( i in 1: length(percentGood2$ErrorRange)){
#  + #What percentage of positions had error < value
#    + percentGood2$goodperc[i] <- ((length(tagDrag$eucDist2[tagDrag$eucDist2 <
#                                                              + percentGood2$ErrorRange[i]])) / (length(tagDrag$eucDist2)))*100
#    + }


#percent good
loc1 <- percentGood1
loc1$type<-c("Stationary Test 1")
#mobile1<- percentGood2
#mobile1$type<-c("Tag Drag Test")
AllperGood <- rbind(loc1)

plot0a2 <- ggplot(data=AllperGood, aes(x=ErrorRange,y=round(goodperc, digits=0)))+
  geom_point(size = 4)+
  geom_line()+
  theme_bw()

#plot0a2 <- ggplot(data=AllperGood, aes(x=ErrorRange,y=round(goodperc, digits=0))) +
 # geom_point(size= 4)+geom_line()+theme_bw() +
  #theme(element_text(size=25),
   #     axis.text.y = element_text(angle=90, hjust= 0.5),
    #    axis.title.y = element_text(vjust= -0.02 ),
     #   axis.title.x = element_text(vjust= -0.02 ),
      #  plot.margin = unit(c(1,1.5,1,1.5), "cm"),
       # panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        #legend.background = element_rect(fill = "transparent", colour = NA),
        #legend.position = c(0.6, 0.5),
        #legend.title=element_blank()) +
#  scale_y_continuous("Percent of all Positions",limits=c(0,100)) +
 # scale_x_continuous("Measured Error (m)",breaks= round(seq(min(
  #  AllperGood$ErrorRange), max(AllperGood$ErrorRange), by = 1),1)) +
  #coord_cartesian( xlim=c(0,20.5))

#1. calculate the average error for each 1 unit HPE Bin.
#Bin by HPE
#create increments to bin by
breaks= seq(0, max(ref2$HPE),by=1)
# create a column with bins
ref2$bin <- cut(ref2$HPE,breaks)
#2. For each 1 m bin, the average HPE is calculated("binmean")
#use taplly to bin the error and calculate a min
binMean <- tapply(ref2$error,ref2$bin,mean)
#count the number of HPEs that fit into each bin value
binNum <- tapply(ref2$error,ref2$bin,length)
binNum[is.na(binNum)] <- 0
#New dataframe for bin data,cacluate xe and ye
bin<- data.frame(binMean, binNum)
#3. Every HPE value in a given 1 m bin has a corresponding HPEm. Each of...
#...these HPEm distances is composed of 2 elements: the error (difference...
#...between the calculated and measured position) in the X direction, and...
#...the error in the Y direction.

#I will refer to these error values as Xe and Ye, respectively.
#error in the x direction
ref2$xe <- sqrt((tag2Easting-(ref2East))^2)
#error in the y direction
ref2$ye <- sqrt((tag2Northing-(ref2North))^2)
#4. For each bin, the standard deviations of Xe and Ye are calculated.
bin$xeSd <- tapply(ref2$xe,ref2$bin,sd)
bin$yeSd <- tapply(ref2$ye,ref2$bin,sd)
#5. To convert the 2-dimensional standard deviations calculated in #4 into a...
#...single measure, the 2DRMS error is calculated from the standard...
#...deviations of Xe and Ye
bin$RMS2d <- 2*sqrt((bin$xeSd)^2 + (bin$yeSd)^2)
#6. Now create a line plot and a dataframe just for the...
#...numbers we need, (ie when we have at least 10 tag...
#...transmissions and an HPE less than 21)
bin$avgHPE <- tapply(ref2$HPE,ref2$bin,mean)
smallBin <- bin[ which(bin$binNum > 10),]
smallBin <- smallBin[which(smallBin$avgHPE < 25 ),]
res3 <- lm(smallBin$RMS2d ~ smallBin$avgHPE)

#PLOT THE POINT SPREAD
#Changing a bunch of theme elements to make it pretty.

xyplot<- data.frame(tag2Easting, tag2Northing)
plot1a2<- qplot(tag2Easting, tag2Northing , alpha= 1/100,data = xyplot, xlab='Easting (m)', ylab= 'Northing (m) ') + theme_bw()
plot1a2<- plot1a2 + geom_point(data=xyplot,aes(ref2East,ref2North),col="white", shape = 18, size= 4)
plot1a2<- plot1a2 + theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 1, face = "plain",family = "serif"),
        axis.text.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = 0, face = "plain",family = "serif"),  
        axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = 1, face = "plain",family = "serif"),
        legend.position = c(0.12,0.98), legend.title = element_text(),legend.text = element_text(size = 18, family = "serif"),
        plot.title = element_text(face = "bold", size = 16, family = "serif", hjust = .5),
        plot.margin = margin(.3,.3,.3,.3, "cm"),
        plot.background = element_rect(colour = "black", fill=NA, size=1.5))
plot1a2<- plot1a2+theme(legend.position="none")#remove legend

plot1a2<- plot1a2 + coord_cartesian(xlim=c(318350 ,318450),
                                  ylim=c(5222425,5222500))
#Graph of HPE to error relationship
plot2a2<- qplot(ref2$HPE, ref2$error , alpha= 1/100,
               data = xyplot, xlab= 'HPE', ylab= 'Measured Error (m)')+
  scale_y_continuous(limits = c(0, 32), breaks = seq(0,32, by = 4)) +
  scale_x_continuous(limits = c(8, 30)) + 
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 1, face = "plain",family = "serif"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = 0, face = "plain",family = "serif"),  
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = 1, face = "plain",family = "serif"),
        legend.position = c(0.12,0.98), legend.title = element_text(),legend.text = element_text(size = 18, family = "serif"),
        plot.title = element_text(face = "bold", size = 16, family = "serif", hjust = .5),
        plot.margin = margin(.2,.2,.2,.2, "cm"),
        plot.background = element_rect(colour = "black", fill=NA, size=1.5))
plot2a2<- plot2a2 +theme(text = element_text(size=16),
                       plot.title = element_text(vjust= 2), axis.title.y = element_text(vjust= -0.02),
                       axis.title.x = element_text(vjust= -0.02))
plot2a2 <- plot2a2 + geom_point(data=smallBin, shape= 19,
                              size= 5,col='black', bg= 'black',alpha= 1, aes(x = avgHPE, y = RMS2d)) +
  theme(legend.position = "none")
plot2a2<- plot2a2 + geom_point(data=smallBin, shape= 19, size= 4.2,col='white',
                             alpha= 1, aes(x = avgHPE, y = RMS2d)) + theme(legend.position = "none")
plot2a2<- plot2a2 + geom_point(data=smallBin, shape= 4, size = 3,color='red',
                             aes(x = avgHPE, y = RMS2d)) +
  theme(legend.position = "none")

plot2a2<- plot2a2 +geom_abline(data=res3, col='white',
                             aes(intercept=res3$coefficients[1], slope=res3$coefficients[2] ), size=1)
plot2a2<- plot2a2 +geom_abline(data=res3, col='black', linetype='dashed',
                             aes(intercept=res3$coefficients[1], slope=res3$coefficients[2]), size= 1)
#plot2a2

#Add model equation and r2 to graph
#function to add linear model to graph
#Setup for a linear model to be added to a ggplot
lm_eqn = function(m) {
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  if (coef(m)[2] >= 0) {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","
                     ~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","
                     ~~italic(r)^2~"="~r2,l)
  }
  as.character(as.expression(eq));
}
#calls in linear model function for linear model res3
plot2a2=plot2a2 + annotate("text",x = 15, y = 31,
                         label = lm_eqn(res3), size= 5, parse=TRUE,family = "serif")
plot2a2

#What type of HPE cutoff should we use if we know what level of error we want
#!#! YOU MUST CHANGE accGoal to your desired error
#now evaluate how HPE cutoffs work if we want a specific error value
HPE<- c(5:HPEHigh)# we want a list of all possible HPE cutoffs
AllHPE <- data.frame(HPE)
AllHPE$ptsReject <- NA# number of points removed
AllHPE$ptsRejectP <- NA# percetage of all ponits
AllHPE$ptsRetain <- NA# number of poitns retained
AllHPE$ptsRetainP <- NA# percetage of all ponits
AllHPE$incorrectReject <- NA# good= error < accGoal and unacceptably...
#...erroneous= error > accGoal
AllHPE$incorrectRejectP <- NA
AllHPE$incorrectRetain<- NA
AllHPE$incorrectRetainP <- NA
AllHPE$incorrectRetainPvsRetain <-NA # incorretly retained/all retained *100
AllHPE$correctReject <- NA
AllHPE$correctRejectP <- NA

AllHPE$correctRetain <- NA
AllHPE$correctRetainP <- NA
AllHPE$goodDataLossP <- NA# incorrectReject/goodpts *100
AllHPE$badDataRetainP <- NA# incorrectRetain/AllHPE$ptsRetain * 100
AllHPE$sdDEV<-NA
AllHPE$avgErr<-NA
AllHPE$maxErr <- NA

# Number of all points that were unacceptably erroneous
badpts <-length(ref2$error[ref2$error > accGoal ])
# % of all points were unacceptably erroneous
badptsP <-(length(ref2$error[ref2$error > accGoal ])/(length(
  + ref2$error)))*100
#Number of all points that are good
goodpts <-length(ref2$error[ref2$error <= accGoal ])
#% of all points that are good
goodptsP <- (length(ref2$error[ref2$error <= accGoal ])/(length(ref2$error)))*100

for( i in 1: length(AllHPE$HPE)){AllHPE$ptsReject[i] <- length( which(ref2$HPE > AllHPE$HPE[i]))
AllHPE$ptsRejectP[i] <- length( which(ref2$HPE > AllHPE$HPE[i]))/(length(ref2$error))*100 
#% of all points dropped
AllHPE$ptsRetain[i] <- length( which(ref2$HPE < AllHPE$HPE[i]))
AllHPE$ptsRetainP[i] <- length( which(ref2$HPE < AllHPE$HPE[i]))/(length(ref2$error))*100 
#% of all points dropped
#incorrect reject
AllHPE$incorrectReject[i] <- length( which( (ref2$HPE[ ref2$error <= accGoal ] > AllHPE$HPE[i]) == TRUE))
#divided by total with acceptable error
AllHPE$incorrectRejectP[i] <-(AllHPE$incorrectReject[i]/length(ref2$error[ref2$error <= accGoal ]) )*100
#incorrect retain
AllHPE$incorrectRetain[i] <-length( which((ref2$HPE[ ref2$error > accGoal ] < AllHPE$HPE[i]) == TRUE))
#divided by accurate points
AllHPE$incorrectRetainP[i] <-(AllHPE$incorrectRetain[i]/length(ref2$error[ref2$error > accGoal ]))*100
#divided by total retained by filter
AllHPE$incorrectRetainPvsRetain[i] <-(AllHPE$incorrectRetain[i]/AllHPE$ptsRetain[i] )*100
#Correctly Rejected
AllHPE$correctReject[i] <-(length(which((ref2$HPE[ref2$error > accGoal] > AllHPE$HPE[i]) == TRUE)))
#divided by sum of all positions
AllHPE$correctRejectP[i]<- (AllHPE$correctReject[i]/AllHPE$ptsRetain[i])*100
#Correctly Retained
AllHPE$correctRetain[i] <- length(which((ref2$HPE[ ref2$error < accGoal] < AllHPE$HPE[i]) == TRUE))
#divided by sum of all positions
AllHPE$correctRetainP[i]<-(AllHPE$correctRetain[i]/length(ref2$error))*100
AllHPE$goodDataLossP[i] <- (AllHPE$incorrectReject[i] / goodpts)*100 
#incorrectReject/goodpts *100
AllHPE$badDataRetainP[i] <- (AllHPE$incorrectRetain[i]/AllHPE$ptsRetain[i])*100 
# incorrectRetain/AllHPE$ptsRetain * 100
AllHPE$avgErr[i] <- mean(ref2$error[ ref2$HPE < AllHPE$HPE[i]]) 
#mean of all error for each hpe
AllHPE$sdDEV [i]<- sd(ref2$error[ ref2$HPE < AllHPE$HPE[i]]) 
#Standard deviation of all error for each HPE
AllHPE$maxErr[i] <- max(ref2$error[ ref2$HPE < AllHPE$HPE[i]]) 
#mean of all error for each hpe
AllHPE$count15[i] <- length(which( (ref2$HPE[ ref2$error > 15] < AllHPE$HPE[i]) == TRUE ))
}
start <- HPELow
end <- HPEHigh
AllHPEperf1 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$correctRejectP[start:end])
AllHPEperf1$var <-"Correct Rejection"
AllHPEperf2 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$incorrectRetainP[start:end])
AllHPEperf2$var <-"Incorrect Retainment"
AllHPEperf3 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$incorrectRejectP[start:end])
AllHPEperf3$var <-"Incorrect Rejection"
#Combine all seconts for graph
AllHPEperf <- rbind(AllHPEperf1,AllHPEperf3,AllHPEperf2)

# PLOT FOR INDICATING HPE VALUES AND THEIR RATE OF INCORRECT REJECTION AND RETAINMENT
#note legend is alphabetical but graph is by dataframe order
plot3a2 <- qplot( AllHPE$incorrectRetainPvsRetain, AllHPE$incorrectRejectP,size=3, xlab= '% Incorrecly Retained (of all retained positions)',ylab= ' % Incorrectly Rejected (of all acceptable positions)' ) + 
  theme_bw() + theme( plot.margin = unit(c(1,1.5,1,1.5), "cm"),panel.grid.major=element_blank(), panel.grid.minor=element_blank())

plot3a2<- plot3a2 +geom_line(data=AllHPE, aes(x=AllHPE$incorrectRetainPvsRetain, y=AllHPE$incorrectRejectP),color= 'black', size=0.25)
plot3a2<- plot3a2 + theme_bw() + theme( text = element_text(size=15), #plot.title = element_text(vjust= 2),
                                      axis.text.y = element_text(angle=90, hjust= 0.4),
                                      axis.title.y = element_text(vjust= .3),
                                      axis.title.x = element_text(vjust= -0.02 ),
                                      plot.margin = unit(c(1,1.5,1,1.5), "cm"),
                                      panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                      legend.position = c(1.5, 0.8),
                                      legend.background = element_rect(fill = "transparent",colour = NA),
                                      legend.title=element_blank() )
plot3a2 <-plot3a2+ geom_text(cex=6,aes(x=(AllHPE$incorrectRetainPvsRetain[3:15]),
                                     y=(AllHPE$incorrectRejectP[3:15])+4 , label=AllHPE$HPE[3:15], group=NULL),)                      

# PLOT FOR AVERAGE ERROR IN METERS AND THE ASSOCIATED HPE VALUE
plot4a2 <- ggplot( data= AllHPE[start:end,],
                  aes(x = AllHPE$HPE[start:end])) +
  geom_point(aes(y = AllHPE$avgErr[start:end],
                 shape= "Mean Error" ), size= 4, color= 'black') +
  geom_point(aes(y = AllHPE$maxErr[start:end],
                 shape="Maximum Error" ), size= 4, color='black') +
  geom_line(aes(y = OutlierGoal, ), size= 1, lty= 'dashed') +
  theme_bw() + theme( text = element_text(size=25),
                      plot.title = element_text(vjust= 2),
                      axis.text.y = element_text(angle=90, hjust= 0.5),
                      axis.title.y = element_text(vjust= -0.02 ),
                      axis.title.x = element_text(vjust= -0.02 ),
                      plot.margin = unit(c(1,1.5,1,1.5), "cm"),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      legend.position = c(0.75, 0.75),
                      legend.background = element_rect(
                        fill = "transparent",colour = NA),
                      legend.title=element_blank() ) 
  xlab('HPE' ) +
  coord_cartesian( ylim=c(0,26)) +
  scale_y_continuous("Error (m)",breaks= round(seq(min(0),
                                                   max(26), by = 1),1)) +
  scale_x_continuous("HPE",breaks=round(seq(min(AllHPE$HPE),
                                            max(AllHPE$HPE), by = 1),1)) +
  geom_text(cex=6,aes(x=AllHPE$HPE[1:25],y= (29.25),
                      label=AllHPE$count15[1:25], group=NULL),)
plot4a2

###############################################################################
########## DO NOT RUN THIS SECTION ############################################
###############################################################################

#with(tagDrag, eucDist-eucDist2)
#badHPE<- tagDrag$HPE[tagDrag$HPE >= HPEfilter ]
#badeucDist3 <-tagDrag$eucDist2[tagDrag$HPE >= HPEfilter]
#badHPENorthing<- tagDrag$Northing[tagDrag$HPE >= HPEfilter ]
#badHPEEasting<- tagDrag$Easting[tagDrag$HPE >= HPEfilter ]
#badErrorNorthing<- tagDrag$Northing[tagDrag$eucDist2 >= accGoal ]
#badErrorEasting<- tagDrag$Easting[tagDrag$eucDist2 >= accGoal ]
#eucDist3<- tagDrag$eucDist2[tagDrag$eucDist2 >= accGoal ]
#incorrectlyRejected<- (tagDrag[tagDrag$eucDist2 <6 & tagDrag$HPE>8,])
#correctlyRejected<- na.omit(tagDrag[tagDrag$eucDist2 >6 & tagDrag$HPE>8,],)
#incorrectlyRetained<- na.omit(tagDrag[tagDrag$eucDist2 >6 & tagDrag$HPE<8,],)

#correctlyRetained <- na.omit(tagDrag[tagDrag$eucDist2 < 6 & tagDrag$HPE<8,],)
#badError<- data.frame("Easting"=badErrorEasting,"Northing"=badErrorNorthing,
#                      + "eucDist3" = eucDist3, "HPE" = tagDrag$HPE[tagDrag$eucDist2 >= accGoal])
#badError <- na.omit(badError)
#New frame here
#tag<-tagDrag[ is.na(tagDrag$gpsEasting)==FALSE,]

# identify unacceptably erroneous points in the plot(must change HPE)
#plot5a2<- qplot(correctlyRetained$gpsEasting2, correctlyRetained$gpsNorthing2,
#               + xlab= ' Easting (m)', ylab= 'Northing (m) ') + theme_bw() + geom_point(data=
#                                                                                          + coast,aes(coast$Easting,coast$Northing),col="black", shape = 18, size= 5)
#plot5a2<- plot5a2 +theme(text = element_text(size=20),
#                       + plot.title = element_text(vjust= 2),
#                       + axis.text.y = element_text(angle=90, hjust= 0.5),
#                       + axis.title.y = element_text(vjust=0.4 ),
#                       + axis.title.x = element_text(vjust= -0.01 ),
#                       + legend.position="none" ) #make the axis
#receivers
#plot5a2<- plot5a2 + geom_point(data=rec2010, aes(x=Easting, y=Northing), shape= 3)
#all locations with a high HPE
#plot5a2 <- plot5a2 + geom_point(data= incorrectlyRejected,aes(x=Easting,
#                                                            + y=Northing), shape = 19, size= 5, col='Orange') # circled orange
#plot5a2 <- plot5a2 + geom_point(data= correctlyRejected,aes(x=Easting,
#                                                          + y=Northing), shape = 19, size= 5, col='Blue') # circled blue
#plot5a2 <- plot5a2 + geom_point(data= incorrectlyRetained,aes(x=Easting,
#                                                            + y=Northing), shape = 19, size= 5, col='red') # circled red
#plot5a2 <- plot5a2 + geom_point(data=correctlyRetained,aes(x=Easting,y=Northing),
#                              + shape = 19, size= 5, col='black') # circled black
#plot5a2<- plot5a2 +coord_cartesian(xlim=c(xAxisMin,xAxisMax ),
#                                 + ylim=c(yAxisMin,yAxisMax))
#plot5b<- qplot(correctlyRetained$HPE, correctlyRetained$eucDist2 , alpha= 1/100,
#               + xlab= 'HPE', ylab= 'Euclidian distance (m)' )+ scale_y_continuous(limits =
#                                                                                     + c(0,21)) + scale_x_continuous(limits = c(0, 21)) + theme_bw() + theme(
#                                                                                       + plot.margin = unit(c(1,1.5,1,1.5), "cm"), panel.grid.major=element_blank(),
#                                                                                       + panel.grid.minor=element_blank())+ theme(legend.position="none") +
#  + scale_colour_identity()
#plot5b<- plot5b +theme(text = element_text(size=20),
#                       + plot.title = element_text(vjust= 2),
#                       + axis.title.y = element_text(vjust= -0.02 ),
#                       + axis.title.x = element_text(vjust= -0.02 ) ) #make the axis
#plot5b<- plot5b + geom_point(data= incorrectlyRejected,aes(x=HPE,y=eucDist2,
#                                                           + col="orange"), alpha=1,shape = 19, size= 4)
#plot5b<- plot5b + geom_point(data= correctlyRejected,aes(x=HPE, y=eucDist2,
#                                                         + col="blue"), shape = 19, size= 4)
#plot5b<- plot5b + geom_point(data= incorrectlyRetained,aes(x=HPE, y= eucDist2,
#                                                           + col="red"), shape = 19, size= 4)

#Area of Receiver Coverage
tagll <- data.frame(X = rec$LONG, Y = rec$LAT)
attr(tagll, "projection") <- "LL"
xyrec <- convUL(tagll) * 1000
rec$EASTING<- xyrec[,1]
rec$NORTHING<- xyrec[,2]

#Converting fish location to UTM
tagll <- data.frame(X = Trial$LON, Y = Trial$LAT)
attr(tagll, "projection") <- "LL"
xyrec <- convUL(tagll) * 1000
Trial$EASTING<- xyrec[,1]
Trial$NORTHING<- xyrec[,2]

#Evaluate Acutal fish trackin data (HPE FILTER)
plot6a2 <- qplot(Trial$EASTING, Trial$NORTHING, xlab= ' Easting (m)', ylab= 'Northing (m)', size=I(1)) + theme_bw()
plot6a2<- plot6a2 + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                  ylim=c(yAxisMin,yAxisMax))
plot6a2<- plot6a2 +theme(text = element_text(size=25),
                       plot.title = element_text(vjust= 2),
                       axis.text.y = element_text(angle=90, hjust= 0.5),
                       axis.title.y = element_text(vjust= -0.02 ),
                       axis.title.x = element_text(vjust= -0.02 ),
                       plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot6a2<- plot6a2 + geom_point(data=Trial[Trial$HPE < 3,], aes(x=EASTING,
                                                             y=NORTHING),col="blue", shape = 20, size= I(1))
plot6a2<- plot6a2 + geom_point(data=Trial[Trial$HPE >= 3 &
                                          Trial$HPE < 8 ,],aes(x=EASTING, y=NORTHING),col="#PINK",
                             shape = 20, size= I(1))
plot6a2<- plot6a2 + geom_point(data=Trial[Trial$HPE >= 8 &
                                          Trial$HPE < 12 ,],aes(x=EASTING, y=NORTHING),col="green",
                             shape = 20, size= I(1))
plot6a2<- plot6a2 + geom_point(data=Trial[Trial$HPE >= 12 &
                                          Trial$HPE < 17 ,],aes(x=EASTING, y=NORTHING),col="yellow",
                             shape = 20,size= I(1))
plot6a2<- plot6a2 + geom_point(data=Trial[Trial$HPE >= 17 &
                                          Trial$HPE < 22 ,],aes(x=EASTING, y=NORTHING),col="#FFCC33",
                             shape = 20,size= I(1))
plot6a2<- plot6a2 + geom_point(data=Trial[Trial$HPE >= 22 ,],aes(x=EASTING,
                                                               y=NORTHING),col="red", shape = 20, size= I(1))
plot6a2<- plot6a2 + geom_point(data=rec,aes(EASTING,NORTHING),fill="black",
                             col="white", shape = 24, size= 3, alpha=1)
plot6a2<- plot6a2 + geom_point(data=bm_outline,aes(e,n),col="black",
                             shape = 18, size= 1)

xydata<-data.frame(X=Trial$EASTING, Y=Trial$NORTHING)

#create a polygon that covers the area of the data
fishPoly<- calcConvexHull(xydata, keepExtra=FALSE)
#Polygon(fishPoly)
#calcArea, calcCentroid are likley both useful
#Create a column in pos data.frame whethere in or outside array
Trial$array<- point.in.polygon(Trial$EASTING,Trial$NORTHING, rec$EASTING,
                               + rec$NORTHING, mode.checked= TRUE)
#plot of array polygon, with fish positions removed, 0 means not in array

#2 Filter data FOR FIRST SET OF PLOTS
cutoff <- HPEfilter
posDrop <- nrow(Trial[(Trial$HPE > cutoff),])
posDropArrayT <- ((Trial$HPE[Trial$array != 0]) > cutoff)
posDropArray <- length(subset(posDropArrayT, posDropArrayT == TRUE))
DropOutT <- ((Trial$HPE[Trial$array == 0]) > cutoff)
posDropOut <- length(subset( DropOutT, DropOutT == TRUE))
posTotal1<- nrow(Trial)
posTotal2<- length(Trial$array[Trial$array == 0])
posTotal3<- length(Trial$array[Trial$array != 0])
percDropFil1 <- (posDrop/posTotal1)*100 # % DroppedTota;
percADropFil1 <- (posDropArray/posTotal2)*100 #% dropped in the array
percODropFil1 <- (posDropOut/posTotal3)*100 #% dropped outside of the array

#Now identify the silhoutte of the point cluster using calcConvexHull
xydata1<-data.frame(X=Trial$EASTING[Trial$HPE < cutoff ], Y=Trial$NORTHING[Trial$HPE <
                                                                             + cutoff ])
#create a polygon that covers the area of the data
fishPoly1<- calcConvexHull(xydata1, keepExtra=FALSE)
 head(xydata1)
 points(fishPoly1$X, fishPoly1$Y, col='red')
 plot(Trial$EASTING,Trial$NORTHING)
 #calcArea, calcCentroid are likley both useful
 areaFil1 <- calcArea(fishPoly1) #calculates area of polygon

#Plot 1- all fish positions
plot7a2 <- qplot(Trial$EASTING, Trial$NORTHING, xlab= 'Easting (m)', ylab=
                  'Northing (m)',alpha=1/100, size=I(1)) + theme_bw()
plot7a2 <- plot7a2 + coord_cartesian(xlim=c(xAxisMin,xAxisMax), ylim=
                                     c(yAxisMin,yAxisMax))
plot7a2 <- plot7a2 + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                              shape = 17, size= 2, alpha=1)
plot7a2 <- plot7a2 + theme(text = element_text(size=13),
                        plot.title = element_text(vjust= 2),
                        axis.text.y = element_text(angle=90, hjust= 0.5),
                        axis.title.y = element_text(vjust= -0.02 ),
                        axis.title.x = element_text(vjust= -0.02 ),
                        plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7a2 <- plot7a2 + theme(legend.position="none") #remove legend
#plot 2- all positions with HPE < cutoff
plot7b2 <- qplot( Trial$EASTING[Trial$HPE<= cutoff], Trial$NORTHING[Trial$HPE<= cutoff
], xlab= ' Easting (m)', ylab= 'Northing (m) ', alpha=1/100, size=I(1)) +
  theme_bw()
plot7b2<- plot7b2 + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                  ylim=c(yAxisMin,yAxisMax))
plot7b2<- plot7b2 + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                             shape = 17, size= 2, alpha=1)
plot7b2<- plot7b2 +theme(text = element_text(size=13),
                       plot.title = element_text(vjust= 2),
                       axis.text.y = element_text(angle=90, hjust= 0.5),
                       axis.title.y = element_text(vjust= -0.02 ),
                       axis.title.x = element_text(vjust= -0.02 ),
                       plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7b2<- plot7b2 + theme(legend.position="none") #remove legend
#plot 3- all positions with HPE > cutoff)
plot7c2 <- qplot( Trial$EASTING[Trial$HPE> cutoff], Trial$NORTHING[Trial$HPE> cutoff],
                 xlab= ' Easting (m)', ylab= 'Northing (m) ',alpha=1/100, size=I(1)) +
  theme_bw()
plot7c2<- plot7c2 + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                  ylim=c(yAxisMin,yAxisMax))
plot7c2<- plot7c2 + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                             shape = 17, size= 2, alpha=1)
plot7c2<- plot7c2 +theme(text = element_text(size=13),
                       plot.title = element_text(vjust= 2),
                       axis.text.y = element_text(angle=90, hjust= 0.5),
                       axis.title.y = element_text(vjust= -0.02 ),
                       axis.title.x = element_text(vjust= -0.02 ),
                       plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7c2<- plot7c2 + theme(legend.position="none") #remove legend 

################################################################################
################################################################################
#REF TAG 3 APRIL DATA
################################################################################
# Reference tag 3 appears to be in poor area of array coverage resulting in a poor accuracy for its tag
# Nothing we can do not of course but filtering data might have to be adjusted appropriately

# RUN THE OBJECT YOU WANT TO EXTRACT DATA BY DATE WHEN YOU WANT TO CALCULATE HPE TO HPEm FOR EACH DATE RANGE
ref3 <- read_csv("E:/BSUGradBurbot/DATA/HPE_Data/TRANSMITTER-BMRef03-CALC-POSITIONS.csv")
str(ref3)
ref3 <- ref3[ref3$DATETIME >= "2019-04-10" & ref3$DATETIME <= "2019-04-30",]
ref3 <- ref3[ref3$DATETIME >= "2019-06-20" & ref3$DATETIME <= "2019-07-10",]
ref3 <- ref3[ref3$DATETIME >= "2019-12-21" & ref3$DATETIME <= "2020-01-11",]

#Known Easting and Northing (UTM) for fixed tag location if known...
#...otherwise leave as NA
ref3East <-NA
ref3North <-NA
#tag transmission rate
transRateStat <- 600 #stationary tag average transmission rate
#transRateMobile <- 980 #tag drag tag test average transmission rate
#Select HPE filter, can be left at default until evaluating the array
HPEfilter <- 15
#Map Ploting Limits , a consistent cord-cartesian cutoff for each grid, look...
#...at study site and identify minimum and maximum easting and Northing
xAxisMin <-316153
xAxisMax <-320248  
yAxisMin <-5219142
yAxisMax <-5224226
#Error Range of interest, just leave this for most applications
ErrorRange <- c(5:30)

#Accuracy Goals: These are deterimined by your reserach questions and
#...analysis needs.
#Oultier Maximum
OutlierGoal <- 15

#Accuracy Goal for all Points
accGoal <- 6

#Goal for average accuracy if studying path lengths (1 order of magnitude...
#...less than max path length)
avgAccGoal<- 2.77

#Setup for Evaluations
HPELow <- 5 #Lowest HPE to consider in plot.
HPEHigh<- 30 #Higest HPE to consider in plot.

#Length of Tag Drag in Seconds
#tagDragLength <- 4860 #recorded during the test

#UTC time zone difference to correct for variable times
#timeZoneDiff=4*60*60

# Convert Latitude longitude to UTM
tagll <- data.frame(X = ref3$LON, Y = ref3$LAT)
attr(tagll, "projection") <- "LL"
xy <- convUL(tagll) * 1000
#Rename Postions
tag3Easting<- xy[,1]
tag3Northing<- xy[,2]
ref3$Easting<- tag3Easting
ref3$Northing<- tag3Northing

# Median Point if a measured point is unavailable
ref3East <- ifelse(is.na(ref3East), median(ref3$Easting),
                   + ref3East)
ref3North <- ifelse(is.na(ref3North), median(ref3$Northing),
                    + ref3North)
#calculate distance between each point and the average
ref3$error <- sqrt( (tag3Easting-(ref3East))^2 +
                         + (tag3Northing-(ref3North))^2)

#What is the fix rate based on transmit rate and test length
#convert date time
tagDateTime<- as.POSIXct((ref3$DATETIME), tz='UTC')
testLength<- as.numeric(max(tagDateTime))-as.numeric(min(tagDateTime))

#number of transmissions that should occur
possibleFixes <- (as.numeric(max(tagDateTime))-
                    + as.numeric(min(tagDateTime)))/transRateStat
#total number of transmisisons heard
actualFixes <- length(ref3$DATETIME)
#Fix Rate
fixRateStat <- (actualFixes/possibleFixes)*100

################################################################################
########## DO NOT RUN ANY OF THIS CODE #########################################
################################################################################
#Process tag drag data
#tagDragLat<-tagDrag$LAT
#tagDragLon<-tagDrag$LON
#tagDragDateTime<-tagDrag$DATETIME
#tagDragDateTime<-as.POSIXct((tagDragDateTime),tz='UTC')-timeZoneDiff
#Note that I subtract sections to standardize time (4hrs)
#tagDragHPE<-tagDrag$HPE
#tagDragdepth<-tagDrag$DEPTH
#tagDragDRX<-tagDrag$DRX
#convert latitude and longitude to UTM for tag drag
#tagll <- data.frame(X = tagDragLon, Y = tagDragLat)
#attr(tagll, "projection") <- "LL"
#xy <- convUL(tagll) * 1000
#Easting and Northing for tags
#tagDragEasting<- xy[,1]
#tagDragNorthing<- xy[,2]
#Format receiver locations
#preparing receiver data for use in plots
#tagll <- data.frame(X = rec2010$Longitude, Y=rec2010$Latitude)
#attr(tagll, "projection") <- "LL"
#xyrec <- convUL(tagll) * 1000
#rec2010$Easting<- xyrec[,1]
#rec2010$Northing<- xyrec[,2]
#Create a new tag dataframe to work from these variables
#tagDrag <-data.frame("Easting"= tagDragEasting, "Northing" =
#                           + tagDragNorthing,"dateTime" = tagDragDateTime, "HPE" =tagDragHPE,
#                         + "depth" = tagDragdepth , "DRX" = tagDragDRX)
#Pull out necessary columns from the GPS file
#gpsLat<- gpsDrag$Latitude
#gpsLon<- gpsDrag$Longitude
#gpsTime<-gpsDrag$Time
#gpsDate<-gpsDrag$Date
#gpsDateTime <- strptime(paste(gpsDate, gpsTime),
#                      + format="%m/%d/%Y %H:%M:%S", tz= 'UTC')
#Convert Gps positions from Lat./Lon. to UTM
#gpsll <- data.frame(X = gpsLon, Y = gpsLat)
#attr(gpsll, "projection") <- "LL"
#xy <- convUL(gpsll) * 1000
#gpsEasting<-xy[,1]
#gpsNorthing<-xy[,2]
#FIX RATE
#number of transmissions that should occur (drag times 13:12:59-13:57:24,
#...14:22:19-15:10:08, 15:35:08-15:59:42) which is 2665 s, 709 s, 1474 s,
#...corresponding to 162 expected transmissions.
#possibleFixesMobile <- tagDragLength/transRateMobile
#total number of transmisisons heard
#actualFixesMobile <- length(tagDrag$dateTime)
#Fix Rate
#fixRateMobile <- (actualFixesMobile/possibleFixesMobile)
#Estimating tag drag accuracy
# Estimate your true position from your GPS data
#tagDrag$gpsEasting <- approx(x=gpsDateTime,
#                          + y=gpsEasting, xout=tagDrag$dateTime)$y
#tagDrag$gpsNorthing <- approx(x=gpsDateTime,
#                            + y=gpsNorthing, xout=tagDrag$dateTime)$y
#5 Calculate the Euclidian distance between esimated True and actual measured...
#...tag positions
#tagDrag$eucDist <- sqrt(((tagDrag$Easting-tagDrag$gpsEasting)^2)+
#                             + ((tagDrag$Northing-tagDrag$gpsNorthing)^2))
#gps <- data.frame(DateTime= gpsDateTime, Easting=gpsEasting,
#                   + Northing=gpsNorthing)
#estimateTimeOffset <- function(tOffset, tagDrag, gps){
+ # Estimate your true position from your GPS data
  #    + tagDrag$gpsEasting <- approx(x=gps$DateTime, y=gps$Easting,
  #                                   + xout=tagDrag$dateTime+tOffset)$y
  #    + tagDrag$gpsNorthing <- approx(x=gps$DateTime,y=gps$Northing,
  #                                    + xout=tagDrag$dateTime+tOffset)$y
  #    +
  #      + tagDrag$eucDist <- sqrt((tagDrag$Easting-tagDrag$Easting)^2 +
  #                                  + ((tagDrag$Northing-tagDrag$gpsNorthing)^2))
  #      +
  #        + return(sum(tagDrag$eucDist^2,na.rm=T))
  #      + }
  #use optim to find tOffset that minimizes the sum of squared error
#offset.optim <- optim(0, fn=estimateTimeOffset, tag=tagDrag, gps=gps,
#                          + method="BFGS")
#tOffset <- offset.optim$par
#6 now calculate gps positions based on adusted tag position timestamps:
# this is because the Vemco times are not a perfect match.
# Estimate your true position from your GPS data
#tagDrag$gpsEasting2 <- approx(x=gpsDateTime, y=gps$Easting,
#                                  + xout=tagDrag$dateTime+tOffset)$y
#tagDrag$gpsNorthing2 <- approx(x=gpsDateTime,y=gps$Northing,
#                                 + xout=tagDrag$dateTime+tOffset)$y
#tagDrag$eucDist2 <- sqrt(((tagDrag$Easting-tagDrag$gpsEasting2)^2)+
#                             + ((tagDrag$Northing-tagDrag$gpsNorthing2)^2))
#Output a list of variables for accessing Array using stationary tags
################################################################################
################################################################################
#For Stationary Test
#Fix Rate %
fixRateStat

#Average Accuracy (m)
avgAccStat <- mean(ref3$error)
avgAccStat

#Median Accuracy (m)
medianAccStat <- median(ref3$error)
medianAccStat

#Proportion with goal accuracy
propAccStat<-length(which( ref3$error < accGoal ))/
  + length(ref3$error)
propAccStat

#Number of points with greater than outlier size
outlierCountStat <- length(which(ref3$error > OutlierGoal))
outlierCountStat

#For tag drag Test
#Fix Rate
#fixRateMobile

#Average Accuracy
#meanAccMobile<- mean(na.omit(tagDrag$eucDist2))
#meanAccMobile

#Median Accuracy
#medianAccMobile <- median(na.omit(tagDrag$eucDist2))
#medianAccMobile

#Proportion with Goal accuracy
#propAccMobile <-length(which( tagDrag$eucDist2 < accGoal ))/
#  + length(tagDrag$eucDist2)
#propAccMobile

#Number of points with greater than outlier size (in our study 15 m)
#outlierCountMobile <- length(which( tagDrag$eucDist2 > OutlierGoal ))
#outlierCountMobile


#For Fixed Tag first
#What percentage of positions had and error below this value
percentGood1<- data.frame(ErrorRange)
percentGood1$goodperc<-NA

for( i in 1: length(percentGood1$ErrorRange)){percentGood1$goodperc [i] <- ((length(ref3$error[ref3$error < percentGood1$ErrorRange[i]])) / (length(ref3$error)))*100}

#tag drag Tag
#percentGood2<- data.frame(ErrorRange)
#percentGood2$goodperc<-NA
#for( i in 1: length(percentGood2$ErrorRange)){
#  + #What percentage of positions had error < value
#    + percentGood2$goodperc[i] <- ((length(tagDrag$eucDist2[tagDrag$eucDist2 <
#                                                              + percentGood2$ErrorRange[i]])) / (length(tagDrag$eucDist2)))*100
#    + }


#percent good
loc1 <- percentGood1
loc1$type<-c("Stationary Test 1")
#mobile1<- percentGood2
#mobile1$type<-c("Tag Drag Test")
AllperGood <- rbind(loc1)

# PLOT FOR SHOWING PERCENTAGE OF GOOD POINTS AND THEIR ASSOCIATED MEASURE OF ERROR IN METERS
plot0a3 <- ggplot(data=AllperGood, aes(x=AllperGood$ErrorRange,y=round(AllperGood$goodperc, digits=0)))+
  geom_point(size = 4)+
  geom_line()+
  theme_bw()

#plot0a3 <- ggplot(data=AllperGood, aes(x=ErrorRange,y=round(goodperc, digits=0))) +
#  geom_point(size= 4)+geom_line()+theme_bw() +
#  theme(element_text(size=25),
#        axis.text.y = element_text(angle=90, hjust= 0.5),
#        axis.title.y = element_text(vjust= -0.02 ),
#        axis.title.x = element_text(vjust= -0.02 ),
#        plot.margin = unit(c(1,1.5,1,1.5), "cm"),
#        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
#        legend.background = element_rect(fill = "transparent", colour = NA),
#        legend.position = c(0.6, 0.5),
#        legend.title=element_blank()) +
#  scale_y_continuous("Percent of all Positions",limits=c(0,100)) +
#  scale_x_continuous("Measured Error (m)",breaks= round(seq(min(
#    AllperGood$ErrorRange), max(AllperGood$ErrorRange), by = 1),1)) +
#  coord_cartesian( xlim=c(0,20.5))

#1. calculate the average error for each 1 unit HPE Bin.
#Bin by HPE
#create increments to bin by
breaks= seq(0, max(ref3$HPE),by=1)
# create a column with bins
ref3$bin <- cut(ref3$HPE,breaks)
#2. For each 1 m bin, the average HPE is calculated("binmean")
#use taplly to bin the error and calculate a min
binMean <- tapply(ref3$error,ref3$bin,mean)
#count the number of HPEs that fit into each bin value
binNum <- tapply(ref3$error,ref3$bin,length)
binNum[is.na(binNum)] <- 0
#New dataframe for bin data,cacluate xe and ye
bin<- data.frame(binMean, binNum)
#3. Every HPE value in a given 1 m bin has a corresponding HPEm. Each of...
#...these HPEm distances is composed of 2 elements: the error (difference...
#...between the calculated and measured position) in the X direction, and...
#...the error in the Y direction.

#I will refer to these error values as Xe and Ye, respectively.
#error in the x direction
ref3$xe <- sqrt((tag3Easting-(ref3East))^2)
#error in the y direction
ref3$ye <- sqrt((tag3Northing-(ref3North))^2)
#4. For each bin, the standard deviations of Xe and Ye are calculated.
bin$xeSd <- tapply(ref3$xe,ref3$bin,sd)
bin$yeSd <- tapply(ref3$ye,ref3$bin,sd)
#5. To convert the 2-dimensional standard deviations calculated in #4 into a...
#...single measure, the 2DRMS error is calculated from the standard...
#...deviations of Xe and Ye
bin$RMS2d <- 2*sqrt((bin$xeSd)^2 + (bin$yeSd)^2)
#6. Now create a line plot and a dataframe just for the...
#...numbers we need, (ie when we have at least 10 tag...
#...transmissions and an HPE less than 21)
bin$avgHPE <- tapply(ref3$HPE,ref3$bin,mean)
smallBin <- bin[ which(bin$binNum >= 200),]
smallBin <- smallBin[which(smallBin$avgHPE < 25 ),]
res3 <- lm(smallBin$RMS2d ~ smallBin$avgHPE)

#PLOT THE POINT SPREAD
#Changing a bunch of theme elements to make it pretty.
xyplot<- data.frame(tag3Easting, tag3Northing)
plot1a3 <- qplot(tag3Easting, tag3Northing , alpha= 1/100,data = xyplot, xlab='Easting (m)', ylab= 'Northing (m) ') + theme_bw()
plot1a3 <- plot1a3 + geom_point(data=xyplot,aes(ref3East,ref3North),col="white", shape = 18, size= 4)
plot1a3 <- plot1a3 + theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = 1, face = "plain"),
        legend.position = c(0.12,0.98), legend.title = element_text(),legend.text = element_text(size = 18, family = "serif"),
        plot.title = element_text(face = "bold", size = 16, family = "serif", hjust = .5),
        plot.margin = margin(.3,.3,.3,.3, "cm"),
        plot.background = element_rect(colour = "black", fill=NA, size=1.5))
plot1a3 <- plot1a3+theme(legend.position="none")#remove legend

plot1a3a <- plot1a3 + coord_cartesian(xlim=c(317275 ,317425),
                                      ylim=c(5220950,5221250))

#Graph of HPE to error relationship
plot2a3 <- qplot(ref3$HPE, ref3$error , alpha= 1/100,
                data = xyplot, xlab= 'HPE', ylab= 'Measured Error (m)')+
  scale_y_continuous(limits = c(0, 50)) +
  scale_x_continuous(limits = c(8, 32)) + 
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = 1, face = "plain"),
        legend.position = c(0.12,0.98), legend.title = element_text(),legend.text = element_text(size = 18, family = "serif"),
        plot.title = element_text(face = "bold", size = 16, family = "serif", hjust = .5),
        plot.margin = margin(.2,.2,.2,.2, "cm"),
        plot.background = element_rect(colour = "black", fill=NA, size=1.5))
plot2a3 <- plot2a3 +theme(text = element_text(size=20),
                         plot.title = element_text(vjust= 2), axis.title.y = element_text(vjust= -0.02),
                         axis.title.x = element_text(vjust= -0.02))
plot2a3 <- plot2a3 + geom_point(data=smallBin, shape= 19,
                                size= 5,col='black', bg= 'black',alpha= 1, aes(x = avgHPE, y = RMS2d)) +
  theme(legend.position = "none")
plot2a3 <- plot2a3 + geom_point(data=smallBin, shape= 19, size= 4.2,col='white',
                               alpha= 1, aes(x = avgHPE, y = RMS2d)) + theme(legend.position = "none")
plot2a3 <- plot2a3 + geom_point(data=smallBin, shape= 4, size = 3,color='red',
                               aes(x = avgHPE, y = RMS2d)) +
  theme(legend.position = "none")

plot2a3 <- plot2a3 +geom_abline(data=res3, col='white',
                               aes(intercept=res3$coefficients[1], slope=res3$coefficients[2] ), size=1)
plot2a3 <- plot2a3 +geom_abline(data=res3, col='black', linetype='dashed',
                               aes(intercept=res3$coefficients[1], slope=res3$coefficients[2]), size= 1)
#Add model equation and r2 to graph
#function to add linear model to graph
#Setup for a linear model to be added to a ggplot
lm_eqn = function(m) {
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  if (coef(m)[2] >= 0) {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","
                     ~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","
                     ~~italic(r)^2~"="~r2,l)
  }
  as.character(as.expression(eq));
}
#calls in linear model function for linear model res3
plot2a3a = plot2a3 + annotate("text",x = 11, y = 31,
                           label = lm_eqn(res3), size= 5, parse=TRUE)

plot2a3a
#What type of HPE cutoff should we use if we know what level of error we want
#!#! YOU MUST CHANGE accGoal to your desired error
#now evaluate how HPE cutoffs work if we want a specific error value
HPE<- c(5:HPEHigh)# we want a list of all possible HPE cutoffs
AllHPE <- data.frame(HPE)
AllHPE$ptsReject <- NA# number of points removed
AllHPE$ptsRejectP <- NA# percetage of all ponits
AllHPE$ptsRetain <- NA# number of poitns retained
AllHPE$ptsRetainP <- NA# percetage of all ponits
AllHPE$incorrectReject <- NA# good= error < accGoal and unacceptably...
#...erroneous= error > accGoal
AllHPE$incorrectRejectP <- NA
AllHPE$incorrectRetain<- NA
AllHPE$incorrectRetainP <- NA
AllHPE$incorrectRetainPvsRetain <-NA # incorretly retained/all retained *100
AllHPE$correctReject <- NA
AllHPE$correctRejectP <- NA

AllHPE$correctRetain <- NA
AllHPE$correctRetainP <- NA
AllHPE$goodDataLossP <- NA# incorrectReject/goodpts *100
AllHPE$badDataRetainP <- NA# incorrectRetain/AllHPE$ptsRetain * 100
AllHPE$sdDEV<-NA
AllHPE$avgErr<-NA
AllHPE$maxErr <- NA

# Number of all points that were unacceptably erroneous
badpts <-length(ref3$error[ref3$error > accGoal ])
# % of all points were unacceptably erroneous
badptsP <-(length(ref3$error[ref3$error > accGoal ])/(length(
  + ref3$error)))*100
#Number of all points that are good
goodpts <-length(ref3$error[ref3$error <= accGoal ])
#% of all points that are good
goodptsP <- (length(ref3$error[ref3$error <= accGoal ])/(length(ref3$error)))*100

for( i in 1: length(AllHPE$HPE)){AllHPE$ptsReject[i] <- length( which(ref3$HPE > AllHPE$HPE[i]))
AllHPE$ptsRejectP[i] <- length( which(ref3$HPE > AllHPE$HPE[i]))/(length(ref3$error))*100 
#% of all points dropped
AllHPE$ptsRetain[i] <- length( which(ref3$HPE < AllHPE$HPE[i]))
AllHPE$ptsRetainP[i] <- length( which(ref3$HPE < AllHPE$HPE[i]))/(length(ref3$error))*100 
#% of all points dropped
#incorrect reject
AllHPE$incorrectReject[i] <- length( which( (ref3$HPE[ ref3$error <= accGoal ] > AllHPE$HPE[i]) == TRUE))
#divided by total with acceptable error
AllHPE$incorrectRejectP[i] <-(AllHPE$incorrectReject[i]/length(ref3$error[ref3$error <= accGoal ]) )*100
#incorrect retain
AllHPE$incorrectRetain[i] <-length( which((ref3$HPE[ ref3$error > accGoal ] < AllHPE$HPE[i]) == TRUE))
#divided by accurate points
AllHPE$incorrectRetainP[i] <-(AllHPE$incorrectRetain[i]/length(ref3$error[ref3$error > accGoal ]))*100
#divided by total retained by filter
AllHPE$incorrectRetainPvsRetain[i] <-(AllHPE$incorrectRetain[i]/AllHPE$ptsRetain[i] )*100
#Correctly Rejected
AllHPE$correctReject[i] <-(length(which((ref3$HPE[ref3$error > accGoal] > AllHPE$HPE[i]) == TRUE)))
#divided by sum of all positions
AllHPE$correctRejectP[i]<- (AllHPE$correctReject[i]/AllHPE$ptsRetain[i])*100
#Correctly Retained
AllHPE$correctRetain[i] <- length(which((ref3$HPE[ ref3$error < accGoal] < AllHPE$HPE[i]) == TRUE))
#divided by sum of all positions
AllHPE$correctRetainP[i]<-(AllHPE$correctRetain[i]/length(ref3$error))*100
AllHPE$goodDataLossP[i] <- (AllHPE$incorrectReject[i] / goodpts)*100 
#incorrectReject/goodpts *100
AllHPE$badDataRetainP[i] <- (AllHPE$incorrectRetain[i]/AllHPE$ptsRetain[i])*100 
# incorrectRetain/AllHPE$ptsRetain * 100
AllHPE$avgErr[i] <- mean(ref3$error[ ref3$HPE < AllHPE$HPE[i]]) 
#mean of all error for each hpe
AllHPE$sdDEV [i]<- sd(ref3$error[ ref3$HPE < AllHPE$HPE[i]]) 
#Standard deviation of all error for each HPE
AllHPE$maxErr[i] <- max(ref3$error[ ref3$HPE < AllHPE$HPE[i]]) 
#mean of all error for each hpe
AllHPE$count15[i] <- length(which( (ref3$HPE[ ref3$error > 15] < AllHPE$HPE[i]) == TRUE ))
}
start <- HPELow
end <- HPEHigh
AllHPEperf1 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$correctRejectP[start:end])
AllHPEperf1$var <-"Correct Rejection"
AllHPEperf2 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$incorrectRetainP[start:end])
AllHPEperf2$var <-"Incorrect Retainment"
AllHPEperf3 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$incorrectRejectP[start:end])
AllHPEperf3$var <-"Incorrect Rejection"
#Combine all seconts for graph
AllHPEperf <- rbind(AllHPEperf1,AllHPEperf3,AllHPEperf2)

# PLOT FOR INDICATING HPE VALUES AND THEIR RATE OF INCORRECT REJECTION AND RETAINMENT
#note legend is alphabetical but graph is by dataframe order
plot3a3 <- qplot( AllHPE$incorrectRetainPvsRetain, AllHPE$incorrectRejectP,size=3, xlab= '% Incorrecly Retained (of all retained positions)',ylab= ' % Incorrectly Rejected (of all acceptable positions)' ) + 
  theme_bw() + theme( plot.margin = unit(c(1,1.5,1,1.5), "cm"),panel.grid.major=element_blank(), panel.grid.minor=element_blank())

plot3a3 <- plot3a3 +geom_line(data=AllHPE, aes(x=AllHPE$incorrectRetainPvsRetain, y=AllHPE$incorrectRejectP),color= 'black', size=0.25)
plot3a3 <- plot3a3 + theme_bw() + theme( text = element_text(size=15), #plot.title = element_text(vjust= 2),
                                        axis.text.y = element_text(angle=90, hjust= 0.4),
                                        axis.title.y = element_text(vjust= .3),
                                        axis.title.x = element_text(vjust= -0.02 ),
                                        plot.margin = unit(c(1,1.5,1,1.5), "cm"),
                                        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                        legend.position = c(1.5, 0.8),
                                        legend.background = element_rect(fill = "transparent",colour = NA),
                                        legend.title=element_blank() )
plot3a3 <-plot3a3+ geom_text(cex=6,aes(x=(AllHPE$incorrectRetainPvsRetain[3:15]),
                                       y=(AllHPE$incorrectRejectP[3:15])+4 , label=AllHPE$HPE[3:15], group=NULL),)                      

# PLOT FOR AVERAGE ERROR IN METERS AND THE ASSOCIATED HPE VALUE
plot4a3 <- ggplot( data= AllHPE[start:end,],
                   aes(x = AllHPE$HPE[start:end])) +
  geom_point(aes(y = AllHPE$avgErr[start:end],
                 shape = "Mean Error" ), size= 4, color= 'black') +
  geom_point(aes(y = AllHPE$maxErr[start:end],
                 shape = "Maximum Error" ), size= 4, color='black') +
  geom_line(aes(y = OutlierGoal, ), size= 1, lty= 'dashed') +
  theme_bw() + theme( text = element_text(size=25),
                      plot.title = element_text(vjust= 2),
                      axis.text.y = element_text(angle=90, hjust= 0.5),
                      axis.title.y = element_text(vjust= -0.02 ),
                      axis.title.x = element_text(vjust= -0.02 ),
                      plot.margin = unit(c(1,1.5,1,1.5), "cm"),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      legend.position = c(0.35, 0.75),
                      legend.background = element_rect(
                        fill = "transparent",colour = NA),
                      legend.title=element_blank() ) 
#  xlab('HPE') +
#  coord_cartesian( ylim=c(0,20)) +
#  scale_y_continuous("Error (m)",breaks= round(seq(min(0),
#                                                   max(20), by = 2),1)) +
#  scale_x_continuous("HPE",breaks=round(seq(min(AllHPE$HPE),
#                                            max(AllHPE$HPE), by = 1),1)) +
#  geom_text(cex=6,aes(x=AllHPE$HPE[1:21],y= (10),
#                     label=AllHPE$count15[1:21], group=NULL),)

plot4a3

#with(tagDrag, eucDist-eucDist2)
#badHPE<- tagDrag$HPE[tagDrag$HPE >= HPEfilter ]
#badeucDist3 <-tagDrag$eucDist2[tagDrag$HPE >= HPEfilter]
#badHPENorthing<- tagDrag$Northing[tagDrag$HPE >= HPEfilter ]
#badHPEEasting<- tagDrag$Easting[tagDrag$HPE >= HPEfilter ]
#badErrorNorthing<- tagDrag$Northing[tagDrag$eucDist2 >= accGoal ]
#badErrorEasting<- tagDrag$Easting[tagDrag$eucDist2 >= accGoal ]
#eucDist3<- tagDrag$eucDist2[tagDrag$eucDist2 >= accGoal ]
#incorrectlyRejected<- (tagDrag[tagDrag$eucDist2 <6 & tagDrag$HPE>8,])
#correctlyRejected<- na.omit(tagDrag[tagDrag$eucDist2 >6 & tagDrag$HPE>8,],)
#incorrectlyRetained<- na.omit(tagDrag[tagDrag$eucDist2 >6 & tagDrag$HPE<8,],)

#correctlyRetained <- na.omit(tagDrag[tagDrag$eucDist2 < 6 & tagDrag$HPE<8,],)
#badError<- data.frame("Easting"=badErrorEasting,"Northing"=badErrorNorthing,
#                      + "eucDist3" = eucDist3, "HPE" = tagDrag$HPE[tagDrag$eucDist2 >= accGoal])
#badError <- na.omit(badError)
#New frame here
#tag<-tagDrag[ is.na(tagDrag$gpsEasting)==FALSE,]


# identify unacceptably erroneous points in the plot(must change HPE)
#plot5a3 <- qplot(correctlyRetained$gpsEasting2, correctlyRetained$gpsNorthing2,
#                + xlab= ' Easting (m)', ylab= 'Northing (m) ') + theme_bw() + geom_point(data=
#                                                                                           + coast,aes(coast$Easting,coast$Northing),col="black", shape = 18, size= 5)
#plot5a3 <- plot5a3 +theme(text = element_text(size=20),
#                         + plot.title = element_text(vjust= 2),
#                         + axis.text.y = element_text(angle=90, hjust= 0.5),
#                         + axis.title.y = element_text(vjust=0.4 ),
#                         + axis.title.x = element_text(vjust= -0.01 ),
#                         + legend.position="none" ) #make the axis
#receivers
#plot5a3 <- plot5a3 + geom_point(data=rec2010, aes(x=Easting, y=Northing), shape= 3)
#all locations with a high HPE
#plot5a3 <- plot5a3 + geom_point(data= incorrectlyRejected,aes(x=Easting,
#                                                              + y=Northing), shape = 19, size= 5, col='Orange') # circled orange
#plot5a3 <- plot5a3 + geom_point(data= correctlyRejected,aes(x=Easting,
#                                                            + y=Northing), shape = 19, size= 5, col='Blue') # circled blue
#plot5a3 <- plot5a3 + geom_point(data= incorrectlyRetained,aes(x=Easting,
#                                                              + y=Northing), shape = 19, size= 5, col='red') # circled red
#plot5a3 <- plot5a3 + geom_point(data=correctlyRetained,aes(x=Easting,y=Northing),
#                                + shape = 19, size= 5, col='black') # circled black
#plot5a3 <- plot5a3 +coord_cartesian(xlim=c(xAxisMin,xAxisMax ),
#                                   + ylim=c(yAxisMin,yAxisMax))
#plot5b3 <- qplot(correctlyRetained$HPE, correctlyRetained$eucDist2 , alpha= 1/100,
#               + xlab= 'HPE', ylab= 'Euclidian distance (m)' )+ scale_y_continuous(limits =
#                                                                                     + c(0,21)) + scale_x_continuous(limits = c(0, 21)) + theme_bw() + theme(
#                                                                                       + plot.margin = unit(c(1,1.5,1,1.5), "cm"), panel.grid.major=element_blank(),
#                                                                                       + panel.grid.minor=element_blank())+ theme(legend.position="none") +
#  + scale_colour_identity()
#plot5b3 <- plot5b3 +theme(text = element_text(size=20),
#                       + plot.title = element_text(vjust= 2),
#                       + axis.title.y = element_text(vjust= -0.02 ),
#                       + axis.title.x = element_text(vjust= -0.02 ) ) #make the axis
#plot5b3 <- plot5b3 + geom_point(data= incorrectlyRejected,aes(x=HPE,y=eucDist2,
#                                                           + col="orange"), alpha=1,shape = 19, size= 4)
#plot5b3 <- plot5b3 + geom_point(data= correctlyRejected,aes(x=HPE, y=eucDist2,
#                                                         + col="blue"), shape = 19, size= 4)
#plot5b3 <- plot5b3 + geom_point(data= incorrectlyRetained,aes(x=HPE, y= eucDist2,
#                                                           + col="red"), shape = 19, size= 4)

#Area of Receiver Coverage
tagll <- data.frame(X = rec$LONG, Y = rec$LAT)
attr(tagll, "projection") <- "LL"
xyrec <- convUL(tagll) * 1000
rec$EASTING<- xyrec[,1]
rec$NORTHING<- xyrec[,2]

#Converting fish location to UTM
tagll <- data.frame(X = Trial$LON, Y = Trial$LAT)
attr(tagll, "projection") <- "LL"
xyrec <- convUL(tagll) * 1000
Trial$EASTING<- xyrec[,1]
Trial$NORTHING<- xyrec[,2]

#Evaluate Acutal fish trackin data (HPE FILTER)
plot6a3 <- qplot(Trial$EASTING, Trial$NORTHING, xlab= ' Easting (m)', ylab= 'Northing (m)', size=I(1)) + theme_bw()
plot6a3 <- plot6a3 + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                    ylim=c(yAxisMin,yAxisMax))
plot6a3 <- plot6a3 +theme(text = element_text(size=25),
                         plot.title = element_text(vjust= 2),
                         axis.text.y = element_text(angle=90, hjust= 0.5),
                         axis.title.y = element_text(vjust= -0.02 ),
                         axis.title.x = element_text(vjust= -0.02 ),
                         plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot6a3 <- plot6a3 + geom_point(data=Trial[Trial$HPE < 3,], aes(x=EASTING,
                                                               y=NORTHING),col="blue", shape = 20, size= I(1))
plot6a3 <- plot6a3 + geom_point(data=Trial[Trial$HPE >= 3 &
                                            Trial$HPE < 5 ,],aes(x=EASTING, y=NORTHING),col="#56B4E9",
                               shape = 20, size= I(1))
plot6a3 <- plot6a3 + geom_point(data=Trial[Trial$HPE >= 5 &
                                            Trial$HPE < 10 ,],aes(x=EASTING, y=NORTHING),col="green",
                               shape = 20, size= I(1))
plot6a3 <- plot6a3 + geom_point(data=Trial[Trial$HPE >= 10 &
                                            Trial$HPE < 15 ,],aes(x=EASTING, y=NORTHING),col="yellow",
                               shape = 20,size= I(1))
plot6a3 <- plot6a3 + geom_point(data=Trial[Trial$HPE >= 15 &
                                            Trial$HPE < 20 ,],aes(x=EASTING, y=NORTHING),col="#FFCC33",
                               shape = 20,size= I(1))
plot6a3 <- plot6a3 + geom_point(data=Trial[Trial$HPE >= 20 ,],aes(x=EASTING,
                                                                 y=NORTHING),col="red", shape = 20, size= I(1))
plot6a3 <- plot6a3 + geom_point(data=rec,aes(EASTING,NORTHING),fill="black",
                               col="white", shape = 24, size= 3, alpha=1)
plot6a3 <- plot6a3 + geom_point(data=bm_outline,aes(e,n),col="black",
                               shape = 18, size= 1)

xydata <- data.frame(X=Trial$EASTING, Y=Trial$NORTHING)

#create a polygon that covers the area of the data
fishPoly<- calcConvexHull(xydata, keepExtra=FALSE)
#Polygon(fishPoly)
#calcArea, calcCentroid are likley both useful
#Create a column in pos data.frame whethere in or outside array
Trial$array<- point.in.polygon(Trial$EASTING,Trial$NORTHING, rec$EASTING,
                               + rec$NORTHING, mode.checked= TRUE)
#plot of array polygon, with fish positions removed, 0 means not in array

#2 Filter data FOR FIRST SET OF PLOTS
cutoff <- HPEfilter
posDrop <- nrow(Trial[(Trial$HPE > cutoff),])
posDropArrayT <- ((Trial$HPE[Trial$array != 0]) > cutoff)
posDropArray <- length(subset(posDropArrayT, posDropArrayT == TRUE))
DropOutT <- ((Trial$HPE[Trial$array == 0]) > cutoff)
posDropOut <- length(subset( DropOutT, DropOutT == TRUE))
posTotal1<- nrow(Trial)
posTotal2<- length(Trial$array[Trial$array == 0])
posTotal3<- length(Trial$array[Trial$array != 0])
percDropFil1 <- (posDrop/posTotal1)*100 # % DroppedTota;
percADropFil1 <- (posDropArray/posTotal2)*100 #% dropped in the array
percODropFil1 <- (posDropOut/posTotal3)*100 #% dropped outside of the array

#Now identify the silhoutte of the point cluster using calcConvexHull
xydata1<-data.frame(X=Trial$EASTING[Trial$HPE < cutoff ], Y=Trial$NORTHING[Trial$HPE <
                                                                             + cutoff ])
#create a polygon that covers the area of the data
fishPoly1<- calcConvexHull(xydata1, keepExtra=FALSE)
head(xydata1)
points(fishPoly1$X, fishPoly1$Y, col='red')
plot(Trial$EASTING,Trial$NORTHING)
#calcArea, calcCentroid are likley both useful
areaFil1 <- calcArea(fishPoly1) #calculates area of polygon

#Plot 1- all fish positions
plot7a3 <- qplot(Trial$EASTING, Trial$NORTHING, xlab= 'Easting (m)', ylab=
                   'Northing (m)',alpha=1/100, size=I(1)) + theme_bw()
plot7a3 <- plot7a3 + coord_cartesian(xlim=c(xAxisMin,xAxisMax), ylim=
                                       c(yAxisMin,yAxisMax))
plot7a3 <- plot7a3 + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                                shape = 17, size= 2, alpha=1)
plot7a3 <- plot7a3 + theme(text = element_text(size=13),
                           plot.title = element_text(vjust= 2),
                           axis.text.y = element_text(angle=90, hjust= 0.5),
                           axis.title.y = element_text(vjust= -0.02 ),
                           axis.title.x = element_text(vjust= -0.02 ),
                           plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7a3 <- plot7a3 + theme(legend.position="none") #remove legend
#plot 2- all positions with HPE < cutoff
plot7b3 <- qplot( Trial$EASTING[Trial$HPE<= cutoff], Trial$NORTHING[Trial$HPE<= cutoff
], xlab= ' Easting (m)', ylab= 'Northing (m) ', alpha=1/100, size=I(1)) +
  theme_bw()
plot7b3 <- plot7b3 + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                    ylim=c(yAxisMin,yAxisMax))
plot7b3 <- plot7b3 + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                               shape = 17, size= 2, alpha=1)
plot7b3 <- plot7b3 +theme(text = element_text(size=13),
                         plot.title = element_text(vjust= 2),
                         axis.text.y = element_text(angle=90, hjust= 0.5),
                         axis.title.y = element_text(vjust= -0.02 ),
                         axis.title.x = element_text(vjust= -0.02 ),
                         plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7b3 <- plot7b3 + theme(legend.position="none") #remove legend
#plot 3- all positions with HPE > cutoff)
plot7c3 <- qplot( Trial$EASTING[Trial$HPE> cutoff], Trial$NORTHING[Trial$HPE> cutoff],
                  xlab= ' Easting (m)', ylab= 'Northing (m) ',alpha=1/100, size=I(1)) +
  theme_bw()
plot7c3 <- plot7c3 + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                    ylim=c(yAxisMin,yAxisMax))
plot7c3 <- plot7c3 + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                               shape = 17, size= 2, alpha=1)
plot7c3 <- plot7c3 +theme(text = element_text(size=13),
                         plot.title = element_text(vjust= 2),
                         axis.text.y = element_text(angle=90, hjust= 0.5),
                         axis.title.y = element_text(vjust= -0.02 ),
                         axis.title.x = element_text(vjust= -0.02 ),
                         plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7c3 <- plot7c3 + theme(legend.position="none") #remove legend

###############################################################################
###############################################################################
################################################################################
#REF TAG 4 APRIL DATA
################################################################################

# RUN THE OBJECT YOU WANT TO EXTRACT DATA BY DATE WHEN YOU WANT TO CALCULATE HPE TO HPEm FOR EACH DATE RANGE
ref4 <- read_csv("E:/BSUGradBurbot/DATA/HPE_Data/TRANSMITTER-BMRef04-CALC-POSITIONS.csv")
str(ref4)
ref4 <- ref4[ref4$DATETIME >= "2019-04-10" & ref4$DATETIME <= "2019-04-30",]
ref4 <- ref4[ref4$DATETIME >= "2019-06-20" & ref4$DATETIME <= "2019-07-10",]
ref4 <- ref4[ref4$DATETIME >= "2019-12-21" & ref4$DATETIME <= "2020-01-11",]

#Known Easting and Northing (UTM) for fixed tag location if known...
#...otherwise leave as NA
ref4East <-NA
ref4North <-NA
#tag transmission rate
transRateStat <- 600 #stationary tag average transmission rate
#transRateMobile <- 980 #tag drag tag test average transmission rate
#Select HPE filter, can be left at default until evaluating the array
HPEfilter <- 15
#Map Ploting Limits , a consistent cord-cartesian cutoff for each grid, look...
#...at study site and identify minimum and maximum easting and Northing
xAxisMin <-316153
xAxisMax <-320248  
yAxisMin <-5219142
yAxisMax <-5224226
#Error Range of interest, just leave this for most applications
ErrorRange <- c(1:30)

#Accuracy Goals: These are deterimined by your reserach questions and
#...analysis needs.
#Oultier Maximum
OutlierGoal <- 10

#Accuracy Goal for all Points
accGoal <- 3

#Goal for average accuracy if studying path lengths (1 order of magnitude...
#...less than max path length)
avgAccGoal<- 1.77

#Setup for Evaluations
HPELow <- 1 #Lowest HPE to consider in plot.
HPEHigh<- 25 #Higest HPE to consider in plot.

#Length of Tag Drag in Seconds
#tagDragLength <- 4860 #recorded during the test

#UTC time zone difference to correct for variable times
#timeZoneDiff=4*60*60

# Convert Latitude longitude to UTM
tagll <- data.frame(X = ref4$LON, Y = ref4$LAT)
attr(tagll, "projection") <- "LL"
xy <- convUL(tagll) * 1000
#Rename Postions
tag4Easting<- xy[,1]
tag4Northing<- xy[,2]
ref4$Easting<- tag4Easting
ref4$Northing<- tag4Northing

# Median Point if a measured point is unavailable
ref4East <- ifelse(is.na(ref4East), median(ref4$Easting),
                   + ref4East)
ref4North <- ifelse(is.na(ref4North), median(ref4$Northing),
                    + ref4North)
#calculate distance between each point and the average
ref4$error <- sqrt( (tag4Easting-(ref4East))^2 +
                         + (tag4Northing-(ref4North))^2)

#What is the fix rate based on transmit rate and test length
#convert date time
tagDateTime<- as.POSIXct((ref4$DATETIME), tz='UTC')
testLength<- as.numeric(max(tagDateTime))-as.numeric(min(tagDateTime))

#number of transmissions that should occur
possibleFixes <- (as.numeric(max(tagDateTime))-
                    + as.numeric(min(tagDateTime)))/transRateStat
#total number of transmisisons heard
actualFixes <- length(ref4$DATETIME)
#Fix Rate
fixRateStat <- (actualFixes/possibleFixes)*100

################################################################################
########## DO NOT RUN ANY OF THIS CODE #########################################
################################################################################
#Process tag drag data
#tagDragLat<-tagDrag$LAT
#tagDragLon<-tagDrag$LON
#tagDragDateTime<-tagDrag$DATETIME
#tagDragDateTime<-as.POSIXct((tagDragDateTime),tz='UTC')-timeZoneDiff
#Note that I subtract sections to standardize time (4hrs)
#tagDragHPE<-tagDrag$HPE
#tagDragdepth<-tagDrag$DEPTH
#tagDragDRX<-tagDrag$DRX
#convert latitude and longitude to UTM for tag drag
#tagll <- data.frame(X = tagDragLon, Y = tagDragLat)
#attr(tagll, "projection") <- "LL"
#xy <- convUL(tagll) * 1000
#Easting and Northing for tags
#tagDragEasting<- xy[,1]
#tagDragNorthing<- xy[,2]
#Format receiver locations
#preparing receiver data for use in plots
#tagll <- data.frame(X = rec2010$Longitude, Y=rec2010$Latitude)
#attr(tagll, "projection") <- "LL"
#xyrec <- convUL(tagll) * 1000
#rec2010$Easting<- xyrec[,1]
#rec2010$Northing<- xyrec[,2]
#Create a new tag dataframe to work from these variables
#tagDrag <-data.frame("Easting"= tagDragEasting, "Northing" =
#                           + tagDragNorthing,"dateTime" = tagDragDateTime, "HPE" =tagDragHPE,
#                         + "depth" = tagDragdepth , "DRX" = tagDragDRX)
#Pull out necessary columns from the GPS file
#gpsLat<- gpsDrag$Latitude
#gpsLon<- gpsDrag$Longitude
#gpsTime<-gpsDrag$Time
#gpsDate<-gpsDrag$Date
#gpsDateTime <- strptime(paste(gpsDate, gpsTime),
#                      + format="%m/%d/%Y %H:%M:%S", tz= 'UTC')
#Convert Gps positions from Lat./Lon. to UTM
#gpsll <- data.frame(X = gpsLon, Y = gpsLat)
#attr(gpsll, "projection") <- "LL"
#xy <- convUL(gpsll) * 1000
#gpsEasting<-xy[,1]
#gpsNorthing<-xy[,2]
#FIX RATE
#number of transmissions that should occur (drag times 13:12:59-13:57:24,
#...14:22:19-15:10:08, 15:35:08-15:59:42) which is 2665 s, 709 s, 1474 s,
#...corresponding to 162 expected transmissions.
#possibleFixesMobile <- tagDragLength/transRateMobile
#total number of transmisisons heard
#actualFixesMobile <- length(tagDrag$dateTime)
#Fix Rate
#fixRateMobile <- (actualFixesMobile/possibleFixesMobile)
#Estimating tag drag accuracy
# Estimate your true position from your GPS data
#tagDrag$gpsEasting <- approx(x=gpsDateTime,
#                          + y=gpsEasting, xout=tagDrag$dateTime)$y
#tagDrag$gpsNorthing <- approx(x=gpsDateTime,
#                            + y=gpsNorthing, xout=tagDrag$dateTime)$y
#5 Calculate the Euclidian distance between esimated True and actual measured...
#...tag positions
#tagDrag$eucDist <- sqrt(((tagDrag$Easting-tagDrag$gpsEasting)^2)+
#                             + ((tagDrag$Northing-tagDrag$gpsNorthing)^2))
#gps <- data.frame(DateTime= gpsDateTime, Easting=gpsEasting,
#                   + Northing=gpsNorthing)
#estimateTimeOffset <- function(tOffset, tagDrag, gps){
+ # Estimate your true position from your GPS data
  #    + tagDrag$gpsEasting <- approx(x=gps$DateTime, y=gps$Easting,
  #                                   + xout=tagDrag$dateTime+tOffset)$y
  #    + tagDrag$gpsNorthing <- approx(x=gps$DateTime,y=gps$Northing,
  #                                    + xout=tagDrag$dateTime+tOffset)$y
  #    +
  #      + tagDrag$eucDist <- sqrt((tagDrag$Easting-tagDrag$Easting)^2 +
  #                                  + ((tagDrag$Northing-tagDrag$gpsNorthing)^2))
  #      +
  #        + return(sum(tagDrag$eucDist^2,na.rm=T))
  #      + }
  #use optim to find tOffset that minimizes the sum of squared error
#offset.optim <- optim(0, fn=estimateTimeOffset, tag=tagDrag, gps=gps,
#                          + method="BFGS")
#tOffset <- offset.optim$par
#6 now calculate gps positions based on adusted tag position timestamps:
# this is because the Vemco times are not a perfect match.
# Estimate your true position from your GPS data
#tagDrag$gpsEasting2 <- approx(x=gpsDateTime, y=gps$Easting,
#                                  + xout=tagDrag$dateTime+tOffset)$y
#tagDrag$gpsNorthing2 <- approx(x=gpsDateTime,y=gps$Northing,
#                                 + xout=tagDrag$dateTime+tOffset)$y
#tagDrag$eucDist2 <- sqrt(((tagDrag$Easting-tagDrag$gpsEasting2)^2)+
#                             + ((tagDrag$Northing-tagDrag$gpsNorthing2)^2))
#Output a list of variables for accessing Array using stationary tags
################################################################################
################################################################################


#For Stationary Test
#Fix Rate %
fixRateStat

#Average Accuracy (m)
avgAccStat <- mean(ref4$error)
avgAccStat

#Median Accuracy (m)
medianAccStat <- median(ref4$error)
medianAccStat

#Proportion with goal accuracy
propAccStat<-length(which( ref4$error < accGoal ))/
  + length(ref4$error)
propAccStat

#Number of points with greater than outlier size
outlierCountStat <- length(which(ref4$error > OutlierGoal))
outlierCountStat

#####################################################
###### DO NOT RUN THIS CODE ##################
#####################################################

#For tag drag Test
#Fix Rate
#fixRateMobile

#Average Accuracy
#meanAccMobile<- mean(na.omit(tagDrag$eucDist2))
#meanAccMobile

#Median Accuracy
#medianAccMobile <- median(na.omit(tagDrag$eucDist2))
#medianAccMobile

#Proportion with Goal accuracy
#propAccMobile <-length(which( tagDrag$eucDist2 < accGoal ))/
#  + length(tagDrag$eucDist2)
#propAccMobile

#Number of points with greater than outlier size (in our study 15 m)
#outlierCountMobile <- length(which( tagDrag$eucDist2 > OutlierGoal ))
#outlierCountMobile
################################################################################

#For Fixed Tag first
#What percentage of positions had and error below this value
percentGood1<- data.frame(ErrorRange)
percentGood1$goodperc<-NA

for( i in 1: length(percentGood1$ErrorRange)){percentGood1$goodperc [i] <- ((length(ref4$error[ref4$error < percentGood1$ErrorRange[i]])) / (length(ref4$error)))*100}

################################################################################
#tag drag Tag
#percentGood2<- data.frame(ErrorRange)
#percentGood2$goodperc<-NA
#for( i in 1: length(percentGood2$ErrorRange)){
#  + #What percentage of positions had error < value
#    + percentGood2$goodperc[i] <- ((length(tagDrag$eucDist2[tagDrag$eucDist2 <
#                                                              + percentGood2$ErrorRange[i]])) / (length(tagDrag$eucDist2)))*100
#    + }
################################################################################

#percent good
loc1 <- percentGood1
loc1$type<-c("Stationary Test 1")
#mobile1<- percentGood2
#mobile1$type<-c("Tag Drag Test")
AllperGood <- rbind(loc1)

plot0a4 <- ggplot(data=AllperGood, aes(x=AllperGood$ErrorRange,y=round(AllperGood$goodperc, digits=0)))+
  geom_point(size = 4)+
  geom_line()+
  theme_bw()

#plot0a4 <- ggplot(data=AllperGood, aes(x=ErrorRange,y=round(goodperc, digits=0))) +
#  geom_point(size= 4)+geom_line()+theme_bw() +
#  theme(element_text(size=25),
#        axis.text.y = element_text(angle=90, hjust= 0.5),
#        axis.title.y = element_text(vjust= -0.02 ),
#        axis.title.x = element_text(vjust= -0.02 ),
#        plot.margin = unit(c(1,1.5,1,1.5), "cm"),
#        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
#        legend.background = element_rect(fill = "transparent", colour = NA),
#        legend.position = c(0.6, 0.5),
#        legend.title=element_blank()) +
#  scale_y_continuous("Percent of all Positions",limits=c(0,100)) +
#  scale_x_continuous("Measured Error (m)",breaks= round(seq(min(
#    AllperGood$ErrorRange), max(AllperGood$ErrorRange), by = 1),1)) +
#  coord_cartesian( xlim=c(0,20.5))

#1. calculate the average error for each 1 unit HPE Bin.
#Bin by HPE
#create increments to bin by
breaks= seq(0, max(ref4$HPE),by=1)
# create a column with bins
ref4$bin <- cut(ref4$HPE,breaks)
#2. For each 1 m bin, the average HPE is calculated("binmean")
#use taplly to bin the error and calculate a min
binMean <- tapply(ref4$error,ref4$bin,mean)
#count the number of HPEs that fit into each bin value
binNum <- tapply(ref4$error,ref4$bin,length)
binNum[is.na(binNum)] <- 0
#New dataframe for bin data,cacluate xe and ye
bin<- data.frame(binMean, binNum)
#3. Every HPE value in a given 1 m bin has a corresponding HPEm. Each of...
#...these HPEm distances is composed of 2 elements: the error (difference...
#...between the calculated and measured position) in the X direction, and...
#...the error in the Y direction.

#I will refer to these error values as Xe and Ye, respectively.
#error in the x direction
ref4$xe <- sqrt((tag4Easting-(ref4East))^2)
#error in the y direction
ref4$ye <- sqrt((tag4Northing-(ref4North))^2)
#4. For each bin, the standard deviations of Xe and Ye are calculated.
bin$xeSd <- tapply(ref4$xe,ref4$bin,sd)
bin$yeSd <- tapply(ref4$ye,ref4$bin,sd)
#5. To convert the 2-dimensional standard deviations calculated in #4 into a...
#...single measure, the 2DRMS error is calculated from the standard...
#...deviations of Xe and Ye
bin$RMS2d <- 2*sqrt((bin$xeSd)^2 + (bin$yeSd)^2)
#6. Now create a line plot and a dataframe just for the...
#...numbers we need, (ie when we have at least 10 tag...
#...transmissions and an HPE less than 21)
bin$avgHPE <- tapply(ref4$HPE,ref4$bin,mean)
smallBin <- bin[ which(bin$binNum > 10),]
smallBin <- smallBin[which(smallBin$avgHPE < 25 ),]
# IRREGULAR CALCULATION OF 2DRMS so i filtered out the weird result
smallBin <- smallBin[which(smallBin$RMS2d < 9 ),]
res4 <- lm(smallBin$RMS2d ~ smallBin$avgHPE)

#PLOT THE POINT SPREAD
#Changing a bunch of theme elements to make it pretty.

xyplot<- data.frame(tag4Easting, tag4Northing)
plot1a4 <- qplot(tag4Easting, tag4Northing , alpha= 1/100,data = xyplot, xlab='Easting (m)', ylab= 'Northing (m) ') + theme_bw()
plot1a4 <- plot1a4 + geom_point(data=xyplot,aes(ref4East,ref4North),col="white", shape = 18, size= 4)
plot1a4 <- plot1a4 +theme(text = element_text(size=15),plot.title = element_text(vjust= 2), 
                          axis.text.y = element_text(angle=90, hjust= 0.5), 
                          axis.title.y = element_text(vjust= -0.02), 
                          axis.title.x = element_text(vjust= -0.02), 
                          plot.margin = unit(c(1,1.5,1,1.5), "cm"))
plot1a4 <- plot1a4+theme(legend.position="none")#remove legend
#plot1a4 <- plot1a4 + coord_cartesian(xlim=c(319500 ,319650),
#                                    ylim=c(5223500,5223660))
#Graph of HPE to error relationship
plot2a4 <- qplot(ref4$HPE, ref4$error , alpha= 1/100,
                 data = xyplot, xlab= 'HPE', ylab= 'Measured Error (m)')+
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_continuous(limits = c(5, 40)) + theme_bw() +
  theme(plot.margin=unit(c(1,1.5,1,1.5), "cm"),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+ theme(legend.position="none")
plot2a4 <- plot2a4 +theme(text = element_text(size=20),
                          plot.title = element_text(vjust= 2), axis.title.y = element_text(vjust= -0.02),
                          axis.title.x = element_text(vjust= -0.02))
plot2a4 <- plot2a4 + geom_point(data=smallBin, shape= 19,
                                size= 5,col='black', bg= 'black',alpha= 1, aes(x = avgHPE, y = RMS2d)) +
  theme(legend.position = "none")
plot2a4 <- plot2a4 + geom_point(data=smallBin, shape= 19, size= 4.2,col='white',
                                alpha= 1, aes(x = avgHPE, y = RMS2d)) + theme(legend.position = "none")
plot2a4 <- plot2a4 + geom_point(data=smallBin, shape= 4, size = 4,color='red',
                                aes(x = avgHPE, y = RMS2d)) +
  theme(legend.position = "none")

plot2a4 <- plot2a4 +geom_abline(data=res4, col='white',
                                aes(intercept=res4$coefficients[1], slope=res4$coefficients[2] ), size=1)
plot2a4 <- plot2a4 +geom_abline(data=res4, col='black', linetype='dashed',
                                aes(intercept=res4$coefficients[1], slope=res4$coefficients[2]), size= 1)
#Add model equation and r2 to graph
#function to add linear model to graph
#Setup for a linear model to be added to a ggplot
lm_eqn = function(m) {
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  if (coef(m)[2] >= 0) {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","
                     ~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","
                     ~~italic(r)^2~"="~r2,l)
  }
  as.character(as.expression(eq));
}
#calls in linear model function for linear model res4
plot2a4 = plot2a4 + annotate("text",x = 10, y = 20,
                             label = lm_eqn(res4), size= 4, parse=TRUE)

#What type of HPE cutoff should we use if we know what level of error we want
#!#! YOU MUST CHANGE accGoal to your desired error
#now evaluate how HPE cutoffs work if we want a specific error value
HPE<- c(5:HPEHigh)# we want a list of all possible HPE cutoffs
AllHPE <- data.frame(HPE)
AllHPE$ptsReject <- NA # number of points removed
AllHPE$ptsRejectP <- NA # percetage of all ponits
AllHPE$ptsRetain <- NA # number of poitns retained
AllHPE$ptsRetainP <- NA # percetage of all ponits
AllHPE$incorrectReject <- NA # good= error < accGoal and unacceptably...
#...erroneous= error > accGoal
AllHPE$incorrectRejectP <- NA
AllHPE$incorrectRetain<- NA
AllHPE$incorrectRetainP <- NA
AllHPE$incorrectRetainPvsRetain <-NA # incorretly retained/all retained *100
AllHPE$correctReject <- NA
AllHPE$correctRejectP <- NA

AllHPE$correctRetain <- NA
AllHPE$correctRetainP <- NA
AllHPE$goodDataLossP <- NA# incorrectReject/goodpts *100
AllHPE$badDataRetainP <- NA# incorrectRetain/AllHPE$ptsRetain * 100
AllHPE$sdDEV<-NA
AllHPE$avgErr<-NA
AllHPE$maxErr <- NA

# Number of all points that were unacceptably erroneous
badpts <-length(ref4$error[ref4$error > accGoal ])
# % of all points were unacceptably erroneous
badptsP <-(length(ref4$error[ref4$error > accGoal ])/(length(
  + ref4$error)))*100
#Number of all points that are good
goodpts <-length(ref4$error[ref4$error <= accGoal ])
#% of all points that are good
goodptsP <- (length(ref4$error[ref4$error <= accGoal ])/(length(ref4$error)))*100

for( i in 1: length(AllHPE$HPE)){AllHPE$ptsReject[i] <- length( which(ref4$HPE > AllHPE$HPE[i]))
AllHPE$ptsRejectP[i] <- length( which(ref4$HPE > AllHPE$HPE[i]))/(length(ref4$error))*100 
#% of all points dropped
AllHPE$ptsRetain[i] <- length( which(ref4$HPE < AllHPE$HPE[i]))
AllHPE$ptsRetainP[i] <- length( which(ref4$HPE < AllHPE$HPE[i]))/(length(ref4$error))*100 
#% of all points dropped
#incorrect reject
AllHPE$incorrectReject[i] <- length( which( (ref4$HPE[ ref4$error <= accGoal ] > AllHPE$HPE[i]) == TRUE))
#divided by total with acceptable error
AllHPE$incorrectRejectP[i] <-(AllHPE$incorrectReject[i]/length(ref4$error[ref4$error <= accGoal ]) )*100
#incorrect retain
AllHPE$incorrectRetain[i] <-length( which((ref4$HPE[ ref4$error > accGoal ] < AllHPE$HPE[i]) == TRUE))
#divided by accurate points
AllHPE$incorrectRetainP[i] <-(AllHPE$incorrectRetain[i]/length(ref4$error[ref4$error > accGoal ]))*100
#divided by total retained by filter
AllHPE$incorrectRetainPvsRetain[i] <-(AllHPE$incorrectRetain[i]/AllHPE$ptsRetain[i] )*100
#Correctly Rejected
AllHPE$correctReject[i] <-(length(which((ref4$HPE[ref4$error > accGoal] > AllHPE$HPE[i]) == TRUE)))
#divided by sum of all positions
AllHPE$correctRejectP[i]<- (AllHPE$correctReject[i]/AllHPE$ptsRetain[i])*100
#Correctly Retained
AllHPE$correctRetain[i] <- length(which((ref4$HPE[ ref4$error < accGoal] < AllHPE$HPE[i]) == TRUE))
#divided by sum of all positions
AllHPE$correctRetainP[i]<-(AllHPE$correctRetain[i]/length(ref4$error))*100
AllHPE$goodDataLossP[i] <- (AllHPE$incorrectReject[i] / goodpts)*100 
#incorrectReject/goodpts *100
AllHPE$badDataRetainP[i] <- (AllHPE$incorrectRetain[i]/AllHPE$ptsRetain[i])*100 
# incorrectRetain/AllHPE$ptsRetain * 100
AllHPE$avgErr[i] <- mean(ref4$error[ ref4$HPE < AllHPE$HPE[i]]) 
#mean of all error for each hpe
AllHPE$sdDEV [i]<- sd(ref4$error[ ref4$HPE < AllHPE$HPE[i]]) 
#Standard deviation of all error for each HPE
AllHPE$maxErr[i] <- max(ref4$error[ ref4$HPE < AllHPE$HPE[i]]) 
#mean of all error for each hpe
AllHPE$count15[i] <- length(which( (ref4$HPE[ ref4$error > 15] < AllHPE$HPE[i]) == TRUE ))
}
start <- HPELow
end <- HPEHigh
AllHPEperf1 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$correctRejectP[start:end])
AllHPEperf1$var <-"Correct Rejection"
AllHPEperf2 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$incorrectRetainP[start:end])
AllHPEperf2$var <-"Incorrect Retainment"
AllHPEperf3 <- data.frame(HPE=AllHPE$HPE[start:end],perc=AllHPE$incorrectRejectP[start:end])
AllHPEperf3$var <-"Incorrect Rejection"
#Combine all seconts for graph
AllHPEperf <- rbind(AllHPEperf1,AllHPEperf3,AllHPEperf2)

#note legend is alphabetical but graph is by dataframe order
plot3a4 <- qplot( AllHPE$incorrectRetainPvsRetain, AllHPE$incorrectRejectP,size=3, xlab= '% Incorrecly Retained (of all retained positions)',ylab= ' % Incorrectly Rejected (of all acceptable positions)' ) + 
  theme_bw() + theme( plot.margin = unit(c(1,1.5,1,1.5), "cm"),panel.grid.major=element_blank(), panel.grid.minor=element_blank())

plot3a4 <- plot3a4 +geom_line(data=AllHPE, aes(x=AllHPE$incorrectRetainPvsRetain, y=AllHPE$incorrectRejectP),color= 'black', size=0.25)
plot3a4 <- plot3a4 + theme_bw() + theme( text = element_text(size=15), #plot.title = element_text(vjust= 2),
                                         axis.text.y = element_text(angle=90, hjust= 0.4),
                                         axis.title.y = element_text(vjust= .3),
                                         axis.title.x = element_text(vjust= -0.02 ),
                                         plot.margin = unit(c(1,1.5,1,1.5), "cm"),
                                         panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                         legend.position = c(1.5, 0.8),
                                         legend.background = element_rect(fill = "transparent",colour = NA),
                                         legend.title=element_blank() )
plot3a4 <-plot3a4+ geom_text(cex=6,aes(x=(AllHPE$incorrectRetainPvsRetain[3:15]),
                                       y=(AllHPE$incorrectRejectP[3:15])+4 , label=AllHPE$HPE[3:15], group=NULL),)                      

plot4a4 <- ggplot( data= AllHPE[start:end,],
                   aes(x = AllHPE$HPE[start:end])) +
  geom_point(aes(y = AllHPE$avgErr[start:end],
                 shape = "Mean Error" ), size= 4, color= 'black') +
  geom_point(aes(y = AllHPE$maxErr[start:end],
                 shape = "Maximum Error" ), size= 4, color='black') +
  geom_line(aes(y = OutlierGoal, ), size= 1, lty= 'dashed') +
  theme_bw() + theme( text = element_text(size=25),
                      plot.title = element_text(vjust= 2),
                      axis.text.y = element_text(angle=90, hjust= 0.5),
                      axis.title.y = element_text(vjust= -0.02 ),
                      axis.title.x = element_text(vjust= -0.02 ),
                      plot.margin = unit(c(1,1.5,1,1.5), "cm"),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      legend.position = c(0.35, 0.75),
                      legend.background = element_rect(
                        fill = "transparent",colour = NA),
                      legend.title=element_blank() ) 
xlab('HPE') +
  coord_cartesian( ylim=c(0,20)) +
  scale_y_continuous("Error (m)",breaks= round(seq(min(0),
                                                   max(20), by = 2),1)) +
  scale_x_continuous("HPE",breaks=round(seq(min(AllHPE$HPE),
                                            max(AllHPE$HPE), by = 1),1)) +
  geom_text(cex=6,aes(x=AllHPE$HPE[1:21],y= (10),
                      label=AllHPE$count15[1:21], group=NULL),)

plot4a4

###############################################################################
##### DO NOT RUN THIS SECTION WE DO NOT HAVE MOVING TAG TEST DATA ############
################################################################################

#with(tagDrag, eucDist-eucDist2)
#badHPE<- tagDrag$HPE[tagDrag$HPE >= HPEfilter ]
#badeucDist4 <-tagDrag$eucDist2[tagDrag$HPE >= HPEfilter]
#badHPENorthing<- tagDrag$Northing[tagDrag$HPE >= HPEfilter ]
#badHPEEasting<- tagDrag$Easting[tagDrag$HPE >= HPEfilter ]
#badErrorNorthing<- tagDrag$Northing[tagDrag$eucDist2 >= accGoal ]
#badErrorEasting<- tagDrag$Easting[tagDrag$eucDist2 >= accGoal ]
#eucDist4<- tagDrag$eucDist2[tagDrag$eucDist2 >= accGoal ]
#incorrectlyRejected<- (tagDrag[tagDrag$eucDist2 <6 & tagDrag$HPE>8,])
#correctlyRejected<- na.omit(tagDrag[tagDrag$eucDist2 >6 & tagDrag$HPE>8,],)
#incorrectlyRetained<- na.omit(tagDrag[tagDrag$eucDist2 >6 & tagDrag$HPE<8,],)

#correctlyRetained <- na.omit(tagDrag[tagDrag$eucDist2 < 6 & tagDrag$HPE<8,],)
#badError<- data.frame("Easting"=badErrorEasting,"Northing"=badErrorNorthing,
#                      + "eucDist4" = eucDist4, "HPE" = tagDrag$HPE[tagDrag$eucDist2 >= accGoal])
#badError <- na.omit(badError)
#New frame here
#tag<-tagDrag[ is.na(tagDrag$gpsEasting)==FALSE,]


# identify unacceptably erroneous points in the plot(must change HPE)
#plot5a4 <- qplot(correctlyRetained$gpsEasting2, correctlyRetained$gpsNorthing2,
#                + xlab= ' Easting (m)', ylab= 'Northing (m) ') + theme_bw() + geom_point(data=
#                                                                                           + coast,aes(coast$Easting,coast$Northing),col="black", shape = 18, size= 5)
#plot5a4 <- plot5a4 +theme(text = element_text(size=20),
#                         + plot.title = element_text(vjust= 2),
#                         + axis.text.y = element_text(angle=90, hjust= 0.5),
#                         + axis.title.y = element_text(vjust=0.4 ),
#                         + axis.title.x = element_text(vjust= -0.01 ),
#                         + legend.position="none" ) #make the axis
#receivers
#plot5a4 <- plot5a4 + geom_point(data=rec2010, aes(x=Easting, y=Northing), shape= 3)
#all locations with a high HPE
#plot5a4 <- plot5a4 + geom_point(data= incorrectlyRejected,aes(x=Easting,
#                                                              + y=Northing), shape = 19, size= 5, col='Orange') # circled orange
#plot5a4 <- plot5a4 + geom_point(data= correctlyRejected,aes(x=Easting,
#                                                            + y=Northing), shape = 19, size= 5, col='Blue') # circled blue
#plot5a4 <- plot5a4 + geom_point(data= incorrectlyRetained,aes(x=Easting,
#                                                              + y=Northing), shape = 19, size= 5, col='red') # circled red
#plot5a4 <- plot5a4 + geom_point(data=correctlyRetained,aes(x=Easting,y=Northing),
#                                + shape = 19, size= 5, col='black') # circled black
#plot5a4 <- plot5a4 +coord_cartesian(xlim=c(xAxisMin,xAxisMax ),
#                                   + ylim=c(yAxisMin,yAxisMax))
#plot5b4 <- qplot(correctlyRetained$HPE, correctlyRetained$eucDist2 , alpha= 1/100,
#               + xlab= 'HPE', ylab= 'Euclidian distance (m)' )+ scale_y_continuous(limits =
#                                                                                     + c(0,21)) + scale_x_continuous(limits = c(0, 21)) + theme_bw() + theme(
#                                                                                       + plot.margin = unit(c(1,1.5,1,1.5), "cm"), panel.grid.major=element_blank(),
#                                                                                       + panel.grid.minor=element_blank())+ theme(legend.position="none") +
#  + scale_colour_identity()
#plot5b4 <- plot5b4 +theme(text = element_text(size=20),
#                       + plot.title = element_text(vjust= 2),
#                       + axis.title.y = element_text(vjust= -0.02 ),
#                       + axis.title.x = element_text(vjust= -0.02 ) ) #make the axis
#plot5b4 <- plot5b4 + geom_point(data= incorrectlyRejected,aes(x=HPE,y=eucDist2,
#                                                           + col="orange"), alpha=1,shape = 19, size= 4)
#plot5b4 <- plot5b4 + geom_point(data= correctlyRejected,aes(x=HPE, y=eucDist2,
#                                                         + col="blue"), shape = 19, size= 4)
#plot5b4 <- plot5b4 + geom_point(data= incorrectlyRetained,aes(x=HPE, y= eucDist2,
#                                                           + col="red"), shape = 19, size= 4)
################################################################################################

#Area of Receiver Coverage
tagll <- data.frame(X = rec$LONG, Y = rec$LAT)
attr(tagll, "projection") <- "LL"
xyrec <- convUL(tagll) * 1000
rec$EASTING<- xyrec[,1]
rec$NORTHING<- xyrec[,2]

#Converting fish location to UTM
tagll <- data.frame(X = Trial$LON, Y = Trial$LAT)
attr(tagll, "projection") <- "LL"
xyrec <- convUL(tagll) * 1000
Trial$EASTING<- xyrec[,1]
Trial$NORTHING<- xyrec[,2]

#Evaluate Actual fish trackin data (HPE FILTER)
plot6a4 <- qplot(Trial$EASTING, Trial$NORTHING, xlab= ' Easting (m)', ylab= 'Northing (m)', size=I(1)) + theme_bw()
plot6a4 <- plot6a4 + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                     ylim=c(yAxisMin,yAxisMax))
plot6a4 <- plot6a4 +theme(text = element_text(size=25),
                          plot.title = element_text(vjust= 2),
                          axis.text.y = element_text(angle=90, hjust= 0.5),
                          axis.title.y = element_text(vjust= -0.02 ),
                          axis.title.x = element_text(vjust= -0.02 ),
                          plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot6a4 <- plot6a4 + geom_point(data=Trial[Trial$HPE < 3,], aes(x=EASTING,
                                                                y=NORTHING),col="blue", shape = 20, size= I(1))
plot6a4 <- plot6a4 + geom_point(data=Trial[Trial$HPE >= 3 &
                                             Trial$HPE < 5 ,],aes(x=EASTING, y=NORTHING),col="#56B4E9",
                                shape = 20, size= I(1))
plot6a4 <- plot6a4 + geom_point(data=Trial[Trial$HPE >= 5 &
                                             Trial$HPE < 10 ,],aes(x=EASTING, y=NORTHING),col="green",
                                shape = 20, size= I(1))
plot6a4 <- plot6a4 + geom_point(data=Trial[Trial$HPE >= 10 &
                                             Trial$HPE < 15 ,],aes(x=EASTING, y=NORTHING),col="yellow",
                                shape = 20,size= I(1))
plot6a4 <- plot6a4 + geom_point(data=Trial[Trial$HPE >= 15 &
                                             Trial$HPE < 20 ,],aes(x=EASTING, y=NORTHING),col="#FFCC33",
                                shape = 20,size= I(1))
plot6a4 <- plot6a4 + geom_point(data=Trial[Trial$HPE >= 20 ,],aes(x=EASTING,
                                                                  y=NORTHING),col="red", shape = 20, size= I(1))
plot6a4 <- plot6a4 + geom_point(data=rec,aes(EASTING,NORTHING),fill="black",
                                col="white", shape = 24, size= 3, alpha=1)
plot6a4 <- plot6a4 + geom_point(data=bm_outline,aes(e,n),col="black",
                                shape = 18, size= 1)

xydata <- data.frame(X=Trial$EASTING, Y=Trial$NORTHING)

#create a polygon that covers the area of the data
fishPoly<- calcConvexHull(xydata, keepExtra=FALSE)
#Polygon(fishPoly)
#calcArea, calcCentroid are likley both useful
#Create a column in pos data.frame whethere in or outside array
Trial$array<- point.in.polygon(Trial$EASTING,Trial$NORTHING, rec$EASTING,
                               + rec$NORTHING, mode.checked= TRUE)
#plot of array polygon, with fish positions removed, 0 means not in array

#2 Filter data FOR FIRST SET OF PLOTS
cutoff <- HPEfilter
posDrop <- nrow(Trial[(Trial$HPE > cutoff),])
posDropArrayT <- ((Trial$HPE[Trial$array != 0]) > cutoff)
posDropArray <- length(subset(posDropArrayT, posDropArrayT == TRUE))
DropOutT <- ((Trial$HPE[Trial$array == 0]) > cutoff)
posDropOut <- length(subset( DropOutT, DropOutT == TRUE))
posTotal1<- nrow(Trial)
posTotal2<- length(Trial$array[Trial$array == 0])
posTotal3<- length(Trial$array[Trial$array != 0])
percDropFil1 <- (posDrop/posTotal1)*100 # % DroppedTota;
percADropFil1 <- (posDropArray/posTotal2)*100 #% dropped in the array
percODropFil1 <- (posDropOut/posTotal3)*100 #% dropped outside of the array

#Now identify the silhoutte of the point cluster using calcConvexHull
xydata1<-data.frame(X=Trial$EASTING[Trial$HPE < cutoff ], Y=Trial$NORTHING[Trial$HPE <
                                                                             + cutoff ])
#create a polygon that covers the area of the data
fishPoly1<- calcConvexHull(xydata1, keepExtra=FALSE)
head(xydata1)
points(fishPoly1$X, fishPoly1$Y, col='red')
plot(Trial$EASTING,Trial$NORTHING)
#calcArea, calcCentroid are likley both useful
areaFil1 <- calcArea(fishPoly1) #calculates area of polygon

#Plot 1- all fish positions
plot7a4 <- qplot(Trial$EASTING, Trial$NORTHING, xlab= 'Easting (m)', ylab=
                   'Northing (m)',alpha=1/100, size=I(1)) + theme_bw()
plot7a4 <- plot7a4 + coord_cartesian(xlim=c(xAxisMin,xAxisMax), ylim=
                                       c(yAxisMin,yAxisMax))
plot7a4 <- plot7a4 + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                                shape = 17, size= 2, alpha=1)
plot7a4 <- plot7a4 + theme(text = element_text(size=13),
                           plot.title = element_text(vjust= 2),
                           axis.text.y = element_text(angle=90, hjust= 0.5),
                           axis.title.y = element_text(vjust= -0.02 ),
                           axis.title.x = element_text(vjust= -0.02 ),
                           plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7a4 <- plot7a4 + theme(legend.position="none") #remove legend
#plot 2- all positions with HPE < cutoff
plot7b4 <- qplot( Trial$EASTING[Trial$HPE<= cutoff], Trial$NORTHING[Trial$HPE<= cutoff
], xlab= ' Easting (m)', ylab= 'Northing (m) ', alpha=1/100, size=I(1)) +
  theme_bw()
plot7b4 <- plot7b4 + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                     ylim=c(yAxisMin,yAxisMax))
plot7b4 <- plot7b4 + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                                shape = 17, size= 2, alpha=1)
plot7b4 <- plot7b4 +theme(text = element_text(size=13),
                          plot.title = element_text(vjust= 2),
                          axis.text.y = element_text(angle=90, hjust= 0.5),
                          axis.title.y = element_text(vjust= -0.02 ),
                          axis.title.x = element_text(vjust= -0.02 ),
                          plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7b4 <- plot7b4 + theme(legend.position="none") #remove legend
#plot 3- all positions with HPE > cutoff)
plot7c4 <- qplot( Trial$EASTING[Trial$HPE> cutoff], Trial$NORTHING[Trial$HPE> cutoff],
                  xlab= ' Easting (m)', ylab= 'Northing (m) ',alpha=1/100, size=I(1)) +
  theme_bw()
plot7c4 <- plot7c4 + coord_cartesian(xlim=c(xAxisMin,xAxisMax),
                                     ylim=c(yAxisMin,yAxisMax))
plot7c4 <- plot7c4 + geom_point(data=rec,aes(EASTING, NORTHING),col="red",
                                shape = 17, size= 2, alpha=1)
plot7c4 <- plot7c4 +theme(text = element_text(size=13),
                          plot.title = element_text(vjust= 2),
                          axis.text.y = element_text(angle=90, hjust= 0.5),
                          axis.title.y = element_text(vjust= -0.02 ),
                          axis.title.x = element_text(vjust= -0.02 ),
                          plot.margin = unit(c(1,1.5,1,1.5), "cm") ) #make the axis
plot7c4 <- plot7c4 + theme(legend.position="none") #remove legend




########################################################
####### GROUP TOGETHER PLOTS THAT WILL BE MERGED #######
########################################################
#PLOT THE POINT SPREAD
#Changing a bunch of theme elements to make it pretty.

xyplot<- data.frame(tag2Easting, tag2Northing)
plot1a2<- qplot(tag2Easting, tag2Northing , alpha= 1/100,data = xyplot, xlab='Easting (m)', ylab= 'Northing (m) ') + theme_bw()
plot1a2<- plot1a2 + geom_point(data=xyplot,aes(ref2East,ref2North),col="white", shape = 18, size=5)
plot1a2<- plot1a2 + theme_classic()+
  scale_y_continuous(breaks = c(5222440,5222470,5222500))+
  scale_x_continuous(breaks = c(318350,318400,318450))+
  theme(axis.text.x = element_text(color = "grey20", size = 19, angle = 0, hjust = .5, vjust = 1, face = "plain",family = "serif"),
        axis.text.y = element_text(color = "grey20", size = 19, angle = 90, hjust = .5, vjust = 0, face = "plain",family = "serif"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 1, face = "plain",family = "serif"),
        legend.position = c(0.12,0.98), legend.title = element_text(),legend.text = element_text(size = 20, family = "serif"),
        plot.title = element_text(face = "bold", size = 20, family = "serif", hjust = .5),
        plot.margin = margin(.8,.5,.3,.3, "cm"),
        plot.background = element_rect(colour = "black", fill=NA, size=1))
plot1a2<- plot1a2+theme(legend.position="none")#remove legend

plot1a2<- plot1a2 + coord_cartesian(xlim=c(318350 ,318450),
                                    ylim=c(5222425,5222500))
#plot1a2
#####################################################################

#Graph of HPE to error relationship
plot2a2<- qplot(ref2$HPE, ref2$error , alpha= 1/100,
                data = xyplot, xlab= 'HPE', ylab= 'Measured Error (m)')+
  scale_y_continuous(limits = c(0, 32), breaks = seq(0,32, by = 6)) +
  scale_x_continuous(limits = c(9, 30)) + 
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 19, angle = 0, hjust = .5, vjust = 1, face = "plain",family = "serif"),
        axis.text.y = element_text(color = "grey20", size = 19, angle = 90, hjust = .5, vjust = 0, face = "plain",family = "serif"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 1, face = "plain",family = "serif"),
        legend.position = c(0.12,0.98), legend.title = element_text(),legend.text = element_text(size = 20, family = "serif"),
        plot.title = element_text(face = "bold", size = 19, family = "serif", hjust = .5),
        plot.margin = margin(.2,.2,.2,.2, "cm"),
        plot.background = element_rect(colour = "black", fill=NA, size=1))
plot2a2<- plot2a2 +theme(text = element_text(size=15),
                         plot.title = element_text(vjust= 2), axis.title.y = element_text(vjust= -0.02),
                         axis.title.x = element_text(vjust= -0.02))
plot2a2 <- plot2a2 + geom_point(data=smallBin, shape= 19,
                                size= 4,col='black', bg= 'black',alpha= 1, aes(x = avgHPE, y = RMS2d)) +
  theme(legend.position = "none")

plot2a2<- plot2a2 +geom_abline(data=res3, col='white',
                               aes(intercept=res3$coefficients[1], slope=res3$coefficients[2] ), size=2)
plot2a2<- plot2a2 +geom_abline(data=res3, col='red', linetype='dashed',
                               aes(intercept=res3$coefficients[1], slope=res3$coefficients[2]), size= 2)

plot2a2<- plot2a2 + geom_point(data=smallBin, shape= 19, size= 5,col='red',
                               alpha= 1, aes(x = avgHPE, y = RMS2d)) + theme(legend.position = "none")
#plot2a2<- plot2a2 + geom_point(data=smallBin, shape= 4, size = 5,color='white',
 #                              aes(x = avgHPE, y = RMS2d)) +
#  theme(legend.position = "none")
#plot2a2

#Add model equation and r2 to graph
#function to add linear model to graph
#Setup for a linear model to be added to a ggplot
lm_eqn = function(m) {
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  if (coef(m)[2] >= 0) {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","
                     ~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","
                     ~~italic(r)^2~"="~r2,l)
  }
  as.character(as.expression(eq));
}

#calls in linear model function for linear model res3
plot2a2=plot2a2 + annotate("text",x = 15, y = 31.5,
                           label = lm_eqn(res3), size= 6, parse=TRUE,family = "serif")
#plot2a2

#################################################################################
library(ggplot2)
library(ggpubr)
#http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/
figure1 <- ggarrange(plot1a2, plot2a2,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure1

citation()
