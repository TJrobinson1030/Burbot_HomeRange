################################################################################
############# HOME RANGE DATA ANALYSIS AND VISUALIZATION #######################
################################################################################
rm(list=ls(all=TRUE))

library(readr)
library(anytime)
library(tidyverse)
library(rstatix)

# Merge all depth data to one file - file includes only unique detections #
Burb_Home_Range = list.files(path = "D:/BSUGradBurbot/DATA/FishData/FinalFilteredData/HomeRange_AllSeasons/LakeTrim",
                             pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                                                      # Store all files in list
  bind_rows                                                                 # Combine data sets into one data set 
unique(Burb_Home_Range$DATE)

Burb_Home_Range$DATE = anytime(as.factor(Burb_Home_Range$DATE))
Burb_Home_Range$DATE = as.Date(Burb_Home_Range$DATE)

#stats summary by sex
###### IMPORTANT TABLE #######
###### TABLE 1.3 ######
##############################
stats_90_sex_all <- Burb_Home_Range %>% 
  group_by(DATE,SEX) %>%
  get_summary_stats(HR90, type = "full")

#palette for colorblind
cbbPaletteG <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#convert to posixCT then to date for graphing purposes
stats_90_sex_all$DATE = anytime(as.factor(stats_90_sex_all$DATE))
stats_90_sex_all$DATE = as.Date(stats_90_sex_all$DATE)

str(stats_90_sex_all)

#interpolated line by group
#MALE
#spline for line
male <- stats_90_sex_all[stats_90_sex_all$SEX == "M",]
spline.d.m <- as.data.frame(spline(male$DATE,male$mean))
spline.d.m$x <- as.Date(spline.d.m$x, origin = "1970-01-01")
spline.d.m$SEX <- "M"

#spline for polygon (error)
hr.ci.m <- male %>% 
  dplyr :: select(DATE,mean,ci)
hr.ci.m <- hr.ci.m %>% mutate(ymin = mean-ci,
                              ymax = mean+ci)
spline.lci.m <- as.data.frame(spline(hr.ci.m$DATE,hr.ci.m$ymin))
spline.uci.m <- as.data.frame(spline(hr.ci.m$DATE,hr.ci.m$ymax))
spline.ci.m <- merge(spline.lci.m,spline.uci.m, by = "x")
spline.ci.m$x <- as.Date(spline.ci.m$x, origin = "1970-01-01")
colnames(spline.ci.m) <- c("DATE","ci.l","ci.u")
spline.ci.m$SEX <- "M"

#FEMALE
#spline for line
female <- stats_90_sex_all[stats_90_sex_all$SEX == "F",]
spline.d.f <- as.data.frame(spline(female$DATE,female$mean))
spline.d.f$x <- as.Date(spline.d.f$x, origin = "1970-01-01")
spline.d.f$SEX <- "F"

#spline for polygon (error)
hr.ci.f <- female %>% 
  dplyr :: select(DATE,mean,ci)
hr.ci.f <- hr.ci.f %>% mutate(ymin = mean-ci,
                              ymax = mean+ci)
spline.lci.f <- as.data.frame(spline(hr.ci.f$DATE,hr.ci.f$ymin))
spline.uci.f <- as.data.frame(spline(hr.ci.f$DATE,hr.ci.f$ymax))
spline.ci.f <- merge(spline.lci.f,spline.uci.f, by = "x")
spline.ci.f$x <- as.Date(spline.ci.f$x, origin = "1970-01-01")
colnames(spline.ci.f) <- c("DATE","ci.l","ci.u")
spline.ci.f$SEX <- "F"

#UNKNOWN
#spline for line
unknown <- stats_90_sex_all[stats_90_sex_all$SEX == "U",]
spline.d.u <- as.data.frame(spline(unknown$DATE,unknown$mean))
spline.d.u$x <- as.Date(spline.d.u$x, origin = "1970-01-01")
spline.d.u$SEX <- "U"

#spline for polygon (error)
hr.ci.u <- unknown %>% 
  dplyr :: select(DATE,mean,ci)
hr.ci.u <- hr.ci.u %>% mutate(ymin = mean-ci,
                              ymax = mean+ci)
spline.lci.u <- as.data.frame(spline(hr.ci.u$DATE,hr.ci.u$ymin))
spline.uci.u <- as.data.frame(spline(hr.ci.u$DATE,hr.ci.u$ymax))
spline.ci.u <- merge(spline.lci.u,spline.uci.u, by = "x")
spline.ci.u$x <- as.Date(spline.ci.u$x, origin = "1970-01-01")
colnames(spline.ci.u) <- c("DATE","ci.l","ci.u")
spline.ci.u$SEX <- "U"

#merge data frames
spline.d <- rbind(spline.d.m,spline.d.f,spline.d.u)
str(spline.d)
colnames(spline.d) <- c("DATE","MEAN","SEX")
spline.d <- spline.d %>% arrange(DATE)

spline.ci <- rbind(spline.ci.m,spline.ci.f,spline.ci.u)
str(spline.ci)
spline.ci <- spline.ci %>% arrange(DATE)

spline <- cbind(spline.d,spline.ci) %>% 
  dplyr :: select(1,2,5,6,7)

#create labels for dates in plot
break.vec <- c(as.Date("2019-04-10"),
               seq(from = as.Date("2019-04-13"), to = as.Date("2020-05-27"),
                   by = "34 days"),
               as.Date("2020-05-29"))

#theme for 4:3 slide in pres
theme_P <-theme(plot.title = element_text(face = "bold", size = 18, family = "serif"),
      axis.text.x = element_text(color = "grey20", size = 18, angle = -45, hjust = .5, vjust = 0, face = "plain",family = "serif"),
      axis.text.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = 0.9, face = "plain",family = "serif"),  
      axis.title.x = element_text(color = "grey20", size = 19, angle = 0, hjust = .5, vjust = 0, face = "plain",family = "serif"),
      axis.title.y = element_text(color = "grey20", size = 19, angle = 90, hjust = .5, vjust = 0, face = "plain",family = "serif"),
      legend.position = c(0.15,0.91), legend.title = element_text(),legend.text = element_text(size = 19, family = "serif"),
      plot.margin = unit(c(0.3,1.5,-0.6,.1),"cm"),
      plot.background = element_rect(colour = "black", fill="white", size=1))

#figure for moving average home range by sex
HR_25day <- ggplot(data = spline, aes(x = DATE, y = MEAN, colour = SEX, group = SEX))+
  geom_rect(xmin = as.Date("2019-12-03"), xmax = as.Date("2020-04-29"), ymin = -0.2,ymax = 2.09,alpha=.02, colour = "grey",fill = "grey")+
  geom_rect(xmin = as.Date("2019-04-13"), xmax = as.Date("2019-04-27"), ymin = -0.2,ymax = 2.09,alpha=.02, colour = "grey",fill = "grey")+
  geom_line(aes(x = DATE, y = MEAN, colour = SEX))+
  geom_ribbon(aes(x = DATE, ymin=ci.l, ymax=ci.u, fill = SEX),alpha = 0.1)+
  #for color coding to match with ribbon
  scale_colour_manual(values=cbbPaletteG, 
                      name = "", 
                      breaks=c("F","M","U"),
                      labels=c("Female", "Male","Unknown"))+
  #for color coding to match with data points
  scale_fill_manual(values=cbbPaletteG, 
                    name = "", 
                    breaks=c("F","M","U"),
                    labels=c("Female", "Male","Unknown"))+
  geom_point(data = stats_90_sex_all, aes(x = DATE, y = mean, colour = SEX), size = 2)+
  scale_y_continuous(limits = c(-0.2,2.1))+
  scale_x_date(date_labels = "%d/%b/%y", breaks = break.vec)+
  labs(title="Estimated Home Range Size of Burbot",y=expression ("90% Kernel Density Range  " (km^2)), fill = "", x = "")+
  theme_classic()+
  theme(plot.title = element_text(face = "bold", size = 16, family = "serif"),
        axis.text.x = element_text(color = "grey20", size = 15, angle = -45, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = 0.9, face = "plain",family = "serif"),  
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        legend.position = c(0.12,0.9), legend.title = element_text(),legend.text = element_text(size = 15, family = "serif"),
        plot.margin = unit(c(0.3,1.5,-0.3,.1),"cm"),
        plot.background = element_rect(colour = "black", fill="white", size=1))+
  coord_cartesian(expand = FALSE)

#https://r-graphics.org/recipe-annotate-rect
HR_25day
HR_25day_small <- HR_25day

#save figure for pres
setwd("E:/BSUGradBurbot/Figures_Plots/Thesis_Pres_Fig")
#ggsave("HR_25day.PNG", width = 12, height = 9, units = "in",dpi= 300)
#ggsave("HR_25day_small.PNG", width = 5, height = 4, units = "in",dpi= 300)

################################################################################
######################### STATS WITH SEASON ############################
################################################################################

Ice19 <- Burb_Home_Range[Burb_Home_Range$DATE >= "2019-04-13" & Burb_Home_Range$DATE <= "2019-04-27",]
Ice20 <- Burb_Home_Range[Burb_Home_Range$DATE >= "2019-12-03" & Burb_Home_Range$DATE <= "2020-04-29",]
Open19 <- Burb_Home_Range[Burb_Home_Range$DATE >= "2019-04-27" & Burb_Home_Range$DATE <= "2019-12-03",]
Open20 <- Burb_Home_Range[Burb_Home_Range$DATE >= "2020-04-29" & Burb_Home_Range$DATE <= "2020-07-29",]

Ice <-rbind(Ice19,Ice20)
Open <-rbind(Open19,Open20)
Ice$SEASON <- "Ice"
Open$SEASON <- "Open"

Burb_Home_Range <- rbind(Ice,Open)


#stats summary by SEASON
library(rstatix)
stats_90_season_all <- Burb_Home_Range %>% 
  group_by(DATE,SEASON) %>%
  get_summary_stats(HR90, type = "full")

#interpolated line by group
#spline for line
spline.l.season <- as.data.frame(spline(stats_90_season_all$DATE,stats_90_season_all$mean))
spline.l.season$x <- as.Date(spline.l.season$x, origin = "1970-01-01")
colnames(spline.l.season) <- c("DATE","MEAN")

#spline for polygon (error)
hr.ci.seas <- stats_90_season_all %>% 
  dplyr ::  select(DATE,mean,ci)
hr.ci.seas <- hr.ci.seas %>% mutate(ymin = mean-ci,
                                    ymax = mean+ci)
spline.lci.seas <- as.data.frame(spline(hr.ci.seas$DATE,hr.ci.seas$ymin))
spline.uci.seas <- as.data.frame(spline(hr.ci.seas$DATE,hr.ci.seas$ymax))
spline.ci.season <- merge(spline.lci.seas,spline.uci.seas, by = "x")
spline.ci.season$x <- as.Date(spline.ci.season$x, origin = "1970-01-01")
colnames(spline.ci.season) <- c("DATE","ci.l","ci.u")

#cbind data frames
spline.seas <- cbind(spline.l.season,spline.ci.season) %>% 
  dplyr ::   select(1,2,4,5)

#create labels for dates in plot
break.vec <- c(as.Date("2019-04-10"),
               seq(from = as.Date("2019-04-13"), to = as.Date("2020-05-25"),
                   by = "34 days"),
               as.Date("2020-05-29"))

#palette for colorblind
cbbPaletteG <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#figure for moving average home range by sex
HR_25day_season <- ggplot(data = spline.seas, aes(x = DATE, y = MEAN))+
  geom_line(aes(x = DATE, y = MEAN))+
  geom_ribbon(aes(x = DATE, ymin=ci.l, ymax=ci.u),alpha = 0.1)+
  #for color coding to match with ribbon
  scale_colour_manual(values=cbbPaletteG, 
                      name = "", 
                      breaks=c("Ice","Open"),
                      labels=c("Ice", "Open"))+
  #for color coding to match with data points
  scale_fill_manual(values=cbbPaletteG, 
                    name = "", 
                    breaks=c("Ice","Open"),
                    labels=c("Ice", "Open"))+
  geom_point(data = stats_90_season_all, aes(x = DATE, y = mean, colour = SEASON), size = 3)+
  scale_y_continuous(limits = c(0.0000001,2))+
  scale_x_date(date_labels = "%d/%b/%y", breaks = break.vec)+
  labs(y=expression ("90% Kernel Density Range " (Km^2)), size = 8, family="serif", fill = "", x = "")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = -45, hjust = .5, vjust = .1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = 0.5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = 0, face = "plain"),
        legend.position = c(0.21,0.88), legend.title = element_text(),legend.text = element_text(size = 18, family = "serif"),
        plot.title = element_text(face = "bold", size = 15, family = "serif"),
        plot.margin = unit(c(0.6,1.2,-0.6,.1),"cm"),
        plot.background = element_rect(colour = "black", fill=NA, size=1))+
  coord_cartesian(expand = FALSE)

HR_25day_season


######## 
# EXTRA TABLES FOR STATS
stats_table_season <- Burb_Home_Range %>% 
  group_by(SEASON) %>% 
  get_summary_stats(HR90,type = "full")

stats_table_season_sex <- Burb_Home_Range %>% 
  group_by(SEASON,SEX) %>% 
  get_summary_stats(HR90,type = "full")

stats_table_date <- Burb_Home_Range %>% 
  group_by(DATE) %>% 
  get_summary_stats(HR90,type = "full")

################################################################################
########### CORE RANGE BURBOT ##################################################
################################################################################

#stats summary by sex
library(rstatix)
stats_50_sex_all <- Burb_Home_Range %>% 
  group_by(DATE,SEX) %>%
  get_summary_stats(HR50, type = "full")

#palette for colorblind
cbbPaletteG <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#convert to posixCT then to date for graphing purposes
stats_50_sex_all$DATE = anytime(as.factor(stats_50_sex_all$DATE))
stats_50_sex_all$DATE = as.Date(stats_50_sex_all$DATE)

str(stats_50_sex_all)

#interpolated line by group
#MALE
#spline for line
male.cr <- stats_50_sex_all[stats_50_sex_all$SEX == "M",]
spline.cr.m <- as.data.frame(spline(male.cr$DATE,male.cr$mean))
spline.cr.m$x <- as.Date(spline.cr.m$x, origin = "1970-01-01")
spline.cr.m$SEX <- "M"

#spline for polygon (error)
cr.ci.m <- male.cr %>% 
  dplyr ::   select(DATE,mean,ci)
cr.ci.m <- cr.ci.m %>% mutate(ymin = mean-ci,
                              ymax = mean+ci)
spline.lci.m.cr <- as.data.frame(spline(cr.ci.m$DATE,cr.ci.m$ymin))
spline.uci.m.cr <- as.data.frame(spline(cr.ci.m$DATE,cr.ci.m$ymax))
spline.ci.m.cr <- merge(spline.lci.m.cr,spline.uci.m.cr, by = "x")
spline.ci.m.cr$x <- as.Date(spline.ci.m.cr$x, origin = "1970-01-01")
colnames(spline.ci.m.cr) <- c("DATE","ci.l","ci.u")
spline.ci.m.cr$SEX <- "M"

#FEMALE
#spline for line
female.cr <- stats_50_sex_all[stats_50_sex_all$SEX == "F",]
spline.cr.f <- as.data.frame(spline(female.cr$DATE,female.cr$mean))
spline.cr.f$x <- as.Date(spline.cr.f$x, origin = "1970-01-01")
spline.cr.f$SEX <- "F"

#spline for polygon (error)
cr.ci.f <- female.cr %>% 
  dplyr ::   select(DATE,mean,ci)
cr.ci.f <- cr.ci.f %>% mutate(ymin = mean-ci,
                              ymax = mean+ci)
spline.lci.f.cr <- as.data.frame(spline(cr.ci.f$DATE,cr.ci.f$ymin))
spline.uci.f.cr <- as.data.frame(spline(cr.ci.f$DATE,cr.ci.f$ymax))
spline.ci.f.cr <- merge(spline.lci.f.cr,spline.uci.f.cr, by = "x")
spline.ci.f.cr$x <- as.Date(spline.ci.f.cr$x, origin = "1970-01-01")
colnames(spline.ci.f.cr) <- c("DATE","ci.l","ci.u")
spline.ci.f.cr$SEX <- "F"

#UNKNOWN
#spline for line
unknown.cr <- stats_50_sex_all[stats_50_sex_all$SEX == "U",]
spline.cr.u <- as.data.frame(spline(unknown.cr$DATE,unknown.cr$mean))
spline.cr.u$x <- as.Date(spline.cr.u$x, origin = "1970-01-01")
spline.cr.u$SEX <- "U"

#spline for polygon (error)
cr.ci.u <- unknown.cr %>% 
  dplyr ::   select(DATE,mean,ci)
cr.ci.u <- cr.ci.u %>% mutate(ymin = mean-ci,
                              ymax = mean+ci)
spline.lci.u.cr <- as.data.frame(spline(cr.ci.u$DATE,cr.ci.u$ymin))
spline.uci.u.cr <- as.data.frame(spline(cr.ci.u$DATE,cr.ci.u$ymax))
spline.ci.u.cr <- merge(spline.lci.u.cr,spline.uci.u.cr, by = "x")
spline.ci.u.cr$x <- as.Date(spline.ci.u.cr$x, origin = "1970-01-01")
colnames(spline.ci.u.cr) <- c("DATE","ci.l","ci.u")
spline.ci.u.cr$SEX <- "U"

#merge data frames
spline.cr <- rbind(spline.cr.m,spline.cr.f,spline.cr.u)
str(spline.cr)
colnames(spline.cr) <- c("DATE","MEAN","SEX")
spline.cr <- spline.cr %>% arrange(DATE)

spline.ci.cr <- rbind(spline.ci.m.cr,spline.ci.f.cr,spline.ci.u.cr)
str(spline.ci.cr)
spline.ci.cr <- spline.ci.cr %>% arrange(DATE)

spline.cr <- cbind(spline.cr,spline.ci.cr) %>% 
  dplyr ::   select(1,2,5,6,7)

#create labels for dates in plot
break.vec <- c(as.Date("2019-04-10"),
               seq(from = as.Date("2019-04-13"), to = as.Date("2020-05-27"),
                   by = "34 days"),
               as.Date("2020-05-29"))
#geom ribbon
#https://stackoverflow.com/questions/54648389/ggplotly-with-geom-ribbon-grouping

#figure for moving average home range by sex
CR_25day <- ggplot(data = spline.cr, aes(x = DATE, y = MEAN, colour = SEX, group = SEX))+
  geom_rect(xmin = as.Date("2019-12-03"), xmax = as.Date("2020-04-29"), ymin = -0.15,ymax = 1,alpha=.02, colour = "grey",fill = "grey")+
  geom_rect(xmin = as.Date("2019-04-13"), xmax = as.Date("2019-04-27"), ymin = -0.15,ymax = 1,alpha=.02, colour = "grey",fill = "grey")+
  geom_line(aes(x = DATE, y = MEAN, colour = SEX))+
  geom_ribbon(aes(x = DATE, ymin=ci.l, ymax=ci.u, fill = SEX),alpha = 0.1)+
  #for color coding to match with ribbon
  scale_colour_manual(values=cbbPaletteG, 
                      name = "", 
                      breaks=c("F","M","U"),
                      labels=c("Female", "Male","Unknown"))+
  #for color coding to match with data points
  scale_fill_manual(values=cbbPaletteG, 
                    name = "", 
                    breaks=c("F","M","U"),
                    labels=c("Female", "Male","Unknown"))+
  geom_point(data = stats_50_sex_all, aes(x = DATE, y = mean, colour = SEX), size = 2)+
  scale_y_continuous(limits = c(-0.15,1))+
  scale_x_date(date_labels = "%d/%b/%y", breaks = break.vec)+
  labs(title="Estimated Core Range Size of Burbot",y=expression ("50% Kernel Density Range  " (km^2)),fill = "", x = "")+
  theme_classic()+
  #theme size is good for 4:3 slide on powerpoint at figure full size of slide
  theme(plot.title = element_text(face = "bold", size = 16, family = "serif"),
        axis.text.x = element_text(color = "grey20", size = 15, angle = -45, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = 0.9, face = "plain",family = "serif"),  
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = 0, face = "plain",family = "serif"),
        legend.position = c(0.12,0.9), legend.title = element_text(),legend.text = element_text(size = 15, family = "serif"),
        plot.margin = unit(c(0.3,1.5,-0.3,.1),"cm"),
        plot.background = element_rect(colour = "black", fill="white", size=1))+
  coord_cartesian(expand = FALSE)

#plot.margin(top,right,bottom,left)
#https://stackoverflow.com/questions/18252827/increasing-area-around-plot-area-in-ggplot2

CR_25day
CR_25day_small <- CR_25day

#save figure for pres
setwd("E:/BSUGradBurbot/Figures_Plots/Thesis_Pres_Fig")
#ggsave("CR_25day.PNG", width = 12, height = 9, units = "in",dpi= 300)
#ggsave("CR_25day_small.PNG", width = 6, height = 4.5, units = "in",dpi= 300)
#
################################################################################

###### 
# EXTRA TABLES
#stats summary by overall season
library(rstatix)
stats_50_seas_sex_all <- Burb_Home_Range %>% 
  group_by(SEASON,SEX) %>%
  get_summary_stats(HR50, type = "full")

################################################################################
######################### STATS WITH SEASON ############################
################################################################################

setwd("E:/BSUGradBurbot/DATA/FishData/FinalFilteredData/HomeRange_AllSeasons")
#write.csv(Burb_Home_Range,"HR_Per_Fish.csv")

#stats summary by SEASON
library(rstatix)
stats_50_season_all <- Burb_Home_Range %>% 
  group_by(DATE,SEASON) %>%
  get_summary_stats(HR50, type = "full")


#convert to posixCT then to date for graphing purposes
stats_50_season_all$DATE = anytime(as.factor(stats_50_season_all$DATE))
stats_50_season_all$DATE = as.Date(stats_50_season_all$DATE)

str(stats_50_season_all)

#interpolated line by group
#spline for line
spline.l.season.cr <- as.data.frame(spline(stats_50_season_all$DATE,stats_50_season_all$mean))
spline.l.season.cr$x <- as.Date(spline.l.season.cr$x, origin = "1970-01-01")
colnames(spline.l.season.cr) <- c("DATE","MEAN")

#spline for polygon (error)
cr.ci.seas <- stats_50_season_all %>% 
  dplyr ::   select(DATE,mean,ci)
cr.ci.seas <- cr.ci.seas %>% mutate(ymin = mean-ci,
                                    ymax = mean+ci)
spline.lci.seas.cr <- as.data.frame(spline(cr.ci.seas$DATE,cr.ci.seas$ymin))
spline.uci.seas.cr <- as.data.frame(spline(cr.ci.seas$DATE,cr.ci.seas$ymax))
spline.ci.season.cr <- merge(spline.lci.seas.cr,spline.uci.seas.cr, by = "x")
spline.ci.season.cr$x <- as.Date(spline.ci.season.cr$x, origin = "1970-01-01")
colnames(spline.ci.season.cr) <- c("DATE","ci.l","ci.u")

#column bind
spline.seas.cr <- cbind(spline.l.season.cr,spline.ci.season.cr) %>% 
  dplyr ::   select(1,2,4,5)

#create labels for dates in plot
#https://stackoverflow.com/questions/14759676/specification-of-first-and-last-tick-marks-with-scale-x-date
break.vec <- c(as.Date("2019-04-10"),
               seq(from = as.Date("2019-04-13"), to = as.Date("2020-05-25"),
                   by = "34 days"),
               as.Date("2020-05-29"))

#palette for colorblind
cbbPaletteG <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#figure for moving average home range by sex
CR_25day_season <- ggplot(data = spline.seas.cr, aes(x = DATE, y = MEAN))+
  geom_line(aes(x = DATE, y = MEAN))+
  geom_ribbon(aes(x = DATE, ymin=ci.l, ymax=ci.u),alpha = 0.1)+
  #for color coding to match with ribbon
  scale_colour_manual(values=cbbPaletteG, 
                      name = "", 
                      breaks=c("Ice","Open"),
                      labels=c("Ice", "Open"))+
  #for color coding to match with data points
  scale_fill_manual(values=cbbPaletteG, 
                    name = "", 
                    breaks=c("Ice","Open"),
                    labels=c("Ice", "Open"))+
  geom_point(data = stats_50_season_all, aes(x = DATE, y = mean, colour = SEASON), size = 3)+
  scale_y_continuous(limits = c(0.0000001,.601))+
  scale_x_date(date_labels = "%d/%b/%y", breaks = break.vec)+
  labs(y=expression ("50% Kernel Density Range " (Km^2)), size = 8, family="serif", fill = "", x = "")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = -45, hjust = .5, vjust = .1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = 0.5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = 0, face = "plain"),
        legend.position = c(0.09,0.88), legend.title = element_text(),legend.text = element_text(size = 18, family = "serif"),
        plot.title = element_text(face = "bold", size = 15, family = "serif"),
        plot.margin = unit(c(0.6,1.2,-0.6,.1),"cm"),
        plot.background = element_rect(colour = "black", fill=NA, size=1))+
  coord_cartesian(expand = FALSE)

CR_25day_season



# EXTRACT LARGEST HR FOR EACH FISH
# dplyr
require(dplyr)
    #https://stackoverflow.com/questions/24237399/how-to-select-the-rows-with-maximum-values-in-each-group-with-dplyr
large_HR <- Burb_Home_Range %>% group_by(TRANSMITTER) %>% filter(HR90 == max(HR90))


################################################################################
################## GENERALIZED LINEAR MIXED MODELS #############################
################################################################################

############# HOME RANGE ####################
stats_90_table_season <- Burb_Home_Range %>% 
  group_by(SEASON) %>% 
  get_summary_stats(HR90,type = "full")

stats_90_table_season_trans <- Burb_Home_Range %>% 
  group_by(SEASON, TRANSMITTER) %>% 
  get_summary_stats(HR90,type = "full")
setwd("C:/Users/tjrobinson/Desktop")
write.csv(stats_90_table_season_trans,"HR_Per_Fish_Season.csv")
ggplot(data = Burb_Home_Range, aes(x=HR90, fill = SEASON, color = SEASON))+
  geom_histogram(alpha=0.5,position = "identity")+
  facet_grid(SEASON~.)

stats_90_table_season_sex <- Burb_Home_Range %>% 
  group_by(SEASON,SEX) %>% 
  get_summary_stats(HR90,type = "full")

stats_50_table_season_sex <- Burb_Home_Range %>% 
  group_by(SEASON,SEX) %>% 
  get_summary_stats(HR50,type = "full")

stats_90_table_date <- Burb_Home_Range %>% 
  group_by(DATE) %>% 
  get_summary_stats(HR90,type = "full")

stats_90_fish <- Burb_Home_Range %>% 
  group_by(TRANSMITTER, SEASON) %>% 
  get_summary_stats(HR90, type = "full")

stats_50_fish <- Burb_Home_Range %>% 
  group_by(TRANSMITTER, SEASON) %>% 
  get_summary_stats(HR50, type = "full")

open_stats <- stats_90_fish[stats_90_fish$SEASON =="Open",]
ice_stats <- stats_90_fish[stats_90_fish$SEASON =="Ice",]

################################################################################

library(lme4)
library(MuMIn)
  #https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples
library(lmerTest)
library(car)

str(Burb_Home_Range)

HR_glmm <- glmer(HR90 ~ (1|TRANSMITTER), data = Burb_Home_Range,
                 family = Gamma(link = "log"))
HR_glmm1 <- glmer(HR90 ~ SEX + SEASON + BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
HR_glmm2 <- glmer(HR90 ~ (SEX * SEASON) + BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
HR_glmm3 <- glmer(HR90 ~ SEX + (SEASON * BASIN) + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
HR_glmm4 <- glmer(HR90 ~ (SEX * BASIN) + SEASON + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
HR_glmm5 <- glmer(HR90 ~ SEX + SEASON + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
HR_glmm6 <- glmer(HR90 ~ (SEX * SEASON) + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
HR_glmm7 <- glmer(HR90 ~ SEASON + BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
HR_glmm8 <- glmer(HR90 ~ (SEASON * BASIN) + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
HR_glmm9 <- glmer(HR90 ~ SEX + BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
HR_glmm10 <- glmer(HR90 ~ (SEX * BASIN) + (1|TRANSMITTER), data = Burb_Home_Range,
                   family = Gamma(link = "log"))
HR_glmm11 <- glmer(HR90 ~ SEX + (1|TRANSMITTER), data = Burb_Home_Range,
                   family = Gamma(link = "log"))
HR_glmm12 <- glmer(HR90 ~ SEASON + (1|TRANSMITTER), data = Burb_Home_Range,
                   family = Gamma(link = "identity"))
HR_glmm13 <- glmer(HR90 ~ BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                   family = Gamma(link = "log"))

# use the mod.sel function to conduct model selection
# and put output into object out.put
out.put.hr <- model.sel(HR_glmm,HR_glmm1,HR_glmm2,HR_glmm3,HR_glmm4,HR_glmm5,
                        HR_glmm6,HR_glmm7,HR_glmm8,HR_glmm9,HR_glmm10,
                        HR_glmm11,HR_glmm12,HR_glmm13)
out.put.hr
a <-sw(out.put.hr)
##############
#remotes::install_github("ddsjoberg/gtsummary")
library(gtsummary)

tbl <- tbl_regression(HR_glmm6,exponentiate = TRUE)

###
Burb_Home_Range1 <- Burb_Home_Range

Burb_Home_Range1$SEX <- as.factor(Burb_Home_Range1$SEX)  
Burb_Home_Range1$SEASON <- as.factor(Burb_Home_Range1$SEASON)  
structure(Burb_Home_Range1)

contrasts(Burb_Home_Range1$SEX) <- "contr.sum"
contrasts(Burb_Home_Range1$SEASON) <- "contr.sum"

Anova(HR_glmm6, type = 3)

###
# df has subjects (S), two factors (X1,X2) each w/levels (a,b), and continuous response (Y)
library(multcomp) # for glht
library(emmeans) # for emm, emmeans

summary(glht(HR_glmm6, emm(pairwise ~ SEASON*SEX)), test=adjusted(type="holm"))
# or, using the Tukey HSD correction instead of Holm's
summary(emmeans(HR_glmm6, pairwise ~ SEASON*SEX, adjust="tukey", mode="linear.predictor", type="Score"))

################################################################################
HR_emm <- emmeans(HR_glmm6, ~ SEASON*SEX)

#WANT THIS TABLE
# TABLE 2 Home Range Paper
pairs(HR_emm, simple = "each")

#WANT THIS TABLE
contrast(HR_emm, "consec", simple = "each", combine = TRUE)

#
str(HR_emm)

contrast(emmeans(HR_glmm6, ~ SEASON*SEX),
         interaction = c("poly", "consec", "consec"))

joint_tests(HR_glmm6)
joint_tests(HR_glmm6, by = "SEASON")

mvcontrast(HR_emm, "pairwise", mult.name = c("SEASON"))

update(mvcontrast(HR_emm, "consec", mult.name = "SEASON", by = "SEASON"), 
       by = NULL)
##########################################

############# CORE RANGE ####################
stats_50_table_season_sex <- Burb_Home_Range %>% 
  group_by(SEASON,SEX) %>% 
  get_summary_stats(HR50,type = "full")

stats_50_table_date_sex <- Burb_Home_Range %>% 
  group_by(DATE,SEX) %>% 
  get_summary_stats(HR50,type = "full")

library(lme4)
library(MuMIn)
library(lmerTest)
library(car)

str(Burb_Home_Range)

CR_glmm <- glmer(HR50 ~ (1|TRANSMITTER), data = Burb_Home_Range,
                 family = Gamma(link = "log"))
CR_glmm1 <- glmer(HR50 ~ SEX + SEASON + BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
CR_glmm2 <- glmer(HR50 ~ (SEX * SEASON) + BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
CR_glmm3 <- glmer(HR50 ~ SEX + (SEASON * BASIN) + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
CR_glmm4 <- glmer(HR50 ~ (SEX * BASIN) + SEASON + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
CR_glmm5 <- glmer(HR50 ~ SEX + SEASON + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
CR_glmm6 <- glmer(HR50 ~ (SEX * SEASON) + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
CR_glmm7 <- glmer(HR50 ~ SEASON + BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
CR_glmm8 <- glmer(HR50 ~ (SEASON * BASIN) + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
CR_glmm9 <- glmer(HR50 ~ SEX + BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
CR_glmm10 <- glmer(HR50 ~ (SEX * BASIN) + (1|TRANSMITTER), data = Burb_Home_Range,
                   family = Gamma(link = "log"))
CR_glmm11 <- glmer(HR50 ~ SEX + (1|TRANSMITTER), data = Burb_Home_Range,
                   family = Gamma(link = "log"))
CR_glmm12 <- glmer(HR50 ~ SEASON + (1|TRANSMITTER), data = Burb_Home_Range,
                   family = Gamma(link = "log"))
CR_glmm13 <- glmer(HR50 ~ BASIN + (1|TRANSMITTER), data = Burb_Home_Range,
                   family = Gamma(link = "log"))

# use the mod.sel function to conduct model selection
# and put output into object out.put
out.put.CR <- model.sel(CR_glmm,CR_glmm1,CR_glmm2,CR_glmm3,CR_glmm4,CR_glmm5,
                        CR_glmm6,CR_glmm7,CR_glmm8,CR_glmm9,CR_glmm10,
                        CR_glmm11,CR_glmm12,CR_glmm13)

out.put.CR
out.put.hr
importance(out.put.CR)

#residuals
hist(resid(CR_glmm6))
qqnorm(resid(CR_glmm6))
qqline(resid(CR_glmm6))
plot(fitted(CR_glmm6),resid(CR_glmm6))
abline(h=0,col="grey")
lines(lowess(fitted(CR_glmm6)[is.finite(fitted(CR_glmm6))],resid(CR_glmm6)[is.finite(fitted(CR_glmm6))]),col="red")
#
Anova(CR_glmm6, type = 3)
###
# df has subjects (S), two factors (X1,X2) each w/levels (a,b), and continuous response (Y)
library(multcomp) # for glht
library(emmeans) # for emm, emmeans
summary(glht(CR_glmm6, emm(pairwise ~ SEASON*SEX)), test=adjusted(type="holm"))
# or, using the Tukey HSD correction instead of Holm's
summary(emmeans(CR_glmm6, pairwise ~ SEASON*SEX, adjust="tukey", mode="linear.predictor", type="Score"))

################################################################################
CR_emm <- emmeans(CR_glmm6, ~ SEASON*SEX)


#https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/#:~:text=The%20emmeans%20package%20has%20helper%20functions%20for%20commonly,or%20trt.vs.ctrlk%2C%20and%20even%20consecutive%20comparisons%20via%20consec.
#WANT THIS TABLE
pairs(CR_emm,simple="each")

#WANT THIS TABLE
contrast(CR_emm, "consec", simple = "each", combine = TRUE)

#
str(CR_emm)

contrast(emmeans(CR_glmm6, ~ SEASON*SEX),
         interaction = c("poly", "consec", "consec"))

joint_tests(CR_glmm6)
joint_tests(CR_glmm6, by = "SEASON")

mvcontrast(CR_emm, "pairwise", mult.name = c("SEASON"))

update(mvcontrast(CR_emm, "consec", mult.name = "SEASON", by = "SEASON"), 
       by = NULL)
summary(CR_glmm1)
summary(CR_glmm2)
summary(CR_glmm3)
summary(CR_glmm4)
summary(CR_glmm5)
summary(CR_glmm6)
summary(CR_glmm7)

#comparison by group
Anova(CR_glmm6,type = 2)
Anova(CR_glmm6,type = 3)
Anova(glmer(HR50 ~ (SEASON * SEX) + (1|TRANSMITTER), data = Burb_Home_Range,
            family = Gamma(link = "log"), contrasts=list(SEASON=contr.sum, SEX=contr.sum)), type=3)


###################################################
######### NOT NEEDED ##############################
###################################################
#comparison with emmeans package
HR_glmm6 <- glmer(HR90 ~ (SEX * SEASON) + (1|TRANSMITTER), data = Burb_Home_Range,
                  family = Gamma(link = "log"))
# USE THIS ONE
library(emmeans)
m_eff <- emmeans(HR_glmm6, spec = '(SEASON*SEX)')
m_eff <- emmeans(HR_glmm6, ~ SEX|SEASON)
summary(m_eff)

m_tukey <- contrast(m_eff, method = "pairwise")

summary(m_tukey)

confint(m_tukey)

plot(m_tukey,xlabs = "Difference in Home Range by Sex")
