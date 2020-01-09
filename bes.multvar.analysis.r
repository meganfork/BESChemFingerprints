### R script for multivariate water chem analysis of BES - adapted from previous NMDS/Radar chart scripts
### Created by AJR
### Created on 19 June 2017
### Last edited by: MLF
### Last edited on: 17 December 2019
### Last edits: Updated calculations for new data (WY2000-WY2018), automated removal of bad data points, updated to dplyr functions, added download of NWIS data

#First set your working directory - this should be changed to reflect wherever you have saved the .csv files
#Note that if you use a project (here 'BaltimoreStreamChem' the directory will be set to the project directory)
# setwd("C:/Users/reisingera/Dropbox/Documents/R/Balt-ChemFingerprints")

library(tidyverse)
library(fmsb)
library(car)
library(vegan)
library(animation)
library(viridis)
library(dataRetrieval)
library(plotrix)
library(lubridate)

#function to extract p-value from a linear model:
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


#### Do some data cleanup from the rawdata emailed by Lisa on 13 Dec 2019 (data from 15-Oct-98 through 7-Jun-19):
bes.waterchem.raw<-read.csv(file="BESFullWaterChemRaw.csv", header=T)

bes.waterchem.working<-bes.waterchem.raw[bes.waterchem.raw$Site=='GFCP'|bes.waterchem.raw$Site=='GFVN'|bes.waterchem.raw$Site=='GFGB'|bes.waterchem.raw$Site=='GFGL'|
    bes.waterchem.raw$Site=='DRKR'|bes.waterchem.raw$Site=='POBR'|bes.waterchem.raw$Site=='BARN'|bes.waterchem.raw$Site=='MCDN',]

bes.waterchem.working<-bes.waterchem.working[,c(1:11,13:15)]

write.csv(bes.waterchem.working, file="bes.waterchem.csv",row.names = F)

### Load working data 
waterchem.data<-read.csv(file="bes.waterchem.csv", header=T)

### clean up obviously bad data points
waterchem.clean<-waterchem.data

waterchem.clean <- waterchem.clean %>% mutate(temperature = replace(temperature, temperature > 100, NA)) #replace temperatures >100 with NA
waterchem.clean <- waterchem.clean %>% mutate(pH = replace(pH, pH > 14 | pH <1, NA)) #replace pH greater than 14 or less than 1
waterchem.clean <- waterchem.clean %>% mutate(dox = replace(dox, dox > 100, NA)) #replace dox greater than 100 or less than 1

### remove duplicated rows
waterchem.clean <- distinct(waterchem.clean)

write.csv(waterchem.clean, file="bes.waterchem.clean.csv",row.names = F) #write out the clean data



#Read in the clean working data for analysis
lulc.data<-read.csv(file="beslulc.csv", header=T)#Updated BES LULC from Claire Welty (see email sent from Claire on 11 Jan 2017)

data<-read.csv(file="bes.waterchem.clean.csv", header=T)
data$Date<-dmy(data$Date)

#add a water year column
data$water.year<-rep(NA, length(data[,1]))

#fill the water year. If the julian date is before 275 (Oct 2 I believe, or Oct 1 for leap year...probably could do this better but it works) then the 
#water year for that row is the same as the actual year, if the julian date is >= 275 then the water year assigned is the next calendar year
for(i in 1:length(data$water.year)){
  if(data$Julian_Date[i]<275){
    data$water.year[i]=data$Year[i]
  }
  else{
    (data$water.year[i]=data$Year[i]+1)
  }
}

###Make a color pallette for 19 classes (19 years of data from water year 2000-2018)
colors<-viridis(19)
#colors<-c("gray20","gray40", "gray60", "gray80", rgb(19,149,186,maxColorValue=255),rgb(17,120,153,maxColorValue=255),rgb(17,120,153,maxColorValue=255),
         # rgb(15,91,120,maxColorValue=255),rgb(192,46,29,maxColorValue=255), rgb(217,78,31,maxColorValue=255), rgb(241,108,32,maxColorValue=255), 
         #  rgb(239,139,44,maxColorValue=255),rgb(236,170,56,maxColorValue=255),rgb(235,200,68,maxColorValue=255),rgb(162,184,108,maxColorValue=255))


#calculate the mean, max, and mean for each variable of interest by water year and site
summary.data<-data %>% 
  group_by(water.year, Site) %>% 
  summarize(        mean.Cl=mean(Cl, na.rm=T), max.Cl=max(Cl, na.rm=T), min.Cl=min(Cl, na.rm=T),
                    mean.NO3=mean(NO3, na.rm=T), max.NO3=max(NO3, na.rm=T), min.NO3=min(NO3, na.rm=T),
                    mean.PO4=mean(PO4, na.rm=T), max.PO4=max(PO4, na.rm=T), min.PO4=min(PO4, na.rm=T),
                    mean.SO4=mean(SO4, na.rm=T), max.SO4=max(SO4, na.rm=T), min.SO4=min(SO4, na.rm=T),
                    mean.TN=mean(TN, na.rm=T), max.TN=max(TN, na.rm=T), min.TN=min(TN, na.rm=T),
                    mean.TP=mean(TP, na.rm=T), max.TP=max(TP, na.rm=T), min.TP=min(TP, na.rm=T))

#Remove 2019 from the records as we don't have TN and TP from 1999 nor a full record from 2019
summary.data<-summary.data[summary.data$water.year!=1999 & summary.data$water.year!=2019,]

range.data<-data %>% 
  group_by(water.year, Site) %>% 
  summarize(        chloride.range<-max(Cl, na.rm=T) - min(Cl, na.rm=T),
                    NO3.range<-max(NO3, na.rm=T)-min(NO3, na.rm=T),
                    PO4.range<-max(PO4, na.rm=T)-min(PO4, na.rm=T),
                    SO4.range<-max(SO4, na.rm=T)-min(SO4, na.rm=T),
                    TN.range<-max(TN, na.rm=T)-min(TN, na.rm=T),
                    TP.range<-max(TP, na.rm=T)-min(TP, na.rm=T))
 
range.data<-range.data[range.data$water.year!=1999 & range.data$water.year!=2019,]

colnames(range.data)<-c("water.year", "Site", "Cl", "NO3", "PO4", "SO4", "TN", "TP")



#Make radar charts showing the annual mean values for each site:
windows(height=3, width=6)
par(mfrow=c(2,4), mar=c(0.6,1,0.6,1))
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$Site=='POBR',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16,cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='POBR')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$Site=='BARN',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='BARN')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$Site=='MCDN',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='MCDN')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$Site=='DRKR',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='DRKR')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$Site=='GFGL',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='GFGL')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$Site=='GFGB',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='GFGB')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$Site=='GFVN',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='GFVN')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                 summary.data[summary.data$Site=='GFCP',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                    plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                    expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='GFCP')

means<-summary.data[,c(1:3,6,9,12,15,18)]

means.lulc<-left_join(means,lulc.data,by=c("Site"="site"))


##############################################
###      Download NWIS discharge data      ###
##############################################

USGS.site.nos<-data.frame(matrix(c(
  "POBR","01583570",
  "BARN","01583580",
  "GFGB","01589197",
  "GFVN","01589300",
  "GFCP","01589352",
  "DRKR","01589330",
  "GFGL","01589180"),ncol=2,byrow=T))
colnames(USGS.site.nos)<-c("siteName","siteNumber")
Q.start<-"1998-10-01"
Q.end<-"2019-09-30"


BESflow<-readNWISdv(siteNumbers = USGS.site.nos$siteNumber,parameterCd = '00060',startDate = Q.start,endDate = Q.end)
colnames(BESflow)<-c("agency_cd","site_no","Date","meanDailyQ_cfs","qual_cd")
BESflow$Date<-ymd(BESflow$Date)

BESflow$meanDailyQ_Ls<-BESflow$meanDailyQ_cfs*28.3168


###########################################
###      Annual C-Q relationships       ###
###########################################


## join flow data and chemistry data (looped over sites, not including MCDN which doesn't have flow):

#create a list of sites
CandQ.bysite<-list()

for (i in 1:nrow(USGS.site.nos)){
CandQ.bysite[[i]]<-right_join(select(filter(BESflow,site_no==USGS.site.nos$siteNumber[i]),site_no,Date,meanDailyQ_Ls),filter(data,Site==as.character(USGS.site.nos$siteName[i])),by="Date")
}
names(CandQ.bysite)<-USGS.site.nos$siteName

## calculate the C-Q relationship for each by looping over year within site

BES.water.years<-seq(from=2000,to=2018,by=1)

annualCQ.bysite<-list() 
#initialize the list with blank dataframes
for (i in 1:length(CandQ.bysite)){
  annualCQ.bysite[[i]]<-data.frame(water.year=BES.water.years,
                                   Cl.slope=rep(NA,times=length(BES.water.years)),
                                   Cl.int=rep(NA,times=length(BES.water.years)),
                                   Cl.p.val=rep(NA,times=length(BES.water.years)),
                                   Cl.r2=rep(NA,times=length(BES.water.years)),
                                   Cl.n.obs=rep(NA,times=length(BES.water.years)),
                                   SO4.slope=rep(NA,times=length(BES.water.years)),
                                   SO4.int=rep(NA,times=length(BES.water.years)),
                                   SO4.p.val=rep(NA,times=length(BES.water.years)),
                                   SO4.r2=rep(NA,times=length(BES.water.years)),
                                   SO4.n.obs=rep(NA,times=length(BES.water.years)),
                                   NO3.slope=rep(NA,times=length(BES.water.years)),
                                   NO3.int=rep(NA,times=length(BES.water.years)),
                                   NO3.p.val=rep(NA,times=length(BES.water.years)),
                                   NO3.r2=rep(NA,times=length(BES.water.years)),
                                   NO3.n.obs=rep(NA,times=length(BES.water.years)),
                                   PO4.slope=rep(NA,times=length(BES.water.years)),
                                   PO4.int=rep(NA,times=length(BES.water.years)),
                                   PO4.p.val=rep(NA,times=length(BES.water.years)),
                                   PO4.r2=rep(NA,times=length(BES.water.years)),
                                   PO4.n.obs=rep(NA,times=length(BES.water.years)),
                                   TN.slope=rep(NA,times=length(BES.water.years)),
                                   TN.int=rep(NA,times=length(BES.water.years)),
                                   TN.p.val=rep(NA,times=length(BES.water.years)),
                                   TN.r2=rep(NA,times=length(BES.water.years)),
                                   TN.n.obs=rep(NA,times=length(BES.water.years)),
                                   TP.slope=rep(NA,times=length(BES.water.years)),
                                   TP.int=rep(NA,times=length(BES.water.years)),
                                   TP.p.val=rep(NA,times=length(BES.water.years)),
                                   TP.r2=rep(NA,times=length(BES.water.years)),
                                   TP.n.obs=rep(NA,times=length(BES.water.years)))
}



#calculate C-Q slopes, p-values, and r2 for each constituent for each water year for each site
for (i in 1:length(CandQ.bysite)){
  for (j in 1:length(BES.water.years)){
    
    wateryear.sub<-filter(CandQ.bysite[[i]],water.year==BES.water.years[j])
    
    #For each solute, an if statement returns NA if the number of annual observations for either flow or conc is less than 40; otherwise C-Q is calculated
    #Cl
    if(length(which(is.na(wateryear.sub$Cl)==F))>40 & length(which(is.na(wateryear.sub$meanDailyQ_Ls)==F))>40){
      Cl.lm<-lm(log10(wateryear.sub$Cl)~log10(wateryear.sub$meanDailyQ_Ls+0.283168/2),na.action = na.omit)
      annualCQ.bysite[[i]]$Cl.slope[j]<-summary(Cl.lm)$coefficients[2]
      annualCQ.bysite[[i]]$Cl.int[j]<-summary(Cl.lm)$coefficients[1]
      annualCQ.bysite[[i]]$Cl.p.val[j]<-lmp(Cl.lm)
      annualCQ.bysite[[i]]$Cl.r2[j]<-summary(Cl.lm)$adj.r.squared
      annualCQ.bysite[[i]]$Cl.n.obs[j]<-length(which(is.na(wateryear.sub$Cl)==F))
    } else {
      annualCQ.bysite[[i]]$Cl.slope[j]<-NA
      annualCQ.bysite[[i]]$Cl.int[j]<-NA
      annualCQ.bysite[[i]]$Cl.p.val[j]<-NA
      annualCQ.bysite[[i]]$Cl.r2[j]<-NA
      annualCQ.bysite[[i]]$Cl.n.obs[j]<-NA
    }
        
    
    #SO4
    if(length(which(is.na(wateryear.sub$SO4)==F))>40 & length(which(is.na(wateryear.sub$meanDailyQ_Ls)==F))>40){
      SO4.lm<-lm(log10(wateryear.sub$SO4)~log10(wateryear.sub$meanDailyQ_Ls+0.283168/2),na.action = na.omit)
      annualCQ.bysite[[i]]$SO4.slope[j]<-summary(SO4.lm)$coefficients[2]
      annualCQ.bysite[[i]]$SO4.int[j]<-summary(SO4.lm)$coefficients[1]
      annualCQ.bysite[[i]]$SO4.p.val[j]<-lmp(SO4.lm)
      annualCQ.bysite[[i]]$SO4.r2[j]<-summary(SO4.lm)$adj.r.squared
      annualCQ.bysite[[i]]$SO4.n.obs[j]<-length(which(is.na(wateryear.sub$SO4)==F))
    } else {
      annualCQ.bysite[[i]]$SO4.slope[j]<-NA
      annualCQ.bysite[[i]]$SO4.int[j]<-NA
      annualCQ.bysite[[i]]$SO4.p.val[j]<-NA
      annualCQ.bysite[[i]]$SO4.r2[j]<-NA
      annualCQ.bysite[[i]]$SO4.n.obs[j]<-NA
    }
    
    #NO3
    if(length(which(is.na(wateryear.sub$NO3)==F))>40 & length(which(is.na(wateryear.sub$meanDailyQ_Ls)==F))>40){
      NO3.lm<-lm(log10(wateryear.sub$NO3+0.005)~log10(wateryear.sub$meanDailyQ_Ls+0.283168/2),na.action = na.omit)
      annualCQ.bysite[[i]]$NO3.slope[j]<-summary(NO3.lm)$coefficients[2]
      annualCQ.bysite[[i]]$NO3.int[j]<-summary(NO3.lm)$coefficients[1]
      annualCQ.bysite[[i]]$NO3.p.val[j]<-lmp(NO3.lm)
      annualCQ.bysite[[i]]$NO3.r2[j]<-summary(NO3.lm)$adj.r.squared
      annualCQ.bysite[[i]]$NO3.n.obs[j]<-length(which(is.na(wateryear.sub$NO3)==F))
    } else {
      annualCQ.bysite[[i]]$NO3.slope[j]<-NA
      annualCQ.bysite[[i]]$NO3.int[j]<-NA
      annualCQ.bysite[[i]]$NO3.p.val[j]<-NA
      annualCQ.bysite[[i]]$NO3.r2[j]<-NA
      annualCQ.bysite[[i]]$NO3.n.obs[j]<-NA
    }
    
    #PO4
    if(length(which(is.na(wateryear.sub$PO4)==F))>40 & length(which(is.na(wateryear.sub$meanDailyQ_Ls)==F))>40){
      PO4.lm<-lm(log10(wateryear.sub$PO4)~log10(wateryear.sub$meanDailyQ_Ls+0.283168/2),na.action = na.omit)
      annualCQ.bysite[[i]]$PO4.slope[j]<-summary(PO4.lm)$coefficients[2]
      annualCQ.bysite[[i]]$PO4.int[j]<-summary(PO4.lm)$coefficients[1]
      annualCQ.bysite[[i]]$PO4.p.val[j]<-lmp(PO4.lm)
      annualCQ.bysite[[i]]$PO4.r2[j]<-summary(PO4.lm)$adj.r.squared
      annualCQ.bysite[[i]]$PO4.n.obs[j]<-length(which(is.na(wateryear.sub$PO4)==F))
    } else {
      annualCQ.bysite[[i]]$PO4.slope[j]<-NA
      annualCQ.bysite[[i]]$PO4.int[j]<-NA
      annualCQ.bysite[[i]]$PO4.p.val[j]<-NA
      annualCQ.bysite[[i]]$PO4.r2[j]<-NA
      annualCQ.bysite[[i]]$PO4.n.obs[j]<-NA
    }
    
   #TN
    if(length(which(is.na(wateryear.sub$TN)==F))>40 & length(which(is.na(wateryear.sub$meanDailyQ_Ls)==F))>40){
    TN.lm<-lm(log10(wateryear.sub$TN)~log10(wateryear.sub$meanDailyQ_Ls+0.283168/2),na.action = na.omit)
    annualCQ.bysite[[i]]$TN.slope[j]<-summary(TN.lm)$coefficients[2]
    annualCQ.bysite[[i]]$TN.int[j]<-summary(TN.lm)$coefficients[1]
    annualCQ.bysite[[i]]$TN.p.val[j]<-lmp(TN.lm)
    annualCQ.bysite[[i]]$TN.r2[j]<-summary(TN.lm)$adj.r.squared
    annualCQ.bysite[[i]]$TN.n.obs[j]<-length(which(is.na(wateryear.sub$TN)==F))
    } else {
      annualCQ.bysite[[i]]$TN.slope[j]<-NA
      annualCQ.bysite[[i]]$TN.int[j]<-NA
      annualCQ.bysite[[i]]$TN.p.val[j]<-NA
      annualCQ.bysite[[i]]$TN.r2[j]<-NA
      annualCQ.bysite[[i]]$TN.n.obs[j]<-NA
    }
    
    #TP
    if(length(which(is.na(wateryear.sub$TP)==F))>40 & length(which(is.na(wateryear.sub$meanDailyQ_Ls)==F))>40){
      TP.lm<-lm(log10(wateryear.sub$TP)~log10(wateryear.sub$meanDailyQ_Ls+0.283168/2),na.action = na.omit)
      annualCQ.bysite[[i]]$TP.slope[j]<-summary(TP.lm)$coefficients[2]
      annualCQ.bysite[[i]]$TP.int[j]<-summary(TP.lm)$coefficients[1]
      annualCQ.bysite[[i]]$TP.p.val[j]<-lmp(TP.lm)
      annualCQ.bysite[[i]]$TP.r2[j]<-summary(TP.lm)$adj.r.squared
      annualCQ.bysite[[i]]$TP.n.obs[j]<-length(which(is.na(wateryear.sub$TP)==F))
      } else {
      annualCQ.bysite[[i]]$TP.slope[j]<-NA
      annualCQ.bysite[[i]]$TP.int[j]<-NA
      annualCQ.bysite[[i]]$TP.p.val[j]<-NA
      annualCQ.bysite[[i]]$TP.r2[j]<-NA
      annualCQ.bysite[[i]]$TP.n.obs[j]<-NA
      }
  }
}
names(annualCQ.bysite)<-USGS.site.nos$siteName

for (i in 1:length(annualCQ.bysite)){
  annualCQ.bysite[[i]]$Site<-names(annualCQ.bysite)[i]
}

annualCQ<-do.call("rbind",annualCQ.bysite) #make "by-site" list into one dataframe


### make figures for each site+solute combination showing all years in viridris
solutes<-c("Cl","NO3","PO4","SO4","TN","TP")

for (i in 1:length(CandQ.bysite)){
  for (j in 1:length(solutes)){
    
    png(paste0("Figures/C-Q visualizations/",names(CandQ.bysite)[i],"_",solutes[j],"annualCQ_00to18.png"),height=10,width=10,units='cm',res=300)
    par(mar = c(2.5,3,0.5,0.5),
        mgp=c(1.3,0.2,0))
    
    plot(log10(CandQ.bysite[[i]][,which(colnames(CandQ.bysite[[i]])==solutes[j])]+0.005)~log10(CandQ.bysite[[i]][,"meanDailyQ_Ls"]+0.283168/2),
         xlab=expression(paste("log(Q) (L s"^-1,")")),ylab=substitute(expression("log("*nn*") (mg L"^-1*")"),list(nn=solutes[j])),
         xlim=c(-1,5),ylim=c(-2.5,3.5),axes=F,type="n")
    axis(1,tck=0.02,cex.axis=0.85)
    axis(2,tck=0.02,cex.axis=0.85)
    box()
    text(2,3.4,names(CandQ.bysite)[i],font=2)
    
    for (k in 1:length(BES.water.years)){
    wateryear.sub<-filter(CandQ.bysite[[i]],water.year==BES.water.years[k])
    
    points(log10(wateryear.sub[,which(colnames(CandQ.bysite[[i]])==solutes[j])]+0.005)~log10(wateryear.sub[,"meanDailyQ_Ls"]+0.283168/2),pch=1,col=colors[k])
    
    
    #if statement based on presence and significance of C-Q relationship for plotting annual ablines
    
    if(is.na(annualCQ.bysite[[i]][k,which(colnames(annualCQ.bysite[[i]])==paste0(solutes[j],".p.val"))])==F & 
       annualCQ.bysite[[i]][k,which(colnames(annualCQ.bysite[[i]])==paste0(solutes[j],".p.val"))]<0.05){
      
          abline(annualCQ.bysite[[i]][k,which(colnames(annualCQ.bysite[[i]])==paste0(solutes[j],".int"))],
                 annualCQ.bysite[[i]][k,which(colnames(annualCQ.bysite[[i]])==paste0(solutes[j],".slope"))],
                 col=colors[k],lwd=2)  
      } else if (is.na(annualCQ.bysite[[i]][k,which(colnames(annualCQ.bysite[[i]])==paste0(solutes[j],".p.val"))])==F & 
               annualCQ.bysite[[i]][k,which(colnames(annualCQ.bysite[[i]])==paste0(solutes[j],".p.val"))]>0.05){
          
          abline(annualCQ.bysite[[i]][k,which(colnames(annualCQ.bysite[[i]])==paste0(solutes[j],".int"))],
                 annualCQ.bysite[[i]][k,which(colnames(annualCQ.bysite[[i]])==paste0(solutes[j],".slope"))],
                 col=colors[k],lty=2)  
      }
    
    }
    dev.off()
  }
}


###################################################
###      Annual CV of chem vs. CV of flow       ###
###################################################

annualCVC.CVQ.bysite<-list()
#initialize the list with blank dataframes
for (i in 1:length(CandQ.bysite)){
  annualCVC.CVQ.bysite[[i]]<-data.frame(water.year=BES.water.years,
                                   flow.CV=rep(NA,times=length(BES.water.years)),
                                   Cl.CV=rep(NA,times=length(BES.water.years)),
                                   Cl.flow.CVratio=rep(NA,times=length(BES.water.years)),
                                   SO4.CV=rep(NA,times=length(BES.water.years)),
                                   SO4.flow.CVratio=rep(NA,times=length(BES.water.years)),
                                   NO3.CV=rep(NA,times=length(BES.water.years)),
                                   NO3.flow.CVratio=rep(NA,times=length(BES.water.years)),
                                   PO4.CV=rep(NA,times=length(BES.water.years)),
                                   PO4.flow.CVratio=rep(NA,times=length(BES.water.years)),
                                   TN.CV=rep(NA,times=length(BES.water.years)),
                                   TN.flow.CVratio=rep(NA,times=length(BES.water.years)),
                                   TP.CV=rep(NA,times=length(BES.water.years)),
                                   TP.flow.CVratio=rep(NA,times=length(BES.water.years)))
}



#calculate CV of flow, CV of concentration, and the ratio of these for each constituent for each water year for each site
for (i in 1:length(CandQ.bysite)){
  for (j in 1:length(BES.water.years)){
    
    wateryear.sub<-filter(CandQ.bysite[[i]],water.year==BES.water.years[j])
    
    if(length(which(is.na(wateryear.sub$meanDailyQ_Ls)==F))>40){
    annualCVC.CVQ.bysite[[i]]$flow.CV[j]<-sd(wateryear.sub$meanDailyQ_Ls)/mean(wateryear.sub$meanDailyQ_Ls)
    } else {
      annualCVC.CVQ.bysite[[i]]$flow.CV[j]<-NA
    }
    
    
    #For each solute, an if statement returns NA if the number of annual observations for conc is less than 40; otherwise CV and ratios are calculated
    #Cl
    if(length(which(is.na(wateryear.sub$Cl)==F))>40){
      annualCVC.CVQ.bysite[[i]]$Cl.CV[j]<-sd(wateryear.sub$Cl,na.rm = T)/mean(wateryear.sub$Cl,na.rm=T)
    } else {
      annualCVC.CVQ.bysite[[i]]$Cl.CV[j]<-NA
    }
    
    if(is.na(annualCVC.CVQ.bysite[[i]]$flow.CV[j])==F & is.na(annualCVC.CVQ.bysite[[i]]$Cl.CV[j])==F){
      annualCVC.CVQ.bysite[[i]]$Cl.flow.CVratio[j]<-annualCVC.CVQ.bysite[[i]]$Cl.CV[j]/annualCVC.CVQ.bysite[[i]]$flow.CV[j]
    }
    
    #SO4
    if(length(which(is.na(wateryear.sub$SO4)==F))>40){
      annualCVC.CVQ.bysite[[i]]$SO4.CV[j]<-sd(wateryear.sub$SO4,na.rm=T)/mean(wateryear.sub$SO4,na.rm=T)
    } else {
      annualCVC.CVQ.bysite[[i]]$SO4.CV[j]<-NA
    }
    
    if(is.na(annualCVC.CVQ.bysite[[i]]$flow.CV[j])==F & is.na(annualCVC.CVQ.bysite[[i]]$SO4.CV[j])==F){
      annualCVC.CVQ.bysite[[i]]$SO4.flow.CVratio[j]<-annualCVC.CVQ.bysite[[i]]$SO4.CV[j]/annualCVC.CVQ.bysite[[i]]$flow.CV[j]
    }
    
    #NO3
    if(length(which(is.na(wateryear.sub$NO3)==F))>40){
      annualCVC.CVQ.bysite[[i]]$NO3.CV[j]<-sd(wateryear.sub$NO3,na.rm=T)/mean(wateryear.sub$NO3,na.rm=T)
    } else {
      annualCVC.CVQ.bysite[[i]]$NO3.CV[j]<-NA
    }
    
    if(is.na(annualCVC.CVQ.bysite[[i]]$flow.CV[j])==F & is.na(annualCVC.CVQ.bysite[[i]]$NO3.CV[j])==F){
      annualCVC.CVQ.bysite[[i]]$NO3.flow.CVratio[j]<-annualCVC.CVQ.bysite[[i]]$NO3.CV[j]/annualCVC.CVQ.bysite[[i]]$flow.CV[j]
    }
    
    #PO4
    if(length(which(is.na(wateryear.sub$PO4)==F))>40){
      annualCVC.CVQ.bysite[[i]]$PO4.CV[j]<-sd(wateryear.sub$PO4,na.rm=T)/mean(wateryear.sub$PO4,na.rm=T)
    } else {
      annualCVC.CVQ.bysite[[i]]$PO4.CV[j]<-NA
    }
    
    if(is.na(annualCVC.CVQ.bysite[[i]]$flow.CV[j])==F & is.na(annualCVC.CVQ.bysite[[i]]$PO4.CV[j])==F){
      annualCVC.CVQ.bysite[[i]]$PO4.flow.CVratio[j]<-annualCVC.CVQ.bysite[[i]]$PO4.CV[j]/annualCVC.CVQ.bysite[[i]]$flow.CV[j]
    }
    
    #TN
    if(length(which(is.na(wateryear.sub$TN)==F))>40){
      annualCVC.CVQ.bysite[[i]]$TN.CV[j]<-sd(wateryear.sub$TN,na.rm=T)/mean(wateryear.sub$TN,na.rm=T)
    } else {
      annualCVC.CVQ.bysite[[i]]$TN.CV[j]<-NA
    }
    
    if(is.na(annualCVC.CVQ.bysite[[i]]$flow.CV[j])==F & is.na(annualCVC.CVQ.bysite[[i]]$TN.CV[j])==F){
      annualCVC.CVQ.bysite[[i]]$TN.flow.CVratio[j]<-annualCVC.CVQ.bysite[[i]]$TN.CV[j]/annualCVC.CVQ.bysite[[i]]$flow.CV[j]
    }
    
    #TP
    if(length(which(is.na(wateryear.sub$TP)==F))>40){
      annualCVC.CVQ.bysite[[i]]$TP.CV[j]<-sd(wateryear.sub$TP,na.rm=T)/mean(wateryear.sub$TP,na.rm=T)
    } else {
      annualCVC.CVQ.bysite[[i]]$TP.CV[j]<-NA
    }
    
    if(is.na(annualCVC.CVQ.bysite[[i]]$flow.CV[j])==F & is.na(annualCVC.CVQ.bysite[[i]]$TP.CV[j])==F){
      annualCVC.CVQ.bysite[[i]]$TP.flow.CVratio[j]<-annualCVC.CVQ.bysite[[i]]$TP.CV[j]/annualCVC.CVQ.bysite[[i]]$flow.CV[j]
    }
  }
}
names(annualCVC.CVQ.bysite)<-USGS.site.nos$siteName

for (i in 1:length(annualCVC.CVQ.bysite)){
  annualCVC.CVQ.bysite[[i]]$Site<-names(annualCVC.CVQ.bysite)[i]
}

annualCVratios<-do.call("rbind",annualCVC.CVQ.bysite) #make "by-site" list into one dataframe





# Multisolute Analysis Plan:
# Run PCAs for the following annual summaries of BES stream chem: annual mean, annual range, annual min, annual max, annual C-Q slope, annual CV of chem:annual CV of flow
# plot each of these using viridis color ramp to differentiate among years (show progression through time, maybe add arrows between successive years?)
# Perhaps also run PCAs for all the points within each site, and compare hull areas or centroids for each year to get a sense of how sampling "snapshots" of water chem are changing
# think about how to compare reactive and conservative solutes (CV ratio?)
# think about some way to summarize multi-solute exceeedence probabilities? 
# Changing flow thresholds for N going from reactive to conservative behavior
# Adding pharma data when we get the data back from Jerker?

################################################################
###              Annual Ordinations (minus MCDN)             ###
################################################################
site.colors<-c('#1b9e77',
               '#d95f02',
               '#7570b3',
               '#e7298a',
               '#66a61e',
               '#e6ab02',
               '#666666')
site.symbols<-c(3,16,25,15,24,18,1,4)



#### 1. Annual means 
means.minusMCDN<-filter(means.lulc,Site!="MCDN") #remove the ag site from analysis

mean.pca<-rda(means.minusMCDN[,3:8], scale=T)
# PC1 explains 49.42% of the variation, PC2 explains 29.46% of the variation

scores.mean<-scores(mean.pca)
row.names(scores.mean$species)<-c("Cl","NO3","PO4","SO4","TN","TP")


## Biplot A: site by color with loadings

png("FIGURES/PCAonAnnualMeans(noMCDN)_A.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(mean.pca,xlab="PC 1 (49.4%)",ylab="PC 2 (29.5%)",type="n",axes=F,xlim=c(-1,2.2),ylim=c(-2.2,1))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

arrows(0,0,scores.mean$species[,1],scores.mean$species[,2],length=0.1,angle=30)
Nudge<-1.05
text(scores.mean$species[,1]*Nudge,scores.mean$species[,2]*Nudge,rownames(scores.mean$species),cex=0.8,font=2)
  
# add sites as colors

for (i in 1:length(unique(means.minusMCDN$Site))){
points(as.vector(scores.mean$sites[which(means.minusMCDN$Site==unique(means.minusMCDN$Site)[i]),1]),as.vector(scores.mean$sites[which(means.minusMCDN$Site==unique(means.minusMCDN$Site)[i]),2]),col=site.colors[i])
}
legend('bottomleft',legend = unique(means.minusMCDN$Site),col=site.colors,pch=1,pt.lwd=2)

dev.off()


## Biplot B: zoomed in with site by symbol and year by color

png("FIGURES/PCAonAnnualMeans(noMCDN)_B.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(mean.pca,xlab="PC 1 (49.4%)",ylab="PC 2 (29.5%)",type="n",axes=F,xlim=c(-1,1),ylim=c(-1.5,1))
axis(1,tck=0.015)
axis(2,tck=0.015)
box()

# add points with site indicated by symbol and year by color (viridis color palette for years is names "colors")

for (i in 1:length(unique(means.minusMCDN$Site))){
  for (j in 1:length(unique(means.minusMCDN$water.year))){
    points(as.vector(scores.mean$sites[which(means.minusMCDN$Site==unique(means.minusMCDN$Site)[i]&means.minusMCDN$water.year==unique(means.minusMCDN$water.year)[j]),1]),as.vector(scores.mean$sites[which(means.minusMCDN$Site==unique(means.minusMCDN$Site)[i]&means.minusMCDN$water.year==unique(means.minusMCDN$water.year)[j]),2]),pch=site.symbols[i],col=colors[j],bg=colors[j])
  }
}

legend('bottomleft',legend = unique(means.minusMCDN$Site),pch=site.symbols,pt.bg = "black",bty='n',cex=0.9)
gradient.rect(-1.3,-0.7,-0.85,-0.6,col=colors)

text(-1.3,-0.57,"2000",adj=c(0,0),cex=0.7)
text(-0.85,-0.57,"2018",adj=c(1,0),cex=0.7)


dev.off()



#### 2. Annual range 
range.minusMCDN<-filter(range.data,Site!="MCDN") #remove the ag site from analysis

range.pca<-rda(range.minusMCDN[,3:8], scale=T)
# PC1 explains 37.2% of the variation, PC2 explains 20.11% of the variation

scores.range<-scores(range.pca)
row.names(scores.range$species)<-c("Cl","NO3","PO4","SO4","TN","TP")


## Biplot A: site by color with loadings

png("FIGURES/PCAonAnnualRanges(noMCDN)_A.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(range.pca,xlab="PC 1 (37.2%)",ylab="PC 2 (20.1%)",type="n",axes=F,xlim=c(-1,2),ylim=c(-2,1.5))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

arrows(0,0,scores.range$species[,1],scores.range$species[,2],length=0.1,angle=30)
Nudge<-1.05
text(scores.range$species[,1]*Nudge,scores.range$species[,2]*Nudge,rownames(scores.range$species),cex=0.8,font=2)

# add sites as colors

for (i in 1:length(unique(range.minusMCDN$Site))){
  points(as.vector(scores.range$sites[which(range.minusMCDN$Site==unique(range.minusMCDN$Site)[i]),1]),as.vector(scores.range$sites[which(range.minusMCDN$Site==unique(range.minusMCDN$Site)[i]),2]),col=site.colors[i])
}
legend('bottomleft',legend = unique(range.minusMCDN$Site),col=site.colors,pch=1,pt.lwd=2)

dev.off()


## Biplot B: zoomed in with site by symbol and year by color

png("FIGURES/PCAonAnnualRanges(noMCDN)_B.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(range.pca,xlab="PC 1 (37.2%)",ylab="PC 2 (20.1%)",type="n",axes=F,xlim=c(-0.9,1.4),ylim=c(-1.8,1))
axis(1,tck=0.015)
axis(2,tck=0.015)
box()

# add points with site indicated by symbol and year by color (viridis color palette for years is names "colors")

for (i in 1:length(unique(range.minusMCDN$Site))){
  for (j in 1:length(unique(range.minusMCDN$water.year))){
    points(as.vector(scores.range$sites[which(range.minusMCDN$Site==unique(range.minusMCDN$Site)[i]&range.minusMCDN$water.year==unique(range.minusMCDN$water.year)[j]),1]),as.vector(scores.range$sites[which(range.minusMCDN$Site==unique(range.minusMCDN$Site)[i]&range.minusMCDN$water.year==unique(range.minusMCDN$water.year)[j]),2]),pch=site.symbols[i],col=colors[j],bg=colors[j])
  }
}

legend('bottomleft',legend = unique(range.minusMCDN$Site),pch=site.symbols,pt.bg = "black",bty='n',cex=0.9)
gradient.rect(-1.2,-0.8,-0.7,-0.7,col=colors)

text(-1.2,-0.67,"2000",adj=c(0,0),cex=0.7)
text(-0.7,-0.67,"2018",adj=c(1,0),cex=0.7)


dev.off()




#### 3. Annual maxima 
summary.minusMCDN<-filter(summary.data,Site!="MCDN") #remove the ag site from analysis

max.pca<-rda(summary.minusMCDN[,seq(from=4,by=3,length.out = 6)], scale=T)
# PC1 explains 36.79% of the variation, PC2 explains 20.54% of the variation

scores.max<-scores(max.pca)
row.names(scores.max$species)<-c("Cl","NO3","PO4","SO4","TN","TP")


## Biplot A: site by color with loadings

png("FIGURES/PCAonAnnualMaxima(noMCDN)_A.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(max.pca,xlab="PC 1 (36.8%)",ylab="PC 2 (20.5%)",type="n",axes=F,xlim=c(-1,1.8),ylim=c(-1.8,2))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

arrows(0,0,scores.max$species[,1],scores.max$species[,2],length=0.1,angle=30)
Nudge<-1.05
text(scores.max$species[,1]*Nudge,scores.max$species[,2]*Nudge,rownames(scores.max$species),cex=0.8,font=2)

# add sites as colors

for (i in 1:length(unique(summary.minusMCDN$Site))){
  points(as.vector(scores.max$sites[which(summary.minusMCDN$Site==unique(summary.minusMCDN$Site)[i]),1]),as.vector(scores.max$sites[which(summary.minusMCDN$Site==unique(summary.minusMCDN$Site)[i]),2]),col=site.colors[i])
}
legend('topleft',legend = unique(summary.minusMCDN$Site),col=site.colors,pch=1,pt.lwd=2)

dev.off()


## Biplot B: zoomed in with site by symbol and year by color

png("FIGURES/PCAonAnnualMaxima(noMCDN)_B.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(max.pca,xlab="PC 1 (36.8%)",ylab="PC 2 (20.5%)",type="n",axes=F,xlim=c(-1,1.8),ylim=c(-1.8,2))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

# add points with site indicated by symbol and year by color (viridis color palette for years is names "colors")

for (i in 1:length(unique(summary.minusMCDN$Site))){
  for (j in 1:length(unique(summary.minusMCDN$water.year))){
    points(as.vector(scores.max$sites[which(summary.minusMCDN$Site==unique(summary.minusMCDN$Site)[i]&summary.minusMCDN$water.year==unique(summary.minusMCDN$water.year)[j]),1]),as.vector(scores.max$sites[which(summary.minusMCDN$Site==unique(summary.minusMCDN$Site)[i]&summary.minusMCDN$water.year==unique(summary.minusMCDN$water.year)[j]),2]),pch=site.symbols[i],col=colors[j],bg=colors[j])
  }
}

legend('topleft',legend = unique(summary.minusMCDN$Site),pch=site.symbols,pt.bg = "black",bty='n',cex=0.9)
gradient.rect(-1.6,0.775,-0.9,0.65,col=colors)

text(-1.6,0.55,"2000",adj=c(0,0),cex=0.7)
text(-0.9,0.55,"2018",adj=c(1,0),cex=0.7)


dev.off()



#### 4. Annual minima 
min.pca<-rda(summary.minusMCDN[,seq(from=5,by=3,length.out = 6)], scale=T)
# PC1 explains 40.73% of the variation, PC2 explains 28.38% of the variation

scores.min<-scores(min.pca)
row.names(scores.min$species)<-c("Cl","NO3","PO4","SO4","TN","TP")


## Biplot A: site by color with loadings

png("FIGURES/PCAonAnnualMinima(noMCDN)_A.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(min.pca,xlab="PC 1 (40.7%)",ylab="PC 2 (28.4%)",type="n",axes=F,xlim=c(-1,2),ylim=c(-1.5,1.75))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

arrows(0,0,scores.min$species[,1],scores.min$species[,2],length=0.1,angle=30)
Nudge<-1.05
text(scores.min$species[,1]*Nudge,scores.min$species[,2]*Nudge,rownames(scores.min$species),cex=0.8,font=2)

# add sites as colors

for (i in 1:length(unique(summary.minusMCDN$Site))){
  points(as.vector(scores.min$sites[which(summary.minusMCDN$Site==unique(summary.minusMCDN$Site)[i]),1]),as.vector(scores.min$sites[which(summary.minusMCDN$Site==unique(summary.minusMCDN$Site)[i]),2]),col=site.colors[i])
}
legend('bottomleft',legend = unique(summary.minusMCDN$Site),col=site.colors,pch=1,pt.lwd=2)

dev.off()


## Biplot B: zoomed in with site by symbol and year by color

png("FIGURES/PCAonAnnualMinima(noMCDN)_B.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(min.pca,xlab="PC 1 (40.7%)",ylab="PC 2 (28.4%)",type="n",axes=F,xlim=c(-1,1.75),ylim=c(-0.8,1.5))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

# add points with site indicated by symbol and year by color (viridis color palette for years is names "colors")

for (i in 1:length(unique(summary.minusMCDN$Site))){
  for (j in 1:length(unique(summary.minusMCDN$water.year))){
    points(as.vector(scores.min$sites[which(summary.minusMCDN$Site==unique(summary.minusMCDN$Site)[i]&summary.minusMCDN$water.year==unique(summary.minusMCDN$water.year)[j]),1]),as.vector(scores.min$sites[which(summary.minusMCDN$Site==unique(summary.minusMCDN$Site)[i]&summary.minusMCDN$water.year==unique(summary.minusMCDN$water.year)[j]),2]),pch=site.symbols[i],col=colors[j],bg=colors[j])
  }
}

legend('bottomleft',legend = unique(summary.minusMCDN$Site),pch=site.symbols,pt.bg = "black",bty='n',cex=0.9)
gradient.rect(-0.55,-1.05,-0.05,-0.95,col=colors)

text(-0.55,-0.9,"2000",adj=c(0,0),cex=0.7)
text(-0.05,-0.9,"2018",adj=c(1,0),cex=0.7)


dev.off()


#### 5. Annual C-Q slopes (not including POBR or MCDN and omitting TP and PO4 data)
annualCQforPCA<-filter(annualCQ,Site!="POBR" & Site!="MCDN")
annualCQforPCA<-filter(annualCQforPCA,is.na(Cl.slope)==F & is.na(NO3.slope)==F & is.na(SO4.slope)==F & is.na(TN.slope)==F)

CQ.pca<-rda(annualCQforPCA[,c(2,6,10,18)],scale=T)
# PC1 explains 51.54 of the variation, PC2 explains 33.50% of the variation

scores.CQ<-scores(CQ.pca)
row.names(scores.CQ$species)<-c("Cl","SO4","NO3","TN")


## Biplot A: site by color with loadings

png("FIGURES/PCAonAnnualCQslopes_A.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(CQ.pca,xlab="PC 1 (51.5%)",ylab="PC 2 (33.5%)",type="n",axes=F,xlim=c(-1.2,2.3),ylim=c(-2.2,0.9))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

arrows(0,0,scores.CQ$species[,1],scores.CQ$species[,2],length=0.1,angle=30)
Nudge<-1.05
text(scores.CQ$species[,1]*Nudge,scores.CQ$species[,2]*Nudge,rownames(scores.CQ$species),cex=0.8,font=2)

# add sites as colors

for (i in 1:length(unique(annualCQforPCA$Site))){
  points(as.vector(scores.CQ$sites[which(annualCQforPCA$Site==unique(annualCQforPCA$Site)[i]),1]),as.vector(scores.CQ$sites[which(annualCQforPCA$Site==unique(annualCQforPCA$Site)[i]),2]),col=site.colors[i])
}
legend('bottomright',legend = unique(annualCQforPCA$Site),col=site.colors,pch=1,pt.lwd=2)

dev.off()


## Biplot B: zoomed in with site by symbol and year by color

png("FIGURES/PCAonAnnualCQslopes_B.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(CQ.pca,xlab="PC 1 (51.5%)",ylab="PC 2 (33.5%)",type="n",axes=F,xlim=c(-1.2,1.9),ylim=c(-1,0.9))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

# add points with site indicated by symbol and year by color (viridis color palette for years is names "colors")

for (i in 1:length(unique(annualCQforPCA$Site))){
  for (j in 1:length(unique(annualCQforPCA$water.year))){
    points(as.vector(scores.CQ$sites[which(annualCQforPCA$Site==unique(annualCQforPCA$Site)[i]&annualCQforPCA$water.year==unique(annualCQforPCA$water.year)[j]),1]),as.vector(scores.CQ$sites[which(annualCQforPCA$Site==unique(annualCQforPCA$Site)[i]&annualCQforPCA$water.year==unique(annualCQforPCA$water.year)[j]),2]),pch=site.symbols[i],col=colors[j],bg=colors[j])
  }
}

legend('bottomright',legend = unique(annualCQforPCA$Site),pch=site.symbols,pt.bg = "black",bty='n',cex=0.9)
gradient.rect(1.4,-0.7,2,-0.6,col=colors)

text(1.4,-0.55,"2000",adj=c(0,0),cex=0.7)
text(2,-0.55,"2018",adj=c(1,0),cex=0.7)


dev.off()


#### 6. Annual CVs of solute concentrations
annualconcCVforPCA<-filter(annualCVratios,is.na(Cl.CV)==F & is.na(NO3.CV)==F & is.na(SO4.CV)==F & is.na(TN.CV)==F & is.na(TP.CV)==F & is.na(PO4.CV)==F)


concCV.pca<-rda(annualconcCVforPCA[,seq(from=3,by=2,length.out = 6)], scale=T)
# PC1 explains 35.95% of the variation, PC2 explains 21.03% of the variation

scores.concCV<-scores(concCV.pca)
row.names(scores.concCV$species)<-c("Cl","SO4","NO3","PO4","TN","TP")


## Biplot A: site by color with loadings

png("FIGURES/PCAonAnnualCVofsoluteconc_A.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(min.pca,xlab="PC 1 (36.0%)",ylab="PC 2 (21.0%)",type="n",axes=F,xlim=c(-2,1),ylim=c(-2.1,1.25))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

arrows(0,0,scores.concCV$species[,1],scores.concCV$species[,2],length=0.1,angle=30)
Nudge<-1.05
text(scores.concCV$species[,1]*Nudge,scores.concCV$species[,2]*Nudge,rownames(scores.concCV$species),cex=0.8,font=2)

# add sites as colors

for (i in 1:length(unique(annualconcCVforPCA$Site))){
  points(as.vector(scores.concCV$sites[which(annualconcCVforPCA$Site==unique(annualconcCVforPCA$Site)[i]),1]),as.vector(scores.concCV$sites[which(annualconcCVforPCA$Site==unique(annualconcCVforPCA$Site)[i]),2]),col=site.colors[i])
}
legend('bottomright',legend = unique(annualconcCVforPCA$Site),col=site.colors,pch=1,pt.lwd=2)

dev.off()


## Biplot B: zoomed in with site by symbol and year by color

png("FIGURES/PCAonAnnualCVofsoluteconc_B.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(min.pca,xlab="PC 1 (36.0%)",ylab="PC 2 (21.0%)",type="n",axes=F,xlim=c(-1.8,0.6),ylim=c(-1.2,1.1))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

# add points with site indicated by symbol and year by color (viridis color palette for years is names "colors")

for (i in 1:length(unique(annualconcCVforPCA$Site))){
  for (j in 1:length(unique(annualconcCVforPCA$water.year))){
    points(as.vector(scores.concCV$sites[which(annualconcCVforPCA$Site==unique(annualconcCVforPCA$Site)[i]&annualconcCVforPCA$water.year==unique(annualconcCVforPCA$water.year)[j]),1]),as.vector(scores.concCV$sites[which(annualconcCVforPCA$Site==unique(annualconcCVforPCA$Site)[i]&annualconcCVforPCA$water.year==unique(annualconcCVforPCA$water.year)[j]),2]),pch=site.symbols[i],col=colors[j],bg=colors[j])
  }
}

legend('bottomleft',legend = unique(annualconcCVforPCA$Site),pch=site.symbols,pt.bg = "black",bty='n',cex=0.9)
gradient.rect(-1.85,-0.48,-1.4,-0.38,col=colors)

text(-1.85,-0.33,"2000",adj=c(0,0),cex=0.7)
text(-1.4,-0.33,"2018",adj=c(1,0),cex=0.7)


dev.off()



#### 7. Annual ratios of CV-chem vs. CV-flow
annualCVratioforPCA<-filter(annualCVratios,is.na(Cl.flow.CVratio)==F & is.na(NO3.flow.CVratio)==F & is.na(SO4.flow.CVratio)==F & is.na(TN.flow.CVratio)==F & is.na(TP.flow.CVratio)==F & is.na(PO4.flow.CVratio)==F)


CVratio.pca<-rda(annualCVratioforPCA[,seq(from=4,by=2,length.out = 6)], scale=T)
# PC1 explains 53.97% of the variation, PC2 explains 17.40% of the variation

scores.CVratio<-scores(CVratio.pca)
row.names(scores.CVratio$species)<-c("Cl","SO4","NO3","PO4","TN","TP")


## Biplot A: site by color with loadings

png("FIGURES/PCAonAnnualCVratio_A.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(min.pca,xlab="PC 1 (54.0%)",ylab="PC 2 (17.4%)",type="n",axes=F,xlim=c(-2.2,0.8),ylim=c(-0.8,2.2))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

arrows(0,0,scores.CVratio$species[,1],scores.CVratio$species[,2],length=0.1,angle=30)
Nudge<-1.05
text(scores.CVratio$species[,1]*Nudge,scores.CVratio$species[,2]*Nudge,rownames(scores.CVratio$species),cex=0.8,font=2)

# add sites as colors

for (i in 1:length(unique(annualCVratioforPCA$Site))){
  points(as.vector(scores.CVratio$sites[which(annualCVratioforPCA$Site==unique(annualCVratioforPCA$Site)[i]),1]),as.vector(scores.CVratio$sites[which(annualCVratioforPCA$Site==unique(annualCVratioforPCA$Site)[i]),2]),col=site.colors[i])
}
legend('topleft',legend = unique(annualCVratioforPCA$Site),col=site.colors,pch=1,pt.lwd=2)

dev.off()


## Biplot B: zoomed in with site by symbol and year by color

png("FIGURES/PCAonAnnualCVratio_B.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(min.pca,xlab="PC 1 (54.0%)",ylab="PC 2 (17.4%)",type="n",axes=F,xlim=c(-1.75,0.5),ylim=c(-0.75,2))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

# add points with site indicated by symbol and year by color (viridis color palette for years is names "colors")

for (i in 1:length(unique(annualCVratioforPCA$Site))){
  for (j in 1:length(unique(annualCVratioforPCA$water.year))){
    points(as.vector(scores.CVratio$sites[which(annualCVratioforPCA$Site==unique(annualCVratioforPCA$Site)[i]&annualCVratioforPCA$water.year==unique(annualCVratioforPCA$water.year)[j]),1]),as.vector(scores.CVratio$sites[which(annualCVratioforPCA$Site==unique(annualCVratioforPCA$Site)[i]&annualCVratioforPCA$water.year==unique(annualCVratioforPCA$water.year)[j]),2]),pch=site.symbols[i],col=colors[j],bg=colors[j])
  }
}

legend('topleft',legend = unique(annualCVratioforPCA$Site),pch=site.symbols,pt.bg = "black",bty='n',cex=0.9)
gradient.rect(-2.05,1.15,-1.55,1.05,col=colors)

text(-2.05,0.96,"2000",adj=c(0,0),cex=0.7)
text(-1.55,0.96,"2018",adj=c(1,0),cex=0.7)


dev.off()





################################
###   AJ's ordination code   ###
################################

### Ordination on annual means

mean.pca<-rda(means.lulc[,3:8], scale=T)
summary(mean.pca)
permanova.mean<-adonis2(means.lulc[,3:8]~means.lulc$water.year*means.lulc$Site)
permanova.mean
mean.pca.scores<-scores(mean.pca)


# Make a PCA plot
pca.plot<-ordiplot(mean.pca, type='n', xlab='PC 1 (55.9%)', ylab='PC 2 (31.6%)', xlim=c(-2,1), ylim=c(-2, 1))
points(pca.plot, "sites", pch=19, col="#fd8d3c", cex=1)
ordihull(mean.pca, groups=means.lulc$site, col="#fd8d3c", draw='polygon')
ordihull(mean.pca, groups=means.lulc$site, col="#fd8d3c", draw='lines')

orditorp(mean.pca, "species", col="#08519c", cex=1, air=0.1)
#text(pca.plot$sites[,1], pca.plot$sites[,2], means$site)#this tells you which site is represented by each point but it's crazy busy so don't use on final plot
text(0.45,-0.7, "DRKR")
text(0.8,0.35, "POBR")
text(0.6, 0.5, "BARN")
text(-0.2, 0.35, "GFGB")
text(0.5,0.1, "GFVN")
text(-0.5, -0.15,"GFGL")
text(0.25,-0.4, "GFCP")
text(-1.2, 0.7, "MCDN")

windows()
plot(means.lulc$isc, mean.pca.scores$sites[,1], pch=16, col='gray50', xlab="ISC (%)", ylab="PC 1 Score", cex.lab=1.5)
summary(lm(mean.pca.scores$sites[,1]~means.lulc$isc))
text(35, -1, "p>0.1", cex=1.5)

plot(means.lulc$isc, mean.pca.scores$sites[,2])
summary(lm(mean.pca.scores$sites[,2]~means.lulc$isc))

plot(means.lulc$tot.ag, mean.pca.scores$sites[,1], pch=16, col='gray50', xlab="Ag in watershed (proportion)", ylab="PC 1 Score", cex.lab=1.5)
model1<-(lm(mean.pca.scores$sites[,1]~means.lulc$tot.ag))
clip(min(means.lulc$tot.ag), max(means.lulc$tot.ag), min(mean.pca.scores$sites[,1]), max(mean.pca.scores$sites[,1]))
summary(model1)
abline(model1, col='red', lwd=2)
text(0.6, 0.5, "p<0.001")
text(0.6, 0.4, "r2 = 0.70")
text(0.75,-0.5, "MCDN", xpd=T)

plot(means.lulc$tot.urb, mean.pca.scores$sites[,2], pch=16, col='gray50', xlab="Urban in watershed (proportion)", ylab="PC 2 Score", cex.lab=1.5)
model1<-(lm(mean.pca.scores$sites[,2]~means.lulc$tot.urb))
clip(min(means.lulc$tot.urb), max(means.lulc$tot.urb), min(mean.pca.scores$sites[,2]), max(mean.pca.scores$sites[,2]))
summary(model1)
abline(model1, col='red', lwd=2)

text(0.1, -1, "p<0.001")
text(0.1, -1.1, "r2 = 0.53")

plot(means.lulc$isc, mean.pca.scores$sites[,2], pch=16, col='gray50', xlab="ISC (%)", ylab="PC 2 Score", cex.lab=1.5)
model1<-(lm(mean.pca.scores$sites[,2]~means.lulc$isc))
clip(min(means.lulc$isc), max(means.lulc$isc), min(mean.pca.scores$sites[,2]), max(mean.pca.scores$sites[,2]))
summary(model1)
abline(model1, col='red', lwd=2)

text(5, -1, "p<0.001")
text(5, -1.1, "r2 = 0.78")


#Now if we drop MCDN

means.nomcdn<-means.lulc[means.lulc$site!='MCDN',]

nomcdn.mean.pca<-rda(means.nomcdn[,3:8], scale=T)
summary(nomcdn.mean.pca)
permanova.nomcdn.mean<-adonis2(means.nomcdn[,3:8]~means.nomcdn$water.year*means.nomcdn$site)
permanova.nomcdn.mean
mean.nomcdn.pca.scores<-scores(nomcdn.mean.pca)

windows()

nomcdn.pca.plot<-ordiplot(nomcdn.mean.pca, type='n', xlab='PC 1 (49.7%)', ylab='PC 2 (29.7%)')
points(nomcdn.pca.plot, "sites", pch=19, col="#fd8d3c", cex=1)
ordihull(nomcdn.mean.pca, groups=means.nomcdn$site, col="#fd8d3c", draw='polygon')
ordihull(nomcdn.mean.pca, groups=means.nomcdn$site, col="#fd8d3c", draw='lines')

orditorp(nomcdn.mean.pca, "species", col="#08519c", cex=1, air=0.1)
#text(nomcdn.pca.plot$sites[,1], nomcdn.pca.plot$sites[,2], means.nomcdn$site)#this tells you which site is represented by each point but it's crazy busy so don't use on final plot
text(1.0,-0.7, "DRKR")
text(-.75,-.6, "POBR")
text(-0.5, 0.5, "BARN")
text(0.2, 1.3, "GFGB")
text(-0.4,0, "GFVN")
text(0.55, 0.55,"GFGL")
text(0.1,-0.3, "GFCP")

windows()
plot(means.nomcdn$isc, mean.nomcdn.pca.scores$sites[,1], pch=16, col='gray50', xlab="ISC (%)", ylab="PC 1 Score", cex.lab=1.5)
model1<-(lm(mean.nomcdn.pca.scores$sites[,1]~means.nomcdn$isc))
clip(min(means.nomcdn$isc), max(means.nomcdn$isc), min(mean.nomcdn.pca.scores$sites[,1]), max(mean.nomcdn.pca.scores$sites[,1]))
summary(model1)
abline(model1, col='red', lwd=2)

text(5,1,"p<0.001")
text(5,0.9, "r^2=0.55")

plot(means.nomcdn$tot.urb, mean.nomcdn.pca.scores$sites[,1], pch=16, col='gray50', xlab="Urban in watershed (proportion)", ylab="PC 1 Score", cex.lab=1.5)
model1<-(lm(mean.nomcdn.pca.scores$sites[,1]~means.nomcdn$tot.urb))
clip(min(means.nomcdn$tot.urb), max(means.nomcdn$tot.urb), min(mean.nomcdn.pca.scores$sites[,1]), max(mean.nomcdn.pca.scores$sites[,1]))
summary(model1)
abline(model1, col='red', lwd=2)

text(0.1,1,"p<0.001")
text(0.1,0.9, "r^2=0.67")

plot(means.nomcdn$isc, mean.nomcdn.pca.scores$sites[,2], pch=16, col='gray50', xlab="ISC (%)", ylab="PC 1 Score", cex.lab=1.5)
model1<-(lm(mean.nomcdn.pca.scores$sites[,2]~means.nomcdn$isc))
clip(min(means.nomcdn$isc), max(means.nomcdn$isc), min(mean.nomcdn.pca.scores$sites[,2]), max(mean.nomcdn.pca.scores$sites[,2]))
summary(model1)
abline(model1, col='red', lwd=2)

text(5,1.2,"p=0.006")
text(5,1.1, "r^2=0.47")
text(0,-0.7,"POBR", xpd=T)
text(1,0.6,"BARN", xpd=T)

plot(means.nomcdn$tot.urb, mean.nomcdn.pca.scores$sites[,2])
summary(lm(mean.nomcdn.pca.scores$sites[,2]~means.nomcdn$tot.urb))

plot(means.nomcdn$tot.ag, mean.nomcdn.pca.scores$sites[,1])
summary(lm(mean.nomcdn.pca.scores$sites[,1]~means.nomcdn$tot.ag))

meanNMDS<-metaMDS(summary.data[,c(3,6,9,12,15,18)], k=2, trymax=100)
windows(height=5, width=5)
par(mar=c(4,4,0.2,0.2))
nmds.plot<-ordiplot(meanNMDS, type='n')
ordihull(meanNMDS, groups=summary.data$site, col='#fdae6b', draw='polygon')
ordihull(meanNMDS, groups=summary.data$site, col='#fdae6b', draw='lines')
points(nmds.plot, "sites", pch=19, col="#fd8d3c", cex=1)
orditorp(meanNMDS, "species", col='#08519c', cex=1, air=0.1)
text(0.38, 0.38, "POBR")
text(-0.18,0.27, "MCDN")
text(-0.25, 0.1, "BARN")
text(-0.2, -0.15, "GFGB")
text(0.30, -0.25, "DRKR")
text(0.13,0.1,"GFCP")
text(-0.1, -0.2, "GFVN")
text(0.02, 0.03, "GFGL")
windows(height=5,width=5)
stressplot(meanNMDS)

#Run a PERMANOVA on the data
mean.permanova<-adonis2(summary.data[,c(3,6,9,12,15,18)]~summary.data$site*summary.data$water.year)

#Repeat with ranges across sites (at least make the plots...):
windows(height=4, width=6)
par(mfrow=c(2,4), mar=c(1,1,1,1))
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$Site=='POBR',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$Site=='BARN',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$Site=='MCDN',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$Site=='DRKR',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$Site=='GFGL',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$Site=='GFGB',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$Site=='GFVN',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$Site=='GFCP',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)

