### R script for multivariate water chem analysis of BES - adapted from previous NMDS/Radar chart scripts
### Created by AJR
### Created on 19 June 2017
### Last edited by: MLF
### Last edited on: 17 December 2019
### Last edits: Updated calculattions for new data (WY2000-WY2018), automated removal of bad data points, updated to dplyr functions, added download of NWIS data

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

write.csv(waterchem.clean, file="bes.waterchem.clean.csv",row.names = F) #write out the clean data



#Read in the clean working data for analysis
lulc.data<-read.csv(file="beslulc.csv", header=T)#Updated BES LULC from Claire Welty (see email sent from Claire on 11 Jan 2017)

data<-read.csv(file="bes.waterchem.clean.csv", header=T)


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

###Make a color pallette for 19 classes (19 years of data)
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
  "POBR",'01583570',
  "BARN",'01583580',
  "GFGB",'01589197',
  "GFVN",'01589300',
  "GFCP",'01589352',
  "DRKR",'01589330',
  "GFGL",'01589180'),ncol=2,byrow=T))
colnames(USGS.site.nos)<-c("siteName","siteNumber")
Q.start<-"1998-10-01"
Q.end<-"2019-09-30"


BESflow<-readNWISdv(siteNumbers = USGS.site.nos$siteNumber,parameterCd = '00060',startDate = Q.start,endDate = Q.end)



# Multisolute Analysis Plan:
# Run PCAs for the following annual summaries of BES stream chem: annual mean, annual range, annual min, annual max, annual C-Q slope, annual CV of chem:annual CV of flow
# plot each of these using viridis color ramp to differentiate among years (show progression through time, maybe add arrows between successive years?)
# Perhaps also run PCAs for all the points within each site, and compare hull areas or centroids for each year to get a sense of how sampling "snapshots" of water chem are changing
# think about how to compare reactive and conservative solutes (CV ratio?)
# think about some way to summarize multi-solute exceeedence probabilities? 
# Changing flow thresholds for N going from reactive to conservative behavior
# Adding pharma data when we get the data back from Jerker?

####################################################
###              Annual Ordinations              ###
####################################################

### 1. Annual means 
mean.pca<-rda(means.lulc[,3:8], scale=T)
# PC1 explains 56.45% of the variation, PC2 explains 31.06% of the variation

scores.mean<-scores(mean.pca,scaling=2)

## Biplot

# png("FIGURES/PCAonAnnualMeans",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(mean.pca,xlab="PC 1 (56.5%)",ylab="PC 2 (31.1%)",type="n",axes=F,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()


arrows(0,0,scores.mean$species[,1],scores.mean$species[,2],length=0.1,angle=30)
Nudge<-1.1
text(scores.mean$species[,1]*Nudge,scores.mean$species[,2]*Nudge,rownames(scores.mean$species),cex=0.8,font=2)



### Fit Initial DOC concentration and proportional DOC loss to the PCA as environmental variables ###
#First, make a dataframe for the DOC variables
DOC.info<-cbind(paste0(bio$Sample.location,'_',bio$Experiment),bio[,c(6,8:11)])
colnames(DOC.info)<-c('name_season',colnames(bio[,c(6,8:11)]))

foo<-as.data.frame(paste0(bio.points.ordered[,2],'_',bio.points.ordered[,5]))
colnames(foo)<-'name_season'
DOC.env<-merge(DOC.info,foo,by="name_season")
#Standardize the data:
sdzd.DOC.env<-matrix(nrow=nrow(DOC.env),ncol=ncol(DOC.env))
for (i in 2:ncol(DOC.env)){
  sdzd.DOC.env[,i]<-(DOC.env[,i]-mean(DOC.env[,i]))/sd(DOC.env[,i])
}
sdzd.DOC.env[,1]<-DOC.env[,1]

DOC.fit<-envfit(pca2,sdzd.DOC.env[,2:6])

scores.2<-scores(pca2,scaling=2)
rownames(scores.2$species)<-colnames(bio.pca2)
foo<-scores(pca2,display="species",scaling=2)
DOC.scores<-scores(DOC.fit,"vectors",scaling=2) #this gives the scaled locations of the arrow heads


# PCA plots

### PANEL A: PCA Biplot with loadings ###

png("D:/DISSERTATION/Urban Streams/BIOAVAILABILITY/BIOAVAIL FIGURES/ES&T figures/OM_PCA_Loadings.png",height=5,width=5,units='in',res=300)
par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(pca2,xlab="PC 1 (47.9%)",ylab="PC 2 (20.4%)",type="n",axes=F,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

rownames(scores.2$species)<-c("C1","C2","C3","HIX","BIX","SUVA","FI")

arrows(0,0,scores.2$species[,1],scores.2$species[,2],length=0.1,angle=30)
Nudge<-1.1
text(scores.2$species[,1]*Nudge,scores.2$species[,2]*Nudge,rownames(scores.2$species),cex=0.8,font=2)


#Add points based on their type (symbol) and incubation time (color)
## "Stream Flows" = blue
### Daylighted Upstream Channel = U
points(scores.2$sites[c(88,94,184,190),1],scores.2$sites[c(88,94,184,190),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(89,95,185,191),1],scores.2$sites[c(89,95,185,191),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(90,96,186,192),1],scores.2$sites[c(90,96,186,192),2],cex=0.9,pch=16,col="#deebf7")
### ISCO storm samples = S
points(scores.2$sites[c(34,37,40,43,46,49,52,55,58,97,136,139,142,145,148,151,154,157,160,163),1],scores.2$sites[c(34,37,40,43,46,49,52,55,58,97,136,139,142,145,148,151,154,157,160,163),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(35,38,41,44,47,50,53,56,59,98,137,140,143,146,149,152,155,158,161,164),1],scores.2$sites[c(35,38,41,44,47,50,53,56,59,98,137,140,143,146,149,152,155,158,161,164),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(36,39,42,45,48,51,54,57,60,99,138,141,144,147,150,153,156,159,162,165),1],scores.2$sites[c(6,39,42,45,48,51,54,57,60,99,138,141,144,147,150,153,156,159,162,165),2],cex=0.9,pch=16,col="#deebf7")

## "Primary sources" = brown
### Biofilm = b
points(scores.2$sites[c(97),1],scores.2$sites[c(97),2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[c(98),1],scores.2$sites[c(98),2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[c(99),1],scores.2$sites[c(99),2],cex=0.9,pch=16,col="#fce2cc")
### Grass = g
points(scores.2$sites[115,1],scores.2$sites[115,2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[116,1],scores.2$sites[116,2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[117,1],scores.2$sites[117,2],cex=0.9,pch=16,col="#fce2cc")
### Leaf (green, senesced, pine needles) = l
points(scores.2$sites[c(61,67,118),1],scores.2$sites[c(61,67,118),2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[c(62,68,119),1],scores.2$sites[c(62,68,119),2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[c(63,69,120),1],scores.2$sites[c(63,69,120),2],cex=0.9,pch=16,col="#fce2cc")
### Soil = s
points(scores.2$sites[85,1],scores.2$sites[85,2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[86,1],scores.2$sites[86,2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[87,1],scores.2$sites[87,2],cex=0.9,pch=16,col="#fce2cc")
### Mulch = m
points(scores.2$sites[64,1],scores.2$sites[64,2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[65,1],scores.2$sites[65,2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[66,1],scores.2$sites[66,2],cex=0.9,pch=16,col="#fce2cc")

## "Leachates" = brown
### CB = C
points(scores.2$sites[c(7,10,13,16,103,106,109,112),1],scores.2$sites[c(7,10,13,16,103,106,109,112),2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[c(8,11,14,17,104,107,110,113),1],scores.2$sites[c(8,11,14,17,104,107,110,113),2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[c(9,12,15,18,105,108,111,114),1],scores.2$sites[c(9,12,15,18,105,108,111,114),2],cex=0.9,pch=16,col="#fce2cc")
### gutter = G
points(scores.2$sites[c(28,31,130,133),1],scores.2$sites[c(28,31,130,133),2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[c(29,32,131,134),1],scores.2$sites[c(29,32,131,134),2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[c(30,33,132,135),1],scores.2$sites[c(30,33,132,135),2],cex=0.9,pch=16,col="#fce2cc")

## "Engineered Headwater Flows" = blue
### CB = C
points(scores.2$sites[c(1,4),1],scores.2$sites[c(1,4),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(2,5),1],scores.2$sites[c(2,5),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(3,6),1],scores.2$sites[c(3,6),2],cex=0.9,pch=16,col="#deebf7")
### gutter = G
points(scores.2$sites[c(19,22,25,121,124,127),1],scores.2$sites[c(19,22,25,121,124,127),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(20,23,26,122,125,128),1],scores.2$sites[c(20,23,26,122,125,128),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(21,24,27,123,126,129),1],scores.2$sites[c(21,24,27,123,126,129),2],cex=0.9,pch=16,col="#deebf7")
### Road runoff = R
points(scores.2$sites[c(70,73,76,166,169,172),1],scores.2$sites[c(70,73,76,166,169,172),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(71,74,77,167,170,173),1],scores.2$sites[c(71,74,77,167,170,173),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(72,75,78,168,171,174),1],scores.2$sites[c(72,75,78,168,171,174),2],cex=0.9,pch=16,col="#deebf7")
### Storm Drain = D
points(scores.2$sites[c(79,82,175,178,181),1],scores.2$sites[c(79,82,175,178,181),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(80,83,176,179,182),1],scores.2$sites[c(80,83,176,179,182),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(81,84,177,180,183),1],scores.2$sites[c(81,84,177,180,183),2],cex=0.9,pch=16,col="#deebf7")

rownames(DOC.scores)<-c("Initial [DOC]","Initial DOC Mass Loss","Initial Proportion DOC Lost","Total DOC Mass Loss","Total Proportion DOC Lost")
arrows(0,0,DOC.scores[c(1,5),1],DOC.scores[c(1,5),2],length = 0.07,angle=30,col="gray50",lwd=2)
text(c(0.5,-0.6),c(0.25,-0.21),c("Proportion DOC Lost","Initial [DOC]"),col="gray50",cex=0.7)

dev.off()

#####################################
#  Panel B - trajectories 0-6  days #
#####################################
png("D:/DISSERTATION/Urban Streams/BIOAVAILABILITY/BIOAVAIL FIGURES/ES&T figures/OM.PCA.BiplotTrajectories0.6.png",height=5,width=5,units='in',res=300)

par(mar=c(3,3,0.2,0.2))
par(mgp=c(1.5,0.4,0))
plot(scores.2$sites[,1],scores.2$sites[,2],xlab="PC 1 (47.9%)",ylab="PC 2 (20.4%)",type="n",axes=F,xlim=c(-2,1))
axis(1,tck=0.02)
axis(2,tck=0.02)
box()

#Add points based on their type (symbol) and incubation time (color)
## "Stream Flows" = blue
### Daylighted Upstream Channel = U
points(scores.2$sites[c(88,94,184,190),1],scores.2$sites[c(88,94,184,190),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(89,95,185,191),1],scores.2$sites[c(89,95,185,191),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(90,96,186,192),1],scores.2$sites[c(90,96,186,192),2],cex=0.9,pch=16,col="#deebf7")
### ISCO storm samples = S
points(scores.2$sites[c(34,37,40,43,46,49,52,55,58,97,136,139,142,145,148,151,154,157,160,163),1],scores.2$sites[c(34,37,40,43,46,49,52,55,58,97,136,139,142,145,148,151,154,157,160,163),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(35,38,41,44,47,50,53,56,59,98,137,140,143,146,149,152,155,158,161,164),1],scores.2$sites[c(35,38,41,44,47,50,53,56,59,98,137,140,143,146,149,152,155,158,161,164),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(36,39,42,45,48,51,54,57,60,99,138,141,144,147,150,153,156,159,162,165),1],scores.2$sites[c(6,39,42,45,48,51,54,57,60,99,138,141,144,147,150,153,156,159,162,165),2],cex=0.9,pch=16,col="#deebf7")

## "Primary sources" = brown
### Biofilm = b
points(scores.2$sites[c(97),1],scores.2$sites[c(97),2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[c(98),1],scores.2$sites[c(98),2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[c(99),1],scores.2$sites[c(99),2],cex=0.9,pch=16,col="#fce2cc")
### Grass = g
points(scores.2$sites[115,1],scores.2$sites[115,2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[116,1],scores.2$sites[116,2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[117,1],scores.2$sites[117,2],cex=0.9,pch=16,col="#fce2cc")
### Leaf (green, senesced, pine needles) = l
points(scores.2$sites[c(61,67,118),1],scores.2$sites[c(61,67,118),2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[c(62,68,119),1],scores.2$sites[c(62,68,119),2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[c(63,69,120),1],scores.2$sites[c(63,69,120),2],cex=0.9,pch=16,col="#fce2cc")
### Soil = s
points(scores.2$sites[85,1],scores.2$sites[85,2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[86,1],scores.2$sites[86,2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[87,1],scores.2$sites[87,2],cex=0.9,pch=16,col="#fce2cc")
### Mulch = m
points(scores.2$sites[64,1],scores.2$sites[64,2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[65,1],scores.2$sites[65,2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[66,1],scores.2$sites[66,2],cex=0.9,pch=16,col="#fce2cc")

## "Leachates" = brown
### CB = C
points(scores.2$sites[c(7,10,13,16,103,106,109,112),1],scores.2$sites[c(7,10,13,16,103,106,109,112),2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[c(8,11,14,17,104,107,110,113),1],scores.2$sites[c(8,11,14,17,104,107,110,113),2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[c(9,12,15,18,105,108,111,114),1],scores.2$sites[c(9,12,15,18,105,108,111,114),2],cex=0.9,pch=16,col="#fce2cc")
### gutter = G
points(scores.2$sites[c(28,31,130,133),1],scores.2$sites[c(28,31,130,133),2],cex=0.9,pch=16,col="#bc5200")
points(scores.2$sites[c(29,32,131,134),1],scores.2$sites[c(29,32,131,134),2],cex=0.9,pch=16,col="#e5a067")
points(scores.2$sites[c(30,33,132,135),1],scores.2$sites[c(30,33,132,135),2],cex=0.9,pch=16,col="#fce2cc")

## "Engineered Headwater Flows" = blue
### CB = C
points(scores.2$sites[c(1,4),1],scores.2$sites[c(1,4),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(2,5),1],scores.2$sites[c(2,5),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(3,6),1],scores.2$sites[c(3,6),2],cex=0.9,pch=16,col="#deebf7")
### gutter = G
points(scores.2$sites[c(19,22,25,121,124,127),1],scores.2$sites[c(19,22,25,121,124,127),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(20,23,26,122,125,128),1],scores.2$sites[c(20,23,26,122,125,128),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(21,24,27,123,126,129),1],scores.2$sites[c(21,24,27,123,126,129),2],cex=0.9,pch=16,col="#deebf7")
### Road runoff = R
points(scores.2$sites[c(70,73,76,166,169,172),1],scores.2$sites[c(70,73,76,166,169,172),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(71,74,77,167,170,173),1],scores.2$sites[c(71,74,77,167,170,173),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(72,75,78,168,171,174),1],scores.2$sites[c(72,75,78,168,171,174),2],cex=0.9,pch=16,col="#deebf7")
### Storm Drain = D
points(scores.2$sites[c(79,82,175,178,181),1],scores.2$sites[c(79,82,175,178,181),2],cex=0.9,pch=16,col="#3182bd")
points(scores.2$sites[c(80,83,176,179,182),1],scores.2$sites[c(80,83,176,179,182),2],cex=0.9,pch=16,col="#9ecae1")
points(scores.2$sites[c(81,84,177,180,183),1],scores.2$sites[c(81,84,177,180,183),2],cex=0.9,pch=16,col="#deebf7")

# arrows for stream flow
for (i in 1:nrow(S.pca2)){
  arrows(S.pca2[i,1],S.pca2[i,2],S.pca2[i,3],S.pca2[i,4],length=0.075,lwd=1,col="#3182bd")
}
# for (i in 1:nrow(S.pca2)){
# arrows(S.pca2[i,3],S.pca2[i,4],S.pca2[i,5],S.pca2[i,6],length=0.075,lwd=1,col="#BB8FCE")
# }

# arrows for EH flow
for (i in 1:nrow(EH.pca2)){
  arrows(EH.pca2[i,1],EH.pca2[i,2],EH.pca2[i,3],EH.pca2[i,4],length=0.075,lwd=1,col="#3182bd")
}
# for (i in 1:nrow(EH.pca2)){
# arrows(EH.pca2[i,3],EH.pca2[i,4],EH.pca2[i,5],EH.pca2[i,6],length=0.075,lwd=1,col="#ABB2B9")
# }

# arrows for EH leachates
for (i in 1:nrow(L.pca2)){
  arrows(L.pca2[i,1],L.pca2[i,2],L.pca2[i,3],L.pca2[i,4],length=0.075,lwd=1,col="#bc5200")
}
# for (i in 1:nrow(L.pca2)){
# arrows(L.pca2[i,3],L.pca2[i,4],L.pca2[i,5],L.pca2[i,6],length=0.075,lwd=1,col="#E6B0AA")
# }

# arrows for primary source leachates
for (i in 1:nrow(P.pca2)){
  arrows(P.pca2[i,1],P.pca2[i,2],P.pca2[i,3],P.pca2[i,4],length=0.075,lwd=1,col="#bc5200")
}
# for (i in 1:nrow(P.pca2)){
# arrows(P.pca2[i,3],P.pca2[i,4],P.pca2[i,5],P.pca2[i,6],length=0.075,lwd=1,col="#A9DFBF")
# }

dev.off()















################################
###   AJ's ordination code   ###
################################

### Ordination on annual means

mean.pca<-rda(means.lulc[,3:8], scale=T)
summary(mean.pca)
permanova.mean<-adonis2(means.lulc[,3:8]~means.lulc$water.year*means.lulc$site)
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

