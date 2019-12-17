### R script for multivariate water chem analysis of BES - adapted from previous NMDS/Radar chart scripts
### Created by AJR
### Created on 19 June 2017
### Last edited by: MLF
### Last edited on: 17 December 2019
### Last edits: Revised how we calculate means, max and mins using ddply, made radar charts, ran PCAs, made figures to put into powerpoint for Peter and Emma
### Next steps: 

#First set your working directory - this should be changed to reflect wherever you have saved the .csv files
#Note that if you use a project (here 'BaltimoreStreamChem' the directory will be set to the project directory)
# setwd("C:/Users/reisingera/Dropbox/Documents/R/Balt-ChemFingerprints")

library(plyr)
library(fmsb)
library(car)
library(vegan)
library(animation)

#Do some data cleanup from the rawdata emailed by Lisa on 13 Dec 2019 (data from 15-Oct-98 through 7-Jun-19):
bes.waterchem.raw<-read.csv(file="BESFullWaterChemRaw.csv", header=T)

bes.waterchem.working<-bes.waterchem.raw[bes.waterchem.raw$site=='GFCP'|bes.waterchem.raw$site=='GFVN'|bes.waterchem.raw$site=='GFGB'|bes.waterchem.raw$site=='GFGL'|
    bes.waterchem.raw$site=='DRKR'|bes.waterchem.raw$site=='POBR'|bes.waterchem.raw$site=='BARN'|bes.waterchem.raw$site=='MCDN',]

bes.waterchem.working<-bes.waterchem.working[,c(1:10,13:15)]

write.csv(bes.waterchem.working, file="bes.waterchem.csv")

#Read in the working data
lulc.data<-read.csv(file="beslulc.csv", header=T)#Updated BES LULC from Claire Welty (see email sent from Claire on 11 Jan 2017)

### Look at the data and then clean up in excel
waterchem.data<-read.csv(file="bes.waterchem.csv", header=T)


#Make a dataframe with the maximum values from the entire range of sampling
max<-as.numeric(apply(waterchem.data, 2, max, na.rm=T))
max # looks like we have some bad data points for temp, do, and ph - will clean up in excel

data<-read.csv(file="bes.waterchem.clean.csv", header=T)

#add a water year column
data$water.year<-rep(NA, length(data[,1]))

#fill the water year. If the julian date is before 275 (Oct 2 I believe, or Oct 1 for leap year...probably could do this better but it works) then the 
#water year for that row is the same as the actual year, if the julian date is >= 275 then the water year assigned is the next calendar year
for(i in 1:length(data$water.year)){
  if(data[i,3]<275){
    data[i,14]=data[i,2]
  }
  else{
    (data[i,14]=data[i,2]+1)
  }
}

###Make a color pallette for 15 classes (15 years of data)
colors<-c("gray20","gray40", "gray60", "gray80", rgb(19,149,186,maxColorValue=255),rgb(17,120,153,maxColorValue=255),rgb(17,120,153,maxColorValue=255),
          rgb(15,91,120,maxColorValue=255),rgb(192,46,29,maxColorValue=255), rgb(217,78,31,maxColorValue=255), rgb(241,108,32,maxColorValue=255), 
          rgb(239,139,44,maxColorValue=255),rgb(236,170,56,maxColorValue=255),rgb(235,200,68,maxColorValue=255),rgb(162,184,108,maxColorValue=255))


#calculate the mean, max, and mean for each variable of interest by water year and site

summary.data<-ddply(data[,c(4:10,14)], .(water.year, site), summarize, 
                    mean.chloride=mean(chloride, na.rm=T), max.chloride=max(chloride, na.rm=T), min.chloride=min(chloride, na.rm=T),
                    mean.no3=mean(no3, na.rm=T), max.no3=max(no3, na.rm=T), min.no3=min(no3, na.rm=T),
                    mean.po4=mean(po4, na.rm=T), max.po4=max(po4, na.rm=T), min.po4=min(po4, na.rm=T),
                    mean.so4=mean(so4, na.rm=T), max.so4=max(so4, na.rm=T), min.so4=min(so4, na.rm=T),
                    mean.tn=mean(tn, na.rm=T), max.tn=max(tn, na.rm=T), min.tn=min(tn, na.rm=T),
                    mean.tp=mean(tp, na.rm=T), max.tp=max(tp, na.rm=T), min.tp=min(tp, na.rm=T))
#Remove 1999 and 2015 from the records as we don't have a full record of those years
summary.data<-summary.data[summary.data$water.year!=1999&summary.data$water.year!=2015,]

range.data<-ddply(data[,c(4:10,14)], .(water.year, site), summarize,
                  chloride.range<-max(chloride, na.rm=T)-min(chloride, na.rm=T),
                  no3.range<-max(no3, na.rm=T)-min(no3, na.rm=T),
                  po4.range<-max(po4, na.rm=T)-min(po4, na.rm=T),
                  so4.range<-max(so4, na.rm=T)-min(so4, na.rm=T),
                  tn.range<-max(tn, na.rm=T)-min(tn, na.rm=T),
                  tp.range<-max(tp, na.rm=T)-min(tp, na.rm=T))
range.data<-range.data[range.data$water.year!=1999&range.data$water.year!=2015,]

colnames(range.data)<-c("water.year", "site", "chloride", "no3", "po4", "so4", "tn", "tp")

#Make radar charts showing the annual mean values for each site:
windows(height=3, width=6)
par(mfrow=c(2,4), mar=c(0.6,1,0.6,1))
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$site=='POBR',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='POBR')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$site=='BARN',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='BARN')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$site=='MCDN',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='MCDN')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$site=='DRKR',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='DRKR')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$site=='GFGL',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='GFGL')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$site=='GFGB',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='GFGB')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                  summary.data[summary.data$site=='GFVN',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                  plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                  expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='GFVN')
radarchart(rbind(apply(summary.data[,c(3,6,9,12,15,18)], 2, max, na.rm=T), apply(summary.data[,c(3,6,9,12,15,18)], 2, min, na.rm=T), 
                 summary.data[summary.data$site=='GFCP',c(3,6,9,12,15,18)]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
                    plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                    expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors, title='GFCP')

means<-summary.data[,c(1:3,6,9,12,15,18)]

means.lulc<-merge(means, lulc.data)


windows()

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
                 range.data[range.data$site=='POBR',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$site=='BARN',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$site=='MCDN',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$site=='DRKR',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$site=='GFGL',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$site=='GFGB',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$site=='GFVN',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)
radarchart(rbind(apply(range.data[,3:8], 2, max, na.rm=T), apply(range.data[,3:8], 2, min, na.rm=T), 
                 range.data[range.data$site=='GFCP',3:8]), axistype=0, axislabcol='gray20', palcex=1.2, vlcex=1, plty=1, 
           plwd=2, pty=16, cglcol='gray20', maxmin=T, vlabels=c(expression(paste("Cl"^" -")), expression(paste("NO"[3]^" -")), 
                                                                expression(paste("PO"[4]^" 3-")), expression(paste("SO"[4]^" 2-")), "TN", "TP"), pcol=colors)

