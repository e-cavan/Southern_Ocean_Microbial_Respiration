#############################################################
# R code to reproduce all statistical analyses and plots from the main manuscript of Cavan & Boyd.

rm(list=ls())

# Set working directory
setwd('~/path/to/folder')


## You will need to download the satellite chl file from https://oceancolor.gsfc.nasa.gov/cgi/l3
# Select Aqua MODIS Chl at 4km resolution

# Read in data
# open in situ pump POC depth profile
isp <- read.csv('/Users/elcavan/Documents/Papers/SOTS/Rev 2/Cavan_Boyd_2018_POC_depth.csv', header=T)
# Read in binned CTD data for salinty and temperature
ctd <- read.csv('/Users/elcavan/Documents/Papers/SOTS/Rev 2/Cavan_Boyd_2018_CTD.csv', header=T)
# open O2 file
resp <- read.csv('/Users/elcavan/Documents/Papers/SOTS/Rev 2/Cavan_Boyd_2018_O2_rates.csv')
# open POC vials file
poc <- read.csv('/Users/elcavan/Documents/Papers/SOTS/Rev 2/Cavan_Boyd_2018_POC_vials.csv', header=T)


# function to calculate standard error of the mean
std <- function(x) sd(x) / sqrt(length(x))

#setwd('/Users/elcavan/Documents/Papers/SOTS/Rev 2/')

# load libraries
library("ncdf4")
library("reshape2")
library(squash)
library(maps)
library(mapdata)
library(mapproj)
library(maptools)
library(lme4)
library(viridis)
library(sp)
library(rgdal)
library(raster)
library(rasterVis)

#***********************************************************
# Plot satellite chl map (Fig. 2a)
#***********************************************************

# Script for extracting data from NetCDF files adapted from Antonio Olinto Avila-da-Silva, Instituto de Pesca, Brasil
# https://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=5954


# test to see if file (.csv) has already been created with satellite data
# (I have included one in the data repository (Cavan_MODISA_chl_old.csv), if you struggle to manipulate the netcdf 
        #file you can use this .csv to build the map)
file.exists("MODISA_chl.csv")     # caution new data will be appended to
# rename this file if it already exists
file.rename("MODISA_chl.csv","Cavan_MODISA_chl_old.csv")


# set the study area. THis is set for the sub-antarctic southern ocean, SW of Tasmania
lonmax<-150
lonmin<-135
latmax<--40
latmin<--50

# Using downloaded satellite Chl file...
# create a list of files and indicate its length - as only plotting one month (ie not average over year) then just one file
f <- list.files(".", pattern="*.L3m_MO_CHL_chlor_a_4km.nc",full.names=F)
lf<-length(f)

# variable of interest
var<-"chlor_a"


# open netCDF file              # can do this as a for loop if have more than one file to process
data<-nc_open(f[1])
# extract data
lon<-ncvar_get(data,"lon")
lat<-ncvar_get(data,"lat")
value<-ncvar_get(data,var)
unit<-ncatt_get(data,var,"units")$value
# matrix to data.frame
dimnames(value)<-list(lon=lon,lat=lat)
dat.var<-melt(value,id="lon")
# select data from the study area taking out missing data
dat.varSAtmp<-subset(dat.var,lon<=lonmax & lon>=lonmin & lat<=latmax &
                       lat>=latmin & value<45)
# extract date information
dateini<-ncatt_get(data,0,"time_coverage_start")$value
dateend<-ncatt_get(data,0,"time_coverage_end")$value
datemean<-mean(c(as.Date(dateend,"%Y-%m-%dT%H:%M:%OSZ"),as.Date(dateini,"%Y-%m-%dT%H:%M:%OSZ")))
year<-substring(datemean,0,4)
month<-substring(datemean,6,7)

# prepare final data set
dat.varSA<-data.frame(rep(as.integer(year,nrow(dat.varSAtmp))),rep(as.integer(month,nrow(dat.varSAtmp))),
                      dat.varSAtmp,rep(unit,nrow(dat.varSAtmp)),rep(var,nrow(dat.varSAtmp)))
names(dat.varSA)<-c("year","month","lon","lat","value","unit","var")
# save csv file
fe<-file.exists("MODISA_chl.csv")
write.table(dat.varSA,"MODISA_chl.csv",row.names=FALSE,col.names=!fe,sep=",",dec=".",append=fe)

# close connection
nc_close(data)





# Convert dataframe to raster for plotting
chl <- read.csv('MODISA_chl.csv', sep=',')

# subset data - lon, lat, value of 1 plot
chl_m<-chl[,c(3,4,5)]
# put value (chl) column first
chl_m <- subset(chl_m, select=c(3,1:2))
# set max chl level
chl_m <- chl_m[chl_m$value<2,]
# data fram of lon and lat
xy<-chl_m[,c(2,3)]

# convert to spatial data frame fram using value df and lon/lat df
spdf <- SpatialPointsDataFrame(coords = xy, data = chl_m,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# convert pixels due to irregular point data, tolerance given from gridden(pixels)=TRUE
pixels <- SpatialPixelsDataFrame(spdf, tolerance = 0.00448656, spdf@data)
# tell r it is gridded
gridded(pixels) = TRUE
# convert to raster using raster package
r<-raster(pixels)


# levelplot 
boundaries<-map("worldHires", ylim=c(-50,-40), xlim=c(135,150), col='gray90', fill=TRUE, myborder=0.01,
                boundary=FALSE)
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                             proj4string=CRS(projection(r)))
mapTheme <- rasterTheme(region=viridis(300))
b <- levelplot(r, contour=FALSE, labels=FALSE, margin=FALSE, par.settings=mapTheme,
               ylim=c(-50,-40), xlim=c(135,150), 
               main=expression(paste('Chl-a (',mu,'g L'^'-1'*')')),
               colorkey=list(at=seq(0, 2, 0.02), labels=list(at=c(0,1,2), 
                                                             labels=c('0','1','2'))))
par(mar=c(1,1,1,1), mfcol=c(2,2), bty='o', cex=1)
# co-ordinates of SOTS site (red point)
t <- 141.5
y<- -46
xy <- data.frame(t,y)
yx <- SpatialPoints(xy[,1:2])
# plot map with tasmania and SOTS site (red point)
b + layer(sp.polygons(bPols, fill='gray')) +layer(sp.points(yx, pch=16, col='red', cex=2), columns = 1)



## find chl at Southern Ocean time series site
#identify lat long
s <- rowColFromCell(r, cellFromXY(r, cbind(141.5,-46)))
#from answer use raster to find value
r[s[1], s[2]]






par(mar=c(5,5,1,1), mfcol=c(1,1),cex=1)

#***********************************************************
# Plot depth profiles of POC and temperature/salinity (Fig. 2b)
#***********************************************************
# In situ pump data

# split isps
isp1 <- isp[isp$Type=='ISP1',]
isp2 <- isp[isp$Type=='ISP2',]

# average across replicate punches of filter (or depth, n = 3)
isp1z <- aggregate(isp1$ug.L.C.1, by=list(isp1$depth), mean)
isp2z <- aggregate(isp2$ug.L.C.1, by=list(isp2$depth), mean)


# std error across replicate punches of filter
isp1z.err <- aggregate(isp1$ug.L.C.1, by=list(isp1$depth), std)
isp2z.err <- aggregate(isp2$ug.L.C.1, by=list(isp2$depth), std)

# rename
names(isp1z) <- c('depth', 'poc')
names(isp2z) <- c('depth', 'poc')



#plot with error bars
par(mar=c(6,5,5,1))
plot(isp1z$poc, isp1z$depth,  pch=16, ylim=c(500,0),
     ylab='Depth (m)', xlab='', xaxt='n',type='n', yaxt='n',
     cex.axis=1.5, cex.lab=1.5)
axis(3, at=seq(0,125,20), cex.axis=1.5)
axis(2, labels=FALSE)
mtext(expression(paste('POC (',mu,'g L'^'-1'*')')), 3, 2, cex=1.5)
text(y=seq(500, 0, -100),cex=1.5,
     par('usr')[1], pos=2,labels=as.vector(c(500 , 400 , 300 , 200 , 100 ,0 )),srt=0, xpd=TRUE)


# error bars isp1
xy.error.bars<- function (x,y,xbar){
  points(isp1z$poc, isp1z$depth,  pch=16, xaxt='n', xlab='',
         cex=1, xpd=FALSE)
  arrows(x-xbar, y, x+xbar, y, code=3, angle=90, length=0.05, lwd=1.2, xpd=FALSE) 
}
x<-isp1z$poc
y<-isp1z$depth
xb<-isp1z.err$x
xy.error.bars(x,y,xb)


# error bars isp2
xy.error.bars<- function (x,y,xbar){
  points(isp2z$poc, isp2z$depth,  pch=16, xaxt='n', xlab='',
         cex=1, xpd=FALSE)
  arrows(x-xbar, y, x+xbar, y, code=3, angle=90, length=0.05, lwd=1.2, xpd=FALSE) 
}
x<-isp2z$poc
y<-isp2z$depth
xb<-isp2z.err$x
xy.error.bars(x,y,xb)

# bind isp1 and isp2 in one dataframe
isp.all <- rbind(isp1z, isp2z)


# calculate z* (remineralisation length scale) for poc ug/l
z <- lm((isp.all$depth)~log(isp.all$poc)); summary(z) # z* = 114, p<0.01, r2=0.63

# plot z* lines
y <- seq(29, 500,2)
x <- 124.719221 * exp((y-0)/-113)
lines(x,y, lty=2, lwd=2)


# power law function (Martin's b)
zlm <- lm(log(isp.all$depth) ~ log(isp.all$poc)); summary(zlm) # b = 0.7
x <- exp((log(y) - coef(zlm)[1])/coef(zlm)[2]) # p<0.001, r2 = 0.88
lines(x,y, lwd=2)
legend('bottomright', c('b', 'z*'), lty=c(1,2), lwd=2, bty='n', cex=1.1)



## plot with salinity and temperature
# temperature
par(new=T)
plot(ctd$Temperature, ctd$Nominal.depth, type='n',  axes=F, xlab='', ylab='', ylim=c(500,0),
     xlim=c(8,14))
lines(ctd$Temperature, ctd$Nominal.depth, col='red', lwd=2)
axis(1, at=seq(8,14,2), cex.axis=1.2, col='red', col.axis='red')
mtext(expression(paste('Temperature (',degree,'C)')), 1, 2.1, cex=1.2, col='red')

# add salinity
par(new=T)
plot(ctd$Salinity, ctd$Nominal.depth, type='n',  axes=F, xlab='', ylab='', ylim=c(500,0),
     xlim=c(34,35.5))
lines(ctd$Salinity, ctd$Nominal.depth, col='blue', lwd=2)
axis(1, at=seq(34,35.5,0.5), cex.axis=1.2, line=3, col='blue', col.axis='blue')
mtext('Salinity', 1, 3.5, cex=1.2, col='blue')




#***********************************************************
# Take average of oxygen concentration at each temperature of each time slot
#***********************************************************

# Code for each row that retains the vials 
# set parameters as factors
resp$Net<-as.factor(resp$Net)
resp$Temp<-as.factor(resp$Temp)
resp$Vial<-as.factor(resp$Vial)
resp$Time<-as.factor(resp$Time)

resp <- within(resp, Code1 <- (Net:Temp:Vial:Time)[drop = TRUE])


# average o2 concentration per time per temp
resp_mean <- aggregate(resp$O2.umol_L, by=list(resp$Code1), FUN=mean)




# Add parameters to new data frame

#net
net_1<-rep('net1', 63)
net_2<-rep('net2', 63)
net_3<-rep('net3', 63)
net_4<-rep('net4', 63)
net_5<-rep('net5', 63)
resp_mean$Nets<-c(net_1, net_2, net_3, net_4, net_5)

#temperature
tempC<-rep(12, 21) # control 12 C
tempM<-rep(17, 21) # mid 17 C
tempH<-rep(19, 21) # high 20 C
temp2C<-rep(12, 21) # control 12 C
temp2M<-rep(16, 21) # mid 17 C
temp2H<-rep(21, 21) # high 20 C
temp3C<-rep(12, 21) # control 12 C
temp3M<-rep(15, 21) # mid 17 C
temp3H<-rep(22, 21) # high 20 C
temp4C<-rep(12, 21) # control 12 C
temp4M<-rep(18, 21) # mid 17 C
temp4H<-rep(20, 21) # high 20 C
temp5C<-rep(12, 21) # control 12 C
temp5M<-rep(14, 21) # mid 17 C
temp5H<-rep(15, 21) # high 20 C

resp_mean$Temp <- c(tempC, tempM, tempH,temp2C, temp2M, temp2H,temp3C, temp3M, temp3H,
                   temp4C, temp4M, temp4H,temp5C, temp5M, temp5H)




#Time
t <- c('T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')
resp_mean$Time <- rep(t, 15)

#Vial
A1<-rep('A1', 7)
A2<-rep('A2', 7)
A3<-rep('A3', 7)
B1<-rep('B1', 7)
B2<-rep('B2', 7)
B3<-rep('B3', 7)
C1<-rep('C1', 7)
C2<-rep('C2', 7)
C3<-rep('C3', 7)
Vial<-c(A1,A2,A3,B1,B2,B3,C1,C2,C3)
resp_mean$Vial<-rep(Vial, 5)


# edit column names for columns 1 and 2
names(resp_mean)[c(1:2)] <- c('Code1', 'O2.umol_L')
# check col names
names(resp_mean)



# add the 'real' time measurements were taken in (hours)
realtime <- c(0, 0.5, 1, 2, 4, 6, 7) # first experiment (Net 1) only ran for 7 hours
net1time<-rep(realtime, 9)

realtime2 <- c(0, 0.5, 1, 2, 4, 6, 8) # rest ran for 8 hours
nettime<-rep(realtime2, 9)

time1<- c(net1time, rep(nettime, 4))
resp_mean$Real_Time <- time1





#***********************************************************
# Test for significant decline in O2 over time
#***********************************************************

resp_mean$Net <- as.factor(resp_mean$Net)
resp_mean$Temp <- as.factor(resp_mean$Temp)
resp_mean$Vial <- as.factor(resp_mean$Vial)

resp_mean <- within(resp_mean, Code2 <- (Net:Temp:Vial)[drop = TRUE])

# test significance of regression with linear mixed effect model, model 1
resp_mod.1 <-lmList(O2.umol_L ~ Real_Time|Code2, na.action = na.omit, data = resp_mean)
summary(resp_mod.1) 
  # All temps significant decrease in o2 with time (over 7-8 hours) apart from rows 1:4 and 19 (p>0.05). Remove later.


# Use model coefficients to get oxygen uptake as umol/L/h
coef.mod.1 <- coef(resp_mod.1)
colnames(coef.mod.1) <- c("c", "slope") # slope is umol/L/h or uM/h and c is the intercept

# make slope (oxygen uptake rates) positive
coef.mod.1$slope <- coef.mod.1$slope * -1
head(coef.mod.1)  # view data frame



# Add parameters back in
# new column temp
tempC<-rep(12, 3) # mid 17 C
tempM<-rep(17, 3) # mid 17 C
tempH<-rep(19, 3) # high 20 C
temp2C<-rep(12, 3) # control 12 C
temp2M<-rep(16, 3) # mid 17 C
temp2H<-rep(21, 3) # high 20 C
temp3C<-rep(12, 3) # control 12 C
temp3M<-rep(15, 3) # mid 17 C
temp3H<-rep(22, 3) # high 20 C
temp4C<-rep(12, 3) # control 12 C
temp4M<-rep(18, 3) # mid 17 C
temp4H<-rep(20, 3) # high 20 C
temp5C<-rep(12, 3) # control 12 C
temp5M<-rep(14, 3) # mid 17 C
temp5H<-rep(15, 3) # high 20 C
coef.mod.1$Temp<-c(tempC, tempM, tempH,temp2C, temp2M, temp2H,temp3C, temp3M, temp3H,
             temp4C, temp4M, temp4H,temp5C, temp5M, temp5H)

# Time of sampling
time1<-rep(20.5, 9) # mid 17 C
time2<-rep(18, 9) # mid 17 C
time3<-rep(12.5, 9) # high 20 C
time4<-rep(12, 9) # high 20 C
time5<-rep(10, 9) # mid 17 C
coef.mod.1$Net_Time<-c(time1,time2,time3,time4,time5)








#***********************************************************
# Make plot of O2 concentration with time (Fig. 2a)
#***********************************************************

# open plotting space
# grdev <- function(...) {get(getOption("device"))(...)}
# quartz(width=7,height=7) # code for mac
# par(mar=c(5,5,1,1), mfrow=c(1,1), bty='o', cex=1)

# make new code column to take mean o2 concentration
resp_mean$Real_Time <- as.factor(resp_mean$Real_Time)
resp_mean <- within(resp_mean, Code3 <- (Temp:Real_Time)[drop = TRUE])

# Take mean of oxygen concentration for each temperature and measurement time
o2 <- aggregate(resp_mean$O2.umol_L, by=list(resp_mean$Code3), mean)
std <- function(x) sd(x)/sqrt(length(x)) # function for std error of the mean
o2.se <- aggregate(resp_mean$O2.umol_L, by=list(resp_mean$Code3), std) # std error o2 per temp


# normalise to first o2 concentration for each temperature to compare on plot
o2$norm <- c((o2$x[1:8]/o2$x[1]), (o2$x[9:15]/o2$x[9]), (o2$x[16:22]/o2$x[16]),
             (o2$x[23:29]/o2$x[23]), (o2$x[30:36]/o2$x[30]), (o2$x[37:43]/o2$x[37]), (o2$x[44:50]/o2$x[44]),
             (o2$x[51:57]/o2$x[51]), (o2$x[58:64]/o2$x[58]), (o2$x[65:71]/o2$x[65]))


# 'normalise' std error
o2.se$norm <- o2.se$x/100

### add other coloumns back into df      
o2 <- o2[-7,]         # At temp 12 the first experiment at 7 hours and the rest had 8 hours so remove this time
o2.se <- o2.se[-7,]   # add for std error
# now add factors back in
o2$time <- rep(c(0,0.5, 1,2,4,6,8), 10)
o2$temp <- c(rep(12, 7), rep(14, 7),rep(15, 7),rep(16, 7),rep(17, 7),rep(18, 7),rep(19, 7),
             rep(20, 7), rep(21, 7), rep(22, 7))
o2.se$time <- rep(c(0,0.5, 1,2,4,6,8), 10)
o2.se$temp <- c(rep(12, 7), rep(14, 7),rep(15, 7),rep(16, 7),rep(17, 7),rep(18, 7),rep(19, 7),
                rep(20, 7), rep(21, 7), rep(22, 7))


# make plot
plot(o2$time, o2$norm, bg=rainbow(10)[as.factor(o2$temp)], pch=21, ylim=c(0.7,1.03),
     xlab= 'Time (hours)', ylab = expression(paste('Normalised ',O[2],
                                                   ' concentration')), type='n', yaxt='n',
     cex.lab=1.4, cex.axis=1.4)
axis(2,  labels=FALSE)
text(y=seq(0.7, 1, 0.1),par('usr')[1], pos=2,labels=as.vector(c(0.70,0.80,0.90,1.00)),srt=0, xpd=TRUE, cex=1.2)

# make colour pallette
cols <- c(colorRampPalette(c( "dodgerblue4","lightblue"))((5)), colorRampPalette(c("rosybrown1", "firebrick"))(5))

# error bars
xy.error.bars<- function (x,y,ybar){
  points(o2$time,o2$norm,pch=21,  xaxt='n', xlab='',
         cex=1, xpd=FALSE, bg=cols[as.factor(o2$temp)])
  arrows(x, y-ybar, x, y+ybar, code=3, angle=90, length=0.05, lwd=1.2, xpd=FALSE,
         col=cols[as.factor(o2$temp)]) 
}
x<-o2$time
y<-o2$norm
yb<-o2.se$norm
xy.error.bars(x,y, yb)


# Get slopes of normalised regressions between time and o2 concentrations
resp_mod.2 <-lmList(norm ~ time|temp, na.action = na.omit, data = o2)
norm.coef <- coef(resp_mod.2)

time <- seq(0, 8, 0.1) # hours (x)
temp_12 <- 1 + (norm.coef$time[1] * time)  # (y)
temp_14 <- 1 + (norm.coef$time[2] * time)
temp_15 <- 1 + (norm.coef$time[3] * time)  
temp_16 <- 1 + (norm.coef$time[4] * time)  
temp_17 <- 1 + (norm.coef$time[5] * time)
temp_18 <- 1 + (norm.coef$time[6] * time)
temp_19 <- 1 + (norm.coef$time[7] * time)
temp_20 <- 1 + (norm.coef$time[8] * time)
temp_21 <- 1 + (norm.coef$time[9] * time)
temp_22 <- 1 + (norm.coef$time[10]  * time)

# lines  
lines(time, temp_12, col=cols[1], lwd=2)
lines(time, temp_14, col=cols[2], lwd=2)
lines(time, temp_15, col=cols[3], lwd=2)
lines(time, temp_16, col=cols[4], lwd=2)
lines(time, temp_17, col=cols[5], lwd=2)
lines(time, temp_18, col=cols[6], lwd=2)
lines(time, temp_19, col=cols[7], lwd=2)
lines(time, temp_20, col=cols[8], lwd=2)
lines(time, temp_21, col=cols[9], lwd=2)
lines(time, temp_22, col=cols[10], lwd=2)
temp<-c(12,14,15,16,17,18,19,20,21,22)
legend('bottomleft', c('12', '14','15','16', '17','18', '19','20', '21', '22'),
       pt.bg=cols[as.factor(temp)], bty='n', pch=21, cex=1.5)









#***********************************************************
# Plot temperature against oxygen consumption (Fig. 2b)
#***********************************************************
# open plotting space (Mac)
# grdev <- function(...) {get(getOption("device"))(...)}
# quartz(width=7,height=7) # code for mac
# par(mar=c(5,5,3,1), mfrow=c(1,1), bty='o', cex=1)

# calculate mean o2 consumption at each temperature
mean.net <- aggregate(coef.mod.1$slope, by=list(coef.mod.1$Temp), mean)
mean.net.sd <- aggregate(coef.mod.1$slope, by=list(coef.mod.1$Temp), std)

# create empty plot
plot((mean.net$Group.1), (mean.net$x), xlab=expression(paste('Temperature (',degree,'C)')), 
     ylab=expression(paste(,O[2], ' consumption (',mu,'mol L'^'-1'*' h'^'-1'*')')),
     pch=16, ylim=c(0,8), type='n', cex.lab=2, cex.axis=2, yaxt='n')
axis(2, at=seq(0,8,2), labels=FALSE)
text(y=seq(0,8,2),par('usr')[1], pos=2,labels=as.vector(c(0,2,4,6,8)),srt=0, xpd=TRUE, cex=2)

# error bars
xy.error.bars<- function (x,y,ybar){
  points(mean.net$Group.1, mean.net$x,  pch=16, xaxt='n', xlab='',
         cex=2, xpd=FALSE)
  arrows(x, y-ybar, x, y+ybar, code=3, angle=90, length=0.05, lwd=1.2, xpd=FALSE) 
}
x<-mean.net$Group.1
y<-mean.net$x
yb<-mean.net.sd$x
xy.error.bars(x,y,yb)

### calculate regression
r<-lm((coef.mod.1$slope)~(coef.mod.1$Temp)); summary(r) 
# set x (temperature)
temp_mod<-seq(0,25,0.1)
# set y (o2 consumption)
met<-(0.359*temp_mod)-2.53
# plot regression lines
lines(temp_mod, met, lwd=2)
legend('topleft', c('p < 0.01', expression(paste('R'^'2'*' = 0.5'))), bty='n', cex=2)

#### add confidence intervals
newx<-seq(10,24, length.out=45)
prd<-predict(r,newdata=data.frame(x=newx),interval = c("confidence"), 
             level = 0.95,type="response")
coef.mod.1$pred.l <-prd[,2]
coef.mod.1$pred.h <-prd[,3]
cf.int.l <- aggregate(coef.mod.1$pred.l, by=list(coef.mod.1$Temp), mean)
cf.int.h <- aggregate(coef.mod.1$pred.h, by=list(coef.mod.1$Temp), mean)
# coordinates
peep <- c(10, 0)
peep2 <- c(10, 1.4)
peep3 <- c(24, 7.2)
peep4 <- c(24, 5.1)
cf.int.h <- rbind( peep2, cf.int.h, peep3)
cf.int.l <- rbind( peep, cf.int.l, peep4)
# plot CI
polygon(c(cf.int.l$Group.1, rev(cf.int.l$Group.1)), c(cf.int.h$x, rev(cf.int.l$x)), col=adjustcolor('gray68', alpha.f=0.5),
        border=NULL, lty=0)













#***********************************************************
# Merge with POC data to calculate k (d-1)
#***********************************************************
#covert o2 consumption (slope) from umol/L/h to umol/h
head(coef.mod.1)
coef.mod.1$O2.umol.h <- coef.mod.1$slope * 0.02 # volume of vials = 20 mls = 0.02 L

# Reorder poc data by nets and vials first
poc <- poc[order(poc$Type, poc$Vial),]

#add to coef.mod1 data
coef.mod.1$poc <- poc$Correct.C..ug.
coef.mod.1$pon <- poc$Correct.N..ug.

# remove non significant slopes identified line *****
coef.mod.1 <- coef.mod.1[-c(1:4, 19),]

# convert POC from ug to umol
coef.mod.1$poc.umol <- coef.mod.1$poc/12


### calculate k (POC turnover rate or mass-normalised respirartion rate)
  # k = oxygen consumption / heterotrophic biomass
  # we assume heterotrophic biomass is 50 % of total mixed layer POC so
  # divide POC by 2. see more in text and sensitivity analysis runs later.

coef.mod.1$k <- (coef.mod.1$O2.umol.h / (coef.mod.1$poc.umol/2))*24   # * 24 converts from /h to /d









# Remove data from net 4 as very high k values - anomalies (k > +- 2 standard devations from the mean)
rm.net4 <- coef.mod.1[-c(23:31),]






#***********************************************************
# Apply linear mixed effect models and make plots (Figs. 4a and b)
#***********************************************************




# linear mixed effect model 
        # with net deployment time as the random effect
fm01 <- lmer(log(k) ~ Temp + (1|Net_Time), rm.net4, REML=F); summary(fm01)

# linear model 
fm03 <- lm(log(k) ~ Temp, rm.net4); summary(fm03); anova(fm01, fm03)





########### Plot Fig. 4a

# calculate mean k per temperature
k_data <- aggregate(rm.net4$k, by=list(rm.net4$Temp), FUN=mean)

# calculate standard error of k per temperature for error bars
k_data_std <- aggregate(rm.net4$k, by=list(rm.net4$Temp), FUN=std)


# open plotting device (mac)
# grdev <- function(...) {get(getOption("device"))(...)}
# quartz(width=7,height=7) # code for mac
# par(mar=c(5,5,1,1), mfrow=c(1,1), bty='o', cex=1)


# make basic plot 
plot(k_data$Group.1, k_data$x, pch=16, ylab=expression(paste('k (d'^'-1'*')')), 
         xlab=expression(paste('Temperature (',degree,'C)')), ylim=c(0,5),
         cex=2, cex.lab=1.5, type='n', xlim=c(11,22), cex.axis=1.4, yaxt='n')
    axis(2, at=seq(0,5,1), labels=FALSE)
    text(y=seq(0,5,1),par('usr')[1], pos=2,labels=as.vector(c(0,1,2,3,4,5)),srt=0, xpd=TRUE, cex=1.5)
    
    
    
    ### plot CI polygon
    #first CI and PI using predict-like method, using code posted here: http://glmm.wikidot.com/faq
    newdat<-data.frame(x=seq(10,24,length=20))    # create new data set using x
    newdat$y <- exp(-1.85847+ ( 0.12178* newdat$x))   # simulate y from lmer regression line
    mm <- as.matrix(data.frame(rep(1,20), newdat$x))  # form matrix of intercept (1s) and x
    
    pvar1 <- diag(mm %*% tcrossprod(vcov(fm01),mm))   # cross multiply matrix with lmer
    newdat <- data.frame(
      newdat
      , plo = newdat$y-1.96*sqrt(pvar1)
      , phi = newdat$y+1.96*sqrt(pvar1)
    )                                               # calc lower and upper CIs
    
    polygon(c(newdat$x, rev(newdat$x)), c(newdat$plo, rev(newdat$phi)), col=adjustcolor('gray68', alpha.f=0.6),
            border=NULL, lty=0)
    
    
    
    ### CI for k's between 12- 15 C temperature, the predicted increase by 2100
    
    # make dataframe with just temps 12 - 15 C
    t.2100 <- rm.net4[rm.net4$Temp<15.1,]
    # linear mixed effect model 
    fm01.2100 <- lmer(log(k) ~ Temp + (1|Net_Time), t.2100, REML=F); summary(fm01)
    
    #first CI and PI using predict-like method, using code posted here: http://glmm.wikidot.com/faq
    newdat<-data.frame(x=seq(12,15,length=9))    # create new data set using x
    newdat$y <- exp(-1.53993+ ( 0.11977* newdat$x))   # simulate y from lmer regression line
    mm <- as.matrix(data.frame(rep(1,9), newdat$x))  # form matrix of intercept (1s) and x
    
    pvar1 <- diag(mm %*% tcrossprod(vcov(fm01.2100),mm))   # cross multiply matrix with lmer
    newdat <- data.frame(
      newdat
      , plo = newdat$y-1.96*sqrt(pvar1)
      , phi = newdat$y+1.96*sqrt(pvar1)
    )                                               # calc lower and upper CIs
    
    polygon(c(newdat$x, rev(newdat$x)), c(newdat$plo, rev(newdat$phi)), col=adjustcolor('gray50', alpha.f=0.6),
            border=NULL, lty=0)
    
    
    
    # plot points and error bars 
    xy.error.bars<- function (x,y,ybar){
      points(k_data$Group.1, k_data$x,pch=16,  xaxt='n', xlab='',
             cex=2, xpd=FALSE)
      arrows(x, y-ybar, x, y+ybar, code=3, angle=90, length=0.05, lwd=2, xpd=FALSE) 
    }
    x<-k_data$Group.1 
    y<-k_data$x
    yb<-k_data_std$x
    xy.error.bars(x,y,yb)
    
    # mean k
    mean(rm.net4$k) # mean k is 1.5 d-1
    
    
    # plot regression lines
    # over 10 C range
    x<- seq(10,24,1)
    y <- exp(-1.85847+ ( 0.12178* x)) # mixed effect model inc. time of sample, p<0.001 r2 = 0.45, fm01
    
    # estimate r2 = 0.45
    
    # SSY:
    SSY <- sum(k_data$x^2)-sum(k_data$x)^2/length(k_data$x)
    print(SSY)
    # SSE (predicted): 
    pred <-exp(-1.85847+ ( 0.12178* k_data$Group.1))  # Predicted y
    SSE <- sum((k_data$x-pred)^2)
    print(SSE)  
    # R2
    print((SSY-SSE)/SSY)
    
    # plot lines  
    
    lines(x,y, lwd=3)
    
    
    
    
    
    
    # do it for 12 - 15 C (temperature by year 2100)
    
    # mean
    k.2100.mean <- aggregate(t.2100$k, by=list(t.2100$Temp), FUN=mean)
    # calculate standard error
    k_2100_std <- aggregate(t.2100$k, by=list(t.2100$Temp), FUN=std)
    
    
    # estimate r2
    
    # SSY:
    SSY <- sum(k.2100.mean$x^2)-sum(k.2100.mean$x)^2/length(k.2100.mean$x)
    print(SSY)
    # SSE (predicted): 
    pred <-exp(-1.53993+ ( 0.11977* k.2100.mean$Group.1))  # Predicted y
    SSE <- sum((k.2100.mean$x-pred)^2)
    print(SSE)  
    # R2
    print((SSY-SSE)/SSY) # 0.49
    
    
    # plot lines  
    x<- seq(12,15,1)
    y <- exp(-1.53993+ ( 0.11977* x)) # mixed effect model inc. time of sample, p<0.001 r2 = 0.45, fm01
    lines(x,y, lwd=3, lty=2)
    
    
    legend('topleft', c('Experimental range', '2100 range'), lty=c(1,2), lwd=2, bty='n', cex=1.5 )
    
    








### Fig. 4b

# open plotting device (mac)
# grdev <- function(...) {get(getOption("device"))(...)}
# quartz(width=7,height=7) # code for mac
# par(mar=c(5,5,1,1), mfrow=c(1,1), bty='o', cex=1)



# calculate parameters to apply mte
rm.net4$kelvin <- rm.net4$Temp + 273.15 # temperature in Kelvin (K)
constant <- 0.0000862 # Boltzmann constant (units = eV / K)
rm.net4$t <- 1 / (constant * rm.net4$kelvin) # x-axis (1/botzmann constant * temp in Kelvin) units = eV

# construct linear models
fm01 <- lmer(log(k) ~ t + (1|Net_Time), rm.net4, REML=F); summary(fm01) # p<0.01 r2=0.4
fm03 <- lm(log(k) ~ t, rm.net4); summary(fm03); anova(fm01 ,fm03)



#### make plot
plot(rm.net4$t, log(rm.net4$k), ylim=c(-1,2), ylab=expression(paste('ln(k)')), xlab='1/cT',
         cex.lab=1.5, cex=2, type='n', xlim=c(39.2,40.7), cex.axis=1.4, yaxt='n')
    axis(2, at=seq(-1, 2, 0.5), labels=FALSE)
    text(y=seq(-1, 2, 0.5),par('usr')[1], pos=2,labels=as.vector(seq(-1, 2, 0.5)),srt=0, xpd=TRUE, cex=1.5)
    
    
    # CI for normal (around averaged data, this is what R2 from)
    k_data$kelvin <- k_data$Group.1 + 273.15
    k_data$t <-1 / (constant * k_data$kelvin)
    
    newdat<-data.frame(x=seq(38,42,length=20))    # create new data set using x
    newdat$y <- 35.3366 + (-0.8784 *  newdat$x)  # simulate y from lmer regression line
    mm <- as.matrix(data.frame(rep(1,20), newdat$x))  # form matrix of intercept (1s) and x
    
    pvar1 <- diag(mm %*% tcrossprod(vcov(fm01),mm))   # cross multiply matrix with lmer
    newdat <- data.frame(
      newdat
      , plo = newdat$y-1.96*sqrt(pvar1)
      , phi = newdat$y+1.96*sqrt(pvar1)
    )                                               # calc lower and upper CIs
    
    polygon(c(newdat$x, rev(newdat$x)), c(newdat$plo, rev(newdat$phi)), col=adjustcolor('gray68', alpha.f=0.6),
            border=NULL, lty=0)
    
    
    
    
    
    # CI for 12 -15 C (temperature by year 2100)
    # compute the means
    t.2100.mean <- aggregate(t.2100$k, by=list(t.2100$Temp), FUN = mean)
    t.2100.std <- aggregate(t.2100$k, by=list(t.2100$Temp), FUN = std)
    
    # LME model from all data (not means)
    t.2100$kelvin <- t.2100$Temp + 273.15
    t.2100$t <- 1 / (constant * t.2100$kelvin)
    fm01 <- lmer(log(k) ~ t + (1|Net_Time), t.2100, REML=F); summary(fm01) # Ea = -0.85
    
    
    newdat<-data.frame(x=seq(40.2,40.7,length=20))    # create new data set using x
    newdat$y <- 34.3539 + (-0.8469 *  newdat$x)  # simulate y from lmer regression line
    mm <- as.matrix(data.frame(rep(1,20), newdat$x))  # form matrix of intercept (1s) and x
    
    pvar1 <- diag(mm %*% tcrossprod(vcov(fm01),mm))   # cross multiply matrix with lmer
    newdat <- data.frame(
      newdat
      , plo = newdat$y-1.96*sqrt(pvar1)
      , phi = newdat$y+1.96*sqrt(pvar1)
    )                                               # calc lower and upper CIs
    
    polygon(c(newdat$x, rev(newdat$x)), c(newdat$plo, rev(newdat$phi)), col=adjustcolor('gray50', alpha.f=0.6),
            border=NULL, lty=0)
    
    
    
    # plot points
    
    # error bars norm
    xy.error.bars<- function (x,y,ybar){
      points(k_data$t, log(k_data$x),pch=16,  xaxt='n', xlab='',
             cex=2, xpd=FALSE)
      arrows(x, y-ybar, x, y+ybar, code=3, angle=90, length=0.05, lwd=2, xpd=FALSE) 
    }
    x<-k_data$t
    y<-log(k_data$x)
    yb<-(k_data_std$x)
    xy.error.bars(x,y,yb)
    
    
    
    x<- seq(39,41, 0.5)
    y<- 35.3366 + (-0.8784 * x) 
    lines(x,y, lwd=2)
    
    # SSY:
    SSY <- sum(log(k_data$x)^2)-sum(log(k_data$x))^2/length(log(k_data$x))
    print(SSY)
    # SSE (predicted):
    pred <-(35.3366+ ( -0.8784* k_data$t))  # Predicted y
    SSE <- sum((log(k_data$x)-(pred))^2)
    print(SSE)
    # R2
    print((SSY-SSE)/SSY) # 0.40
    
    
    
    
    
    ####### do for 12 -15 C
    
    # SSY:
    SSY <- sum(log(t.2100.mean$x)^2)-sum(log(t.2100.mean$x))^2/length(log(t.2100.mean$x))
    print(SSY)
    # SSE (predicted):
    pred <-(34.3539+ ( -0.8469* t.2100.mean$t))  # Predicted y
    SSE <- sum((log(t.2100.mean$x)-pred)^2)
    print(SSE)
    # R2
    print((SSY-SSE)/SSY) # 0.41
    
    x<- seq(40.2,40.7, 0.05)
    y<- 34.3539 + (-0.8469 * x) # slope = 1.3, higher than mte theory r2=0.68 p<0.001
    lines(x, y, lwd=3, lty=2)
    
    
    legend('topright', c('Ea = 0.88 eV', 'Ea = 0.85 eV'), lty=c(1,2), lwd=2, bty='n', cex=1.5 )
    
    

    
    
    
    
    
    
#***********************************************************
# Test effect of assuming a constant heterotrophic biomass in all vials
#***********************************************************

# Strzepek 2005 - hetero = 46 - 56 % of total poc (avg = 49 %) Plot in supplementary materials
  
# 8 different temps, randomly associate number from 0.4 - 0.6 and then 0.46 - 0.56 to each temp, 
    # then recalculate k




#40 - 60 %

a <- vector()
ae <- vector()
hea.c <- vector()

for (i in 1:1000){    # run 1000 times
  rm.net4$HB <- runif(nrow(rm.net4), min=0.4, max=0.6) # assign random number to each row
  rm.net4$k_vary_hb <- (rm.net4$O2.umol.h / (rm.net4$poc.umol * rm.net4$HB)) * 24 # recalculate k using random number
  lm_MTE_01 <- lmer(log(k_vary_hb) ~ t + (1|Net_Time), rm.net4, REML=F) # make linear model
  e<- coef(summary(lm_MTE_01))[2,1] # get the slope
  a <- c(a, e) # create vector with all 1000 slopes
  ce <- coef(summary(lm_MTE_01))[1,1] # get intersect
  ae<- c(ae, ce) # creat vector of all 1000 interects
}

# comibine to make datafram (for supplementary plot)
hea.c <- cbind(a, ae)
hea.c <- as.data.frame(hea.c)

# Mean Ea using biomass between 40 - 60 %
mean(hea.c$a)
std(hea.c$a) # std error for the mean



# open plotting space (Mac)
# grdev <- function(...) {get(getOption("device"))(...)}
# quartz(width=7,height=7) # code for mac
# par(mar=c(5,5,3,1), mfrow=c(1,1), bty='o', cex=1)

# plot all 1000 lines on Fig. 4b
plot(k_data$t,log(k_data$x), pch=16, ylab=expression(paste('ln(k)')), xlab='1/cT')
# add 1000 lines (40 - 60 %)
for (i in 1:nrow(hea.c)){
  x <- seq(39, 41, length.out=1000)
  y<- hea.c$ae[i] + (hea.c$a[i] * x)
  lines(x,y, lwd=0.4, col='grey')
}





# Repeat for 46 - 56 % (actual range in data from Strzepek 2005)
ea <- vector()
c <- vector()
ea.c <- vector()

for (i in 1:1000){
  rm.net4$HB <- runif(nrow(rm.net4), min=0.46, max=0.56)
  rm.net4$k_vary_hb <- (rm.net4$O2.umol.h / (rm.net4$poc.umol * rm.net4$HB)) * 24
  lm_MTE_01 <- lmer(log(k_vary_hb) ~ t + (1|Net_Time), rm.net4, REML=F)
  e<- coef(summary(lm_MTE_01))[2,1]
  ea <- c(ea, e)
  ce <- coef(summary(lm_MTE_01))[1,1]
  c<- c(c, ce)
}

# Mean Ea using biomass between 46 - 56 %
mean(ea)
std(ea)
# put in df
ea.c <- cbind(ea, c)
ea.c <- as.data.frame(ea.c)

# plot 1000 lines where biomass is 46 - 56 % of POC
for (i in 1:nrow(ea.c)){
  x <- seq(39, 41, length.out=1000)
  y<- ea.c$c[i] + (ea.c$ea[i] * x)
  lines(x,y, lwd=0.4)
}

# Add line used in main analysis, where biomass = 50 % of POC
x<- seq(39,41, 0.5)
y<- 35.3366 + (-0.8784 * x) 
lines(x, y, col='orange')

legend('topright', c('HB = 40 - 60 % of POC','HB = 46 - 56 % of POC', 'HB = 50 % POC'), col = c('grey', 'black', 'orange'), lwd=2, bty='n')
points(k_data$t,log(k_data$x), pch=16) # add points again as some masked by lines
    





#***********************************************************
# Reaction norm, Fig. 5.
#***********************************************************
temp <- seq(7, 22, 1)
temp.b  <- seq(7, 17, 1)
#boyd <-c(0,0.2, 0.4,0.6, 0.8, 1, 0.99, 0.98, 0.97,0.96, 0.95)
boyd <- c(0.2, 0.5, 0.8, 0.93, 0.98, 1, 0.99, 0.98, 0.97, 0.96, 0.95)
#plot(temp.b,boyd, type='l')
spi <- spline(temp.b,boyd)

par(mfrow=c(1,1), mfcol=c(1,1), mar=c(5,5,2,2))

plot(spi$x, spi$y, type='l', xlim=c(7,24), ylim=c(0.2,2.4), ylab='', 
     xlab=expression(paste('Temperature (',degree,'C)')),lty=1, lwd=2, col='red', yaxt='n')
axis(2, at=seq(0.2,2.4,0.4), labels=FALSE)
text(y=seq(0.2,2.4,0.4),par('usr')[1], pos=2,labels=as.vector(sprintf("%.1f", round(seq(0.2,2.4,0.4),2))),srt=0, xpd=TRUE, cex=1)


mtext(expression(paste('k (d'^'-1'*')')), 2, 3.5)
mtext(expression(paste(,mu,'/',mu,'max')), 2, 2)

x<- seq(7,24,1)
y <- exp(-1.85847+ ( 0.12178* x)) # mixed effect model inc. time of sample, p<0.001 r2 = 0.45
lines(x, y, lwd=2, lty=1)

x <- c(12, 14, 15, 16,17,19,21,22) # temp exp done at
y <- exp(-1.85847+ ( 0.12178* x))
points(x,y, pch=16)

# control temperature of this study (12 C)
abline(0,1 ,v=12, col='black', lty=2, lwd=2)
arrows(13, 2.3, 14.9, 2.3, length=0.1, lwd=2, col='black', lty=1)

# in situ temperature of Boyd et al study
abline(0,1 ,v=10.6, col='red', lty=1, lwd=2)
arrows(10.6, 2.2, 13.7, 2.2, length=0.1, lwd=2, col='red', lty=1)

arrows(8, 2.6, 11.8, 2.6, length=0.1, lwd=2, col='black', xpd=TRUE, code=3)
arrows(12.2, 2.6, 22, 2.6, length=0.1, lwd=2, col='black', xpd=TRUE, code=3)

legend('bottomright', c('Boyd et al.', 'This study'),  lwd=2, bty='n', cex=0.8, 
       col=c('red', 'black'))

