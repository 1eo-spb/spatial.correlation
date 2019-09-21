setwd("C:/Users/surov/OneDrive/MTU/01. Research/06. Spatial correlation/r")
# C:\Users\surov\OneDrive\MTU\01. Research\06. Spatial correlation\r

# suscep -> param1 (Pentane)
# gracit -> param2 (Propane)
# x is for Longitude
# y is for Latitude

library(rworldmap)
library(ggplot2)
library(ggsn)
library(gstat)
library(sp)
library(maptools)
library(raster)
library(rgdal)
library(lattice)
library(DAAG)
library(spdep)
library(ggmap)
register_google(key="xxx",account_type="standard", write = TRUE) # https://console.cloud.google.com

################################################################################
#
# Data setting functions
# R1DS - Read 1st Data Set from a file Parameter-1.csv
# R2DS - Read 2nd Data Set from a file Parameter-2.csv
#
################################################################################
R1DS <- function(file)
{
  param1.raw <- read.table(file=file, header=TRUE, sep=",", na.strings=T);
  param1.dat <- param1.raw[,c(1,2,3)];
  #names(param1.dat) <-c("Lat","Long",'Param1')
}
#####
R2DS <- function(file)
{
  param2.raw <- read.table(file=file, header=TRUE, sep=",", na.strings=T);
  param2.dat <- param2.raw[,c(1,2,3)];
  #names(param2.dat) <-c("Lat","Long",'Param2')
}
################################################################################
#
# The function calculates the data limits in data sets 
# 0 - find - Just fild all the values to use it further
# x - calc - calculate the limits to be set automatically where 'x' is  persentage over the limits
#
################################################################################
LIMITS <- function(what)
{
  x = what;
  if (x == 0) {
    max_lat_p1 <- max(param1.dat[,1]);
    min_lat_p1 <- min(param1.dat[,1]);
    max_lat_p2 <- max(param2.dat[,1]);
    min_lat_p2 <- min(param2.dat[,1]);
    
    max_lon_p1 <- max(param1.dat[,2]);
    min_lon_p1 <- min(param1.dat[,2]);
    max_lon_p2 <- max(param2.dat[,2]);
    min_lon_p2 <- min(param2.dat[,2]);
    
    mean_lat_p1 <- mean(param1.dat[,1]);
    mean_lat_p2 <- mean(param2.dat[,1]);
    mean_lon_p1 <- mean(param1.dat[,2]);
    mean_lon_p2 <- mean(param2.dat[,2]);
    
  }  else {
     max_lat_p1 <- max(param1.dat[,1]);
     min_lat_p1 <- min(param1.dat[,1]);
     max_lat_p2 <- max(param2.dat[,1]);
     min_lat_p2 <- min(param2.dat[,1]);
     
     max_lon_p1 <- max(param1.dat[,2]);
     min_lon_p1 <- min(param1.dat[,2]);
     max_lon_p2 <- max(param2.dat[,2]);
     min_lon_p2 <- min(param2.dat[,2]);
     
     mean_lat_p1 <- mean(param1.dat[,1]);
     mean_lat_p2 <- mean(param2.dat[,1]);
     mean_lon_p1 <- mean(param1.dat[,2]);
     mean_lon_p2 <- mean(param2.dat[,2]);
     
     if (max_lat_p1 >= max_lat_p2) {
       lim_lat_max <- max_lat_p1; 
     } else {lim_lat_max <- max_lat_p2};
     
     if (min_lat_p1 <= min_lat_p2) {
       lim_lat_min <- min_lat_p1; 
     } else {lim_lat_min <- min_lat_p2};
     
     if (max_lon_p1 >= max_lon_p2) {
       lim_lon_max <- max_lon_p1; 
     } else {lim_lon_max <- max_lon_p2};
     
     if (min_lon_p1 <= min_lon_p2) {
       lim_lon_min <- min_lon_p1; 
     } else {lim_lon_min <- min_lon_p2};
     
     eax <- lim_lat_max + (x * ((lim_lat_max - lim_lat_min) / 100));
     ean <- lim_lat_min - (x * ((lim_lat_max - lim_lat_min) / 100));
     eox <- lim_lon_max + (x * ((lim_lon_max - lim_lon_min) / 100));
     eon <- lim_lon_min - (x * ((lim_lon_max - lim_lon_min) / 100));
    }
}
################################################################################
#
# Data printing
# SDP - Simple Data Plot
# xmin/xmax - Longitude limits
# ymin/ymax - Latitude limits
#
################################################################################
SDP <- function(xmin,xmax,ymin,ymax)
{
  par(mfrow = c(1,1), xaxs = 'i', yaxs = 'i', mar = c(4,4,3,2.5),cex.lab = 1, cex.axis = 1.0)
  plot(1, 1, type = 'n', xlim = c(xmin=xmin, xmax=xmax), ylim = c(ymin=ymin,ymax=ymax), ylab ='', axes = FALSE, xlab='' )
  box(col = 'grey')
  grid(col = 'grey', lty = 1)
  points(param1.dat$Long, param1.dat$Lat, type="p", pch=1, col="green")
  points(param2.dat$Long, param2.dat$Lat, type="p", pch=4, col="red")
  title(main="Sampling location")
  axis(side = 1, tick = FALSE, col = 'grey', line = -0.8)
  axis(side = 2, tick = FALSE, col = 'grey', line = 0, las = 0)
  mtext(expression("Latitude"), side=2, line=2.5)
  mtext(expression("Longitude"), side=1, line=2.5)
}
################################################################################
#
# Printing on a global world map
#
################################################################################
# dev.off() 
# newmap <- getMap(resolution = "coarse") # different resolutions available
# plot(newmap)
# points(Long, Lat, data=param1.dat, col="red")
################################################################################
#
# Printing on a Google map
# GMP - Google Map Plot
# type  - "terrain", "terrain-background", "satellite", "roadmap", "hybrid", "terrain", "watercolor", and "toner"
# Zoom scale - 1-21
# data  - 1 for Graviti, 2 for susceptibility, 3 for both
#
################################################################################
GMP <- function(type,zscale,data)
{
 n = data;
 googlat.s <- mean(c(ean,eax)) # coordinate of the center of param1.area (latitude)
 googlon.s <- mean(c(eon,eox)) # coordinate of the center of param1.area (longitude)
 d.mich.s <- get_map(location=c(googlon.s,googlat.s), maptype=(type=type),color="color", scale=2, zoom=(zscale=zscale)) # create map for param1 
 if (n == 1) {
    ggmap(d.mich.s)+geom_point(aes(param1.dat$Long, param1.dat$Lat, size=param1.dat[,3]), colour="yellow", data=param1.dat, na.rm=T)
  } else if (n == 2) {
    ggmap(d.mich.s)+geom_point(aes(param2.dat$Long, param2.dat$Lat, size=param2.dat[,3]), colour="red", data=param2.dat, na.rm=T) #+geom_point(aes(param1.dat$Long, param1.dat$Lat, size=param1.dat$Pentane), colour="red", data=param1.dat, na.rm=T);
  } else if (n == 3) {
    ggmap(d.mich.s)+geom_point(aes(param2.dat$Long, param2.dat$Lat, size=param2.dat[,3]), colour="red", data=param2.dat, na.rm=T)+geom_point(aes(param1.dat$Long, param1.dat$Lat, size=param1.dat[,3]), colour="yellow", data=param1.dat, na.rm=T);
  } else if (n == 4) {
    ggmap(d.mich.s)+geom_point(aes(param1.dat$Long, param1.dat$Lat, size=param1.dat[,3]), colour="yellow", data=param1.dat, na.rm=T)+geom_point(aes(param2.dat$Long, param2.dat$Lat, size=param2.dat[,3]), colour="red", data=param2.dat, na.rm=T);
  } else if (n == 5) {
    ggmap(d.mich.s)+geom_point(aes(param2.dat$Long, param2.dat$Lat, size=0.5), colour="red", data=param2.dat, na.rm=T)+geom_point(aes(param1.dat$Long, param1.dat$Lat, size=0.5), colour="yellow", data=param1.dat, na.rm=T);
  } else if (n == 6) {
    ggmap(d.mich.s)+geom_point(aes(param1.dat$Long, param1.dat$Lat, size=0.5), colour="yellow", data=param1.dat, na.rm=T)+geom_point(aes(param2.dat$Long, param2.dat$Lat, size=0.5), colour="red", data=param2.dat, na.rm=T);
  } else {print ("Wrong argument for data type! :P")}
}
# dev.off()
# googlat.s <- mean(param1.dat$Lat)  # coordinate of the center of susc.area (latitude)
# googlon.s <- mean(param1.dat$Long) # coordinate of the center of susc.area (longitude)
# googlat.g <- mean(param2.dat$Lat)  # coordinate of the center of grav.area (latitude)
# googlon.g <- mean(param2.dat$Long) # coordinate of the center of grav.area (longitude)
# d.mich.s <- get_map(location=c(googlon.s,googlat.s),zoom=14, scale=2, maptype="hybrid", source ="google", crop=TRUE, color="color", language="en-EN") # create map for param1 
# d.mich.g <- get_map(location=c(googlon.s,googlat.s),zoom=14, scale=2, maptype="hybrid", source ="google", crop=TRUE, color="color", language="en-EN") # create map for param2
# ggmap(d.mich.s)+geom_point(aes(param1.dat$Long, param1.dat$Lat, size=param1.dat$Pentane), colour="red", data=param1.dat, na.rm=T)
# ggmap(d.mich.g)+geom_point(aes(param2.dat$Long, param2.dat$Lat, size=param2.dat$Propane), colour="red", data=param2.dat, na.rm=T)
# ggmap(d.mich.s)+geom_point(aes(param1.dat$Long, param1.dat$Lat, size=param1.dat$Pentane), colour="red", data=param1.dat, na.rm=T)+geom_point(aes(param2.dat$Long, param2.dat$Lat, size=param2.dat$Propane), colour="yellow", data=param2.dat, na.rm=T)
################################################################################
#
#
#
################################################################################
# param1.sp <- param1.dat
# param2.sp <- param2.dat
# bbox(param1.sp)
# bbox(param2.sp)
# coordinates(param1.sp) <- ~ Long + Lat
# coordinates(param2.sp) <- ~ Long + Lat
# concentration using the bubble function
# bubble(param1.sp, "Pentane",col= "red", main = "Paramener #1 (Pentane)", axes=T)
# bubble(param2.sp, "Propane",col= "red", main = "Paramener #2 (Propane)", axes=T)
# 
# display.brewer.all()
#
# param1.pred <- param2.idw[["var1.pred"]]
# param2.pred <- param2.idw[["var1.pred"]]
################################################################################
#
# Transform to coordinate & Calculating of coefficient of correlation
# CCS - Calculation of coefficients of correlation for Pentane
# CCG - Calculation of coefficients of correlation for Propane
# idp - weight
#
################################################################################
CCS <- function(ipd)
{
  param1.sp <- param1.dat
  coordinates(param1.sp) <- ~ Long + Lat
  ns <- 0;
  ss <- 1;
  while(ss <= (ipd=ipd)){
    param1.idw.valid5 <- krige.cv(Pentane~1, param1.sp, nfold=5, set=list(idp=ss))
    ns[ss] <- cor(param1.idw.valid5$observed, param1.idw.valid5$var1.pred)
    ss <- ss + 1
  }
  ccs <- which.max(ns);
}
#####
CCG <- function(ipd)
{
  param2.sp <- param2.dat
  coordinates(param2.sp) <- ~ Long + Lat
  ng <- 0;
  sg <- 1;
  while(sg <= (ipd=ipd)){
    param2.idw.valid5 <- krige.cv(Propane~1, param2.sp, nfold=5, set=list(idp=sg))
    ng[sg] <- cor(param2.idw.valid5$observed, param2.idw.valid5$var1.pred)
    sg <- sg + 1
  }
  ccg <- which.max(ng);
}
################################################################################
#
# Printing of spatial distribution
# SPG - Spatial Plot of Gravity
# SPS - Spatial Plot of Susceptibility
# grd - dimention of grid
# mxd - only observations within a distance of maxdist from the prediction location are used for prediction or simulation
#
################################################################################
SPS <- function(grd,mxd)
{
  param1.sp <- param1.dat;
  coordinates(param1.sp) <- ~ Long + Lat;
  param1.grid <- Sobj_SpatialGrid(param1.sp,maxDim=(grd=grd))$SG;
  gridded(param1.grid) = TRUE
  param1.idw <- krige(Pentane~1, param1.sp, param1.grid, maxdist=(mxd=mxd), set=list(idp=cs));
  spplot(param1.idw["var1.pred"], main = "Pentane distribution inverse distance weighted interpolations", xlab="Longitude",ylab="Latitude", scales = list(draw = TRUE), panel = panel.gridplot);
  }
#####
SPG <- function(grd,mxd)
{
  param2.sp <- param2.dat
  coordinates(param2.sp) <- ~ Long + Lat
  param2.grid <- Sobj_SpatialGrid(param2.sp,maxDim=(grd=grd))$SG
  gridded(param2.grid) = TRUE
  param2.idw <- krige(Propane~1, param2.sp, maxdist=(mxd=mxd), param2.grid, set=list(idp=cg))
  spplot(param2.idw["var1.pred"], main = "Propane distribution inverse distance weighted interpolations", xlab="Longitude",ylab="Latitude", scales = list(draw = TRUE), panel = panel.gridplot);
}

################################################################################
#
# PGS - Polygon for Grav & Sucs
#
################################################################################

PGS <- function(pol1,pol2)
{
# names(param1.dat) <- c("Lat","Long","Pentane")
  param1.spt <- param1.dat;
  coordinates(param1.spt) <- ~ Long + Lat;
  param2.spt <- param2.dat;
  coordinates(param2.spt) <- ~ Long + Lat;
  
  r1 <- pol1
  r2 <- pol2
  
  poly1 <- Polygon(r1)
  poly2 <- Polygon(r2)
  
  poly.list.g <- Polygons(list(poly1),1)
  poly.list.s <- Polygons(list(poly2),1)
  
  surv <- SpatialPolygons(list(poly.list.g))
  samp <- SpatialPolygons(list(poly.list.s))

  par(mfrow = c(1, 2));
  plot(param1.spt, axes=T, xlim=lolim, ylim=lalim, col="blue", xlab="Longitude", ylab="Latitude", main="Pentane")
  plot(samp, border="red", add=T)
  plot(surv, border="green", add=T)
  plot(param2.spt, axes=T, xlim=lolim, ylim=lalim, col="orange", xlab="Longitude", ylab="Latitude", main="Propane")
  plot(samp, border="red", add=T)
  plot(surv, border="green", add=T)
}

################################################################################
#
# RCV - Convert to raster & cross validation
# grd - dimention of grid
# mxd - only observations within a distance of maxdist from the prediction location are used for prediction or simulation
#
################################################################################
RCV <- function(grd,mxd)
{
# names(param1.dat) <-c("Lat","Long","Pentane")
  param1.spt <- param1.dat;
  coordinates(param1.spt) <- ~ Long + Lat;
  param2.spt <- param2.dat;
  coordinates(param2.spt) <- ~ Long + Lat;
  
  param1.grd <- Sobj_SpatialGrid(param1.spt, maxDim=(grd=grd))$SG;
  param2.grd <- Sobj_SpatialGrid(param2.spt, maxDim=(grd=grd))$SG;
  
  gridded(param1.grd) = TRUE
  gridded(param2.grd) = TRUE
  
  param1.idw <- krige(Pentane~1, param1.spt, param1.grd,  maxdist=(mxd=mxd), set=list(idp=cs));
  param2.idw <- krige(Propane~1, param2.spt, param2.grd,  maxdist=(mxd=mxd), set=list(idp=cg));
  
  param1.ras <- raster(param1.idw, layer=1,  values=TRUE);
  param2.ras <- raster(param2.idw, layer=1,  values=TRUE);
  
  param1.val <- param1.spt[["Pentane"]]
  param2.val <- extract(param2.ras, param1.spt, method='bilinear')
  
  val <- data.frame(cbind(param1.val,param2.val));
  
  par(mfrow = c(1, 2))
  cc1 <- CVlm(data=val, form.lm=formula(param1.val~param2.val), plotit="Observed")
  cc2 <- CVlm(data=val, form.lm=formula(param1.val~param2.val), plotit="Residual")

  # param2.val <- param2.spt[["Propane"]]             #temp
  # param1.val <- extract(param1.ras, param2.spt)     #temp
  # val <- data.frame(cbind(param2.val,param1.val));  #temp
  # cc1 <- CVlm(data=val, form.lm=formula(param2.val~param1.val), plotit="Observed")  #temp
  # cc2 <- CVlm(data=val, form.lm=formula(param2.val~param1.val), plotit="Residual")  #temp
  
  # Not run: The basic syntax for a regression analysis in R is lm(Y~X)
  # where Y is the dependent variable to be predicted and X is the is the independent variable measured.
  # cc1[c("ms","df")]
  # df - data.frame
  # ms - Mean square
}
  
################################################################################
#
# TMC - Two Maps Comparison
# grd - dimention of grid
# mxd - only observations within a distance of maxdist from the prediction location are used for prediction or simulation
#
################################################################################

TMC <- function(grd,mxd)
{
  dat <- cbind.data.frame(param2.dat, param1.val)
  dat <- dat[,c(1,2,3,4)]
  names(dat) <- c("Lat","Long","Pentane","Propane")
  
  par1sh.spt <- dat[,-4];
  par2sh.spt <- dat[,-3];
  
  par1sh.spt[3] <- rnorm(par1sh.spt$Pentane, mean = 0, sd = 0.5)
  par2sh.spt[3] <- rnorm(par2sh.spt$Propane, mean = 0, sd = 0.5)
  
  
  coordinates(par2sh.spt) <- ~ Long + Lat;
  coordinates(par1sh.spt) <- ~ Long + Lat;
  
  par2sh.grd <- Sobj_SpatialGrid(par2sh.spt, maxDim=(grd=grd))$SG;
  par1sh.grd <- Sobj_SpatialGrid(par1sh.spt, maxDim=(grd=grd))$SG;
  
  gridded(par2sh.grd) = TRUE
  gridded(par1sh.grd) = TRUE
  
  par2sh.idw <- krige(Propane~1, par2sh.spt, par2sh.grd, maxdist=(mxd=mxd), set=list(idp=cg));
  par1sh.idw <- krige(Pentane~1, par1sh.spt, par1sh.grd, maxdist=(mxd=mxd), set=list(idp=cs));
  
  par2sh.ras <- raster(par2sh.idw, layer=1,  values=TRUE);
  par1sh.ras <- raster(par1sh.idw, layer=1,  values=TRUE);
  
  # y - par1sh.ras
  # x - par2sh.ras
  
  subt <- overlay(par1sh.ras, par2sh.ras, fun=function(x,y){return(x-y)})
  summ <- overlay(par1sh.ras, par2sh.ras, fun=function(x,y){return(x+y)})
  #urx <- overlay(par1sh.ras, par2sh.ras, fun=function(x,y){return(x*y)})  
  #cellStats(uc, stat = "sum",na.rm = TRUE) 
  #unique(getValues(uc)) 
  
  par(mfrow = c(1, 2), xaxs = 'i', yaxs = 'i', mar = c(4,4,3,2.5),cex.lab = 1, cex.axis = 1.0)
  plot(summ, xlim = c(-86.152, -86.098), ylim = c(44.4, 44.43), ylab ='', axes = FALSE, xlab='' )
  box(col = 'grey')
  grid(col = 'grey', lty = 1)
  title(main="The sums of normalized values")
  axis(side = 1, tick = FALSE, col = 'grey', line = -0.8)
  axis(side = 2, tick = FALSE, col = 'grey', line = 0, las = 0)
  mtext(expression("Latitude"), side=2, line=2.5)
  mtext(expression("Longitude"), side=1, line=2.5)
  
  plot(subt, xlim = c(-86.152, -86.098), ylim = c(44.4, 44.43), ylab ='', axes = FALSE, xlab='' )
  box(col = 'grey')
  grid(col = 'grey', lty = 1)
  title(main="Differences of normalized values")
  axis(side = 1, tick = FALSE, col = 'grey', line = -0.8)
  axis(side = 2, tick = FALSE, col = 'grey', line = 0, las = 0)
  mtext(expression("Latitude"), side=2, line=2.5)
  mtext(expression("Longitude"), side=1, line=2.5)
 }
#####################################################
# lamin <- 44.40;
# lamax <- 44.43;
# lomin <- -86.147;
# lomax <- -86.098;
# lolim <- c(lomin, lomax);
# lalim <- c(lamin, lamax); 
################################################################################
#
# Creation of Pacage "Magnalys"
#
################################################################################

# package.skeleton(list=c("CCG", "CCS", "GMP", "PGS", "RCV", "R2DS", "R1DS", "SDP", "SPG",  "SPS","TMC"), name="Magnalys")

################################################################################
#
#
#
################################################################################
# xlim=c(-86.147, -86.098), ylim=c(44.4, 44.43)


param1.dat <- R1DS('Parameter-1.csv')
param2.dat <- R2DS('Parameter-2.csv')
LIMITS(3)
x <-3
SDP(eon, eox, ean, eax)
GMP("hybrid", 14, 1) #plot param #1
GMP("hybrid", 14, 2) #plot param #2
GMP("hybrid", 14, 3) #plot param #1 over param #2
GMP("hybrid", 14, 4) #plot param #2 over param #1
GMP("hybrid", 14, 5) #plot unit size spots #1 over param #2
GMP("hybrid", 14, 6) #plot unit size spots #2 over param #1


lamin <- 44.40;
lamax <- 44.43;
lomin <- -86.147;
lomax <- -86.098;
lolim <- c(lomin, lomax);
lalim <- c(lamin, lamax);
grid.dim <- 500





#cs <- CCS(10)
#cg <- CCG(10)
#SPS(500,10)
#SPG(500,10)


cs <- CCS(10)

cg <- CCG(10)

SPS(grid.dim,5)

SPG(grid.dim,10)

r1 <- cbind(c(-86.142, -86.131, -86.1, -86.1, -86.144,-86.144),
           c(44.42, 44.424, 44.421, 44.406, 44.407, 44.415))


r2 <- cbind(c(-86.142, -86.131, -86.1, -86.1, -86.144,-86.144),
            c(44.42, 44.424, 44.421, 44.406, 44.407, 44.415))

PGS(r1,r2)

param2.val <- RCV(grid.dim,10)

TMC(grid.dim,20)

