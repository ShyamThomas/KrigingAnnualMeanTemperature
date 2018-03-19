#####################################################################################################
##   Script for interpolating and mapping annual mean temperature for Ontario                      ##             
##   based on Envirnment Canada weather station data (in progress)                                 ##
##                                                                                                 ##
##   Date: 18th January 2018                                                                       ##
#####################################################################################################

library(ggplot2)
library("sp", lib.loc="~/R/win-library/3.3")
library(rgdal)
library(maptools)
library("gstat", lib.loc="~/R/win-library/3.3")
library("automap", lib.loc="~/R/win-library/3.3")
library(raster)
library(automap)

setwd("C:/Users/thoma/Rprojects/OntarioClimateInterpolation")

### Start by reading and mapping all station for a given year with daily weather data (at least 334 days of data was set asthreshold)###

OntMap = readShapePoly("1mprovnc.shp")
ont.stns2010=read.csv("OntarioTempData2010.csv")
head(ont.stns2010)
map=ggplot()+geom_polygon(data=OntMap, aes(x=long, y=lat, group=group))+
geom_point(data=ont.stns2010, aes(x=ont.stns2010$Max.of.Longitude, y=ont.stns2010$Max.of.Latitude, colour=ont.stns2010$Average.of.MeanTemp.C.2))
map+scale_color_gradient2(low="red", mid= "yellow", high="blue", midpoint=4.5)+theme(legend.position = c(0.9,0.85))+theme(legend.title=element_blank())

png("Ont2010_StnsMeanTemp.png", width=5, height = 5, units = "in", res = 900)




ont.lamb <- "+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

studyarea = readShapePoly("1mprovnc.shp")
proj4string(studyarea) = CRS(ont.lamb)
data2010.temp=read.csv("OntarioTempData2010.csv")
head(data2010.temp)
coordinates(data2010.temp) = cbind(data2010.temp$Max.of.Longitude, data2010.temp$Max.of.Latitude)
class(data2010.temp)
proj4string(data2010.temp) = CRS(ont.lamb)
crs(data2010.temp)


SA.ext=extent(studyarea)
gridsize=0.05

r = raster(SA.ext, res=0.05)
SA.r = rasterize(studyarea, r, fun="first")
plot(SA.r)
SA.grid = as(SA.r,"SpatialGridDataFrame")
plot(SA.grid)
points(data2010.temp,pch=20, col="red")

data2010.temp2=data2010.temp [-zerodist(data2010.temp)[,1],] ## remove duplicate x,y coordinates

### Fit a variogram and find the best variogram model to simulate the spatial variance###
var.temp=variogram(Average.of.MeanTemp.C.2~1, data2010.temp2)
summary(var.temp)
plot(var.temp)
exp.mod=vgm(psill=10, model="Exp", range=8, nugget=0)
vgm.fit = fit.variogram(object = var.temp, model = exp.mod)
plot(var.temp, model=vgm.fit)

### Using package automap###
temp.autoKrig=autoKrige(Average.of.MeanTemp.C.2~1, data2010.temp2, SA.grid)

### Try a universal variogram###
var.temp2=autofitVariogram(data2010.temp2$Average.of.MeanTemp.C.2 ~ data2010.temp2$Max.of.Longitude + data2010.temp2$Max.of.Latitude, data2010.temp2)
plot(var.temp2)
summary(var.temp2)



### Use the best variogram model to krig unsampled locations across Ontario ###
temp.krig=krige(Average.of.MeanTemp.C.2~1, data2010.temp2, SA.grid, model = exp.fit)
spplot(temp.krig,"var1.pred")
spplot(temp.krig,"var1.var")

##### Convert to raster and save in ascii format#######
temp.raster.pred=raster(temp.krig)

plot(temp.raster.pred)
temp.raster.var=raster(temp.krig, layer=2)
plot(temp.raster.var)

writeRaster(temp.raster, "TempPredictions", format = "ascii")


##################################################################################################################################################################
####################################################Extract Average Temperature of Ontario BsM Lakes##############################################################
OntLakes = readShapePoly("BsM_CleanedLakeBoundary.shp")
proj4string(OntLakes) = CRS(ont.lamb)
class(OntLakes)
plot(OntLakes)
##Overlay temperature raster with lake polygons##
plot(temp.raster.pred)
plot(studyarea, add=TRUE, border="black")
plot(OntLakes, border="white", add=TRUE)

meantmp.lakes=extract(temp.raster.pred,OntLakes, fun = mean, na.rm = TRUE, sp=TRUE)
head(meantmp.lakes)
