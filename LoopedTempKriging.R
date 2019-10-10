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

###################Creating a single interpolated surface for Ontario####

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


########################################################################################################################################
########################################################################################################################################

####Looping temperature datafiles to create multiple Kriged raasters######

library(sp)
library(rgdal)
library(maptools)
library(gstat)
library(raster)
library(automap)


setwd("C:/Users/thoma/Rprojects/OntarioClimateInterpolation")

ont.lamb <- "+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

ontborder = readShapePoly("1mprovnc.shp")
proj4string(ontborder) = CRS(ont.lamb)
class(ontborder)
plot(ontborder)

temp.files=list() #create an empty list
listcsv = dir(path= "./TempFiles", pattern = "*.csv") #set the path and file extension
listcsv # check!

newnames=c("LakeName", "Lat", "Long", "TotalTemp", "SumofDays", "AverageTemp") #assigning a better set of names for columns

## A simple for loop to read all the temperature csv files
for (k in 1:length(listcsv)){
temp.files[[k]] = read.csv(listcsv[k])
colnames(temp.files[[k]])=newnames
}
class(temp.files)
str(temp.files) #Check!


png("Ont2008_2012StnsMeanTemp.png", width=5, height = 5, units = "in", res = 900)

par(mfrow=c(2,3),mar=c(3,2.5,2,2.5))

temp.labels=seq(2008,2012,1)
rbPal = colorRampPalette(c('red','yellow','blue'))

## Now another loop to assign the new names and more...

for (i in 1:length(temp.files)){

coordinates(temp.files[[i]])=cbind(temp.files[[i]]$Long,temp.files[[i]]$Lat)
proj4string(temp.files[[i]]) = CRS(ont.lamb)
plot(ontborder, main=temp.labels[i], axes=TRUE)
temp.files[[i]]$Col= rbPal(12)[as.numeric(cut(temp.files[[i]]$AverageTemp, breaks=12))]
plot(temp.files[[i]], add=TRUE, pch=16, cex=0.75, col=temp.files[[i]]$Col, title=temp.labels[i])
}

dev.off()


###Now lets make an empty raster based on Ontario polygon and then grid it to make a kriging tempelate

ont.ext=extent(ontborder)
summary(ont.ext)

ont.raster.null=raster(ont.ext, res=0.05)
ont.raster.final=rasterize(ontborder,ont.raster.null, fun="first")
summary(ont.raster.final)
plot(ont.raster.final)

ont.grid=as(ont.raster.final, "SpatialGridDataFrame")
summary(ont.grid)
plot(ont.grid)


###For loop to reiterate the kriging process across the different data samples using the spatial grid of Ontario

temp.autoKrig=list()
temp.raster.pred=list()


png("Ont2008_2012AnnMeanTemp_rasters.png", width=5, height = 5, units = "in", res = 900)
par(mfrow=c(2,3),mar=c(3,2.5,2,2.5))

      for (i in 1:length(temp.files)){
          temp.autoKrig[[i]]=autoKrige(AverageTemp~1, temp.files[[i]], ont.grid)
          #plot(temp.autoKrig[[i]])
          temp.raster.pred[[i]]=raster(temp.autoKrig[[i]]$krige_output)
          plot(temp.raster.pred[[i]], main=temp.labels[i])
          }


###Now another loop to capture mean lake temperature for each year

OntLakes = readShapePoly("BsM_CleanedLakeBoundary.shp")
proj4string(OntLakes) = CRS(ont.lamb)
class(OntLakes)
plot(OntLakes)

meantmp.lakes=list()
par(mfrow=c(2,3),mar=c(3,2.5,2,2.5))

     for (i in 1:length(temp.raster.pred)){
          #plot(temp.raster.pred[[i]])
          #plot(studyarea, add=TRUE, border="black")
          #plot(OntLakes, border="white", add=TRUE)
          meantmp.lakes[[i]]=extract(temp.raster.pred[[i]],OntLakes, fun = mean, na.rm = TRUE, sp=TRUE)
     }

### final set of codes to extract the mean temp of all lakes for the 5 years

laketemps = matrix(ncol=length(meantmp.lakes), nrow=length(meantmp.lakes[[1]]$var1.pred))
for(i in 1:5){
laketemps[,i] =meantmp.lakes[[i]]$var1.pred
}
laketemps.df=data.frame(laketemps)
colnames(laketemps.df)=temp.labels
class(laketemps.df)
head(laketemps.df)

meantmp.lakes.df=data.frame(meantmp.lakes[[1]][,c(1,6:10,14:21)])
meantmp.lakes.final.df=cbind(meantmp.lakes.final.df,output.df)
class(meantmp.lakes.final.df)
head(meantmp.lakes.final.df)

meantmp.lakes.final.df$MeanTemp = rowMeans(meantmp.lakes.final.df[,17:19])
head(meantmp.lakes.final.df)
write.csv(meantmp.lakes.final.df, "meantemp.ontlakes.2008_2012.csv")
