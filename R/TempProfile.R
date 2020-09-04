setwd("C:\\Users\\Martijn\\Dropbox\\PhD\\Panama2011_share\\PAR&Tmonitoring")

library(plyr)
library(scatterplot3d)
library(rgl)

data<-read.csv("ProfileForJMP.csv",sep=",", dec=".",header=TRUE)

palette(heat.colors(10))
plot3d(data$Temp,data$height, data$Time_24, col=data$Temp, size=3)
lines(data$Temp,data$height)