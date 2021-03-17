library(tidyverse)
# Load data & clean dataframes-------------------------------------------


params.all<-read.csv("data/params.all.28032018.csv")
# remove duplicate columns that will mess up merging with other df that has the same columns
pars.core <- params.all[,-c(2,6,7,11:21),drop=FALSE] 

P400<-read.csv("data/P400.150321.csv")
# remove duplicate columns that will mess up merging with other df that has the same columns
P400.core <- P400[,-c(7,9,10,13,15,18,20,26:31,41:71),drop=FALSE] 

P800<-read.csv("data/P800.150321.csv")
# remove duplicate columns that will mess up merging with other df that has the same columns
P800.core <- P800[,-c(5,17,19,23:29,39:70),drop=FALSE] 

stom.lim.400<-merge(P400.core, pars.core, by=c("Treat", "Plant", "Leaf", "Rep"))  
stom.lim.800<-merge(P800.core, pars.core, by=c("Treat", "Plant", "Leaf", "Rep"))  


# Calculate stomatal limitation at 400 ppm --------------------------------

# calculate Kc following Bernacci et al 2001 (PCE). R=8.314 J mol-1K-1
stom.lim.400$Kc<-404.9*(exp(1))^((79403*(stom.lim.400$Tk-298.15))/(298.15*8.314*stom.lim.400$Tk))

# calculate K0 following Bernacci et al 2001 (PCE). R=8.314 J mol-1K-1
stom.lim.400$K0<-278.4*(exp(1))^((36380*(stom.lim.400$Tk-298.15))/(298.15*8.314*stom.lim.400$Tk))

### calculate gamma* following Bernacci et al 2001 (PCE).
stom.lim.400$g.star<-42.75*(exp(1))^((37830*(stom.lim.400$Tk-298.15))/(298.15*8.314*stom.lim.400$Tk))

# from this calculate Km. One imput here is Oi, or the intercellular oxygen concentration. This is assumed to be 210 at sea level (de Kauwe et al 2016)
stom.lim.400$Km<-stom.lim.400$Kc*(1+(210/stom.lim.400$K0))

# Calculate what Anet would have been at each point had CO2 (or gs) not been limiting by setting ci to ca. 
# and then using Farquhar et al 1980 to calculate A.gross 
stom.lim.400$Ag<-(stom.lim.400$Vcmax*(stom.lim.400$CO2R-stom.lim.400$g.star))/(stom.lim.400$CO2R+stom.lim.400$Km)

#### now turn that into Anet using the fact that Rday was assumed to be 1.5% of Vcmax
stom.lim.400$An.mod<-stom.lim.400$Ag-0.015*stom.lim.400$Vcmax

# now determine the extent of stomatal restriction or "l"
stom.lim.400$l<-1-(stom.lim.400$Photo/stom.lim.400$An.mod)




# Calculate stomatal limitation at 800 ppm --------------------------------

stom.lim.800$Kc<-404.9*(exp(1))^((79403*(stom.lim.800$Tk-298.15))/(298.15*8.314*stom.lim.800$Tk))
stom.lim.800$K0<-278.4*(exp(1))^((36380*(stom.lim.800$Tk-298.15))/(298.15*8.314*stom.lim.800$Tk))
stom.lim.800$g.star<-42.75*(exp(1))^((37830*(stom.lim.800$Tk-298.15))/(298.15*8.314*stom.lim.800$Tk))
stom.lim.800$Km<-stom.lim.800$Kc*(1+(210/stom.lim.800$K0))
stom.lim.800$Ag<-(stom.lim.800$Vcmax*(stom.lim.800$CO2R-stom.lim.800$g.star))/(stom.lim.800$CO2R+stom.lim.800$Km)
stom.lim.800$An.mod<-stom.lim.800$Ag-0.015*stom.lim.800$Vcmax
stom.lim.800$l<-1-(stom.lim.800$Photo/stom.lim.800$An.mod)


# plot stomatal limitation factor 'l' against temperature -----------------
# plot l at 400 ppm for plants at control conditions and at 800 ppm at treatment CO2
png(filename="figures/Fig6_stomatalLimitation.png", width = 600, height =600)
par(mfrow=c(2,2), mar=c(0.5,0.5,0,0), oma=c(4,4,0,1))
par(tcl=0.5) #TO PUT TICK MARKS ON THE INSIDE
mgp=c(2, 0.5, 0)

gray.trans <- "#55555555"
plot(stom.lim.400$l[stom.lim.400$Treat=="c.c"]~stom.lim.400$Tleaf[stom.lim.400$Treat=="c.c"], las=1, 
     ylim=c(0.1,0.9), xlim=c(25,45), pch=19,cex=2, xaxt='n',
     col=gray.trans,
     cex.axis=1.5, ylab = 'Stomatal Limitation')
axis(1, at=c(25, 30, 35, 40, 45), labels=F)
text(33,0.85, "(a) Control", font=3,family="serif", cex=1.7)
text(33,0.75, "400 ppm", font=3,family="serif", cex=1.7)


plot(stom.lim.800$l[stom.lim.800$Treat=="c.w"]~stom.lim.800$Tleaf[stom.lim.800$Treat=="c.w"], 
     las=1, ylim=c(0.1,0.9), xlim=c(25,45), pch=19,cex=2, xaxt='n',yaxt='n',col=gray.trans, cex.axis=1.5)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8), labels=F)
axis(1, at=c(25, 30, 35, 40, 45), labels=F)
text(33,0.85, expression(paste("(b) ",Control %->% Treatment)), font=3,family="serif", cex=1.7)
text(33,0.75, "800 ppm", font=3,family="serif", cex=1.7)


plot(stom.lim.800$l[stom.lim.800$Treat=="w.w"]~stom.lim.800$Tleaf[stom.lim.800$Treat=="w.w"], 
     las=1,
     ylim=c(0.1,0.9), xlim=c(25,45), 
     pch=19,cex=2, col=gray.trans, cex.axis=1.5, 
     xlab=expression(paste("Leaf temperature"~degree*C)))
text(33,0.85, "(c) Treatment", font=3,family="serif", cex=1.7)
text(33,0.75, "800 ppm", font=3,family="serif", cex=1.7)

plot(stom.lim.400$l[stom.lim.400$Treat=="w.c"]~stom.lim.400$Tleaf[stom.lim.400$Treat=="w.c"], 
     las=1, ylim=c(0.1,0.9), xlim=c(25,45), pch=19,cex=2, col=gray.trans, yaxt='n', cex.axis=1.5)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8), labels=F)
text(33,0.85, expression(paste("(d) ",Treatment %->% Control)), font=3,family="serif", cex=1.7)
text(33,0.75, "400 ppm", font=3,family="serif", cex=1.7)
dev.off()
