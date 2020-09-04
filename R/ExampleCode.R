devtools::install_github('padpadpadpad/nls.multstart')
library(nls.multstart) 


params.all<-read.csv("params.all.28032018.csv")

params.all$Tk<-params.all$Tleaf+273.15

# create empty frame for Medlyn parameters
pars.med<-NULL

# try to do so in a for loop by treatment (c.c, w.w, c.w, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)

treat<-as.character(unique(params.all$Treat))

for (i in 1:length(treat))
{
  dat2<-subset(params.all, Treat==treat[i])
  
  fit_med <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                        (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                        data = dat2,
                        iter = 1000,
                        start_lower = c(kopt = 160, Hd = 100, Ha = 50, Topt = 300),
                        start_upper = c(kopt = 300, Hd = 900, Ha = 300, Topt = 330),
                        supp_errors = 'Y',
                        na.action = na.omit,
                        #convergence_count = 500,
                        lower = c(kopt = 150, Hd = 50, Ha = 10, Topt = 300))
  results1<-summary(fit_med)
  results2<-as.data.frame(results1$parameters)
  pars.med<-rbind(pars.med, c(treat[i], results2[1,1], results2[2,1], results2[3,1], results2[4,1],
                                        results2[1,2], results2[2,2], results2[3,2], results2[4,2]))
  
  plot(dat2$Vcmax~dat2$Tk, ylim=c(50,325), xlim=c(300,320))
  x1<-seq(300,317, by=0.5)
  y1<-results2[1,1] * ((results2[2,1] * (2.718282^((results2[3,1]*(x1-results2[4,1]))/(x1*0.008314*results2[4,1])))) / 
                (results2[2,1] - (results2[3,1]*(1-(2.718282^((results2[2,1]*(x1-results2[4,1]))/(x1*0.008314*results2[4,1])))))))
  points(x1,y1, type="l", col="red")
}

pars.med.f<-as.data.frame(pars.med)
colnames(pars.med.f)<-list("Treat","Vcmax.opt","Hd.Vcmax","Ha.Vcmax","Topt.Vcmax",
                           "SEM.Vcmax.opt","SEM.Hd.Vcmax","SEM.Ha.Vcmax","SEM.Topt.Vcmax")

for (i in 2:length(pars.med.f))
{
  pars.med.f[,c(i)]<-as.numeric(as.character(pars.med.f[,c(i)]))
}
