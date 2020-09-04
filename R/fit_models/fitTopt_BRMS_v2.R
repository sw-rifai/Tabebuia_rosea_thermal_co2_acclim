library(nls.multstart);
# library(rstan); 
library(tidyverse); library(bayesplot);
# library(rstantools)
library(rstan)
library(xfun); library(brms); 
options(mc.cores = parallel::detectCores()-1)
set.seed(1); 

params.all<-read_csv("data/params.all.28032018.csv")
params.all$Tk<-params.all$Tleaf+273.15

# create empty frame for Medlyn parameters
pars.med<-NULL

# treatments: (c.c, w.w, c.w, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
treat<-as.character(unique(params.all$Treat))

#---BRMS fit nonlinear VCmax
bfit_c.c <- brm(
  bf(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                       (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))), 
     kopt~1, 
     Ha~1, 
     # Hd~1, # Hd now fixed at 200 
     Topt~1,
     nl=T),
  data = params.all %>% filter(Treat=="c.c"),
  family = gaussian(link='identity'),
  chains = 3, 
  iter = 2000, 
  prior = c(prior(normal(60,5),nlpar="Ha",lb=25,ub=200), #25,150 
            prior(normal(250, 5), nlpar="kopt",lb=160,ub=300), # 230,300 
            prior(normal(200,1),nlpar="Hd",lb=190,ub=210), # 150,250
            prior(normal(315,5),nlpar="Topt",lb=300,ub=330)), 
  control=list(adapt_delta = 0.995, max_treedepth = 15))


bfit_c.w <- update(bfit_c.c, newdata=params.all %>% filter(Treat=='c.w'))
bfit_w.c <- update(bfit_c.c, newdata=params.all %>% filter(Treat=='w.c'))
bfit_w.w <- update(bfit_c.c, newdata=params.all %>% filter(Treat=='w.w'))


plot(Vcmax~Tk, data=params.all, col=as.factor(params.all$Treat),pch=20, xlim=c(299,322))

abline(v=fixef(bfit_c.c, robust = T)["Topt_Intercept","Estimate"],col='blue')
abline(v=fixef(bfit_c.c, robust = T)["Topt_Intercept","Q2.5"],col='blue',lty=3)
abline(v=fixef(bfit_c.c, robust = T)["Topt_Intercept","Q97.5"],col='blue',lty=3)

abline(v=fixef(bfit_w.c, robust = T)["Topt_Intercept","Estimate"],col='orchid')
abline(v=fixef(bfit_w.c, robust = T)["Topt_Intercept","Q2.5"],col='orchid',lty=3)
abline(v=fixef(bfit_w.c, robust = T)["Topt_Intercept","Q97.5"],col='orchid',lty=3)

abline(v=fixef(bfit_c.w, robust = T)["Topt_Intercept","Estimate"],col='darkgreen')
abline(v=fixef(bfit_c.w, robust = T)["Topt_Intercept","Q2.5"],col='darkgreen',lty=3)
abline(v=fixef(bfit_c.w, robust = T)["Topt_Intercept","Q97.5"],col='darkgreen',lty=3)

abline(v=fixef(bfit_w.w, robust = T)["Topt_Intercept","Estimate"],col='red')
abline(v=fixef(bfit_w.w, robust = T)["Topt_Intercept","Q2.5"],col='red',lty=3)
abline(v=fixef(bfit_w.w, robust = T)["Topt_Intercept","Q97.5"],col='red',lty=3)



p_w.w <- marginal_effects(bfit_w.w,spaghetti = 10)
p_w.w+geom_point(data=params.all, aes(Tk, Vcmax),col='blue')
