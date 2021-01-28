library(nls.multstart); library(rstan); library(tidyverse); library(bayesplot);
library(rstantools)
options(mc.cores = parallel::detectCores())
set.seed(1); 

params.all<-read_csv("data/params.all.28032018.csv")

params.all$Tk<-params.all$Tleaf+273.15

# create empty frame for Medlyn parameters
pars.med<-NULL

# removed for loop 
# treatments: (c.c, w.w, c.w, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
treat<-as.character(unique(params.all$Treat))


#########################################################################################
# FIT TREATMENT 1
#########################################################################################
idx <- 1
tmp_dat <- subset(params.all, Treat==treat[idx])
fit_med <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                           (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                         data = tmp_dat,
                         iter = 1000,
                         start_lower = c(kopt = 160, Hd = 100, Ha = 50, Topt = 300),
                         start_upper = c(kopt = 300, Hd = 900, Ha = 300, Topt = 330),
                         supp_errors = 'Y',
                         na.action = na.omit,
                         #convergence_count = 500,
                         lower = c(kopt = 150, Hd = 50, Ha = 10, Topt = 300))
fit_med %>% summary()
tmp_dat %>% plot(Vcmax~Tk, data=.)

file_path <- "R/FitTopt.stan";
lines <- readLines(file_path, encoding="ASCII");
for (n in 1:length(lines)) cat(lines[n],'\n');

m_stan1 <- stan(model_code = lines, #"R/StanSim/meanOfPreds.stan",
               data = list(Vcmax=tmp_dat$Vcmax,
                           Tk=tmp_dat$Tk,
                           N = dim(tmp_dat)[1]),
               warmup = 2000,
               iter=4000, thin=1, chains=3, verbose=T, 
               control=list(adapt_delta=0.9), 
               algorithm = "HMC")
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(m_stan1, digits=3, pars=c("kopt","Hd","Ha","Topt","sigma"))
mcmc_trace(as.array(m_stan1), c("sigma","Topt","Hd","Ha","kopt"))
mcmc_dens(as.array(m_stan1), c("sigma","Topt","Hd","Ha","kopt"))


# # alpha=shape; rate=beta
# # u=alpha/beta
# # variance alpha/beta^2
# curve(dgamma(x, shape = 30, rate = 0.1),0,1000)


post_df <- rstan::extract(m_stan1)

r_Vcmax <- function(kopt, Hd, Ha, Tk, Topt){
  Vcmax <- kopt * ( (Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                      (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))    ); 
  return(Vcmax); 
}

df_med <- lapply(post_df,"median") %>% as.data.frame()
df_5 <- lapply(post_df,"quantile",0.05) %>% as.data.frame();
df_95 <- lapply(post_df,"quantile",0.95) %>% as.data.frame();
df_ppd <- tibble(Tk=290:320,Vcmax05=NA,Vcmax50=NA,Vcmax95=NA);

for(i in 1:dim(df_ppd)[1]){
  df_ppd$Vcmax05[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.05, na.rm=T)
  df_ppd$Vcmax50[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.5, na.rm=T)
  df_ppd$Vcmax95[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.95, na.rm=T)
}


png(paste0("figures/",treat[idx],".png"),res = 300, pointsize=16, width = 20, height=15, units="cm")
col_pal <- viridis::inferno(5, begin = 0.2, end = 0.8)
curve(r_Vcmax(post_df$kopt[1], post_df$Hd[1], post_df$Ha[1], Tk, post_df$Topt[1]), 
      xname = "Tk",from = 290,to=320, 
      ylim=c(40,max(tmp_dat$Vcmax)),
      col=rgb(0,0,0,alpha=0.01), 
      ylab="Vcmax",xlab="Temperature [K]",main=treat[idx])
for(i in 1:1000){
  curve(r_Vcmax(post_df$kopt[i], post_df$Hd[i], post_df$Ha[i], Tk, post_df$Topt[i]), 
        xname = "Tk",from = 290,to=320,add=T,col=rgb(0,0,0,alpha=0.01))
}
lines(Vcmax50~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=3)
lines(Vcmax05~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=1, lty=3)
lines(Vcmax95~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=1, lty=3)
points(Vcmax~Tk, data=tmp_dat, col=col_pal[4],pch=20)
abline(v=post_df$Topt %>% median(), col=col_pal[5], lwd=2)
abline(v=post_df$Topt %>% quantile(0.05), col=col_pal[5], lwd=1, lty=3)
abline(v=post_df$Topt %>% quantile(0.95), col=col_pal[5], lwd=1, lty=3)
dev.off()
#########################################################################################
# END TREATMENT 1
#########################################################################################


#########################################################################################
# FIT TREATMENT 2
#########################################################################################
idx <- 2
tmp_dat <- subset(params.all, Treat==treat[idx])
fit_med <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                           (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                         data = tmp_dat,
                         iter = 1000,
                         start_lower = c(kopt = 160, Hd = 100, Ha = 50, Topt = 300),
                         start_upper = c(kopt = 300, Hd = 900, Ha = 300, Topt = 330),
                         supp_errors = 'Y',
                         na.action = na.omit,
                         #convergence_count = 500,
                         lower = c(kopt = 150, Hd = 50, Ha = 10, Topt = 300))
fit_med %>% summary()
tmp_dat %>% plot(Vcmax~Tk, data=.)

file_path <- "R/FitTopt.stan";
lines <- readLines(file_path, encoding="ASCII");
for (n in 1:length(lines)) cat(lines[n],'\n');

m_stan2 <- stan(model_code = lines, #"R/StanSim/meanOfPreds.stan",
               data = list(Vcmax=tmp_dat$Vcmax,
                           Tk=tmp_dat$Tk,
                           N = dim(tmp_dat)[1]),
               warmup = 2000,
               iter=4000, thin=1, chains=3, verbose=T, 
               control=list(adapt_delta=0.9), 
               algorithm = "NUTS")
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(m_stan2, digits=3, pars=c("kopt","Hd","Ha","Topt","sigma"))
mcmc_trace(as.array(m_stan2), c("sigma","Topt","Hd","Ha","kopt"))
mcmc_dens(as.array(m_stan2), c("sigma","Topt","Hd","Ha","kopt"))

post_df <- rstan::extract(m_stan2)

r_Vcmax <- function(kopt, Hd, Ha, Tk, Topt){
  Vcmax <- kopt * ( (Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                      (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))    ); 
  return(Vcmax); 
}

df_med <- lapply(post_df,"median") %>% as.data.frame()
df_5 <- lapply(post_df,"quantile",0.05) %>% as.data.frame();
df_95 <- lapply(post_df,"quantile",0.95) %>% as.data.frame();
df_ppd <- tibble(Tk=290:320,Vcmax05=NA,Vcmax50=NA,Vcmax95=NA);

for(i in 1:dim(df_ppd)[1]){
  df_ppd$Vcmax05[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.05, na.rm=T)
  df_ppd$Vcmax50[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.5, na.rm=T)
  df_ppd$Vcmax95[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.95, na.rm=T)
}


png(paste0("figures/",treat[idx],".png"),res = 300, pointsize=16, width = 20, height=15, units="cm")
col_pal <- viridis::inferno(5, begin = 0.2, end = 0.8)
curve(r_Vcmax(post_df$kopt[1], post_df$Hd[1], post_df$Ha[1], Tk, post_df$Topt[1]), 
      xname = "Tk",from = 290,to=320, 
      ylim=c(40,max(tmp_dat$Vcmax)),
      col=rgb(0,0,0,alpha=0.01), 
      ylab="Vcmax",xlab="Temperature [K]",main=treat[idx])
for(i in 1:1000){
  curve(r_Vcmax(post_df$kopt[i], post_df$Hd[i], post_df$Ha[i], Tk, post_df$Topt[i]), 
        xname = "Tk",from = 290,to=320,add=T,col=rgb(0,0,0,alpha=0.01))
}
lines(Vcmax50~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=3)
lines(Vcmax05~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=1, lty=3)
lines(Vcmax95~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=1, lty=3)
points(Vcmax~Tk, data=tmp_dat, col=col_pal[4],pch=20)
abline(v=post_df$Topt %>% median(), col=col_pal[5], lwd=2)
abline(v=post_df$Topt %>% quantile(0.05), col=col_pal[5], lwd=1, lty=3)
abline(v=post_df$Topt %>% quantile(0.95), col=col_pal[5], lwd=1, lty=3)
dev.off()
#########################################################################################
# END TREATMENT 2
#########################################################################################


#########################################################################################
# FIT TREATMENT 3 - THIS ONE MIGHT NEED SPECIAL PRIORS
#########################################################################################
idx <- 3
tmp_dat <- subset(params.all, Treat==treat[idx])
fit_med <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                           (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                         data = tmp_dat,
                         iter = 1000,
                         start_lower = c(kopt = 160, Hd = 100, Ha = 50, Topt = 300),
                         start_upper = c(kopt = 300, Hd = 900, Ha = 300, Topt = 330),
                         supp_errors = 'Y',
                         na.action = na.omit,
                         #convergence_count = 500,
                         lower = c(kopt = 150, Hd = 50, Ha = 10, Topt = 300))
fit_med %>% summary()
tmp_dat %>% plot(Vcmax~Tk, data=.)

file_path <- "R/FitTopt.stan";
lines <- readLines(file_path, encoding="ASCII");
for (n in 1:length(lines)) cat(lines[n],'\n');

set.seed(2)
m_stan3 <- stan(model_code = lines, #"R/StanSim/meanOfPreds.stan",
               data = list(Vcmax=tmp_dat$Vcmax,
                           Tk=tmp_dat$Tk,
                           N = dim(tmp_dat)[1]),
               warmup = 2000,
               iter=4000, thin=1, chains=3, verbose=T, 
               control=list(adapt_delta=0.95), 
               algorithm = "NUTS")
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(m_stan3, digits=3, pars=c("kopt","Hd","Ha","Topt","sigma"))
mcmc_trace(as.array(m_stan3), c("sigma","Topt","Hd","Ha","kopt"))
mcmc_dens(as.array(m_stan3), c("sigma","Topt","Hd","Ha","kopt"))

post_df <- rstan::extract(m_stan3)

r_Vcmax <- function(kopt, Hd, Ha, Tk, Topt){
  Vcmax <- kopt * ( (Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                      (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))    ); 
  return(Vcmax); 
}

df_med <- lapply(post_df,"median") %>% as.data.frame()
df_5 <- lapply(post_df,"quantile",0.05) %>% as.data.frame();
df_95 <- lapply(post_df,"quantile",0.95) %>% as.data.frame();
df_ppd <- tibble(Tk=290:320,Vcmax05=NA,Vcmax50=NA,Vcmax95=NA);

for(i in 1:dim(df_ppd)[1]){
  df_ppd$Vcmax05[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.05, na.rm=T)
  df_ppd$Vcmax50[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.5, na.rm=T)
  df_ppd$Vcmax95[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.95, na.rm=T)
}


png(paste0("figures/",treat[idx],".png"),res = 300, pointsize=16, width = 20, height=15, units="cm")
col_pal <- viridis::inferno(5, begin = 0.2, end = 0.8)
curve(r_Vcmax(post_df$kopt[1], post_df$Hd[1], post_df$Ha[1], Tk, post_df$Topt[1]), 
      xname = "Tk",from = 290,to=320, 
      ylim=c(40,max(tmp_dat$Vcmax)),
      col=rgb(0,0,0,alpha=0.01), 
      ylab="Vcmax",xlab="Temperature [K]",main=treat[idx])
for(i in 1:1000){
  curve(r_Vcmax(post_df$kopt[i], post_df$Hd[i], post_df$Ha[i], Tk, post_df$Topt[i]), 
        xname = "Tk",from = 290,to=320,add=T,col=rgb(0,0,0,alpha=0.01))
}
lines(Vcmax50~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=3)
lines(Vcmax05~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=1, lty=3)
lines(Vcmax95~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=1, lty=3)
points(Vcmax~Tk, data=tmp_dat, col=col_pal[4],pch=20)
abline(v=post_df$Topt %>% median(), col=col_pal[5], lwd=2)
abline(v=post_df$Topt %>% quantile(0.05), col=col_pal[5], lwd=1, lty=3)
abline(v=post_df$Topt %>% quantile(0.95), col=col_pal[5], lwd=1, lty=3)
dev.off()
#########################################################################################
# END TREATMENT 3
#########################################################################################


#########################################################################################
# FIT TREATMENT 4
#########################################################################################
idx <- 4
tmp_dat <- subset(params.all, Treat==treat[idx])
fit_med <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                           (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                         data = tmp_dat,
                         iter = 1000,
                         start_lower = c(kopt = 160, Hd = 100, Ha = 50, Topt = 300),
                         start_upper = c(kopt = 300, Hd = 900, Ha = 300, Topt = 330),
                         supp_errors = 'Y',
                         na.action = na.omit,
                         #convergence_count = 500,
                         lower = c(kopt = 150, Hd = 50, Ha = 10, Topt = 300))
fit_med %>% summary()
tmp_dat %>% plot(Vcmax~Tk, data=.)

file_path <- "R/FitTopt.stan";
lines <- readLines(file_path, encoding="ASCII");
for (n in 1:length(lines)) cat(lines[n],'\n');

m_stan4 <- stan(model_code = lines, #"R/StanSim/meanOfPreds.stan",
               data = list(Vcmax=tmp_dat$Vcmax,
                           Tk=tmp_dat$Tk,
                           N = dim(tmp_dat)[1]),
               warmup = 2000,
               iter=4000, thin=1, chains=3, verbose=T, 
               control=list(adapt_delta=0.9), 
               algorithm = "NUTS")
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(m_stan4, digits=3, pars=c("kopt","Hd","Ha","Topt","sigma"))
mcmc_trace(as.array(m_stan4), c("sigma","Topt","Hd","Ha","kopt"))
mcmc_dens(as.array(m_stan4), c("sigma","Topt","Hd","Ha","kopt"))

post_df <- rstan::extract(m_stan4)

r_Vcmax <- function(kopt, Hd, Ha, Tk, Topt){
  Vcmax <- kopt * ( (Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                      (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))    ); 
  return(Vcmax); 
}

df_med <- lapply(post_df,"median") %>% as.data.frame()
df_5 <- lapply(post_df,"quantile",0.05) %>% as.data.frame();
df_95 <- lapply(post_df,"quantile",0.95) %>% as.data.frame();
df_ppd <- tibble(Tk=290:320,Vcmax05=NA,Vcmax50=NA,Vcmax95=NA);

for(i in 1:dim(df_ppd)[1]){
  df_ppd$Vcmax05[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.05, na.rm=T)
  df_ppd$Vcmax50[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.5, na.rm=T)
  df_ppd$Vcmax95[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.95, na.rm=T)
}


png(paste0("figures/",treat[idx],".png"),res = 300, pointsize=16, width = 20, height=15, units="cm")
col_pal <- viridis::inferno(5, begin = 0.2, end = 0.8)
curve(r_Vcmax(post_df$kopt[1], post_df$Hd[1], post_df$Ha[1], Tk, post_df$Topt[1]), 
      xname = "Tk",from = 290,to=320, 
      ylim=c(40,max(tmp_dat$Vcmax)),
      col=rgb(0,0,0,alpha=0.01), 
      ylab="Vcmax",xlab="Temperature [K]",main=treat[idx])
for(i in 1:1000){
  curve(r_Vcmax(post_df$kopt[i], post_df$Hd[i], post_df$Ha[i], Tk, post_df$Topt[i]), 
        xname = "Tk",from = 290,to=320,add=T,col=rgb(0,0,0,alpha=0.01))
}
lines(Vcmax50~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=3)
lines(Vcmax05~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=1, lty=3)
lines(Vcmax95~Tk, data=df_ppd, type="l",col=col_pal[2],lwd=1, lty=3)
points(Vcmax~Tk, data=tmp_dat, col=col_pal[4],pch=20)
abline(v=post_df$Topt %>% median(), col=col_pal[5], lwd=2)
abline(v=post_df$Topt %>% quantile(0.05), col=col_pal[5], lwd=1, lty=3)
abline(v=post_df$Topt %>% quantile(0.95), col=col_pal[5], lwd=1, lty=3)
dev.off()
#########################################################################################
# END TREATMENT 4
#########################################################################################


m_stan1  %>% print(par="Topt")
m_stan2  %>% print(par="Topt") 
m_stan3  %>% print(par="Topt") 
m_stan4  %>% print(par="Topt")

m_stan1@inits
# pars.med.f<-as.data.frame(pars.med)
# colnames(pars.med.f)<-list("Treat","Vcmax.opt","Hd.Vcmax","Ha.Vcmax","Topt.Vcmax",
#                            "SEM.Vcmax.opt","SEM.Hd.Vcmax","SEM.Ha.Vcmax","SEM.Topt.Vcmax")
# 
# for (i in 2:length(pars.med.f))
# {
#   pars.med.f[,c(i)]<-as.numeric(as.character(pars.med.f[,c(i)]))
# }
