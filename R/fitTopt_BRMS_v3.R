library(nls.multstart); library(rstan); library(tidyverse)
options(mc.cores = parallel::detectCores()-1)
set.seed(1); 

params.all<-read_csv("data/params.all.28032018.csv")
params.all$Tk<-params.all$Tleaf+273.15

# create empty frame for Medlyn parameters
pars.med<-NULL

# treatments: (c.c, w.w, c.w, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
treat<-as.character(unique(params.all$Treat))

#NLS multistart
nfit_c.c <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                           (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                         data = params.all %>% filter(Treat=='c.c'),
                         iter = 1000,
                         start_lower = c(kopt = 160, Hd = 200, Ha = 50, Topt = 300),
                         start_upper = c(kopt = 300, Hd = 200, Ha = 300, Topt = 330),
                         supp_errors = 'Y',
                         na.action = na.omit,
                         #convergence_count = 500,
                         upper = c(kopt=300, Hd=200, Ha=150, Topt=330), 
                         lower = c(kopt = 150, Hd = 200, Ha = 10, Topt = 273.15))
nfit_c.w <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                            (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                          data = params.all %>% filter(Treat=='c.w'),
                          iter = 1000,
                          start_lower = c(kopt = 160, Hd = 200, Ha = 50, Topt = 300),
                          start_upper = c(kopt = 300, Hd = 200, Ha = 300, Topt = 330),
                          supp_errors = 'Y',
                          na.action = na.omit,
                          #convergence_count = 500,
                          upper = c(kopt=300, Hd=200, Ha=150, Topt=330), 
                          lower = c(kopt = 150, Hd = 200, Ha = 10, Topt = 273.15))
nfit_w.c <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                            (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                          data = params.all %>% filter(Treat=='w.c'),
                          iter = 1000,
                          start_lower = c(kopt = 160, Hd = 200, Ha = 50, Topt = 300),
                          start_upper = c(kopt = 300, Hd = 200, Ha = 300, Topt = 330),
                          supp_errors = 'Y',
                          na.action = na.omit,
                          #convergence_count = 500,
                          upper = c(kopt=300, Hd=200, Ha=150, Topt=330), 
                          lower = c(kopt = 150, Hd = 200, Ha = 10, Topt = 273.15))
nfit_w.w <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                            (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                          data = params.all %>% filter(Treat=='w.w'),
                          iter = 1000,
                          start_lower = c(kopt = 160, Hd = 200, Ha = 50, Topt = 300),
                          start_upper = c(kopt = 300, Hd = 200, Ha = 300, Topt = 330),
                          supp_errors = 'Y',
                          na.action = na.omit,
                          #convergence_count = 500,
                          upper = c(kopt=300, Hd=200, Ha=150, Topt=330), 
                          lower = c(kopt = 150, Hd = 200, Ha = 10, Topt = 273.15))


# bfit_c.c_Hd200 <- brm(
#   bf(Vcmax ~ kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
#                        (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))), 
#      kopt~1, 
#      Ha~1, 
#      # Hd~1, # Hd now fixed at 200 
#      Topt~1,
#      nl=T),
#   data = params.all %>% filter(Treat=="c.c"),
#   family = gaussian(link='identity'),
#   chains = 3, 
#   iter = 20000, thin=20,
#   prior = c(prior(normal(60,5),nlpar="Ha",lb=25,ub=200), #25,150 
#             prior(normal(250, 5), nlpar="kopt",lb=50,ub=500), # 230,300 
#             # prior(normal(200,1),nlpar="Hd",lb=190,ub=210), # 150,250
#             prior(normal(315,5),nlpar="Topt",lb=290,ub=330)), 
#   control=list(adapt_delta = 0.995, max_treedepth = 15))
# bfit_c.w_Hd200 <- brm(
#   bf(Vcmax ~ kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
#                        (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))), 
#      kopt~1, 
#      Ha~1, 
#      # Hd~1, # Hd now fixed at 200 
#      Topt~1,
#      nl=T),
#   data = params.all %>% filter(Treat=="c.w"),
#   family = gaussian(link='identity'),
#   chains = 3, 
#   iter = 20000, thin=20,
#   prior = c(prior(normal(60,5),nlpar="Ha",lb=25,ub=200), #25,150 
#             prior(normal(250, 5), nlpar="kopt",lb=50,ub=500), # 230,300 
#             # prior(normal(200,1),nlpar="Hd",lb=190,ub=210), # 150,250
#             prior(normal(315,5),nlpar="Topt",lb=290,ub=330)), 
#   control=list(adapt_delta = 0.995, max_treedepth = 15))
# bfit_w.c_Hd200 <- brm(
#   bf(Vcmax ~ kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
#                        (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))), 
#      kopt~1, 
#      Ha~1, 
#      # Hd~1, # Hd now fixed at 200 
#      Topt~1,
#      nl=T),
#   data = params.all %>% filter(Treat=="w.c"),
#   family = gaussian(link='identity'),
#   chains = 3, 
#   iter = 20000, thin=20,
#   prior = c(prior(normal(60,5),nlpar="Ha",lb=25,ub=200), #25,150 
#             prior(normal(250, 5), nlpar="kopt",lb=50,ub=500), # 230,300 
#             # prior(normal(200,1),nlpar="Hd",lb=190,ub=210), # 150,250
#             prior(normal(315,5),nlpar="Topt",lb=290,ub=330)), 
#   control=list(adapt_delta = 0.995, max_treedepth = 15))
# bfit_w.w_Hd200 <- brm(
#   bf(Vcmax ~ kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
#                        (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))), 
#      kopt~1, 
#      Ha~1, 
#      # Hd~1, # Hd now fixed at 200 
#      Topt~1,
#      nl=T),
#   data = params.all %>% filter(Treat=="w.w"),
#   family = gaussian(link='identity'),
#   chains = 3, 
#   iter = 20000, thin=20,
#   prior = c(prior(normal(60,5),nlpar="Ha",lb=25,ub=200), #25,150 
#             prior(normal(250, 5), nlpar="kopt",lb=50,ub=500), # 230,300 
#             # prior(normal(200,1),nlpar="Hd",lb=190,ub=210), # 150,250
#             prior(normal(315,5),nlpar="Topt",lb=290,ub=330)), 
#   control=list(adapt_delta = 0.995, max_treedepth = 15))
# 
# # bfit_c.w_Hd200 <- update(bfit_c.c_Hd200, newdata=params.all %>% filter(Treat=="c.w"))
# # bfit_w.c_Hd200 <- update(bfit_c.c_Hd200, newdata=params.all %>% filter(Treat=="w.c"))
# # bfit_w.w_Hd200 <- update(bfit_c.c_Hd200, newdata=params.all %>% filter(Treat=="w.w"))
# 
# 
# # Plot Predicted Fits Across Treatments -----------------------------------
# png("figures/CompareToptEst_NLS_BRMS.png", height = 800, width=1000, units='px', pointsize = 18)
# sim_Tk <- seq(298,327,0.1)
# params.all %>% 
#   plot(Vcmax~Tk, data=., xlim=c(299,322), main="Topt");
# lines(posterior_predict(bfit_c.c_Hd200, transform = colMeans, 
#                         newdata = data.frame(Tk=sim_Tk))~sim_Tk, 
#       col='black')
# lines(posterior_predict(bfit_c.w_Hd200, transform = colMeans, 
#                         newdata = data.frame(Tk=sim_Tk))~sim_Tk, 
#       col='orange')
# lines(posterior_predict(bfit_w.c_Hd200, transform = colMeans, 
#                         newdata = data.frame(Tk=sim_Tk))~sim_Tk, 
#       col='orchid')
# lines(posterior_predict(bfit_w.w_Hd200, transform = colMeans, 
#                         newdata = data.frame(Tk=sim_Tk))~sim_Tk, 
#       col='red')
# abline(v=fixef(bfit_c.c_Hd200, robust = T)["Topt_Intercept","Estimate"],col='black',lwd=2)
# abline(v=fixef(bfit_c.w_Hd200, robust = T)["Topt_Intercept","Estimate"],col='orange',lwd=2)
# abline(v=fixef(bfit_w.c_Hd200, robust = T)["Topt_Intercept","Estimate"],col='orchid',lwd=2)
# abline(v=fixef(bfit_w.w_Hd200, robust = T)["Topt_Intercept","Estimate"],col='red',lwd=2)
# 
# 
# abline(v=summary(nfit_c.c)$coefficients["Topt","Estimate"],col='black',lty=3,lwd=2)
# abline(v=summary(nfit_c.w)$coefficients["Topt","Estimate"],col='orange',lty=3,lwd=2)
# abline(v=summary(nfit_w.c)$coefficients["Topt","Estimate"],col='orchid',lty=3,lwd=2)
# abline(v=summary(nfit_w.w)$coefficients["Topt","Estimate"],col='red',lty=3,lwd=2)
# 
# legend("topleft", col=c("black","orange" ,"orchid","red"), lty=c(1,1,1,1), 
#        legend=c("c.c","c.w","w.c","w.w"))
# legend('bottomright', legend=c('Bayesian Fit', 'NLS Fit'), lty=c(1,3))
# dev.off()


# Rstan fit ---------------------------------------------------------------
library(rstan)
stan_c.c <- stan(file = "R/FitTopt.stan",
                data = list(Vcmax=params.all %>% filter(Treat=="c.c") %>% pull(Vcmax),
                            Tk=params.all %>% filter(Treat=="c.c") %>% pull(Tk),
                            N = params.all %>% filter(Treat=="c.c") %>% dim() %>% .[1]),
                warmup = 2000,
                iter=4000, thin=1, chains=3, verbose=T,
                cores=3,
                control=list(adapt_delta=0.9))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(stan_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(stan_c.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(stan_c.c), c("sigma","Topt","Ha","kopt"))

stan_c.w <- stan(file = "R/FitTopt.stan",
                 data = list(Vcmax=params.all %>% filter(Treat=="c.w") %>% pull(Vcmax),
                             Tk=params.all %>% filter(Treat=="c.w") %>% pull(Tk),
                             N = params.all %>% filter(Treat=="c.w") %>% dim() %>% .[1]),
                 warmup = 2000,
                 iter=4000, thin=1, chains=3, verbose=T,
                 cores=3,save_dso = T,
                 control=list(adapt_delta=0.9))
print(stan_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(stan_c.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(stan_c.w), c("sigma","Topt","Ha","kopt"))

stan_w.c <- stan(file = "R/FitTopt.stan",
                 data = list(Vcmax=params.all %>% filter(Treat=="w.c") %>% pull(Vcmax),
                             Tk=params.all %>% filter(Treat=="w.c") %>% pull(Tk),
                             N = params.all %>% filter(Treat=="w.c") %>% dim() %>% .[1]),
                 warmup = 2000,
                 iter=4000, thin=1, chains=3, verbose=T,
                 cores=3,save_dso = T,
                 control=list(adapt_delta=0.9))
print(stan_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(stan_w.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(stan_w.c), c("sigma","Topt","Ha","kopt"))

stan_w.w <- stan(file = "R/FitTopt.stan",
                 data = list(Vcmax=params.all %>% filter(Treat=="w.w") %>% pull(Vcmax),
                             Tk=params.all %>% filter(Treat=="w.w") %>% pull(Tk),
                             N = params.all %>% filter(Treat=="w.w") %>% dim() %>% .[1]),
                 warmup = 2000,
                 iter=4000, thin=1, chains=3, verbose=T,
                 cores=3,save_dso = T,
                 control=list(adapt_delta=0.9))
print(stan_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(stan_w.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(stan_w.w), c("sigma","Topt","Ha","kopt"))


# Plot Predicted Fits Across Treatments -----------------------------------
r_Vcmax <- function(kopt, Hd=200, Ha, Tk, Topt){
  Vcmax <- kopt * ( (Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                      (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))    ); 
  return(Vcmax); 
}

post_df <- stan_c.c %>% as.data.frame() %>% as_tibble() %>% mutate(Treat="c.c")
post_df <- bind_rows(post_df, stan_c.w %>% as.data.frame() %>% as_tibble() %>% mutate(Treat="c.w"))
post_df <- bind_rows(post_df, stan_w.c %>% as.data.frame() %>% as_tibble() %>% mutate(Treat="w.c"))
post_df <- bind_rows(post_df, stan_w.w %>% as.data.frame() %>% as_tibble() %>% mutate(Treat="w.w"))
post_df <- post_df %>% 
  sample_n(100) %>% 
  left_join(., 
            as_tibble(expand.grid(Tk=seq(298,318,by=0.25), Treat=c("c.c","c.w","w.c","w.w")))) %>% 
  mutate(Vcmax=r_Vcmax(kopt=kopt,Hd=200,Ha=Ha,Tk=Tk,Topt=Topt))

topt_df <- tibble(Tk = summary(stan_c.c)$summary["Topt","50%"], Treat="c.c") %>% 
          bind_rows(.,
          tibble(Tk = summary(stan_c.w)$summary["Topt","50%"], Treat="c.w"), 
          tibble(Tk = summary(stan_w.c)$summary["Topt","50%"], Treat="w.c"), 
          tibble(Tk = summary(stan_w.w)$summary["Topt","50%"], Treat="w.w"))

post_df %>% 
  group_by(Treat, Tk) %>% 
  summarize(Vcmax_0.5 = median(Vcmax), 
            Vcmax_0.95 = quantile(Vcmax, 0.95), 
            Vcmax_0.05 = quantile(Vcmax, 0.05)) %>% 
  ggplot(data=., aes(Tk, Vcmax_0.5, color=Treat))+
  geom_ribbon(aes(x=Tk, ymax=Vcmax_0.95, ymin=Vcmax_0.05, fill=Treat), 
              alpha=0.25, lty=0)+
  geom_line()+
  geom_point(data=params.all, aes(Tk, Vcmax,color=Treat))+
  geom_vline(data=topt_df, aes(xintercept=Tk, color=Treat))+
  labs(x='Temperature [K]',y="Vcmax [ ]")+
  theme_bw()+
  scale_color_viridis_d(end=0.8)+
  scale_fill_viridis_d(end=0.8)



df_med <- lapply(post_df,"median") %>% as.data.frame()
df_5 <- lapply(post_df,"quantile",0.05) %>% as.data.frame();
df_95 <- lapply(post_df,"quantile",0.95) %>% as.data.frame();
df_ppd <- tibble(Tk=290:320,Vcmax05=NA,Vcmax50=NA,Vcmax95=NA);

for(i in 1:dim(df_ppd)[1]){
  df_ppd$Vcmax05[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.05, na.rm=T)
  df_ppd$Vcmax50[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.5, na.rm=T)
  df_ppd$Vcmax95[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.95, na.rm=T)
}



png("figures/CompareToptEst_NLS_BRMS.png", height = 800, width=1000, units='px', pointsize = 18)
sim_Tk <- seq(298,327,0.1)
params.all %>% 
  plot(Vcmax~Tk, data=., xlim=c(299,322), main="Topt");

lines(posterior_predict(bfit_c.c_Hd200, transform = colMeans, 
                        newdata = data.frame(Tk=sim_Tk))~sim_Tk, 
      col='black')
lines(posterior_predict(bfit_c.w_Hd200, transform = colMeans, 
                        newdata = data.frame(Tk=sim_Tk))~sim_Tk, 
      col='orange')
lines(posterior_predict(bfit_w.c_Hd200, transform = colMeans, 
                        newdata = data.frame(Tk=sim_Tk))~sim_Tk, 
      col='orchid')
lines(posterior_predict(bfit_w.w_Hd200, transform = colMeans, 
                        newdata = data.frame(Tk=sim_Tk))~sim_Tk, 
      col='red')

summary(stan_c.c)$summary["Topt","50%"]

abline(v=summary(stan_c.c)$summary["Topt","50%"],col='black',lwd=2)
abline(v=summary(stan_c.w)$summary["Topt","50%"],col='orange',lwd=2)
abline(v=summary(stan_w.c)$summary["Topt","50%"],col='orchid',lwd=2)
abline(v=summary(stan_w.w)$summary["Topt","50%"],col='red',lwd=2)


abline(v=summary(nfit_c.c)$coefficients["Topt","Estimate"],col='black',lty=3,lwd=2)
abline(v=summary(nfit_c.w)$coefficients["Topt","Estimate"],col='orange',lty=3,lwd=2)
abline(v=summary(nfit_w.c)$coefficients["Topt","Estimate"],col='orchid',lty=3,lwd=2)
abline(v=summary(nfit_w.w)$coefficients["Topt","Estimate"],col='red',lty=3,lwd=2)

legend("topleft", col=c("black","orange" ,"orchid","red"), lty=c(1,1,1,1), 
       legend=c("c.c","c.w","w.c","w.w"))
legend('bottomright', legend=c('Bayesian Fit', 'NLS Fit'), lty=c(1,3))
dev.off()
