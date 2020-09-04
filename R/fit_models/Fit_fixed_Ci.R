# Hi Sami, 
# Sorry to give you more work to do, but I wonder if you'd be able to fit some more
# curves for the Tabebuia data. From each A-ci curve I calculated what photosynthesis
# would be at a fixed Ci, roughly corresponding to what Ci is at current ambient CO2. 
# By fitting curves of photosynthesis controlled for Ci the effect of stomatal 
# conductance is taken out of the equation, so the comparison of the parameters of 
# these curves with those of the P400, should give an indication of the role of 
# stomatal conductance on driving the temperature response. It looks like the curves
# will be practically the same as those for P400 (so probably don't need different
# priors etc), suggesting that in contrast to the tree work I've done before, these
# saplings are not strongly controlled by conductance at all. 
# 
# There's a column with P270, and one for P505, which is the fixed Ci corresponding
# to Ci under 800 ppm, "ambient" in the treatment dome.
# 
# Cheers,
# Martijn

### Plot Vcmax CO2 400 & 800

library(tidyverse); library(bayesplot); library(rstan)
options(mc.cores = parallel::detectCores()-1); 
rstan_options(auto_write = TRUE);
set.seed(1); 

dat_400<-read_csv("data/P400.2904.csv")
dat_all<-read_csv("data/params.all.28032018.csv")
dat_fci<-read_csv("data/params_20190124.csv")
dat_fci <- dat_fci %>% 
  mutate(Tk = Tleaf+273.15)
dat400 <- read_csv("data/P400.2904.csv")
dat800 <- read_csv("data/P800.2904.csv")

# treatments: (w.w, w.w, w.c, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
# treat<-as.character(unique(params.all$Treat))

# Fit photo T opt ---------------------------------------------------------------------
# CC - P270
photo270_cc_fci <- stan(file = "Stan/FitPhoto_fixed_Ci.stan",
                  data = list(A=dat_fci %>% filter(Treat=='c.c') %>% pull(P270),
                              Tk=dat_fci %>% filter(Treat=='c.c') %>% pull(Tk),
                              N = dat_fci %>% filter(Treat=='c.c') %>% dim() %>% .[1]),
                  warmup = 4000, save_dso = T,
                  iter=6000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo270_cc_fci, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo270_cc_fci), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo270_cc_fci), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

# CC - P505
photo505_cc_fci <- stan(file = "Stan/FitPhoto_fixed_Ci.stan",
                        data = list(A=dat_fci %>% filter(Treat=='c.c') %>% pull(P505),
                                    Tk=dat_fci %>% filter(Treat=='c.c') %>% pull(Tk),
                                    N = dat_fci %>% filter(Treat=='c.c') %>% dim() %>% .[1]),
                        warmup = 2000, save_dso = T,
                        iter=4000, thin=2, chains=3, verbose=T,
                        cores=3,
                        control=list(adapt_delta=0.95))
print(photo505_cc_fci, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo505_cc_fci), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo505_cc_fci), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

#------------------------------------------------------------------------------------------
# cw - P270
photo270_cw_fci <- stan(file = "Stan/FitPhoto_fixed_Ci.stan",
                        data = list(A=dat_fci %>% filter(Treat=='c.w') %>% pull(P270),
                                    Tk=dat_fci %>% filter(Treat=='c.w') %>% pull(Tk),
                                    N = dat_fci %>% filter(Treat=='c.w') %>% dim() %>% .[1]),
                        warmup = 4000, save_dso = T,
                        iter=6000, thin=2, chains=3, verbose=T,
                        cores=3,
                        control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo270_cw_fci, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo270_cw_fci), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo270_cw_fci), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

# cw - P505
photo505_cw_fci <- stan(file = "Stan/FitPhoto_fixed_Ci.stan",
                        data = list(A=dat_fci %>% filter(Treat=='c.w') %>% pull(P505),
                                    Tk=dat_fci %>% filter(Treat=='c.w') %>% pull(Tk),
                                    N = dat_fci %>% filter(Treat=='c.w') %>% dim() %>% .[1]),
                        warmup = 2000, save_dso = T,
                        iter=4000, thin=2, chains=3, verbose=T,
                        cores=3,
                        control=list(adapt_delta=0.95))
print(photo505_cw_fci, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo505_cw_fci), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo505_cw_fci), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


#---------------------------------------------------------------------------------
# wc - P270
photo270_wc_fci <- stan(file = "Stan/FitPhoto_fixed_Ci.stan",
                        data = list(A=dat_fci %>% filter(Treat=='w.c') %>% pull(P270),
                                    Tk=dat_fci %>% filter(Treat=='w.c') %>% pull(Tk),
                                    N = dat_fci %>% filter(Treat=='w.c') %>% dim() %>% .[1]),
                        warmup = 4000, save_dso = T,
                        iter=6000, thin=2, chains=3, verbose=T,
                        cores=3,
                        control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo270_wc_fci, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo270_wc_fci), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo270_wc_fci), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

# wc - P505
photo505_wc_fci <- stan(file = "Stan/FitPhoto_fixed_Ci.stan",
                        data = list(A=dat_fci %>% filter(Treat=='w.c') %>% pull(P505),
                                    Tk=dat_fci %>% filter(Treat=='w.c') %>% pull(Tk),
                                    N = dat_fci %>% filter(Treat=='w.c') %>% dim() %>% .[1]),
                        warmup = 2000, save_dso = T,
                        iter=4000, thin=2, chains=3, verbose=T,
                        cores=3,
                        control=list(adapt_delta=0.95))
print(photo505_wc_fci, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo505_wc_fci), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo505_wc_fci), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

#---------------------------------------------------------------------------------
# ww - P270
photo270_ww_fci <- stan(file = "Stan/FitPhoto_fixed_Ci.stan",
                        data = list(A=dat_fci %>% filter(Treat=='w.w') %>% pull(P270),
                                    Tk=dat_fci %>% filter(Treat=='w.w') %>% pull(Tk),
                                    N = dat_fci %>% filter(Treat=='w.w') %>% dim() %>% .[1]),
                        warmup = 4000, save_dso = T,
                        iter=6000, thin=2, chains=3, verbose=T,
                        cores=3,
                        control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo270_ww_fci, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo270_ww_fci), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo270_ww_fci), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

# ww - P505
photo505_ww_fci <- stan(file = "Stan/FitPhoto_fixed_Ci.stan",
                        data = list(A=dat_fci %>% filter(Treat=='w.w') %>% pull(P505),
                                    Tk=dat_fci %>% filter(Treat=='w.w') %>% pull(Tk),
                                    N = dat_fci %>% filter(Treat=='w.w') %>% dim() %>% .[1]),
                        warmup = 2000, save_dso = T,
                        iter=4000, thin=2, chains=3, verbose=T,
                        cores=3,
                        control=list(adapt_delta=0.95))
print(photo505_ww_fci, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo505_ww_fci), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo505_ww_fci), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

as.data.frame(photo270_cc_fci)


res_fci <- bind_rows(
summary(photo270_cc_fci, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
  as_tibble() %>% 
  mutate(param = c("kopt","Ha","Topt","sigma","lp__"), 
         Treat='c.c', Ci = 'fixed 207'), 
summary(photo270_cw_fci, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
  as_tibble() %>% 
  mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
         Treat='c.w', Ci = 'fixed 207'), 
summary(photo270_wc_fci, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
  as_tibble() %>% 
  mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
         Treat='w.c', Ci = 'fixed 207'), 
summary(photo270_ww_fci, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
  as_tibble() %>% 
  mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
         Treat='w.w', Ci = 'fixed 207'),
summary(photo505_cc_fci, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
  as_tibble() %>% 
  mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
         Treat='c.c', Ci = 'fixed 505'), 
summary(photo505_cw_fci, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
  as_tibble() %>% 
  mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
         Treat='c.w', Ci = 'fixed 505'), 
summary(photo505_wc_fci, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
  as_tibble() %>% 
  mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
         Treat='w.c', Ci = 'fixed 505'), 
summary(photo505_ww_fci, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
  as_tibble() %>% 
  mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
         Treat='w.w', Ci = 'fixed 505'))

res_fci %>% 
  filter(!param %in% c('lp__','sigma')) %>%
  ggplot(data=., aes(x=Treat, y=`50%`))+
  geom_linerange(aes(ymin = `10%`, ymax = `90%`),col='darkolivegreen3',lwd=2)+
  geom_linerange(aes(ymin = `25%`, ymax = `75%`),col='darkolivegreen4',lwd=3)+
  geom_point(col='black')+
  theme_bw()+
  facet_wrap(~param+Ci, scales = 'free', ncol = 2)

res_fci %>% write_csv(x = ., path = paste0("outputs/params_photosyn_fixedCI_",
                                           Sys.Date(),".csv"))



# Plot fixed Ci -----------------------------------------------------------
f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}

object <- photo270_cc_fci
tmp1 <- dat_fci %>% filter(Treat=='c.c')
p1 <- object %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,319,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='blue')+
  theme_bw()+
  labs(y=bquote(paste('V'['cmax'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  # ylim(0,325)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p1 <- p1 + annotate("text", 
                    x = 21, 
                    y = c(summary(object)$summary["kopt","50%"]*1.25,
                          summary(object)$summary["kopt","50%"]*1.25 - 5), 
                    hjust=0,size=4,
                    label = c("Control", 
                              "Treatment"),parse=T)+
  annotate("point", x = 26.5, y = c(summary(object)$summary["kopt","50%"]*1.25,
                                    summary(object)$summary["kopt","50%"]*1.25 - 5), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,16))
tmp1 <- (tmp1 %>% filter(Treat=="c.c")) %>% select(Tk,P270)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,P270), size=4, alpha=0.5,
                      col='blue', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(object)$summary["Topt","50%"]-273.15)), 
                      col='blue')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(object)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(object)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(object, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(object, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)
p1
