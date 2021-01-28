library(tidyverse); library(bayesplot);
library(rstan)
# library(xfun); library(brms); 
options(mc.cores = parallel::detectCores()-1)
set.seed(1); 

dat400 <- read_csv("data/P400.2904.csv")
dat800 <- read_csv("data/P800.2904.csv")

#*****************************************************************************
# Photo 400 FIT TREATMENTS c.c & c.w w.c w.w                 ----------------
#*****************************************************************************
library(rstan)
set.seed(3)
photo400_c.c <- stan(file = "Stan/FitPhoto.stan",
                   data = list(A=dat400 %>% filter(Treat=="c.c") %>% pull(Photo),
                               Tk=dat400 %>% filter(Treat=="c.c") %>% pull(Tk),
                               N = dat400 %>% filter(Treat=="c.c") %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=4000, thin=2, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo400_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo400_c.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo400_c.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

shinystan::launch_shinystan(photo400_c.c)

photo400_c.w <- stan(file = "Stan/FitPhoto.stan",
                   data=list(A=dat400 %>% filter(Treat=="c.w") %>% pull(Photo), 
                             Tk=dat400 %>% filter(Treat=="c.w") %>% pull(Tk),
                             N=dat400 %>% filter(Treat=="c.w") %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=4000, thin=2, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.99))
print(photo400_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo400_c.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo400_c.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()
#*******************************************************************************

#*****************************************************************************
# FIT TREATMENTS w.w & w.c ----------------
#*****************************************************************************
library(rstan)
set.seed(3)
photo400_w.w <- stan(file = "Stan/FitPhoto.stan",
                   data = list(A=dat400 %>% filter(Treat=="w.w") %>% pull(Photo),
                               Tk=dat400 %>% filter(Treat=="w.w") %>% pull(Tk),
                               N = dat400 %>% filter(Treat=="w.w") %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=4000, thin=2, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo400_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo400_w.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo400_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

photo400_w.c <- stan(file = "Stan/FitPhoto.stan",
                   data=list(A=dat400 %>% filter(Treat=="w.c") %>% pull(Photo), 
                             Tk=dat400 %>% filter(Treat=="w.c") %>% pull(Tk),
                             N=dat400 %>% filter(Treat=="w.c") %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=6000, thin=2, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.95))
print(photo400_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo400_w.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo400_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()
#*******************************************************************************

#*****************************************************************************
# Photo 800 FIT TREATMENTS c.c & c.w w.c w.w                 ----------------
#*****************************************************************************
library(rstan)
set.seed(3)
photo800_c.c <- stan(file = "Stan/FitPhoto.stan",
                     data = list(A=dat800 %>% filter(Treat=="c.c") %>% pull(Photo),
                                 Tk=dat800 %>% filter(Treat=="c.c") %>% pull(Tk),
                                 N = dat800 %>% filter(Treat=="c.c") %>% dim() %>% .[1]),
                     warmup = 2000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo800_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo800_c.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo800_c.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


photo800_c.w <- stan(file = "Stan/FitPhoto.stan",
                     data=list(A=dat800 %>% filter(Treat=="c.w") %>% pull(Photo), 
                               Tk=dat800 %>% filter(Treat=="c.w") %>% pull(Tk),
                               N=dat800 %>% filter(Treat=="c.w") %>% dim() %>% .[1]),
                     warmup = 2000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.99))
print(photo800_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo800_c.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo800_c.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()
#*******************************************************************************

#*****************************************************************************
# FIT TREATMENTS w.w & w.c ----------------
#*****************************************************************************
library(rstan)
set.seed(3)
photo800_w.w <- stan(file = "Stan/FitPhoto.stan",
                     data = list(A=dat800 %>% filter(Treat=="w.w") %>% pull(Photo),
                                 Tk=dat800 %>% filter(Treat=="w.w") %>% pull(Tk),
                                 N = dat800 %>% filter(Treat=="w.w") %>% dim() %>% .[1]),
                     warmup = 2000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo800_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo800_w.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo800_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

photo800_w.c <- stan(file = "Stan/FitPhoto.stan",
                     data=list(A=dat800 %>% filter(Treat=="w.c") %>% pull(Photo), 
                               Tk=dat800 %>% filter(Treat=="w.c") %>% pull(Tk),
                               N=dat800 %>% filter(Treat=="w.c") %>% dim() %>% .[1]),
                     warmup = 2000, save_dso = T,
                     iter=6000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.95))
print(photo800_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(photo800_w.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(photo800_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()
#*******************************************************************************

res_p400_p800 <- bind_rows(
  summary(photo400_c.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","lp__"), 
           Treat='c.c', CO2 = 400), 
  summary(photo400_c.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
           Treat='c.w', CO2 = 400), 
  summary(photo400_w.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
           Treat='w.c', CO2 = 400), 
  summary(photo400_w.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
           Treat='w.w', CO2 = 400),
  summary(photo800_c.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","lp__"), 
           Treat='c.c', CO2 = 800), 
  summary(photo800_c.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
           Treat='c.w', CO2 = 800), 
  summary(photo800_w.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
           Treat='w.c', CO2 = 800), 
  summary(photo800_w.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","lp__"),
           Treat='w.w', CO2 = 800))
  

  
res_p400_p800 %>% 
  filter(!param %in% c('lp__','sigma')) %>%
  ggplot(data=., aes(x=Treat, y=`50%`-273.15))+
  # geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`),col='darkolivegreen3',lwd=2)+
  # geom_linerange(aes(ymin = `25%`, ymax = `75%`),col='darkolivegreen4',lwd=3)+
  geom_point(col='black')+
  theme_bw()+
  facet_wrap(~param+CO2, 
             # labeller = 'label_both',
             scales = 'free',
             ncol = 2, as.table = T)

write_csv(res_p400_p800, path = paste0("outputs/params_photosyn_Ca_400_800_",Sys.Date(),".csv"))




