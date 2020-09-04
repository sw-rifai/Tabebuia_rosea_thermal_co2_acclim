### Plot Vcmax CO2 400 & 800
library(tidyverse); library(bayesplot); library(rstan)
options(mc.cores = parallel::detectCores()-1); 
rstan_options(auto_write = TRUE);
set.seed(1); 

dat_all<-read_csv("data/params.all.28032018.csv")
dat_all <- dat_all %>% 
  mutate(Tk = Tleaf+273.15)
dat400 <- read_csv("data/P400.2904.csv")
dat800 <- read_csv("data/P800.2904.csv")

# treatments: (w.w, w.w, w.c, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
# treat<-as.character(unique(params.all$Treat))

# Fit Vcmax ---------------------------------------------------------------------
vc400_c.c <- stan(file = "Stan/FitTopt_sim_Vcmax.stan",
                   data = list(Vcmax=dat_all %>% filter(Treat=='c.c') %>% pull(Vcmax),
                               Tk=dat_all %>% filter(Treat=='c.c') %>% pull(Tk),
                               N = dat_all %>% filter(Treat=='c.c') %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=4000, thin=2, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(vc400_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma","Vcmax_Topt"))
mcmc_trace(as.array(vc400_c.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(vc400_c.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


vc400_w.w <- stan(file = "Stan/FitTopt_sim_Vcmax.stan",
                  data = list(Vcmax=dat_all %>% filter(Treat=='w.w') %>% pull(Vcmax),
                              Tk=dat_all %>% filter(Treat=='w.w') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='w.w') %>% dim() %>% .[1]),
                  warmup = 2000, save_dso = T,
                  iter=4000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(vc400_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(vc400_w.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(vc400_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

vc400_c.w <- stan(file = "Stan/FitTopt_sim_Vcmax.stan",
                  data = list(Vcmax=dat_all %>% filter(Treat=='c.w') %>% pull(Vcmax),
                              Tk=dat_all %>% filter(Treat=='c.w') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='c.w') %>% dim() %>% .[1]),
                  warmup = 2000, save_dso = T,
                  iter=4000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(vc400_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(vc400_c.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(vc400_c.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

vc400_w.c <- stan(file = "Stan/FitTopt_sim_Vcmax.stan",
                  data = list(Vcmax=dat_all %>% filter(Treat=='w.c') %>% pull(Vcmax),
                              Tk=dat_all %>% filter(Treat=='w.c') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='w.c') %>% dim() %>% .[1]),
                  warmup = 2000, save_dso = T,
                  iter=4000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(vc400_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(vc400_w.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(vc400_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


# Fit Jmax ---------------------------------------------------------------------
jm400_c.c <- stan(file = "Stan/FitTopt_Jmax.stan",
                  data = list(Jmax=dat_all %>% filter(Treat=='c.c') %>% pull(Jmax),
                              Tk=dat_all %>% filter(Treat=='c.c') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='c.c') %>% dim() %>% .[1]),
                  warmup = 2000, save_dso = T,
                  iter=4000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(jm400_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma","Jmax_Topt"))
mcmc_trace(as.array(jm400_c.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(jm400_c.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


jm400_w.w <- stan(file = "Stan/FitTopt_Jmax.stan",
                  data = list(Jmax=dat_all %>% filter(Treat=='w.w') %>% pull(Jmax),
                              Tk=dat_all %>% filter(Treat=='w.w') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='w.w') %>% dim() %>% .[1]),
                  warmup = 2000, save_dso = T,
                  iter=4000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(jm400_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(jm400_w.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(jm400_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

jm400_c.w <- stan(file = "Stan/FitTopt_Jmax.stan",
                  data = list(Jmax=dat_all %>% filter(Treat=='c.w') %>% pull(Jmax),
                              Tk=dat_all %>% filter(Treat=='c.w') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='c.w') %>% dim() %>% .[1]),
                  warmup = 2000, save_dso = T,
                  iter=4000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(jm400_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(jm400_c.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(jm400_c.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

jm400_w.c <- stan(file = "Stan/FitTopt_Jmax.stan",
                  data = list(Jmax=dat_all %>% filter(Treat=='w.c') %>% pull(Jmax),
                              Tk=dat_all %>% filter(Treat=='w.c') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='w.c') %>% dim() %>% .[1]),
                  warmup = 2000, save_dso = T,
                  iter=4000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(jm400_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(jm400_w.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(jm400_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()





res_vc_j <- bind_rows(
  summary(vc400_c.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","Vcmax_Topt","lp__"), 
           Treat='c.c', rate = 'Vcmax'), 
  summary(vc400_c.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","Vcmax_Topt","lp__"),
           Treat='c.w', rate = 'Vcmax'), 
  summary(vc400_w.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","Vcmax_Topt","lp__"),
           Treat='w.c', rate = 'Vcmax'), 
  summary(vc400_w.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","Vcmax_Topt","lp__"),
           Treat='w.w', rate = 'Vcmax'),
  summary(jm400_c.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","Jmax_Topt","lp__"),
           Treat='c.c', rate = 'Jmax'), 
  summary(jm400_c.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","Jmax_Topt","lp__"),
           Treat='c.w', rate = 'Jmax'), 
  summary(jm400_w.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","Jmax_Topt","lp__"),
           Treat='w.c', rate = 'Jmax'), 
  summary(jm400_w.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("kopt","Ha","Topt","sigma","Jmax_Topt","lp__"),
           Treat='w.w', rate = 'Jmax'))

res_vc_j %>% 
  filter(!param %in% c('lp__','sigma')) %>%
  ggplot(data=., aes(x=Treat, y=`50%`))+
  geom_linerange(aes(ymin = `10%`, ymax = `90%`),col='darkolivegreen3',lwd=2)+
  geom_linerange(aes(ymin = `25%`, ymax = `75%`),col='darkolivegreen4',lwd=3)+
  geom_point(col='black')+
  theme_bw()+
  facet_wrap(~param+rate, 
             labeller = 'label_both',
             scales = 'free', ncol = 2)

write_csv(res_vc_j, path = paste0("outputs/params_vcmax_jmax_",Sys.Date(),".csv"))






