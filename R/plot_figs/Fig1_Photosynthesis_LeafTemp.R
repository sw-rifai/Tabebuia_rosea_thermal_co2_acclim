library(nls.multstart);
library(tidyverse); library(bayesplot);
library(rstan)
# library(xfun); library(brms); 
options(mc.cores = parallel::detectCores()-1)
set.seed(1); 

params.all<-read_csv("data/params.all.28032018.csv")
dat400 <- read_csv("data/P400.2904.csv")
dat800 <- read_csv("data/P800.2904.csv")

params.all$Tk<-params.all$Tleaf+273.15

# create empty frame for Medlyn parameters
pars.med<-NULL

# treatments: (c.c, w.w, c.w, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
treat<-as.character(unique(params.all$Treat))

dat400 %>% glimpse
dat400 %>% ggplot(data=., aes(Tleaf, Photo))+geom_point()
params.all %>% filter(Treat=="c.c") %>% glimpse
params.all %>% filter(Treat=="c.c") %>% pull(Plant) %>% table
params.all %>% filter(Treat=="c.c") %>% pull(Leaf) %>% table
dat400 %>% 
  ggplot(data=., aes(Tk, Photo, color=Treat))+geom_point()+
  geom_smooth()


#*****************************************************************************
# FIT TREATMENT c.c ----------------
#*****************************************************************************
library(rstan)
sphoto_c.c <- stan(file = "R/FitPhoto.stan",
                 data = list(A=dat400 %>% filter(Treat=="c.c") %>% pull(Photo),
                             Tk=dat400 %>% filter(Treat=="c.c") %>% pull(Tk),
                             N = dat400 %>% filter(Treat=="c.c") %>% dim() %>% .[1]),
                 warmup = 2000, save_dso = T,
                 iter=4000, thin=1, chains=3, verbose=T,
                 cores=3,
                 control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(sphoto_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(sphoto_c.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(sphoto_c.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

sphoto_c.w <- stan(file = "R/FitPhoto.stan",
                   data=list(A=dat400 %>% filter(Treat=="c.w") %>% pull(Photo), 
                        Tk=dat400 %>% filter(Treat=="c.w") %>% pull(Tk),
                        N=dat400 %>% filter(Treat=="c.w") %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=4000, thin=1, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.95))
print(sphoto_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(sphoto_c.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(sphoto_c.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

sphoto_w.c <- stan(file = "R/FitPhoto.stan",
                   data=list(A=dat400 %>% filter(Treat=="w.c") %>% pull(Photo), 
                             Tk=dat400 %>% filter(Treat=="w.c") %>% pull(Tk),
                             N=dat400 %>% filter(Treat=="w.c") %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=4000, thin=1, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.95))
print(sphoto_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(sphoto_w.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(sphoto_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

sphoto_w.w <- stan(file = "R/FitPhoto.stan",
                   data=list(A=dat400 %>% filter(Treat=="w.w") %>% pull(Photo), 
                             Tk=dat400 %>% filter(Treat=="w.w") %>% pull(Tk),
                             N=dat400 %>% filter(Treat=="w.w") %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=4000, thin=1, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.95))
print(sphoto_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(sphoto_w.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(sphoto_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

library(rethinking)
rethinking::HPDI(as.matrix(sphoto_c.c, pars="Topt"), prob = 0.9)
rethinking::HPDI(as.matrix(sphoto_c.w, pars="Topt"), prob = 0.9)
rethinking::HPDI(as.matrix(sphoto_w.c, pars="Topt"), prob = 0.9)
rethinking::HPDI(as.matrix(sphoto_w.w, pars="Topt"), prob = 0.9)

f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
          (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}

print(sphoto_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))
sphoto_c.c %>% print

# Plot Stan Fit c.c and c.w -----------------------------------------------
p1 <- sphoto_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='blue')+
  theme_bw()+
  labs(y=bquote(paste('P'[400],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,40)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p1 <- p1 + annotate("text", x = 36, y = c(37,33), hjust=0,size=4,
            label = c("Control", 
                      "paste(Control %->% \" Treatment\")"),parse=T)+
  annotate("point", x = 35, y = c(37,33), 
           colour = c('blue',"red"), 
           size = 5, shape=c(16,1))
tmp1 <- (dat400 %>% filter(Treat=="c.c")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=4, alpha=0.5,
                col='blue', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(sphoto_c.c)$summary["Topt","50%"]-273.15)), 
                      col='blue')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(sphoto_c.c)$summary["Topt","2.5%"]-273.15), 
                      col='blue',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(sphoto_c.c)$summary["Topt","97.5%"]-273.15), 
                          col='blue',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(sphoto_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
                          col='blue',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(sphoto_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)


p2_lines <- sphoto_c.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='red')
tmp2 <- (dat400 %>% filter(Treat=="c.w")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, shape=1,
                      col='red', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(sphoto_c.w)$summary["Topt","50%"]-273.15), 
                      col='red')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(sphoto_c.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(sphoto_c.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(sphoto_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(sphoto_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)

p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave("figures/p400_c.c_c.w.png", width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


# Plot Stan Fit w.w and w.c -----------------------------------------------
p1 <- sphoto_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  # arrange(desc(lp)) %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  # scale_color_viridis_c("B", direction=1, end=0.975)+
  theme_bw()+
  labs(y=bquote(paste('P'[400],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,40)
p1
tmp1 <- (dat400 %>% filter(Treat=="w.w")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=4, alpha=0.6,
                      col='red', show.legend = F, inherit.aes = F)
p1  
p1_topt <- geom_vline(aes(xintercept=(summary(sphoto_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(sphoto_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(sphoto_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)

p2_lines <- sphoto_w.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  # arrange(desc(lp)) %>% 
  # ggplot(data=., aes(Tk, Photo, group=lp, color=lp))+
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
# scale_color_viridis_c("B", direction=1, end=0.975)+
# theme_bw()
tmp2 <- (dat400 %>% filter(Treat=="w.c")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, shape=1,
                        col='blue', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(sphoto_c.w)$summary["Topt","50%"]-273.15), 
                      col='blue')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(sphoto_c.w)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(sphoto_c.w)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_2.5+p1_topt_97.5+
  p2_topt_2.5+p2_topt_97.5
ggsave("figures/p400_w.w_w.c.png", width = 120, height=100, units='mm')

p1+p2_lines+p2_points+p1_topt+p2_topt
# p1_topt_2.5+p1_topt_97.5+
# p2_topt_2.5+p2_topt_97.5
ggsave("figures/p400_w.w_w.c_noToptCI.png", width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


sphoto_c.c %>% 
  as.data.frame() %>% 
  head() %>%
  as_tibble() %>% 
  complete(., Tk=full_seq(285,310,1))

  # ggplot(data=., aes(x=Tk))+
  # stat_function(fun=f_photo, args=list(kopt=kopt,Ha=Ha,Tk=Tk))

library(tidybayes)
sphoto_c.c %>% recover_types()
sphoto_c.c %>% spread_draws("Topt") %>% head()











nls_c.c <- nls_multstart(Photo ~ kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                         data = dat400 %>% filter(Treat=="c.w"),
                         iter = 1000,
                         start_lower = c(kopt = 160, Ha = 50, Topt = 300),
                         start_upper = c(kopt = 300, Ha = 150, Topt = 330),
                         supp_errors = 'Y',
                         na.action = na.omit,
                         #convergence_count = 500,
                         lower = c(kopt = 10, Ha = 10, Topt = 283))
summary(nls_c.c)
tibble(Photo=predict(nls_c.c, newdata=tibble(Tk=seq(275,330))),
  Tk=seq(275,330), 
  Tc=Tk-273.15) %>% 
  ggplot(data=., aes(Tc,Photo))+geom_point()+
  geom_point(data=dat400 %>% filter(Treat=="c.w"), aes(Tleaf, Photo),col='blue')









# BRMS approach -----------------------------------------------------------
# parameters {
#   real<lower=10,upper=500> kopt;
#   // real<lower=150,upper=250> Hd;
#   real<lower=10,upper=250> Ha; 
#   real<lower=283,upper=335> Topt;
#   real<lower=0> sigma;
# }
# kopt ~ normal(25,100);
# // Hd ~ normal(200,10);
# Ha ~ normal(100,25); 
# Topt ~ normal(310,10); 

library(brms)
bphoto_c.c <- brm(
  bf(Photo ~ kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                       (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))),
     kopt~1,
     Ha~1,
     # Hd~1, # Hd now fixed at 200
     Topt~1,
     nl=T),
  data = dat400 %>% filter(Treat=="c.c"),
  family = Gamma(link='identity'),
  chains = 3,
  iter = 20000, thin=20,
  prior = c(prior(normal(50,5),nlpar="Ha",lb=10,ub=250), #25,150
            prior(normal(50, 5), nlpar="kopt",lb=10,ub=500), # 230,300
            # prior(normal(200,1),nlpar="Hd",lb=190,ub=210), # 150,250
            prior(normal(300,10),nlpar="Topt",lb=283,ub=325)),
  control=list(adapt_delta = 0.995, max_treedepth = 15))
plot(bphoto_c.c, type='mcmc_trace')
sphoto_c.c
bphoto_c.c %>% summary()
bphoto_c.c %>% marginal_effects(points=T)
bphoto_c.c %>% pairs()

bphoto_c.w <- brm(
  bf(Photo ~ kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                       (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))),
     kopt~1,
     Ha~1,
     # Hd~1, # Hd now fixed at 200
     Topt~1,
     nl=T),
  data = dat400 %>% filter(Treat=="c.w"),
  family = gaussian(link='identity'),
  chains = 3,
  iter = 20000, thin=20,
  prior = c(prior(normal(50,10),nlpar="Ha",lb=10,ub=250), #25,150
            prior(normal(15, 5), nlpar="kopt",lb=7,ub=30), # 230,300
            # prior(normal(200,1),nlpar="Hd",lb=190,ub=210), # 150,250
            prior(normal(300,15),nlpar="Topt",lb=283,ub=325)),
  control=list(adapt_delta = 0.995, max_treedepth = 15))
bphoto_c.w
bphoto_c.w %>% marginal_effects()
bphoto_c.w %>% stanplot(., type='dens')


bphoto_w.c <- brm(
  bf(Photo ~ kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                       (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))),
     kopt~1,
     Ha~1,
     # Hd~1, # Hd now fixed at 200
     Topt~1,
     nl=T),
  data = dat400 %>% filter(Treat=="w.c"),
  family = gaussian(link='identity'),
  chains = 3,
  iter = 20000, thin=20,
  prior = c(prior(normal(50,10),nlpar="Ha",lb=10,ub=250), #25,150
            prior(normal(15, 5), nlpar="kopt",lb=7,ub=30), # 230,300
            # prior(normal(200,1),nlpar="Hd",lb=190,ub=210), # 150,250
            prior(normal(300,15),nlpar="Topt",lb=283,ub=325)),
  control=list(adapt_delta = 0.995, max_treedepth = 15))
bphoto_w.c
bphoto_w.c %>% marginal_effects()
bphoto_w.c %>% stanplot(., type='dens')



brms::stanplot(bphoto_c.w, type='trace')

f_photo <- function(kopt,Ha,Tk,Topt,Hd=200){
  photo <- kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}
curve(f_photo(kopt = 20, Ha=100, Topt=305, Tk=x),290,320, ylim=c(0,30))
curve(f_photo(kopt = 20, Ha=100, Topt=305, Tk=x, Hd=250),290,320, add=T,col='red')
curve(f_photo(kopt = 20, Ha=100, Topt=305, Tk=x, Hd=300),290,320, add=T,col='red')

dat800 %>% 
  ggplot(data=., aes(Tk, Photo, color=Treat))+
  geom_point()+
  geom_smooth(se=F)
params.all %>% 
  ggplot(data=., aes(Vcmax, Tk, color=Treat))+
  geom_point()+
  geom_smooth(se=F)
params.all %>% 
  ggplot(data=., aes(Jmax, Tk, color=Treat))+
  geom_point()+
  geom_smooth(se=F)
