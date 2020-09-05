library(tidyverse); library(bayesplot);
library(rstan); library(brms); library(nls.multstart)
# library(xfun); library(brms); 
options(mc.cores = parallel::detectCores()-1)
set.seed(1); 

dat400 <- read_csv("data/P400.2904.csv")
dat800 <- read_csv("data/P800.2904.csv")

dat400 <- dat400 %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0)) %>% 
  mutate(p1_w = ifelse(Treat %in% c("w.c","w.w"), 1,0)) %>% 
  mutate(p2_w = ifelse(Treat %in% c("c.w","w.w"), 1,0))


#*****************************************************************************
# Compare functional forms using Photo 400 
#*****************************************************************************

# Peaked Arrhenius 
np400_1 <- nls_multstart(Photo ~ 
    (kopt+k_cw*cw+k_wc*wc+k_ww*ww) *
      ((200 * (2.718282^(((Ha+Ha_cw*cw+Ha_wc*wc+Ha_ww*ww)*(Tk-((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))))/(Tk*0.008314*((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww)))))) /
        (200 - ((Ha+Ha_cw*cw+Ha_wc*wc+Ha_ww*ww)*(1-(2.718282^((200*(Tk-((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))))/(Tk*0.008314*((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))))))))),
            data = dat400,
            iter = 100,
            start_lower = c(kopt = 20, k_cw=-2,k_wc=-2,k_ww=-2, 
                            Ha = 2,Ha_cw=0,Ha_wc=0, Ha_ww=0, 
                            Topt = 273+25, Topt_cw=0, Topt_wc=0, Topt_ww=0),
            start_upper = c(kopt = 36,k_cw=2,k_wc=2,k_ww=2, 
                            Ha = 100, Ha_cw=2, Ha_wc=2, Ha_ww=2, 
                            Topt = 273+40, Topt_cw=1, Topt_wc=1, Topt_ww=1),
            # supp_errors = 'Y',
            na.action = na.omit,
            convergence_count = 500,
            lower = c(kopt = 19.5,k_cw=-10,k_wc=-10,k_ww=-10, 
                      Ha = 2,Ha_cw=-20, Ha_wc=-20, Ha_ww=-20, 
                      Topt = 273+20, Topt_cw=-3, Topt_wc=-3, Topt_ww=-3),
            upper = c(kopt=40,k_cw=10,k_wc=10,k_ww=10, 
                      Ha=120,Ha_cw=200, Ha_wc=200, Ha_ww=200, 
                      Topt=273+40, Topt_cw=5, Topt_wc=5, Topt_ww=10))
summary(np400_1)
yardstick::rsq_trad_vec(dat400$Photo, predict(np400_1, newdata=dat400))
dat400 %>% 
  mutate(pred=predict(np400_1, newdata=.)) %>% 
  ggplot(data=.,aes(pred,Photo,color=Treat))+
  geom_abline(aes(intercept=0,slope=1),col='#cf0000')+
  geom_point()+
  geom_smooth(method='lm',se=F)+
  coord_equal()+
  scale_color_viridis_d(option='A',end=0.8)

p_pa <- expand_grid(Treat=unique(dat400$Treat), 
            Tc=seq(25,47,length.out = 100)) %>%
  mutate(Tk=Tc+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0)) %>% 
  mutate(p1_w = ifelse(Treat %in% c("w.c","w.w"), 1,0)) %>% 
  mutate(p2_w = ifelse(Treat %in% c("c.w","w.w"), 1,0)) %>% 
  mutate(pred=predict(nfit_1, newdata=.)) %>% 
  ggplot(data=., aes(Tc, pred, color=Treat))+
  geom_line(lwd=1)+
  geom_point(data=dat400, aes(Tk-273.15, Photo, color=Treat))+
  scale_color_viridis_d(end=0.9)+
  labs(title='Peaked Arrhenius'); p_pa


#--- concave parabola ---------------------------------------------
np400_2 <- nls_multstart(Photo~ (kopt+k_cw*cw+k_wc_*wc+k_ww*ww) - 
                           b*(Tk-(Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))**2, 
                         data = dat400, 
                         iter = 100, 
                         supp_errors = 'Y',
                         start_lower = list(kopt = 14,k_cw=-1,k_wc=-1,k_ww=-1, b=0, 
                                            Topt=273+10,Topt_cw=0,Topt_wc=0,Topt_ww=0), 
                         start_upper = list(kopt = 30,k_cw=0,k_wc=0,k_ww=0, b=0.5, 
                                            Topt=273+50,Topt_cw=1,Topt_wc=1,Topt_ww=1))

summary(np400_2)
dat400 %>% 
  mutate(pred=predict(np400_2, newdata=.)) %>% 
  ggplot(data=.,aes(pred,Photo))+
  geom_abline(aes(intercept=0,slope=1),col='#cf0000')+
  geom_point()+
  geom_smooth(method='lm',se=F)+
  coord_equal()
yardstick::rsq_trad_vec(dat400$Photo, predict(np400_2))

p_parab <- expand_grid(Treat=unique(dat400$Treat), 
            Tc=seq(25,47,length.out = 100)) %>%
  mutate(Tk=Tc+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0)) %>% 
  mutate(p1_w = ifelse(Treat %in% c("w.c","w.w"), 1,0)) %>% 
  mutate(p2_w = ifelse(Treat %in% c("c.w","w.w"), 1,0)) %>% 
  mutate(pred=predict(np400_2, newdata=.)) %>% 
  ggplot(data=., aes(Tc, pred, color=Treat))+
  geom_line(lwd=1)+
  geom_point(data=dat400, aes(Tk-273.15, Photo, color=Treat))+
  scale_color_viridis_d(end=0.9)+
  geom_hline(aes(yintercept=coef(np400_2)['kopt']))+
  geom_hline(aes(yintercept=coef(np400_2)['kopt']+coef(np400_2)['k_wc']))+
  geom_hline(aes(yintercept=coef(np400_2)['kopt']+coef(np400_2)['k_cw']))+
  geom_hline(aes(yintercept=coef(np400_2)['kopt']+coef(np400_2)['k_ww']))+
  labs(title='Concave Parabola'); p_parab



# Ricker Function ---------------------------------------------------------
np400_3 <- nls_multstart(Photo~(a+a_wc*wc+a_cw*cw+a_ww*ww)*
                           Tk*exp(-(b+b_wc*wc+b_cw*cw+b_ww*ww)*Tk), 
                         data=dat400, 
                         iter=1000, 
                         supp_errors = 'Y', 
                         control=list(maxiter=200),
                         # control = nls.lm.control(list(maxiter=50)),
                         start_lower=c(a=0,a_wc=0,a_cw=0,a_ww=0,
                                       b=0.003,b_wc=0,b_cw=0,b_ww=0), 
                         start_upper=c(a=0.25,a_wc=0.01,a_cw=0.01,a_ww=0.01,
                                       b=0.0033,b_wc=0,b_cw=0,b_ww=0), 
                         lower=c(a=0,a_wc=-0.05, a_cw=-0.05, a_ww=-0.05, 
                                 b=0.003, b_wc=-0.0005, b_cw=-0.0005, b_ww=-0.0005), 
                         upper=c(a=0.3,a_wc=0.05, a_cw=0.05, a_ww=0.05, 
                                 b=0.00339, b_wc=0.0005, b_cw=0.0005, b_ww=0.0005)
                         )
summary(np400_3)              
1/coef(np400_3)['b'] # Topt
exp(-1)*coef(np400_3)['a']/coef(np400_3)['b'] # kopt
dat400 %>% 
  mutate(pred=predict(np400_3, newdata=.)) %>% 
  ggplot(data=.,aes(pred,Photo))+
  geom_abline(aes(intercept=0,slope=1),col='#cf0000')+
  geom_point()+
  geom_smooth(method='lm',se=F)+
  coord_equal()
yardstick::rsq_trad_vec(dat400$Photo, predict(np400_3))

p_ricker <- expand_grid(Treat=unique(dat400$Treat), 
            Tc=seq(25,47,length.out = 100)) %>%
  mutate(Tk=Tc+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0)) %>% 
  mutate(p1_w = ifelse(Treat %in% c("w.c","w.w"), 1,0)) %>% 
  mutate(p2_w = ifelse(Treat %in% c("c.w","w.w"), 1,0)) %>% 
  mutate(pred=predict(np400_3, newdata=.)) %>% 
  ggplot(data=., aes(Tc, pred, color=Treat))+
  geom_line(lwd=1)+
  geom_point(data=dat400, aes(Tk-273.15, Photo, color=Treat))+
  scale_color_viridis_d(end=0.9)+
  labs(title='Ricker'); p_ricker


bbmle::AICtab(np400_1, np400_2, np400_3)

library(patchwork)
p_pa+p_parab+p_ricker+plot_layout(guides = 'collect')
ggsave(filename = "figures/compare_photo400_functional_forms.png", 
       width=20, height = 8, units='cm')
################################################################################
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

np400_cc <- nls_multstart(Photo~ (kopt+k_p1*p1_w+k_p2*p2_w) - 
                            b*(Tk-(Topt+Topt_p1*p1_w+Topt_p2*p2_w))**2, 
                          data = dat400, 
                          iter = 10000, 
                          supp_errors = 'Y',
                          start_lower = list(kopt = 14,k_p1=-1,k_p2=-1, b=0, Topt=273+10, 
                                             Topt_p1=0,Topt_p2=0), 
                          start_upper = list(kopt = 30,k_p1=0,k_p2=0, b=0.5, Topt=273+50, 
                                             Topt_p1=1,Topt_p2=1))
summary(np400_cc)
dat400 %>% 
  mutate(pred=predict(np400_cc, newdata=.)) %>% 
  ggplot(data=.,aes(pred,Photo))+
  geom_abline(aes(intercept=0,slope=1),col='#cf0000')+
  geom_point()+
  geom_smooth(method='lm',se=F)+
  coord_equal()
yardstick::rsq_trad_vec(dat400$Photo, predict(np400_cc))






np400_cc_2 <- nls.multstart::nls_multstart(Photo~ b*(Tk-Topt)*(1-exp(c*(Tk-Topt))), 
                data = dat400, 
                iter=1000,
                supp_errors = 'Y',
                start_lower = list(b=1,  Topt=300, c=-30), 
                start_upper = list(b=10, Topt=320, c=30))
np400_cc_2

summary(np400_cc_2)
dat400 %>% 
  mutate(pred=predict(np400_cc_2, newdata=.)) %>% 
  ggplot(data=.,aes(pred,Photo))+
  geom_abline(aes(intercept=0,slope=1),col='#cf0000')+
  geom_point()+
  geom_smooth(method='lm',se=F)+
  coord_equal()
yardstick::rsq_trad_vec(dat400$Photo, predict(np400_cc_2))


set.seed(3)
bprior <- prior(normal(20,5), nlpar=Aopt, lb=10)+
          prior(normal(0.5,0.25), nlpar=b, lb=0)+
          prior(normal(273+25,5), nlpar=Topt, lb=283)
f <- bf(Photo~ Aopt - b*(Tk-Topt)**2, 
        Aopt + b + Topt~1,
        nl=TRUE)
make_stancode(f, prior=bprior, family=gaussian(),
              data=dat400)

p400_cc <- brm(f,
             data = dat400, 
             prior = bprior, 
             # sample_prior = 'only'
             algorithm = 'sampling',
             control = list(adapt_delta=0.99,
                            max_treedepth=13),
             # chains = 3,
             # iter = 250
)
# p400_cc <- update(p400_cc, recompile = F, chains=4, iter=6000)
p400_cc
plot(p400_cc)
bayes_R2(p400_cc)


set.seed(3)
bprior_2 <- prior(normal(20,5), nlpar=Aopt, lb=10)+
  prior(normal(0.5,0.25), nlpar=b, lb=0)+
  prior(normal(273+25,5), nlpar=Topt, lb=283)
f_2 <- bf(Photo~ Aopt - b*(Tk-Topt)**2, 
        Aopt + b + Topt~1,
        nl=TRUE)
make_stancode(f_2, prior=bprior_2, family=gaussian(),
              data=dat400)

p400_cc_2 <- brm(f,
               data = dat400, 
               prior = bprior, 
               # sample_prior = 'only'
               algorithm = 'sampling',
               control = list(adapt_delta=0.99,
                              max_treedepth=13),
               # chains = 3,
               # iter = 250
)















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
dat400_cc <- dat400 %>% filter(Treat=='c.c')




f_pa <- function(kopt,k_cc,k_cw,k_wc,k_ww, 
                 Ha, Topt){
  (kopt+k_cc*cc+k_cw*cw+k_wc*wc+k_ww*ww) * 
    ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
       (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
}
mod_grad_pa <- deriv(
  body(f_pa)[[2]], 
  namevec = c("kopt","k_cc","k_cw","k_wc","k_ww", 
              "Ha","Topt"), 
  function.arg = f_pa
)






# Peaked Arrhenius with Hd estimate ---------------------------------------
nfit_2 <- nls_multstart(Photo ~ 
         (kopt+k_cw*cw+k_wc*wc+k_ww*ww) *
         ((200 * (2.718282^(((Ha+Ha_cw*cw+Ha_wc*wc+Ha_ww*ww)*(Tk-((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))))/(Tk*0.008314*((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww)))))) /
            (200 - ((Ha+Ha_cw*cw+Ha_wc*wc+Ha_ww*ww)*(1-(2.718282^((200*(Tk-((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))))/(Tk*0.008314*((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))))))))),
       data = dat400,
       iter = 100,
       supp_errors = 'Y',
       start_lower = c(kopt = 20, k_cw=-2,k_wc=-2,k_ww=-2, 
                       Ha = 2,Ha_cw=0,Ha_wc=0, Ha_ww=0, 
                       Topt = 273+25, Topt_cw=0, Topt_wc=0, Topt_ww=0),
       start_upper = c(kopt = 36,k_cw=2,k_wc=2,k_ww=2, 
                       Ha = 100, Ha_cw=2, Ha_wc=2, Ha_ww=2, 
                       Topt = 273+40, Topt_cw=1, Topt_wc=1, Topt_ww=1),
       na.action = na.omit,
       convergence_count = 500,
       lower = c(kopt = 19.5,k_cw=-10,k_wc=-10,k_ww=-10, 
                 Ha = 2,Ha_cw=-20, Ha_wc=-20, Ha_ww=-20, 
                 Topt = 273+20, Topt_cw=-3, Topt_wc=-3, Topt_ww=-3),
       upper = c(kopt=40,k_cw=10,k_wc=10,k_ww=10, 
                 Ha=120,Ha_cw=200, Ha_wc=200, Ha_ww=200, 
                 Topt=273+40, Topt_cw=5, Topt_wc=5, Topt_ww=10))
summary(nfit_2)


nfit_2 <- nlme(Photo ~ 
          (kopt+k_cw*cw+k_wc*wc+k_ww*ww) *
          ((200 * (2.718282^(((Ha+Ha_cw*cw+Ha_wc*wc+Ha_ww*ww)*(Tk-((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))))/(Tk*0.008314*((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww)))))) /
             (200 - ((Ha+Ha_cw*cw+Ha_wc*wc+Ha_ww*ww)*(1-(2.718282^((200*(Tk-((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))))/(Tk*0.008314*((Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))))))))),
          fixed=kopt+k_cw+k_wc+k_ww+Ha+Ha_cw+Ha_wc+Ha_ww+
            Topt+Topt_cw+Topt_wc+Topt_ww~1,
          random = Topt~1,
        data = dat400,
        # iter = 100,
        # supp_errors = 'Y',
        start = c(kopt = 20, k_cw=-2,k_wc=-2,k_ww=-2, 
                        Ha = 2,Ha_cw=0,Ha_wc=0, Ha_ww=0, 
                        Topt = 273+25, Topt_cw=0, Topt_wc=0, Topt_ww=0))
        # start_upper = c(kopt = 36,k_cw=2,k_wc=2,k_ww=2, 
        #                 Ha = 100, Ha_cw=2, Ha_wc=2, Ha_ww=2, 
        #                 Topt = 273+40, Topt_cw=1, Topt_wc=1, Topt_ww=1),
        # na.action = na.omit,
        # convergence_count = 500,
        # lower = c(kopt = 19.5,k_cw=-10,k_wc=-10,k_ww=-10, 
        #           Ha = 2,Ha_cw=-20, Ha_wc=-20, Ha_ww=-20, 
        #           Topt = 273+20, Topt_cw=-3, Topt_wc=-3, Topt_ww=-3),
        # upper = c(kopt=40,k_cw=10,k_wc=10,k_ww=10, 
        #           Ha=120,Ha_cw=200, Ha_wc=200, Ha_ww=200, 
        #           Topt=273+40, Topt_cw=5, Topt_wc=5, Topt_ww=10))


install.packages('plantecophys')
plantecophys::Photosyn()


library(mgcv)
fg <- gam(Photo~s(Tk, bs='cs')+s(Tk,by=Treat,bs=c('bs','fs'))+s(Plant,bs='re'), 
    data=dat400 %>% 
      mutate(Treat=as_factor(Treat), 
             Plant=as_factor(paste(Treat,Plant))), 
    method='REML')
summary(fg)
plot(fg)
mgcViz::getViz(fg) %>% plot


expand_grid(Treat=unique(dat400$Treat), 
            Tc=seq(25,47,length.out = 100), 
            Plant=1:5) %>%
  mutate(Tk=Tc+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0)) %>% 
  mutate(p1_w = ifelse(Treat %in% c("w.c","w.w"), 1,0)) %>% 
  mutate(p2_w = ifelse(Treat %in% c("c.w","w.w"), 1,0)) %>% 
  mutate(Treat=as_factor(Treat), 
         Plant=as_factor(paste(Treat,Plant))) %>%  
  mutate(pred=predict(fg, newdata=., exclude="s(Plant)")) %>% 
  ggplot(data=., aes(Tc, pred, color=Treat))+
  geom_line(lwd=1)+
  geom_point(data=dat400, aes(Tk-273.15, Photo, color=Treat))+
  scale_color_viridis_d(end=0.9)

