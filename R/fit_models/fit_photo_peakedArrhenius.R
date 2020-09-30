library(tidyverse); library(bayesplot);
library(rstan); library(brms); library(nls.multstart)
options(mc.cores = parallel::detectCores()-1)
set.seed(1); 

dat_fci<-read_csv("data/params_20190124.csv")
dat_fci <- dat_fci %>% 
  mutate(Tk = Tleaf+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0))




bprior <- 
  prior(normal(237 , 10), nlpar = kopt, lb=100, ub=400)+
  prior(normal( 0,2 ), nlpar = kwc,  lb=-30, ub=30)+
  prior(normal( 0,2 ), nlpar = kcw, lb=-30, ub=30)+
  prior(normal(0,2), nlpar = kww, lb=-30,ub=30)+
  prior(normal(310,10), nlpar = Topt,  lb=300, ub=320)+
  prior(normal( 0,2 ), nlpar = Tcw,  lb=-5, ub=5)+
  prior(normal( 0,2 ), nlpar = Twc, lb=-5, ub=5)+
  prior(normal( 0,2 ), nlpar = Tww, lb=-5, ub=5)+
  prior(normal( 80,10 ), nlpar = Ha,  lb=10, ub=250)+
  prior(normal( 0,2 ), nlpar = Hawc,  lb=-20, ub=20)+
  prior(normal( 0,2 ), nlpar = Hacw, lb=-20, ub=20)+
  prior(normal( 0,2 ), nlpar = Haww, lb=-20, ub=20)

f <- bf(Vcmax ~ 
      (kopt + kcw*cw + kwc*wc + kww*ww) * 
      ((200 * (2.718282^(((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
      (Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
        (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))) / 
         (200 - ((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
    (1-(2.718282^((200*(Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
      (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww)))))))),
    kopt + kcw + kwc + kww + 
      Topt + Tcw + Twc + Tww + 
      Ha + Hacw + Hawc + Haww ~ 1, 
    nl=TRUE)

make_stancode(f, prior=bprior, family=gaussian(),
              data=dat_fci)
parr_p270 <- brm(f,
           data = dat_fci, 
           prior = bprior, 
           # sample_prior = 'only')
           algorithm = 'sampling',
           control = list(adapt_delta=0.999)
)

