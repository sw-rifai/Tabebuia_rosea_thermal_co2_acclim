library(nls.multstart); library(rstan); library(tidyverse)
library(brms)
options(mc.cores = parallel::detectCores()-1)
set.seed(1); 

params.all<-read_csv("data/params.all.28032018.csv")
params.all$Tk<-params.all$Tleaf+273.15

# create empty frame for Medlyn parameters
pars.med<-NULL

# treatments: (c.c, w.w, c.w, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
treat<-as.character(unique(params.all$Treat))

params.all <- params.all %>% 
          mutate(phase1 = ifelse(Treat %in% c("c.c","c.w"), 0,1),
                 phase2 = ifelse(Treat %in% c("c.c","w.c"), 0,1), 
                 phase3 = ifelse(Treat == "w.w",1,0)) %>% 
  mutate(id = row_number()) %>% 
  mutate(plant_id = paste0("plant_",Treat,"_",Plant))


# Vcmax Stan fits - vary 3params by Treatment --------------------------------------------------

bprior <- 
  prior(normal(237 , 10), nlpar = kopt, lb=100, ub=400)+
  prior(normal( 0,2 ), nlpar = kopt1,  lb=-30, ub=30)+
  prior(normal( 0,2 ), nlpar = kopt2, lb=-30, ub=30)+
  prior(normal(0,2), nlpar = kopt3, lb=-30,ub=30)+
  prior(normal(310,10), nlpar = Topt,  lb=300, ub=320)+
  prior(normal( 0,2 ), nlpar = Topt1,  lb=-5, ub=5)+
  prior(normal( 0,2 ), nlpar = Topt2, lb=-5, ub=5)+
  prior(normal( 0,2 ), nlpar = Topt3, lb=-5, ub=5)+
  prior(normal( 80,10 ), nlpar = Ha,  lb=10, ub=250)+
  prior(normal( 0,2 ), nlpar = Ha1,  lb=-20, ub=20)+
  prior(normal( 0,2 ), nlpar = Ha2, lb=-20, ub=20)

# f <- bf(Vcmax ~ 
#           (kopt + kopt1*phase1 + kopt2*phase2) * 
#           ((200 * (2.718282^(((Ha + Ha1*phase1 + Ha2*phase2)*(Tk-(Topt + Topt1*phase1 + Topt2*phase2)))/(Tk*0.008314*(Topt + Topt1*phase1 + Topt2*phase2))))) / 
#              (200 - ((Ha + Ha1*phase1 + Ha2*phase2)*(1-(2.718282^((200*(Tk-(Topt + Topt1*phase1 + Topt2*phase2)))/(Tk*0.008314*(Topt + Topt1*phase1 + Topt2*phase2)))))))),
#         kopt + kopt1 + kopt2 + Topt + Topt1 + Topt2 + Ha + Ha1 + Ha2 ~ 1, 
#         nl=TRUE)
f <- bf(Vcmax ~ 
          (kopt + kopt1*phase1 + kopt2*phase2 + kopt3*phase3) * 
          ((200 * (2.718282^(((Ha + Ha1*phase1 + Ha2*phase2)*(Tk-(Topt + Topt1*phase1 + Topt2*phase2 +Topt3*phase3)))/(Tk*0.008314*(Topt + Topt1*phase1 + Topt2*phase2+Topt3*phase3))))) / 
             (200 - ((Ha + Ha1*phase1 + Ha2*phase2)*(1-(2.718282^((200*(Tk-(Topt + Topt1*phase1 + Topt2*phase2+Topt3*phase3)))/(Tk*0.008314*(Topt + Topt1*phase1 + Topt2*phase2+Topt3*phase3)))))))),
        kopt + kopt1 + kopt2 + kopt3 +
          Topt + Topt1 + Topt2 + Topt3 + Ha + Ha1 + Ha2 ~ 1, 
        nl=TRUE)

make_stancode(f, prior=bprior, family=gaussian(),
              data=params.all)
fit <- brm(f,
           data = params.all, 
           prior = bprior, 
           # sample_prior = 'only')
           algorithm = 'sampling',
           control = list(adapt_delta=0.999),
           chains = 4,
           iter = 1000
)
summary(fit)
plot(fit)
bayes_R2(fit)
fit3 <- update(fit, 
               algorithm='sampling', 
               chains=4, 
               iter = 1000,
               control=list(adapt_delta=0.999)
)
bayes_R2(fit3)
summary(fit3)
plot(fit3)



expand_grid(Tk=seq(min(params.all$Tk),max(params.all$Tk),length.out = 1000), 
            Treat = treat) %>% 
  mutate(phase1 = ifelse(Treat %in% c("c.c","c.w"), 0,1),
         phase2 = ifelse(Treat %in% c("c.c","w.c"), 0,1), 
         phase3 = ifelse(Treat == "w.w",1,0)) %>% 
  mutate(pred = predict(fit, newdata=.)[,'Estimate']) %>% 
  ggplot(data=., aes(Tk,pred,color=Treat))+
  geom_line()+
  geom_point(data=params.all, aes(Tk, Vcmax,color=Treat))+
  geom_smooth(data=params.all, aes(Tk, Vcmax,color=Treat), 
              se=F)+
  scale_color_viridis_d(option='B', end=0.9)+
  scale_y_continuous(limits=c(0,500))+
  theme_linedraw()


# Vcmax Stan fits - only vary Tmax --------------------------------------------------

prior_Topt <- 
  prior(normal(237 , 10), nlpar = kopt, lb=100, ub=400)+
  prior(normal( 0,2 ), nlpar = kopt1,  lb=-30, ub=30)+
  prior(normal( 0,2 ), nlpar = kopt2, lb=-30, ub=30)+
  prior(normal(0,2), nlpar = kopt3, lb=-30,ub=30)+
  prior(normal(310,10), nlpar = Topt,  lb=300, ub=320)+
  prior(normal( 0,2 ), nlpar = Topt1,  lb=-5, ub=5)+
  prior(normal( 0,2 ), nlpar = Topt2, lb=-5, ub=5)+
  prior(normal( 0,2 ), nlpar = Topt3, lb=-5, ub=5)+
  prior(normal( 80,10 ), nlpar = Ha,  lb=10, ub=250)+
  prior(normal(200,10), nlpar=Hd, lb=150, ub=250)

f_Topt <- bf(Vcmax ~ 
          (kopt + kopt1*phase1 + kopt2*phase2 + kopt3*phase3) * 
          ((Hd * (2.718282^(((Ha)*(Tk-(Topt + Topt1*phase1 + Topt2*phase2 +Topt3*phase3)))/(Tk*0.008314*(Topt + Topt1*phase1 + Topt2*phase2+Topt3*phase3))))) / 
             (Hd - ((Ha)*(1-(2.718282^((Hd*(Tk-(Topt + Topt1*phase1 + Topt2*phase2+Topt3*phase3)))/(Tk*0.008314*(Topt + Topt1*phase1 + Topt2*phase2+Topt3*phase3)))))))),
        kopt + kopt1 + kopt2 + kopt3 +
          Topt + Topt1 + Topt2 + Topt3 ~ 1, 
          Hd + Ha ~ (1|plant_id),
        nl=TRUE)

make_stancode(f_Topt, prior=prior_Topt, family=gaussian(),
              data=params.all)
fit_Topt <- brm(f_Topt,
           data = params.all, 
           prior = prior_Topt, 
           # sample_prior = 'only')
           algorithm = 'sampling',
           control = list(adapt_delta=0.999, 
                          max_treedepth=15),
           chains = 4,
           iter = 1000
)
summary(fit_Topt)
plot(fit_Topt)
bayes_R2(fit_Topt)

params.all %>% 
  group_by(Treat) %>% 
  filter(Vcmax==max(Vcmax)) %>% 
  select(Treat, id)
fit_Topt <- update(fit_Topt, 
                   newdata = params.all %>% 
                     filter(id != 167), 
                   control = list(adapt_delta=0.999, 
                                  max_treedepth=15))
summary(fit_Topt)

expand_grid(Tk=seq(min(params.all$Tk),max(params.all$Tk),length.out = 1000), 
            Treat = treat) %>% 
  mutate(phase1 = ifelse(Treat %in% c("c.c","c.w"), 0,1),
         phase2 = ifelse(Treat %in% c("c.c","w.c"), 0,1), 
         phase3 = ifelse(Treat == "w.w",1,0)) %>% 
  mutate(pred = predict(fit_Topt, newdata=.,re_formula = NA)[,'Estimate']) %>% 
  ggplot(data=., aes(Tk,pred,color=Treat))+
  geom_line()+
  geom_point(data=params.all, aes(Tk, Vcmax,color=Treat))+
  # geom_smooth(data=params.all, aes(Tk, Vcmax,color=Treat), 
  #             se=F)+
  scale_color_viridis_d(option='B', end=0.9)+
  scale_y_continuous(limits=c(0,500))+
  theme_linedraw()
#**************************************************************************
# END SECTION *************************************************************
#**************************************************************************










(kopt + kopt1*phase1 + kopt2*phase2)
(Topt + Topt1*phase1 + Topt2*phase2)
(Ha + Ha1*phase1 + Ha2*phase2)

fn_opt <- function(Tk, phase1, phase2, 
                   kopt,kopt1,kopt2, 
                   Ha, Ha1,Ha2, 
                   Topt, Topt1, Topt2){
  (kopt + kopt1*phase1 + kopt2*phase2) * 
    ((200 * (2.718282^(((Ha + Ha1*phase1 + Ha2*phase2)*(Tk-(Topt + Topt1*phase1 + Topt2*phase2)))/(Tk*0.008314*(Topt + Topt1*phase1 + Topt2*phase2))))) / 
            (200 - ((Ha + Ha1*phase1 + Ha2*phase2)*(1-(2.718282^((200*(Tk-(Topt + Topt1*phase1 + Topt2*phase2)))/(Tk*0.008314*(Topt + Topt1*phase1 + Topt2*phase2))))))))
}

mod_grad_arrh <- deriv(
  body(fn_opt)[[2]], 
  namevec = c("kopt","kopt1","kopt2",
              "Ha","Ha1","Ha2",
              "Topt","Topt1","Topt2"), 
  function.arg = fn_opt
)

#NLS multistart
nfit <- nls_multstart(Vcmax ~ mod_grad_arrh(Tk, phase1, phase2,
  kopt,kopt1,kopt2,Ha,Ha1,Ha2,Topt,Topt1,Topt2),
          data = params.all,
          iter = 1000,
          start_lower = c(kopt = 160,kopt1=0,kopt2=0,
                          Ha = 50, Ha1=0,Ha2=0,
                          Topt = 300, Topt1=0,Topt2=0),
          start_upper = c(kopt = 300, kopt1=1, kopt2=1, 
                          Ha = 300, Ha1=1, Ha2=1, 
                          Topt = 330, Topt1=1,Topt2=1),
          supp_errors = 'Y',
          na.action = na.omit,
          #convergence_count = 500,
                # "kopt""kopt1""kopt2" "Ha""Ha1""Ha2" "Topt""Topt1""Topt2"
          lower = c(100, -100, -100,    50, -5, -50,     300, -10, -10), 
          upper = c(500, 100, 100,      100, 50, 50,       330, 10, 10) 
)

summary(nfit)

expand_grid(Tk=seq(min(params.all$Tk),max(params.all$Tk),length.out = 100), 
            Treat = treat) %>% 
      mutate(phase1 = ifelse(Treat %in% c("c.c","c.w"), 0,1),
       phase2 = ifelse(Treat %in% c("c.c","w.c"), 0,1)) %>% 
  mutate(pred = predict(nfit, newdata=.)) %>% 
  ggplot(data=., aes(Tk,pred,color=Treat))+
  geom_line()+
  scale_fill_viridis_d()




library(lme4)

fit2 <- nlmer(Vcmax ~ mod_grad_arrh(Tk, phase1, phase2,
                     kopt,kopt1,kopt2,Ha,Ha1,Ha2,Topt,Topt1,Topt2) ~
                kopt|plant_id, 
      data=params.all %>% mutate(plant_id = paste0(Treat,Plant)), 
      start=coef(nfit))
summary(fit2)
      # start = c(kopt = 160,kopt1=0,kopt2=0,
      #                 Ha = 50, Ha1=0,Ha2=0,
      #                 Topt = 300, Topt1=0,Topt2=0))

predict.merMod(fit2, newdata=params.all, re.form=~0, 
        type='link')

expand_grid(Tk=seq(min(params.all$Tk),max(params.all$Tk),length.out = 100), 
            Treat = treat) %>% 
  mutate(phase1 = ifelse(Treat %in% c("c.c","c.w"), 0,1),
         phase2 = ifelse(Treat %in% c("c.c","w.c"), 0,1)) %>% 
  mutate(pred = predict(fit2, newdata=., 
                        re.form=NA)) %>% 
  ggplot(data=., aes(Tk,pred,color=Treat))+
  geom_line()+
  scale_fill_viridis_d()




