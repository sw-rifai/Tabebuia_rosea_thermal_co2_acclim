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


#*****************************************************************************
# Fit Parabolic function Photo 270
#*****************************************************************************
np270_2 <- nls_multstart(P270~ (kopt+k_cw*cw+k_wc*wc+k_ww*ww) - 
                           b*(Tk-(Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))**2, 
                         data = dat_fci, 
                         iter = 100, 
                         supp_errors = 'Y',
                         start_lower = list(kopt = 14,k_cw=-1,k_wc=-1,k_ww=-1, b=0, 
                                            Topt=273+10,Topt_cw=0,Topt_wc=0,Topt_ww=0), 
                         start_upper = list(kopt = 30,k_cw=0,k_wc=0,k_ww=0, b=0.5, 
                                            Topt=273+50,Topt_cw=1,Topt_wc=1,Topt_ww=1))

summary(np270_2, c(0.1,0.5,0.9))

set.seed(3)
bprior <- prior(normal(20,5), nlpar=kopt)+
  prior(normal(0,5), nlpar=kwc)+
  prior(normal(0,5), nlpar=kcw)+
  prior(normal(0,5), nlpar=kww)+
  prior(normal(0.25,0.25), nlpar=b, lb=0)+
  prior(normal(300,10), nlpar=Topt, lb=273)+
  prior(normal(0,5), nlpar=Twc)+
  prior(normal(0,5), nlpar=Tcw)+
  prior(normal(0,5), nlpar=Tww)

f <- bf(P270~kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2, 
        kopt+kwc+kcw+kww+b+Topt+Twc+Tcw+Tww ~ 1,
        nl=TRUE)
make_stancode(f, prior=bprior, family=gaussian(),
              data=dat_fci)

fit_p270 <- brm(f,
                data = dat_fci, 
                prior = bprior, 
                # sample_prior = 'only'
                algorithm = 'sampling',
                control = list(adapt_delta=0.99,
                               max_treedepth=13),
                # chains = 3,
                # iter = 250
)

summary(fit_p270, prob=0.8)$fixed
plot(fit_p270,ask = F)
bayes_R2(fit_p270)
sm_p270 <- broom.mixed::tidy(fit_p270, conf.level=0.8, conf.method="HPDinterval")

#*****************************************************************************
# Fit Parabolic function Photo 505 
#*****************************************************************************
np505_2 <- nls_multstart(P505~ (kopt+k_cw*cw+k_wc*wc+k_ww*ww) - 
                           b*(Tk-(Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))**2, 
                         data = dat505, 
                         iter = 100, 
                         supp_errors = 'Y',
                         start_lower = list(kopt = 14,k_cw=-1,k_wc=-1,k_ww=-1, b=0, 
                                            Topt=273+10,Topt_cw=0,Topt_wc=0,Topt_ww=0), 
                         start_upper = list(kopt = 30,k_cw=0,k_wc=0,k_ww=0, b=0.5, 
                                            Topt=273+50,Topt_cw=1,Topt_wc=1,Topt_ww=1))

summary(np505_2)

set.seed(3)
bprior <- prior(normal(30,5), nlpar=kopt)+
  prior(normal(0,5), nlpar=kwc)+
  prior(normal(0,5), nlpar=kcw)+
  prior(normal(0,5), nlpar=kww)+
  prior(normal(0.25,0.25), nlpar=b, lb=0)+
  prior(normal(305,10), nlpar=Topt, lb=273)+
  prior(normal(0,5), nlpar=Twc)+
  prior(normal(0,5), nlpar=Tcw)+
  prior(normal(0,5), nlpar=Tww)

f <- bf(P505~kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2, 
        kopt+kwc+kcw+kww+b+Topt+Twc+Tcw+Tww ~ 1,
        nl=TRUE)
make_stancode(f, prior=bprior, family=gaussian(),
              data=dat_fci)

fit_p505 <- brm(f,
                data = dat_fci, 
                prior = bprior, 
                # sample_prior = 'only'
                algorithm = 'sampling',
                control = list(adapt_delta=0.99,
                               max_treedepth=13),
                # chains = 3,
                # iter = 250
)
summary(fit_p505,prob=c(0.8))$fixed
plot(fit_p505, ask=F)
bayes_R2(fit_p505)
prior_summary(fit_p505)

sm_p505 <- broom.mixed::tidy(fit_p505, conf.level=0.8, conf.method="HPDinterval")



# Plotting ----------------------------------------------------------------


fn <- function(x){
  x <- gsub("_Intercept","",x)
  x <- gsub("b_","",x)
  x
}
fn("b")
vec_labels <- c("c.c"="Control", 
                "c.w"="paste(Control %->% \" Treatment\")",
                "w.c"="paste(Treatment %->% \" Control\")", 
                "w.w"="Treatment")

#************************************************************************
# FOR THE TOP ROW of P270 -------------------
#************************************************************************

p1 <- fit_p270 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kcw,kwc,kww,Topt,Tcw,Twc,Tww,b),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","w.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Photo, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("c.c","w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P270,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Tww_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Tww_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Tww_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  scale_color_manual("",
                     values=c("c.c"='blue', "w.w"='#cf0000'), 
                     labels=c("c.c"="Control",
                              "c.w"="paste(Control %->% \" Treatment\")",
                              "w.c"="paste(Treatment %->% \" Control\")", 
                              "w.w"="Treatment"))+
  theme_bw()+
  labs(y=bquote(paste('P'[270],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,40)+
  scale_x_continuous(limits=c(25,46), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = c(0.95,0.95), 
        legend.justification = c(0.95,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p1


#************************************************************************
# FOR THE MIDDLE ROW of P270 ---------------------
#************************************************************************
p2 <- fit_p270 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kcw,kwc,kww,Topt,Tcw,Twc,Tww,b),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","c.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Photo, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("c.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P270,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("c.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P270,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F,shape=1)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Tcw_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Tcw_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Tcw_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  scale_color_manual("", 
                     values=c("blue","red"), 
                     breaks=c("c.c","c.w"),
                     labels=c("Control",
                              sprintf('Control\u2192Treatment')
                     ),
                     guide=guide_legend(override.aes = list(shape=c(20,1), 
                                                            size=c(6,3)))
  )+
  theme_bw()+
  labs(y=bquote(paste('P'[270],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,40)+
  scale_x_continuous(limits=c(25,46), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = c(0.95,0.95), 
        legend.justification = c(0.95,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p2


#************************************************************************
# FOR THE BOTTOM ROW of P270 -------------------------
#************************************************************************
p3 <- fit_p270 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kcw,kwc,kww,Topt,Tcw,Twc,Tww,b),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("w.w","w.c")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Photo, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P270,color=Treat),
             size=5,alpha=0.5,
             inherit.aes = F, pch=20)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("w.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P270,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F, pch=1)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Tww_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Tww_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Tww_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Twc_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Twc_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p270,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p270,term=="Twc_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  scale_color_manual("", 
                     values=c("red","blue"), 
                     breaks=c("w.w","w.c"),
                     labels=c("Treatment",
                              sprintf('Treatment\u2192Control')
                     ), 
                     guide=guide_legend(override.aes = list(shape=c(20,1), 
                                                            size=c(6,3)))
                     # parse("paste(Control %->% \" Treatment\")"))
  )+
  theme_bw()+
  labs(y=bquote(paste('P'[270],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,40)+
  scale_x_continuous(limits=c(25,46), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = c(0.95,0.95), 
        legend.justification = c(0.95,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.title.x = element_text(size=18),
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20)); p3




library(patchwork)
p1/p2/p3

ggsave(p1/p2/p3,
       filename = paste0("figures/photo_Topt/p270_parabolicFit_CI80_",Sys.Date(),".png"),
       width = 120, height=3*100, units='mm')


#************************************************************************
# FOR THE TOP ROW of P505 --------------------------------------
#************************************************************************

p4 <- fit_p505 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kcw,kwc,kww,Topt,Tcw,Twc,Tww,b),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","w.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Photo, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("c.c","w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P505,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Tww_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Tww_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Tww_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  scale_color_manual("",
                     values=c("c.c"='blue', "w.w"='#cf0000'), 
                     labels=c("c.c"="Control",
                              "c.w"="paste(Control %->% \" Treatment\")",
                              "w.c"="paste(Treatment %->% \" Control\")", 
                              "w.w"="Treatment"))+
  theme_bw()+
  labs(y=bquote(paste('P'[505],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(position = 'right', limits = c(0,40))+
  scale_x_continuous(limits=c(25,46), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = 'none',#c(0.95,0.95), 
        legend.justification = c(0.95,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title.y = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p4


#************************************************************************
# FOR THE MIDDLE ROW of P505 ----------------------
#************************************************************************
p5 <- fit_p505 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kcw,kwc,kww,Topt,Tcw,Twc,Tww,b),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","c.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Photo, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("c.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P505,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("c.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P505,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F,shape=1)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Tcw_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Tcw_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Tcw_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  scale_color_manual("", 
                     values=c("blue","red"), 
                     breaks=c("c.c","c.w"),
                     labels=c("Control",
                              sprintf('Control\u2192Treatment')
                     ),
                     guide=guide_legend(override.aes = list(shape=c(20,1), 
                                                            size=c(6,3)))
  )+
  theme_bw()+
  labs(y=bquote(paste('P'[505],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(position = 'right', limits = c(0,40))+
  scale_x_continuous(limits=c(25,46), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = 'none',#c(0.95,0.95), 
        legend.justification = c(0.95,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title.y = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p5


#************************************************************************
# FOR THE BOTTOM ROW of P505 -----------------------------
#************************************************************************
p6 <- fit_p505 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kcw,kwc,kww,Topt,Tcw,Twc,Tww,b),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("w.w","w.c")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Photo, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P505,color=Treat),
             size=5,alpha=0.5,
             inherit.aes = F, pch=20)+
  geom_point(data=dat_fci %>% 
               filter(Treat%in%c("w.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P505,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F, pch=1)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Tww_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Tww_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Tww_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Twc_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Twc_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p505,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p505,term=="Twc_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  scale_color_manual("", 
                     values=c("red","blue"), 
                     breaks=c("w.w","w.c"),
                     labels=c("Treatment",
                              sprintf('Treatment\u2192Control')
                     ), 
                     guide=guide_legend(override.aes = list(shape=c(20,1), 
                                                            size=c(6,3)))
                     # parse("paste(Control %->% \" Treatment\")"))
  )+
  theme_bw()+
  labs(y=bquote(paste('P'[505],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(position = 'right', limits = c(0,40))+
  scale_x_continuous(limits=c(25,46), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = 'none',# c(0.95,0.95), 
        legend.justification = c(0.95,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18), 
        axis.title.y = element_text(size=20), 
        plot.margin = margin(l=20)); p6




library(patchwork)
p4/p5/p6

ggsave(p4/p5/p6,
       filename = paste0("figures/photo_Topt/p505_parabolicFit_CI80_",Sys.Date(),".png"),
       width = 120, height=3*100, units='mm')




ggsave((p1/p2/p3)|(p4/p5/p6),
       filename = paste0("figures/photo_Topt/p270_p505_parabolicFit_CI80_",Sys.Date(),".png"),
       width = 2*120, height=3*100, units='mm')
ggsave((p1/p2/p3)|(p4/p5/p6),
       filename = paste0("figures/photo_Topt/p270_p505_parabolicFit_CI80_",Sys.Date(),".pdf"),
       width = 2*120, height=3*100, units='mm')
