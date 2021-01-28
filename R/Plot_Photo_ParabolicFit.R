library(tidyverse); library(bayesplot);
library(rstan); library(brms); library(nls.multstart)
options(mc.cores = parallel::detectCores()-1, 
        pillar.sigfig=3)
options(pillar.sigfig = 5)
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

dat800 <- dat800 %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0)) %>% 
  mutate(p1_w = ifelse(Treat %in% c("w.c","w.w"), 1,0)) %>% 
  mutate(p2_w = ifelse(Treat %in% c("c.w","w.w"), 1,0))

#*****************************************************************************
# Fit Parabolic function Photo 400 ------------------------------------------
#*****************************************************************************
np400_2 <- nls_multstart(Photo~ (kopt+k_cw*cw+k_wc*wc+k_ww*ww) - 
                           b*(Tk-(Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))**2, 
                         data = dat400, 
                         iter = 100, 
                         supp_errors = 'Y',
                         start_lower = list(kopt = 14,k_cw=-1,k_wc=-1,k_ww=-1, b=0, 
                                            Topt=273+10,Topt_cw=0,Topt_wc=0,Topt_ww=0), 
                         start_upper = list(kopt = 30,k_cw=0,k_wc=0,k_ww=0, b=0.5, 
                                            Topt=273+50,Topt_cw=1,Topt_wc=1,Topt_ww=1))

summary(np400_2, c(0.1,0.5,0.9))

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

f <- bf(Photo~kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2, 
        kopt+kwc+kcw+kww+b+Topt+Twc+Tcw+Tww ~ 1,
        nl=TRUE)
make_stancode(f, prior=bprior, family=gaussian(),
              data=dat400)

parabolic_p400 <- brm(f,
               data = dat400, 
               prior = bprior, 
               # sample_prior = 'only'
               algorithm = 'sampling',
               control = list(adapt_delta=0.99,
                              max_treedepth=13),
               # chains = 3,
               # iter = 250
)
write_rds(parabolic_p400, file = paste0("outputs/parr_parabolic_p400_",Sys.Date(),".rds"))
parabolic_p400
summary(parabolic_p400, prob=0.8)$fixed
plot(parabolic_p400,ask = F)
bayes_R2(parabolic_p400)
sm_p400 <- broom.mixed::tidy(parabolic_p400, conf.level=0.8, conf.method="HPDinterval")
sm_p400
brms::pp_check(parabolic_p400, nsamples=100)


#*****************************************************************************
# Fit Parabolic function Photo 800 -------------------------------------------
#*****************************************************************************
np800_2 <- nls_multstart(Photo~ (kopt+k_cw*cw+k_wc*wc+k_ww*ww) - 
                           b*(Tk-(Topt+Topt_cw*cw+Topt_wc*wc+Topt_ww*ww))**2, 
                         data = dat800, 
                         iter = 100, 
                         supp_errors = 'Y',
                         start_lower = list(kopt = 14,k_cw=-1,k_wc=-1,k_ww=-1, b=0, 
                                            Topt=273+10,Topt_cw=0,Topt_wc=0,Topt_ww=0), 
                         start_upper = list(kopt = 30,k_cw=0,k_wc=0,k_ww=0, b=0.5, 
                                            Topt=273+50,Topt_cw=1,Topt_wc=1,Topt_ww=1))

summary(np800_2)

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

f <- bf(Photo~kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2, 
        kopt+kwc+kcw+kww+b+Topt+Twc+Tcw+Tww ~ 1,
        nl=TRUE)
make_stancode(f, prior=bprior, family=gaussian(),
              data=dat800)

parabolic_p800 <- brm(f,
                data = dat800, 
                prior = bprior, 
                # sample_prior = 'only'
                algorithm = 'sampling',
                control = list(adapt_delta=0.99,
                               max_treedepth=13),
                # chains = 3,
                # iter = 250
)
write_rds(parabolic_p800, file = paste0("outputs/parr_parabolic_p800_",Sys.Date(),".rds"))
parabolic_p800
summary(parabolic_p800,prob=c(0.8))$fixed
plot(parabolic_p800, ask=F)
bayes_R2(parabolic_p800)
prior_summary(parabolic_p800)

sm_p800 <- broom.mixed::tidy(parabolic_p800, conf.level=0.8, conf.method="HPDinterval")



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
# FOR THE TOP ROW of P400 -------------------
#************************************************************************

p1 <- parabolic_p400 %>% 
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
  geom_point(data=dat400 %>% 
               filter(Treat%in%c("c.c","w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Tww_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Tww_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Tww_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  scale_color_manual("",
                     values=c("c.c"='blue', "w.w"='#cf0000'), 
                     labels=c("c.c"="Control",
                              "c.w"="paste(Control %->% \" Treatment\")",
                              "w.c"="paste(Treatment %->% \" Control\")", 
                              "w.w"="Treatment"))+
  theme_bw()+
  labs(y=bquote(paste('P'[400],'(',mu,'mol m'^'-2','s'^'-1',')')),
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
# FOR THE MIDDLE ROW of P400 ---------------------
#************************************************************************
p2 <- parabolic_p400 %>% 
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
  geom_point(data=dat400 %>% 
               filter(Treat%in%c("c.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat400 %>% 
               filter(Treat%in%c("c.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F,shape=1)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Tcw_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Tcw_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Tcw_(Intercept)")$conf.low-273.15), 
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
  labs(y=bquote(paste('P'[400],'(',mu,'mol m'^'-2','s'^'-1',')')),
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
# FOR THE BOTTOM ROW of P400 -------------------------
#************************************************************************
p3 <- parabolic_p400 %>% 
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
  geom_point(data=dat400 %>% 
               filter(Treat%in%c("w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=5,alpha=0.5,
             inherit.aes = F, pch=20)+
  geom_point(data=dat400 %>% 
               filter(Treat%in%c("w.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F, pch=1)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Tww_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Tww_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Tww_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Twc_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Twc_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p400,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p400,term=="Twc_(Intercept)")$conf.low-273.15), 
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
  labs(y=bquote(paste('P'[400],'(',mu,'mol m'^'-2','s'^'-1',')')),
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
  filename = paste0("figures/photo_Topt/p400_parabolicFit_CI80_",Sys.Date(),".png"),
  width = 120, height=3*100, units='mm')


#************************************************************************
# FOR THE TOP ROW of P800 --------------------------------------
#************************************************************************

p4 <- parabolic_p800 %>% 
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
  geom_point(data=dat800 %>% 
               filter(Treat%in%c("c.c","w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Tww_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Tww_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Tww_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  scale_color_manual("",
                     values=c("c.c"='blue', "w.w"='#cf0000'), 
                     labels=c("c.c"="Control",
                              "c.w"="paste(Control %->% \" Treatment\")",
                              "w.c"="paste(Treatment %->% \" Control\")", 
                              "w.w"="Treatment"))+
  theme_bw()+
  labs(y=bquote(paste('P'[800],'(',mu,'mol m'^'-2','s'^'-1',')')),
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
# FOR THE MIDDLE ROW of P800 ----------------------
#************************************************************************
p5 <- parabolic_p800 %>% 
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
  geom_point(data=dat800 %>% 
               filter(Treat%in%c("c.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat800 %>% 
               filter(Treat%in%c("c.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F,shape=1)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$conf.low-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Tcw_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Tcw_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Tcw_(Intercept)")$conf.low-273.15), 
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
  labs(y=bquote(paste('P'[800],'(',mu,'mol m'^'-2','s'^'-1',')')),
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
# FOR THE BOTTOM ROW of P800 -----------------------------
#************************************************************************
p6 <- parabolic_p800 %>% 
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
  geom_point(data=dat800 %>% 
               filter(Treat%in%c("w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=5,alpha=0.5,
             inherit.aes = F, pch=20)+
  geom_point(data=dat800 %>% 
               filter(Treat%in%c("w.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Photo,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F, pch=1)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Tww_(Intercept)")$estimate-273.15), 
             color='red')+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Tww_(Intercept)")$conf.high-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Tww_(Intercept)")$conf.low-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Twc_(Intercept)")$estimate-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Twc_(Intercept)")$conf.high-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
                   filter(sm_p800,term=="Twc_(Intercept)")$conf.low-273.15), 
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
  labs(y=bquote(paste('P'[800],'(',mu,'mol m'^'-2','s'^'-1',')')),
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
       filename = paste0("figures/photo_Topt/p800_parabolicFit_CI80_",Sys.Date(),".png"),
       width = 120, height=3*100, units='mm')




ggsave((p1/p2/p3)|(p4/p5/p6),
       filename = paste0("figures/photo_Topt/p400_p800_parabolicFit_CI80_",Sys.Date(),".png"),
       width = 2*120, height=3*100, units='mm')
ggsave((p1/p2/p3)|(p4/p5/p6),
       filename = paste0("figures/photo_Topt/p400_p800_parabolicFit_CI80_",Sys.Date(),".pdf"),
       width = 2*120, height=3*100, units='mm')
