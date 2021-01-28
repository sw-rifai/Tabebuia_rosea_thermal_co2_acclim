library(tidyverse); library(bayesplot);
library(rstan); library(brms); library(nls.multstart)
options(mc.cores = parallel::detectCores()-1)
set.seed(1); 

dat400 <- read_csv("data/P400.2904.csv") %>% rename(P400=Photo)
dat800 <- read_csv("data/P800.2904.csv") %>% rename(P800=Photo)
dat400 <- dat400 %>% 
  mutate(Tk = Tleaf+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0))
dat800 <- dat800 %>% 
  mutate(Tk = Tleaf+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0))


# Ricker function: photo ~ a*Tk*exp(-b*Tk)
#*****************************************************************************
# Fit ricker function Photo 400
#*****************************************************************************
bprior <- prior(normal(13,2), nlpar=a, lb=5)+
  prior(normal(-5,5), nlpar=awc)+
  prior(normal(-5,5), nlpar=acw)+
  prior(normal(-5,5), nlpar=aww)+
  prior(normal(45,3), nlpar=b, lb=20)+
  prior(normal(0,5), nlpar=bwc)+
  prior(normal(0,5), nlpar=bcw)+
  prior(normal(0,5), nlpar=bww)+
  prior(normal(-200,5), nlpar=c)+
  prior(normal(-5,5), nlpar=cwc)+
  prior(normal(-5,5), nlpar=ccw)+
  prior(normal(-5,5), nlpar=cww)+
  prior(normal(-259,5), nlpar=d, ub=-100)

f <- bf(P400~(a+acw*cw+awc*wc+aww*ww)*(Tk+d)*
          exp(-(Tk+d)/(b+bcw*cw+bwc*wc+bww*ww)) + 
          (c+ccw*cw+cwc*wc+cww*ww),
        a+awc+acw+aww+b+bcw+bwc+bww+c+ccw+cwc+cww+d~1, 
        nl=TRUE)
make_stancode(f, prior=bprior, family=gaussian(),
              data=dat400)
ricker_p400 <- brm(f,
                   data = dat400, 
                   prior = bprior, 
                   # sample_prior = 'only'
                   algorithm = 'sampling',
                   control = list(adapt_delta=0.995,
                                  max_treedepth=15),
                   chains = 4,
                   iter = 2000
)
summary(ricker_p400)
summary(ricker_p400, prob=0.8)$fixed
plot(ricker_p400,ask = F)
bayes_R2(ricker_p400)
sm_p400 <- broom.mixed::tidy(ricker_p400, conf.level=0.8, conf.method="HPDinterval")

topt_p400_cc <- brms::posterior_summary(fixef(ricker_p400, summary = F)[,'b_Intercept']-
                                          fixef(ricker_p400, summary = F)[,'d_Intercept']-
                                          273.15,
                                        prob=c(0.1,0.9))
topt_p400_cw <- brms::posterior_summary(fixef(ricker_p400, summary = F)[,'b_Intercept']-
                                          fixef(ricker_p400, summary = F)[,'d_Intercept']-
                                          273.15+
                                          fixef(ricker_p400, summary = F)[,'bcw_Intercept'], 
                                        prob=c(0.1,0.9))
topt_p400_wc <- brms::posterior_summary(fixef(ricker_p400, summary = F)[,'b_Intercept']-
                                          fixef(ricker_p400, summary = F)[,'d_Intercept']-
                                          273.15+
                                          fixef(ricker_p400, summary = F)[,'bwc_Intercept'], 
                                        prob=c(0.1,0.9))
topt_p400_ww <- brms::posterior_summary(fixef(ricker_p400, summary = F)[,'b_Intercept']-
                                          fixef(ricker_p400, summary = F)[,'d_Intercept']-
                                          273.15+
                                          fixef(ricker_p400, summary = F)[,'bww_Intercept'], 
                                        prob=c(0.1,0.9))
brms::loo(ricker_p400)

#*****************************************************************************
# Fit ricker function Photo 800 
#*****************************************************************************
bprior <- prior(normal(13,2), nlpar=a, lb=5)+
  prior(normal(-5,5), nlpar=awc)+
  prior(normal(-5,5), nlpar=acw)+
  prior(normal(-5,5), nlpar=aww)+
  prior(normal(45,3), nlpar=b, lb=20)+
  prior(normal(0,5), nlpar=bwc)+
  prior(normal(0,5), nlpar=bcw)+
  prior(normal(0,5), nlpar=bww)+
  prior(normal(-200,5), nlpar=c)+
  prior(normal(-5,5), nlpar=cwc)+
  prior(normal(-5,5), nlpar=ccw)+
  prior(normal(-5,5), nlpar=cww)+
  prior(normal(-259,5), nlpar=d, ub=-100)

f <- bf(P800~(a+acw*cw+awc*wc+aww*ww)*(Tk+d)*
          exp(-(Tk+d)/(b+bcw*cw+bwc*wc+bww*ww)) + 
          (c+ccw*cw+cwc*wc+cww*ww), 
        a+awc+acw+aww+b+bcw+bwc+bww+c+ccw+cwc+cww+d~1, 
        nl=TRUE)
make_stancode(f, prior=bprior, family=gaussian(),
              data=dat800)
ricker_p800 <- brm(f,
                   data = dat800, 
                   prior = bprior, 
                   algorithm = 'sampling',
                   control = list(adapt_delta=0.99,
                                  max_treedepth=13),
                   # chains = 3,
                   # iter = 250
)
summary(ricker_p800,prob=c(0.8))$fixed
plot(ricker_p800, ask=F)
bayes_R2(ricker_p800)
prior_summary(ricker_p800)

sm_p800 <- broom.mixed::tidy(ricker_p800, conf.level=0.8, conf.method="HPDinterval")

topt_p800_cc <- brms::posterior_summary(fixef(ricker_p800, summary = F)[,'b_Intercept']-
                                          fixef(ricker_p800, summary = F)[,'d_Intercept']-
                                          273.15,
                                        prob=c(0.1,0.9))
topt_p800_cw <- brms::posterior_summary(fixef(ricker_p800, summary = F)[,'b_Intercept']-
                                          fixef(ricker_p800, summary = F)[,'d_Intercept']-
                                          273.15+
                                          fixef(ricker_p800, summary = F)[,'bcw_Intercept'], 
                                        prob=c(0.1,0.9))
topt_p800_wc <- brms::posterior_summary(fixef(ricker_p800, summary = F)[,'b_Intercept']-
                                          fixef(ricker_p800, summary = F)[,'d_Intercept']-
                                          273.15+
                                          fixef(ricker_p800, summary = F)[,'bwc_Intercept'], 
                                        prob=c(0.1,0.9))
topt_p800_ww <- brms::posterior_summary(fixef(ricker_p800, summary = F)[,'b_Intercept']-
                                          fixef(ricker_p800, summary = F)[,'d_Intercept']-
                                          273.15+
                                          fixef(ricker_p800, summary = F)[,'bww_Intercept'], 
                                        prob=c(0.1,0.9))
brms::loo(ricker_p800)


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
p1 <-
  ricker_p400 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,a,awc,acw,aww,b,bcw,bwc,bww,c,ccw,cwc,cww,d),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=(a+acw*cw+awc*wc+aww*ww)*(Tk+d)*
           exp(-(Tk+d)/(b+bcw*cw+bwc*wc+bww*ww)) + 
           (c+ccw*cw+cwc*wc+cww*ww)) %>% 
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
             aes(Tc, P400,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_p400_cc[1,"Estimate"]), 
             color='blue')+
  geom_vline(aes(xintercept=topt_p400_cc[1,"Q10"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p400_cc[1,"Q90"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p400_ww[1,"Estimate"]), 
             color='red')+
  geom_vline(aes(xintercept=topt_p400_ww[1,"Q10"]), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_p400_ww[1,"Q90"]), 
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
p2 <-
  ricker_p400 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,a,awc,acw,aww,b,bcw,bwc,bww,c,ccw,cwc,cww,d),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=(a+acw*cw+awc*wc+aww*ww)*(Tk+d)*
           exp(-(Tk+d)/(b+bcw*cw+bwc*wc+bww*ww)) + 
           (c+ccw*cw+cwc*wc+cww*ww)) %>% 
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
               filter(Treat%in%c("c.c","c.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P400,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_p400_cc[1,"Estimate"]), 
             color='blue')+
  geom_vline(aes(xintercept=topt_p400_cc[1,"Q10"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p400_cc[1,"Q90"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p400_cw[1,"Estimate"]), 
             color='red')+
  geom_vline(aes(xintercept=topt_p400_cw[1,"Q10"]), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_p400_cw[1,"Q90"]), 
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
p3 <-   ricker_p400 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,a,awc,acw,aww,b,bcw,bwc,bww,c,ccw,cwc,cww,d),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=(a+acw*cw+awc*wc+aww*ww)*(Tk+d)*
           exp(-(Tk+d)/(b+bcw*cw+bwc*wc+bww*ww)) + 
           (c+ccw*cw+cwc*wc+cww*ww)) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("w.c","w.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Photo, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat400 %>% 
               filter(Treat%in%c("w.c","w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P400,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_p400_wc[1,"Estimate"]), 
             color='blue')+
  geom_vline(aes(xintercept=topt_p400_wc[1,"Q10"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p400_wc[1,"Q90"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p400_ww[1,"Estimate"]), 
             color='red')+
  geom_vline(aes(xintercept=topt_p400_ww[1,"Q10"]), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_p400_ww[1,"Q90"]), 
             color='red',lty=3)+
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
       filename = paste0("figures/photo_Topt/p400_rickerFit_CI80_",Sys.Date(),".png"),
       width = 120, height=3*100, units='mm')


#************************************************************************
# FOR THE TOP ROW of P800 --------------------------------------
#************************************************************************
p4 <-
  ricker_p800 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,a,awc,acw,aww,b,bcw,bwc,bww,c,ccw,cwc,cww,d),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=(a+acw*cw+awc*wc+aww*ww)*(Tk+d)*
           exp(-(Tk+d)/(b+bcw*cw+bwc*wc+bww*ww)) + 
           (c+ccw*cw+cwc*wc+cww*ww)) %>% 
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
             aes(Tc, P800,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_p800_cc[1,"Estimate"]), 
             color='blue')+
  geom_vline(aes(xintercept=topt_p800_cc[1,"Q10"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p800_cc[1,"Q90"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p800_ww[1,"Estimate"]), 
             color='red')+
  geom_vline(aes(xintercept=topt_p800_ww[1,"Q10"]), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_p800_ww[1,"Q90"]), 
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
p5 <-
  ricker_p800 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,a,awc,acw,aww,b,bcw,bwc,bww,c,ccw,cwc,cww,d),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=(a+acw*cw+awc*wc+aww*ww)*(Tk+d)*
           exp(-(Tk+d)/(b+bcw*cw+bwc*wc+bww*ww)) + 
           (c+ccw*cw+cwc*wc+cww*ww)) %>% 
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
               filter(Treat%in%c("c.c","c.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P800,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_p800_cc[1,"Estimate"]), 
             color='blue')+
  geom_vline(aes(xintercept=topt_p800_cc[1,"Q10"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p800_cc[1,"Q90"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p800_cw[1,"Estimate"]), 
             color='red')+
  geom_vline(aes(xintercept=topt_p800_cw[1,"Q10"]), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_p800_cw[1,"Q90"]), 
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
p6 <- ricker_p800 %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__, 
         b=b_b_Intercept) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,a,awc,acw,aww,b,bcw,bwc,bww,c,ccw,cwc,cww,d),
         Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Photo=(a+acw*cw+awc*wc+aww*ww)*(Tk+d)*
           exp(-(Tk+d)/(b+bcw*cw+bwc*wc+bww*ww)) + 
           (c+ccw*cw+cwc*wc+cww*ww)) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("w.c","w.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Photo, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat800 %>% 
               filter(Treat%in%c("w.c","w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, P800,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_p800_wc[1,"Estimate"]), 
             color='blue')+
  geom_vline(aes(xintercept=topt_p800_wc[1,"Q10"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p800_wc[1,"Q90"]), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_p800_ww[1,"Estimate"]), 
             color='red')+
  geom_vline(aes(xintercept=topt_p800_ww[1,"Q10"]), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_p800_ww[1,"Q90"]), 
             color='red',lty=3)+
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
  # 
  # ricker_p800 %>% 
  # as.data.frame() %>% 
  # sample_n(100) %>% 
  # as_tibble() %>% 
  # rename(lp=lp__, 
  #        b=b_b_Intercept) %>%
  # rename_all(fn) %>% 
  # expand(nesting(lp,kopt,kcw,kwc,kww,Topt,Tcw,Twc,Tww,b),
  #        Tk=seq(295,318,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  # mutate(Photo=kopt+kcw*cw+kwc*wc+kww*ww-b*(Tk-(Topt+Tcw*cw+Twc*wc+Tww*ww))**2) %>% 
# mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
#                        wc==1 & cw==0 & ww==0 ~"w.c",
#                        ww==1 & wc==0 & cw==0 ~"w.w",
#                        ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
# filter(is.na(Treat)==F) %>% 
# mutate(lp_Treat=paste(lp,Treat)) %>% 
# filter(Treat %in% c("w.w","w.c")) %>%
# mutate(Tc=Tk-273.15) %>% 
# ggplot(data=.,aes(Tc, Photo, color=Treat, group=lp_Treat))+
# geom_line(alpha=0.05)+
# geom_point(data=dat400 %>% 
#              filter(Treat%in%c("w.w")) %>% 
#              mutate(Tc=Tk-273.15), 
#            aes(Tc, P800,color=Treat),
#            size=5,alpha=0.5,
#            inherit.aes = F, pch=20)+
# geom_point(data=dat400 %>% 
#              filter(Treat%in%c("w.c")) %>% 
#              mutate(Tc=Tk-273.15), 
#            aes(Tc, P800,color=Treat),
#            size=4,alpha=0.5,
#            inherit.aes = F, pch=1)+
# geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
#                  filter(sm_p800,term=="Tww_(Intercept)")$estimate-273.15), 
#            color='red')+
# geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
#                  filter(sm_p800,term=="Tww_(Intercept)")$conf.high-273.15), 
#            color='red',lty=3)+
# geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
#                  filter(sm_p800,term=="Tww_(Intercept)")$conf.low-273.15), 
#            color='red',lty=3)+
# geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
#                  filter(sm_p800,term=="Twc_(Intercept)")$estimate-273.15), 
#            color='blue')+
# geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
#                  filter(sm_p800,term=="Twc_(Intercept)")$conf.high-273.15), 
#            color='blue',lty=3)+
# geom_vline(aes(xintercept=filter(sm_p800,term=="Topt_(Intercept)")$estimate+
#                  filter(sm_p800,term=="Twc_(Intercept)")$conf.low-273.15), 
#            color='blue',lty=3)+
# scale_color_manual("", 
#                    values=c("red","blue"), 
#                    breaks=c("w.w","w.c"),
#                    labels=c("Treatment",
#                             sprintf('Treatment\u2192Control')
#                    ), 
#                    guide=guide_legend(override.aes = list(shape=c(20,1), 
#                                                           size=c(6,3)))
#                    # parse("paste(Control %->% \" Treatment\")"))
# )+
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
       filename = paste0("figures/photo_Topt/p800_rickerFit_CI80_",Sys.Date(),".png"),
       width = 120, height=3*100, units='mm')




ggsave((p1/p2/p3)|(p4/p5/p6),
       filename = paste0("figures/photo_Topt/p400_p800_rickerFit_CI80_",Sys.Date(),".png"),
       width = 2*120, height=3*100, units='mm')
ggsave((p1/p2/p3)|(p4/p5/p6),
       filename = paste0("figures/photo_Topt/p400_p800_rickerFit_CI80_",Sys.Date(),".pdf"),
       width = 2*120, height=3*100, units='mm')






































