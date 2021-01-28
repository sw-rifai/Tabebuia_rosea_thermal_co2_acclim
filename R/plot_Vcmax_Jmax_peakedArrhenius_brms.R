library(tidyverse); library(bayesplot);
library(rstan); library(brms); library(nls.multstart)
library(patchwork)
options(mc.cores = parallel::detectCores()-1)
set.seed(1); 

dat_all<-read_csv("data/params.all.28032018.csv")
# dat_fci<-read_csv("data/params_20190124.csv")
dat_all <- dat_all %>% 
  mutate(Tk = Tleaf+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0))

# Vcmax ~ kopt*(Hd*exp( (Ha*(Tk-Topt))/(Tk*R*Topt) ))

# Fit Vcmax ---------------------------------------------------------------
bprior_v <- 
  prior(normal(250 , 100), nlpar = kopt, lb=25, ub=500)+
  prior(normal( 0,100), nlpar = kwc,  lb=-250, ub=250)+
  prior(normal( 0,100), nlpar = kcw, lb=-250, ub=250)+
  prior(normal(0,100), nlpar = kww, lb=-250, ub=250)+
  prior(normal(310,10), nlpar = Topt,  lb=290, ub=335)+
  prior(normal( 0,2 ), nlpar = Tcw,  lb=-10, ub=10)+
  prior(normal( 0,2 ), nlpar = Twc, lb=-10, ub=10)+
  prior(normal( 0,2 ), nlpar = Tww, lb=-10, ub=10)+
  prior(normal( 75,20 ), nlpar = Ha,  lb=25, ub=150)+
  prior(normal( 0,20 ), nlpar = Hawc,  lb=-25, ub=25)+
  prior(normal( 0,20 ), nlpar = Hacw, lb=-25, ub=25)+
  prior(normal( 0,20 ), nlpar = Haww, lb=-25, ub=25)

f_v <- bf(Vcmax ~ 
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

make_stancode(f_v, prior=bprior_v, family=gaussian(),
              data=dat_all)
parr_vcmax <- brm(f_v,
           data = dat_all, 
           prior = bprior_v, 
           # sample_prior = 'only')
           algorithm = 'sampling',
           control = list(adapt_delta=0.999)
)
write_rds(parr_vcmax, file = paste0("outputs/parr_vcmax_",Sys.Date(),".rds"))
parr_vcmax
plot(parr_vcmax)
pp_check(parr_vcmax, nsamples = 30)


sm_vcmax <- broom.mixed::tidy(parr_vcmax, conf.level=0.8, conf.method="HPDinterval")

topt_vcmax_cc <- brms::posterior_summary(fixef(parr_vcmax, summary = F)[,'Topt_Intercept'],
                                        prob=c(0.1,0.9))
topt_vcmax_wc <- brms::posterior_summary(fixef(parr_vcmax, summary = F)[,'Twc_Intercept'],
                                         prob=c(0.1,0.9))
topt_vcmax_cw <- brms::posterior_summary(fixef(parr_vcmax, summary = F)[,'Tcw_Intercept'],
                                         prob=c(0.1,0.9))
topt_vcmax_ww <- brms::posterior_summary(fixef(parr_vcmax, summary = F)[,'Tww_Intercept'],
                                         prob=c(0.1,0.9))

# Vcmax: Calculate delta S entropy estimates per treatment -------------------------------------
# Control 
fixef(parr_vcmax, summary = F) %>% as_tibble() %>%
  mutate(entropy = 1000*200/(Topt_Intercept) + 8.3145*log(1000*Ha_Intercept/(1000*200-1000*Ha_Intercept))) %>%
  pull(entropy) %>% quantile(., c(0.1,0.5,0.9))

# Control -> Treatment
fixef(parr_vcmax, summary = F) %>% as_tibble() %>%
  mutate(entropy = 1000*200/(Topt_Intercept+Tcw_Intercept) + 
           8.3145*log((1000*Ha_Intercept+Hacw_Intercept)/(1000*200-1000*(Ha_Intercept+Hacw_Intercept)))) %>%
  pull(entropy) %>% quantile(., c(0.1,0.5,0.9))

# Treatment 
fixef(parr_vcmax, summary = F) %>% as_tibble() %>%
  mutate(entropy = 1000*200/(Topt_Intercept+Tww_Intercept) + 
           8.3145*log((1000*Ha_Intercept+Haww_Intercept)/(1000*200-1000*(Ha_Intercept+Haww_Intercept)))) %>%
  pull(entropy) %>% quantile(., c(0.1,0.5,0.9))

# Treatment 
fixef(parr_vcmax, summary = F) %>% as_tibble() %>%
  mutate(entropy = 1000*200/(Topt_Intercept+Twc_Intercept) + 
           8.3145*log((1000*Ha_Intercept+Hawc_Intercept)/(1000*200-1000*(Ha_Intercept+Hawc_Intercept)))) %>%
  pull(entropy) %>% quantile(., c(0.1,0.5,0.9))


# Fit Jmax ---------------------------------------------------------------
bprior_j <- 
  prior(normal(250 , 100), nlpar = kopt, lb=25, ub=500)+
  prior(normal( 0,100), nlpar = kwc,  lb=-250, ub=250)+
  prior(normal( 0,100), nlpar = kcw, lb=-250, ub=250)+
  prior(normal(0,100), nlpar = kww, lb=-250, ub=250)+
  prior(normal(310,10), nlpar = Topt,  lb=290, ub=335)+
  prior(normal( 0,2 ), nlpar = Tcw,  lb=-10, ub=10)+
  prior(normal( 0,2 ), nlpar = Twc, lb=-10, ub=10)+
  prior(normal( 0,2 ), nlpar = Tww, lb=-10, ub=10)+
  prior(normal( 75,20 ), nlpar = Ha,  lb=25, ub=150)+
  prior(normal( 0,20 ), nlpar = Hawc,  lb=-25, ub=25)+
  prior(normal( 0,20 ), nlpar = Hacw, lb=-25, ub=25)+
  prior(normal( 0,20 ), nlpar = Haww, lb=-25, ub=25)

f_j <- bf(Jmax ~ 
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

make_stancode(f_j, prior=bprior_j, family=gaussian(),
              data=dat_all)
parr_jmax <- brm(f_j,
                  data = dat_all, 
                  prior = bprior_j, 
                  # sample_prior = 'only')
                  algorithm = 'sampling',
                  control = list(adapt_delta=0.999)
)
write_rds(parr_jmax, file = paste0("outputs/parr_jmax_",Sys.Date(),".rds"))

parr_jmax
plot(parr_jmax)
pp_check(parr_jmax, nsamples = 30)

sm_jmax <- broom.mixed::tidy(parr_jmax, conf.level=0.8, conf.method="HPDinterval")

topt_jmax_cc <- brms::posterior_summary(fixef(parr_jmax, summary = F)[,'Topt_Intercept'],
                                         prob=c(0.1,0.9))
topt_jmax_wc <- brms::posterior_summary(fixef(parr_jmax, summary = F)[,'Twc_Intercept'],
                                         prob=c(0.1,0.9))
topt_jmax_cw <- brms::posterior_summary(fixef(parr_jmax, summary = F)[,'Tcw_Intercept'],
                                         prob=c(0.1,0.9))
topt_jmax_ww <- brms::posterior_summary(fixef(parr_jmax, summary = F)[,'Tww_Intercept'],
                                         prob=c(0.1,0.9))

# Jmax: Calculate delta S entropy estimates per treatment -------------------------------------
# Control 
fixef(parr_jmax, summary = F) %>% as_tibble() %>%
  mutate(entropy = 1000*200/(Topt_Intercept) + 8.3145*log(1000*Ha_Intercept/(1000*200-1000*Ha_Intercept))) %>%
  pull(entropy) %>% quantile(., c(0.1,0.5,0.9))

# Control -> Treatment
fixef(parr_jmax, summary = F) %>% as_tibble() %>%
  mutate(entropy = 1000*200/(Topt_Intercept+Tcw_Intercept) + 
           8.3145*log((1000*Ha_Intercept+Hacw_Intercept)/(1000*200-1000*(Ha_Intercept+Hacw_Intercept)))) %>%
  pull(entropy) %>% quantile(., c(0.1,0.5,0.9))

# Treatment 
fixef(parr_jmax, summary = F) %>% as_tibble() %>%
  mutate(entropy = 1000*200/(Topt_Intercept+Tww_Intercept) + 
           8.3145*log((1000*Ha_Intercept+Haww_Intercept)/(1000*200-1000*(Ha_Intercept+Haww_Intercept)))) %>%
  pull(entropy) %>% quantile(., c(0.1,0.5,0.9))

# Treatment 
fixef(parr_jmax, summary = F) %>% as_tibble() %>%
  mutate(entropy = 1000*200/(Topt_Intercept+Twc_Intercept) + 
           8.3145*log((1000*Ha_Intercept+Hawc_Intercept)/(1000*200-1000*(Ha_Intercept+Hawc_Intercept)))) %>%
  pull(entropy) %>% quantile(., c(0.1,0.5,0.9))


# Plotting helpers --------------------------------------------------------
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



# Plot Vcmax --------------------------------------------------------------
#************************************************************************
# FOR THE TOP ROW of Vcmax -------------------
#************************************************************************
p1 <-
  parr_vcmax %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kwc,kcw,kww, 
                 Topt,Twc,Tcw,Tww, 
                 Ha,Hacw,Hawc,Haww),
         Tk=seq(295,320,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Vcmax = 
           (kopt + kcw*cw + kwc*wc + kww*ww) * 
           ((200 * (2.718282^(((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                                 (Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))) / 
              (200 - ((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                        (1-(2.718282^((200*(Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                        (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))))))
         ) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","w.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Vcmax, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("c.c","w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Vcmax,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=topt_vcmax_cc[1,"Q10"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_vcmax_cc[1,"Q90"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_vcmax_ww[1,"Estimate"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='red')+
  geom_vline(aes(xintercept=topt_vcmax_ww[1,"Q10"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_vcmax_ww[1,"Q90"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  scale_color_manual("",
                     values=c("c.c"='blue', "w.w"='#cf0000'), 
                     labels=c("c.c"="Control",
                              "c.w"="paste(Control %->% \" Treatment\")",
                              "w.c"="paste(Treatment %->% \" Control\")", 
                              "w.w"="Treatment"))+
  theme_bw()+
  labs(y=bquote(paste('V'['cmax'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,400)+
  scale_x_continuous(limits=c(25,46.5), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = c(0.05,0.95), 
        legend.justification = c(0.05,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p1

#************************************************************************
# FOR THE MIDDLE ROW of Vcmax -------------------
#************************************************************************
p2 <-
  parr_vcmax %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kwc,kcw,kww, 
                 Topt,Twc,Tcw,Tww, 
                 Ha,Hacw,Hawc,Haww),
         Tk=seq(295,320,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Vcmax = 
           (kopt + kcw*cw + kwc*wc + kww*ww) * 
           ((200 * (2.718282^(((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                                 (Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))) / 
              (200 - ((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                        (1-(2.718282^((200*(Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                        (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))))))
  ) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","c.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Vcmax, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("c.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Vcmax,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("c.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Vcmax,color=Treat),
             size=4,alpha=0.5,pch=1,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=topt_vcmax_cc[1,"Q10"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_vcmax_cc[1,"Q90"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_vcmax_cw[1,"Estimate"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='red')+
  geom_vline(aes(xintercept=topt_vcmax_cw[1,"Q10"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_vcmax_cw[1,"Q90"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
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
  labs(y=bquote(paste('V'['cmax'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,400)+
  scale_x_continuous(limits=c(25,46.5), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = c(0.05,0.95), 
        legend.justification = c(0.05,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p2

#************************************************************************
# FOR THE Bottom ROW of Vcmax -------------------
#************************************************************************
p3 <-
  parr_vcmax %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kwc,kcw,kww, 
                 Topt,Twc,Tcw,Tww, 
                 Ha,Hacw,Hawc,Haww),
         Tk=seq(295,320,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Vcmax = 
           (kopt + kcw*cw + kwc*wc + kww*ww) * 
           ((200 * (2.718282^(((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                                 (Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))) / 
              (200 - ((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                        (1-(2.718282^((200*(Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                        (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))))))
  ) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("w.c","w.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Vcmax, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Vcmax,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("w.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Vcmax,color=Treat),
             size=4,alpha=0.5,pch=1,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_vcmax_wc[1,"Estimate"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=topt_vcmax_wc[1,"Q10"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_vcmax_wc[1,"Q90"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_vcmax_ww[1,"Estimate"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='red')+
  geom_vline(aes(xintercept=topt_vcmax_ww[1,"Q10"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_vcmax_ww[1,"Q90"]+topt_vcmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  scale_color_manual("", 
                     values=c("red","blue"), 
                     breaks=c("w.w","w.c"),
                     labels=c("Treatment",
                              sprintf('Treatment\u2192Control')
                     ),
                     guide=guide_legend(override.aes = list(shape=c(20,1), 
                                                            size=c(6,3)))
  )+
  theme_bw()+
  labs(y=bquote(paste('V'['cmax'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,400)+
  scale_x_continuous(limits=c(25,46.5), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = c(0.05,0.95), 
        legend.justification = c(0.05,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p3


# Plot Jcmax --------------------------------------------------------------
#************************************************************************
# FOR THE TOP ROW of Jcmax -------------------
#************************************************************************
p4 <-
  parr_jmax %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kwc,kcw,kww, 
                 Topt,Twc,Tcw,Tww, 
                 Ha,Hacw,Hawc,Haww),
         Tk=seq(295,320,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Jmax = 
           (kopt + kcw*cw + kwc*wc + kww*ww) * 
           ((200 * (2.718282^(((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                                 (Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))) / 
              (200 - ((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                        (1-(2.718282^((200*(Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                        (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))))))
  ) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","w.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Jmax, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("c.c","w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Jmax,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_jmax_cc[1,"Estimate"]-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=topt_jmax_cc[1,"Q10"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_jmax_cc[1,"Q90"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_jmax_ww[1,"Estimate"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='red')+
  geom_vline(aes(xintercept=topt_jmax_ww[1,"Q10"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_jmax_ww[1,"Q90"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  scale_color_manual("",
                     values=c("c.c"='blue', "w.w"='#cf0000'),
                     labels=c("c.c"="Control",
                              "c.w"="paste(Control %->% \" Treatment\")",
                              "w.c"="paste(Treatment %->% \" Control\")",
                              "w.w"="Treatment"))+
  theme_bw()+
  labs(y=bquote(paste('J'['max'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(position='right', limits=c(0,400))+
  scale_x_continuous(limits=c(25,46.5), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = 'none', #c(0.05,0.95), 
        legend.justification = c(0.05,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title.y = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p4

#************************************************************************
# FOR THE MIDDLE ROW of Jmax -------------------
#************************************************************************
p5 <-
  parr_jmax %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kwc,kcw,kww, 
                 Topt,Twc,Tcw,Tww, 
                 Ha,Hacw,Hawc,Haww),
         Tk=seq(295,320,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Jmax = 
           (kopt + kcw*cw + kwc*wc + kww*ww) * 
           ((200 * (2.718282^(((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                                 (Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))) / 
              (200 - ((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                        (1-(2.718282^((200*(Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                        (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))))))
  ) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","c.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Jmax, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("c.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Jmax,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("c.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Jmax,color=Treat),
             size=4,alpha=0.5,pch=1,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_jmax_cc[1,"Estimate"]-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=topt_jmax_cc[1,"Q10"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_jmax_cc[1,"Q90"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_jmax_cw[1,"Estimate"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='red')+
  geom_vline(aes(xintercept=topt_jmax_cw[1,"Q10"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_jmax_cw[1,"Q90"]+topt_jmax_cc[1,"Estimate"]-273.15), 
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
  labs(y=bquote(paste('J'['max'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(limits=c(0,400), position='right')+
  scale_x_continuous(limits=c(25,46.5), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = 'none',#c(0.05,0.95), 
        legend.justification = c(0.05,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title.y = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p5

#************************************************************************
# FOR THE Bottom ROW of Jmax -------------------
#************************************************************************
p6 <-
  parr_jmax %>% 
  as.data.frame() %>% 
  sample_n(100) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>% 
  expand(nesting(lp,kopt,kwc,kcw,kww, 
                 Topt,Twc,Tcw,Tww, 
                 Ha,Hacw,Hawc,Haww),
         Tk=seq(295,320,0.1), cw=c(0,1), wc=c(0,1), ww=c(0,1)) %>% 
  mutate(Jmax = 
           (kopt + kcw*cw + kwc*wc + kww*ww) * 
           ((200 * (2.718282^(((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                                 (Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))) / 
              (200 - ((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                        (1-(2.718282^((200*(Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                        (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))))))
  ) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("w.c","w.w")) %>%
  mutate(Tc=Tk-273.15) %>% 
  ggplot(data=.,aes(Tc, Jmax, color=Treat, group=lp_Treat))+
  geom_line(alpha=0.05)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("w.w")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Jmax,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat_all %>% 
               filter(Treat%in%c("w.c")) %>% 
               mutate(Tc=Tk-273.15), 
             aes(Tc, Jmax,color=Treat),
             size=4,alpha=0.5,pch=1,
             inherit.aes = F)+
  geom_vline(aes(xintercept=topt_jmax_wc[1,"Estimate"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='blue')+
  geom_vline(aes(xintercept=topt_jmax_wc[1,"Q10"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_jmax_wc[1,"Q90"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='blue',lty=3)+
  geom_vline(aes(xintercept=topt_jmax_ww[1,"Estimate"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='red')+
  geom_vline(aes(xintercept=topt_jmax_ww[1,"Q10"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  geom_vline(aes(xintercept=topt_jmax_ww[1,"Q90"]+topt_jmax_cc[1,"Estimate"]-273.15), 
             color='red',lty=3)+
  scale_color_manual("", 
                     values=c("red","blue"), 
                     breaks=c("w.w","w.c"),
                     labels=c("Treatment",
                              sprintf('Treatment\u2192Control')
                     ),
                     guide=guide_legend(override.aes = list(shape=c(20,1), 
                                                            size=c(6,3)))
  )+
  theme_bw()+
  labs(y=bquote(paste('J'['max'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(position='right',limits=c(0,400))+
  scale_x_continuous(limits=c(25,46.5), breaks = seq(25,45,by=5), expand = expansion(0,0.1))+
  theme(legend.position = 'none',#c(0.05,0.95), 
        legend.justification = c(0.05,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title.y = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL); p6



# (p1/p2/p3)|(p4/p5/p6)
ggsave((p1/p2/p3)|(p4/p5/p6),
       filename = paste0("figures/vcmax_jmax_brms_",Sys.Date(),".png"),
       width = 2*120, height=3*100, units='mm')
ggsave((p1/p2/p3)|(p4/p5/p6),
       filename = paste0("figures/vcmax_jmax_brms_",Sys.Date(),".pdf"),
       width = 2*120, height=3*100, units='mm')
