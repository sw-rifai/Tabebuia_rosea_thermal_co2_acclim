

library(tidyverse)
library(ggpubr)
theme_set(theme_pubr())

P400<-read.csv("P400.150321.csv")
P800<-read.csv("P800.150321.csv")

x.cc<-P400$Tleaf[P400$Treat=="c.c"]
y.cc<-P400$Cond[P400$Treat=="c.c"]

x.ww<-P800$Tleaf[P800$Treat=="w.w"]
y.ww<-P800$Cond[P800$Treat=="w.w"]

x.cw<-P800$Tleaf[P800$Treat=="c.w"]
y.cw<-P800$Cond[P800$Treat=="c.w"]

x.wc<-P400$Tleaf[P400$Treat=="w.c"]
y.wc<-P400$Cond[P400$Treat=="w.c"]

(p1 <- tibble(x.cc,y.cc) %>%
  ggplot(data=.,aes(x.cc,y.cc))+
  geom_point()+ylim(0,0.8)+xlim(25,45)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihoo
                               select=TRUE),  # penalizes the wiggliness of the spline
              color='blue',fill='dodgerblue')+theme_classic())




(p2 <- tibble(x.ww,y.ww) %>%
  ggplot(data=.,aes(x.ww,y.ww))+
  geom_point()+ylim(0,0.8)+xlim(25,45)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihoo
                               select=TRUE),  # penalizes the wiggliness of the spline
              color='red',fill='red')+theme_classic())


(p3 <- tibble(x.cw,y.cw) %>%
  ggplot(data=.,aes(x.cw,y.cw))+
  geom_point()+ylim(0,0.8)+xlim(25,45)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihoo
                               select=TRUE),  # penalizes the wiggliness of the spline
              color='pink',fill='pink')+theme_classic())



(p4 <- tibble(x.wc,y.wc) %>%
  ggplot(data=.,aes(x.wc,y.wc))+
  geom_point()+ylim(0,0.8)+xlim(25,45)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihoo
                               select=TRUE),  # penalizes the wiggliness of the spline
              color='lightblue',fill='lightblue')+theme_classic())




# Now same for gs vs VPD
  

# top panel; treatment at 800 ppm, control at 400 ppm 

x.cc.v<-P400$VpdL[P400$Treat=="c.c"]
y.cc.g<-P400$Cond[P400$Treat=="c.c"]

x.ww.v<-P800$VpdL[P800$Treat=="w.w"]
y.ww.g<-P800$Cond[P800$Treat=="w.w"]

x.wc.v<-P400$VpdL[P400$Treat=="w.c"]
y.wc.g<-P400$Cond[P400$Treat=="w.c"]

x.cw.v<-P800$VpdL[P800$Treat=="c.w"]
y.cw.g<-P800$Cond[P800$Treat=="c.w"]



(p5<-tibble(x.cc.v,y.cc.g) %>%
  ggplot(data=.,aes(x.cc.v,y.cc.g))+
  geom_point()+ylim(0,0.8)+xlim(0,4)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihoo
                               select=TRUE),  # penalizes the wiggliness of the spline
              color='blue',fill='dodgerblue')+theme_classic())


(p6<-tibble(x.ww.v,y.ww.g) %>%
  ggplot(data=.,aes(x.ww.v,y.ww.g))+
  geom_point()+ylim(0,0.8)+xlim(0,4)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihoo
                               select=TRUE),  # penalizes the wiggliness of the spline
              color='red',fill='red')+theme_classic())


(p7<-tibble(x.cw.v,y.cw.g) %>%
  ggplot(data=.,aes(x.cw.v,y.cw.g))+
  geom_point()+ylim(0,0.8)+xlim(0,4)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihoo
                               select=TRUE),  # penalizes the wiggliness of the spline
              color='pink',fill='pink')+theme_classic())


(p8<-tibble(x.wc.v,y.wc.g) %>%
  ggplot(data=.,aes(x.wc.v,y.wc.g))+
  geom_point()+ylim(0,0.8)+xlim(0,4)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihoo
                               select=TRUE),  # penalizes the wiggliness of the spline
              color='lightblue',fill='lightblue')+theme_classic())



figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8,
                    labels = c("(a)", "(b)", "(c)","(d)","(e)", "(f)", "(g)","(h)"),
                    ncol = 2, nrow = 4)
figure




