library(tidyverse)
library(ggpubr)
theme_set(theme_pubr())

P400<-read.csv("data/P400.150321.csv")
P800<-read.csv("data/P800.150321.csv")

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



################################################################################
# REVAMP ------------------------------------
library(tidyverse)
library(ggpubr)

p400<-read_csv("data/P400.150321.csv") %>% 
  mutate(Tk = Tleaf+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0)) %>% 
  mutate(ca=400)
p800<-read_csv("data/P800.150321.csv") %>% 
  mutate(Tk = Tleaf+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0)) %>% 
  mutate(ca=800)


# By Tleaf ****************************************************
p_a <- bind_rows(
  p400 %>% filter(Treat %in% c('c.c')),
  p800 %>% filter(Treat %in% c('w.w')),
   ) %>% 
  ggplot(data=.,aes(Tleaf,Cond,color=Treat,fill=Treat))+
  geom_point(size=3.5)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihood
                               select=TRUE),  # penalizes the wiggliness of the spline
                    )  +
  annotate("text", x = 39-5, y = c(0.8, 0.7), hjust=0,size=4,
           # label = c("paste(Treatment %->% \" Treatment\")", 
           #           "paste(Treatment %->% \" Control\")"),
           label = c("Control | 400 ppm", 
                     "Treatment | 800 ppm"),
           parse=F)+
  annotate("point", x = 38.5-5, y = c(0.8,0.7), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,16))+
  scale_y_continuous(limits=c(0,0.85), 
                     breaks=c(0,0.2,0.4,0.6,0.8),
                     expand=c(0,0))+
  scale_x_continuous(limits=c(24.5,45.5), 
                     expand=c(0,0),
                     breaks=c(25,30,35,40,45))+
  scale_color_manual(values=c("c.c"="blue","w.w"="red"))+
  scale_fill_manual(values=c("c.c"="blue","w.w"="red"))+
  labs(x=NULL,y=expression(paste(g[s]~(mol~m**-2~s**-1))))+
  theme_linedraw()+
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()); p_a


p_c <- bind_rows(
  p400 %>% filter(Treat %in% c('c.c')),
  p800 %>% filter(Treat %in% c('c.w')),
) %>% 
  ggplot(data=.,aes(Tleaf,Cond,color=Treat,fill=Treat))+
  geom_point(data=. %>% filter(Treat=='c.c'),pch=20,size=3.5)+
  geom_point(data=. %>% filter(Treat=='c.w'),pch=1,size=3.5)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihood
                               select=TRUE),  # penalizes the wiggliness of the spline
  )  +
  annotate("text", x = 39-5, y = c(0.8, 0.7), hjust=0,size=4,
           # label = c("paste(Treatment %->% \" Treatment\")", 
           #           "paste(Treatment %->% \" Control\")"),
           label = c("Control | 400ppm", 
                     "Control \U2192 Treatment | 800 ppm"),
           parse=F)+
  annotate("point", x = 38.5-5, y = c(0.8,0.7), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,16))+
  scale_y_continuous(limits=c(0,0.85), 
                     breaks=c(0,0.2,0.4,0.6,0.8),
                     expand=c(0,0))+
  scale_x_continuous(limits=c(24.5,45.5), 
                     expand=c(0,0),
                     breaks=c(25,30,35,40,45))+
  scale_color_manual(values=c("c.c"="blue","c.w"="red"))+
  scale_fill_manual(values=c("c.c"="blue","c.w"="red"))+
  labs(x=NULL,y=expression(paste(g[s]~(mol~m**-2~s**-1))))+
  theme_linedraw()+
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()); p_c


p_e <- bind_rows(
  p800 %>% filter(Treat %in% c('w.w')),
  p400 %>% filter(Treat %in% c('w.c')),
) %>% 
  ggplot(data=.,aes(Tleaf,Cond,color=Treat,fill=Treat))+
  geom_point(data=. %>% filter(Treat=='w.w'),pch=16,size=3.5)+
  geom_point(data=. %>% filter(Treat=='w.c'),pch=1,size=3.5)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihood
                               select=TRUE),  # penalizes the wiggliness of the spline
  )  +
  annotate("text", x = 39-5, y = c(0.8, 0.7), hjust=0,size=4,
           # label = c("paste(Treatment %->% \" Treatment\")", 
           #           "paste(Treatment %->% \" Control\")"),
           label = c("Treatment | 800ppm", 
                     "Treatment \U2192 Control | 400 ppm"),
           parse=F)+
  annotate("point", x = 38.5-5, y = c(0.8,0.7), 
           colour = c('red',"blue"), 
           size = 5, alpha=0.5, shape=c(16,16))+
  scale_y_continuous(limits=c(0,0.85), 
                     breaks=c(0,0.2,0.4,0.6,0.8),
                     expand=c(0,0))+
  scale_x_continuous(limits=c(24.5,45.5), 
                     expand=c(0,0),
                     breaks=c(25,30,35,40,45))+
  scale_color_manual(values=c("w.c"="blue","w.w"="red"))+
  scale_fill_manual(values=c("w.c"="blue","w.w"="red"))+
  labs(x=expression(paste(Leaf~temperature~(degree*C))),
       y=expression(paste(g[s]~(mol~m**-2~s**-1))))+
  theme_linedraw()+
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        panel.grid = element_blank()); p_e







# By VPD ***************************************************** 
p_b <- bind_rows(
  p400 %>% filter(Treat %in% c('c.c')),
  p800 %>% filter(Treat %in% c('w.w')),
) %>% 
  ggplot(data=.,aes(VpdL,Cond,color=Treat,fill=Treat))+
  geom_point(size=3.5)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihood
                               select=TRUE),  # penalizes the wiggliness of the spline
  )  +
  # annotate("text", x = 39, y = c(0.8, 0.7), hjust=0,size=4,
  #          # label = c("paste(Treatment %->% \" Treatment\")", 
  #          #           "paste(Treatment %->% \" Control\")"),
  #          label = c("Control | 400 ppm", 
  #                    "Treatment | 800 ppm"),
  #          parse=F)+
  # annotate("point", x = 38.5, y = c(0.8,0.7), 
  #          colour = c('blue',"red"), 
  #          size = 5, alpha=0.5, shape=c(16,16))+
  scale_y_continuous(limits=c(0,0.85), 
                     breaks=c(0,0.2,0.4,0.6,0.8),
                     expand=c(0,0))+
  scale_x_continuous(limits=c(0,4), 
                     expand=c(0,0),
                     breaks=c(1:4))+
  scale_color_manual(values=c("c.c"="blue","w.w"="red"))+
  scale_fill_manual(values=c("c.c"="blue","w.w"="red"))+
  labs(x=NULL,y=expression(paste(g[s]~(mol~m**-2~s**-1))))+
  theme_linedraw()+
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()); p_b


p_d <- bind_rows(
  p400 %>% filter(Treat %in% c('c.c')),
  p800 %>% filter(Treat %in% c('c.w')),
) %>% 
  ggplot(data=.,aes(VpdL,Cond,color=Treat,fill=Treat))+
  geom_point(data=. %>% filter(Treat=='c.c'),pch=20,size=3.5)+
  geom_point(data=. %>% filter(Treat=='c.w'),pch=1,size=3.5)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihood
                               select=TRUE),  # penalizes the wiggliness of the spline
  )  +
  # annotate("text", x = 39, y = c(0.8, 0.7), hjust=0,size=4,
  #          # label = c("paste(Treatment %->% \" Treatment\")", 
  #          #           "paste(Treatment %->% \" Control\")"),
  #          label = c("Control | 400ppm", 
  #                    "Control \U2192 Treatment | 800 ppm"),
  #          parse=F)+
  # annotate("point", x = 38.5, y = c(0.8,0.7), 
  #          colour = c('blue',"red"), 
  #          size = 5, alpha=0.5, shape=c(16,16))+
  scale_y_continuous(limits=c(0,0.85), 
                     breaks=c(0,0.2,0.4,0.6,0.8),
                     expand=c(0,0))+
  scale_x_continuous(limits=c(0,4), 
                     expand=c(0,0),
                     breaks=c(0:4))+
  scale_color_manual(values=c("c.c"="blue","c.w"="red"))+
  scale_fill_manual(values=c("c.c"="blue","c.w"="red"))+
  labs(x=NULL,y=expression(paste(g[s]~(mol~m**-2~s**-1))))+
  theme_linedraw()+
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()); p_d


p_f <- bind_rows(
  p800 %>% filter(Treat %in% c('w.w')),
  p400 %>% filter(Treat %in% c('w.c')),
) %>% 
  ggplot(data=.,aes(VpdL,Cond,color=Treat,fill=Treat))+
  geom_point(data=. %>% filter(Treat=='w.w'),pch=16,size=3.5)+
  geom_point(data=. %>% filter(Treat=='w.c'),pch=1,size=3.5)+
  geom_smooth(method='gam', # generalized additive model
              formula=y~s(x,bs='cs'), # cubic regression spline; good at not getting too wiggly
              method.args=list(method='REML', # restricted maximum likelihood
                               select=TRUE),  # penalizes the wiggliness of the spline
  )  +
  # annotate("text", x = 40-1, y = c(0.8, 0.7), hjust=0,size=4,
  #          # label = c("paste(Treatment %->% \" Treatment\")", 
  #          #           "paste(Treatment %->% \" Control\")"),
  #          label = c("Treatment | 800ppm", 
  #                    "Treatment \U2192 Control | 400 ppm"),
  #          parse=F)+
  # annotate("point", x = 39.5-1, y = c(0.8,0.7), 
  #          colour = c('red',"blue"), 
  #          size = 5, alpha=0.5, shape=c(16,16))+
  scale_y_continuous(limits=c(0,0.85), 
                     breaks=c(0,0.2,0.4,0.6,0.8),
                     expand=c(0,0))+
  scale_x_continuous(limits=c(0,4), 
                     expand=c(0,0),
                     breaks=c(0:4))+
  scale_color_manual(values=c("w.c"="blue","w.w"="red"))+
  scale_fill_manual(values=c("w.c"="blue","w.w"="red"))+
  labs(x=expression(paste(Leaf~temperature~(degree*C))),
       y=expression(paste(g[s]~(mol~m**-2~s**-1))))+
  theme_linedraw()+
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        panel.grid = element_blank()); p_f



ggarrange(p_a,p_b,p_c,p_d,p_e,p_f, 
          labels=c("(a)","(b)","(c)","(d)","(e)","(f)"),
          hjust=-3.5,
          vjust = 2,
          ncol=2, nrow=3)
ggsave(filename = 'figures/Fig5_Gs_Tleaf_VPD.png',
       width=250, 
       height=250,
       dpi=300,
       units='mm')
