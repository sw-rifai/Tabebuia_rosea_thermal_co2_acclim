library(tidyverse); 
library(tidybayes); 

dat_all<-read_csv("data/params.all.28032018.csv")
dat_all <- dat_all %>% 
  mutate(Tk = Tleaf+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0))


fn <- function(x){
  x <- gsub("_Intercept","",x)
  x <- gsub("b_","",x)
  x
}

f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}
f_photo_mp <- function(kopt,kcw,kwc,kww,
                       Ha,Hacw,Hawc,Haww,
                       Topt,Tcw,Twc,Tww,
                       cw,wc,ww,
                       Tk,...){
  photo <- (kopt + kcw*cw + kwc*wc + kww*ww) * 
   ((200 * (2.718282^(((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                          (Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                         (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))) / 
       (200 - ((Ha + Hacw*cw + Hawc*wc + Haww*ww)*
                 (1-(2.718282^((200*(Tk-(Topt + Tcw*cw + Twc*wc +Tww*ww)))/
                                 (Tk*0.008314*(Topt + Tcw*cw + Twc*wc +Tww*ww))))))))
                                  
  # photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
  #                    (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}

vc_fit <- read_rds("outputs/parr_vcmax_2021-03-15.rds")
jm_fit <- read_rds("outputs/parr_jmax_2021-03-15.rds")

# panel (a) control & treatment -------------------------------------------------
dat_pa_cc <- inner_join(
  {vc_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=1,cw=0,wc=0,ww=0) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Vcmax=Photo) %>% 
      select(Tk,Vcmax)},
  {jm_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=1,cw=0,wc=0,ww=0) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Jmax=Photo) %>% 
      select(Tk,Jmax)}, by='Tk')

dat_pa_ww <- inner_join(
  {vc_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=0,cw=0,wc=0,ww=1) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Vcmax=Photo) %>% 
      select(Tk,Vcmax)},
  {jm_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=0,cw=0,wc=0,ww=1) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Jmax=Photo) %>% 
      select(Tk,Jmax)}, by='Tk')


bind_rows(dat_pa_cc %>% mutate(treat="c.c"), 
          dat_pa_ww %>% mutate(treat="w.w")) %>% 
  ggplot(data=., aes(Tk-273.15, Jmax/Vcmax, color=treat))+
  geom_line(alpha=1)+
  theme_bw()+
  geom_point(data=(dat_all %>% select(Treat,Tk,Vcmax,Jmax) %>% 
                     rename(treat=Treat) %>% filter(treat %in% c("c.c","w.w"))), 
             aes(Tk-273.15,Jmax/Vcmax), pch=20, 
             size=4, alpha=0.4,
             show.legend = F)+
  scale_color_manual(values=c("c.c"='blue',"w.w"='red'))+
  labs(y=bquote(paste('J'['max']/'V'['Cmax'])), 
       xlab=NULL
       # x=bquote(paste('Leaf Temperature ('^degree,'C)'))
  )+
  # scale_x_continuous(limits=c(25,46))+
  scale_y_continuous(limits=c(0,2),
                     breaks = c(0,0.5,1,1.5),
                     expand = c(0,0))+
  annotate("text", x = 25, y = c(0.4, 0.3), hjust=0,size=4,
           # label = c("paste(Treatment %->% \" Treatment\")", 
           #           "paste(Treatment %->% \" Control\")"),
           label = c("Control", 
                     "Treatment"),
           parse=T)+
  annotate("point", x = 24, y = c(0.4,0.3), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,16))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = 'none')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expansion(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        axis.title.x.bottom = element_text(size=20))+
  xlab(label = NULL)
ggsave(paste0("figures/jmax_vmax_ratio_c.c_w.w_",Sys.Date(),".png"),
       width = 120, height=100, units='mm')

# END **************************************************************************





# panel (b) control & control->treatment -------------------------------------------------
dat_pb_cc <- inner_join(
  {vc_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=1,cw=0,wc=0,ww=0) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Vcmax=Photo) %>% 
      select(Tk,Vcmax)},
  {jm_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=1,cw=0,wc=0,ww=0) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Jmax=Photo) %>% 
      select(Tk,Jmax)}, by='Tk')

dat_pb_cw <- inner_join(
  {vc_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=0,cw=1,wc=0,ww=0) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Vcmax=Photo) %>% 
      select(Tk,Vcmax)},
  {jm_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=0,cw=1,wc=0,ww=0) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Jmax=Photo) %>% 
      select(Tk,Jmax)}, by='Tk')

bind_rows(dat_pb_cc %>% mutate(treat='c.c'), 
          dat_pb_cw %>% mutate(treat="c.w")) %>% 
  ggplot(data=., aes(Tk-273.15, Jmax/Vcmax, color=treat))+
  geom_line(alpha=1)+
  theme_bw()+
  geom_point(data=(dat_all %>% select(Treat,Tk,Vcmax,Jmax) %>% 
                     rename(treat=Treat) %>% filter(treat %in% c("c.c","c.w"))), 
             aes(Tk-273.15,Jmax/Vcmax,shape=treat), 
             # shape=1, 
             size=4, alpha=0.4,
             show.legend = F)+
  scale_shape_manual(values=c("c.c"=16,"c.w"=1))+
  scale_color_manual(values=c("c.c"='blue',"c.w"='red'))+
  labs(y=bquote(paste('J'['max']/'V'['Cmax'])), 
       xlab=NULL
       # x=bquote(paste('Leaf Temperature ('^degree,'C)'))
  )+
  scale_x_continuous(limits=c(25,46))+
  scale_y_continuous(limits=c(0,2), 
                     breaks=c(0,0.5,1,1.5),
                     expand = c(0,0))+
  annotate("text", x = 25, y = c(0.4, 0.3), hjust=0,size=4,
           label = c("Control",
                     "paste(Control %->% \" Treatment\")"),
           # label = c("Control", 
           #           "Treatment"),
           parse=T)+
  annotate("point", x = 24, y = c(0.4,0.3), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = 'none')+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_blank(), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        axis.title.x.bottom = element_text(size=20))+
  xlab(label = NULL)
ggsave(paste0("figures/jmax_vmax_ratio_c.c_c.w_",Sys.Date(),".png"),
       width = 120, height=100, units='mm')
# END **************************************************************************


# panel (c) treatment & treatment->control -------------------------------------------------
dat_pc_ww <- inner_join(
  {vc_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=0,cw=0,wc=0,ww=1) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Vcmax=Photo) %>% 
      select(Tk,Vcmax)},
  {jm_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=0,cw=0,wc=0,ww=1) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Jmax=Photo) %>% 
      select(Tk,Jmax)}, by='Tk')

dat_pc_wc <- inner_join(
  {vc_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=0,cw=0,wc=1,ww=0) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Vcmax=Photo) %>% 
      select(Tk,Vcmax)},
  {jm_fit %>% as.data.frame() %>% 
      as.data.frame() %>% 
      mode_hdi(.width = 0.5) %>% 
      rename_all(fn) %>% 
      expand(nesting(kopt,kcw,kwc,kww,
                     Ha,Hacw,Hawc,Haww,
                     Topt,Tcw,Twc,Tww),
             Tk=seq(295,321,0.1), 
             cc=0,cw=0,wc=1,ww=0) %>% 
      mutate(Photo=pmap(., f_photo_mp)) %>% 
      unnest(cols='Photo') %>% 
      ungroup() %>% 
      rename(Jmax=Photo) %>% 
      select(Tk,Jmax)}, by='Tk')

bind_rows(dat_pc_ww %>% mutate(treat='w.w'), 
          dat_pc_wc %>% mutate(treat='w.c')) %>% 
  ggplot(data=., aes(Tk-273.15, Jmax/Vcmax, color=treat))+
  geom_line(alpha=1)+
  theme_bw()+
  geom_point(data=(dat_all %>% select(Treat,Tk,Vcmax,Jmax) %>% 
                     rename(treat=Treat) %>% filter(treat %in% c("w.w","w.c"))), 
             aes(Tk-273.15,Jmax/Vcmax,shape=treat), 
             # shape=1, 
             size=4, alpha=0.4,
             show.legend = F)+
  scale_shape_manual(values=c("w.w"=16,"w.c"=1))+
  scale_color_manual(values=c("w.w"='red',"w.c"='blue'))+
  labs(y=bquote(paste('J'['max']/'V'['Cmax'])),
       x=bquote(paste('Leaf Temperature ('^degree,'C)'))
  )+
  scale_x_continuous(limits=c(25,46))+
  scale_y_continuous(limits=c(0,2),
                     breaks=c(0,0.5,1,1.5),
                     expand = c(0,0))+
  annotate("text", x = 25, y = c(0.4, 0.3), hjust=0,size=4,
           label = c("Treatment",
                     "paste(Treatment %->% \" Control\")"),
           # label = c("Control", 
           #           "Treatment"),
           parse=T)+
  annotate("point", x = 24, y = c(0.4,0.3), 
           colour = c('red',"blue"), 
           size = 5, alpha=0.5, shape=c(16,1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = 'none')+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expansion(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        axis.title.x.bottom = element_text(size=20))
ggsave(paste0("figures/jmax_vmax_ratio_w.w_w.c_",Sys.Date(),".png"),
       width = 120, height=100, units='mm')
# END **************************************************************************



library(magick)
i1 <- image_read("figures/jmax_vmax_ratio_c.c_w.w_2021-03-16.png")
i2 <- image_read("figures/jmax_vmax_ratio_c.c_c.w_2021-03-16.png")
i3 <- image_read("figures/jmax_vmax_ratio_w.w_w.c_2021-03-16.png")
i1 <- image_annotate(i1, "(a)", size = 80, color = "black",
                     location = "+1+1")
i2 <- image_annotate(i2, "(b)", size = 80, color = "black",
                     location = "+1+1")
i3 <- image_annotate(i3, "(c)", size = 80, color = "black",
                     location = "+1+1")
ijoin_ratio <- image_append(c(i1,i2,i3), stack=T)
ijoin_ratio
# img_top <- c(i1,i2)
# img_mid <- c(i3,i4)
# img_bottom <- c(i5,i6)
# ijoin_top <- image_append(img_top)
# ijoin_mid <- image_append(img_mid)
# ijoin_bottom <- image_append(img_bottom)
# ijoin_P400P800 <- image_append(c(ijoin_top,ijoin_mid, ijoin_bottom), stack=T)

image_write(ijoin_ratio, 
            path = paste0("figures/Fig4_Ratio_Jmax_Vcmax_",Sys.Date(),".png"), 
            format = "png")







