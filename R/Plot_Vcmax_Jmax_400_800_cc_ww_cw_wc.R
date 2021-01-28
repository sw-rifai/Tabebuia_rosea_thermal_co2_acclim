
f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
# Plot Stan Vcmax fit cc and ww -----------------------------------------------
p1 <- vc400_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,319,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='blue')+
  theme_bw()+
  labs(y=bquote(paste('V'['cmax'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,325)+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20))+
  labs(x=NULL)

p1 <- p1 + annotate("text", x = 22.5, y = c(295+30,275+30), hjust=0,size=4,
                    label = c("Control", 
                              "Treatment"),parse=T)+
  annotate("point", x = 28.5, y = c(295+30,275+30), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,16))
tmp1 <- (dat_all %>% filter(Treat=="c.c")) %>% select(Tk,Vcmax)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Vcmax), size=5, alpha=0.5,
                      col='blue', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(vc400_c.c)$summary["Topt","50%"]-273.15)), 
                      col='blue')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(vc400_c.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(vc400_c.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(vc400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(vc400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)


p2_lines <- vc400_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,319,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='red')
  # scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))

tmp2 <- (dat_all %>% filter(Treat=="w.w")) %>% select(Tk,Vcmax)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Vcmax), size=5, shape=16,alpha=0.5,
                        col='red', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(vc400_w.w)$summary["Topt","50%"]-273.15), 
                      col='red')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(vc400_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(vc400_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(vc400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(vc400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)

vc400_c.c_w.w <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90

ggsave(paste0("figures/vcmax400_c.c_w.w_",Sys.Date(),".png"),plot=vc400_c.c_w.w, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/vcmax400_c.c_w.w_",Sys.Date(),".svg"),
       plot=vc400_c.c_w.w, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/vcmax400_c.c_w.w_",Sys.Date(),".pdf"),
       plot=vc400_c.c_w.w, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
# Plot Stan Vcmax fit cc and cw -----------------------------------------------
p1 <- vc400_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='blue')+
  theme_bw()+
  labs(y=bquote(paste('V'['cmax'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,325)+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20))+
  labs(x=NULL)


p1 <- p1 + annotate("text", x = 22.5, y = c(295+30,275+30), hjust=0,size=4,
                    label = c("paste(Control %->% \" Control\")", 
                              "paste(Control %->% \" Treatment\")"),parse=T)+
  annotate("point", x = 33.5, y = c(295+30,275+30), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,1))
tmp1 <- (dat_all %>% filter(Treat=="c.c")) %>% select(Tk,Vcmax)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Vcmax), size=5, alpha=0.5,
                      col='blue', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(vc400_c.c)$summary["Topt","50%"]-273.15)), 
                      col='blue')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(vc400_c.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(vc400_c.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(vc400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(vc400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)
p1 <- p1+ theme(axis.text.x.bottom = element_text(size=18), 
                axis.text.y.left = element_text(size=18), 
                axis.title.y.left = element_text(size=20))+
  labs(x=NULL)


p2_lines <- vc400_c.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='red')
tmp2 <- (dat_all %>% filter(Treat=="c.w")) %>% select(Tk,Vcmax)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Vcmax), size=5, shape=1,alpha=0.5,
                        col='red', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(vc400_c.w)$summary["Topt","50%"]-273.15), 
                      col='red')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(vc400_c.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(vc400_c.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(vc400_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(vc400_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)

vc400_c.c_c.w <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/vcmax400_c.c_c.w_",Sys.Date(),".png"),
       plot=vc400_c.c_c.w, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/vcmax400_c.c_c.w_",Sys.Date(),".svg"),
       plot=vc400_c.c_c.w, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/vcmax400_c.c_c.w_",Sys.Date(),".pdf"),
       plot=vc400_c.c_c.w, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
# Plot Stan Vcmax fit ww and wc -----------------------------------------------
p1 <- vc400_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  labs(y=bquote(paste('V'['cmax'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,325)+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20),
        axis.title.x.bottom = element_text(size=20))

p1 <- p1 + annotate("text", x = 22.5, y = c(295+30,275+30), hjust=0,size=4,
                    label = c("paste(Treatment %->% \" Treatment\")", 
                              "paste(Treatment %->% \" Control\")"),parse=T)+
  annotate("point", x = 35.5, y = c(295+30,275+30), 
           colour = c('red',"blue"), 
           size = 5, alpha=0.5, shape=c(16,1))
tmp1 <- (dat_all %>% filter(Treat=="w.w")) %>% select(Tk,Vcmax)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Vcmax), size=5, alpha=0.5,
                      col='red', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(vc400_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(vc400_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(vc400_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(vc400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(vc400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p2_lines <- vc400_w.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp2 <- (dat_all %>% filter(Treat=="w.c")) %>% select(Tk,Vcmax)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Vcmax), size=5, shape=1,alpha=0.5,
                        col='blue', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(vc400_w.c)$summary["Topt","50%"]-273.15), 
                      col='blue')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(vc400_w.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(vc400_w.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(vc400_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(vc400_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

vc400_w.w_w.c <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/vcmax400_w.w_w.c_",Sys.Date(),".png"),plot=vc400_w.w_w.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/vcmax400_w.w_w.c_",Sys.Date(),".svg"),plot=vc400_w.w_w.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/vcmax400_w.w_w.c_",Sys.Date(),".pdf"),plot=vc400_w.w_w.c, width = 120, height=100, units='mm')

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


###############################################################################
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
# Plot Stan Jmax fit cc and ww -----------------------------------------------
p1 <- jm400_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,319,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='blue')+
  theme_bw()+
  labs(y=bquote(paste('J'['max'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,325)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20))+
  labs(x=NULL)


p1 <- p1 + annotate("text", x = 23, y = c(295+30,275+30), hjust=0,size=4,
                    label = c("Control", 
                              "Treatment"),parse=T)+
  annotate("point", x = 28.5, y = c(295+30,275+30), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,16))
tmp1 <- (dat_all %>% filter(Treat=="c.c")) %>% select(Tk,Jmax)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Jmax), size=5, alpha=0.5,
                      col='blue', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(jm400_c.c)$summary["Topt","50%"]-273.15)), 
                      col='blue')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(jm400_c.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(jm400_c.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(jm400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(jm400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)


p2_lines <- jm400_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,319,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='red')
tmp2 <- (dat_all %>% filter(Treat=="w.w")) %>% select(Tk,Jmax)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Jmax), size=5, shape=16,alpha=0.5,
                        col='red', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(jm400_w.w)$summary["Topt","50%"]-273.15), 
                      col='red')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(jm400_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(jm400_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(jm400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(jm400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)

jm400_c.c_w.w <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/jmax400_c.c_w.w_",Sys.Date(),".png"),plot=jm400_c.c_w.w, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/jmax400_c.c_w.w_",Sys.Date(),".pdf"),plot=jm400_c.c_w.w, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/jmax400_c.c_w.w_",Sys.Date(),".svg"),plot=jm400_c.c_w.w, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

# Plot Stan Jmax fit cc and cw -----------------------------------------------
p1 <- jm400_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='blue')+
  theme_bw()+
  labs(y=bquote(paste('J'['max'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,325)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20))+
  labs(x=NULL)

p1 <- p1 + annotate("text", x = 22.5, y = c(295+30,275+30), hjust=0,size=4,
                    label = c("paste(Control %->% \" Control\")", 
                              "paste(Control %->% \" Treatment\")"),parse=T)+
  annotate("point", x = 34.5, y = c(295+30,275+30), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,1))
tmp1 <- (dat_all %>% filter(Treat=="c.c")) %>% select(Tk,Jmax)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Jmax), size=5, alpha=0.5,
                      col='blue', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(jm400_c.c)$summary["Topt","50%"]-273.15)), 
                      col='blue')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(jm400_c.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(jm400_c.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(jm400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(jm400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)


p2_lines <- jm400_c.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='red')
tmp2 <- (dat_all %>% filter(Treat=="c.w")) %>% select(Tk,Jmax)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Jmax), size=5, shape=1,alpha=0.5,
                        col='red', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(jm400_c.w)$summary["Topt","50%"]-273.15), 
                      col='red')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(jm400_c.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(jm400_c.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(jm400_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(jm400_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)

jm400_c.c_c.w <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/jmax400_c.c_c.w_",Sys.Date(),".png"),plot=jm400_c.c_c.w, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/jmax400_c.c_c.w_",Sys.Date(),".svg"),plot=jm400_c.c_c.w, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/jmax400_c.c_c.w_",Sys.Date(),".pdf"),plot=jm400_c.c_c.w, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


# Plot Stan Jmax fit ww and wc -----------------------------------------------
p1 <- jm400_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  labs(y=bquote(paste('J'['max'],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,325)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        axis.title.x.bottom = element_text(size=20))

p1 <- p1 + annotate("text", x = 21+1.5, y = c(295+30,275+30), hjust=0,size=4,
                    label = c("paste(Treatment %->% \" Treatment\")", 
                              "paste(Treatment %->% \" Control\")"),parse=T)+
  annotate("point", x = 34+1.5, y = c(295+30,275+30), 
           colour = c('red',"blue"), 
           size = 5, alpha=0.5, shape=c(16,1))
tmp1 <- (dat_all %>% filter(Treat=="w.w")) %>% select(Tk,Jmax)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Jmax), size=5, alpha=0.5,
                      col='red', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(jm400_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(jm400_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(jm400_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(jm400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(jm400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p2_lines <- jm400_w.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp2 <- (dat_all %>% filter(Treat=="w.c")) %>% select(Tk,Jmax)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Jmax), size=5, shape=1,alpha=0.5,
                        col='blue', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(jm400_w.c)$summary["Topt","50%"]-273.15), 
                      col='blue')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(jm400_w.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(jm400_w.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(jm400_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(jm400_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

jm400_w.w_w.c <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/jmax400_w.w_w.c_",Sys.Date(),".png"),plot=jm400_w.w_w.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/jmax400_w.w_w.c_",Sys.Date(),".svg"),plot=jm400_w.w_w.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/vector_format/jmax400_w.w_w.c_",Sys.Date(),".pdf"),plot=jm400_w.w_w.c, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


library(magick)
i1 <- image_read("figures/vcmax400_c.c_w.w_2019-06-06.png")
i2 <- image_read("figures/jmax400_c.c_w.w_2019-06-06.png")
i3 <- image_read("figures/vcmax400_c.c_c.w_2019-06-06.png")
i4 <- image_read("figures/jmax400_c.c_c.w_2019-06-06.png")
i5 <- image_read("figures/vcmax400_w.w_w.c_2019-06-06.png")
i6 <- image_read("figures/jmax400_w.w_w.c_2019-06-06.png")
i1 <- image_annotate(i1, "(a)", size = 90, color = "black",
                     location = "+150+5")
i2 <- image_annotate(i2, "(b)", size = 90, color = "black",
                     location = "+160+5")
i3 <- image_annotate(i3, "(c)", size = 90, color = "black",
                     location = "+150+5")
i4 <- image_annotate(i4, "(d)", size = 90, color = "black",
                     location = "+160+5")
i5 <- image_annotate(i5, "(e)", size = 90, color = "black",
                     location = "+150+5")
i6 <- image_annotate(i6, "(f)", size = 90, color = "black",
                     location = "+160+5")

img_top <- c(i1,i2)
img_mid <- c(i3,i4)
img_bottom <- c(i5,i6)
ijoin_top <- image_append(img_top)
ijoin_mid <- image_append(img_mid)
ijoin_bottom <- image_append(img_bottom)
ijoin_P400P800 <- image_append(c(ijoin_top,ijoin_mid, ijoin_bottom), stack=T)

image_write(ijoin_P400P800, 
            path = paste0("figures/Joined_Vcmax_Jmax_",Sys.Date(),".png"), format = "png")
