library(tidyverse); library(bayesplot);
library(rstan)

params.all<-read_csv("data/params.all.28032018.csv")
dat400 <- read_csv("data/P400.2904.csv")
dat800 <- read_csv("data/P800.2904.csv")
params.all$Tk<-params.all$Tleaf+273.15

# Plot Stan Fit P400 c.c and w.w -----------------------------------------------
f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}
nsamps <- 1000

p1 <- photo400_w.w %>%  
  as.data.frame() %>% 
  sample_n(nsamps) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  labs(y=bquote(paste('P'[400],'(',mu,'mol m'^'-2','s'^'-1',')')), 
       x=NULL,
       # x=bquote(paste('Leaf Temperature ('^degree,'C)'))
  )+
  ylim(0,40)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18),
        axis.title.x.bottom = element_text(size=18),
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20))


p1 <- p1 + annotate("text", x = 36.5, y = c(37,33), hjust=0,size=4,
                    label = c("Control", 
                              "Treatment"),parse=T)+
  annotate("point", x = 35.5, y = c(37,33), 
           colour = c("blue","red"), 
           size = 5, shape=c(16,16))
tmp1 <- (dat400 %>% filter(Treat=="w.w")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=5, alpha=0.5,
                      col='red', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(photo400_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(photo400_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(photo400_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(photo400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(photo400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p2_lines <- photo400_c.c %>%  
  as.data.frame() %>% 
  sample_n(nsamps) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp2 <- (dat400 %>% filter(Treat=="c.c")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, 
                        col='blue', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(photo400_c.c)$summary["Topt","50%"]-273.15), 
                      col='blue')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(photo400_c.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(photo400_c.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(photo400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(photo400_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

p400_w.w_c.c <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave("figures/p400_c.c_w.w_.png",plot=p400_w.w_c.c, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


# Plot Stan Fit w.w and w.c -----------------------------------------------
f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}


p1 <- photo400_w.w %>%  
  as.data.frame() %>% 
  sample_n(nsamps) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  labs(y=bquote(paste('P'[400],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)'))
       )+
  ylim(0,40)+
  scale_x_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18),
        axis.title.x.bottom = element_text(size=18),
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20))
  
p1 <- p1 + annotate("text", x = 36.5, y = c(37,33), hjust=0,size=4,
                    label = c("Treatment", 
                              "paste(Treatment %->% \" Control\")"),parse=T)+
  annotate("point", x = 35.5, y = c(37,33), 
           colour = c('red',"blue"), 
           size = 5, shape=c(16,1))
tmp1 <- (dat400 %>% filter(Treat=="w.w")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=5, alpha=0.5,
                      col='red', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(photo400_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(photo400_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(photo400_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(photo400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(photo400_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p2_lines <- photo400_w.c %>%  
  as.data.frame() %>% 
  sample_n(nsamps) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp2 <- (dat400 %>% filter(Treat=="w.c")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, shape=1,
                        col='blue', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(photo400_w.c)$summary["Topt","50%"]-273.15), 
                      col='blue')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(photo400_w.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(photo400_w.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(photo400_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(photo400_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

p400_w.w_w.c <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/p400_w.w_w.c_",Sys.Date(),".png"),plot=p400_w.w_w.c, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-



