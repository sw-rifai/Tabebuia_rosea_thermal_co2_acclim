library(tidyverse); library(bayesplot);
library(rstan)
options(mc.cores = parallel::detectCores()-1)
set.seed(3); 

# params.all<-read_csv("data/params.all.28032018.csv")
# params.all$Tk<-params.all$Tleaf+273.15
dat270 <- read_csv("data/params_20190124.csv") %>% 
  rename(Photo = P270, Tk=Tleaf) %>% 
  select(-P505) %>% 
  mutate(Tk=Tk+273.15)

dat505 <- read_csv("data/params_20190124.csv") %>% 
  rename(Photo = P505, Tk=Tleaf) %>% 
  select(-P270) %>% 
  mutate(Tk=Tk+273.15)

#*****************************************************************************
# Photo 270 FIT TREATMENTS c.c & c.w w.c w.w                 ----------------
#*****************************************************************************
library(rstan)
set.seed(3)
photo270_c.c <- stan(file = "Stan/FitPhoto.stan",
                     data = list(A=dat270 %>% filter(Treat=="c.c") %>% pull(Photo),
                                 Tk=dat270 %>% filter(Treat=="c.c") %>% pull(Tk),
                                 N = dat270 %>% filter(Treat=="c.c") %>% dim() %>% .[1]),
                     warmup = 4000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo270_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
# mcmc_trace(as.array(photo270_c.c), c("sigma","Topt","Ha","kopt"))
# mcmc_dens(as.array(photo270_c.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

photo270_c.w <- stan(file = "Stan/FitPhoto.stan",
                     data=list(A=dat270 %>% filter(Treat=="c.w") %>% pull(Photo), 
                               Tk=dat270 %>% filter(Treat=="c.w") %>% pull(Tk),
                               N=dat270 %>% filter(Treat=="c.w") %>% dim() %>% .[1]),
                     warmup = 4000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.99))
print(photo270_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
# mcmc_trace(as.array(photo270_c.w), c("sigma","Topt","Ha","kopt"))
# mcmc_dens(as.array(photo270_c.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()
#*******************************************************************************

#*****************************************************************************
# FIT TREATMENTS w.w & w.c ----------------
#*****************************************************************************
library(rstan)
set.seed(3)
photo270_w.w <- stan(file = "Stan/FitPhoto.stan",
                     data = list(A=dat270 %>% filter(Treat=="w.w") %>% pull(Photo),
                                 Tk=dat270 %>% filter(Treat=="w.w") %>% pull(Tk),
                                 N = dat270 %>% filter(Treat=="w.w") %>% dim() %>% .[1]),
                     warmup = 4000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo270_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
# mcmc_trace(as.array(photo270_w.w), c("sigma","Topt","Ha","kopt"))
# mcmc_dens(as.array(photo270_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

photo270_w.c <- stan(file = "Stan/FitPhoto.stan",
                     data=list(A=dat270 %>% filter(Treat=="w.c") %>% pull(Photo), 
                               Tk=dat270 %>% filter(Treat=="w.c") %>% pull(Tk),
                               N=dat270 %>% filter(Treat=="w.c") %>% dim() %>% .[1]),
                     warmup = 4000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.95))
print(photo270_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
# mcmc_trace(as.array(photo270_w.c), c("sigma","Topt","Ha","kopt"))
# mcmc_dens(as.array(photo270_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()
#*******************************************************************************

#*****************************************************************************
# Photo 505 FIT TREATMENTS c.c & c.w w.c w.w                 ----------------
#*****************************************************************************
library(rstan)
set.seed(3)
photo505_c.c <- stan(file = "Stan/FitPhoto.stan",
                     data = list(A=dat505 %>% filter(Treat=="c.c") %>% pull(Photo),
                                 Tk=dat505 %>% filter(Treat=="c.c") %>% pull(Tk),
                                 N = dat505 %>% filter(Treat=="c.c") %>% dim() %>% .[1]),
                     warmup = 4000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo505_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
# mcmc_trace(as.array(photo505_c.c), c("sigma","Topt","Ha","kopt"))
# mcmc_dens(as.array(photo505_c.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


photo505_c.w <- stan(file = "Stan/FitPhoto.stan",
                     data=list(A=dat505 %>% filter(Treat=="c.w") %>% pull(Photo), 
                               Tk=dat505 %>% filter(Treat=="c.w") %>% pull(Tk),
                               N=dat505 %>% filter(Treat=="c.w") %>% dim() %>% .[1]),
                     warmup = 4000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.99))
print(photo505_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
# mcmc_trace(as.array(photo505_c.w), c("sigma","Topt","Ha","kopt"))
# mcmc_dens(as.array(photo505_c.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()
#*******************************************************************************

#*****************************************************************************
# FIT TREATMENTS w.w & w.c ----------------
#*****************************************************************************
library(rstan)
set.seed(3)
photo505_w.w <- stan(file = "Stan/FitPhoto.stan",
                     data = list(A=dat505 %>% filter(Treat=="w.w") %>% pull(Photo),
                                 Tk=dat505 %>% filter(Treat=="w.w") %>% pull(Tk),
                                 N = dat505 %>% filter(Treat=="w.w") %>% dim() %>% .[1]),
                     warmup = 4000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(photo505_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
# mcmc_trace(as.array(photo505_w.w), c("sigma","Topt","Ha","kopt"))
# mcmc_dens(as.array(photo505_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

photo505_w.c <- stan(file = "Stan/FitPhoto.stan",
                     data=list(A=dat505 %>% filter(Treat=="w.c") %>% pull(Photo), 
                               Tk=dat505 %>% filter(Treat=="w.c") %>% pull(Tk),
                               N=dat505 %>% filter(Treat=="w.c") %>% dim() %>% .[1]),
                     warmup = 4000, save_dso = T,
                     iter=8000, thin=2, chains=3, verbose=T,
                     cores=3,
                     control=list(adapt_delta=0.95))
print(photo505_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
# mcmc_trace(as.array(photo505_w.c), c("sigma","Topt","Ha","kopt"))
# mcmc_dens(as.array(photo505_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()
#*******************************************************************************



# Functions ---------------------------------------------------------------
pow <- function(alpha,beta){return(alpha^beta)}
# f_photo <- function(kopt,Ha,Tk,Topt, ...){
#   photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
#                      (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
#   return(photo)}
f_photo <- function(kopt,Ha,Tk,Topt, ...){
  pow <- function(alpha,beta){return(alpha^beta)}
  out <- kopt * ((200*pow(2.718282,((Ha*(Tk-Topt))/(Tk*0.008314*Topt))))/ 
                   (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(out)
}






#########################################################################################
# FOR THE TOP ROW
#-########################################################################################
#
# Plot Stan Fit P270 c.c and w.w -----------------------------------------------
p1 <- photo270_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  labs(y=bquote(paste('P'[270],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,40)+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))+
  labs(x=NULL)
p1 <- p1 + annotate("text", x = 36.5+2, y = c(37,33), hjust=0,size=4,
                    label = c("Control", 
                              "Treatment"),parse=T)+
  annotate("point", x = 35.5+2, y = c(37,33), 
           colour = c("blue","red"), 
           size = 5, shape=c(16,16), alpha=0.5)
tmp1 <- (dat270 %>% filter(Treat=="w.w")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=5, alpha=0.5,
                      col='red', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(photo270_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(photo270_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(photo270_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(photo270_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(photo270_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p2_lines <- photo270_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp2 <- (dat270 %>% filter(Treat=="c.c")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, alpha=0.5, 
                        col='blue', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(photo270_c.c)$summary["Topt","50%"]-273.15), 
                      col='blue')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(photo270_c.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(photo270_c.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(photo270_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(photo270_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

p270_w.w_c.c <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/photo_Topt/p270_c.c_w.w_",Sys.Date(),".png"),plot=p270_w.w_c.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/photo_Topt/p270_c.c_w.w_",Sys.Date(),".pdf"),plot=p270_w.w_c.c, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


# Plot Stan Fit P505 c.c and w.w -----------------------------------------------
p1 <- photo505_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  labs(y=bquote(paste('P'[505],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=NULL
       # x=bquote(paste('Leaf Temperature ('^degree,'C)'))
  )+
  scale_y_continuous(position = "right", limits = c(0,40))+         ### ONLY FOR P505!
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=20),
        axis.title.x.bottom = element_text(size=20),
        axis.text.y.right = element_text(size=20), 
        axis.title.y.right = element_text(size=20), 
        plot.margin = margin(l=20))

tmp1 <- (dat505 %>% filter(Treat=="w.w")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=5, alpha=0.5,
                      col='red', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(photo505_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(photo505_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(photo505_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(photo505_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(photo505_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p2_lines <- photo505_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp2 <- (dat505 %>% filter(Treat=="c.c")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, alpha=0.5, 
                        col='blue', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(photo505_c.c)$summary["Topt","50%"]-273.15), 
                      col='blue')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(photo505_c.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(photo505_c.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(photo505_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(photo505_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

p505_w.w_c.c <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/photo_Topt/p505_c.c_w.w_",Sys.Date(),".png"),plot=p505_w.w_c.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/photo_Topt/p505_c.c_w.w_",Sys.Date(),".pdf"),plot=p505_w.w_c.c, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-















##########################################################################################
# FOR THE MIDDLE ROW
##########################################################################################
#
# Plot Stan Fit P270 c.c and c.w -----------------------------------------------
p1 <- photo270_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='blue')+
  theme_bw()+
  labs(y=bquote(paste('P'[270],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=NULL)+
  # x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  ylim(0,40)+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))
p1 <- p1 + annotate("text", x = 36+2, y = c(37,33), hjust=0,size=4,
                    label = c("Control", 
                              "paste(Control %->% \" Treatment\")"),parse=T)+
  annotate("point", x = 35+2, y = c(37,33), 
           colour = c('blue',"red"), 
           size = 5, shape=c(16,1), alpha=0.5)
tmp1 <- (dat270 %>% filter(Treat=="c.c")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=5, alpha=0.5,
                      col='blue', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(photo270_c.c)$summary["Topt","50%"]-273.15)), 
                      col='blue')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(photo270_c.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(photo270_c.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(photo270_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(photo270_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)


p2_lines <- photo270_c.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='red')
tmp2 <- (dat270 %>% filter(Treat=="c.w")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, shape=1,
                        col='red', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(photo270_c.w)$summary["Topt","50%"]-273.15), 
                      col='red')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(photo270_c.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(photo270_c.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(photo270_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(photo270_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)

p270_cc_cw <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(plot = p270_cc_cw, 
       paste0("figures/photo_Topt/p270_c.c_c.w_",Sys.Date(),".png"), 
       width = 120, height=100, units='mm')
ggsave(plot = p270_cc_cw, 
       paste0("figures/photo_Topt/p270_c.c_c.w_",Sys.Date(),".pdf"), 
       width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


# Plot Stan Fit P505 c.c and c.w -----------------------------------------------
p1 <- photo505_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='blue')+
  theme_bw()+
  labs(y=bquote(paste('P'[505],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=NULL)+
  # x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(position = 'right', limits = c(0,40))+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.right = element_text(size=18), 
        axis.title.y.right = element_text(size=20), 
        plot.margin = margin(l=20))
# p1 <- p1 + annotate("text", x = 36, y = c(37,33), hjust=0,size=4,
#                     label = c("Control", 
#                               "paste(Control %->% \" Treatment\")"),parse=T)+
#   annotate("point", x = 35, y = c(37,33), 
#            colour = c('blue',"red"), 
#            size = 5, shape=c(16,1))
tmp1 <- (dat505 %>% filter(Treat=="c.c")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=5, alpha=0.5,
                      col='blue', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(photo505_c.c)$summary["Topt","50%"]-273.15)), 
                      col='blue')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(photo505_c.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(photo505_c.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(photo505_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(photo505_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)


p2_lines <- photo505_c.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='red')
tmp2 <- (dat505 %>% filter(Treat=="c.w")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, shape=1,
                        col='red', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(photo505_c.w)$summary["Topt","50%"]-273.15), 
                      col='red')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(photo505_c.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(photo505_c.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(photo505_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(photo505_c.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)

p505_cc_cw <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(plot = p505_cc_cw, 
       paste0("figures/photo_Topt/p505_c.c_c.w_",Sys.Date(),".png"), 
       width = 120, height=100, units='mm')
ggsave(plot = p505_cc_cw, 
       paste0("figures/photo_Topt/p505_c.c_c.w_",Sys.Date(),".pdf"), 
       width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-








##########################################################################################
# FOR THE BOTTOM ROW
##########################################################################################
# Plot Stan Fit P270 w.w and w.c -----------------------------------------------
p1 <- photo270_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  labs(y=bquote(paste('P'[270],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(position = 'left', limits = c(0,40))+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18),
        axis.title.x.bottom = element_text(size=20),
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))

p1 <- p1 + annotate("text", 
                    x = 36.5+1, 
                    y = c(37,33), hjust=0,size=4,
                    label = c("Treatment", 
                              "paste(Treatment %->% \"Control\")"),parse=T)+
  annotate("point", x = 35.5+1, y = c(37,33), 
           colour = c('red',"blue"), 
           size = 5,alpha=0.5, shape=c(16,1))
tmp1 <- (dat270 %>% filter(Treat=="w.w")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=5, alpha=0.5,
                      col='red', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(photo270_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(photo270_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(photo270_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(photo270_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(photo270_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p2_lines <- photo270_w.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp2 <- (dat270 %>% filter(Treat=="w.c")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, shape=1, 
                        col='blue', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(photo270_w.c)$summary["Topt","50%"]-273.15), 
                      col='blue')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(photo270_w.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(photo270_w.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(photo270_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(photo270_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

p270_w.w_w.c <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/photo_Topt/p270_w.w_w.c_",Sys.Date(),".png"),plot=p270_w.w_w.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/photo_Topt/p270_w.w_w.c_",Sys.Date(),".pdf"),plot=p270_w.w_w.c, width = 120, height=100, units='mm')
gc(reset = T, full = T)
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-



# Plot Stan Fit P505 w.w and w.c -----------------------------------------------
p1 <- photo505_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  labs(y=bquote(paste('P'[505],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(position = 'right', limits = c(0,40))+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.title.x.bottom = element_text(size=20),
        axis.text.y.right = element_text(size=18), 
        axis.title.y.right = element_text(size=20), 
        plot.margin = margin(l=20))

# p1 <- p1 + annotate("text", 
#                     x = 36.5-3, 
#                     y = c(37,33), hjust=0,size=5,
#                     label = c("Treatment", 
#                               "paste(Treatment %->% \"Control\")"),parse=T)+
#   annotate("point", x = 35.5-3, y = c(37,33), 
#            colour = c('red',"blue"), 
#            size = 5,alpha=0.5, shape=c(16,1))
tmp1 <- (dat505 %>% filter(Treat=="w.w")) %>% select(Tk,Photo)
p1 <- p1 + geom_point(data=tmp1, aes(Tk-273.15,Photo), size=5, alpha=0.5,
                      col='red', show.legend = F, inherit.aes = F)
p1_topt <- geom_vline(aes(xintercept=(summary(photo505_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p1_topt_2.5 <- geom_vline(aes(xintercept=summary(photo505_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p1_topt_97.5 <- geom_vline(aes(xintercept=summary(photo505_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p1_topt_10 <- geom_vline(aes(
  xintercept=summary(photo505_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p1_topt_90 <- geom_vline(aes(
  xintercept=summary(photo505_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p2_lines <- photo505_w.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(298.15,320,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp2 <- (dat505 %>% filter(Treat=="w.c")) %>% select(Tk,Photo)
p2_points <- geom_point(data=tmp2, aes(Tk-273.15,Photo), size=5, shape=1, 
                        col='blue', show.legend = F, inherit.aes = F)
p2_topt <- geom_vline(aes(xintercept=summary(photo505_w.c)$summary["Topt","50%"]-273.15), 
                      col='blue')
p2_topt_2.5 <- geom_vline(aes(xintercept=summary(photo505_w.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p2_topt_97.5 <- geom_vline(aes(xintercept=summary(photo505_w.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p2_topt_10 <- geom_vline(aes(
  xintercept=summary(photo505_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p2_topt_90 <- geom_vline(aes(
  xintercept=summary(photo505_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

p505_w.w_w.c <- p1+p2_lines+p2_points+p1_topt+p2_topt+p1_topt_10+p1_topt_90+p2_topt_10+p2_topt_90
ggsave(paste0("figures/photo_Topt/p505_w.w_w.c_",Sys.Date(),".png"),plot=p505_w.w_w.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/photo_Topt/p505_w.w_w.c_",Sys.Date(),".pdf"),plot=p505_w.w_w.c, width = 120, height=100, units='mm')
gc(reset = T, full = T)

#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-









library(magick)
i1 <- image_read("figures/photo_Topt/p270_c.c_w.w_2019-06-11.png")
i2 <- image_read("figures/photo_Topt/p505_c.c_w.w_2019-06-11.png")
i3 <- image_read("figures/photo_Topt/p270_c.c_c.w_2019-06-11.png")
i4 <- image_read("figures/photo_Topt/p505_c.c_c.w_2019-06-11.png")
i5 <- image_read("figures/photo_Topt/p270_w.w_w.c_2019-06-11.png")
i6 <- image_read("figures/photo_Topt/p505_w.w_w.c_2019-06-11.png")
i1 <- image_annotate(i1, "(a)", size = 100, color = "black",
                     location = "+320+20")
i2 <- image_annotate(i2, "(b)", size = 100, color = "black",
                     location = "+85+20")
i3 <- image_annotate(i3, "(c)", size = 100, color = "black",
                     location = "+320+20")
i4 <- image_annotate(i4, "(d)", size = 100, color = "black",
                     location = "+85+20")
i5 <- image_annotate(i5, "(e)", size = 100, color = "black",
                     location = "+320+20")
i6 <- image_annotate(i6, "(f)", size = 100, color = "black",
                     location = "+85+20")
# i2 <- image_border(i2, color='white',geometry = '100x1')
# i4 <- image_border(i4, color='white',geometry = '100x1')
# i6 <- image_border(i6, color='white',geometry = '100x1')

img_top <- c(i1,i2)
img_mid <- c(i3,i4)
img_bottom <- c(i5,i6)
ijoin_top <- image_append(img_top)
ijoin_mid <- image_append(img_mid)
ijoin_bottom <- image_append(img_bottom)
ijoin_P270P505 <- image_append(c(ijoin_top,ijoin_mid, ijoin_bottom), stack=T)

image_write(ijoin_P270P505, path = paste0("figures/photo_Topt/Joined_P270P505_",Sys.Date(),".png"),
            format = "png")













