library(tidyverse); library(bayesplot);
library(rstan)
# library(xfun); library(brms); 
options(mc.cores = parallel::detectCores()-1)
set.seed(3); 

params.all<-read_csv("data/params.all.28032018.csv")
dat800 <- read_csv("data/P800.2904.csv")

params.all$Tk<-params.all$Tleaf+273.15

# create empty frame for Medlyn parameters
pars.med<-NULL

# treatments: (w.w, w.w, w.c, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
treat<-as.character(unique(params.all$Treat))


#*****************************************************************************
# FIT TREATMENTS w.w & w.c ----------------
#*****************************************************************************
library(rstan)
set.seed(3)
sphoto800_w.w <- stan(file = "R/FitPhoto800.stan",
                   data = list(A=dat800 %>% filter(Treat=="w.w") %>% pull(Photo),
                               Tk=dat800 %>% filter(Treat=="w.w") %>% pull(Tk),
                               N = dat800 %>% filter(Treat=="w.w") %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=4000, thin=1, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.99))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(sphoto800_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(sphoto800_w.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(sphoto800_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

sphoto800_w.c <- stan(file = "R/FitPhoto800.stan",
                   data=list(A=dat800 %>% filter(Treat=="w.c") %>% pull(Photo), 
                             Tk=dat800 %>% filter(Treat=="w.c") %>% pull(Tk),
                             N=dat800 %>% filter(Treat=="w.c") %>% dim() %>% .[1]),
                   warmup = 2000, save_dso = T,
                   iter=4000, thin=1, chains=3, verbose=T,
                   cores=3,
                   control=list(adapt_delta=0.99))
print(sphoto800_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(sphoto800_w.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(sphoto800_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()
#*******************************************************************************

# Plot Stan Fit w.w and w.c -----------------------------------------------
f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}


p800_1 <- sphoto800_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  scale_y_continuous(position = 'right', limits = c(0,40))+
  labs(y=bquote(paste('P'[800],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  # ylim(0,40)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# p800_1 <- p800_1 + annotate("text", x = 36.5, y = c(37,33), hjust=0,size=4,
#                     label = c("Treatment", 
#                               "paste(Treatment %->% \" Control\")"),parse=T)+
#   annotate("point", x = 35.5, y = c(37,33), 
#            colour = c('red',"blue"), 
#            size = 5, shape=c(16,1))
tmp800_1 <- (dat800 %>% filter(Treat=="w.w")) %>% select(Tk,Photo)
p800_1 <- p800_1 + geom_point(data=tmp800_1, aes(Tk-273.15,Photo), size=4, alpha=0.5,
                      col='red', show.legend = F, inherit.aes = F)
p800_1_topt <- geom_vline(aes(xintercept=(summary(sphoto800_w.w)$summary["Topt","50%"]-273.15)), 
                      col='red')
p800_1_topt_2.5 <- geom_vline(aes(xintercept=summary(sphoto800_w.w)$summary["Topt","2.5%"]-273.15), 
                          col='red',lty=3)
p800_1_topt_97.5 <- geom_vline(aes(xintercept=summary(sphoto800_w.w)$summary["Topt","97.5%"]-273.15), 
                           col='red',lty=3)
p800_1_topt_10 <- geom_vline(aes(
  xintercept=summary(sphoto800_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p800_1_topt_90 <- geom_vline(aes(
  xintercept=summary(sphoto800_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p800_2_lines <- sphoto800_w.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp800_2 <- (dat800 %>% filter(Treat=="w.c")) %>% select(Tk,Photo)
p800_2_points <- geom_point(data=tmp800_2, aes(Tk-273.15,Photo), size=5, shape=1,
                        col='blue', show.legend = F, inherit.aes = F)
p800_2_topt <- geom_vline(aes(xintercept=summary(sphoto800_w.c)$summary["Topt","50%"]-273.15), 
                      col='blue')
p800_2_topt_2.5 <- geom_vline(aes(xintercept=summary(sphoto800_w.c)$summary["Topt","2.5%"]-273.15), 
                          col='blue',lty=3)
p800_2_topt_97.5 <- geom_vline(aes(xintercept=summary(sphoto800_w.c)$summary["Topt","97.5%"]-273.15), 
                           col='blue',lty=3)
p800_2_topt_10 <- geom_vline(aes(
  xintercept=summary(sphoto800_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p800_2_topt_90 <- geom_vline(aes(
  xintercept=summary(sphoto800_w.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

p800_w.w_w.c <- p800_1+p800_2_lines+p800_2_points+p800_1_topt+p800_2_topt+p800_1_topt_10+p800_1_topt_90+p800_2_topt_10+p800_2_topt_90
ggsave("figures/p800_w.w_w.c.png",plot=p800_w.w_w.c, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-



# Plot Stan Fit c.c and w.w -----------------------------------------------
f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}


p800_1 <- sphoto800_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(Tk-273.15, Photo, group=lp))+
  geom_line(alpha=0.00975,col='red')+
  theme_bw()+
  scale_y_continuous(position = 'right', limits = c(0,40))+
  labs(y=bquote(paste('P'[800],'(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  # ylim(0,40)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# p800_1 <- p800_1 + annotate("text", x = 36.5, y = c(37,33), hjust=0,size=4,
#                     label = c("Treatment", 
#                               "paste(Treatment %->% \" Control\")"),parse=T)+
#   annotate("point", x = 35.5, y = c(37,33), 
#            colour = c('red',"blue"), 
#            size = 5, shape=c(16,1))
tmp800_1 <- (dat800 %>% filter(Treat=="w.w")) %>% select(Tk,Photo)
p800_1 <- p800_1 + geom_point(data=tmp800_1, aes(Tk-273.15,Photo), size=4, alpha=0.5,
                              col='red', show.legend = F, inherit.aes = F)
p800_1_topt <- geom_vline(aes(xintercept=(summary(sphoto800_w.w)$summary["Topt","50%"]-273.15)), 
                          col='red')
p800_1_topt_2.5 <- geom_vline(aes(xintercept=summary(sphoto800_w.w)$summary["Topt","2.5%"]-273.15), 
                              col='red',lty=3)
p800_1_topt_97.5 <- geom_vline(aes(xintercept=summary(sphoto800_w.w)$summary["Topt","97.5%"]-273.15), 
                               col='red',lty=3)
p800_1_topt_10 <- geom_vline(aes(
  xintercept=summary(sphoto800_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='red',lty=3)
p800_1_topt_90 <- geom_vline(aes(
  xintercept=summary(sphoto800_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='red',lty=3)


p800_2_lines <- sphoto800_c.c %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  as_tibble() %>% 
  rename(lp=lp__) %>%
  expand(nesting(lp,kopt,Ha,Topt),Tk=seq(295,318,0.1)) %>% 
  mutate(Photo=pmap(., f_photo)) %>% 
  unnest() %>% 
  ungroup() %>% 
  geom_line(data=., aes(Tk-273.15,Photo,group=lp),alpha=0.01,col='blue')
tmp800_2 <- (dat800 %>% filter(Treat=="c.c")) %>% select(Tk,Photo)
p800_2_points <- geom_point(data=tmp800_2, aes(Tk-273.15,Photo), size=4, 
                            col='blue', show.legend = F, inherit.aes = F)
p800_2_topt <- geom_vline(aes(xintercept=summary(sphoto800_c.c)$summary["Topt","50%"]-273.15), 
                          col='blue')
p800_2_topt_2.5 <- geom_vline(aes(xintercept=summary(sphoto800_c.c)$summary["Topt","2.5%"]-273.15), 
                              col='blue',lty=3)
p800_2_topt_97.5 <- geom_vline(aes(xintercept=summary(sphoto800_c.c)$summary["Topt","97.5%"]-273.15), 
                               col='blue',lty=3)
p800_2_topt_10 <- geom_vline(aes(
  xintercept=summary(sphoto800_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","10%"]-273.15), 
  col='blue',lty=3)
p800_2_topt_90 <- geom_vline(aes(
  xintercept=summary(sphoto800_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Topt","90%"]-273.15), 
  col='blue',lty=3)

p800_c.c_w.w <- p800_1+p800_2_lines+p800_2_points+p800_1_topt+p800_2_topt+p800_1_topt_10+p800_1_topt_90+p800_2_topt_10+p800_2_topt_90
ggsave("figures/p800_c.c_w.w_.png",plot=p800_c.c_w.w, width = 120, height=100, units='mm')
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
