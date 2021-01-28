library(tidybayes); 
### Plot Vcmax CO2 400 & 800
library(tidyverse); library(bayesplot); library(rstan)
options(mc.cores = parallel::detectCores()-1); 
rstan_options(auto_write = TRUE);
set.seed(1); 

dat_all<-read_csv("data/params.all.28032018.csv")
dat_all <- dat_all %>% 
  mutate(Tk = Tleaf+273.15)

# treatments: (w.w, w.w, w.c, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
# treat<-as.character(unique(params.all$Treat))

# Fit Vcmax ---------------------------------------------------------------------
vc400_c.c <- stan(file = "Stan/FitTopt_sim_Vcmax.stan",
                  data = list(Vcmax=dat_all %>% filter(Treat=='c.c') %>% pull(Vcmax),
                              Tk=dat_all %>% filter(Treat=='c.c') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='c.c') %>% dim() %>% .[1]),
                  warmup = 4000, save_dso = T,
                  iter=8000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(vc400_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma","Vcmax_Topt"))
mcmc_trace(as.array(vc400_c.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(vc400_c.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


vc400_w.w <- stan(file = "Stan/FitTopt_sim_Vcmax.stan",
                  data = list(Vcmax=dat_all %>% filter(Treat=='w.w') %>% pull(Vcmax),
                              Tk=dat_all %>% filter(Treat=='w.w') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='w.w') %>% dim() %>% .[1]),
                  warmup = 4000, save_dso = T,
                  iter=8000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(vc400_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(vc400_w.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(vc400_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

vc400_c.w <- stan(file = "Stan/FitTopt_sim_Vcmax.stan",
                  data = list(Vcmax=dat_all %>% filter(Treat=='c.w') %>% pull(Vcmax),
                              Tk=dat_all %>% filter(Treat=='c.w') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='c.w') %>% dim() %>% .[1]),
                  warmup = 4000, save_dso = T,
                  iter=8000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(vc400_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(vc400_c.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(vc400_c.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

vc400_w.c <- stan(file = "Stan/FitTopt_sim_Vcmax.stan",
                  data = list(Vcmax=dat_all %>% filter(Treat=='w.c') %>% pull(Vcmax),
                              Tk=dat_all %>% filter(Treat=='w.c') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='w.c') %>% dim() %>% .[1]),
                  warmup = 4000, save_dso = T,
                  iter=8000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(vc400_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(vc400_w.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(vc400_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


# Fit Jmax ---------------------------------------------------------------------
jm400_c.c <- stan(file = "Stan/FitTopt_Jmax.stan",
                  data = list(Jmax=dat_all %>% filter(Treat=='c.c') %>% pull(Jmax),
                              Tk=dat_all %>% filter(Treat=='c.c') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='c.c') %>% dim() %>% .[1]),
                  warmup = 4000, save_dso = T,
                  iter=8000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(jm400_c.c, digits=3, pars=c("kopt","Ha","Topt","sigma","Jmax_Topt"))
mcmc_trace(as.array(jm400_c.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(jm400_c.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


jm400_w.w <- stan(file = "Stan/FitTopt_Jmax.stan",
                  data = list(Jmax=dat_all %>% filter(Treat=='w.w') %>% pull(Jmax),
                              Tk=dat_all %>% filter(Treat=='w.w') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='w.w') %>% dim() %>% .[1]),
                  warmup = 4000, save_dso = T,
                  iter=8000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(jm400_w.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(jm400_w.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(jm400_w.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

jm400_c.w <- stan(file = "Stan/FitTopt_Jmax.stan",
                  data = list(Jmax=dat_all %>% filter(Treat=='c.w') %>% pull(Jmax),
                              Tk=dat_all %>% filter(Treat=='c.w') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='c.w') %>% dim() %>% .[1]),
                  warmup = 4000, save_dso = T,
                  iter=8000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(jm400_c.w, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(jm400_c.w), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(jm400_c.w), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()

jm400_w.c <- stan(file = "Stan/FitTopt_Jmax.stan",
                  data = list(Jmax=dat_all %>% filter(Treat=='w.c') %>% pull(Jmax),
                              Tk=dat_all %>% filter(Treat=='w.c') %>% pull(Tk),
                              N = dat_all %>% filter(Treat=='w.c') %>% dim() %>% .[1]),
                  warmup = 4000, save_dso = T,
                  iter=8000, thin=2, chains=3, verbose=T,
                  cores=3,
                  control=list(adapt_delta=0.95))
# ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(jm400_w.c, digits=3, pars=c("kopt","Ha","Topt","sigma"))
mcmc_trace(as.array(jm400_w.c), c("sigma","Topt","Ha","kopt"))
mcmc_dens(as.array(jm400_w.c), c("sigma","Topt","Ha","kopt"))+bayesplot::theme_default()


f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}


dat_all<-read_csv("data/params.all.28032018.csv")
dat_all <- dat_all %>% 
  mutate(Tk = Tleaf+273.15)

# treatments: (w.w, w.w, w.c, w.c, for control plants, warmed plants, 
# control plants transferred to warmed, and warmed plants transferred to controll)
# treat<-as.character(unique(params.all$Treat))



# Functions ---------------------------------------------------------------
f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}



#########################################################################################
# FOR THE TOP ROW
#-########################################################################################
# Plot Jmax/Vcmax c.c & w.w -----------------------------------------------
tmp_c.c <- full_join({vc400_c.c %>% 
    as.data.frame() %>% 
    sample_n(1000) %>% 
    mode_hdi(.width = 0.5) %>% 
    expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
    mutate(Photo=pmap(., f_photo)) %>% 
    unnest() %>% 
    ungroup() %>% 
    rename(Vcmax=Photo) %>% 
    select(Tk,Vcmax)},
    {jm400_c.c %>% 
        as.data.frame() %>% 
        mode_hdi(.width = 0.5) %>% 
        expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
        mutate(Photo=pmap(., f_photo)) %>% 
        unnest() %>% 
        ungroup() %>% 
        rename(Jmax=Photo) %>% 
        select(Tk,Jmax)}) %>% 
  mutate(treat='c.c')
tmp_w.w <- full_join({vc400_w.w %>% 
    as.data.frame() %>% 
    mode_hdi(.width = 0.5) %>% 
    expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
    mutate(Photo=pmap(., f_photo)) %>% 
    unnest() %>% 
    ungroup() %>% 
    rename(Vcmax=Photo) %>% 
    select(Tk,Vcmax)},
    {jm400_w.w %>% 
        as.data.frame() %>% 
        mode_hdi(.width = 0.5) %>% 
        expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
        mutate(Photo=pmap(., f_photo)) %>% 
        unnest() %>% 
        ungroup() %>% 
        rename(Jmax=Photo) %>% 
        select(Tk,Jmax)}) %>% 
  mutate(treat='w.w')

p_jmvc_cc_ww <- bind_rows(tmp_c.c, tmp_w.w) %>% 
  ggplot(data=., aes(Tk-273.15, Jmax/Vcmax, color=treat))+
  geom_line(alpha=1)+
  theme_bw()+
  geom_point(data=(dat_all %>% select(Treat,Tk,Vcmax,Jmax) %>% 
                     rename(treat=Treat) %>% filter(treat %in% c("c.c","w.w"))), 
             aes(Tk-273.15,Jmax/Vcmax), 
             # pch=20, 
             size=5, alpha=0.5,
             show.legend = F)+
  scale_color_manual(values=c("c.c"='blue',"w.w"='red'))+
  labs(y=bquote(paste('J'['max']/'V'['Cmax'])), 
       xlab=NA
       # x=bquote(paste('Leaf Temperature ('^degree,'C)'))
  )+
  scale_y_continuous(limits=c(0,2), expand = c(0,0))+
  annotate("point", x = 25, y = c(0.4,0.3), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,16))+
  annotate("text", x = 25+1, y = c(0.4, 0.3), hjust=0,size=5,
           # label = c("paste(Treatment %->% \" Treatment\")", 
           #           "paste(Treatment %->% \" Control\")"),
           label = c("Control", 
                     "Treatment"),
           parse=T)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = 'none')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_blank(), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        axis.title.x.bottom = element_blank(), 
        plot.margin = margin(10,10,10,10))+
  xlab(label = NULL)
p_jmvc_cc_ww

ggsave(plot = p_jmvc_cc_ww, 
       paste0("figures/ratio_jmax_vcmax/jmax_vmax_ratio_c.c_w.w_",Sys.Date(),".png"),
       width = 120, height=100, units='mm')
ggsave(plot = p_jmvc_cc_ww, 
       paste0("figures/ratio_jmax_vcmax/jmax_vmax_ratio_c.c_w.w_",Sys.Date(),".pdf"),
       width = 120, height=100, units='mm')



##########################################################################################
# FOR THE MIDDLE ROW
##########################################################################################
#
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
# Plot Jmax/Vcmax c.c & c.w -----------------------------------------------
tmp_c.c <- full_join({vc400_c.c %>% 
    as.data.frame() %>% 
    mode_hdi(.width = 0.5) %>% 
    expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
    mutate(Photo=pmap(., f_photo)) %>% 
    unnest() %>% 
    ungroup() %>% 
    rename(Vcmax=Photo) %>% 
    select(Tk,Vcmax)},
    {jm400_c.c %>% 
        as.data.frame() %>% 
        mode_hdi(.width = 0.5) %>% 
        expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
        mutate(Photo=pmap(., f_photo)) %>% 
        unnest() %>% 
        ungroup() %>% 
        rename(Jmax=Photo) %>% 
        select(Tk,Jmax)}) %>% 
  mutate(treat='c.c')
tmp_c.w <- full_join({vc400_c.w %>% 
    as.data.frame() %>% 
    mode_hdi(.width = 0.5) %>% 
    expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
    mutate(Photo=pmap(., f_photo)) %>% 
    unnest() %>% 
    ungroup() %>% 
    rename(Vcmax=Photo) %>% 
    select(Tk,Vcmax)},
    {jm400_c.w %>% 
        as.data.frame() %>% 
        mode_hdi(.width = 0.5) %>% 
        expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
        mutate(Photo=pmap(., f_photo)) %>% 
        unnest() %>% 
        ungroup() %>% 
        rename(Jmax=Photo) %>% 
        select(Tk,Jmax)}) %>% 
  mutate(treat='c.w')

p_jmvc_cc_cw <- bind_rows(tmp_c.c, tmp_c.w) %>% 
  ggplot(data=., aes(Tk-273.15, Jmax/Vcmax, color=treat))+
  geom_line(alpha=1)+
  theme_bw()+
  geom_point(data=(dat_all %>% select(Treat,Tk,Vcmax,Jmax) %>% 
                     rename(treat=Treat) %>% filter(treat %in% c("c.c","c.w"))), 
             aes(Tk-273.15,Jmax/Vcmax,shape=treat), 
             # shape=1, 
             size=5, alpha=0.5,
             show.legend = F)+
  scale_shape_manual(values=c("c.c"=16,"c.w"=1))+
  scale_color_manual(values=c("c.c"='blue',"c.w"='red'))+
  labs(y=bquote(paste('J'['max']/'V'['Cmax'])), 
       xlab=NULL
       # x=bquote(paste('Leaf Temperature ('^degree,'C)'))
  )+
  scale_y_continuous(limits=c(0,2), expand = c(0,0))+
  annotate("text", x = 25+1, y = c(0.4, 0.3), hjust=0,size=5,
           label = c("paste(Control %->% \" Control\")",
                     "paste(Control %->% \" Treatment\")"),
           # label = c("Control", 
           #           "Treatment"),
           parse=T)+
  annotate("point", x = 25, y = c(0.4,0.3), 
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
        axis.title.x.bottom = element_blank(), 
        plot.margin = margin(10,10,10,10))+
  xlab(label = NULL)
p_jmvc_cc_cw
ggsave(plot = p_jmvc_cc_cw, 
       paste0("figures/ratio_jmax_vcmax/jmax_vmax_ratio_c.c_c.w_",Sys.Date(),".png"),
       width = 120, height=100, units='mm')
ggsave(plot = p_jmvc_cc_cw, 
       paste0("figures/ratio_jmax_vcmax/jmax_vmax_ratio_c.c_c.w_",Sys.Date(),".pdf"),
       width = 120, height=100, units='mm')

##########################################################################################
# FOR THE BOTTOM ROW
##########################################################################################
# Plot Jmax/Vcmax w.w & w.c -----------------------------------------------
tmp_w.w <- full_join({vc400_w.w %>% 
    as.data.frame() %>% 
    mode_hdi(.width = 0.5) %>% 
    expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
    mutate(Photo=pmap(., f_photo)) %>% 
    unnest() %>% 
    ungroup() %>% 
    rename(Vcmax=Photo) %>% 
    select(Tk,Vcmax)},
    {jm400_w.w %>% 
        as.data.frame() %>% 
        mode_hdi(.width = 0.5) %>% 
        expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
        mutate(Photo=pmap(., f_photo)) %>% 
        unnest() %>% 
        ungroup() %>% 
        rename(Jmax=Photo) %>% 
        select(Tk,Jmax)}) %>% 
  mutate(treat='w.w')
tmp_w.c <- full_join({vc400_w.c %>% 
    as.data.frame() %>% 
    mode_hdi(.width = 0.5) %>% 
    expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
    mutate(Photo=pmap(., f_photo)) %>% 
    unnest() %>% 
    ungroup() %>% 
    rename(Vcmax=Photo) %>% 
    select(Tk,Vcmax)},
    {jm400_w.c %>% 
        as.data.frame() %>% 
        mode_hdi(.width = 0.5) %>% 
        expand(nesting(kopt,Ha,Topt),Tk=seq(295,321,0.1)) %>% 
        mutate(Photo=pmap(., f_photo)) %>% 
        unnest() %>% 
        ungroup() %>% 
        rename(Jmax=Photo) %>% 
        select(Tk,Jmax)}) %>% 
  mutate(treat='w.c')

p_jmvc_ww_wc <- bind_rows(tmp_w.w, tmp_w.c) %>% 
  ggplot(data=., aes(Tk-273.15, Jmax/Vcmax, color=treat))+
  geom_line(alpha=1)+
  theme_bw()+
  geom_point(data=(dat_all %>% select(Treat,Tk,Vcmax,Jmax) %>% 
                     rename(treat=Treat) %>% filter(treat %in% c("w.w","w.c"))), 
             aes(Tk-273.15,Jmax/Vcmax,shape=treat), 
             # shape=1, 
             size=5, alpha=0.5,
             show.legend = F)+
  scale_shape_manual(values=c("w.w"=16,"w.c"=1))+
  scale_color_manual(values=c("w.w"='red',"w.c"='blue'))+
  labs(y=bquote(paste('J'['max']/'V'['Cmax'])),
       x=bquote(paste('Leaf Temperature ('^degree,'C)'))
  )+
  scale_y_continuous(limits=c(0,2), expand = c(0,0))+
  annotate("text", x = 25+1, y = c(0.4, 0.3), hjust=0,size=4,
           label = c("paste(Treatment %->% \" Treatment\")",
                     "paste(Treatment %->% \" Control\")"),
           # label = c("Control", 
           #           "Treatment"),
           parse=T)+
  annotate("point", x = 25, y = c(0.4,0.3), 
           colour = c('red',"blue"), 
           size = 5, alpha=0.5, shape=c(16,1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = 'none')+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        axis.title.x.bottom = element_text(size=20), 
        plot.margin = margin(10,10,10,10))
p_jmvc_ww_wc

ggsave(plot = p_jmvc_ww_wc, 
       filename = paste0("figures/ratio_jmax_vcmax/jmax_vmax_ratio_w.w_w.c_",Sys.Date(),".png"),
       width = 120, height=100, units='mm')
ggsave(plot = p_jmvc_ww_wc, 
       filename = paste0("figures/ratio_jmax_vcmax/jmax_vmax_ratio_w.w_w.c_",Sys.Date(),".pdf"),
       width = 120, height=100, units='mm')
# end *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-



library(magick)
i1 <- image_read("figures/ratio_jmax_vcmax/jmax_vmax_ratio_c.c_w.w_2019-06-10.png")
i2 <- image_read("figures/ratio_jmax_vcmax/jmax_vmax_ratio_c.c_c.w_2019-06-10.png")
i3 <- image_read("figures/ratio_jmax_vcmax/jmax_vmax_ratio_w.w_w.c_2019-06-10.png")
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
            path = paste0("figures/ratio_jmax_vcmax/Joined_Ratio_Jmax_Vcmax_",Sys.Date(),".png"), 
            format = "png")



# # Join Figures ------------------------------------------------------------
# library(magick)
# i1 <- image_read("figures/vcmax400_c.c_w.w_2019-06-10.png")
# i2 <- image_read("figures/jmax400_c.c_w.w_2019-06-10.png")
# i3 <- image_read("figures/vcmax400_c.c_c.w_2019-06-10.png")
# i4 <- image_read("figures/jmax400_c.c_c.w_2019-06-10.png")
# i5 <- image_read("figures/vcmax400_w.w_w.c_2019-06-10.png")
# i6 <- image_read("figures/jmax400_w.w_w.c_2019-06-10.png")
# i1 <- image_annotate(i1, "(a)", size = 90, color = "black",
#                      location = "+300+40")
# i2 <- image_annotate(i2, "(b)", size = 90, color = "black",
#                      location = "+50+40")
# i3 <- image_annotate(i3, "(c)", size = 90, color = "black",
#                      location = "+300+40")
# i4 <- image_annotate(i4, "(d)", size = 90, color = "black",
#                      location = "+50+40")
# i5 <- image_annotate(i5, "(e)", size = 90, color = "black",
#                      location = "+300+40")
# i6 <- image_annotate(i6, "(f)", size = 90, color = "black",
#                      location = "+50+40")
# 
# img_top <- c(i1,i2)
# img_mid <- c(i3,i4)
# img_bottom <- c(i5,i6)
# ijoin_top <- image_append(img_top)
# ijoin_mid <- image_append(img_mid)
# ijoin_bottom <- image_append(img_bottom)
# ijoin_VcmaxJmax <- image_append(c(ijoin_top,ijoin_mid, ijoin_bottom), stack=T)
# 
# image_write(ijoin_VcmaxJmax, 
#             path = paste0("figures/Joined_Vcmax_Jmax_",Sys.Date(),".png"), format = "png")
