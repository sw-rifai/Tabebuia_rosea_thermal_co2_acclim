library(tidyverse); 
library(rstan)

dat_rd <- read_csv("data/Rdark_by_Plant.csv")
names(dat_rd) <- tolower(names(dat_rd))
# dat_rd %>% 
#   group_by(treat) %>% 
#   mutate(mod = lm(rdark~tleaf, data=.))
#   mutate(r_30 = predict(lm(rdark~tleaf, data=.), 
#                       newdata=data.frame(tleaf=30))) 

#*****************************************************************************
# Rdark                 ----------------
#*****************************************************************************
set.seed(3)
rd_c.c <- stan(file = "Stan/Fit_Q10.stan",
               data = list(rdark=dat_rd %>% filter(treat=="c.c") %>% pull(rdark),
                           tleaf=dat_rd %>% filter(treat=="c.c") %>% pull(tleaf),
                           N = dat_rd %>% filter(treat=="c.c") %>% dim() %>% .[1]),
               warmup = 4000, save_dso = T,
               iter=8000, thin=2, chains=3, verbose=T,
               cores=3,
               control=list(adapt_delta=0.95))
# # ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(rd_c.c, digits=3, pars=c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_trace(as.array(rd_c.c), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_dens(as.array(rd_c.c), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))+bayesplot::theme_default()

rd_c.w <- stan(file = "Stan/Fit_Q10.stan",
               data = list(rdark=dat_rd %>% filter(treat=="c.w") %>% pull(rdark),
                           tleaf=dat_rd %>% filter(treat=="c.w") %>% pull(tleaf),
                           N = dat_rd %>% filter(treat=="c.w") %>% dim() %>% .[1]),
               warmup = 4000, save_dso = T,
               iter=8000, thin=2, chains=3, verbose=T,
               cores=3,
               control=list(adapt_delta=0.95))
# # ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(rd_c.w, digits=3, pars=c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_trace(as.array(rd_c.w), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_dens(as.array(rd_c.w), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))+bayesplot::theme_default()

rd_w.c <- stan(file = "Stan/Fit_Q10.stan",
               data = list(rdark=dat_rd %>% filter(treat=="w.c") %>% pull(rdark),
                           tleaf=dat_rd %>% filter(treat=="w.c") %>% pull(tleaf),
                           N = dat_rd %>% filter(treat=="w.c") %>% dim() %>% .[1]),
               warmup = 4000, save_dso = T,
               iter=8000, thin=2, chains=3, verbose=T,
               cores=3,
               control=list(adapt_delta=0.95))
# # ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(rd_w.c, digits=3, pars=c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_trace(as.array(rd_w.c), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_dens(as.array(rd_w.c), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))+bayesplot::theme_default()

rd_w.w <- stan(file = "Stan/Fit_Q10.stan",
               data = list(rdark=dat_rd %>% filter(treat=="w.w") %>% pull(rdark),
                           tleaf=dat_rd %>% filter(treat=="w.w") %>% pull(tleaf),
                           N = dat_rd %>% filter(treat=="w.w") %>% dim() %>% .[1]),
               warmup = 4000, save_dso = T,
               iter=8000, thin=2, chains=3, verbose=T,
               cores=3,
               control=list(adapt_delta=0.95))
# # ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(rd_w.w, digits=3, pars=c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_trace(as.array(rd_w.w), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_dens(as.array(rd_w.w), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))+bayesplot::theme_default()



res_rd <- bind_rows(
  summary(rd_c.c)$summary %>% 
    as_tibble() %>% 
    mutate(param = c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30","lp__"), 
           Treat='c.c'), 
  summary(rd_c.w)$summary %>% 
    as_tibble() %>% 
    mutate(param = c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30","lp__"), 
           Treat='c.w'), 
  summary(rd_w.c)$summary %>% 
    as_tibble() %>% 
    mutate(param = c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30","lp__"), 
           Treat='w.c'), 
  summary(rd_w.w)$summary %>% 
    as_tibble() %>% 
    mutate(param = c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30","lp__"), 
           Treat='w.w')) 

# c.c, c.w, w.c, w.w
# 2.7, 1.9, 1.4, 1.1

res_rd %>% 
  # filter(!param %in% c('lp__','sigma')) %>%
  filter(param %in% c('Q10')) %>%
  ggplot(data=., aes(x=Treat, y=`50%`))+
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`),col='darkolivegreen3',lwd=2)+
  geom_linerange(aes(ymin = `25%`, ymax = `75%`),col='darkolivegreen4',lwd=3)+
  geom_point(col='black')+
  theme_bw()+
  facet_wrap(~param, 
             # labeller = 'label_both',
             scales = 'free',
             ncol = 2, as.table = T)

write_csv(res_rd, path = paste0("outputs/params_rdark30_Q10_",Sys.Date(),".csv"))


  # rdark[n] ~ normal(alpha + beta*tleaf[n], sigma);
  # rdark[n] ~ normal(rdark_30*pow(Q10,((0.1*tleaf[n]-3.0))), sigma_Q10);
fn <- function(rdark_30, Q10, tleaf){
  pow <- function(a,b){return(a**b)}
  out <- rdark_30*pow(Q10,((0.1*tleaf-3.0)))
  return(out)
}

dat_rd %>% filter(treat=="c.c") %>% 
  plot(rdark~tleaf, data=.,col='blue',pch=20, xlim=c(25,45))
curve(fn(
   Q10 = summary(rd_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Q10","50%"], 
   rdark_30 = summary(rd_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["rdark_30","50%"], 
   tleaf=x), add=T,col='blue')
curve(fn(
  Q10 = summary(rd_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Q10","97.5%"], 
  rdark_30 = summary(rd_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["rdark_30","97.5%"], 
  tleaf=x), add=T, lty=3,col='blue')
curve(fn(
  Q10 = summary(rd_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Q10","2.5%"], 
  rdark_30 = summary(rd_c.c, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["rdark_30","2.5%"], 
  tleaf=x), add=T, lty=3,col='blue')

dat_rd %>% filter(treat=="w.w") %>% 
  points(rdark~tleaf, data=.,col='red',pch=20)
curve(fn(
  Q10 = summary(rd_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Q10","50%"], 
  rdark_30 = summary(rd_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["rdark_30","50%"], 
  tleaf=x), add=T,col='red')
curve(fn(
  Q10 = summary(rd_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Q10","97.5%"], 
  rdark_30 = summary(rd_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["rdark_30","97.5%"], 
  tleaf=x), add=T, lty=3,col='red')
curve(fn(
  Q10 = summary(rd_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["Q10","2.5%"], 
  rdark_30 = summary(rd_w.w, probs = c(0.025,0.1,0.5,0.9,0.975))$summary["rdark_30","2.5%"], 
  tleaf=x), add=T, lty=3,col='red')




# Functions ---------------------------------------------------------------
fn <- function(rdark_30, Q10, tleaf, ...){
  pow <- function(a,b){return(a**b)}
  out <- rdark_30*pow(Q10,((0.1*tleaf-3.0)))
  return(out)
}


#########################################################################################
# FOR THE TOP ROW
#-########################################################################################
#
# Defs for all figs
tmin <- 27.5; tmax <- 41.5

# Plot Stan Fit Rd c.c and w.w -----------------------------------------------
tmp_dat1 <- (dat_rd %>% filter(treat=="c.c")) %>% select(tleaf,rdark)
tmp_dat2 <- (dat_rd %>% filter(treat=="w.w")) %>% select(tleaf,rdark)

p_rd_w.w_c.c <- rd_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  rename(lp=lp__) %>%
  expand(nesting(lp,Q10,rdark_30),tleaf=seq(tmin,tmax,1)) %>% 
  mutate(Rd=pmap(., fn)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(tleaf, Rd, group=lp))+
  geom_line(alpha= 0.00975*1,
            col='red')+
  geom_line(data= rd_c.c %>%  
                as.data.frame() %>% 
                sample_n(1000) %>%
                rename(lp=lp__) %>%
                expand(nesting(lp,Q10,rdark_30),tleaf=seq(tmin,tmax,1)) %>% 
                mutate(Rd=pmap(., fn)) %>% 
                unnest() %>% 
                ungroup(),
              aes(tleaf, Rd, group=lp),
      alpha= 0.00975,
              col='blue') + 
  theme_bw()+
  labs(y=bquote(paste('Dark respiration ','(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(breaks = seq(0,3.5,by=0.5), expand=expand_scale(0,0), limits = c(0,3.5))+
  scale_x_continuous(breaks = seq(ceiling(tmin),floor(tmax),by=2), expand = expand_scale(0,0))+
  labs(x=NULL)+
  annotate("point", x = tmin+1, y = c(3.3,3*0.9), 
           colour = c("blue","red"), 
           size = 5, shape=c(16,16), alpha=0.5)+
  annotate("text", x = tmin+2, y = c(3.3,3*0.9), hjust=0,size=4,
                    label = c("Control", 
                              "Treatment"),parse=T)+
    geom_point(data=tmp_dat1, aes(tleaf, rdark), size=5, alpha=0.5,
                      col='blue', show.legend = F, inherit.aes = F) +
  geom_point(data=tmp_dat2, aes(tleaf, rdark), size=5, alpha=0.5,
             col='red', show.legend = F, inherit.aes = F)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_blank(), #element_text(size=18), 
        axis.text.y.left = element_text(size=17), 
        axis.title.y.left = element_text(size=20),
        plot.margin = margin(15,15,15,15)
  )
p_rd_w.w_c.c    
  
ggsave(paste0("figures/ratio_jmax_vcmax/rd_c.c_w.w_",Sys.Date(),".png"),plot=p_rd_w.w_c.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/ratio_jmax_vcmax/rd_c.c_w.w_",Sys.Date(),".pdf"),plot=p_rd_w.w_c.c, width = 120, height=100, units='mm')
# geom_hline(aes(yintercept=(summary(rd_w.w)$summary["rdark_30","97.5%"])), col='red',lty=3)+
# geom_hline(aes(yintercept=(summary(rd_w.w)$summary["rdark_30","50%"])), col='red')+
# geom_hline(aes(yintercept=(summary(rd_w.w)$summary["rdark_30","2.5%"])), col='red',lty=3)+
# geom_hline(aes(yintercept=(summary(rd_c.c)$summary["rdark_30","97.5%"])), col='blue',lty=3)+
# geom_hline(aes(yintercept=(summary(rd_c.c)$summary["rdark_30","50%"])), col='blue')+
# geom_hline(aes(yintercept=(summary(rd_c.c)$summary["rdark_30","2.5%"])), col='blue',lty=3)
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-



# Plot Stan Fit Rd c.c and c.w -----------------------------------------------
# p1 <- 
tmp_dat1 <- (dat_rd %>% filter(treat=="c.c")) %>% select(tleaf,rdark)
tmp_dat2 <- (dat_rd %>% filter(treat=="c.w")) %>% select(tleaf,rdark)

p_rd_c.c_c.w <- rd_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  rename(lp=lp__) %>%
  expand(nesting(lp,Q10,rdark_30),tleaf=seq(tmin,tmax,1)) %>% 
  mutate(Rd=pmap(., fn)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(tleaf, Rd, group=lp))+
  geom_line(alpha= 0.00975*1,
            col='red')+
  geom_line(data= rd_c.c %>%  
              as.data.frame() %>% 
              sample_n(1000) %>%
              rename(lp=lp__) %>%
              expand(nesting(lp,Q10,rdark_30),tleaf=seq(tmin,tmax,1)) %>% 
              mutate(Rd=pmap(., fn)) %>% 
              unnest() %>% 
              ungroup(),
            aes(tleaf, Rd, group=lp),
            alpha= 0.00975,
            col='blue') + 
  theme_bw()+
  labs(y=bquote(paste('Dark respiration ','(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(breaks = seq(0,3.5,by=0.5), expand=expand_scale(0,0), limits = c(0,3.5))+
  scale_x_continuous(breaks = seq(ceiling(tmin),floor(tmax),by=2), expand = expand_scale(0,0))+
  labs(x=NULL)+
  annotate("point", x = tmin+1, y = c(3.3,3*0.9), 
           colour = c("blue","red"), 
           size = 5, shape=c(16,1))+
  annotate("text", x = tmin+2, y = c(3.3,3*0.9), hjust=0,size=4,
           label = c("Control", 
                     "paste(Control %->% \" Treatment\")"),parse=T)+
  geom_point(data=tmp_dat1, aes(tleaf, rdark), size=5, alpha=0.5,
             col='blue', show.legend = F, inherit.aes = F) +
  geom_point(data=tmp_dat2, aes(tleaf, rdark), size=5, shape=1,
             col='red', show.legend = F, inherit.aes = F)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_blank(), #element_text(size=18), 
        axis.text.y.left = element_text(size=17), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(15,15,15,15)
  )
p_rd_c.c_c.w    

ggsave(paste0("figures/ratio_jmax_vcmax/rd_c.c_c.w_",Sys.Date(),".png"),plot=p_rd_c.c_c.w, width = 120, height=100, units='mm')
ggsave(paste0("figures/ratio_jmax_vcmax/rd_c.c_c.w_",Sys.Date(),".pdf"),plot=p_rd_c.c_c.w, width = 120, height=100, units='mm')
# geom_hline(aes(yintercept=(summary(rd_w.w)$summary["rdark_30","97.5%"])), col='red',lty=3)+
# geom_hline(aes(yintercept=(summary(rd_w.w)$summary["rdark_30","50%"])), col='red')+
# geom_hline(aes(yintercept=(summary(rd_w.w)$summary["rdark_30","2.5%"])), col='red',lty=3)+
# geom_hline(aes(yintercept=(summary(rd_c.c)$summary["rdark_30","97.5%"])), col='blue',lty=3)+
# geom_hline(aes(yintercept=(summary(rd_c.c)$summary["rdark_30","50%"])), col='blue')+
# geom_hline(aes(yintercept=(summary(rd_c.c)$summary["rdark_30","2.5%"])), col='blue',lty=3)
#_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

# Plot Stan Fit Rd c.c and c.w -----------------------------------------------
tmp_dat1 <- (dat_rd %>% filter(treat=="w.w")) %>% select(tleaf,rdark)
tmp_dat2 <- (dat_rd %>% filter(treat=="w.c")) %>% select(tleaf,rdark)

p_rd_w.w_w.c <- rd_w.w %>%  
  as.data.frame() %>% 
  sample_n(1000) %>%
  rename(lp=lp__) %>%
  expand(nesting(lp,Q10,rdark_30),tleaf=seq(tmin,tmax,1)) %>% 
  mutate(Rd=pmap(., fn)) %>% 
  unnest() %>% 
  ungroup() %>% 
  ggplot(data=., aes(tleaf, Rd, group=lp))+
  geom_line(alpha= 0.00975*1,
            col='red')+
  geom_line(data= rd_w.c %>%  
              as.data.frame() %>% 
              sample_n(1000) %>%
              rename(lp=lp__) %>%
              expand(nesting(lp,Q10,rdark_30),tleaf=seq(tmin,tmax,1)) %>% 
              mutate(Rd=pmap(., fn)) %>% 
              unnest() %>% 
              ungroup(),
            aes(tleaf, Rd, group=lp),
            alpha= 0.00975,
            col='blue') + 
  theme_bw()+
  labs(y=bquote(paste('Dark respiration ','(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(breaks = seq(0,3.5,by=0.5), expand=expand_scale(0,0), limits = c(0,3.5))+
  scale_x_continuous(breaks = seq(ceiling(tmin),floor(tmax),by=2), expand = expand_scale(0,0))+
  annotate("point", x = tmin+1, y = c(3.3,3*0.9), 
           colour = c("red","blue"), 
           size = 5, shape=c(16,1))+
  annotate("text", x = tmin+2, y = c(3.3,3*0.9), hjust=0,size=4,
           label = c("Treatment", 
                     "paste(Treatment %->% \" Control\")"),parse=T)+
  geom_point(data=tmp_dat1, aes(tleaf, rdark), size=5, alpha=0.5,
             col='red', show.legend = F, inherit.aes = F) +
  geom_point(data=tmp_dat2, aes(tleaf, rdark), size=5, shape=1,
             col='blue', show.legend = F, inherit.aes = F)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=17), 
        axis.title.x.bottom = element_text(size=20), 
        axis.text.y.left = element_text(size=17), 
        axis.title.y.left = element_text(size=18),
        plot.margin = margin(15,15,15,15)
  )
p_rd_w.w_w.c    

ggsave(paste0("figures/ratio_jmax_vcmax/rd_w.w_w.c_",Sys.Date(),".png"),plot=p_rd_w.w_w.c, width = 120, height=100, units='mm')
ggsave(paste0("figures/ratio_jmax_vcmax/rd_w.w_w.c_",Sys.Date(),".pdf"),plot=p_rd_w.w_w.c, width = 120, height=100, units='mm')




# Join images -------------------------------------------------------------
library(magick)
i1 <- image_read("figures/ratio_jmax_vcmax/rd_c.c_w.w_2019-06-10.png")
i2 <- image_read("figures/ratio_jmax_vcmax/rd_c.c_c.w_2019-06-10.png")
i3 <- image_read("figures/ratio_jmax_vcmax/rd_w.w_w.c_2019-06-10.png")
# i1 <- image_annotate(i1, "(a)", size = 100, color = "black",
#                      location = "+320+20")
# i2 <- image_annotate(i2, "(b)", size = 100, color = "black",
#                      location = "+85+20")
# i3 <- image_annotate(i3, "(c)", size = 100, color = "black",
#                      location = "+320+20")
# i1 <- image_border(i1, color='white',geometry = '1x5')
# i2 <- image_border(i2, color='white',geometry = '1x5')
# i3 <- image_border(i3, color='white',geometry = '1x5')

img_join <- image_append(c(i1,i2,i3), stack = T)
# img_top <- c(i1,i2)
# img_mid <- c(i3,i4)
# img_bottom <- c(i5,i6)
# ijoin_top <- image_append(img_top)
# ijoin_mid <- image_append(img_mid)
# ijoin_bottom <- image_append(img_bottom)
# ijoin_P400P800 <- image_append(c(ijoin_top,ijoin_mid, ijoin_bottom), stack=T)

image_write(img_join, path = paste0("figures/ratio_jmax_vcmax/Joined_Dark_Resp_",Sys.Date(),".png"),
            format = "png")

