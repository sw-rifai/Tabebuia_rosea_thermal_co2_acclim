library(tidyverse); library(bayesplot);
library(brms); library(nls.multstart)
library(patchwork)
options(mc.cores = parallel::detectCores()-1)
set.seed(321)

dat_rd <- read_csv("data/Rdark_by_Plant.csv")
dat_rd <- dat_rd %>% mutate(Tk = Tleaf+273.15) %>% 
  mutate(cc=ifelse(Treat=='c.c',1,0), 
         cw=ifelse(Treat=='c.w',1,0), 
         wc=ifelse(Treat=='w.c',1,0), 
         ww=ifelse(Treat=='w.w',1,0))
names(dat_rd) <- tolower(names(dat_rd))



lm(log(rdark)~cw+wc+ww + tleaf + tleaf*wc + tleaf*wc + tleaf*ww, 
   data=dat_rd) %>% summary
parr_rd <- brm(log(rdark)~cw+wc+ww + tleaf + tleaf*cw + tleaf*wc + tleaf*ww, 
               data=dat_rd, family=gaussian(), 
               control=list(adapt_delta=0.99))
bayes_R2(parr_rd)
pp_check(parr_rd, nsamples = 30)
# plot(parr_rd)

# Q10 = exp(10*b)
# ln(rdark) = a + b*tleaf

fixef(parr_rd, summary=F) %>% 
  as_tibble() %>% 
  rename(b_tleaf = tleaf) %>% 
  mutate(Q10_cc = exp(10*b_tleaf)) %>% 
  pull(Q10_cc) %>% 
  summary


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

fn_q10 <- function(rdark_30, Q10, tleaf, ...){
  pow <- function(a,b){return(a**b)}
  out <- rdark_30*pow(Q10,((0.1*tleaf-3.0)))
  return(out)
}

#########################################################################################
# FOR THE TOP ROW of Rdark
#-########################################################################################
# Defs for all figs
tmin <- 27.5; tmax <- 41.5


prd1 <- parr_rd %>% 
  as.data.frame() %>% 
  sample_n(1000) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>%
  rename(a=b, acw=cw, awc=wc,aww=ww) %>% 
  rename(b=tleaf, bwc=`wc:tleaf`, bcw=`cw:tleaf`,bww=`ww:tleaf`) %>% 
  expand(nesting(lp,
                 a,acw,awc,aww,
                 b,bcw,bwc,bww), 
         cw=c(0,1), 
         wc=c(0,1),
         ww=c(0,1), 
         leaf_temps=seq(tmin,tmax,length.out=50)) %>% 
  mutate(Rd = exp((a + acw*cw + awc*wc + aww*ww)+(b + bcw*cw + bwc*wc + bww*ww)*leaf_temps)) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","w.w")) %>% 
  ggplot(data=.,aes(leaf_temps,Rd,color=Treat,group=lp_Treat))+
  geom_line(alpha=0.00975)+
  geom_point(data=dat_rd %>% 
               rename(Treat=treat, 
                      Rd=rdark,
                      leaf_temps=tleaf) %>% 
               filter(Treat%in%c("c.c","w.w")), 
             aes(leaf_temps, Rd,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  theme_bw()+
  scale_color_manual("",
                     values=c("c.c"='blue', "w.w"='red'), 
                     labels=c("c.c"="Control",
                              "c.w"="paste(Control %->% \" Treatment\")",
                              "w.c"="paste(Treatment %->% \" Control\")", 
                              "w.w"="Treatment"))+
  labs(y=bquote(paste('Dark respiration ','(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(breaks = seq(0,3.5,by=0.5), expand=expansion(0,0), limits = c(0,3.5))+
  scale_x_continuous(breaks = seq(ceiling(tmin),floor(tmax),by=2), expand = expansion(0,0))+
  labs(x=NULL)+
  # annotate("point", x = tmin+1, y = c(3.3,3*0.9), 
  #          colour = c("blue","red"), 
  #          size = 5, shape=c(16,16), alpha=0.5)+
  # annotate("text", x = tmin+2, y = c(3.3,3*0.9), hjust=0,size=4,
  #          label = c("Control", 
  #                    "Treatment"),parse=T)+
  theme(legend.position = c(0.05,0.95), 
        legend.justification = c(0.05,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))
  

#************************************************************************
# FOR THE MIDDLE ROW of Rdark -------------------
#************************************************************************
prd2 <- parr_rd %>% 
  as.data.frame() %>% 
  sample_n(1000) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>%
  rename(a=b, acw=cw, awc=wc,aww=ww) %>% 
  rename(b=tleaf, bwc=`wc:tleaf`, bcw=`cw:tleaf`,bww=`ww:tleaf`) %>% 
  expand(nesting(lp,
                 a,acw,awc,aww,
                 b,bcw,bwc,bww), 
         cw=c(0,1), 
         wc=c(0,1),
         ww=c(0,1), 
         leaf_temps=seq(tmin,tmax,length.out=50)) %>% 
  mutate(Rd = exp((a + acw*cw + awc*wc + aww*ww)+(b + bcw*cw + bwc*wc + bww*ww)*leaf_temps)) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("c.c","c.w")) %>% 
  ggplot(data=.,aes(leaf_temps,Rd,color=Treat,group=lp_Treat))+
  geom_line(alpha=0.00975)+
  geom_point(data=dat_rd %>% 
               rename(Treat=treat, 
                      Rd=rdark,
                      leaf_temps=tleaf) %>% 
               filter(Treat%in%c("c.c")), 
             aes(leaf_temps, Rd,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat_rd %>% 
               rename(Treat=treat, 
                      Rd=rdark,
                      leaf_temps=tleaf) %>% 
               filter(Treat%in%c("c.w")), 
             aes(leaf_temps, Rd,color=Treat),
             size=4,alpha=0.5,pch=1,
             inherit.aes = F)+
  theme_bw()+
  scale_color_manual("", 
                     values=c("blue","red"), 
                     breaks=c("c.c","c.w"),
                     labels=c("Control",
                              sprintf('Control\u2192Treatment')
                     ),
                     guide=guide_legend(override.aes = list(shape=c(20,1), 
                                                            size=c(6,3)))
  )+
  labs(y=bquote(paste('Dark respiration ','(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(breaks = seq(0,3.5,by=0.5), expand=expansion(0,0), limits = c(0,3.5))+
  scale_x_continuous(breaks = seq(ceiling(tmin),floor(tmax),by=2), expand = expansion(0,0))+
  labs(x=NULL)+
  theme(legend.position = c(0.05,0.95), 
        legend.justification = c(0.05,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))

#************************************************************************
# FOR THE Bottom ROW of Rdark -------------------
#************************************************************************

prd3 <- parr_rd %>% 
  as.data.frame() %>% 
  sample_n(1000) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>%
  rename(a=b, acw=cw, awc=wc,aww=ww) %>% 
  rename(b=tleaf, bwc=`wc:tleaf`, bcw=`cw:tleaf`,bww=`ww:tleaf`) %>% 
  expand(nesting(lp,
                 a,acw,awc,aww,
                 b,bcw,bwc,bww), 
         cw=c(0,1), 
         wc=c(0,1),
         ww=c(0,1), 
         leaf_temps=seq(tmin,tmax,length.out=50)) %>% 
  mutate(Rd = exp((a + acw*cw + awc*wc + aww*ww)+(b + bcw*cw + bwc*wc + bww*ww)*leaf_temps)) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  filter(Treat %in% c("w.c","w.w")) %>% 
  ggplot(data=.,aes(leaf_temps,Rd,color=Treat,group=lp_Treat))+
  geom_line(alpha=0.00975)+
  geom_point(data=dat_rd %>% 
               rename(Treat=treat, 
                      Rd=rdark,
                      leaf_temps=tleaf) %>% 
               filter(Treat%in%c("w.w")), 
             aes(leaf_temps, Rd,color=Treat),
             size=4,alpha=0.5,
             inherit.aes = F)+
  geom_point(data=dat_rd %>% 
               rename(Treat=treat, 
                      Rd=rdark,
                      leaf_temps=tleaf) %>% 
               filter(Treat%in%c("w.c")), 
             aes(leaf_temps, Rd,color=Treat),
             size=4,alpha=0.5,pch=1,
             inherit.aes = F)+
  theme_bw()+
  scale_color_manual("", 
                     values=c("red","blue"), 
                     breaks=c("w.w","w.c"),
                     labels=c("Treatment",
                              sprintf('Treatment\u2192Control')
                     ),
                     guide=guide_legend(override.aes = list(shape=c(20,1), 
                                                            size=c(6,3)))
  )+
  labs(y=bquote(paste('Dark respiration ','(',mu,'mol m'^'-2','s'^'-1',')')),
       x=bquote(paste('Leaf Temperature ('^degree,'C)')))+
  scale_y_continuous(breaks = seq(0,3.5,by=0.5), expand=expansion(0,0), limits = c(0,3.5))+
  scale_x_continuous(breaks = seq(ceiling(tmin),floor(tmax),by=2), expand = expansion(0,0))+
  labs(x=NULL)+
  theme(legend.position = c(0.05,0.95), 
        legend.justification = c(0.05,0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        plot.margin = margin(l=20))



ggsave(prd1/prd2/prd3,
       filename = paste0("figures/Fig7_rdark_brms_",Sys.Date(),".png"),
       width = 1*120, height=3*100, units='mm')

library(data.table)
out_RdQ10 <- parr_rd %>% 
  as.data.frame() %>% 
  # sample_n(1000) %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>%
  rename(a=b, acw=cw, awc=wc,aww=ww) %>% 
  rename(b=tleaf, bwc=`wc:tleaf`, bcw=`cw:tleaf`,bww=`ww:tleaf`) %>% 
  expand(nesting(lp,
                 a,acw,awc,aww,
                 b,bcw,bwc,bww), 
         cw=c(0,1), 
         wc=c(0,1),
         ww=c(0,1) 
         ) %>% 
  mutate(R30 = exp((a + acw*cw + awc*wc + aww*ww)+(b + bcw*cw + bwc*wc + bww*ww)*30), 
         Q10 = exp(10*((b + bcw*cw + bwc*wc + bww*ww)))) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% #pull(Treat) %>% table
  filter(is.na(Treat)==F) %>% 
  mutate(lp_Treat=paste(lp,Treat)) %>% 
  # filter(Treat %in% c("c.c","w.w"))
  select( Treat, R30, Q10) %>%
  as.data.table() %>% 
  .[,as.data.table(posterior_summary(.SD, probs = c(0.05,0.95)),
                   keep.rownames = TRUE), by='Treat'] %>% 
  .[order(rn)]
write_csv(out_RdQ10, file = paste0('outputs/params_rdark30_Q10_',Sys.Date(),".csv"))


# Acclim_setTemp = Rd[c.c]@30 / Rd[w.w]@30
# Acclim_Homeo = Rd[c.c]@30 / Rd[w.w]@40
out_acclim <- parr_rd %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  rename(lp=lp__) %>%
  rename_all(fn) %>%
  rename(a=b, acw=cw, awc=wc,aww=ww) %>% 
  rename(b=tleaf, bwc=`wc:tleaf`, bcw=`cw:tleaf`,bww=`ww:tleaf`) %>% 
  expand(nesting(lp,
                 a,acw,awc,aww,
                 b,bcw,bwc,bww), 
         cw=c(0,1), 
         wc=c(0,1),
         ww=c(0,1) 
  ) %>% 
  mutate(Treat=case_when(cw==1 & wc==0 & ww==0 ~"c.w", 
                         wc==1 & cw==0 & ww==0 ~"w.c",
                         ww==1 & wc==0 & cw==0 ~"w.w",
                         ww==0 & wc==0 & cw==0 ~"c.c")) %>% 
  filter(is.na(Treat)==FALSE) %>%
  mutate(R26.5 = exp((a + acw*cw + awc*wc + aww*ww)+(b + bcw*cw + bwc*wc + bww*ww)*26.5)) %>% 
  mutate(R29.5 = exp((a + acw*cw + awc*wc + aww*ww)+(b + bcw*cw + bwc*wc + bww*ww)*29.5)) %>% 
  mutate(R30.5 = exp((a + acw*cw + awc*wc + aww*ww)+(b + bcw*cw + bwc*wc + bww*ww)*30.5)) %>% 
  mutate(R35.5 = exp((a + acw*cw + awc*wc + aww*ww)+(b + bcw*cw + bwc*wc + bww*ww)*35.5)) %>% 
  select(Treat,R26.5,R29.5,R30.5,R35.5) %>%
  as.data.table() %>%
  .[,as.data.table(posterior_summary(.SD, probs = c(0.05,0.95)),
                   keep.rownames = TRUE), by='Treat'] %>%
  .[order(rn)]
write_csv(out_acclim, file = paste0('outputs/params_rdark_acclim_',Sys.Date(),".csv"))


