library(tidyverse)
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
library(rstan)
set.seed(3)
rd_c.c <- stan(file = "Stan/Fit_Q10.stan",
                     data = list(rdark=dat_rd %>% filter(treat=="c.c") %>% pull(rdark),
                                 tleaf=dat_rd %>% filter(treat=="c.c") %>% pull(tleaf),
                                 N = dat_rd %>% filter(treat=="c.c") %>% dim() %>% .[1]),
                     warmup = 2000, save_dso = T,
                     iter=4000, thin=2, chains=3, verbose=T,
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
               warmup = 2000, save_dso = T,
               iter=4000, thin=2, chains=3, verbose=T,
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
               warmup = 2000, save_dso = T,
               iter=4000, thin=2, chains=3, verbose=T,
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
               warmup = 2000, save_dso = T,
               iter=4000, thin=2, chains=3, verbose=T,
               cores=3,
               control=list(adapt_delta=0.95))
# # ! NEED TO ENSURE THE MODEL CONVERGED !!!
print(rd_w.w, digits=3, pars=c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_trace(as.array(rd_w.w), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))
mcmc_dens(as.array(rd_w.w), c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30"))+bayesplot::theme_default()



res_rd <- bind_rows(
  summary(rd_c.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30","lp__"), 
           Treat='c.c'), 
  summary(rd_c.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30","lp__"), 
           Treat='c.w'), 
  summary(rd_w.c, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30","lp__"), 
           Treat='w.c'), 
  summary(rd_w.w, probs=c(0.1,0.25,0.5,0.75,0.9))$summary %>% 
    as_tibble() %>% 
    mutate(param = c("alpha","beta","Q10","sigma","sigma_Q10","rdark_30","lp__"), 
           Treat='w.w')) 
  
# c.c, c.w, w.c, w.w
# 2.7, 1.9, 1.4, 1.1

res_rd %>% 
  # filter(!param %in% c('lp__','sigma')) %>%
  filter(param %in% c('Q10')) %>%
  ggplot(data=., aes(x=Treat, y=`50%`))+
  geom_linerange(aes(ymin = `10%`, ymax = `90%`),col='darkolivegreen3',lwd=2)+
  geom_linerange(aes(ymin = `25%`, ymax = `75%`),col='darkolivegreen4',lwd=3)+
  geom_point(col='black')+
  theme_bw()+
  facet_wrap(~param, 
             # labeller = 'label_both',
             scales = 'free',
             ncol = 2, as.table = T)

write_csv(res_rd, path = paste0("outputs/params_rdark30_Q10_",Sys.Date(),".csv"))



# lm(rdark~log(tleaf), data=dat_rd %>% filter(treat=='c.c')) %>% summary
# 
# 
# 
# r_30 <- (tmp %>% filter(treat=='w.w') %>% lm(rdark~tleaf, data=.) %>% 
#            predict(., newdata=data.frame(tleaf=30)))
# nls(rdark ~ r_30*Q10^((tleaf-30)/10), 
#     data=tmp %>% filter(treat=='w.w'), 
#     start=list(Q10=1), 
#     algorithm='default') %>% summary
# nls(rdark ~ r_30*Q10^((tleaf-30)/10), 
#     data=tmp %>% filter(treat=='c.c'), 
#     start=list(Q10=1), 
#     algorithm='plinear') %>% summary
# nls(rdark ~ r_30*Q10^((tleaf-30)/10), 
#     data=tmp %>% filter(treat=='c.c'), 
#     start=list(Q10=1), 
#     algorithm='port') %>% summary
# 
# 
# dat_all %>% glimpse
# 
# dat_fci %>% 
#   ggplot(data=., aes(Tleaf, Rday.25, color=Treat))+
#   geom_point()+
#   geom_smooth(method='lm')
# 
# tmp <- read_csv("data/P400.2904.csv")
# 
# 
# x <- runif(100, 1,3)
# y <- rnorm(100, mean=0.5*exp(x**1.1), sd=sqrt(x))
# plot(y~x)
# 
# LL <- function(b0, b1,sigma){
#   -sum(dnorm(x, mean=b0*exp(x**b1), sd=sigma, log = T))
# }
# LL(1,1,1)
# bbmle::mle2(LL, method="L-BFGS-B", 
#             start=list(b0=0.1, b1=1, sigma=1), 
#             data=list(x=x))
# 
# 
# rgamma(1000, shape = 0.01, rate = 0.01) %>% hist(100)
# 
# dat_rd <- read_csv("data/Rdark_by_Plant.csv")
# names(dat_rd) <- tolower(names(dat_rd))
# tmp <- dat_rd %>% 
#   group_by(treat,plant) %>% 
#   arrange(tleaf) %>% 
#   mutate(tdiff = tleaf-lag(tleaf), 
#          rdiff = rdark-lag(rdark), 
#          r_lag = lag(rdark)) %>%
#   ungroup() %>% 
#   filter(is.na(rdiff)==F) %>% 
#   mutate(tdiff10=tdiff/10)
# 
# 
# 
# dat_rd <- read_csv("data/Rdark_by_Plant.csv")
# names(dat_rd) <- tolower(names(dat_rd))
# # c.c, c.w, w.c, w.w
# # 2.7, 1.9, 1.4, 1.1
# r_30 <- (tmp %>% filter(treat=='w.w') %>% lm(rdark~tleaf, data=.) %>% 
#   predict(., newdata=data.frame(tleaf=30)))
# nls(rdark ~ r_30*Q10^((tleaf-30)/10), 
#     data=tmp %>% filter(treat=='w.w'), 
#     start=list(Q10=1), 
#     algorithm='default') %>% summary
# nls(rdark ~ r_30*Q10^((tleaf-30)/10), 
#     data=tmp %>% filter(treat=='c.c'), 
#     start=list(Q10=1), 
#     algorithm='plinear') %>% summary
# nls(rdark ~ r_30*Q10^((tleaf-30)/10), 
#     data=tmp %>% filter(treat=='c.c'), 
#     start=list(Q10=1), 
#     algorithm='port') %>% summary
# 
# 
# fit <- lm(rdark~exp(I(0.1*tleaf-3)), data=tmp)
# plot(rdark~tleaf, data=tmp %>% filter(treat=='c.c')); abline(fit)
# plot(c(NA),c(NA), ylim=c(0,10), xlim=c(30,40))
# abline(fit)
# 
# r_30
# fn <- function(r_30, Q10, x){ (r_30*Q10**((x-30)/10))}
# curve(fn(1, 3.0474, x), 30,40)
# curve(fn(1, 2.65, x), 30,40,col='red', add=T)
# points(rdark~tleaf, data=dat_rd %>% filter(treat=='c.c'))
# curve(fn(1, 2.5, x), 30,45,col='red',add=T)
# 
# 
# 
# 
# r_30 <- (tmp %>% filter(treat=='c.w') %>% lm(rdark~tleaf, data=.) %>% 
#            predict(., newdata=data.frame(tleaf=30)))
# nls(rdark ~ 0.925*Q10^((tleaf-30)/10), 
#     data=tmp %>% filter(treat=='c.w'), 
#     start=list(Q10=1.5), 
#     algorithm = 'port') %>% summary
# 
# r_30 <- (tmp %>% filter(treat=='w.c') %>% lm(rdark~tleaf, data=.) %>% 
#            predict(., newdata=data.frame(tleaf=30)))
# nls(rdark ~ r_30*Q10^((tleaf-30)/10), 
#     data=tmp %>% filter(treat=='c.c'), 
#     start=list(Q10=1.5), 
#     algorithm = 'port') %>% summary
# 
# r_30 <- (tmp %>% filter(treat=='c.c') %>% lm(rdark~tleaf, data=.) %>% 
#            predict(., newdata=data.frame(tleaf=30)))
# nls(rdark ~ 0.925*Q10^((tleaf-30)/10), 
#     data=tmp %>% filter(treat=='c.c'), 
#     start=list(Q10=1.5), 
#     algorithm = 'plinear') %>% summary
# 
# 
# 
# 
# lm(rdark ~ r_lag**(tdiff/10), data=tmp %>% filter(treat=='c.c')) %>% summary
# lm(rdark ~ r_lag*exp(tdiff/10), data=tmp %>% filter(treat=='c.w'))
# lm(rdark ~ r_lag*exp(tdiff/10), data=tmp %>% filter(treat=='w.c'))
# lm(rdark ~ r_lag*exp(tdiff/10), data=tmp %>% filter(treat=='w.w'))
# 
# lm(rdark ~ r_lag**(tdiff/10), data=tmp %>% filter(treat=='c.c')) %>% summary
# 
# 
# nls(Q10 ~ (rdark/r_lag)^(tdiff10), 
#     data=tmp %>% filter(treat=='c.c'), 
#     start=list(Q10=1.5), 
#     algorithm = 'plinear') %>% summary
# 
# with(tmp %>% filter(treat=='c.c'), (rdark/r_lag)^(tdiff10)) %>% summary
# 
# nls(rdark ~ (r_lag)*(Q10)^(tdiff10*0.1), 
#     data=tmp %>% filter(treat=='c.c'), 
#     start=list(Q10=0.5)) %>% summary
# nls(rdark ~ (r_lag*Q10)**(tdiff/10), data=tmp %>% filter(treat=='c.w'), 
#     start=list(Q10=1)) %>% summary
# nls(rdark ~ (r_lag*Q10)**(tdiff/10), data=tmp %>% filter(treat=='w.c'), 
#     start=list(Q10=1)) %>% summary
# nls(rdark ~ (r_lag*Q10)**(tdiff/10), data=tmp %>% filter(treat=='w.w'), 
#     start=list(Q10=1)) %>% summary
# 
# fn <- function(r1, Q10, x){ (r1*Q10**((x-r1)/10))}
# curve(fn(1, 2, x), 30,45, ylim=c(0,50))
# curve(fn(1, 2.5, x), 30,45,col='red',add=T)
# 
# nls(rdark ~ alpha*exp(0.1*tleaf*beta), 
#     data=dat_rd %>% filter(treat=='c.c'), 
#     start=list(alpha=1, beta=0.001)) %>% 
#   coef %>% exp
# 
# nls(rdark ~ alpha*exp(0.1*tleaf*beta), 
#     data=dat_rd %>% filter(treat=='c.w'), 
#     start=list(alpha=1, beta=0.001)) %>% 
#   coef %>% exp
# 
# nls(rdark ~ alpha*exp(0.1*tleaf*beta), 
#     data=dat_rd %>% filter(treat=='w.c'), 
#     start=list(alpha=1, beta=0.001)) %>% 
#   coef %>% exp
# 
# nls(rdark ~ alpha*exp(0.1*tleaf*beta), 
#     data=dat_rd %>% filter(treat=='w.w'), 
#     start=list(alpha=1, beta=0.001)) %>% 
#   coef
# 
# nls(rdark ~ alpha*exp(0.1*tleaf*beta), 
#     data=dat_rd %>% filter(treat=='w.w'), 
#     start=list(alpha=1, beta=0.001)) %>% 
#   coef %>% exp
# 
# 
# fn <- function(r1, Q10, x){ (r1*Q10**((x-r1)/10))}
# curve(fn(1, 1.8, x), 30,45)
# plot(rdark~tleaf, data=dat_rd %>% filter(treat=='c.c'))
# curve(fn(1, 2.5, x), 30,45,col='red',add=T)
# 
# 
# 
# # Q10=exp(10*beta)
# lm(log(rdark)~tleaf, data=dat_rd %>% filter(treat=="c.c")) %>% 
#   coef %>% .[2]*10 %>% as.numeric %>% .[1]*exp(1)
# lm(log(rdark)~tleaf, data=dat_rd %>% filter(treat=="c.w")) %>% 
#   coef %>% .[2]*10 %>% as.numeric %>% .[1]*exp(1)
# lm(log(rdark)~tleaf, data=dat_rd %>% filter(treat=="w.c")) %>% 
#   coef %>% .[2]*10 %>% as.numeric %>% .[1]*exp(1)
# lm(log(rdark)~tleaf, data=dat_rd %>% filter(treat=="w.w")) %>% 
#   coef %>% .[2]*10 %>% as.numeric
# # c.c, c.w, w.c, w.w
# # 2.7, 1.9, 1.4, 1.1
# 
# 
# Q10 <- function(data, temp, ...){ 
#   var <- data[,sens.of] 
#   temp <- data [,temp] 
#   exponent <- nls(var ~ a*exp(1)^(b*temp),...) 
#   exponent <- coef(exponent)[2] 
#   Q10 <- exp(10*exponent) 
#   names(Q10) <- "Q10" 
#   Q10 <- cbind(exponent, Q10) 
#   row.names(Q10) <- NULL 
#   return(Q10) 
# }
# 
# Q10(tmp %>% filter(treat=='c.c'), temp=tmp %>% filter(treat=='c.c') %>% pull(tleaf))
# dat_rd %>% ggplot(data=., aes(rdark, tleaf, group=plant))+
#   geom_point()+
#   geom_line()+
#   # geom_smooth(se=F, method='lm')+
#   facet_wrap(~treat)
# 
# 
# 
# 
