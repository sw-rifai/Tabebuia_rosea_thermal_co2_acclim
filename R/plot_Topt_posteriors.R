library(tidyverse); 
library(ggdist); # easier ways to plot the posterior
library(colorspace); # many color palettes
library(patchwork) # Join the plots into a composite
# Load the MCMC fits -------------------------------------------------
parabolic_p400 <- read_rds("outputs/parr_parabolic_p400_2021-01-16.rds")
parabolic_p800 <- read_rds("outputs/parr_parabolic_p800_2021-01-16.rds")
parr_vcmax <- read_rds("outputs/parr_vcmax_2021-01-16.rds")
parr_jmax <- read_rds("outputs/parr_jmax_2021-01-16.rds")


(p1 <- parabolic_p400 %>%
  as_tibble() %>% 
  mutate(Topt_cc = b_Topt_Intercept, # CC is the base
         Topt_cw = b_Topt_Intercept+b_Tcw_Intercept, # Treatment effects added to the base
         Topt_wc = b_Topt_Intercept+b_Twc_Intercept,
         Topt_ww = b_Topt_Intercept+b_Tww_Intercept) %>% 
  select(Topt_cc, Topt_cw, Topt_wc, Topt_ww) %>% # Recasting the DF from wide to long
  gather(key=group, value=value) %>% 
  filter(is.na(group)==F) %>% 
  ggplot(aes(y = group, x = value)) +
  stat_halfeye(aes(fill = stat(cut_cdf_qi(
    cdf, 
    .width = c(0.5, 0.9, 0.99), # This controls the number of fill shades in the posterior
    labels = scales::percent_format()
  ))), 
  .width=c(0.5,0.9,0.99)) + # This controls horizontal line increments
  scale_fill_discrete_sequential(palette = "Heat", nmax = 6, order = 3:6, na.translate=F)+ # many palettes to choose from
  labs(
    title = expression(paste(T[opt]~of~Photosynthesis[400])), 
    # subtitle = "aes(fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95))))",
    fill = "Interval", 
    x=expression(paste(T[opt]~(K)))
  )+
  # expand_limits(x=c(295,310))+ # Way to specifiy the x-limits, but it won't cut into the range of the data
  theme_linedraw())

(p2 <- parabolic_p800 %>%
  as_tibble() %>% 
  mutate(Topt_cc = b_Topt_Intercept, 
         Topt_cw = b_Topt_Intercept+b_Tcw_Intercept, 
         Topt_wc = b_Topt_Intercept+b_Twc_Intercept,
         Topt_ww = b_Topt_Intercept+b_Tww_Intercept) %>% 
  select(Topt_cc, Topt_cw, Topt_wc, Topt_ww) %>% 
  gather(key=group, value=value) %>% 
  filter(is.na(group)==F) %>% 
  ggplot(aes(y = group, x = value)) +
  stat_halfeye(aes(fill = stat(cut_cdf_qi(
    cdf, 
    .width = c(0.5, 0.9, 0.99),
    labels = scales::percent_format()
  )))) +
  scale_fill_discrete_sequential(palette = "Heat", nmax = 6, order = 3:6, na.translate=F)+
  labs(
    title = expression(paste(T[opt]~of~Photosynthesis[800])), 
    # subtitle = "aes(fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95))))",
    fill = "Interval", 
    x=expression(paste(T[opt]~(K)))
  )+
    # expand_limits(x=c(295,310))+
  theme_linedraw())

(p3 <- parr_vcmax %>%
    as_tibble() %>% 
    mutate(Topt_cc = b_Topt_Intercept, 
           Topt_cw = b_Topt_Intercept+b_Tcw_Intercept, 
           Topt_wc = b_Topt_Intercept+b_Twc_Intercept,
           Topt_ww = b_Topt_Intercept+b_Tww_Intercept) %>% 
    select(Topt_cc, Topt_cw, Topt_wc, Topt_ww) %>% 
    gather(key=group, value=value) %>% 
    filter(is.na(group)==F) %>% 
    ggplot(aes(y = group, x = value)) +
    stat_halfeye(aes(fill = stat(cut_cdf_qi(
      cdf, 
      .width = c(0.5, 0.9, 0.99),
      labels = scales::percent_format()
    )))) +
    scale_fill_discrete_sequential(palette = "Heat", nmax = 6, order = 3:6, na.translate=F)+
    labs(
      title = expression(paste(T[opt]~of~V[CMax])), 
      # subtitle = "aes(fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95))))",
      fill = "Interval", 
      x=expression(paste(T[opt]~(K)))
    )+
    # expand_limits(x=c(295,310))+
    theme_linedraw())

(p4 <- parr_jmax %>%
    as_tibble() %>% 
    mutate(Topt_cc = b_Topt_Intercept, 
           Topt_cw = b_Topt_Intercept+b_Tcw_Intercept, 
           Topt_wc = b_Topt_Intercept+b_Twc_Intercept,
           Topt_ww = b_Topt_Intercept+b_Tww_Intercept) %>% 
    select(Topt_cc, Topt_cw, Topt_wc, Topt_ww) %>% 
    gather(key=group, value=value) %>% 
    filter(is.na(group)==F) %>% 
    ggplot(aes(y = group, x = value)) +
    stat_halfeye(aes(fill = stat(cut_cdf_qi(
      cdf, 
      .width = c(0.5, 0.9, 0.99),
      labels = scales::percent_format()
    )))) +
    scale_fill_discrete_sequential(palette = "Heat", nmax = 6, order = 3:6, na.translate=F)+
    labs(
      title = expression(paste(T[opt]~of~J[Max])), 
      # subtitle = "aes(fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95))))",
      fill = "Interval", 
      x=expression(paste(T[opt]~(K)))
    )+
    # expand_limits(x=c(295,310))+
    theme_linedraw())





p_out <- (p1|p2)/(p3|p4)+plot_layout(guides='collect')
ggsave(p_out, 
       filename = paste0('figures/Topt_posteriors_P400_P800_VCMax_JMax_',Sys.Date(),'.png'),
       width=20, height=15, units='cm'
       )
