library(tidybayes); 

f_photo <- function(kopt,Ha,Tk,Topt, ...){
  photo <- kopt * ((200 * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                     (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt)))))))
  return(photo)}

# Plot Jmax/Vcmax c.c & w.w -----------------------------------------------
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

bind_rows(tmp_c.c, tmp_w.w) %>% 
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
       xlab=NA
       # x=bquote(paste('Leaf Temperature ('^degree,'C)'))
       )+
  scale_x_continuous(limits=c(25,46))+
  scale_y_continuous(limits=c(0,2), expand = c(0,0))+
  annotate("text", x = 25, y = c(0.4, 0.3), hjust=0,size=4,
           # label = c("paste(Treatment %->% \" Treatment\")", 
           #           "paste(Treatment %->% \" Control\")"),
           label = c("Control", 
                     "Treatment"),
           parse=T)+
  annotate("point", x = 31, y = c(0.4,0.3), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,16))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = 'none')+
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        axis.title.x.bottom = element_text(size=20))+
  xlab(label = NULL)

ggsave(paste0("figures/jmax_vmax_ratio_c.c_w.w_",Sys.Date(),".png"),
       width = 120, height=100, units='mm')
# end *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

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

bind_rows(tmp_c.c, tmp_c.w) %>% 
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
  scale_y_continuous(limits=c(0,2), expand = c(0,0))+
  annotate("text", x = 25, y = c(0.4, 0.3), hjust=0,size=4,
           label = c("paste(Control %->% \" Control\")",
                     "paste(Control %->% \" Treatment\")"),
           # label = c("Control", 
           #           "Treatment"),
           parse=T)+
  annotate("point", x = 36, y = c(0.4,0.3), 
           colour = c('blue',"red"), 
           size = 5, alpha=0.5, shape=c(16,1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = 'none')+
  scale_x_continuous(breaks = seq(25,45,by=5), expand = expand_scale(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x.bottom = element_text(size=18), 
        axis.text.y.left = element_text(size=18), 
        axis.title.y.left = element_text(size=20), 
        axis.title.x.bottom = element_text(size=20))+
  xlab(label = NULL)

ggsave(paste0("figures/jmax_vmax_ratio_c.c_c.w_",Sys.Date(),".png"),
       width = 120, height=100, units='mm')
# end *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

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

bind_rows(tmp_w.w, tmp_w.c) %>% 
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
  scale_y_continuous(limits=c(0,2), expand = c(0,0))+
  annotate("text", x = 25, y = c(0.4, 0.3), hjust=0,size=4,
           label = c("paste(Treatment %->% \" Treatment\")",
                     "paste(Treatment %->% \" Control\")"),
           # label = c("Control", 
           #           "Treatment"),
           parse=T)+
  annotate("point", x = 37.5, y = c(0.4,0.3), 
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
        axis.title.x.bottom = element_text(size=20))

ggsave(paste0("figures/jmax_vmax_ratio_w.w_w.c_",Sys.Date(),".png"),
       width = 120, height=100, units='mm')
# end *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


library(magick)
i1 <- image_read("figures/jmax_vmax_ratio_c.c_w.w_2019-02-20.png")
i2 <- image_read("figures/jmax_vmax_ratio_c.c_c.w_2019-02-20.png")
i3 <- image_read("figures/jmax_vmax_ratio_w.w_w.c_2019-02-20.png")
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
            path = paste0("figures/Ratio_Jmax_Vcmax_",Sys.Date(),".png"), 
            format = "png")

