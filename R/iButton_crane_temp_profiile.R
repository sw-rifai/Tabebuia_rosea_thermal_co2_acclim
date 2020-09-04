library(tidyverse); library(hms); library(viridis)
# start time: 6:30
# end time: 8:50
dat <- read_csv("data/ibutton_profile_crane.csv", skip = 1, col_names = F)

dtime <- seq(parse_date_time("20110101 063000", orders = "ymd HMS"), 
    parse_date_time("20110108 085000", orders = "ymd HMS"),
    by="10 min")
heights <- c(40,38,36,34,30,28,26,24,20,18,16,12,8,6,4,2)
dat[,1] <- heights

colnames(dat) <- c("height",dtime)
dat <- dat %>% gather(dateTime, temp, -height)
dat <- dat %>% mutate(dateTime=as.POSIXct(as.numeric(dateTime),origin="1970-01-01"))
dat <- dat %>% mutate(hour=hour(dateTime))

lfit1 <- dat %>% lm(temp~scale(height)*scale(hour), data=.)
lfit2 <- dat %>% lm(temp~scale(height)+scale(hour), data=.)
pfit2 <- dat %>% lm(temp~scale(height)+scale(hour) + I(hour**2), data=.)
pfit3 <- dat %>% lm(temp~scale(height)+scale(hour) + I(hour**2) + I(hour**3), data=.)
efit <- dat %>% lm(temp~log(height)+log(hour+1), data=.)
jfit <- dat %>% lm(temp ~ height+log(height)+hour+I(hour**2), data=.)
jfit2 <- dat %>% lm(temp ~ height+log(height)+hour+I(hour**2)+I(hour**3), data=.)

summary(lfit2)
summary(pfit3)
summary(efit)
summary(jfit2)

bbmle::BICtab(lfit1,lfit2,pfit2,pfit3,efit,jfit,jfit2)

curve(predict(pfit3, newdata=data.frame(height=30,hour=x)), 1,24,col="red", xlab="Hour",ylab="Temp. [C]")
curve(predict(pfit3, newdata=data.frame(height=2,hour=x)), 1,24,add=T)
curve(predict(pfit3, newdata=data.frame(height=10,hour=x)), 1,24,add=T,col="blue")
curve(predict(pfit3, newdata=data.frame(height=20,hour=x)), 1,24,add=T,col="darkgreen")
legend("topleft",col=c("red","darkgreen","blue","black"), legend=c("30 m","20 m", "10 m", "2 m"),lwd=c(1,1,1,1))

curve(predict(pfit3, newdata=data.frame(height=x,hour=6)), 2,40,col="black", xlab="Height",ylab="Temp. [C]",ylim=c(20,30))
curve(predict(pfit3, newdata=data.frame(height=x,hour=12)), 2,40,add=T,col="red")
curve(predict(pfit3, newdata=data.frame(height=x,hour=18)), 2,40,add=T,col="darkgreen")
curve(predict(pfit3, newdata=data.frame(height=x,hour=24)), 2,40,add=T,col="blue")
legend("topleft",col=c("black","red","dark green","blue"), legend=c("6:00","12:00","18:00","24:00"),lwd=c(1,1,1,1), horiz = T)

curve(predict(jfit2, newdata=data.frame(height=40,hour=x)), 1,24,col="red", xlab="Hour",ylab="Temp. [C]")
curve(predict(jfit2, newdata=data.frame(height=2,hour=x)), 1,24,add=T)
curve(predict(jfit2, newdata=data.frame(height=10,hour=x)), 1,24,add=T,col="blue")
curve(predict(jfit2, newdata=data.frame(height=25,hour=x)), 1,24,add=T,col="darkgreen")
legend("topleft",col=c("red","darkgreen","blue","black"), legend=c("40 m","25 m", "10 m", "2 m"),lwd=c(1,1,1,1))

curve(predict(jfit2, newdata=data.frame(height=x,hour=6)), 2,40,col="black", xlab="Height",ylab="Temp. [C]",ylim=c(20,30))
curve(predict(jfit2, newdata=data.frame(height=x,hour=12)), 2,40,add=T,col="red")
curve(predict(jfit2, newdata=data.frame(height=x,hour=18)), 2,40,add=T,col="darkgreen")
curve(predict(jfit2, newdata=data.frame(height=x,hour=24)), 2,40,add=T,col="blue")
legend("topleft",col=c("black","red","dark green","blue"), legend=c("6:00","12:00","18:00","24:00"),lwd=c(1,1,1,1), horiz = T)





dat %>% ggplot(data=., aes(dateTime,temp,color=height))+geom_path()+
  scale_color_viridis()

dat %>% ggplot(data=., aes(height,temp,color=as.factor(hour)))+
  # geom_point()+
  geom_smooth(se=F)+
  scale_color_viridis(discrete = T,option="A",end=0.9)+
  theme_bw()

dat %>%
  filter(hour<=12) %>% 
  ggplot(data=., aes(height,temp,color=as.factor(hour)))+
  # geom_point()+
  geom_smooth(se=F)+
  scale_color_viridis(discrete = T,option="A",end=0.9)+
  theme_bw()

dat %>%
  filter(hour>12) %>% 
  ggplot(data=., aes(height,temp,color=as.factor(hour)))+
  # geom_point()+
  geom_smooth(se=F)+
  scale_color_viridis(discrete = T,option="A",end=0.9)+
  theme_bw()
