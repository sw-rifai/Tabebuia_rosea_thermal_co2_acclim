devtools::install_github('padpadpadpad/nls.multstart')
library(nls.multstart) 
library(rstan); 
options(mc.cores = parallel::detectCores())
nobs <- 100
kopt_t <- 170; 
Hd_t <- 105; 
Ha_t <- 55; 
Topt_t <- 301.5; 
Tk <- rnorm(nobs, mean = 307, sd = 4)
eps_ind <- rnorm(nobs, mean=0, sd=10)
Vcmax_t <- kopt_t * ((Hd_t * (2.718282^((Ha_t*(Tk-Topt_t))/(Tk*0.008314*Topt_t)))) / 
                       (Hd_t - (Ha_t*(1-(2.718282^((Hd_t*(Tk-Topt_t))/(Tk*0.008314*Topt_t))))))) + eps_ind 

plot(Vcmax_t~Tk)  
sim_dat <- data.frame(Tk=Tk, Vcmax=Vcmax_t)

fit_med <- nls_multstart(Vcmax ~ kopt * ((Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                                           (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))),
                         data = sim_dat,
                         iter = 1000,
                         start_lower = c(kopt = 160, Hd = 100, Ha = 50, Topt = 300),
                         start_upper = c(kopt = 300, Hd = 900, Ha = 300, Topt = 330),
                         supp_errors = 'Y',
                         na.action = na.omit,
                         #convergence_count = 500,
                         lower = c(kopt = 150, Hd = 50, Ha = 10, Topt = 300))

fit_med %>% summary


file_path <- "R/FitTopt.stan";
lines <- readLines(file_path, encoding="ASCII");
for (n in 1:length(lines)) cat(lines[n],'\n');

m_stan <- stan(model_code = lines, #"R/StanSim/meanOfPreds.stan",
               data = list(Vcmax=Vcmax_t,
                           Tk=Tk,
                           N = nobs),
               iter=5000, thin=5, chains=4, verbose=T, 
               control=list(adapt_delta=0.9))

print(m_stan, digits=3, pars=c("kopt","Hd","Ha","Topt","sigma"))
mcmc_trace(as.array(m_stan), c("sigma","Topt","Hd","Ha","kopt"))
mcmc_dens(as.array(m_stan), c("sigma","Topt","Hd","Ha","kopt"))

post_df <- rstan::extract(m_stan)

r_Vcmax <- function(kopt, Hd, Ha, Tk, Topt){
  Vcmax <- kopt * ( (Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                      (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))    ); 
  return(Vcmax); 
}

df_med <- lapply(post_df,"median") %>% as.data.frame()
df_5 <- lapply(post_df,"quantile",0.05) %>% as.data.frame();
df_95 <- lapply(post_df,"quantile",0.95) %>% as.data.frame();
df_ppd <- tibble(Tk=290:320,Vcmax05=NA,Vcmax50=NA,Vcmax95=NA);

for(i in 1:dim(df_ppd)[1]){
 df_ppd$Vcmax05[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.05, na.rm=T)
 df_ppd$Vcmax50[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.5, na.rm=T)
 df_ppd$Vcmax95[i] <- r_Vcmax(kopt=post_df$kopt, Hd=post_df$Hd, Ha=post_df$Ha, Tk=df_ppd$Tk[i], Topt=post_df$Topt) %>% quantile(., 0.95, na.rm=T)
}

curve(r_Vcmax(post_df$kopt[1], post_df$Hd[1], post_df$Ha[1], Tk, post_df$Topt[1]), 
      xname = "Tk",from = 290,to=320, ylim=c(90,200), col=rgb(0,0,0,alpha=0.1), 
      ylab="Vcmax",xlab="Temperature [K]")
for(i in 1:1000){
  curve(r_Vcmax(post_df$kopt[i], post_df$Hd[i], post_df$Ha[i], Tk, post_df$Topt[i]), 
        xname = "Tk",from = 290,to=320,add=T,col=rgb(0,0,0,alpha=0.1))
}
lines(Vcmax50~Tk, data=df_ppd, type="l",col="green",lwd=3)
lines(Vcmax05~Tk, data=df_ppd, type="l",col="green",lwd=2, lty=3)
lines(Vcmax95~Tk, data=df_ppd, type="l",col="green",lwd=2, lty=3)



# Vcmax ~ normal(kopt * ((Hd * (2.718282^((Ha*(Tk[n]-Topt))/(Tk*0.008314*Topt)))) / 
#                       (Hd - (Ha*(1-(2.718282^((Hd*(Tk[n]-Topt))/(Tk[n]*0.008314*Topt)))))))

(Hd - (Ha*(1-(2.718282^((Hd*(Tk[n]-Topt))/(Tk[n]*0.008314*Topt))))))
model_code <-
  '
functions {
real calc_Vcmax(real kopt, real Hd, real Ha, real Tk, real Topt) {
real a; 
real b; 
real d; 
real out; 

// b = (Tk*0.008314*Topt);
// a = pow(2.718282, ((Ha*(Tk-Topt)))/b);
// d = Hd - Ha*(1-(pow(2.718282, ((Hd*(Tk-Topt))/(Tk*0.008314*Topt)))));
// ((Hd*(Tk-Topt))/(Tk*0.008314*Topt))
// d = pow(2.718282,);
out = kopt * ( (Hd * pow(2.718282,((Ha*(Tk-Topt))/(Tk*0.008314*Topt))) ) / 
(Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))    ); 

// out = kopt * (Hd*(a/b))/d;

return out;
}
}
model {}
'
expose_stan_functions(stanc(model_code = model_code))
calc_Vcmax(1, 200, 400, 304, 305)

r_Vcmax <- function(kopt, Hd, Ha, Tk, Topt){
  Vcmax <- kopt * ( (Hd * (2.718282^((Ha*(Tk-Topt))/(Tk*0.008314*Topt)))) / 
                  (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))    ); 
  return(Vcmax); 
}
r_Vcmax2 <- function(kopt, Hd, Ha, Tk, Topt){
  Vcmax <- kopt * ( (Hd * pow(2.718282,((Ha*(Tk-Topt))/(Tk*0.008314*Topt))) ) / 
                      (Hd - (Ha*(1-(2.718282^((Hd*(Tk-Topt))/(Tk*0.008314*Topt))))))    ); 
  return(Vcmax); 
}

r_Vcmax(1, 200, 400, 304, 305)
r_Vcmax2(1, 200, 400, 304, 305)

# a = shape, b=rate
expose_stan_functions("gamma(1,1)")
dgamma(1, shape = 5, rate = 1)

