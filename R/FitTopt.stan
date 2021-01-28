functions {
 real calc_Vcmax(real kopt, real Ha, real Tk, real Topt) {
 real out; 

 out = kopt * ((200*pow(2.718282,((Ha*(Tk-Topt))/(Tk*0.008314*Topt))))/ 
                      (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))); 
 return out;
 }
}

data {
  int<lower=0> N;
  vector[N] Vcmax;
  vector[N] Tk;
}

parameters {
  real<lower=25,upper=500> kopt;
  // real<lower=150,upper=250> Hd;
  real<lower=25,upper=150> Ha; 
  real<lower=290,upper=335> Topt;
  real<lower=0> sigma;
}


model {
  kopt ~ normal(250,100);
  // Hd ~ normal(200,10);
  Ha ~ normal(75,25); 
  Topt ~ normal(310,10); 
  sigma ~ cauchy(0,7.5);
  for (n in 1:N){
    Vcmax[n] ~ normal(calc_Vcmax(kopt, Ha, Tk[n], Topt),sigma);
  }
}