functions {
 real calc_photosynthesis(real kopt, real Ha, real Tk, real Topt) {
 real out; 

 out = kopt * ((200*pow(2.718282,((Ha*(Tk-Topt))/(Tk*0.008314*Topt))))/ 
                      (200 - (Ha*(1-(2.718282^((200*(Tk-Topt))/(Tk*0.008314*Topt))))))); 
 return out;
 }
}

data {
  int<lower=0> N;
  vector[N] A;
  vector[N] Tk;
}

parameters {
  real<lower=1,upper=500> kopt;
  real<lower=1,upper=120> Ha; 
  real<lower=283,upper=335> Topt;
  real<lower=0> sigma;
}


model {
  kopt ~ normal(25,100);
  Ha ~ normal(50,20); 
  Topt ~ normal(310,10); 
  sigma ~ cauchy(0,7.5);
  for (n in 1:N){
    A[n] ~ normal(calc_photosynthesis(kopt, Ha, Tk[n], Topt),sigma);
  }
}