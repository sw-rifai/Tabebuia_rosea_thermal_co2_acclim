library(plantecophys)
# assumption: fixed amount of water in the air
curve(VPDleafToAir(VPD=x,Tleaf = 31, Tair = 30), 0.1,3, 
      xlim=c(0,3), ylim=c(0,3)); 
curve(VPDairToLeaf(VPD=x, Tleaf = 31, Tair = 30), 0.1,3,col='red',add=T)
abline(0,1,col='purple')

?plantecophys::FindTleaf
find_tleaf <- Vectorize(FindTleaf)
find_tleaf(gs=0.1, Tair=25, Wleaf=0.001)
curve(find_tleaf(gs=0.1, Tair=25, Wleaf=x),0.001,0.2)
curve(plantecophys::RHtoVPD(RH = 0.75, TdegC = x), 0,45)
curve(plantecophys::VPDairToLeaf(VPD = 1, Tleaf = x, Tair = 25), 20,30); abline(v=25,h=1,lty=3)

curve(plantecophys::VPDleafToAir(VPD = 1, 
                                 Tleaf = 25, 
                                 Tair = x), 20,30); abline(v=25,h=1,lty=3)
curve(plantecophys::VPDleafToAir(VPD = 1, Tleaf = x, Tair = 25), 20,30); abline(v=25,h=1,lty=3)

curve(VPDleafToAir(VPD = RHtoVPD(RH = 75, TdegC = 30), 
                   Tleaf = x, 
                   Tair = 30), 20,30); abline(v=25,h=1,lty=3)

VPDleafToAir(VPD = RHtoVPD(RH = 75, TdegC = 30), 
             Tleaf = x, 
             Tair = 30)

(esat(35)/1000)-RHtoVPD(100,35)
RHtoVPD(90,35)-RHtoVPD(100,35)

p_t <- function(p_opt, t_opt, t_leaf, omega){
  return(p_opt*exp(-( (t_opt-t_leaf)/omega )**2))
}
med <- function(g0, g1, D, A, Ca){
  g0 + 1.6*(1 + g1/sqrt(D))*A/Ca
}

curve(p_t(p_opt = 20, t_opt = 30, t_leaf=x, omega=9), 0,45)
abline(v=30)

med <- function(g0, g1, D, A, Ca){
  g0 + 1.6*(1 + g1/sqrt(D))*A/Ca
}
curve(med(0,10,x,20,400),0.1,5)
curve(med(g0 = 0,g1 = 10,D=1,A = 20,Ca = x),200,800)



# Co-limitation using VPD_air
curve(med(g0=0, 
          g1=10, 
          D=plantecophys::RHtoVPD(RH = 0.75, TdegC = x), 
          # A=20,
          A=p_t(20, 30, t_leaf = x, omega=9),
          Ca=400), 
      10,45, ylim=c(0,1))

curve(med(g0=0, 
          g1=10, 
          D=plantecophys::RHtoVPD(RH = 0.75, TdegC = x), 
          A=p_t(20, 30, t_leaf = x, omega=9),
          Ca=400), 
      10,45, ylim=c(0,1))
curve(med(0,10,RHtoVPD(0.75,x),A=20,400),10,45,col='darkgreen',add=T)


# Co-limitation using VPD_air
curve(med(g0=0, 
          g1=10, 
          D=plantecophys::RHtoVPD(RH = 75, TdegC = x),   # leaf and air temp perfectly coupled
          # D=VPDleafToAir(VPD = 1, Tleaf = x, Tair = 30),
          # D=VPDleafToAir(VPD = 1, Tleaf = 30, Tair = x),
          # D = VPDleafToAir(VPD = RHtoVPD(RH = 75, TdegC = x), 
          #                  Tleaf = 25, 
          #                  Tair = x),
          A=p_t(p_opt = 20, 
                t_opt = 30, 
                t_leaf = x, 
                omega=9),
          Ca=400), 
      20,35, 
      # ylim=c(0,10)
      ylab='g_s', add = F
      )
abline(v=30)
curve(med(0,10,RHtoVPD(0.75,x),A=20,400),10,45,col='darkgreen',add=T)



curve(med(1,10,x,20,400),0.1,5)


fn <- function(Ca, gs, g0, g1, D){
  return(Ca*(gs - g0)/(1.6*(1 + g1/sqrt(D))))
}
curve(fn(400, med(0.1,10,x,20,400), 0.1, 10, x),0.1,5)

med(0.1,10,1.5,20,400)



p_t <- function(p_opt, t_opt, t_leaf, omega){
  return(p_opt*exp(-( (t_opt-t_leaf)/omega )**2))
}
med <- function(g0, g1, D, A, Ca){
  g0 + 1.6*(1 + g1/sqrt(D))*A/Ca
}
find_tleaf <- Vectorize(plantecophys::FindTleaf)
find_tleaf(gs=0.1, Tair=25, Wleaf=0.04)

curve(med(g0=0, 
          g1=10, 
          # D=plantecophys::RHtoVPD(RH = 0.75, TdegC = x), # D_a
          D=VPDleafToAir(VPD = 2, Tleaf = x, Tair = 30), # D_l
          A=p_t(p_opt = 20, 
                t_opt = 30, 
                t_leaf = x, 
                omega=9),
          Ca=400), 
      25,35, 
      # ylim=c(0,1), 
      ylab='g_s', xlab='leaf_temp')
curve(med(g0 = 0,g1 = 10,
          D=VPDleafToAir(VPD = 2, Tleaf = x, Tair = 30),  # D_l
          # D = RHtoVPD(0.75,x),A=20,Ca = 400), # D_a
          A = 20,
          Ca = 400),
      10,45,col='darkgreen',add=T)
abline(v=30, lty=3)
      


# Photo Topt response
p_t <- function(p_opt, t_opt, t_leaf, omega){
  return(p_opt*exp(-( (t_opt-t_leaf)/omega )**2))
}

# Medlyn gs model 
med <- function(g0, g1, D, A, Ca){
  g0 + 1.6*(1 + g1/sqrt(D))*A/Ca
}
find_tleaf <- Vectorize(plantecophys::FindTleaf)

# Atmospheric VPD with constant RH and with Topt control on A
curve(med(g0=0, 
          g1=10, 
          D=plantecophys::RHtoVPD(RH = 0.75, TdegC = x), # D_a
          # D=VPDleafToAir(VPD = 2, Tleaf = x, Tair = 30), # D_l
          A=p_t(p_opt = 20, 
                t_opt = 30, 
                t_leaf = x, 
                omega=9),
          Ca=400), 
      25,35, 
      ylim=c(0,1),
      ylab='g_s', xlab='leaf_temp'); abline(v=30, lty=3,col='darkgreen')

#--- leaf-to-air VPD (?) with Topt control on A ---
# function: Photo Topt response
p_t <- function(p_opt, t_opt, t_leaf, omega){
  return(p_opt*exp(-( (t_opt-t_leaf)/omega )**2))
}
# function: Medlyn gs model 
med <- function(g0, g1, D, A, Ca){
  g0 + 1.6*(1 + g1/sqrt(D))*A/Ca
}
curve(med(g0 = 0,
          g1 = 10,
          D=VPDairToLeaf(VPD = 2, Tleaf = x, Tair = 30),  # D_l
          # D = RHtoVPD(0.75,x), # D_a
          A=p_t(p_opt = 20, 
                t_opt = 30, 
                t_leaf = x, 
                omega=9),
          Ca = 400),
      25,35,col='darkgreen',add=F, 
      ylab='g_s', xlab='leaf_temp', ylim=c(0.25, 1))
curve(med(g0=0,g1=10,D=VPDairToLeaf(VPD = 2, Tleaf = x, Tair = 30),A=20, Ca=400), 
      25,35,col='blue',add=T)
abline(v=30, lty=3); 
legend('bottomleft',
        legend=c('Medlyn with A=20','Medlyn with A ~ f(leaf_temp)'),
       col=c('blue','darkgreen'), lty=c(1,1))
