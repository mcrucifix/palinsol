
# miscallenous constants

sideralyear = 365*24*3600 + 6*3600 + 9*60 + 10
tropicalyear = 365.2422*24*60*60
n=2*pi/sideralyear
# n=2*pi/tropicalyear
CAC_ref <- 0.003267


alpha_ref = 1/(81.30)/(332958)*(149600e3/384400)^3.
c_ref = 5.15*pi/180  # inclination of lunar orbit

# two definitions for omega: one for sideral day and one for solar day 

omega_sideral=2*pi/(23*60*60+56*60+04)
omega_solar=2*pi/(24*60*60)

# apo is the variable called P0 in Berger's thesis  and
# is referred to as a Newcomb constant
# 
# 

apo <- function(CAC=CAC_ref) {
  apo <- n*n/omega_solar*CAC*1.5*(1-0.017^2)^(-1.5)
  apo 
}

# i find his value for e_ref = 0.044. 
# another possibility, to be very close to the value he has, is 
# to use, for omega, the solar day and not the sideral day. 
# not entirely clear which one to use. For me it should be the sideral day. 
# so there may be an error, here

ala <- function(alpha = alpha_ref, eref = 0.017, ece=0.055, c=c_ref) {

  # eref = 0.017
  # ece = 0.055
 

  # this formula is  the formula for "l_0"  that comes 
  # from his formula for "l"  given en p. 112
  # it uses the fact that on an ellipse  
  # 
  ala <-  1.5*n*n/omega_sideral * CAC * (alpha*(1-eref^2)^(-1.5) + (1-ece^2)^(-1.5))
  
  moonfactor = 1-1.5*sin(c)*sin(c)
  ala <- ala*60*60* sideralyear *(180/pi)  * moonfactor
  ala
}

# tests

print ('apo : Berger reference is 17.3919')
print ( apo()*60*60* sideralyear *(180/pi) ) 
print (apo()*60*60* tropicalyear *(180/pi) )

print ('ala : Berger reference is 54.9066')
print ( ala() )




# je trouve ici 54.67759 au lieu de ses 54.9066. Environ 0.5% de difference.


# the constants he provided are given p. 114 of his thesis, eq 67. SO I did not quite manage to
# reproduce them here, but yet very close, indeed. 
