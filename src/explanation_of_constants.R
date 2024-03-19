
# miscallenous constants

# we take information from 

sideralyear = 365*24*3600 + 6*3600 + 9*60 + 10
tropicalyear = 365.2422*24*60*60
#n=2*pi/sideralyear
n=2*pi/tropicalyear
CAC_ref <- 0.00326742
# e_ref <- 0.017
e_ref <- 0.0
ece_ref <- 0.055

# i think we kind find arguments for using solar day rather than sideral day
# because the rotation of the Earth around the sun is indeed part of its angular momentum
# So it is indeed correct to say that, in a TROPICAL year, the earth is doing 364.24 revolution, 
# sugesting taht the 

# notes: p. 107 : L est le moment angulaire de rotation
# L = par définition C*omega, donc C/L =  1/omega (p. 107 apres eq. 52)

# alpha p. 112 est ml/mtheta * (atheta_al)^3

alpha_ref = 1/(81.30)/(332958)*(149600e3/384400)^3.
c_ref = 5.15*pi/180  # inclination of lunar orbit

# two definitions for omega: one for sideral day and one for solar day 
# thinking of it I guess we should use the omega_solar

omega_sideral=2*pi/(23*60*60+56*60+04)
omega_solar=2*pi/(24*60*60)

# apo is the variable called P0 in Berger's thesis  and
# is referred to as a Newcomb constant
# 
# 

apo <- function(CAC=CAC_ref, e=e_ref) {
  apo <- n*n/omega_solar*CAC*1.5*(1-e^2)^(-1.5)
  apo <- apo*60*60* tropicalyear *(180/pi)
  apo 
}

# another possibility, to be very close to the value he has, is 
# to use, for omega, the solar day and not the sideral day. 
# not entirely clear which one to use. For me it should be the sideral day. 
# so there may be an error, here

ala <- function(CAC=CAC_ref, alpha = alpha_ref, e = e_ref, ece=ece_ref , c=c_ref){
  # e_ref = 0.017
  # ece = 0.055
  # I am not sure these are the values used by Berger in his thesis
  # they are taken here from Berger et al., Terra Research 1989

  # this formula is  the formula for "l_0"  that comes 
  # from his formula for "l"  given en p. 112
  # it uses the fact that on an Keplerian orbit, 
  # the average of 1/r^3 = a*(1-e^2)^(-3/2)
  # see the attached maxima nootebok in scientific doc
  # it comes out that, p. 112, d\ell/de^2 = -3/2*a*(-5/2)
  # 
  moonfactor = 1-1.5*sin(c)*sin(c)
  ala <-  1.5*n*n/omega_solar * CAC * (alpha*(1-ece^2)^(-1.5)*moonfactor + (1-e^2)^(-1.5) )
  ala <- ala*60*60* tropicalyear *(180/pi)
  ala
}

# tests

print ('apo : Berger reference is 17.3919')
print ( apo())

# the closest value ot Berger is found for usig sideralyear, solar day (!) 
# and e_ref = 0. 
# the difference I have with him is one precession cycle over 20e6 se we are
# clearly within the error bar. 

print ('ala : Berger reference is 54.9066')
print ( ala() )




# je trouve ici 54.96646 au lieu de ses 54.9066. Environ 0.1% de differene. C'est très sensible au choix de l'eccentricité lunaire et de 'c' et je ne retrouve pas les valeurs qu'il a utilisées 


# the constants he provided are given p. 114 of his thesis, eq 67. SO I did not quite manage to
# reproduce them here, but yet very close, indeed. 
