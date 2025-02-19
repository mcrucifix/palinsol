
# miscallenous constants

# we take information from 

# sideralyear = 365*24*3600 + 6*3600 + 9*60 + 10
#n=2*pi/sideralyear
tropicalyear = 365.2422*24*60*60
n=2*pi/tropicalyear
CAC_ref <- 0.00326742   # Earth dynamical ellipticity

# note: see e.g. Laskar 2004 and others
# "for a fast rotating planet like the Earth, the 
# "dynamical ellipticity Ed = (2C-A-B)/C  can be 
# considered to be proporotinal to omega^2. 
# 
# e_ref <- 0.017
e_ref <- 0.0      # earth orbit reference excentricity, used as a basis for 
                  # the Taylor development in Berger's codes
ece_ref <- 0.055  # moon orbit excentricity

# i think we kind find arguments for using solar day rather than sideral day
# because the rotation of the Earth around the sun is indeed part of its angular momentum
# So it is indeed correct to say that, in a TROPICAL year, the earth is doing 364.24 revolution, 
# sugesting taht the 

# notes: p. 107 : L est le moment angulaire de rotation
# L = par définition C*omega, donc C/L =  1/omega (p. 107 apres eq. 52)

# alpha p. 112 est ml/mtheta * (atheta_al)^3

earth_moon_distance_ref = 384400
c_ref = 5.15*pi/180  # inclination of lunar orbit

# two definitions for omega: one for sideral day and one for solar day 
# thinking of it I guess we should use the omega_solar

omega_sideral=2*pi/(23*60*60+56*60+04)
omega_solar=2*pi/(24*60*60)

# apo is the variable called P0 in Berger's thesis  and
# is referred to as a Newcomb constant
# 
# alpha is the mass ratio multiplied by the semi-great axes ratios of moon and earth. 
# 

# - this is the public handle 

#' Compute the Newcomb constant and its sensitivity to e^, according to the Sharav-Budnikova model featured by Berger 1973. 
#'
#' -- 
#' 
#' Implementation of expressions for "ell" and "P0" p. 112 of Berger 1973. These variables
#' are fundamental in the computation of the luni-solar precession (also called axial precession)
#' In first approximation, the precession frequency (the term k) appearing in astronomical developments
#' is equal to "ell"*cos(e_bar), where e_bar is the reference value for obliquity. 
#' d(ell)/d(e^2) = 3/2 P0, where "e" is the obliquity. All parameters are given default values
#' coresponding to the present-day state. 
#' 
#' @param CAC is (C-A)/C where A and C are Earth's moments of inertia aroud the plora and equatioral axis of inertia. Non-dimensional
#' @param omega is the rotational angular velocity of the Earth, measured in a geocentric framework, in radians per second
#' @param c is the inclination of the  lunar orbit on the ecliptic, in radians
#' @param ece : moon orbit excentricity
#' @param earth_moon_distance is expressed in km
#' @return a list with "ell" and "P0", ready for use in `compute_tables`, in arc seconds per year. 
#' @author Michel Crucifix, UC Louvain, Belgium
#' @references Berger, A. L. (1973)
#' Berger A. L. (1973), The astronomical theory of paleoclimates, Theses, UCLouvain, 1973. Currently not available in digital version. 
#' Berger A., M. Loutre and V. Dehant (1989), Astronomical frequencies for pre-Quaternary palaeoclimate studies, Terra Nova, (1) 474–479 doi:10.1111/j.1365-3121.1989.tb00413.x
#' @examples
#' # we consider the AstroGeo model of Farhat et al., Astron. Asttroph. 665
#' # length of day at -370 Ma : 21.81
#' # semi-major axis of moon: 361000km
#' N = newcomb_parameters (omega = 2*pi/(21.81*3600), earth_moon_distance = 361000)
#' print (N)
#' # assuming, still according to the AstroGeo22 model, an obliquity of 22.2 degrees, 
#' # this simple model suggests a precession rate of
#' N$ala * cos(22.2 * pi/180)
#' # MORE WORK IS NEEDED HERE BECAUSE  C-A/C also varies with OMEGA
#'
#' @export newcomb_parameters
newcomb_parameters <- function(CAC=CAC_ref, omega = omega_solar, c=c_ref, ece = ece_ref, earth_moon_distance = earth_moon_distance_ref){
  apo_parameter <- ala(CAC, 0., omega)
  ala_parameter <- ala(CAC, earth_moon_distance,  0., ece, c, omega)
  return (list(ala = ala_parameter, apo = apo_parameter))
}

# - the latter two functions remain internal. 

apo <- function(CAC=CAC_ref, e=e_ref, omega = omega_solar) {
  apo <- n*n/omega*CAC*1.5*(1-e^2)^(-1.5)
  apo <- apo*60*60* tropicalyear *(180/pi)
  apo 
}

# another possibility, to be very close to the value he has, is 
# to use, for omega, the solar day and not the sideral day. 
# not entirely clear which one to use. For me it should be the sideral day. 
# so there may be an error, here

ala <- function(CAC=CAC_ref, earth_moon_distance = earth_moon_distance_ref, e = e_ref, ece=ece_ref , c=c_ref, omega = omega_solar){
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

  alpha = 1/(81.30)/(332958)*(149600e3/earth_moon_distance)^3.
  moonfactor = 1-1.5*sin(c)*sin(c)
  ala <-  1.5*n*n/omega * CAC * (alpha*(1-ece^2)^(-1.5)*moonfactor + (1-e^2)^(-1.5) )
  ala <- ala*60*60* tropicalyear *(180/pi)
  ala
}

## # tests
## 
## print ('apo : Berger reference is 17.3919')
## print ( apo())
## 
## # the closest value ot Berger is found for usig sideralyear, solar day (!) 
## # and e_ref = 0. 
## # the difference I have with him is one precession cycle over 20e6 se we are
## # clearly within the error bar. 
## 
## print ('ala : Berger reference is 54.9066')
## print ( ala() )
## 



# je trouve ici 54.96646 au lieu de ses 54.9066. Environ 0.1% de differene. C'est très sensible au choix de l'eccentricité lunaire et de 'c' et je ne retrouve pas les valeurs qu'il a utilisées 


# the constants he provided are given p. 114 of his thesis, eq 67. SO I did not quite manage to
# reproduce them here, but yet very close, indeed. 
