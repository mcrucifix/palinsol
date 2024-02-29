sideralyear = 365*24*3600 + 6*3600 + 9*60 + 10
tropcialyear = 365.2422*24*60*60

n=2*pi/sideralyear
omega=2*pi/(23*60*60+56*60+04)
omega=2*pi/(24*60*60)
CAC = 0.003267

apo = n*n/omega*CAC*1.5*(1-0.017^2)^(-1.5)

# i find his value for e_ref = 0.044. 
# another possibility, to be very close to the value he has, is 
# to use, for omega, the solar day and not the sideral day. 
# not entirely clear which one to use. For me it should be the sideral day. 
# so there may be an error, here

apo*60*60* sideralyear *(180/pi) 
apo*60*60* tropicalyear *(180/pi) 

alpha = 1/(81.30)/(332958)*(149600e3/384400)^3.
eref = 0.017
ece = 0.055

ala = 1.5*n*n/omega * CAC * (alpha*(1-eref^2)^(-1.5) + (1-ece^2)^(-1.5))

c = 5.15*pi/180  # inclination of lunar orbit

moonfactor = 1-1.5*sin(c)*sin(c)

ala*60*60* sideralyear *(180/pi)  * moonfactor

# je trouve ici 54.67759 au lieu de ses 54.9066. Environ 0.5% de difference.


# the constants he provided are given p. 114 of his thesis, eq 67. SO I did not quite manage to
# reproduce them here, but yet very close, indeed. 
