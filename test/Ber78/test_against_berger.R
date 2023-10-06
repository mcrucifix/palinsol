dyn.load(file.path('berger78','insol.so'))


PI=pi
PIR=pi/180.
PIRR=PIR/3600. 
STEP = 360 / 365.25
TEST = 0.0001
TAU = 86.4


berger_DAYINS <- function(ecc, xl, so, tls, S0, lat) {
	SF=S0*TAU/pi
	WW=as.numeric(0.)
	DAYL=as.numeric(0.)
	OUT = .Fortran("DAYINS", ecc=ecc, xl=xl, so=so, tls=tls, lat=lat, pir=PIR, pi=PI, test=TEST, SF=SF, WW=WW, DAYL=DAYL)
	return(c(ww=OUT$WW, W = OUT$WW/86.400, dayl = OUT$DAYL))
}

# tls : true solar longitude
# xl  : true solar longitude of perihelion + 180 [making it heliocentric]

ecc = 0.07
xl = pi
eps = 23*pi/180
so = sin(eps)       
tls = 0. 
S0 = 1368

print ( berger_DAYINS (ecc, xl, so, tls + pi, S0,  lat=60*pi/180) )


require(palinsol)

palinsol::Insol(orbit=c(ecc=ecc, varpi=xl, eps=eps), long=tls,  lat=60*pi/180, S0=S0)


