require(palinsol)
dyn.load(file.path('berger78','insol.so'))


PI=pi
PIR=pi/180.
PIRR=PIR/3600. 
STEP = 360 / 365.25
TEST = 0.0001
TAU = 1.
# TAU = 86.


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


palinsol::Insol(orbit=c(ecc=ecc, varpi=xl, eps=eps), long=tls,  lat=60*pi/180, S0=S0)


for (ecc in seq(7)*0.01)
  for (xl in seq(0,10)*2*pi/10)
    for (eps in seq(15,30,2)*pi/180)
       for (lat in seq(-90,90,10)){
         so = sin(eps)       
         S0 = 1368
berger_tls <- sapply(seq(120)*3, function(tls) {
        berger_DAYINS (ecc, xl*180/pi, so, tls , S0,  lat=lat)
})

palinsol_tls <- sapply(seq(120)*3, function(tls) {
         palinsol::Insol(orbit=c(ecc=ecc, varpi=xl, eps=eps), long=tls*pi/180,  lat=lat*pi/180, S0=S0)
})

tol = 1.e-12

error = max(abs(berger_tls['ww',] - palinsol_tls))

if (max(abs(berger_tls['ww',] - palinsol_tls)) < tol) {
   print(sprintf("PASSED ecc=%.2f lat=%.2f eps=%.2f xl=%.2f ", ecc, lat, eps, xl))
} else {
  print(sprintf("ecc %.2f lat %.2f eps %.2f xl %.2f : FAILED by %.4f  ", ecc, lat, eps, xl, error))

}}

