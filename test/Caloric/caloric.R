load_all()


# supply difference of insolation between a given time (mean longitude + 80.5) and
# 6 months later
diffInsol <- function(orbit, day, lat,  ... ){
  l1 <- day2l (orbit, day )
  l2 <- day2l (orbit, day + 180.)
  return ( 
   Insol( orbit = orbit, long = l2, lat=lat,  ...) - 
   Insol( orbit = orbit, long = l1, lat=lat,  ...) )
}

O <- c(ecc=0.05, varpi=330 * pi / 180, eps = 22.5 * pi / 180)

plot ( Milankovitch(O) , levels = seq(360,570,10 ))

abline(h=13)

diffInsol(O, 80, lat = 60 * pi / 180. )

firstguess <- function(orbit, lat){
  modvarpi <- (orbit['varpi']  + pi/2)  %% ( 2 * pi )  
  ifelse ( modvarpi  > pi, 
  (170 - 90 * tanh(1.7*lat)   + 360 ) %% 360  , 
  (-10 + 90 * tanh(1.7*lat)   + 360 ) %% 360  )
}

lat <- 62.2 * pi/180. 

fg <- firstguess(O, lat)

OUT <-  uniroot(function(day) diffInsol (O, day, lat), c(fg-30, fg+30) )


if (OUT$root){
  d1 <- OUT$root
  d2 <- OUT$root + 180. 
  sol <- Insol_d1d2(O, d1=d1, d2=d2, lat=lat, avg=FALSE)
}

