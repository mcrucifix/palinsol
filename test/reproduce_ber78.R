
require(palinsol)
data(Bre73)
data(La88)

La88_shifted <- La88
shifted_io <- La88$io$Phases - 50*La88$io$Freq 
shifted_epi <- La88$epi$Phases - 50*La88$epi$Freq 

La88_shifted$io$Phases <- shifted_io
La88_shifted$epi$Phases <- shifted_epi

myBER78 <- SB67 (Bre73, 
                  P0 = 17.3919 *  pi/(3600*180.) ,
                  ell = 54.9087 * pi/(3600*180.) , 
                  epsilonbar = 23.394901 * pi/180., 
                  zeta = 3.392506 * pi / 180.,
                  aggregating = TRUE, ber78_strict = TRUE)

myBER90 <- SB67 (La88_shifted, 
                  P0 = 17.3919 *  pi/(3600*180.) ,
                  ell = 54.9066 * pi/(3600*180.) , 
                  epsilonbar = 23.399935 * pi/180., 
                  aggregating = TRUE, ber78_strict = FALSE)


myBER902 <- SB67 (La88, 
                  P0 = 17.3919 *  pi/(3600*180.) ,
                  ell = 54.9066 * pi/(3600*180.) , 
                  epsilonbar = 23.399935 * pi/180., 
                  aggregating = TRUE, ber78_strict = FALSE)

La88_shifted$epi$Phase[seq(10)]
La88$epi$Phase[seq(10)]
myBER90$OBLIQUITY$Phase[seq(10)]
myBER902$OBLIQUITY$Phase[seq(10)]

# ber90
ell = 54.9066 * pi/(3600*180.) 
epsilonbar = 23.399935 * pi/180.
# 
ell * cos(epsilonbar) * 3600*180/pi
# ell * cos(epsilonbar) * 3600*180/pi * 1.002621

times <- seq(-2e6, 0, 1e3)

O <- sapply(times, ber78)
eccpalinsol <-  ts ( O['ecc',], start = -2e6, deltat = 1e3)
cppalinsol <-  ts ( O['esinw',], start = -2e6, deltat = 1e3)
obpalinsol <-  ts ( O['eps',], start = -2e6, deltat = 1e3)


O90 <- sapply(times, ber78)
eccpalinsol90 <-  ts ( O['ecc',], start = -2e6, deltat = 1e3)
cppalinsol90 <-  ts ( O['esinw',], start = -2e6, deltat = 1e3)
obpalinsol90 <-  ts ( O['eps',], start = -2e6, deltat = 1e3)
# eccentricity


e <- Mod ( develop(Bre73$epi, -2e6, 0, 1e3, dfunction=cis) )
Pi  <- Arg ( develop(Bre73$epi, -2e6, 0, 1e3, dfunction=cis) )
Psi  <-  ( develop(myBER78$PSI, -2e6, 0, 1e3) ) 
Obl  <-  ( develop(myBER78$OBLIQUITY, -2e6, 0, 1e3, maxfreq=240) ) 

e90 <- Mod ( develop(La88$epi, -2e6, 0, 1e3, dfunction=cis) )
Pi90  <- Arg ( develop(La88$epi, -2e6, 0, 1e3, dfunction=cis) )
Psi90  <-  ( develop(myBER90$PSI, -2e6, 0, 1e3, maxfreq=1000 ) )
Obl90  <-  ( develop(myBER90$OBLIQUITY, -2e6, 0, 1e3))




omegatilde <- (Pi + Psi + pi) %% (2*pi)
climprecess <- ts( e*sin(omegatilde), start=-2e6, deltat=1e3)


omegatild90 <- (Pi90 + Psi90 + pi) %% (2*pi)
climprecess90 <- ts( e*sin(omegatild90), start=-2e6, deltat=1e3)

plot(eccpalinsol)
lines(e - eccpalinsol, col='red')

# e sin (omegatilde) 


plot(cppalinsol, xlim=c(-2e6, -1.8e5))
lines(climprecess, col='red')
lines(cppalinsol - climprecess, col='blue')


plot(obpalinsol, xlim=c(-2e6, -1.8e5))
lines(Obl, col='red')
plot(Obl - obpalinsol, col='red')


plot(obpalinsol90)
lines(Obl90 , col='red')


plot(cppalinsol90, xlim=c(-2e6, -1.8e5))
lines(climprecess90, xlim=c(-2e6, -1.8e5), col='red')
plot(cppalinsol90 - climprecess90) 

myBER90$OBLIQUITY$Phase[seq(10)] * 180/pi 
myBER902$OBLIQUITY$Phase[seq(10)] * 180/pi 


myBER90$PSI$Phase[seq(10)] * 180/pi  
myBER902$PSI$Phase[seq(10)] * 180/pi 

myBER90$OBLIQUITY$Freq[seq(10)] * 180/pi 
myBER902$OBLIQUITY$Phase[seq(10)] * 180/pi 

plot ( myBER90$OBLIQUITY$Freq[seq(180)] * 3600 * 180/pi -  BER90$Table1$Freq[seq(180)] )

plot ( myBER90$OBLIQUITY$Amp[seq(180)] * 3600 * 180/pi ) 


plot ( myBER90$OBLIQUITY$Phase[seq(80)] *  180/pi - BER90$Table1$Phase[seq(80)] )

