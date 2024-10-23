# require(palinsol)
load_all()
data(Bre73)
data(La88)

La88_shifted <- La88
shifted_io <- La88$io$Phases - 50*La88$io$Freq 
shifted_epi <- La88$epi$Phases - 50*La88$epi$Freq 

La88_shifted$io$Phases <- shifted_io
La88_shifted$epi$Phases <- shifted_epi

myBER78 <- SB67 (Bre73, 
                  P0 = 17.3919 *  pi/(3600*180.) ,
                  ell = 54.9066 * pi/(3600*180.) , 
                #  epsilonbar = 23.394901 * pi/180., 
                  epsilonbar = 23.394901 * pi/180., 
                  zeta = 3.39250 * pi / 180.,
                  aggregating = FALSE, ber78_strict = TRUE, ber90_reproduce = FALSE)

myBER78a <- SB67 (Bre73, 
                  P0 = 17.3919 *  pi/(3600*180.) ,
                  ell = 54.9066 * pi/(3600*180.) , 
                #  epsilonbar = 23.394901 * pi/180., 
                  epsilonbar = 23.394852 * pi/180., 
                  zeta = 3.39250 * pi / 180.,
                  aggregating = TRUE, ber78_strict = TRUE, ber90_reproduce = FALSE)



myBER90 <- SB67 (La88_shifted, 
                  P0 = 17.3919 *  pi/(3600*180.) ,
                  ell = 54.9066 * pi/(3600*180.) , 
                  epsilonbar = 23.399935 * pi/180., 
                  aggregating = TRUE, ber78_strict = FALSE, ber90_reproduce = TRUE)


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


O90 <- sapply(times, ber90)
eccpalinsol90 <-  ts ( O90['ecc',], start = -2e6, deltat = 1e3)
cppalinsol90 <-  ts ( O90['esinw',], start = -2e6, deltat = 1e3)
obpalinsol90 <-  ts ( O90['eps',], start = -2e6, deltat = 1e3)
# eccentricity


e <- Mod ( develop(Bre73$epi, -2e6, 0, 1e3, dfunction=cis) )
Pi  <- Arg ( develop(Bre73$epi, -2e6, 0, 1e3, dfunction=cis) )
Psi  <-  ( develop(myBER78a$PSI, -2e6, 0, 1e3, maxfreq = 80 ) )
Obl  <-  ( develop(myBER78$OBLIQUITY, -2e6, 0, 1e3, maxfreq=240) ) 
Obl_t  <-  ( develop(myBER78$OBLIQUITY, -2e6, 0, 1e3, maxfreq=47) )  # 47 over 240 to reproduce BER78

e90 <- Mod ( develop(La88$epi, -2e6, 0, 1e3, dfunction=cis) )
Pi90  <- Arg ( develop(La88$epi, -2e6, 0, 1e3, dfunction=cis) )
Psi90  <-  ( develop(myBER90$PSI, -2e6, 0, 1e3, maxfreq=1000 ) )
Obl90  <-  ( develop(myBER90$OBLIQUITY, -2e6, 0, 1e3))


plot(Obl, xlim=c(-2e6, -1.9e6))
lines(Obl_t, col='red')
lines(obpalinsol, col='green')

plot(Obl-Obl_t)
plot(Obl-obpalinsol)


myBER78a$PSI$Amp[1] * 3600 * 180/pi
BER78$Table5$Amp[1]


plot (myBER78a$PSI$Amp[seq(78)] * 3600 * 180/pi)
points(abs(BER78$Table5$Amp), col='red')
plot (
       myBER78a$PSI$Amp[seq(78)] * 3600 * 180/pi -         
       abs(BER78$Table5$Amp))



# BER90 is reproduce down to machine limite
# BER78 obliquity is reproduced within 10.e-5
# BER78 clim precession: still some (very) small differences that do not seem to trend (check over 2 Ma) so this must be related to the number of terms retained in intermediary steps  // after a couple of hours I couldn't find out. 
