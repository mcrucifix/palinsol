# #########################################################
# 
# WARNING NOTICE
#
# #########################################################
# 
# This code is a R-implementation of code developed by Andre Berger for his
# thesis (Berger, 1973, UCLouvain), and further developed until publication of
# the Berger and Loutre 1991 article (BL91). This code by Berger (and Loutre)
# includes two Fortran subprograms, called  p7505ff.f and  p75133f.f. To
# the best of my knowledge, it has never been officially released and so cannot
# be explicitly transcribed here. However, this code was used to check the
# validity of the R code below. The present code can be viewed as the legacy of many
# years is of  Berger and Loutre's work and now made available to the public by means
# of this package. I commented places were I could not exactly will produce
# their output and take full responsibility for any mistake or omission. 
# 
# The code takes a trigonometric development of the planetary solution
# (typically the h,q,p,k development of Bretagnon or Laskar) and, using the
# personal equation to the sick and ordered to the masses and eccentricity,
# generates trigonometric developments for climatic precession and obliquity.
# This is original work by Berger,  though he recognises that a similar work
# had been previously undertaken by Sharaf and Budnikova (1967). Unfortunately,
# that original reference published in the former USSR is not available to me. 
# 
# The thesis of André Berger has been scanned but it is not publicly released.
# It is however accessible from this author. 
#
# Below: 
# BL91 : Berger A. and M. F. Loutre (1991), Insolation values for the climate
# of the last 10 million years, Quaternary Science Reviews, (10) 297 - 317
# doi:10.1016/0277-3791(91)90033-Q
# 
# SB67 : Sharaf S. G. and N. A. Boudnikova (1967), Secular variations of
# elements of the Earth's orbit which influences the climates of the geological
# past, Tr. Inst. Theor. Astron. Leningrad, (11) 231-261 
# (this paper remained unaccessible to us)
#
# 
# Michel Crucifix, 2025-02-19
#
# #########################################################

# The following code does all the dirty work but is not to be called by the user. 
# 
# P0 = 17.3919
# ell = 54.9066 
# 
# ellradiants <- ell * (2*pi/3600/360)
# P0radiants  <- P0 * (2*pi/3600/360)

SB67_Internal <- function(ell, P0, epsilonbar, zeta = 1.600753*pi/180,   
                          PlanetarySolution, aggregating = FALSE, 
                          ber78_strict=FALSE, ber90_reproduce = FALSE) {

# in the Berger and Loutre procedure: 
# the inputs are "h" and "alpha", so in this notation
# epsilonbar and zeta.
# psi, psibar and epsilonbarprime are computed internally

# We consider that this routine is the one that will provide psi, psibar and epsilonbarprime on this basis

# However, epsilon_zero and psi_zero are constrained by observations
# this constraints epsilonbarstar and alpha via the equations 25 and 26 in BL91
# this seems easy to do.....
# but I am less sure how to obtain "h" (epsilonbar). That is equation 24
# maybe it is self converging, i.e., first guess 'h' by epsilonbar = h - sum
# (... epsilon bar), then replace epsilonbar etc. It should work


EPI <- as.data.frame(PlanetarySolution$epi)
# fragile !!!! 

I2OM <- with(PlanetarySolution$epi, list(
            N = Amp, 
            sigma = Freq, 
            delta = Phases))

colnames(EPI) <- c('g','M','beta')
if (!is.null(PlanetarySolution$i2o))  {

I2OM <- as.data.frame(PlanetarySolution$i2o)

# EPI <- read.table('ber90/epi.dat')

colnames(I2OM) <- c('sigmaprime', 'Nprime', 'deltaprime')

n <- length(I2OM$Nprime)

# I will privilege here transparent writing, even if less efficient. Fortran compiler would be better

# line 1

TERM1 <- list()
TERM1$sigma <- I2OM$sigmaprime
TERM1$N <- I2OM$Nprime
TERM1$delta <- I2OM$deltaprime

for (i in seq(n)) TERM1$N[i] <- TERM1$N[i] * ( 2 - I2OM$Nprime[i]^2 - 2 * sum(I2OM$Nprime[-i]^2) ) 

# this corresponds to lines p7505810 in his code 
# line 2 (TERM2)

TERM2 <- list()

TERM2$N <- - outer(I2OM$Nprime, I2OM$Nprime, function(Ni,Nj) Nj^2*Ni);
TERM2$sigma <- outer(I2OM$sigmaprime, I2OM$sigmaprime, function(Si,Sj) (2*Sj-Si))
TERM2$delta <- outer(I2OM$deltaprime, I2OM$deltaprime, function(Si,Sj) (2*Sj-Si))


for (h in seq(ncol(TERM2$N))) TERM2$N[h,h] = NA
for (h in seq(ncol(TERM2$sigma))) TERM2$sigma[h,h] = NA
for (h in seq(ncol(TERM2$delta))) TERM2$delta[h,h] = NA

TERM2$N <- as.numeric(TERM2$N)
TERM2$sigma <- as.numeric(TERM2$sigma)
TERM2$delta <- as.numeric(TERM2$delta)

keep <- which(!is.na(TERM2$sigma))

TERM2$N <- TERM2$N[keep]
TERM2$sigma <- TERM2$sigma[keep]
TERM2$delta <- TERM2$delta[keep]


# line 2 (TERM3)

TERM3 <- list()
TERM3$N <- -2 * outer ( outer(I2OM$Nprime, I2OM$Nprime[-n]), I2OM$Nprime)

TERM3$sigma<- outer( outer(I2OM$sigma, I2OM$sigma[-n], function(si, sj) sj-si), I2OM$sigma, function(sj,sk) sj+sk)

TERM3$delta<- outer( outer(I2OM$delta, I2OM$delta[-n], function(si, sj) sj-si), I2OM$delta, function(sj,sk) sj+sk)

for (i in seq(n)) {
  for (j in seq(n-1)) {
    for (k in seq(n)) {
      if ((k >= j) || (k == i) || (j == i) ){
        TERM3$sigma[i, j, k ] = NA;
        TERM3$delta[i, j, k ] = NA;
        TERM3$N[i, j, k ] = NA;
      }
    }}}

TERM3$sigma <- as.numeric(TERM3$sigma)
TERM3$delta <- as.numeric(TERM3$delta)
TERM3$N <- as.numeric(TERM3$N)


keep <- which(!is.na(TERM2$sigma))

TERM3$N <- TERM3$N[keep]
TERM3$sigma <- TERM3$sigma[keep]
TERM3$delta <- TERM3$delta[keep]


IOM <- list(
            N = c(TERM1$N, TERM2$N, TERM3$N), 
            sigma = c(TERM1$sigma, TERM2$sigma, TERM3$sigma), 
            delta = c(TERM1$delta, TERM2$delta, TERM3$delta))


# THE IOM LIST IS VALIDATED. AT LEAST IN AMPLITUDES, AND FOR THE FIRST SERIES OF 80 TERMS
tol <- 1.e-7

nkeep <- which (IOM$N > tol)

IOM <- list( N = IOM$N[nkeep], 
             sigma = IOM$sigma[nkeep], 
             delta = IOM$delta[nkeep]) 




} else  if (!is.null(PlanetarySolution$io))  {

# the planetary solution was supplied in terms of sini , not sin (i/2)
I2OM <- PlanetarySolution$io
IOM <- list(
            N = I2OM$Amp, 
            sigma = I2OM$Freq, 
            delta = I2OM$Phases)
} else {stop ("missing either io or i2o in planetary solution") }

#### development of obliquity

#psi     = 0.0002445400207  + 0.0002443015 - 0.0002445400 

# epsilonbar is the "h" in Berger and Loutre 1991
# epsilonbarprime is the "e_star" in Berger and Loutre
# psi is the "k" in Berger and Loutre
# psibar is the "k_tilde" in Berger and Loutre
# "ell" is the "P_bar" in Berger and Loutre 1991
# zeta is the "alpha" in Berger and Loutre 1991

# note to self: 
# still need to work out whether the psi estimated this way
# is a firt guess and should then be recopmuted with esplinobarstar (iterating, then)
# or whether this is intended
# if need be recomputed, perhaps better to supply the right value
# straightaway if, for example, supplied by Laskar ? 

psi <- ell * cos(epsilonbar)  # the "k" in Berger and Loutre 1991


# so in the following we will use "psi" in all the developments
# this is as BER90. Berger and Loutre, p. 303 have a note about
# this, and explain that it differs a bit from the SB procedure. 
# They also refer to Berger 1988 for more details, I haven't checked yet. 


# print (sprintf('Iteration %d; psi = %.9f', i, psi * 3600 * 180 / pi))

### ATTENTION: 
### 1. there is an extra ai^2 in eq. 21 or BER90 that is not in his code
### 2. his code does something weird and does not include
###    the 3 p0/p term. Need to understand and check this. 


# so for anyuse of epsilonbar below, see remark above

ng <- length(EPI$g)
tane = tan(epsilonbar)
cote = 1./tane

###########

itermax = 1;


if (ber78_strict) itermax = 2 

temp =  outer(EPI$M, EPI$M)
ctemp = outer(EPI$beta, EPI$beta,  function(bi,bk) { cos(bi-bk) })
dij = outer(seq(ng), seq(ng), function(i,k) {i < k})


for (iter in seq(itermax)){

# thesis Berger p. 110 = Sharaf Budnikova

deltaprime <- IOM$delta + zeta

fi <-  psi + IOM$sigma
ci <-  psi / fi 


if (ber90_reproduce){

Constraint <- (1 - sum(IOM$N^2*(1.5 + 0.75 * ci^2 - 2.5 * ci - 0.5 * (ci*(ci-1))*tane^2)) )

psibar <- ell * cos(epsilonbar) * Constraint

print (sprintf("psibar without P0 correction = %.5f", psibar * 180/pi * 3600))

# this is what you need to do to actually reproduce the Berger and Loutre paper, but this contradicts
# what they say!!! 

} else {
Constraint <- (1 - sum(IOM$N^2*(1.5 + 0.75 * ci^2 - 2.5 * ci - 0.5 * (ci*(ci-1))*tane^2)) - 3 * P0 / ell *  sum(temp*ctemp*dij, na.rm=TRUE))


psibar <- ell * cos(epsilonbar) * Constraint
if(ber78_strict)  psibar <- 50.43927 * pi/180 / 3600
print (sprintf("psibar with P0 correction = %.5f", psibar * 180/pi * 3600))

}

if(ber78_strict)  psi <- psibar
# dirty hack

}


######

# the first part of the Constraint (until -3po/ell) is validated
#  in his code, he actually loops to satisfy the initial condition constrain
# you must look at line 214, where the loop begins
# h (epsilonbar here) is first given a first guess. He computes all terms
# with that one, including the precession (!)
# then he only computes the first part of the constrain and tests a new
# psibar, whele the second part of the constrain (-3 P0/ell etc.) is 
# replaced by the difference you should have to satisfy the observation 
# of the current rate. 
# By newton raphson he infers a new 'h', which generates another psi
# and the iteration goes on. 
# His approach is beautiful because he has analytical integrals. 
# We could also do it brute force if we had constraints to satisfy. 





nf <- length(ci)

C1f <- -ci
C2fik <- -0.5 * psi * outer (seq(nf), seq(nf), function(i,k) {
                       1./(fi[i]+fi[k])*((ci[i]^2 + ci[k]^2 - 2) * tane +
                              (ci[i]+ci[k]) * cote)  })


C2fii = -.25 *(ci * (ci^2 - 1) * tane + ci^2 * cote )

C2fikprime <- -0.5 * psi * outer (seq(nf), seq(nf), function(i,k) {
                       1./(fi[i]-fi[k])*((ci[i]^2 - ci[k]^2 + 2*ci[k] - 2*ci[i]) * tane +
                              (ci[i]-ci[k]) * cote)  })



#### DEVELOPPEMENT FOR PSI         

# P0 et ell ont ete repris de sa these, sauf que pour ell il y a peut etre un facteur (1-e_0)3/2




# =====================================
# PRECESSION
# =====================================

#


D1f <- ci * (cote + (ci - 1) * tane)

D2f <- 0.25 * ci * (ci^2 + ci - 1) + 0.125 * ci^2*(ci-1)^2*(tane^2) + .5 * ci^2 * cote^2 

# ligne p7507950

# attention au signe moins de C2fik: corriger dans la these de berger c'est une erreur

D2fik <- psi * outer(fi,fi, function(fi, fk){1/(fi + fk)}) * (
          0.5 * outer (ci, ci,  function(ci,ck) {ci^2+ck^2+ci+ck-ci*ck-2} )
         - C2fik * tane
         - 0.5 *   outer (ci, ci,  function(ci,ck) {ci*(ci-1)+ck*(ck-1)}) * tane^2
         + outer (ci,ci, '+')*cote^2 )


D2fikprime <- psi * outer(fi,fi, function(fi, fk){1/(fi - fk)}) * (
          0.5 * outer (ci, ci,  function(ci,ck) { 5 * (ci+ck) - (ci^2+ck^2) - ci*ck - 6} )
         - C2fikprime * tane
         + 0.5 *  outer (ci*(ci-1), ci*(ci-1),'+') * tane^2) 







# must define here P0, ell but that's another part of the code. 

Dfseconde = 3 * psi * P0 / ell * outer(EPI$g, EPI$g,  function(gi,gk) { 1./(gi-gk)})

# validated

D1 <- D1f - cote   # p. 113 of thesis, first equation

D2 <- D2f - (ci - 0.5)*cote^2 - 0.5 * (ci^2 - 0.5)
D2ik <- D2fik - (outer(ci, ci, '+')-1) * cote^2 - 0.5 *(outer(ci^2, ci^2,'+')-1)
D2ikprime <- D2fikprime - 0.5 * ( outer(ci^2, ci^2,'-')) + outer(ci, ci, '-') 


C1 = C1f + 1
C2ii = C2fii + 0.5*D1f - 0.25 * cote
C2ik = C2fik + 0.5 * outer(D1f, D1f, '+') - 0.50 * cote
C2ikprime = C2fikprime - 0.5 * outer (D1f, D1f, '+') + 0.50 * cote


ETERM1 <- list(  N = C1 * IOM$N, sigma = fi, delta = deltaprime ) 

# ETERM1 est validé

ETERM2 <- list ( N = C2ii * IOM$N^2 , sigma = 2*fi, delta = 2 * deltaprime ) 

# ETREM2 est validé

ETERM3 <- list ( N = C2ik * outer(IOM$N, IOM$N), sigma = outer(fi, fi, '+'), delta = outer(deltaprime, deltaprime, '+'))

for (i in seq(length(IOM$N))) {
  for (k in seq(length(IOM$N))) {
     if ( i >= k ) {
       ETERM3$N[i,k] = NA; 
       ETERM3$sigma[i,k] = NA;
       ETERM3$delta[i,k] = NA}}}


keep <- which(!is.na(ETERM3$sigma))

ETERM3$N <- ETERM3$N[keep]
ETERM3$sigma <- ETERM3$sigma[keep]
ETERM3$delta <- ETERM3$delta[keep]


ETERM4 <- list ( N = C2ikprime * outer(IOM$N, IOM$N), sigma = outer(fi, fi, '-'), delta = outer(deltaprime, deltaprime, '-'))

for (i in seq(length(IOM$N))) {
  for (k in seq(length(IOM$N))) {
     if ( i >= k ) {
       ETERM4$N[i,k] = NA; 
       ETERM4$sigma[i,k] = NA;
       ETERM4$delta[i,k] = NA}}}


keep <- which(!is.na(ETERM4$sigma))

ETERM4$N <- ETERM4$N[keep]
ETERM4$sigma <- ETERM4$sigma[keep]
ETERM4$delta <- ETERM4$delta[keep]




OBLIQUITY <- list ( N = c(ETERM1$N, ETERM2$N, ETERM3$N, ETERM4$N), 
                    sigma = c(ETERM1$sigma, ETERM2$sigma, ETERM3$sigma, ETERM4$sigma), 
                    delta = c(ETERM1$delta, ETERM2$delta, ETERM3$delta, ETERM4$delta))


# so if we read Berger and Loutre 1991: 
# the "constant" of integration are epsilonbar and zeta

epsilonbarstar <- epsilonbar - sum(IOM$N*IOM$N*(0.5*ci*(ci-1)*tane + 0.25*(2*ci-1)*cote) ) 


# keep the first 200 terms

# we could also transform the negative $N$ into positive by changing the Phases (all cosines)

# again one step at a time. ... .

nkeep = 2000
OR <- order(abs(OBLIQUITY$N), decreasing=TRUE)
OR <- OR[seq(min(length(OR), nkeep))]
OBLIQUITY <- with(OBLIQUITY, data.frame(Amp = N[OR], Freq=sigma[OR], Phases=delta[OR]))



# with this, there seems be no duplicates (but maybe there will be once we have recitfied the N's)
# however we have a problem becase the sigma0 generates periods = psi

OBLIQUITY_ORDERED <- OBLIQUITY # %>% group_by(Freq) %>% summarise (Amp=sum(Amp), Phases)
#OBLIQUITY_ORDERED # <- OBLIQUITY %>% group_by(Freq) %>% summarise (Amp=sum(Amp), Phases)
ORDER <- order(abs(OBLIQUITY_ORDERED$Amp), decreasing=TRUE)
OBLIQUITY <- as.data.frame(OBLIQUITY_ORDERED[ORDER, ])

attr(OBLIQUITY, "class") <- 'spectral_decomp'
attr(OBLIQUITY, "trend") <- 0. 

# this needs to be updated
attr(OBLIQUITY, "shift") <- epsilonbarstar

which(duplicated(OBLIQUITY$sigma))

2*pi/ OBLIQUITY$sigma[seq(100)]

require(palinsol)
data(BER90)

if (1 == 0) { # disabled plotting : was for testing
  par(mfrow=c(1,2))
  plot(abs(OBLIQUITY$sigma)/(2*pi), OBLIQUITY$Amp, type='n', ylim=c(-1e-2, 0), xlim=c(0e-5, 4e-5))
  lines(ETERM1$sigma/(2*pi), ETERM1$N, col='red', type='h')
  lines(ETERM2$sigma/(2*pi), ETERM2$N, col='blue' , type='h')
  lines(ETERM3$sigma/(2*pi), ETERM3$N, col='gray', type='h')
  lines(ETERM4$sigma/(2*pi), ETERM4$N, col='green', type='h')
  grid()
  with(BER90$Table1,  plot(Freq/360/3600 , Amp*2*pi/60/60/360, type='h', col='black', ylim=c(-1e-2, 0), xlim=c(0e-5, 4e-5)))
  grid()
}
# ok this is more or less validated. 

#plot(ETERM1$sigma/(2*pi), C1f, type='h')

# attention at this point we haven't tried to reorganise and chase doublons. 
# there will be many of those, but one step at a time. 
# make all frequencies positive for obliquity


# this would be the command to group by frequency
# note you must also make an equivalence for sigma -> abs(sigma) but also correct for te Phases accordingly
# > OBLIQUITY_S <- OBLIQUITY %>% group_by(sigma) %>% summarise (N=sum(N))


####################################
# END OF THE OBLIQUITY 
###################################

################################################
# START  PRECESSION WITH RESPECT TO FIX PLANE
##############################################


PTERM1 <- list(  Amp = Dfseconde * outer(EPI$M,EPI$M), 
                       Freq = outer(EPI$g, EPI$g, '-'), 
                       Phases = outer(EPI$beta, EPI$beta, '-'))


PTERM2 <- list ( Amp = D1*IOM$N, Freq = fi, Phases = deltaprime)
PTERM3 <- list ( Amp = D2*IOM$N^2, Freq = 2*fi, Phases = 2*deltaprime)
PTERM4 <- list ( Amp = D2ik*outer(IOM$N, IOM$N), 
                       Freq = outer(fi, fi, '+'), 
                       Phases = outer(deltaprime, deltaprime, '+'))


PTERM5 <- list ( Amp = D2ikprime*outer(IOM$N, IOM$N), 
                       Freq = outer(fi, fi , '-'), 
                       Phases = outer(deltaprime, deltaprime, '-'))

# print('//////')
# print (PTERM5$Freq[26,61])
# print ('a')
# print ( which (abs(abs(PTERM5$Freq) - abs(2.437144e-6)) < 1.e-12))
# print ('b')
# print (PTERM5$Freq[27,62])
# print('//////')


for (i in seq(ng)) {
  for (k in seq(ng)) {
      if (k >= i)  {
        PTERM1$Amp[i, k] = NA;
        PTERM1$Freq[i, k ] = NA;
        PTERM1$Phases[i, k ] = NA;
      }
  }}

for (i in seq(nf)) {
  for (k in seq(nf)) {
      if (k >= i)  {
        PTERM4$Amp[i, k] = NA;
        PTERM4$Freq[i, k ] = NA;
        PTERM4$Phases[i, k ] = NA;
        PTERM5$Amp[i, k] = NA;
        PTERM5$Freq[i, k ] = NA;
        PTERM5$Phases[i, k ] = NA;
      }
  }}


keep <- which(!is.na(PTERM1$Freq))

PTERM1$Amp <- PTERM1$Amp[keep]
PTERM1$Freq <- PTERM1$Freq[keep]
PTERM1$Phases <- PTERM1$Phases[keep]


keep <- which(!is.na(PTERM4$Freq))


PTERM4$Amp <- PTERM4$Amp[keep]
PTERM4$Freq <- PTERM4$Freq[keep]
PTERM4$Phases <- PTERM4$Phases[keep]


keep <- which(!is.na(PTERM5$Freq))


PTERM5$Amp <- PTERM5$Amp[keep]
PTERM5$Freq <- PTERM5$Freq[keep]
PTERM5$Phases <- PTERM5$Phases[keep]

if (ber78_strict) {
 keep1 <- seq(min(200, length(PTERM1$Amp)))
 keep2 <- seq(min(200, length(PTERM2$Amp)))
 keep3 <- seq(min(200, length(PTERM3$Amp)))
 keep4 <- seq(min(200, length(PTERM4$Amp)))
 keep5 <- seq(min(200, length(PTERM5$Amp)))

PSI <- data.frame (
        Amp = c(PTERM1$Amp[keep1], PTERM2$Amp[keep2], PTERM3$Amp[keep3], PTERM4$Amp[keep4], PTERM5$Amp[keep5]), 
        Freq = c(PTERM1$Freq[keep1], PTERM2$Freq[keep2], PTERM3$Freq[keep3], PTERM4$Freq[keep4], PTERM5$Freq[keep5]), 
        Phases = c(PTERM1$Phases[keep1], PTERM2$Phases[keep2], PTERM3$Phases[keep3], PTERM4$Phases[keep4], PTERM5$Phases[keep5]), 
        Group = c(PTERM1$Amp[keep1]*0+1, PTERM2$Amp[keep2]*0+2, PTERM3$Amp[keep3]*0+3, PTERM4$Amp[keep4]*0+4, PTERM5$Amp[keep5]*0+5))

} else {



# print ('icicic')
# print (PTERM5$Freq[1764])
# print ('icicic 00000')

PSI <- data.frame (
        Amp = c(PTERM1$Amp, PTERM2$Amp, PTERM3$Amp, PTERM4$Amp, PTERM5$Amp), 
        Freq = c(PTERM1$Freq, PTERM2$Freq, PTERM3$Freq, PTERM4$Freq, PTERM5$Freq), 
        Phases = c(PTERM1$Phases, PTERM2$Phases, PTERM3$Phases, PTERM4$Phases, PTERM5$Phases), 
        Group = c(PTERM1$Amp*0+1, PTERM2$Amp*0+2, PTERM3$Amp*0+3, PTERM4$Amp*0+4, PTERM5$Amp*0+5))
}

negatives <- which(PSI$Freq < 0);

PSI$Amp[negatives] = -PSI$Amp[negatives]
PSI$Freq[negatives] = -PSI$Freq[negatives]
PSI$Phases[negatives] = -PSI$Phases[negatives] 
PSI$Phases = PSI$Phases  %% (2*pi)

# we proceed in two steps: first we order by frequencies, and check whether
# grouping is needed. 


# if aggregate (this can be time consuming, but this is the easiest to read, and
# may introduce some small numeric errors) And loose group information
# put efficient after all
 
if (aggregating){
COMPLEX_AMP <- PSI$Amp * cis(PSI$Phases)
tmp <- aggregate( Amp ~ Freq, data = data.frame(Amp= COMPLEX_AMP, Freq=signif(PSI$Freq, 12)), FUN=sum)
PSI <- data.frame(Amp = Mod(tmp$Amp), Freq=tmp$Freq, Phases=Arg(tmp$Amp))
}

ORDER <- order(abs(PSI$Amp), decreasing=TRUE)
PSI <- as.data.frame(PSI[ORDER, ])


# thing to look at
# With La88, the group1[205] and group5[1764] have exactly the same amplitudes
# and frequencies, but different phases. why ? 

# the 205 corresponds to the combination 3 / 51 (205=79+78+48, and 48+3 = 51)
# and tis corresponds to gi terms in the epi development

# the 1764 corresponds to the combination 26 / 61 
# and tis corresponds to gi terms in the IOM development

# it seems impossible to get them identical but yet this seems to be

# print ('<<<<<-')
# a1 = EPI$g[3] - EPI$g[51]
# a2 = IOM$sigma[27] - IOM$sigma[62]
# print (c(a1, a2, a1-a2))
# print ('<<<<<-')




attr(PSI, "class") <- 'spectral_decomp'
attr(PSI, "trend") <- psibar
attr(PSI, "shift") <- zeta

# note: there are 7 occurrences where same Freqs give different Phasess. 
# they come from the merging of 'g' terms with 's' terms


################################################
# END OF PRECESSION WITH RESPECT TO FIXED PLANE
##############################################

################################################
# START OF CLIMATIC PRECESSION WITH RESPECT TO FIXED PLANE
##############################################

CTERM1 <- list( Amp = EPI$M, Freq = EPI$g + psi, Phases = EPI$beta + zeta)
CTERM2 <- list( Amp = as.numeric(outer(EPI$M, D1 * IOM$N)), 
                Freq = outer (EPI$g, IOM$sigma, '+' ) + 2 * psi , 
                Phases = outer (EPI$beta, IOM$delta, '+' ) + 2 * zeta )

CTERM3 <- list( Amp = as.numeric(outer(EPI$M, D1 * IOM$N)), 
                Freq = outer (IOM$sigma, EPI$g, '-' ) + 2 * psi , 
                Phases = outer (IOM$delta, EPI$beta, '-' ) )



# require(palinsol)
# data(BER90)

# par(mfrow=c(1,2))
# BFreq <- BER90$Table5$Freq*(2*pi)/360/3600
# BAmp <- BER90$Table5$Amp*(2*pi)/360/3600

# plot(PSI_ORDERED$Freq, PSI_ORDERED$N, type='h')
# plot(BFreq, BAmp, type='h')



# ellradiants <- ell * (2*pi/3600/360)

CLIMPRECESS <- data.frame ( 
       Amp = c(CTERM1$Amp, CTERM2$Amp, CTERM3$Amp), 
       Freq = c(CTERM1$Freq, CTERM2$Freq, CTERM3$Freq), 
       Phases = c(CTERM1$Phases, CTERM2$Phases, CTERM3$Phases))


attr(CLIMPRECESS, "class") <- 'spectral_decomp'
attr(CLIMPRECESS, "trend") <- 0. 
attr(CLIMPRECESS, "shift") <- 0.


################################
# END OF CLIMATIC PRECESSION
#################################

### THE FOLLOWING CONSTRAINT SHOULD BE RESPECTED
### WHCIH REQUIRES ITERATING OVER either PSIBAR OR AL0
### psi <- ellradiants * cos(epsilonbar) * Constraint

# voir aussi ce que ca implique pour epsilon0 et zeta0

### BERGER USED A NEWTON RAPHSON METHOD WITH ANALYTICAL DERIVATIVE

# if we consider epsilonbar and ell to be the given, then loop by providing epsilonbar and ell until convergence
# etc. 

# Attention this is important because in the end this will impact the precession frequency

# depending on what is provided, we can iterate to make this relationship satisfied. 
                                              
# still need to clarify epsilonbar and the determination of initial conditions but we can leave it for later

class(CLIMPRECESS) <- "discreteSpectrum"
class(PSI) <- "discreteSpectrum"
class(OBLIQUITY) <- "discreteSpectrum"


 return(list(CLIMPRECESS=CLIMPRECESS, 
             PSI = PSI, 
             OBLIQUITY = OBLIQUITY))

}



#' Sharaf-Budnikova model, based on Berger's implementation
#'
#' This is a R-recoding of André Berger's implementation of the Sharak and Budnikova precession
#' model. Using, as inputs, a trigonometric expansion of (e,pi) and (i/2, omega) (e.g. `data(La88)`)
#' newcomb canstants (as, e.g., supplied by `newcomp_parameters`), mean obliquity and mean precession phase at
#' the reference year corresponding to the orbital solution (usually 1950.0 or 2000.0). 
#' Supplies trigonometrical expansions for obliquity (OBLIQUITY), longitude of the vernal equinox with respect to the fixed
#' plane (PSI) and climatic precession. 
#'
#' 
#' @param orbital solution (e.g., `data(La88)`)
#' @param P0  : partial derivative of `P0` with respect to eccentricity
#' @param ell : "newcomb" constant
#' @param epsilonbar : reference obliquity
#' @param zeta : precession phase at year zero (note toself: I need to be a bit more accurate about the meaning here"
#' @param aggregating : should different terms with same frequencies be aggregated ? (default : `TRUE`)
#' @param ber78_strict : original implementation of Berger 1973 thesis, used for the BER78 solution
#' @param ber90_reproduce : hack to actually reproduce the Berger and Loutre 1991 paper, which includes a change with respect to the description in the paper. 
#' @return a list with `discreteSpectrum` for OBLIQUITY, PSI (longitude of perihelion on fixed plane) and CLIMPRECESS
#'         a `discreteSpectrum` is a list with `Amp`, `Freq` (angluar velocity) and `Phases`, with attributes given the
#'         `shift` with respect to zero, and `trend` (for PSI). 
#' @note  Running this routine typically takes one second on a modern computer, and precession and obliquity may then
#'        be computed very efficiently for any period of time where the supplied parameters are considered as valid
#'        However, given the nature of orbital solutions and also the range of validity of the SB development (accurate
#'        up to order 2 in eccentricity), this generally would not span more than 2 or 3 million years. This would be 
#'        enough for the recent past, or for illustrative purposes in the deep past. Modern codes for precession, more
#'        accurate, would however not take much more computing time in practice. 
#'        if aggregating is "FALSE", then the group, as per the arithmetic expansion of Berger (1973) and Berger and Loutre (1991)
#'        belong, is given. 
#' @references Berger, A. L. (1973), Berger and Loutre (1991), Sharaf S. G. and N. A. Boudnikova (1967), Secular variations of elements of the Earth's orbit which influences the climates of the geological past, Tr. Inst. Theor. Astron. Leningrad, (11) 231-261 
#' @examples
#' 
#' data(La88)
#' S <- SB67(La88, aggregating = TRUE)
#' S2 <- SB67(La88, aggregating = FALSE)
#' 
#' psi1 <- develop(S$PSI, -1e6, 0, 1e3)
#' psi2 <- develop(S2$PSI, -1e6, 0, 1e3)
#' plot(psi1)
#' lines(psi1 - psi2)
#' @export
SB67 <- function (PlanetarySolution, 
                  P0 = 17.3919 *  pi/(3600*180.) ,
                  ell = 54.9066 * pi/(3600*180.) , 
                  epsilonbar = 23.399935 * pi/180., 
                  zeta = 1.600753 * pi / 180., aggregating = TRUE, ber78_strict = FALSE, ber90_reproduce = FALSE){

OUT <- SB67_Internal(ell, P0, epsilonbar, zeta, PlanetarySolution, aggregating = aggregating, ber78_strict = ber78_strict, ber90_reproduce = ber90_reproduce)

# le code commenté est pour deux types de reconstructions de la precession climatique. 1. Ec*sin(Pi+psi) ; 2. Par  le developpement. 

### times <- seq(1000)*1000
### esinpi <- reconstruct_spectral(EPPI, times, sin)
### ecospi <- reconstruct_spectral(EPPI, times, cos)
### Pi <- Arg(ecospi + 1i*esinpi)
### Ec <- Mod(ecospi + 1i*esinpi)
### 
### 
### psi <- reconstruct_spectral(OUT$PSI, times, sin) %% (2*pi)
### 
### plot(Pi+psi)
### 
### esinomega <- Ec * sin(Pi+psi)
### 
### plot(esinomega)
### 
### esinomega_method <- reconstruct_spectral(OUT$CLIMPRECESS, times, sin) 
### 
### 
### plot(esinomega, type='l')
### lines(esinomega_method, col='red')


# le mathch est correct mais pas parfait. 
# revoir avec les rotines originales de Berger pour voir si le match etait cense etre parfait. 
# on peut le faire avec palinsol, version deepast (pas deepast_test)

return(OUT)
}



# the  match is perfect down to 1.e-14 (for two points in a time series); so aggregation seems to work well, here. 




