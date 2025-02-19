# thi actually provides the eccentricity spectrum of 
# m(1+A/2-A^2/8+A^3/16) where A= e^2/m^2 - 1 , and
# the definition of m given below. 
# this confirms the accuracy of the development given in Berger an Loutre 1990 
# and the fact that the square root has been developed up to the third order in (1-e^2/m). 
# errors are thus of the order of 5/128A^4, which perhaps surprisingly is still of the order of 0.015 / 0.020 in eccentricity. Errors are seem to be always positive (5/128*A^4 may reach 2.0, as A may reach 2.0, but not easy to prove because the next term, in A^5, can be large as well. 

# in fact using higher order terms helps a bit for capturing the 'low eccentricity' bit, but it actually makes things much worse above 0.04, so that the third order is actually the right deal. 


#' Trigonometrical development of eccentricity spectrum, based on an erbital solution
#'
#' This is a fully transparent R-implementation of the trigonometric epansion of excentricity,
#' based on a trigonometrical epansio of (e,pi) (e.g., `data(La88)`).
#' order  by decreasing frequency
#'
#' @param  G : `discreteSpectrum` object giving the 
#'             trigonometric development of (e * cos (Pi) + \sqrt(1) * e * sin (Pi)) 
#'             with complex amplitudes and frequencies as radians / year
#' @return `discreteSpectrum` object containing  the trigonometric development of eccentricity
#' @references Berger A. and M. F. Loutre (1990), Origine des frequences des elements astronomiques intervenant dans le calcul de l'insolation, Bulletin de la Classe des sciences, (1) 45-106 doi:10.3406/barb.1990.38523
#' @importFrom gtools combinations
#' @importFrom dplyr  bind_rows
#' @importFrom stats aggregate
#' @note # this actually provides the (exact) spectrum of m(1+A/2-A^2/8+A/16) where A= e^2/m^2 - 1 , and m =  sqrt(sum(epi$Amp^2))
#
eccentricity_spectrum <- function(G) {

require(gtools)

# the file EPI containts the development of "e sin Pi" which is the 
# only file needed to produce the development of eccentricity
# This particular file comes from La93.

## Il y asans doute une erreur dans la page 54 de Berger et Loutre, mais ou ???

## Revoir la these de Berger ? 

print('init tables')

EPI <- G

EPI$omega <- as.data.frame(G)[, 'Freq']
EPI$M <- as.data.frame(G)[, 'Amp']
EPI$phi <- as.data.frame(G)[, 'Phases']



combine <- function(a,b, om1, om2, ph1, ph2, factor=1) {
  A <- factor * outer(a,b, "*")
  O <- outer(om1,om2, "+")
  P <- outer(ph1,ph2, "+")
  A <- A[upper.tri(A)]
  O <- O[upper.tri(O)]
  P <- P[upper.tri(P)]
  return(data.frame(A=A,O=O,P=P))
}

combine2 <- function(a,b, om1, om2, ph1, ph2) {
IJ = gtools::combinations(length(a),2)
  A <- apply(IJ, 1, function(V) a[V[1]]*b[V[2]])
  O <- apply(IJ, 1, function(V) om1[V[1]]+om2[V[2]])
  P <- apply(IJ, 1, function(V) ph1[V[1]]+ph2[V[2]])
  return(data.frame(A=A,O=O,P=P))
}




# order  by decreasing frequency

O1 <- order(EPI$omega, decreasing=TRUE)

EPI <- list(omega=EPI$omega[O1], M=EPI$M[O1], phi=EPI$phi[O1])

m <- sqrt(sum(EPI$M^2))
a = EPI$M / m


n=0


B <- combine(a,a, EPI$omega, -EPI$omega,  EPI$phi, -EPI$phi)

bk <- B$A
gk <- B$O
pk <- B$P

# bk <- EPI$M[IJB[,1]] * EPI$M[IJB[,2]]
# gk <- EPI$omega[IJB[,1]] - EPI$omega[IJB[,2]]
# pk <- EPI$phi[IJB[,1]] - EPI$phi[IJB[,2]]

ek1 <- bk 
ek2 <- 0.375*bk*bk*bk 
sumbl2 <- rep(sum(bk^2 ), length(bk))
sumbl2 <- sumbl2 - bk*bk 

ek3 <- 0.75 * bk*sumbl2

ek = ek1+ek2+ek3
print('group 1 : termes d ordre 1')

group <- list()

group[[1]] <- data.frame(A=1-0.25*sum(bk^2 ), O=0, P=0)



group[[2]] <- data.frame(A=ek, O=gk, P=pk)
#  group2 <- -0.25 * data.frame(A=bk^2, O=2*gk, P=2*pk)
#  group5 <- 0.125 * data.frame(A=bk^3, O=3*gk, P=3*pk)

 odecreasing <- order(bk, decreasing=TRUE)
 bk <- bk[odecreasing]#[seq(100)]
 gk <- gk[odecreasing]#[seq(100)]
 pk <- pk[odecreasing]#[seq(100)]

# print('groupes 3 a 9 : termes d ordre 2')

# n=0
# IJ = matrix(0, 124750,2)
# for (i in seq(499)) for (j in seq(i+1,500))  {n=n+1 ; IJ[n,]=c(i,j)}

# bk3 = -0.5 * (bk[IJ[,1]] * bk[IJ[,2]] )
# gk3a = (gk[IJ[,1]] + gk[IJ[,2]] )
# gk3b = (gk[IJ[,1]] - gk[IJ[,2]] )

# pk3a = (pk[IJ[,1]] + pk[IJ[,2]] )
# pk3b = (pk[IJ[,1]] - pk[IJ[,2]] )

# group3 =  data.frame(A=c(bk3,bk3), O=c(gk3a,gk3b), P=c(pk3a,pk3b))



# bk6 = 0.375 * (bk[IJ[,1]] * bk[IJ[,2]]^2 )
# bk7 = 0.375 * (bk[IJ[,1]]^2 * bk[IJ[,2]] )

# gk6a = (gk[IJ[,1]] + 2*gk[IJ[,2]] )
# gk6b = (-gk[IJ[,1]] + 2*gk[IJ[,2]] )

# pk6a = (pk[IJ[,1]] + 2*pk[IJ[,2]] )
# pk6b = (-pk[IJ[,1]] + 2*pk[IJ[,2]] )

# gk7a = (2*gk[IJ[,1]] + gk[IJ[,2]] )
# gk7b = (2*gk[IJ[,1]] - gk[IJ[,2]] )

# pk7a = (2*pk[IJ[,1]] + pk[IJ[,2]] )
# pk7b = (2*pk[IJ[,1]] - pk[IJ[,2]] )


# group6 =  data.frame(A=c(bk6,bk6), O=c(gk6a,gk6b), P=c(pk6a,pk6b))
# group7 =  data.frame(A=c(bk7,bk7), O=c(gk7a,gk7b), P=c(pk7a,pk7b))

group[[3]] <- data.frame(A=-0.25 * (bk^2 ), O=2*gk, P=2*pk)
group[[4]] <-  combine(bk, bk, gk, gk, pk, pk, -0.5)
group[[5]] <-  combine(bk,bk, gk, -gk, pk, -pk, -0.5)
group[[6]] <-  data.frame(A=0.125 * bk^3, O=3*gk, P=3*pk)

group[[7]] <-   combine(bk^2, bk, 2*gk, gk, 2*pk, pk, 0.375)
group[[8]] <-  combine(bk, bk^2, gk, 2*gk, pk, 2*pk, 0.375)
group[[9]] <-  combine(bk^2, bk, 2*gk, -gk, 2*pk, -pk, 0.375)
group[[10]] <-  combine(bk, bk^2, -gk, 2*gk, -pk, 2*pk, 0.375)


maxterms = 121

#IJK = gtools::combinations(50,3)
IJK = gtools::combinations(min(maxterms,length(bk)),3)

# attention le signe moins n'est pas dans Berger

bklm =  0.75 * as.numeric(apply(IJK, 1, function(V)  {bk[V[1]] *  bk[V[2]] * bk[V[3]]}))
gklm1 = as.numeric(apply(IJK, 1, function(V) {gk[V[1]] +  gk[V[2]] + gk[V[3]]}))
gklm2 = as.numeric(apply(IJK, 1, function(V) {gk[V[1]] +  gk[V[2]] - gk[V[3]]}))
gklm3 = as.numeric(apply(IJK, 1, function(V) {gk[V[1]] -  gk[V[2]] + gk[V[3]]}))
gklm4 = as.numeric(apply(IJK, 1, function(V) {gk[V[1]] -  gk[V[2]] - gk[V[3]]}))
pklm1 = as.numeric(apply(IJK, 1, function(V) {pk[V[1]] +  pk[V[2]] + pk[V[3]]}))
pklm2 = as.numeric(apply(IJK, 1, function(V) {pk[V[1]] +  pk[V[2]] - pk[V[3]]}))
pklm3 = as.numeric(apply(IJK, 1, function(V) {pk[V[1]] -  pk[V[2]] + pk[V[3]]}))
pklm4 = as.numeric(apply(IJK, 1, function(V) {pk[V[1]] -  pk[V[2]] - pk[V[3]]}))

group[[11]] =  data.frame(A=c(bklm,bklm,bklm,bklm), 
                      O=c(gklm1,gklm2,gklm3,gklm4),
                      P=c(pklm1,pklm2,pklm3,pklm4))

print('developpement final')

# developpement <- rbind(group0, group1, group2, group3, group4,  group5, group6, group7, group8, group9,  group10)

# print(length(developpement$A))

# developpement <- group10
# ici on prend un petit risque, on ne va garder que les termes dont l'amplitude est > 5.e-7, soit environ 1000 termes. 

print('ecremage')

#group <- lapply(group, function(developpement) developpement[which(abs(developpement$A) > 1.e-6),])

print ('generate partial reconstructions')

times <- seq(0,5e3,2)*1e3
X <- t(sapply(group, function(developpement)
   sapply(times, function(t) sum(m * developpement$A * cos(developpement$O*t+developpement$P))) ))

esinpi <- sapply(times, function(t) sum(EPI$M * sin(EPI$omega*t+EPI$phi)))
ecospi <- sapply(times, function(t) sum(EPI$M * cos(EPI$omega*t+EPI$phi)))
Y = sqrt(esinpi*esinpi+ecospi*ecospi)



XX <- apply(X,2,sum)

Xts <- ts(XX, deltat=times[2]-times[1])
Yts <- ts(Y, deltat=times[2]-times[1])


# thi actually provides the eccentricity spectrum of 
# m(1+A/2-A^2/8+A/16) where A= e^2/m^2 - 1 , and
# the definition of m given above

A = (Yts/m)^2 - 1;
Yts2 <- m*(1 + A/2 - A*A/8 + A*A*A/16)

Yts3 <- m*(1 + A/2 - A*A/8 + A*A*A/16 - 5*A^4/128)
Yts4 <- m*(1 + A/2 - A*A/8 + A*A*A/16 - 5*A^4/128 + 7*A^5/256)

# now you can check that Xts and Yts2 match perfectly, to 1e-16


print ('chasse aux doublons')


developpement <- bind_rows(group) 

developpement_1 <- data.frame(A=developpement$A,  O=developpement$O)
developpement_2 <- data.frame(P=developpement$P,  O=developpement$O)

AFF1 <- aggregate(.~O, data=developpement_1, FUN=sum)
AFF2 <- aggregate(.~O, data=developpement_2, FUN=function(x) x[1])

AFF <- cbind(AFF1, P=AFF2$P)


print ('liste finale:')

negative_amplitudes <- which(AFF$A < 0)

AFF$A[negative_amplitudes] <-      - AFF$A[negative_amplitudes]
AFF$O[negative_amplitudes] <-      - AFF$O[negative_amplitudes]
AFF$P[negative_amplitudes] <-   pi - AFF$P[negative_amplitudes]


o = order(AFF$A, decreasing=TRUE)



out <- list (Freq= AFF$O[o], Amp = AFF$A[o]*m, Phases = AFF$P[o])
class(out) <- "discreteSpectrum"

return(out)
}
#print(length(AFF$A))


# AFF$A = AFF$A*m
# developpement$A = developpement$A*m


##

#### print('ici')
#### print(length(AFF$A))
#### 
#### print('reconstruction...')
#Xts2 <- ts (  sapply(times, function(t) sum(AFF$A * cos(AFF$O*t+AFF$P))), deltat = times[2]-times[1]) * m



# Xts and Xts2 also match perfectly 




#### print('palinsol.')
#### ecc_ber90 <- sapply(times, ber90)['ecc',]
#### 
#### 
#### 
#### 
#### 
#### plot(ecc_reconstruct, type='l')
#### # lines(ecc_reconstruct_non_aggregated, type='l', col='blue')
#### lines(ecc_ber90-0.02, col='red')
#### 
#### plot(ecc_reconstruct - ecc_ber90, type='l')
#### 
#### print('true variances')
#### # print(var(ecc_reconstruct_non_aggregated))
#### print(var(ecc_reconstruct))
#### print(var(ecc_ber90))
#### print(var(ecc_ber90 - ecc_reconstruct))
#### 
#### # Order1 <- order(AFF$O, decreasing=FALSE)
#### # Order2 <- order(AFF$A, decreasing=TRUE)
#### # 
#### # plot(2*pi/AFF$O[Order1], AFF$A[Order1], type='l', xlim=c(0,3e6))
#### # AFF$period = 2*pi/AFF$O
#### # 
#### # AFF[Order2[seq(10)],]
#### # 
#### 
