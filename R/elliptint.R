 SIDERAL_YEAR = 365.25636 
 TROPIC_YEAR  = 365.24219876
 YEAR = SIDERAL_YEAR


W  <- function (phi, eps, ecc, lambda,S0=1365,n=3)
   {
 
    require(gsl) ## gnu scientific library ported by Robin Hankin. Thank you Robin !!! 
    pi2=pi/2
    H00   = pi2
    seps = sin(eps)
    ceps = cos(eps)
    sphi = sin(phi)
    cphi = cos(phi)
    tphi = sphi/cphi

    eq34 <- function (lambda)
    {
    slambda = sin(lambda)
    clambda = cos(lambda)
    k   = seps/cphi
    sesl   = seps*slambda
    tdelta = sesl / sqrt(1-sesl*sesl)
    H0    = acos(-tphi*tdelta)
    Flk   = ellint_F(lambda,k)
    Elk   = ellint_E(lambda,k)

    sphi*seps*(H00-clambda*H0) + cphi*Elk +
        sphi*tphi*Flk - sphi*tphi*ceps*ceps*
        ellint_P(lambda, k, -seps*seps)
    }
#

     eq40  <- function(lambda)
    {
    ## the max(-1, min(1, 
    ## is to account for a numerical artefact when lambda = lambda1,2,3,4

    slambda = sin(lambda)
    clambda = cos(lambda)

    k =  seps/cphi
    k1 = cphi/seps
    psi  = asin(max(-1,min(1,k*slambda)))

    if (clambda < 0) psi = pi-psi
     
    sesl   = seps*slambda
    tdelta = sesl / sqrt(1-sesl*sesl)
    H0    = acos(max(-1,min(1,-tphi*tdelta)))
    Fpk   = ellint_F(psi,k1)
    Epk   = ellint_E(psi,k1)
    Pipk  = ellint_P(psi,k1, -cphi*cphi)

    ( sphi*seps*(H00-clambda*H0) + seps*Epk +
        ceps* ceps/seps * Fpk - sphi*sphi*ceps*ceps/seps*
        Pipk)
    }
#

    eq38  <- function(lambda) { - pi * sphi*seps*cos(lambda) }


    T   = YEAR * 0.086400 * 1000
    xes = sqrt(1-ecc*ecc)
    W0  = S0*T/(2*pi*pi*xes)

    if (phi >= (pi2-eps) | phi <= -(pi2-eps) )
    {
    ## above polar circle 
    lambda1 = (asin(cphi/seps) )
    lambda2 = pi - lambda1
    lambda3 = pi + lambda1
    lambda4 = 2*pi - lambda1
  
    WW=0

    if (lambda > 0)         WW = WW + eq40(min(lambda,lambda1))
    if (lambda > lambda2)   WW = WW + eq40(min(lambda,lambda3)) - eq40(lambda2)
    if (lambda > lambda4)   WW = WW + eq40(lambda) - eq40(lambda4)

    if (phi >= (pi2-eps) ) {  ## northern hemisphere
      if (lambda > lambda1)   WW = WW + eq38(min(lambda,lambda2)) - eq38(lambda1)
      } else 
    { ## Southern hemisphere
      if (lambda > lambda3)   WW = WW + eq38(min(lambda,lambda4)) - eq38(lambda3)
    }

    WW = W0*WW 
    } else ## outside polar circle
    { WW = W0 * eq34(lambda) }
   WW
   }

