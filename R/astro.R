# Copyright (c) 2022 Michel Crucifix <michel.crucifix@uclouvain.be>

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject

# the following conditions:

# The above copyright notice and this permission notice shall be
# incluudedin all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND INFRINGEMENT
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR

# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# When using this package for actual applications, always
# cite the authors of the original insolation solutions
# Berger, Loutre and/or Laskar, see details in man pages


# ------------------------------------------------------------------
# R Code developed for R version 2.15.2 (2012-10-26) -- "Trick or Treat"
# ------------------------------------------------------------------

 ## calculates astronomical elements according to the solutions given below


#' Compute astronomical parameters in the past or in the future
#'
#' --
#'
#'
#' Both \code{ber78} and \code{ber90} compute astronomical elements based on a
#' spectral decomposition (sum of sines and cosines) of obliquity and planetary
#' precession parameters.  \code{ber78} uses the Berger (1978) algorithm and is
#' accurate for +/- 1e6 years about the present. \code{ber90} uses the Berger
#' and Loutre (1991) algorithm and is accurate for +/- 3e6 years about the
#' present (but with a tiny accuracy over the last 50 kyr, usually negligible
#' for any palaeo application, see example below).
#'
#' \code{la04} interpolates tables provided by Laskar (2004), obtained by a
#' simplectic numerical integration of the planetary system, in which the Moon
#' is considered as a planet.  This solution is valid for about 50 Myr around
#' the present.
#'
#' \code{precession}, \code{coprecession} and \code{obliquity} do as
#' \code{astro}, but only return precession (e sin varpi), coprecession (e cos
#' varpi) and obliquity, respectively.
#'
#' @aliases astro ber78 ber90 la04 precession coprecession obliquity
#' @param t Time, years after 1950
#' @param solution solution used. One of \code{ber78}, \code{ber90} or
#' \code{la04}
#' @param degree returns angles in degrees if \code{TRUE}
#' @return A vector of 3 (la04) or 4 (\code{ber78} and \code{ber90})
#' astronomical elements \tabular{ll}{ \code{eps} \tab obliquity, \cr
#' \code{ecc}\tab eccentricity and \cr \code{varpi} true solar longitude of the
#' perihelion. \cr } \code{ber78} and \code{ber90} also return \code{epsp}, the
#' Hilbert transform of obliquity (sines changed in cosines in the spectral
#' decomposition).
#'
#' Angles are returned in radians unless \code{degree=TRUE}
#' @author Michel Crucifix, U. catholique de Louvain, Belgium.
#' @references Berger, A. L. (1978).  Long-term variations of daily insolation
#' and Quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367,
#' doi:10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2
#'
#' Berger and M.F. Loutre (1991), Insolation values for the climate of the last
#' 10 million years, Quaternary Science Reviews, 10, 297 - 317,
#' doi:10.1016/0277-3791(91)90033-Q
#'
#' J. Laskar et al. (2004), A long-term numerical solution for the insolation
#' quantities of the Earth, Astron. Astroph., 428, 261-285,
#' doi:10.1051/0004-6361:20041335
#' @keywords misc
#' @examples
#'
#'
#' ## compare the obliquity over the last 2 Myr with the
#' ## three solutions
#'
#' times <- seq(-2e6,0,1e3)
#' Obl <- function(t) {c(time=t,ber78=ber78(t)['eps'],
#'        ber90=ber90(t)['eps'], la04=la04(t)['eps'])}
#'
#' Obls <- data.frame(t(sapply(times,Obl)))
#' ## may take about 10 seconds to run
#' with(Obls, {
#'   plot(times/1e3, ber78.eps, type='l', xlab='time (kyr)',
#'                                        ylab='Obliquity (radians)')
#'   lines(times/1e3, ber90.eps, type='l', col='red')
#'   lines(times/1e3, la04.eps, type='l', col='green')
#'   })
#'
#' legend('topright', c('ber78','ber90','la04'), col=c('black','red','green'),
#'        lty=1)
#'
#' ## same but with a zoom over the last 300 000 years:
#'
#' T <- which (times > -3e5)
#' with(Obls, {
#'   plot(times[T]/1e3, ber78.eps[T], type='l', xlab='time (kyr)',
#'                                        ylab='Obliquity (radians)')
#'   lines(times[T]/1e3, ber90.eps[T], type='l', col='red')
#'   lines(times[T]/1e3, la04.eps[T], type='l', col='green')
#'   })
#'
#' legend('topright', c('ber78','ber90','la04'), col=c('black','red','green'), 
#'        lty=1)
#'

#' @export astro
astro <- function(t, solution = ber78, degree = FALSE) {
  solution(t, degree)
}


#' @export ber78
ber78 <- function(t, degree = FALSE) {

  psibar <- 50.439273 / 60. / 60. * pi / 180
  estar <- 23.320556
  ##  e0    <- 0.028707
  zeta  <- 3.392506 * pi / 180.
  twopi <- 2*pi

  sectorad <- pi / (180*60.*60.)

  M <-  palinsol::BER78$Table4$Amp
  g <-  palinsol::BER78$Table4$Rate * sectorad
  b <-  palinsol::BER78$Table4$Phase * pi /180
  F <-  palinsol::BER78$Table5$Amp * sectorad
  fp <- palinsol::BER78$Table5$Rate * sectorad
  d <-  palinsol::BER78$Table5$Phase * pi /180.
  A <-  palinsol::BER78$Table1$Amp/60./60.
  f <-  palinsol::BER78$Table1$Rate * sectorad
  phi <- palinsol::BER78$Table1$Phase * pi/180.

  ## Obliquity

  eps <- estar + sum(A * cos(f * t + phi))
  epsp <- estar + sum(A * sin(f * t + phi))

  esinpi <- sum(M * sin(g * t+b))
  ecospi <- sum(M * cos(g * t+b))
  psi    <- psibar*t + zeta + sum(F * sin(fp * t + d))

  e <- sqrt(esinpi^2 + ecospi^2)
  Pi <- atan(esinpi/ecospi) + pi *( ecospi < 0)
  eps <- eps * pi / 180.
  epsp <- epsp * pi / 180.
  varpi <- (Pi+psi+pi) %% (twopi)

  if (degree) {
    rad2deg <- 180 / pi
    eps <- eps * rad2deg
    varpi <- varpi * rad2deg
  }

  c(eps = eps, ecc = e, varpi = varpi, epsp = epsp)
}
 ## Calculates climate orbital elements according to the algorithm given in A. Berger (1978)
 # Berger, A. L. (1978).  Long-term variations of daily insolation and Quaternary climatic changes,
 # J. Atmos. Sci., 35(), 2362-2367

 # This solution is valid for + / - 3e6 years.

## attach table to ber90 function
## Input :  t = time expressed in yr after 1950.0 (reference epoch)

#' @export ber90
ber90 <- function(t,degree=FALSE)
{

  psibar<- 50.41726176/60./60. * pi/180
  estar <- 23.33340950
  zeta  <- 1.60075265 * pi/180.
  twopi <- 2*pi

  sectorad <- pi/(180*60.*60.)

  M <-  palinsol::BER90$Table4$Amp
  g <-  palinsol::BER90$Table4$Rate*sectorad
  b <-  palinsol::BER90$Table4$Phase*pi/180
  F <-  palinsol::BER90$Table5$Amp*sectorad
  fp <- palinsol::BER90$Table5$Rate*sectorad
  d <-  palinsol::BER90$Table5$Phase*pi/180.
  A <-  palinsol::BER90$Table1$Amp/60./60.
  f <-  palinsol::BER90$Table1$Rate*sectorad
  phi<- palinsol::BER90$Table1$Phase*pi/180.

  ## Obliquity

  eps <- estar + sum(A*cos(f*t+phi))
  epsp <- estar + sum(A*sin(f*t+phi))

  esinpi <- sum(M*sin(g*t+b))
  ecospi <- sum(M*cos(g*t+b))
  psi    <- psibar*t + zeta + sum(F*sin(fp*t +d))

  e <- sqrt(esinpi^2+ecospi^2)
  Pi <-atan(esinpi/ecospi)+pi*(ecospi<0)
  eps <- eps * pi/180.
  epsp <- epsp * pi/180.
  varpi <- (Pi+psi+pi) %% (twopi)

  if (degree) {rad2deg <- 180/pi
               eps <- eps*rad2deg
               varpi <- varpi*rad2deg}

  c(eps=eps,ecc=e,varpi=varpi,epsp=epsp)
}

 # This solution is valid for 50e6 years

## Input :  t = time expressed in yr after 1950.0 (reference epoch)

#' @export la04
la04 <- function(t, degree = FALSE)
{

  tka = t/1000. - 0.050  # time elapsed since 1950.0
  if (tka>0)
  {
   local(
   {
    F <-  floor(tka)
    ORB <<- palinsol::LA04$la04future[F+1, ]
    if (! (tka == F)) {
      D  <- tka - floor(tka)
      diff <- palinsol::LA04$la04future[F+2, ] - ORB
      # note : if the diff in varpi is greater than pi,
      # this probably means that we have skipped 2*pi,
      # so we need to correct accordingly
      if (diff$varpi > pi) diff$varpi = diff$varpi - 2*pi
      if (diff$varpi < -pi) diff$varpi = diff$varpi + 2*pi
      #
      ORB <<- ORB + D*diff
    }
   }
   )
  } else {
   local(
   {
    F <-  floor(tka)
    ORB <<- palinsol::LA04$la04past[-F+1, ]
    if (! (tka == F)) {
      D  <- tka - F

      diff <- palinsol::LA04$la04past[-F, ] - ORB
      # note : if the diff in varpi is greater than pi,
      # this probably means that we have skipped 2*pi,
      # so we need to correct accordingly
      if (diff$varpi > pi) diff$varpi = diff$varpi - 2*pi
      if (diff$varpi < -pi) diff$varpi = diff$varpi + 2*pi
      #
      ORB <<- ORB + D*diff
    }
   }
  )
  }

  if (degree) {
    rad2deg <- 180 / pi
    ORB['eps'] <- ORB['eps'] * rad2deg
    ORB['varpi'] <- ORB['varpi'] * rad2deg
  }


   # must return a array (0.92 -> 0.93)
   names <- c('eps','ecc','varpi')
   OUT = as.numeric(ORB[names]) ; names(OUT) <- names
   OUT
}

#' @export precession
precession <- function(t,solution=ber78)
##  as astro, but returns only precession parameter e sin (varpi)
{
  O <- astro(t,solution, degree=FALSE)
  as.numeric(O['ecc'] * sin (O['varpi']))
}

#' @export coprecession
coprecession <- function(t,solution=ber78)
##  as astro, but returns only precession parameter e sin (varpi)
{
  O <- astro(t,solution, degree=FALSE)
  as.numeric(O['ecc'] * cos (O['varpi']))
}

#' @export obliquity
obliquity <- function(t,solution=ber78,degree=FALSE)
##  as astro, but returns only obliquity
{ as.numeric(astro(t,solution, degree=degree)['eps']) }
