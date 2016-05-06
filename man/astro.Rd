\name{astro}
\alias{astro}
\alias{ber78}
\alias{ber90}
\alias{la04}
\alias{precession}
\alias{obliquity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute astronomical parameters in the past or in the future}

\description{
--
}
\usage{
astro(t,solution,degree=FALSE)
ber78(t,degree=FALSE)
ber90(t,degree=FALSE)
la04(t,degree=FALSE)
precession(t,solution)
obliquity(t,solution,degree=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{Time, years after 1950}
  \item{solution}{solution used. One of \code{ber78}, \code{ber90} or \code{la04}}
  \item{degree}{returns angles in degrees if \code{TRUE}}
}
\details{
  
  Both \code{ber78} and  \code{ber90}  compute astronomical elements based on a spectral
  decomposition (sum of sines and cosines) of obliquity and planetary precession parameters.
  \code{ber78} uses the Berger (1978) algorithm and is accurate for +/- 1e6 years about the
  present. \code{ber90} uses the Berger and Loutre (1991) algorithm and is accurate
  for  +/- 3e6 years   about the present (but with a tiny accuracy over the last 50 kyr, 
  usually negligible for any palaeo application, see example below). 

  \code{la04}  interpolates tables provided by Laskar (2004), obtained by a simplectic
  numerical integration of the planetary system, in which the Moon is considered as a planet. 
  This solution is valid for about 50 Myr around the present. 

  \code{precession} and \code{obliquity} do as \code{astro}, but only return
  precession (e sin varpi) parameter and obliquity, respectively.
}

\value{
 A vector of 3 (la04) or 4 (\code{ber78} and \code{ber90}) astronomical elements
 \tabular{ll}{
   \code{eps} \tab obliquity, \cr
   \code{ecc}\tab  eccentricity and \cr \code{varpi} true solar longitude of the perihelion. \cr
  }
  \code{ber78} and \code{ber90} also return 
  \code{epsp}, the Hilbert transform of obliquity (sines changed
  in cosines in the spectral decomposition). 

  Angles are returned in radians unless \code{degree=TRUE}
}

\references{
  Berger, A. L. (1978).  Long-term variations of daily
insolation and Quaternary climatic changes, J. Atmos. Sci., 35,
2362-2367, doi:10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2

  Berger and M.F. Loutre (1991), Insolation values for the climate of the last 10 million years,  Quaternary Science Reviews, 10, 297 - 317, doi:10.1016/0277-3791(91)90033-Q

 J. Laskar et al. (2004), A long-term numerical solution for the insolation quantities of the Earth,  Astron. Astroph., 428, 261-285, doi:10.1051/0004-6361:20041335}
  
\author{Michel Crucifix, U. catholique de Louvain, Belgium.}

\examples{

## compare the obliquity over the last 2 Myr with the three solutions

times <- seq(-2e6,0,1e3)
Obl <- function(t) {c(time=t,ber78=ber78(t)['eps'],
       ber90=ber90(t)['eps'], la04=la04(t)['eps'])}

Obls <- data.frame(t(sapply(times,Obl)))
## may take about 10 seconds to run
with(Obls, {
  plot(times/1e3, ber78.eps, type='l', xlab='time (kyr)', 
                                       ylab='Obliquity (radians)')  
  lines(times/1e3, ber90.eps, type='l', col='red')  
  lines(times/1e3, la04.eps, type='l', col='green')  
  })

legend('topright', c('ber78','ber90','la04'), col=c('black','red','green'), lty=1)

## same but with a zoom over the last 300 000 years:

T <- which (times > -3e5)
with(Obls, {
  plot(times[T]/1e3, ber78.eps[T], type='l', xlab='time (kyr)', 
                                       ylab='Obliquity (radians)')  
  lines(times[T]/1e3, ber90.eps[T], type='l', col='red')  
  lines(times[T]/1e3, la04.eps[T], type='l', col='green')  
  })

legend('topright', c('ber78','ber90','la04'), col=c('black','red','green'), lty=1)




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
