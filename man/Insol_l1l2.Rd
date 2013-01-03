\name{Insol_l1l2}
\alias{Insol_l1l2}
\alias{Insol_d1d2}
 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Time-integrated insolation}

\description{
Computes time-integrated incoming solar radiation (Insol) either between 
given true solar longitudes (\code{Insol_l1l2}) or days of year (\code{Insol_d1d2})
for a given orbit and latitude}

\usage{
Insol_l1l2 (orbit,l1=0,l2=2*pi,lat=65*pi/180,avg=FALSE,ell=TRUE,...)
Insol_d1d2 (orbit,d1,d2,lat=65*pi/180,avg=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orbit}{Output from a solution, such as \code{ber78}, \code{ber90} or \code{la04}}
  \item{lat}{latitude}
  \item{l1}{lower true solar longitude bound of the time-integral}
  \item{l2}{upper true solar longitude bound of the time-integral}
  \item{d1}{lower calendar day (360-day-year) of the time-integral}
  \item{d2}{upper calendar day (360-day-year) of the time-integral}
  \item{avg}{performs a time-average.}
  \item{ell}{uses elliptic integrals for the calculation (much faster)}
  \item{...}{other arguments to be passed to \code{Insol}}
}
\details{
  All angles input measured in radiants.  

  Note that in contrast to Berger (2010)
  we consider the tropic year as the reference, rather than the sideral year, 
  which partly explains some of the small differences with the original publication}

\value{
  Time-integrated insolation in kJ/m2 if \code{avg=TRUE}, else time-average in W/m2
  }

\references{
Berger, A., Loutre, M.F. and Yin Q. (2010), Total irradiation during any time interval of the year using elliptic integrals,  Quaternary Science Reviews, 29, 1968 - 1982.
}
  
\author{Michel Crucifix, U. catholique de Louvain, Belgium.}

\examples{
## reproduces Table 1a of Berger et al. 2010:
lat <- seq(85, -85, -10) * pi/180. ## angles in radiants. 
orbit=c(eps=  23.446 * pi/180., ecc= 0.016724, varpi= (102.04 - 180)* pi/180. )
T <-  sapply(lat, function(x) c(lat = x * 180/pi, 
        m1 =  Insol_l1l2(orbit, 0, 70 * pi/180, lat=x, ell= TRUE, S0=1368) / 1e3,
        m2 =  Insol_l1l2(orbit, 0, 70 * pi/180, lat=x, ell=FALSE, S0=1368) / 1e3) ) 
data.frame(t(T))
 ## reproduces Table 1b of Berger et al. 2010:
lat <- c(85, 55, 0, -55, -85) * pi/180. ## angles in radiants. 
T <-  sapply(lat, function(x) c(lat = x * 180/pi, 
         m1 =  Insol_l1l2(orbit, 30 * pi/180. , 75 * pi/180, 
               lat=x, ell= TRUE, S0=1368) / 1e3,
         m2 =  Insol_l1l2(orbit, 30 * pi/180. , 75 * pi/180, 
               lat=x, ell=FALSE, S0=1368) / 1e3) ) 
 ## reproduces Table 2a of Berger et al. 2010:
lat <- seq(85, -85, -10) * pi/180. ## angles in radiants. 

## 21 march in a 360-d year. By definition : day 80 = 21 march at 12u
d1 = 79.5 
d2 = 79.5 + (10 + 30 + 30 ) * 360/365.2425 ## 30th May in a 360-d year

T <-  sapply(lat, function(x) c(lat = x * 180/pi, 
        m1 =  Insol_d1d2(orbit, d1,d2, lat=x, ell= TRUE, S0=1368) / 1e3,
        m2 =  Insol_d1d2(orbit, d1,d2, lat=x, ell= FALSE, S0=1368) / 1e3))
                          
## I did not quite get the same results as on the table 
## on this one; probably a matter of calendar
## note : the authors in fact used S0=1368 (pers. comm.) 
## 1366 in the paper is a misprint

data.frame(t(T))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

