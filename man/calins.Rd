\name{calins}
\alias{calins}
 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Caloric insolation}

\description{
Computes caloric summer insolation for a given astronomical configuration and 
latitude. }

\usage{
 calins (orbit,lat=65*pi/180,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orbit}{Output from a solution, such as \code{ber78}, \code{ber90} or \code{la04}}
  \item{lat}{latitude}
  \item{...}{Other arguments passed to Insol}
}
\details{ The caloric summer is a notion introduced by M. Milankovitch. It is defined
as the  halve of the tropical year
during for which daily mean insolation are greater than all days of the other halves.
The algorithm is an original algorithm by M. Crucifix, but consistent with 
earlier definitions and algorithms by A. Berger (see examples). Do not confuse
this Berger (1978) reference with the Berger (1978), J. Atm. Sci. of the 
astronomical solution. }

\value{
  Time-integrated insolation in kJ/m2 during the caloric summer. 
  }

\references{
Berger (1978) Long-term variations of caloric insolation resulting from the earth's orbital elements,  Quaternary Research, 9, 139 - 167. 
}
  
\author{Michel Crucifix, U. catholique de Louvain, Belgium.}

\examples{
## reproduces Table 2 of Berger 1978
lat <- seq(90, 0, -10) * pi/180. ## angles in radiants. 
orbit_1 = ber78(0)
orbit_2 = orbit_1
orbit_2 ['eps'] = orbit_2['eps'] + 1*pi/180.

T <-  sapply(lat, function(x) c(lat = x * 180/pi, 
                          calins(orbit_2, lat=x, S0=1365) / (4.18 * 1e1)
                        - calins(orbit_1, lat=x, S0=1365) / (4.18 * 1e1) ) )
data.frame(t(T))
# there are still some differences, of the order of 0.3 %, that are probably related to
# the slightly different methods. 
# 41.8 is the factor from cal/cm2 to  kJ/m2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

