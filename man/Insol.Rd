\name{Insol}
\alias{Insol}
 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes incoming solar radiation (insolation)}

\description{
Computes incoming solar radiation (insolation) for a given astronomical configuration, 
 true solar longitude and latitude}

\usage{
Insol (orbit,long=pi/2, lat=65*pi/180,S0=1365)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orbit}{Output from a solution, such as \code{ber78}, \code{ber90} or \code{la04}}
  \item{long}{true solar longitude}
  \item{lat}{latitude}
  \item{S0}{Total solar irradiance}
}
\details{
  True solar longitude is measured in radiants:
  \tabular{ll}{
    pi/2 \tab for summer solstice\cr
    pi \tab   for autumn equinox\cr
    3 * pi/2 \tab for winter solstice\cr
    0 \tab for spring equinox\cr
  }
  It may be obtained for a given day in the year using the function \code{day2l}.
}

\value{
  Daily-mean insolation (assuming fixed astronomical parameters during a true solar day)
  }

\references{
  Berger, A. L. (1978).  Long-term variations of daily
insolation and Quaternary climatic changes, J. Atmos. Sci., 35,
2362-2367.
}
  
\author{Michel Crucifix, U. catholique de Louvain, Belgium.}

\examples{

## make a little wrapper, with all default values

insolation <- function(times, astrosol=ber78,...)
  sapply(times, function(tt) Insol(orbit=astrosol(tt)))

tts <- seq(from = -400e3, to = 0, by = 1e3)
isl <- insolation(tts, ber78)
plot(tts, isl, typ='l')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
