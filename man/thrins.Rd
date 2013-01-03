\name{thrins}
\alias{thrins}
 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Integrated insolation for all days exceeding a threshold}

\description{
Integrated insolation over the part during which daily-mean insolation 
exceeds a threshold, expressed in W/m2}

\usage{
 thrins  (lat=65*pi/180,orbit,threshold=400,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lat}{latitude}
  \item{orbit}{Output from a solution, such as \code{ber78}, \code{ber90} or \code{la04}}
  \item{threshold}{threshold insolation ,in W/m2}
  \item{...}{other arguments to be passed to Insol}
}
\details{ Algorithm is by M. Crucifix, but the idea of thresholded insolation is due to Huybers and 
Tziperman (2008), reference below.}

\value{
  Time-integrated insolation in kJ/m2 . The quantity is calculated by brute-force integration
  with a 1-degree time-step in true solar longitude and this can be quite slow if long series
  are to be calculated. 
  }

\references{
P. Huybers and E. Tziperman (2008), Integrated summer insolation forcing and 40,000-year glacial cycles: The perspective from an ice-sheet/energy-balance model,  Paleoceanography, 23.
}
  
\author{Michel Crucifix, U. catholique de Louvain, Belgium.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

