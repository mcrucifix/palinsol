\name{LA04}
\alias{LA04}

\docType{data}
\title{Astronomical elements supplied by Laskar et al. 2004}
\description{
   Astronomical elements (longitude of perihelion, 
   obliquity and eccentricity) supplied by Laskar et al. 2004
   by stey of 1ka, from -51Ma to present (la04past) and from present to
   + 21Ma.
   }
\usage{ data(LA04) }
\format{\code{list}, named \code{LA04}, and containing two \code{data.frames}, 
        with the contents of original sources: \code{la04past} and 
        \code{la04future} for past and future astronomical elements, 
        respectively.}
\source{http://www.imcce.fr/Equipes/ASD/insola/earth/earth.ht.}
\references{
  J. Laskar et al., A long-term numerical solution for the 
  insolation quantities of the Earth,  Astron. Astroph., 428, 261-285  200 }

\author{Michel Crucifix, U. catholique de Louvain, Belgium.}

\examples{ 
data(LA04)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
