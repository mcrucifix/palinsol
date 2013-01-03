\name{PSeries}
\alias{PSeries}
\alias{Table2}
 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes incoming solar radiation (insolation)}

\description{
Reconstructs Table2 of Berger (1978), containing the explicit d'Alembert 
decomposotion of e sin w, where e is eccentricity and w, the longitude
of the perihelion, according to a procedure described in 
Berger and Loutre (1990) }

\usage{
data(PSeries)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\details{
  }

\value{
  Results are stored in Table2
  }

\references{
  Berger, A. L. (1978).  Long-term variations of daily
insolation and Quaternary climatic changes, J. Atmos. Sci., 35,
2362-2367.
  A. Berger and M. F. Loutre (1990), Origine des fr\'equences des \'el\'ements astronomiques intervenant dans l'insolation,  Bull. Classe des Sciences, 1-3, 45-106 
}
  
\author{Michel Crucifix, U. catholique de Louvain, Belgium.}

\examples{

## make a little wrapper, with all default values

data(PSeries)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
