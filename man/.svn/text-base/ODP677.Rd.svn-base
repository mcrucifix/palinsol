\name{ODP677}
\alias{ODP677}
\docType{data}
\title{Marine core ODP677}
\description{
   Eastern Tropical Pacific core ODP-677
PI: Nicholas J. Shackleton
Core/Site: ODP-677  
Latitude: 1.12.00N ( 1.2)
Longitude: 83.44.00W ( -83.73333)
Water depth. m: 3450
   
	}		

\usage{data(ODP677)}

\format{
  A data frame on the following variables.
  \tabular{cl}{
 Variables:			\tab \cr
             depth        \tab  Depth (cm)			\cr
             yrBP         \tab radiocarbon years B.P. (Libby half-life; B.P. = 1950)			\cr
             d18Oforams-b \tab delta O18 PDB (Benthic Foraminifers)			\cr
             d18Oforams-p \tab delta O18 PDB (Planctic Foraminifers)		\cr
             notes        \tab Species		\cr
 }
 }

\details{
The time scale is astronmically-tuned, and motivated a revision of the date of the Brunhes-Matayama magnetic reversal.}
\note{The following attributes have been defined, corrsesponding to benthic species, uvigerinaes, and gruber. 
\code{ $benthics
 "CIB     " "CIBKULL " "CIBMELAT" "GLOBCASB" "GLOBOBUL" "GLOBOCAS" "GLOBOGYR" "GLOBORID"
 "GLOBORIG" "GYRGLOB " "GYROID  " "HELEGANS" "LATIC   " "ORIBULGL" "ORID    " "ORIDBUL "
 "ORIDBULL" "ORIDGLOB" "ORIDGYR " "ORIDMIX " "ORIDSBUL" "ORIDUVIG" "ORIGLOBU" "ORIGYRGL"
 "ORIGYRPY" "ORIGYRST" "PWUELCIB" "PWUELCIL" "PWUELL  " "PWUELLAT" "PWUELMEL" "PYRGGLOB"
 "PYRGO   " "SBULGBUL" "SBULL   " "SBULLGBU" "SBULLGLO" "SBULPYRG" "STYLOS  "}

\code{$uvig
 "UVIBULOR" "UVIG    " "UVIGBUL " "UVIGBULL" "UVIGGLOB" "UVIGLOB " "UVIGLOBO" "UVIGLORI"
 "UVIGORBU" "UVIGORGL" "UVIGORID" "UVIGORIG" "UVIGPYRG" "UVIGSBUL" "UVIORBUL" "UVIORGLO"
 "UVIPORID"}

\code{$gruber
 "GRUBER C"}
}
\source{
  http://www.geography.wisc.edu/classes/geog527/odp-677-fwc.txt
  See the file \file{PKG/orig/ODP/ODP677_benthics.txt} for details.
  Original version at 
\url{http://www.geography.wisc.edu/classes/geog527/odp-677-fwc.txt}
}

\references{
  N. J. Shackleton, A. Berger, and W. R. Peltier. An alternative astronomical calibration of the lower Pleistocene timescale based on  ODP  site 677. Transactions of the Royal Society of Edinburgh-Earth Sciences, 81:251-261, 1990.
     }

\examples{
data(ODP677)

## select and plot uvigerinae species
ODP677_uvig <- subset(ODP677,notes \%in\% attr(ODP677,"uvig"))
with(ODP677_uvig,plot(yrBP,d18Oforams.b,type="b"))
}
\keyword{datasets}
