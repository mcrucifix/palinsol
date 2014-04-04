\name{BER78}
\alias{BER78}
\alias{BER90}

\docType{data}
\title{Tables supplied by BER78 and BER90}
\description{
   Tables necessary for computation of astronomical solutions by
   Berger (1978) and
   Beger and Loutre (1990)
   }
\usage{
 data(BER78) 
 data(BER90) 
 }
 \format{\code{Lists},  named \code{BER78} and \code{BER90}, respectively, containing as elements the 
         original published Tables as \code{data.frames}}
% \format{text files reproducing the published tables with original units.}

\references{
 A. L.  Berger (1978).  Long-term variations of daily
insolation and Quaternary climatic changes, J. Atmos. Sci., 35,
2362-2367.
 
 A. Berger and M. F. Loutre (1990), Origine des frequences 
  des elements astronomiques intervenant dans l'insolation,  
  Bull. Classe des Sciences, 1-3, 45-106 

 A. Berger and M.F. Loutre (1991), Insolation values for the climate of 
 the last 10 million years,  Quaternary Science Reviews, 10, 297 - 317. }


\note{Table2 of BER78 was reconstructed based on details given by Berger and Loutre 1990, the code
developed to this end is available in the \code{orig} directory.}

\author{Michel Crucifix, U. catholique de Louvain, Belgium.}

\source{ftp://ftp.elic.ucl.ac.be/berger/berger78/ 
        ftp://ftp.elic.ucl.ac.be/loutre/QSR/ }

\examples{ 
data(BER78)
data(BER90)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
