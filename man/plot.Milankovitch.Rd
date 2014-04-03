\name{plot.Milankovitch}
\alias{plot.Milankovitch}
 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot Milankovitch graph}

\description{plot Milankovitch object}
\usage{
\method{plot}{Milankovitch} (x, months=TRUE, plot_function=image,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Milankovitch object }
  \item{months}{if true : x-axis of the plot indicates months
  conventionnally defined with the true solar longitude; x-axis is simply the true solar longitude
  otherwise}
  \item{...}{Other arguments passed to plotting function}
 \item{plot_function}{function used to plot the matrix. Typically
  \code{contour} or \code{image} but may also be \code{image.plot} if using the \code{fields} package.}
}

  
\author{Michel Crucifix, U. catholique de Louvain, Belgium.}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
