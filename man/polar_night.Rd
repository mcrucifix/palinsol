% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Insol.R
\name{polar_night}
\alias{polar_night}
\title{Begin and end of the polar night}
\usage{
polar_night(lat, orbit)
}
\arguments{
\item{lat}{latitude}

\item{orbit}{Output from a solution, such as \code{ber78}, \code{ber90} or
\code{la04}}
}
\value{
Either a message about the absence of polar night (for specified
reasons), or the true solar longitude, in degrees, of the beginning and end
of the polar night.
}
\description{
Provides the true solar longitude (in degrees) of the beginning and end of
the polar night for a given latitude and orbit
}
\examples{


current_orbit <- la04(0)

# polar night at the equator ? 
polar_night (0, current_orbit)

# polar night at 80 N ? 
polar_night (80*pi/180, current_orbit)

# polar nights expressed as day of year 
l2day(current_orbit, polar_night (80*pi/180, current_orbit))


}
\references{
any standard text book of spherical astronomy
}
\author{
Michel Crucifix, U. catholique de Louvain, Belgium.
}
