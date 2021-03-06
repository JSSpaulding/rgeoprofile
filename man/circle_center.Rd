% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circle_center.R
\name{circle_center}
\alias{circle_center}
\title{Center of the Circle Calculation for Geographic Profiling}
\usage{
circle_center(lat, lon)
}
\arguments{
\item{lat}{a vector of latitudes for the crime incident series}

\item{lon}{a vector of latitudes for the crime incident series}
}
\value{
A latitude and longitude coordinate for the center of the circle
    which encompasses the incidents. This point can be used to prioritize
    the area which contains the offender's anchor point.
}
\description{
A calculation for the center of the circle which bounds the
    incident coordinates for serial crime analysis. This function is among
    the centrographic methods which have been used for geographic profiling.
    The model assumes that the serial perpetrator's home base is relatively
    central among the crime incidents.
}
\examples{
#Using provided dataset for the Boston Strangler Incidents:
data(desalvo)
circle_center(desalvo$lat, desalvo$lon)
}
\author{
Jamie Spaulding, Keith Morris
}
\keyword{methods}
\keyword{spatial}
