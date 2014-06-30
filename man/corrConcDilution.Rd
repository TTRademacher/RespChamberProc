\name{corrConcDilution}
\alias{corrConcDilution}
\title{corrConcDilution}
\description{Calculate concentration corrected for dilution with water vapor		}
\usage{corrConcDilution(ds, colConc = "CO2_Avg", colVapour = "H20_Avg")}
\arguments{
  \item{ds}{data frame with each row one observations, and respective columns}
  \item{colConc}{column name of CO2 concentration [ppm]}
  \item{colVapour}{column name of CO2 concentration [ppt]}
}
\details{If CO2 concentration is measured per moist air, this function will calculate the concentration\
per dry air.}
\value{numeric vector (nrow ds):  concentration of CO2 per dry air [ppm]}
\references{LI-COR, Application Note 129. The Importance of Water Vapor Measurements and Corrections. LI-COR, Inc., 4421 Superior Street, Lincoln, NE 68504, USA.}
\author{Thomas Wutzler, Oscar Perez Priego}




\examples{
data(chamberLoggerEx1s)
ds <- chamberLoggerEx1s
ds$CO2_dry <- corrConcDilution(ds)	
}
