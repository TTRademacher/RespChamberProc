\name{corrFluxDensity}
\alias{corrFluxDensity}
\title{corrFluxDensity}
\description{Calculate total flux inside the chamber from flux per amount of air}
\usage{corrFluxDensity(CO2_molarFlux, volume = 1, temp = 20, pressure = 101325)}
\arguments{
  \item{CO2_molarFlux}{numeric vector of rate of changes in CO2 dry molar fraction [ppm/sec]}
  \item{volume}{numeric scalar: volume of the chamber in [m3]}
  \item{temp}{numeric vector: temperature inside chamber [degC]}
  \item{pressure}{numeric vector: pressure inside chamber [Pa]}
}
\details{The amount of air is determined by the volumen, pressure and temperature (universal gas law)}
\value{numeric vector (nrow ds):  CO2 flux [mumol CO2 /s]}

\author{Thomas Wutzler, Oscar Perez Priego}




\examples{
data(chamberLoggerEx1s)
ds <- chamberLoggerEx1s
  CO2_molarFlux <- ds$CO2_Avg   # just to have a column (here concentrations), the flux should be calculated properly before
ds$CO2_flux <- corrFluxDensity(CO2_molarFlux, pressure=101*1000) #kPa converted to Pa	
  plot(ds$CO2_flux)
}
