\name{convert_mumolPers_to_gPerday}
\alias{convert_mumolPers_to_gPerday}
\title{convert mumolPers to gPerday}
\description{convert vom mumol CO2 / second to g C / day}
\usage{convert_mumolPers_to_gPerday(flux)}\arguments{
  \item{flux}{numeric vector or array: flux}
}
\details{Concentration measures are usually given by micromol CO2 across several seconds, and
the flux, i.e. its slope hence given in micromol CO2 per second.
To compare carbon balances, the units of gC per day are more convenient.  

mumol are converted to mol by /1e6
mol are converted to gC by *12
per second are converted to per day by *3600*24}
\value{numeric vector of flux in other units}

\author{Thomas Wutzler, Oscar Perez Priego}





