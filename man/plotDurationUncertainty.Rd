\name{plotDurationUncertainty}
\alias{plotDurationUncertainty}
\title{plotDurationUncertainty}
\description{plot the increase of uncertainty with decreaseing measurement duration}
\usage{plotDurationUncertainty(ds, fRegress = c(exp = regressFluxExp, 
    lin = regressFluxLinear, tanh = regressFluxTanh), ..., nDur = 20, 
    maxSdFlux = 1)}
\arguments{
  \item{ds}{
}
  \item{fRegress}{list of functions to yield
<< a single flux estimate, see details of \code{\link{calcClosedChamberFlux}}}
  \item{\dots}{further arguments to \code{\link{calcClosedChamberFlux}}}
  \item{nDur}{number of durations to check}
  \item{maxSdFlux}{maxium allowed standard deviation of flux in [mumol / s]}
}
\details{Produces a plot with standard deviation of the flux estimate versus the duration of the measurment.
The lines correspond to the given maxium acceptable standard deviation
and the duration that matches this criterion.}
\value{result of \code{\link{calcClosedChamberFlux}} for the minimum duration, with addition components 
\item{duration}{minimum duration in seconds, with sdFlux < maxSdFlux (or maximum duration if criterion not met)}
\item{statAll}{component stat of the fits for each duration
}}

\author{Thomas Wutzler, Oscar Perez Priego}




\examples{
data(chamberLoggerEx2)
ds <- subset(chamberLoggerEx2, iChunk==99)	# very strong (and therefore precise) uptake
#plot( CO2_dry ~ TIMESTAMP, ds )
resDur <- plotDurationUncertainty( ds, colTemp="AirTemp", volume = 0.6*0.6*0.6, maxSdFlux = 0.8, nDur=10 )
resDur$duration
#plot( flux ~ duration, resDur$statAll )
#plot( sdFlux ~ duration, resDur$statAll )
}
