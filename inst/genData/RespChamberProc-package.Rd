\name{RespChamberProc-package}
\alias{RespChamberProc-package}
\alias{RespChamberProc}
\title{Processing data from respiration chambers .}
\description{
Processing data from respiration chambers.
}
\details{

A first set of functions correct measured CO2 concentrations
\itemize{
\item{ Changing air density with pressure and temperature: \code{\link{corrConcDensity}}  }
\item{ Dilution with water vapor: \code{\link{corrConcDilution}}  }
\item{ Leakage of the chamber: \code{\link{corrConcLeakage}}  }
}%itemize


A second set of functions estimates the flux and its uncertainty from the concentrations 
\itemize{
\item{ Aggregating estimates: \code{\link{calcOpenChamberFlux}}  }
\item{ Select the subset of data after lagTime: \code{\link{selectDataAfterLag}}  }
\item{ Estimate the initial flux by polynomial regression: \code{\link{regressFlux}}  }
\item{  Estimate uncertainty of regrssion due to leverage of starting and end points: \code{\link{sigmaBootLeverage}}  }
}%itemize

Further misc functionality: \itemize{
}%itemize

}%details

\author{Thomas Wutzler, Oscar Pérez Priego}

\keyword{ package }


