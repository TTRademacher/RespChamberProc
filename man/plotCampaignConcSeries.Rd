\name{plotCampaignConcSeries}
\alias{plotCampaignConcSeries}
\title{plotCampaignConcSeries}
\description{get a series of ggplots of the time series and its fits}
\usage{plotCampaignConcSeries(ds, resL = NULL, varName = "CO2_dry", 
    idCol = "iChunk", timeCol = "TIMESTAMP", plotsPerPage = 64, 
    fileName = "", colIds = c(), ggplotList = c())}
\arguments{
  \item{ds}{data frame to plot, with collumns \code{idCol}, \code{timeCol} and \code{varName}}
  \item{resL}{list with results of \code{\link{calcClosedChamberFlux}} for each id-subset in ds}
  \item{varName}{variable to plot}
  \item{idCol}{collumn name of identifier of one time series}
  \item{timeCol}{collumn name of the time collumn}
  \item{plotsPerPage}{number of plots per page}
  \item{fileName}{if non-zero length string, the fileName where all plots are printed to  #paste0(varName,".pdf")}
  \item{colIds}{
}
  \item{ggplotList}{list added to each ggplot.}
}



\author{Thomas Wutzler, Oscar Perez Priego}





