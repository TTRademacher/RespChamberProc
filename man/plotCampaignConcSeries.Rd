\name{plotCampaignConcSeries}
\alias{plotCampaignConcSeries}
\title{plotCampaignConcSeries}
\description{get a series of ggplots of the time series and its fits}
\usage{plotCampaignConcSeries(ds, resL = NULL, varName = "CO2_dry", 
    idCol = "iChunk", timeCol = "TIMESTAMP", qualityFlag = 0, 
    plotsPerPage = 64, fileName = "", colIds = c(), 
    ggplotList = c(), isVerbose = TRUE)}
\arguments{
  \item{ds}{data frame to plot, with collumns \code{idCol}, \code{timeCol} and \code{varName}}
  \item{resL}{list with results of \code{\link{calcClosedChamberFlux}} for each id-subset in ds}
  \item{varName}{variable to plot}
  \item{idCol}{collumn name of identifier of one time series}
  \item{timeCol}{collumn name of the time collumn}
  \item{qualityFlag}{vector of length nrow(ds) of a quality flag. For chunks where
this flag is not 0, subplots are dimmed.}
  \item{plotsPerPage}{number of plots per page}
  \item{fileName}{if non-zero length string, the fileName where all plots are printed to  #paste0(varName,".pdf")}
  \item{colIds}{
}
  \item{ggplotList}{list added to each ggplot.}
  \item{isVerbose}{
}
}



\author{Thomas Wutzler, Oscar Perez Priego}





