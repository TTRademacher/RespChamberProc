\name{plotCampaignConcSeries}
\alias{plotCampaignConcSeries}
\title{plotCampaignConcSeries}
\description{plot all the time series into a series of pdf-pages}
\usage{plotCampaignConcSeries(ds, resL = NULL, varName = "CO2_dry", 
    idCol = "iChunk", timeCol = "TIMESTAMP", plotsPerPage = 64, 
    fileName = paste0(varName, ".pdf"), colIds = c())}\arguments{
  \item{ds}{data frame to plot, with collumns "id", "TIMESTAMP" and \code{varName}}
  \item{resL}{list with results of \code{\link{calcClosedChamberFlux}} for each id-subset in ds}
  \item{varName}{variable to plot}
  \item{idCol}{collumn name of identifier of one time series}
  \item{timeCol}{collumn name of the time collumn}
  \item{plotsPerPage}{number of plots per page}
  \item{fileName}{fileName of the output pdf}
  \item{colIds}{
}
}



\author{Thomas Wutzler, Oscar Perez Priego}




\examples{
if( FALSE ){		# depends on data not available with the package, but demonstrates the way
	load(file.path(projectDir,"SMANIEconc.Rd"))
	ds0 <- ds <- SMANIEconc
	.createresCalcClosedChamberFlux <- function(){
		# note how results of calcClosedChamberFlux are assembled in a list
		resCalcClosedChamberFlux <- resL <- dlply( ds, .(id), function(dss){
					cat( ",", levels(dss$id)[dss$id[1]] )
					res2 <- res <- try( 
							res <- calcClosedChamberFlux( dss, colTemp="AirTemp", volume = 0.6*0.6*0.6
									, debugInfo = list(useOscarsLagDectect = TRUE, omitAutoCorrFit=TRUE, omitEstimateLeverage = TRUE)
							)
							,silent=TRUE)
				})
		save( resCalcClosedChamberFlux, file=file.path(projectDir,"resCalcClosedChamberFlux_SMANIEconc.Rd"))
	}
	load( file.path(projectDir,"resCalcClosedChamberFlux_SMANIEconc.Rd"))
	resL <- resCalcClosedChamberFlux
	names(resL)[ sapply( resL, inherits,  "try-error") ]	# "1.2.19"  "1.3.2"   "1.3.23"  "2.2.270" "2.2.277"
	#
	# plot to file
	fileName=file.path(projectDir,paste0("SMANIEconc_CO2_dry.pdf"))	##<< fileName
	plotCampaignConcSeries(ds, resL=resL, fileName=fileName)
}
}
