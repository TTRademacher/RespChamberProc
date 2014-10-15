plotDurationUncertainty <- function(
	### plot the increase of uncertainty with decreaseing measurement duration
	ds
	,colTime="TIMESTAMP"	##<< column name of time [s]
	,fRegress = c(exp=regressFluxExp, lin=regressFluxLinear, tanh=regressFluxTanh)	##<< list of functions to yield 
			##<< a single flux estimate, see details of \code{\link{calcClosedChamberFlux}}
	,...	##<< further arguments to \code{\link{calcClosedChamberFlux}}
	,nDur = 20		##<< number of durations to check
	,maxSdFlux = 1	##<< maxium allowed standard deviation of flux in [mumol / s]
){
	times <- ds[[colTime]]
	times0 <- as.numeric(times) - as.numeric(times[1])
	resFit0 <- calcClosedChamberFlux(ds, colTime=colTime, fRegress=fRegress, ...)
	#resFit0$stat
	durations <- seq( max(65,resFit0$stat["tLag"]), max(times0), length.out=nDur+1)
	duration <- durations[1]
	#plot( CO2_dry ~ times0, ds)
	resFitsO <- lapply( durations[-c(nDur+1) ], function(duration){
				dss <- subset(ds, times0 <= duration )
				times0s <- times0[times0 <= duration]
				resFit <- calcClosedChamberFlux(dss, tLagFixed=resFit0$stat["tLag"]
					,colTime=colTime, fRegress=fRegress[resFit0$stat["iFRegress"]],...) 
				#plot( CO2_dry ~ times0s, dss)
				#lines( fitted(resFit$model) ~ times0s[times0s >= resFit0$stat["tLag"]], col="red")
				c(resFit, duration = max(times0s) )
			})
	resFits <- c(resFitsO, list(c(resFit0,duration=max(times0)) ))
	#resFit <- resFitsO[[1]]
	durationsR <- sapply( resFits, function(resFit){
				resFit$duration
			}) 
	tmp2 <- cbind( duration=durationsR, t(sapply( resFits, function(resFit){resFit$stat})))
	iMinTime <- if( min(tmp2[,"sdFlux"], na.rm=TRUE) <= maxSdFlux ) min(which( tmp2[,"sdFlux"] <= maxSdFlux )) else nrow(tmp2)
	minDuration <- tmp2[iMinTime,]
	##details<< 
	## Produces a plot with standard deviation of the flux estimate versus the duration of the measurment.
	## The lines correspond to the given maxium acceptable standard deviation
	## and the duration that matches this criterion.
	plot( sdFlux ~ duration, tmp2, xlab="Duration of measurement (s)" , ylab="sd(fluxEstimate)")
	abline(h = maxSdFlux, col="grey", lty="dashed" )
	abline(v = minDuration["duration"], , col="grey", lty="dashed" )
	#plot( flux ~ duration, tmp2 )
	#
	##value<< result of \code{\link{calcClosedChamberFlux}} for the minimum duration, with addition components 
	c(resFits[[ iMinTime ]][1:2]
		, duration=as.numeric(minDuration[1])	##<< minimum duration in seconds, with sdFlux < maxSdFlux (or maximum duration if criterion not met)
		, statAll= list(tmp2)					##<< component stat of the fits for each duration
	)
}
attr(plotDurationUncertainty,"ex") <- function(){
	data(chamberLoggerEx2)
	ds <- subset(chamberLoggerEx2, iChunk==99)	# very strong (and therefore precise) uptake
	#plot( CO2_dry ~ TIMESTAMP, ds )
	resDur <- plotDurationUncertainty( ds, colTemp="AirTemp", volume = 0.6*0.6*0.6, maxSdFlux = 0.8, nDur=10 )
	resDur$duration
	#plot( flux ~ duration, resDur$statAll )
	#plot( sdFlux ~ duration, resDur$statAll )
}
