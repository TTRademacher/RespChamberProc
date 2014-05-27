
calcClosedChamberFlux <- function(
	### Calculate CO2 flux and its uncertainties for a non steady-state canopy chamber.
	ds						##<< data.frame with concentration and time column of a chamber measurement of one replicate
	,colConc="CO2_dry"		##<< column name of CO2 concentration [ppm]
	,colTime="TIMESTAMP"	##<< column name of time [s]
	,colTemp="TA_Avg"       ##<< column name of air temperature inside chamber [°C]
    ,colPressure="Pa"       ##<< column name of air pressure inside chamber [Pa]
	,fRegress = c(lin=regressFluxLinear, tanh=regressFluxTanh, exp=regressFluxExp)	##<< list of functions to yield a single flux estimate, see details  
  	,volume=1               ##<< volume inside the chamber im [m3]
	,isEstimateLeverage	= TRUE	##<< set to FALSE to omit the time consuming bootstrap for uncertainty due to leverage
    ,isStopOnError = TRUE   ##<< set to FALSE to not stop execution on errors, but report NAs
){
	##seealso<< \code{\link{RespChamberProc}}
	
  ##details<< 
  ## The function \code{fRegress} must conform to \code{\link{regressFluxSquare}}, i.e.
  ## return a vector of length 2: the flux estimate and its standard deviation.
  ## Optionally, it may return the model fit object in attribute "model"
  ## If several functions are given, then the best fit is selected according to AIC criterion.
  #
  #plot( ds[,colConc] ~ ds[,colTime] )
  if( !length(names(fRegress)) ) names(fRegress) <- 1:length(fRegress)
  dslRes <- selectDataAfterLag(ds, colConc=colConc, colTime=colTime)
  dsl <- dslRes$ds
  timesOrig <- ds[,colTime]
  times <- dsl[,colTime]
  times0 <- as.numeric(times) - as.numeric(times[1])
  tLag <- as.numeric(timesOrig[ dslRes$lagIndex ]) - as.numeric(timesOrig[1])
  #abline(v=tLag+ds[1,colTime] )
  conc <- dsl[,colConc]
  #
  # removed check for linear fit
  # plot( conc ~ times )
  res <- #try({
	#fReg <- fRegress[[1]]
	#trace(fReg, recover)
    fluxEstL <- lapply( fRegress, function(fReg){	
    	fluxEst <- fReg( conc ,  times)
    })
    iBest <- which.min( do.call(rbind,fluxEstL)[,3] )
    fReg <- fRegress[[iBest]]
    fluxEst <- fluxEstL[[iBest]]
  	#lines( fitted(attr(fluxEst,"model")) ~  times0 , col="maroon")
  	#abline( coefficients(attr(fluxEst,"model"))[1], fluxEst[1], col="blue" )
	coefEst <- coefficients(attr(fluxEst,"model"))
    leverageEst <- if( isTRUE(isEstimateLeverage) ) sigmaBootLeverage(conc, times, fRegress=fReg
		, coefStart=coefEst ) else NA
    ##details<<
    ## There are two kinds of uncertainty associated with the flux.
    ## The first comes from the uncertainty of the slope of concentration increase.
    ## The second comes from the leverage of starting and end points of the regression (estimated by a bootstrap)
    ## return value sdFlux is the maximum of those two components
    #
    #correct fluxes for density and express per chamber instead of per mol air
    fluxEstTotal = corrFluxDensity(fluxEst, volume = volume
                    , temp = ds[ dslRes$lagIndex,colTemp ]
                    , pressure = ds[ dslRes$lagIndex,colPressure ])
    leverageEstTotal = corrFluxDensity(leverageEst, volume = volume
                                       , temp = ds[ dslRes$lagIndex,colTemp ]
                                       , pressure = ds[ dslRes$lagIndex,colPressure ])
    #
    #corrFluxDensity( dsl, vol=2)
  	##value<< list with entries
  	res <- c(
  		flux = as.numeric(fluxEstTotal[1])			        ##<< the estimate of the CO2 flux [mumol / s]
		,fluxMedian = as.numeric(leverageEstTotal[3])
  		,sdFlux = max(fluxEstTotal["sdFlux"],leverageEstTotal["sd"], na.rm=TRUE)	##<< the standard deviation of the CO2 flux
  		,tLag = tLag		                                ##<< time of lag phase in seconds
  		,lagIndex = dslRes$lagIndex 					##<< index of the row at the end of lag-time
  		,sdFluxRegression = as.numeric(fluxEstTotal["sdFlux"]) ##<< the standard deviation of the flux by a single regression of CO2 flux
  		,sdFluxLeverage = as.numeric(leverageEstTotal["sd"])	    ##<< the standard deviation of the flux by leverage of starting or end values of the time series
		,iFRegress=as.numeric(iBest)					##<< index of the best (lowest AIC) regression function
  	)
    attr(res,"model") <- attr(fluxEst, "model")
    res
  #}, silent=TRUE)    
  if( inherits(res,"try-error") ) res <- structure(rep(NA_real_,6),names=c("flux","sdFlux","tLag","lagIndex","sdFluxRegression","sdFluxLeverage"))
  res
}
attr(calcClosedChamberFlux,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
    ds$Pa <- chamberLoggerEx1s$Pa * 1000  # convert kPa to Pa
	
	conc <- ds$CO2_dry <- corrConcDilution(ds)
    resLin <- calcClosedChamberFlux(ds, fRegress=list(regressFluxLinear))
	resSquare <- calcClosedChamberFlux(ds, fRegress=list(regressFluxSquare))
	#resExp <- calcClosedChamberFlux(ds, fRegress=regressFluxExp )
	resTanh <- calcClosedChamberFlux(ds, fRegress=list(regressFluxTanh ))
	
	times <- ds$TIMESTAMP
	times0 <- as.numeric(times) - as.numeric(times[1])
	times0Fit <- times0[times0>resLin["tLag"] ]
	#length(times0Fit)
	plot( ds$CO2_dry ~ times0, xlab="time (s)", ylab="" ); mtext("CO2_dry (ppm)",2,2,las=0)
	abline(v=resLin["tLag"], col="grey", lty="dotted")
	lines( fitted(attributes(resLin)$model) ~ times0Fit , col="grey" )
	lines( fitted(attributes(resSquare)$model) ~ times0Fit , col="blue" )
	#lines( fitted(attributes(resExp)$model) ~ times0Fit , col="purple" )
	if( resLin[1] > 0 ){
		lines( fitted(attributes(resTanh)$model) ~ times0Fit , col="red" )
		legend( "bottomright", inset=c(0.01,0.01), legend=c("Linear","Polynomial","Tanh"), col=c("grey","blue","red"), lty=1)
	}else{
		lines( -fitted(attributes(resTanh)$model) ~ times0Fit , col="red" )
		legend( "topright", inset=c(0.01,0.01), legend=c("Linear","Polynomial","Tanh"), col=c("grey","blue","purple","red"), lty=1)
	}
	
     res <- rbind(resLin, resSquare, resTanh)
	 res
	
}



selectDataAfterLag <- function(
  ### Omit the data within lag-time and normalize times to start after lag
  ds   ##<< data.frame with time and concentration columns
  ,colConc="CO2_dry"		##<< column name of CO2 concentration per dry air [ppm]
  ,colTime="TIMESTAMP"  	##<< column name of time column [s]
  
){
	##seealso<< \code{\link{RespChamberProc}}
	
  # TODO: account for different times
  iBreak <- min(cpt.mean(ds[,colConc],penalty="SIC",method="PELT",class=FALSE)) 
  ##value<< A list with entries
  list(
    lagIndex = iBreak    ##<< the index of the end of the lag period
    ,ds = ds[ (iBreak+1):nrow(ds), ]   ##<< the dataset ds without the lag-period
  )
}
attr(selectDataAfterLag,"ex") <- function(){
  data(chamberLoggerEx1s)
  ds <- chamberLoggerEx1s
  ds$CO2_dry <- corrConcDilution(ds)
  ret <- selectDataAfterLag(ds)
  plot( ds$CO2_dry)
  abline(v=ret$lagIndex, col="red")
}

regressFluxSquare <- function(
		### Estimate the initial flux and its standard deviation by polynomial regression
		conc	  ##<< numeric vector of CO2 concentrations [ppm]
		,times	##<< times of conc measurements	[seconds]
){
  ##seealso<< \code{\link{RespChamberProc}}
  
  ##details<< 
  ## The flux is calculated at the slop of the concnetration change. By
  ## changing the concentration gradient, however, the flux is disturbed. In effect the 
  ## flux will decline over time and concentrations go towards a saturation.
  ##
  ## This method fits a polynomial regression to the concentrations and infers the slope at reports
  ## the slope at the initial time.
  ##
  ## Other functional forms can be fitted to estimate the initial slope:
  ## \itemize{
  ##  \item{ Linear: \code{\link{regressFluxLinear}}  }
  ##  \item{ Hyperbolic tangent saturating function: \code{\link{regressFluxTanh}}  }
  ##  \item{ Exponential: \code{\link{regressFluxExp}}  }
  ##  \item{ Michaelis-Menten type saturating function: \code{\link{regressFluxMenten}}  }
  ## }
  ##
  ## The hyperbolic tangent form (\code{\link{regressFluxTanh}}) has the advantage that 
  ## initially the flux is changing only very slowly. In contrast, whith the exponential form the
  ## slope changes much at the beginning.
		  
  timesSec <- as.numeric(times) - as.numeric(times[1])
  lm1 <- lm( conc ~ poly(timesSec,2,raw = TRUE) )
  res <- c(
    flux = as.vector(coefficients(lm1)[2])
    ,sdFlux = summary(lm1)$coefficients[2,2]
	, AIC = AIC(lm1)
    )
	attr(res,"model") <- lm1
	res
	### numeric vector (2): estimate and its standard deviation of the initial flux [ppm / s]
}
attr(regressFluxSquare,"ex") <- function(){
  data(chamberLoggerEx1s)
  ds <- chamberLoggerEx1s
  conc <- ds$CO2_dry <- corrConcDilution(ds)  
  times <- ds$TIMESTAMP
  regressFluxSquare( conc, times  )
  plot( conc ~ times)
  #lines( fitted(lm1) ~ times )
}

regressFluxExp <- function(
		### Estimate the initial flux by fitting an exponentially saturating function
		conc	  ##<< numeric vector of CO2 concentrations []
		,times	##<< times of conc measurements	[seconds]
		,cSatFac=2  ##<< Position of the initial saturation (0 start, 1 end, >1 outside measured range) 
){
	##seealso<< \code{\link{regressFluxSquare}}
	##seealso<< \code{\link{RespChamberProc}}
	
	timesSec <- as.numeric(times) - as.numeric(times[1])
	#plot( conc ~ timesSec )
	fluxLin <- coefficients(lm(conc ~ timesSec ))[2]
	concP <- if( fluxLin < 0 ) conc else -conc		# invert to decreasing concentrations
	#plot( concP ~ timesSec )
	# increasing concentration
	# set the saturation twice above the range
	cRange <- quantile( concP , probs=c(0.05,0.95))
	cSat0 <- cRange[1] - cSatFac*diff(cRange)
	#abline(h=cSat0)
	#plot( I(concP-cSat0) ~ timesSec )
	cDiff <- concP-cSat0; cDiff[cDiff <= 0] <- NA
	lm1 <- suppressWarnings( lm( log(cDiff) ~ timesSec ) )
	#plot( log(conc-cSat0) ~ timesSec )
	a0 <- exp(coefficients(lm1)[1])
	b0 <- -coefficients(lm1)[2]
	r0 <- a0/-b0
	nlm1 <- if( fluxLin < 0 ){
		nlm1 <- try(nls( conc ~ r/-b * exp(-b*timesSec) + c
			,start = list(r = r0, b =b0, c=cSat0)
			,control=nls.control(tol = 1e-03, minFactor=1/1024/4)
		), silent=TRUE)
	} else {
		nlm1 <- try(nls( conc ~ r/b * exp(b*timesSec) - c		# concB -> -conc; b->-b; r->-r
						,start = list(r = -r0, b =-b0, c=cSat0)
						,control=nls.control(tol = 1e-03, minFactor=1/1024/4)
				), silent=TRUE)
	}
	if( inherits(nlm1,"try-error") ){
		nlm1 <- if( fluxLin < 0 ){
				nlm1 <- nls( conc ~ r/-b * exp(-b*timesSec) + cSat0
						,start = list(r = r0, b =b0)
						,control=nls.control(tol = 1e-03, maxiter=200, minFactor=1/1024/4)
				)
			} else {
				nlm1 <- nls( conc ~ r/b * exp(b*timesSec) - cSat0
						,start = list(r = -r0, b =-b0)
						,control=nls.control(tol = 1e-03, maxiter=200, minFactor=1/1024/4)
				)
			}
	}
	#plot( conc ~ timesSec )
	#lines( I(a0*exp(-b0*timesSec) + cSat0) ~ timesSec, col="maroon"  ) 
	#lines( fitted(nlm1) ~ timesSec, col="blue"  ) 
	#lines( fitted(nlm1) ~ timesSec, col="orange"  ) 
	#plot(resid(nlm1) ~ timesSec )
	#qqnorm( resid(nlm1) ); abline(0,1)
	res <- c(
			flux = coefficients(nlm1)[1]
			,sdFlux = sqrt(vcov(nlm1)[1,1])
			, AIC = AIC(nlm1)
	)
	attr(res,"model") <- nlm1
	res
	### numeric vector (2): estimate and its standard deviation of the initial flux [ppm / s]
}
attr(regressFluxExp,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
	conc <- ds$CO2_dry <- corrConcDilution(ds)  
	times <- ds$TIMESTAMP
	regressFluxExp( conc, times  )
	plot( conc ~ times)
	#lines( fitted(lm1) ~ times )
}


.regressFluxMenten <- function(
		### Estimate the initial flux by fitting a Michaelis-Menten type saturating function
		conc	  ##<< numeric vector of CO2 concentrations []
		,times	##<< times of conc measurements	[seconds]
){
	##seealso<< \code{\link{regressFluxSquare}}
	##seealso<< \code{\link{RespChamberProc}}
	
	timesSec <- as.numeric(times) - as.numeric(times[1])
	#plot( conc ~ timesSec )
	c0 <- quantile( tail(conc, max(10,length(conc)%/%5) ), probs=0.4)
	#abline(h=c0)
	plot( conc-c0 ~ timesSec )
	lm1 <- lm(head(conc,30)-c0 ~ head(timesSec,30) )
	# initial slope in Michaelis menten is f=-r/m, substitute r = -m*f
	# r is the range to cover (intercept from the linear model)
	f0 <- coefficients(lm1)[2]
	m0 <- -coefficients(lm1)[1] / f0
	nlm1 <- nls( conc ~ c + -m*f * (1- timesSec/(m+timesSec)) 
			,start = list(f=f0, m = m0, c = c0)
	)
	#lines( I(c0 + -m0*f0*(1-timesSec/(m0+timesSec))) ~ timesSec, col="maroon"  ) 
	#lines( fitted(nlm1) ~ timesSec, col="blue"  ) 
	#plot(resid(nlm1) ~ timesSec )
	#qqnorm( resid(nlm1) ); abline(0,1)
	res <- c(
			flux = coefficients(nlm1)[1]
			,sdFlux = sqrt(vcov(nlm1)[1,1])
			, AIC = AIC(nlm1)
	)
	attr(res,"model") <- nlm1
	res
	### numeric vector (2): estimate and its standard deviation of the initial flux [ppm / s]
}
attr(regressFluxMenten,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
	conc <- ds$CO2_dry <- corrConcDilution(ds)  
	times <- ds$TIMESTAMP
	regressFluxMenten( conc, times  )
	plot( conc ~ times)
	#lines( fitted(lm1) ~ times )
}

regressFluxTanh <- function(
		### Estimate the initial flux by fitting a hyperbolic tangent saturating function
		conc	  ##<< numeric vector of CO2 concentrations []
		,times	  ##<< times of conc measurements	[seconds]
		,cSatFac=2  ##<< Position of the initial saturation (0 start, 1 end, >1 outside measured range)
		,corAR1=TRUE	##<< try to account for autocorrelation
		,start=c()	##<< numeric vector of starting parameters. May provide from last bootstrap to speed up fitting  
){
	##seealso<< \code{\link{regressFluxSquare}}
	##seealso<< \code{\link{RespChamberProc}}
	#
	timesSec <- as.numeric(times) - as.numeric(times[1])
	fluxLin <- coefficients(lm(conc ~ timesSec ))[2]
	if( !length(start) ){
		concP <- if( fluxLin < 0 ) -conc else conc		# invert, do not forget to invert again when flux calculated
		# increasing concentration
		# set the saturation twice above the range
		cRange <- quantile( concP , probs=c(0.05,0.95))
		cSat0 <- cRange[1] + cSatFac*diff(cRange)
		#abline(h=cSat0)
		#plot( concP-cSat0 ~ timesSec )
		lm1 <- lm(head(concP,30)-cSat0 ~ head(timesSec,30) )
		#abline(coefficients(lm1))
		# c0 is the range to cover (= -intercept from the linear model)
		s0 <- coefficients(lm1)[2]	     # initial slope
		c00 <- -coefficients(lm1)[1]     # range to cover
	}
	#plot( (concP - cSat0)/c00 +1  ~ timesSec )		# that matches tanh def between 0 and 1
	#lines( tanh(timesSec*s0/c00) ~ timesSec)	 	# tanh function  - set equal to above and solve for concP
	#plot( concP ~ timesSec )
	#lines( (tanh(timesSec*s0/c00)-1)*c00 + cSat0)  # initial in concP range
	#lines( (tanh(timesSec*coefficients(nlm1)[1]/coefficients(nlm1)[3])-1)*coefficients(nlm1)[3] + coefficients(nlm1)[2])  # initial in concP range
	#plot( conc ~ timesSec )
	#lines( (tanh(timesSec*-s0/c00)+1)*c00 -cSat0)  # initial in concP range
	#lines( (tanh(timesSec*-coefficients(nlm1)[1]/coefficients(nlm1)[3])+1)*coefficients(nlm1)[3] - coefficients(nlm1)[2])  # initial in concP range
	nlm1 <- try( 
			# deprecated use nlsLM from minpack.lm to switch between Newten and Levenberg, avoid singular gradient in near-linear cases
			# http://www.r-bloggers.com/a-better-nls/
			nlm1 <- if( fluxLin < 0 ){
					gc1 <- gnls( conc ~   (tanh(timesSec*s/c0)-1)*c0 + cSat 
						,start = if(length(start)) start else list(s=s0, cSat = cSat0, c0 = c00)
						,params= c(s+cSat+c0~1)
						,correlation=corAR1( 0.3, form = ~ timesSec)
					)
				} else {
					# increasing concentrations
					gc1 <-gnls( conc ~  (tanh(timesSec*s/c0)+1)*c0 - cSat
						,start = if(length(start)) start else list(s=-s0, cSat = cSat0, c0 = c00)
						,params= c(s+cSat+c0~1)
						,correlation=corAR1( 0.3, form = ~ timesSec)
					)
				}
	,silent=TRUE)
	# lines(fitted(nlm1) ~ timesSec, col="orange" )
	#
	#lines( I(cSat0 + c00*(1-tanh(timesSec*(-s0)))) ~ timesSec, col="maroon"  ) 
	#lines( fitted(nlm1) ~ timesSec, col="purple"  )
	#plot(resid(nlm1) ~ timesSec )
	#qqnorm( resid(nlm1) ); abline(0,1)
	res <- if( inherits(nlm1,"try-error") ) c(
				flux = NA
				,sdFlux = NA
				, AIC = NA
				,autoCorr =  NA 
	) else c(
				flux = coefficients(nlm1)[1]
				,sdFlux = sqrt(vcov(nlm1)[1,1])
				, AIC = AIC(nlm1)
				,autoCorr =  as.numeric(nlme:::coef.corAR1(nlm1$modelStruct$corStruct, unconstrained=FALSE)) 
		)		
	attr(res,"model") <- nlm1
	res
	### numeric vector (2): estimate and its standard deviation of the initial flux [ppm / s]
}
attr(regressFluxTanh,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s[-(1:16),]
	conc <- ds$CO2_dry <- corrConcDilution(ds)  
	times <- ds$TIMESTAMP
	times0 <- as.numeric(times) - as.numeric(times[1])
	regressFluxTanh( conc, times  )
	plot( conc ~ times0)
	#lines( fitted(lm1) ~ times )
	
	
}




regressFluxLinear <- function(
  ### Estimate the initial flux by polynomial regression
  conc	  ##<< numeric vector of CO2 concentrations []
  ,times	##<< times of conc measurements	[seconds]
){
	##seealso<< \code{\link{regressFluxSquare}}
	##seealso<< \code{\link{RespChamberProc}}
	
  timesSec <- as.numeric(times) - as.numeric(times[1])
  lm1 <- lm( conc ~ timesSec )
  res <- c(
    flux = as.vector(coefficients(lm1)[2])
    ,sdFlux = summary(lm1)$coefficients[2,2]
	,AIC=AIC(lm1)
  )
  attr(res,"model") <- lm1
  res
  ### numeric vector (2): estimate and its standard deviation of the initial flux [ppm / s]
}
attr(regressFluxLinear,"ex") <- function(){
  data(chamberLoggerEx1s)
  ds <- chamberLoggerEx1s
  conc <- ds$CO2_Avg
  times <- ds$TIMESTAMP
  plot( conc ~ times )
  regressFluxLinear( conc, times  )
}


sigmaBootLeverage <- function(
		### Estimate uncertainty of regrssion due to leverage of starting and end points. 
		conc	  ##<< numeric vector of CO2 concentrations 
		,times	##<< times of conc measurements
		,fRegress = regressFluxSquare	##<< function to yield a single flux estimate
		,probs=c(0.05,0.5,0.95)			##<< percentiles for which the flux quantiels shoudl be returned
		,nSample=80
		,coefStart=c()	##<< numeric vector of starting parameters. May provide from last bootstrap to speed up fitting  
){
	##seealso<< \code{\link{RespChamberProc}}
	
  periodLength <- diff( as.numeric(times[c(1,length(times))]) )
  if( periodLength < 60 ) warning(paste("Time series of only",periodLength," seconds is too short. Recommended are at least 60 seconds"))
  start <- seq (0, 10)   # indices of starting the time series
  close <- seq(max(15, length(conc)-40), length(conc), 1) # indices of the end (deployment) of the duration
  
  ##defining the function to be bootstrapped based on starting and deployment time
  nSample=80
  starts <-  sample(start, nSample, replace = TRUE)
  closes <-  sample(close, nSample, replace = TRUE)
  zz <- sapply(1:nSample, function(i){
    subIndices <- starts[i]:closes[i]
    fRegress( conc[subIndices], times[subIndices], start=coefStart )[1]
  })
  #plot(density(zz, na.rm=TRUE))
  ## a vector with first entry the standard deviation of the flux, and following entries quantiles for given argument \code{probs}
  c( sd=sd(zz, na.rm=TRUE), quantile(zz, probs, na.rm=TRUE) ) 
}
attr(sigmaBootLeverage,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
  sigmaBootLeverage( ds$CO2_Avg, ds$TIMESTAMP )
}

