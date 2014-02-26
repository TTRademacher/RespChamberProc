
calcClosedChamberFlux <- function(
	### Calculate CO2 flux and its uncertainties for a non steady-state canopy chamber.
	ds						##<< data.frame with concentration and time column of a chamber measurement of one replicate
	,colConc="CO2_dry"		##<< column name of CO2 concentration [ppm]
	,colTime="TIMESTAMP"	##<< column name of time [s]
	,colTemp="TA_Avg"     ##<< column name of air temperature inside chamber [°C]
  ,colPressure="Pa"     ##<< column name of air pressure inside chamber [Pa]
	,fRegress = regressFluxTanh	##<< function to yield a single flux estimate, see details  
  ,volume=1            ##<< volume inside the chamber im [m3]
){
	##details<< 
	## The function \code{fRegress} must conform to \code{\link{regressFluxTanh}}, i.e.
  ## return a vector of length 2: the flux estimate and its standard deviation.
  ## Optionally, it may return the model fit object in attribute "model"

  dslRes <- selectDataAfterLag(ds, colConc=colConc, colTime=colTime)
  dsl <- dslRes$ds
  timesOrig <- ds[,colTime]
  times <- dsl[,colTime]
  times0 <- as.numeric(times) - as.numeric(times[1])
  tLag <- as.numeric(timesOrig[ dslRes$lagIndex ]) - as.numeric(timesOrig[1])
  conc <- dsl[,colConc]
  
  times30 <- head(times,30)	#first 30 seconds for linear estimate
  conc30 <- head(conc,30)
  #plot( conc ~ times0 )
  resLinFlux <- regressFluxLinear( conc30, times30 )
  linFlux <- resLinFlux[1]
  sdResid <- sd(resid(attr(resLinFlux,"model")))
  #abline( coefficients(attr(resLinFlux,"model")) )
  if( abs(linFlux*diff(as.numeric(times30[c(1,length(times30))]))) < 1.5*sdResid ){
    warning("Flux magnitude smaller than noise, using linear flux estimte of the first 30 seconds")
    fRegress <- regressFluxLinear
	times <- times30
	conc <- conc30
  }
    
  	fluxEst <- fRegress( conc ,  times)
	#lines( fitted(attr(fluxEst,"model")) ~  times0 , col="maroon")
	#abline( coefficients(attr(fluxEst,"model"))[1], fluxEst[1], col="blue" )
  leverageEst <- sigmaBootLeverage(conc, times, fRegress=fRegress)
  ##details<<
  ## There are two kinds of uncertainty associated with the flux.
  ## The first comes from the uncertainty of the slope of concentration increase.
  ## The second comes from the leverage of starting and end points of the regression (estimated by a bootstrap)
  ## return value sdFlux is the maximum of those two components
  
  #correct fluxes for density and express per chamber instead of per mol air
  fluxEstTotal = corrFluxDensity(fluxEst, volume = volume
                  , temp = ds[ dslRes$lagIndex,colTemp ]
                  , pressure = ds[ dslRes$lagIndex,colPressure ])
  leverageEstTotal = corrFluxDensity(leverageEst, volume = volume
                                     , temp = ds[ dslRes$lagIndex,colTemp ]
                                     , pressure = ds[ dslRes$lagIndex,colPressure ])
  
  #corrFluxDensity( dsl, vol=2)
	##value<< numceric vector with entries
	res <- c(
		flux = as.numeric(fluxEstTotal[1])			        ##<< the estimate of the CO2 flux [mumol / s]
		,sdFlux = max(fluxEstTotal[2],leverageEstTotal)	##<< the standard deviation of the CO2 flux
		,tLag = tLag		                                ##<< time of lag phase in seconds
		,sdFluxRegression = as.numeric(fluxEstTotal[2]) ##<< the standard deviation of the flux by a single regression of CO2 flux
		,sdFluxLeverage = leverageEstTotal	            ##<< the standard deviation of the flux by leverage of starting or end values of the time series
	)
	attr(res,"model") <- attr(fluxEst, "model")
	res
}
attr(calcClosedChamberFlux,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
  ds$Pa <- chamberLoggerEx1s$Pa * 1000  # convert kPa to Pa
	conc <- ds$CO2_dry <- corrConcDilution(ds)  
    res1 <- calcClosedChamberFlux(ds, fRegress=regressFluxLinear)
	res2 <- calcClosedChamberFlux(ds, fRegress=regressFluxSquare)
	res3 <- calcClosedChamberFlux(ds, fRegress=regressFluxExp )
	res4 <- calcClosedChamberFlux(ds, fRegress=regressFluxTanh )
	
	times <- ds$TIMESTAMP
	times0 <- as.numeric(times) - as.numeric(times[1])
	times0Fit <- times0[times0>res1["tLag"] ]
	#length(times0Fit)
	plot( conc ~ times0)
	abline(v=res1["tLag"])
	lines( fitted(attributes(res1)$model) ~ times0Fit , col="grey" )
	lines( fitted(attributes(res2)$model) ~ times0Fit , col="blue" )
	lines( fitted(attributes(res3)$model) ~ times0Fit , col="purple" )
	lines( fitted(attributes(res4)$model) ~ times0Fit , col="red" )
  
  res <- rbind(res1, res2, res3, res4)
	
}



selectDataAfterLag <- function(
  ### Omit the data within lag-time and normalize times to start after lag
  ds   ##<< data.frame with time and concentration columns
  ,colConc="CO2_dry"		##<< column name of CO2 concentration per dry air [ppm]
  ,colTime="TIMESTAMP"  	##<< column name of time column [s]
  
){
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
		### Estimate the initial flux by polynomial regression
		conc	  ##<< numeric vector of CO2 concentrations []
		,times	##<< times of conc measurements	[seconds]
){
  timesSec <- as.numeric(times) - as.numeric(times[1])
  lm1 <- lm( conc ~ poly(timesSec,2,raw = TRUE) )
  res <- c(
    flux = as.vector(coefficients(lm1)[2])
    ,sdFlux = summary(lm1)$coefficients[2,2]
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
){
	timesSec <- as.numeric(times) - as.numeric(times[1])
	#plot( conc ~ timesSec )
	c0 <- quantile( tail(conc, max(10,length(conc)%/%5) ), probs=0.25)
	#abline(h=c0)
	cDiff <- conc-c0; cDiff[cDiff <= 0] <- NA
	lm1 <- suppressWarnings( lm( log(cDiff) ~ timesSec ) )
	#plot( log(conc-c0) ~ timesSec )
	a0 <- exp(coefficients(lm1)[1])
	b0 <- -coefficients(lm1)[2]
	r0 <- a0/-b0
	nlm1 <- try(nls( conc ~ r/-b * exp(-b*timesSec) + c
		,start = list(r = r, b =b, c=c0)
		,control=nls.control(tol = 1e-03, maxiter=200, minFactor=1/1024/4)
	), silent=TRUE)
	if( inherits(nlm1,"try-error") ){
		nlm1 <- nls( conc ~ r/-b * exp(-b*timesSec) + c0
						,start = list(r = r0, b =b0)
						,control=nls.control(tol = 1e-03, maxiter=200, minFactor=1/1024/4)
				)
	}
	#lines( I(a0*exp(-b0*timesSec) + c0) ~ timesSec, col="maroon"  ) 
	#lines( fitted(nlm1) ~ timesSec, col="blue"  ) 
	#plot(resid(nlm1) ~ timesSec )
	#qqnorm( resid(nlm1) ); abline(0,1)
	res <- c(
			flux = coefficients(nlm1)[1]
			,sdFlux = sqrt(vcov(nlm1)[1,1])
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


regressFluxMenten <- function(
		### Estimate the initial flux by fitting a Michaelis-Menten type saturating function
		conc	  ##<< numeric vector of CO2 concentrations []
		,times	##<< times of conc measurements	[seconds]
){
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
		### Estimate the initial flux by fitting a Michaelis-Menten type saturating function
		conc	  ##<< numeric vector of CO2 concentrations []
		,times	##<< times of conc measurements	[seconds]
){
	timesSec <- as.numeric(times) - as.numeric(times[1])
	#plot( conc ~ timesSec )
	cSat0 <- quantile( tail(conc, max(10,length(conc)%/%5) ), probs=0.4)
	#abline(h=c0)
	#plot( conc-c0 ~ timesSec )
	lm1 <- lm(head(conc,30)-cSat0 ~ head(timesSec,30) )
	# initial slope in Michaelis menten is f=-r/m, substitute r = -m*f
	# r is the range to cover (intercept from the linear model)
	s0 <- coefficients(lm1)[2]
	c00 <- coefficients(lm1)[1]
	#plot( (c00 - (conc - cSat0))/c00  ~ timesSec )
	#lines( tanh(timesSec*-f0) ~ timesSec)
	#nlm0 <- nls( conc ~  cSat + c0*(1-tanh(timesSec*(-s)))
	#		,start = list(s=s0, cSat = cSat0, c0 = c00)
	#)
	#cf0 <- coefficients(nlm0)
	# f = s*c0; subsitute s = f/c0
	nlm1 <- try(nls( conc ~  cSat + c0*(1-tanh(timesSec*(-f/c0)))
			,start = list(f=s0*c00, cSat = cSat0, c0 = c00)
	), silent=TRUE)
	if( inherits(nlm1,"try-error") ){
		# fix cSat
		nlm1 <- nls( conc ~  cSat0 + c0*(1-tanh(timesSec*(-f/c0)))
				,start = list(f=s0*c00, c0 = c00)
		)
	}
	
	#lines( I(cSat0 + c00*(1-tanh(timesSec*(-f0)))) ~ timesSec, col="maroon"  ) 
	#lines( fitted(nlm1) ~ timesSec, col="purple"  )
	#plot(resid(nlm1) ~ timesSec )
	#qqnorm( resid(nlm1) ); abline(0,1)
	res <- c(
			flux = coefficients(nlm1)[1]
			,sdFlux = sqrt(vcov(nlm1)[1,1])
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
  timesSec <- as.numeric(times) - as.numeric(times[1])
  lm1 <- lm( conc ~ timesSec )
  res <- c(
    flux = as.vector(coefficients(lm1)[2])
    ,sdFlux = summary(lm1)$coefficients[2,2]
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
){
	##description<< 
	## 
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
    fRegress( conc[subIndices], times[subIndices] )[1]
  })
  sd(zz)
}
attr(sigmaBootLeverage,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
  sigmaBootLeverage( ds$CO2_Avg, ds$TIMESTAMP )
}

