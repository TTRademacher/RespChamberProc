
calcOpenChamberFlux <- function(
	### Calculate CO2 flux and its uncertainties for an open top chamber.
	ds						##<< data.frame with concentration and time column of a chamber measurement of one replicate
	,colConc="CO2_dry"		##<< column name of CO2 concentration [Amount of substance]
	,colTime="TIMESTAMP"	##<< column name of time
	,fRegress = regressFluxSquare	##<< function to yield a single flux estimate   
){
	##details<< 
	## details
  
  dslRes <- selectDataAfterLag(ds, colConc=colConc, colTime=colTime)
  times0 <- ds[,colTime]
  tLag <- times0[ dslRes$lagIndex ]
  dsl <- dslRes$ds
  times <- dsl[,colTime]
  times0 <- times - times[1]
  conc <- dsl[,colConc]
  
  times30 <- head(times,30)	#first 30 seconds for linear estimate
  conc30 <- head(conc,30)
  #plot( conc ~ times0 )
  resLinFlux <- regressFluxLinear( conc30, times30 )
  linFlux <- resLinFlux[1]
  sdResid <- sd(resid(attr(resLinFlux,"lm1")))
  #abline( coefficients(attr(resLinFlux,"lm1")) )
  if( abs(linFlux*diff(as.numeric(times30[c(1,length(times30))]))) < 1.5*sdResid ){
    warning("Flux magnitude smaller than noise, using linear flux estimte of the first 30 seconds")
    fRegress <- regressFluxLinear
	times <- times30
	conc <- conc30
  }
    
  	fluxEst <- fRegress( conc ,  times)
	#lines( fitted(attr(fluxEst,"lm1")) ~  times0 , col="maroon")
	#abline( coefficients(attr(fluxEst,"lm1"))[1], fluxEst[1], col="blue" )
  leverageEst <- sigmaBootLeverage(conc, times, fRegress=fRegress)
  ##details<<
  ## There are two kinds of uncertainty associated with the flux.
  ## The first comes from the uncertainty of the slope of concentration increase.
  ## The second comes from the leverage of starting and end points of the regression (estimated by a bootstrap)
  ## return value sdFlux is the maximum of those two components
  
	##value<< numceric vector with entries
	c(
		flux = fluxEst[1]			##<< the estimate of the CO2 flux [Amount per time]
		,sdFlux = max(fluxEst[2],leverageEst)	##<< the standard deviation of the CO2 flux
		,tLag = tLag		##<< lag-time between CO2 concentration in chamber and measurement at the sensor[XX]
		,sdFluxRegression = as.numeric(fluxEst[2]) ##<< the standard deviation of the flux by a single regression of CO2 flux
		,sdFluxLeverage = leverageEst	##<< the standard deviation of the flux by leverage of starting or end values of the time series
	)
	
}
attr(calcOpenChamberFlux,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
	conc <- ds$CO2_dry <- corrConcDilution(ds)  
  calcOpenChamberFlux(ds)
	
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
	attr(res,"lm1") <- lm1
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
  attr(res,"lm1") <- lm1
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

