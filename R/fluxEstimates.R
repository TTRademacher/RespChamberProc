
calcOpenChamberFlux <- function(
	### Calculate CO2 flux and its uncertainties for an open top chamber.
	ds						##<< data.frame with concentration and time column of a chamber measurement of one replicate
	,colConc="CO2_corr"		##<< column name of CO2 concentration [Amount of substance]
	,colTime="TIMESTAMP"	##<< column name of time
){
	##details<< 
	## details
	
	##value<< numceric vector with entries
	c(
		flux = flux			##<< the estimate of the CO2 flux [Amount per time]
		,sdFlux = sdFlux	##<< the standard deviation of the CO2 flux
		,tLag = tLag		##<< lag-time between CO2 concentration in chamber and measurement at the sensor[XX]
		,sdFluxRegression = ##<< the standard deviation of the flux by a single regression of CO2 flux
		,sdFluxLeverage = 	##<< the standard deviation of the flux by leverage of starting or end values of the time series
	)
	
}
attr(calcOpenChamberFlux,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
}

estimateLagTime <- function(
		### Estimate the lag-time between CO2 concentration in chamber and measurement at the sensor
		conc	##<< numeric vector of CO2 concentrations 
		,times	##<< times of conc measurements	
){
	##description<< 
	## 
}
attr(estimateLagTime,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
}

regressFlux <- function(
		### Estimate the initial flux by polynomial regression
		conc	##<< numeric vector of CO2 concentrations 
		,times	##<< times of conc measurements	
){
	### numeric vector (2): estimate and its standard deviation
}

sigmaBootLeverage <- function(
		### Estimate uncertainty of regrssion due to leverage of starting and end points. 
		conc	##<< numeric vector of CO2 concentrations 
		,times	##<< times of conc measurements
		,fRegress = regressFlux	##<< function to yield a single flux estimate 
){
	##description<< 
	## 
}
attr(sigmaBootLeverage,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
}

