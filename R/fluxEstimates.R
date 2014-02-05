
calcOpenChamberFlux <- function(
	### Calculate flux and its uncertainties.
	ds						##<< data.frame with concentration and time column
	,colConc="CO2_corr"		##<< column name of CO2 concentration [Amount of substance]
	,colTime="TIMESTAMP"	##<< column name of time
	,rescolFlux="flux"		##<< result column name for calculated flux 
	,rescolSdFlux=paste("sd",rescolFlux,sep="")		##<< result column name for calculated standard deviation of flux 
){
	##details<< 
	## 
	### data.frame with columns flux and sdFlux
	
}
attr(calcOpenChamberFlux,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
}

estimateLagTime <- function(
		### Estimate the lag-time between CO2 concentration and IRGA measurement
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
		### estimate uncertainty of regrssion due to leverage of starting and end points 
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





