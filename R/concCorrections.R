
corrConcDilution <- function(
  ### Calculate concentration corrected for dilution with water vapor		
  ds 					##<< data frame with each row one observations, and respective columns 
  ,colConc="CO2_Avg"	##<< column name of CO2 concentration [ppm]
  ,colVapour="H20_Avg"	##<< column name of CO2 concentration [ppt]
){
  ##details<< 
  ##  If CO2 concentration is measured per moist air, this function will calculate the concentration\
  ## per dry air.
  
  ##references<<
  ## LI-COR, Application Note 129. The Importance of Water Vapor Measurements and Corrections. LI-COR, Inc., 4421 Superior Street, Lincoln, NE 68504, USA.
  ds[,colConc]/(1- (ds[,colVapour])/1000)
  ### numeric vector (nrow ds):  concentration of CO2 per dry air [ppm]
} 
attr(corrConcDilution,"ex") <- function(){
  data(chamberLoggerEx1s)
  ds <- chamberLoggerEx1s
  ds$CO2_dry <- corrConcDilution(ds)	
}



corrFluxDensity <- function(
	### Calculate density corrected concentration		
	ds 					##<< data frame with each row one observations, and respective columns 
	,volume=1			##<< volume of the chamber in m3
	,colConc="CO2_Dry"##<< column name of CO2 concentration per dry air [Amount of substance]
	,colTemp="TA_Avg"	##<< column name of air temperature inside chamber [°C]  
	,colPressure="Pa"	##<< column name of pressure inside chamber [kPa]
){
	##details<< 
	## XX
	R <- 8.3144621	# universal gas constant in [J/K/mol]
	ds[,colConc] * ds[,colTemp] * ds[,colPressure] * ds[,colPressure] * volume * R
	### numeric vector (nrow ds):  corrected concentration [Amount of substance]
} 
attr(corrFluxDensity,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
	ds$CO2_denC <- corrConcDensity(ds)	
}


corrConcLeakage <- function(
	### Calculate concentration corrected for dilution with water vapour		
	conc 				##<< numeric vecgor: concentration [Amount of substance]
	,pLeakage=1			##<< numeric scalar: XX
){
	##details<< 
	## XX
	colc * pLeakage
	### numeric vector (length conc): corrected concentration [Amount of substance]
} 
attr(corrConcLeakage,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
	ds$CO2_leakC <- corrConcLeakage(ds)	
}

