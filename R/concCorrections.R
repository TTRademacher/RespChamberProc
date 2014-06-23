
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
	### Calculate total flux inside the chamber from flux per amount of air
  	CO2_molarFlux 		##<< numeric vector of rate of changes in CO2 dry molar fraction [ppm/sec]
	,volume=1			##<< numeric scalar: volume of the chamber in [m3]
	,temp=20	     	##<< numeric vector: temperature inside chamber [degC]  
	,pressure=101325 	##<< numeric vector: pressure inside chamber [Pa]
){
	##details<< 
	## The amount of air is determined by the volumen, pressure and temperature (universal gas law)
	R <- 8.3144621	# universal gas constant in [J/K/mol]  (J=Pa * m3)  
	CO2_molarFlux * pressure * volume  / (temp+273.15) / R 
	### numeric vector (nrow ds):  CO2 flux [mumol CO2 /s]
} 
attr(corrFluxDensity,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
  	CO2_molarFlux <- ds$CO2_Avg   # just to have a column (here concentrations), the flux should be calculated properly before
	ds$CO2_flux <- corrFluxDensity(CO2_molarFlux, pressure=101*1000) #kPa converted to Pa	
  	plot(ds$CO2_flux)
}


corrFluxLeakage <- function(
	### Calculate flux corrected for leakage of the chamber	
	conc 				##<< numeric vector: concentration [Amount of substance]
	,pLeakage=1			##<< numeric scalar: fraction of concentration lost by leakage
){
	##details<< 
	## XX
	conc * pLeakage
	### numeric vector (length conc): corrected concentration [Amount of substance]
} 
attr(corrFluxLeakage,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
	ds$CO2_leakC <- corrFluxLeakage(ds$CO2_Avg)	
}

