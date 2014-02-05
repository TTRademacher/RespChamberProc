densityCorrFlux <- function(
	### Calculate density corrected concentration		
	ds 					##<< data frame with each row one observations, and respective columns 
	,volume=1			##<< volume of the chamber in m3
	,colConc="CO2_Avg"	##<< column name of CO2 concentration [Amount of substance]
	,colTemp="TA_Avg"	##<< column name of air temperature inside chamber [°C]  
	,colPressure="Pa"	##<< column name of pressure inside chamber [kPa]
){
	##details<< 
	## XX
	R <- 8.3144621	# universal gas constant in [J/K/mol]
	ds[,colConc] * ds[,colTemp] * ds[,colPressure] * ds[,colPressure] * volume * R
	### numeric vector (nrow ds):  corrected concentration [Amount of substance]
} 
attr(densityCorrFlux,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
	ds$CO2_denC <- densityCorrFlux(ds)	
}


dilutionCorrFlux <- function(
		### Calculate concentration corrected for dilution with water vapour		
		ds 					##<< data frame with each row one observations, and respective columns 
		,colConc="CO2_Avg"	##<< column name of CO2 concentration [Amount of substance]
		,colVapour="H20_Avg"	##<< column name of CO2 concentration [Amount of substance]
){
	##details<< 
	## XX
	ds[,colConc] * ds[,colVapour] 
	### numeric vector (nrow ds):  corrected concentration [Amount of substance]
} 
attr(dilutionCorrFlux,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
	ds$CO2_dilC <- dilutionCorrFlux(ds)	
}


leakageCorrFlux <- function(
	### Calculate concentration corrected for dilution with water vapour		
	conc 				##<< numeric vecgor: concentration [Amount of substance]
	,pLeakage=1			##<< numeric scalar: XX
){
	##details<< 
	## XX
	colc * pLeakage
	### numeric vector (length conc): corrected concentration [Amount of substance]
} 
attr(leakageCorrFlux,"ex") <- function(){
	data(chamberLoggerEx1s)
	ds <- chamberLoggerEx1s
	ds$CO2_leakC <- leakageCorrFlux(ds)	
}

