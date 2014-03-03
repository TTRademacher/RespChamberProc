ds0 <- readDat("tmp/MANIP_Ch1_1.dat", tz="CET")

summary(ds)
ds <- subset(ds0, as.numeric(TIMESTAMP) >= as.numeric(as.POSIXct("2014-03-03 00:00:01 CET")))
ds$CO2_dry <- corrConcDilution(ds, colConc = "CO2_LI840", colVapour = "H2O_LI840")
ds$Pa <- ds$AirPres * 100   # convert hPa to Pa

plot( CO2_LI840 ~ TIMESTAMP, ds)

dsi <- subset( ds, Collar == 1)
timeSec <- as.numeric(dsi$TIMESTAMP) - as.numeric( dsi$TIMESTAMP[1] )
dsi <- subset(dsi, timeSec < 4*60 )

plot( CO2_LI840 ~ TIMESTAMP, dsi)


dsi$CO2_dry <- corrConcDilution(dsi, colConc = "CO2_LI840", colVapour = "H2O_LI840")
dsi$Pa <- dsi$AirPres * 100   # convert hPa to Pa
res <- calcClosedChamberFlux(dsi, colTemp = "AirTemp", colPressure = "Pa", fRegress = c(regressFluxLinear, 
				regressFluxTanh), volume = 0.6*0.6*0.6, isEstimateLeverage = TRUE)

#plotting
conc <- ds$CO2_dry <- corrConcDilution(ds)
times <- ds$TIMESTAMP
times0 <- as.numeric(times) - as.numeric(times[1])
times0Fit <- times0[times0>res["tLag"] ]
#length(times0Fit)
plot( ds$CO2_dry ~ times0, xlab="time (s)", ylab="" ); mtext("CO2_dry (ppm)",2,2,las=0)
abline(v=res["tLag"], col="grey", lty="dotted")
lines( -fitted(attributes(res)$model) ~ times0Fit , col="red" )

ds$TA_Avg <- ds$AirTemp

library(plyr)
#?dlply

dsCollar <- subset( ds, Collar==9)
resL <- dlply( ds, .(Collar), function(dsCollar){
			print( dsCollar$Collar[1] )
			timeSec <- as.numeric(dsCollar$TIMESTAMP) - as.numeric( dsCollar$TIMESTAMP[1] )
			plot( CO2_dry ~ TIMESTAMP, dsCollar )
			dsi <- subset(dsCollar, timeSec < 4*60 )
			if( max(timeSec) < 1*60 ){
				#warning("time series too short")
				return( c( Collar=dsi$Collar[1], rep(NA, 5) ) )
			}
			res <- calcClosedChamberFlux(dsi, colTemp = "AirTemp", colPressure = "Pa", fRegress = c(regressFluxLinear, 
							regressFluxTanh), volume = 0.6*0.6*0.6, isEstimateLeverage = TRUE)
			c( Collar=dsi$Collar[1], res )
		})
