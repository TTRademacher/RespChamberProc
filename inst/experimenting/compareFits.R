tmp <- read.csv("tmp/Flux_R.csv")
i <- 1
dsAll <- do.call( rbind, lapply( 1:24, function(i){
			data.frame( series=i, time=tmp[,1], CO2 = tmp[,i+1], water = tmp[,i+1+24] ) 
		}))
dsAll$CO2_dry <- corrConcDilution(dsAll, "CO2", "water")
dsAll$Pa <- 101*1000	# sea level pressure
dsAll$TA_Avg <- 21		# room temperature
dsAll$TIMESTAMP <- ds$time

.tmp.f <- function(){
	library(ggplot2)
	p1 <- ggplot( dsAll, aes(x=time, y=CO2_dry) ) + geom_point() + facet_wrap( ~ series); p1
}

ds <- subset(dsAll, series==17)		# saturating consumption
ds <- subset(dsAll, series==15)		# saturating production


dss <- subset(dsAll, series != 18)
library(plyr)
resFluxesL <- dlply( dss, .(series), function(ds){
			cat(ds$series[1],",")
			resi <- calcClosedChamberFlux(ds, fRegress=regressFluxTanh)
			c( series = ds$series[1], resi )
		})
resFluxesTanh <- do.call(rbind, resFluxesL)

resFluxesL <- dlply( dss, .(series), function(ds){
			resi <- calcClosedChamberFlux(ds, fRegress=regressFluxSquare)
			c( series = ds$series[1], resi )
		})
resFluxesPoly <- do.call(rbind, resFluxesL)

(resFluxesPoly - resFluxesTanh ) / resFluxesTanh 


ds <- subset(dss, series==2)		# differences between tanh and poly
ds <- subset(dss, series==5)		# differences between tanh and poly
ds <- subset(dss, series==13)		# differences between tanh and poly
#ds <- subset(dss, series==16)		# differences between tanh and poly
ds <- subset(dss, series==23)		# differences between tanh and poly
