tmp <- read.csv("tmp/Flux_R.csv")
i <- 1
dsAll <- do.call( rbind, lapply( 1:24, function(i){
			data.frame( series=i, time=tmp[,1], CO2 = tmp[,i+1], water = tmp[,i+1+24] ) 
		}))
dsAll$CO2_dry <- corrConcDilution(dsAll, "CO2", "water")
dsAll$Pa <- 101*1000	# sea level pressure
dsAll$TA_Avg <- 21		# room temperature
dsAll$TIMESTAMP <- dsAll$time

.tmp.f <- function(){
	library(ggplot2)
	p1 <- ggplot( dsAll, aes(x=time, y=CO2_dry) ) + geom_point() + facet_wrap( ~ series, scales="free"); p1
}





dss <- subset(dsAll, series != 18)

ds <- subset(dss, series==17)		# saturating consumption
ds <- subset(dss, series==15)		# saturating production
res <- calcClosedChamberFlux(ds )
# plot( CO2~time, ds )

ds <- subset(dss, series==3)		# saturating production


library(plyr)
resFluxesL <- dlply( dss, .(series), function(ds){
			cat(ds$series[1],",")
			#trace(calcClosedChamberFlux, recover )	#untrace(calcClosedChamberFlux)
			#trace(regressFluxTanh, recover )	#untrace(regressFluxTanh)
			resi <- calcClosedChamberFlux(ds, fRegress=c(regressFluxLinear))
			c( series = ds$series[1], resi )
		})
resFluxesLin <- do.call(rbind, resFluxesL)

resFluxesL <- dlply( dss, .(series), function(ds){
			cat(ds$series[1],",")
			#trace(calcClosedChamberFlux, recover )	#untrace(calcClosedChamberFlux)
			#trace(regressFluxTanh, recover )	#untrace(regressFluxTanh)
			resi <- calcClosedChamberFlux(ds, fRegress=c(regressFluxTanh))
			c( series = ds$series[1], resi )
		})
resFluxesTanh <- do.call(rbind, resFluxesL)

resFluxesL <- dlply( dss, .(series), function(ds){
			cat(ds$series[1],",")
			resi <- calcClosedChamberFlux(ds, fRegress=c(regressFluxSquare))
			c( series = ds$series[1], resi )
		})
resFluxesPoly <- do.call(rbind, resFluxesL)

resFluxesL <- dlply( dss, .(series), function(ds){
			cat(ds$series[1],",")
			#trace(regressFluxExp, recover )	#untrace(regressFluxExp)
			resi <- calcClosedChamberFlux(ds, fRegress=c(regressFluxExp))
			c( series = ds$series[1], resi )
		})
resFluxesExp <- do.call(rbind, resFluxesL)

.tmp.f <- function(){
	resFluxesL <- dlply( dss, .(series), function(ds){
				resi <- calcClosedChamberFlux(ds)
				c( series = ds$series[1], resi )
			})
	resFluxesAIC <- do.call(rbind, resFluxesL)
}

AIC <- cbind(resFluxesLin[,"AIC"], resFluxesPoly[,"AIC"], resFluxesExp[,"AIC"], resFluxesTanh[,"AIC"] )
apply( AIC, 1, which.min )



(resFluxesPoly - resFluxesTanh ) / resFluxesTanh
(resFluxesExp - resFluxesTanh ) / resFluxesTanh

plot( abs(flux) ~ series, resFluxesAIC, pch="x" )
points( abs(flux) ~ series, resFluxesLin, col="grey" )
points( abs(flux) ~ series, resFluxesExp, col="red" )
points( abs(flux) ~ series, resFluxesTanh, col="blue" )
points( abs(flux) ~ series, resFluxesPoly, col="maroon" )

ds <- subset(dss, series==2)		# differences between tanh and poly
ds <- subset(dss, series==5)		# differences between tanh and poly
ds <- subset(dss, series==13)		# differences between tanh and poly
#ds <- subset(dss, series==16)		# differences between tanh and poly
ds <- subset(dss, series==23)		# differences between tanh and poly

ds <- subset(dss, series==15)		# differences between tanh and exp
ds <- subset(dss, series==17)		# differences between tanh and exp # but strange points after lag-time
ds <- subset(dss, series==21)		# differences between tanh and exp

ds <- subset(dss, series==6)		# tanh fitted but not exp

#trace(calcClosedChamberFlux, recover )	#untrace(calcClosedChamberFlux)
#trace(regressFluxExp, recover )	#untrace(regressFluxExp)
rExp <- calcClosedChamberFlux(ds, fRegress=c(regressFluxExp)); mExp <- attr(rExp,"model")
#trace(regressFluxTanh, recover )	#untrace(regressFluxTanh)
rTanh <- calcClosedChamberFlux(ds, fRegress=c(regressFluxTanh)); mTanh <- attr(rTanh,"model")
rLin <- calcClosedChamberFlux(ds, fRegress=c(regressFluxLinear)); mLin <- attr(rLin,"model")
rPoly <- calcClosedChamberFlux(ds, fRegress=c(regressFluxSquare)); mPoly <- attr(rPoly,"model")
#rAIC <- calcClosedChamberFlux(ds ); mAIC <- attr(rAIC,"model")
rbind(rLin, rExp, rTanh, rPoly)
qqnorm(resid(mTanh, type="normalized")); abline(0,1)
qqnorm(resid(mExp, type="normalized")); abline(0,1)
qqnorm(resid(mPoly, type="pearson")); abline(0,1)
qqnorm(resid(mLin, type="pearson")); abline(0,1)
plot( CO2_dry ~ time, ds)
points( CO2_dry ~ time, ds[ds$time <= rExp["tLag"], ], col="lightgrey")
points( CO2_dry ~ time, ds[ds$time > 5*rExp["tLag"], ], col="lightgrey")
abline(v=rExp["tLag"])
lines(fitted(mExp) ~ I(ds$time[ds$time > rExp["tLag"] ]))
lines(fitted(mTanh) ~ I(ds$time[ds$time > rExp["tLag"] ]), col="blue")
#lines(fitted(mAIC) ~ I(ds$time[ds$time > rExp["tLag"] ]), col="maroon")
lines(fitted(mPoly) ~ I(ds$time[ds$time > rExp["tLag"] ]), col="green")

plot( resid(mTanh) ~ I(ds$time[ds$time > rExp["tLag"] ]))
points( resid(mExp) ~ I(ds$time[ds$time > rExp["tLag"] ]), col="blue")
qqnorm(resid(mTanh)); abline(0,1)

windows()
qqnorm(resid(mExp)); abline(0,1)
qqline( resid(mExp), distribution=function(p){ qexp(p)}, col="blue" )
plot(density(resid(mExp)))
x <- seq(-1.5,1.5,length.out=30); 
lines( dnorm( x, sd=sd(resid(mExp)) ) ~ x, col="red")
lines( dexp( x ) ~ x, col="red")
lines(density(resid(mTanh)), col="blue")

acf(resid(mTanh))
acf(resid(mExp))


ds <- subset(dss, series==17)		# large difference between methods
