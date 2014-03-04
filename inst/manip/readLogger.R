ds0 <- readDat("tmp/MANIP_Ch1_1.dat", tz="CET")

ds <- subset(ds0, as.numeric(TIMESTAMP) >= as.numeric(as.POSIXct("2014-03-03 00:00:01 CET")))
ds$CO2_dry <- corrConcDilution(ds, colConc = "CO2_LI840", colVapour = "H2O_LI840")
ds$Pa <- ds$AirPres * 100   # convert hPa to Pa

#dsChunksClean <- subsetContiguous(ds)
dsChunksClean <- subsetContiguous(ds, fIsBadChunk=function(dsi){ var(dsi$CO2_dry)==0})
length(unique(dsChunksClean$iChunk))

# plot the time series
library(ggplot2)
p1 <- ggplot( dsChunksClean, aes(x=TIMESTAMP, y=CO2_dry) ) + geom_point() + facet_wrap( ~ iChunk, scales = "free")

#-- calculate flux and extract environmental conditions, may be parallelized
library(plyr)
library(doSNOW)
nNode = 2	# number of processors
cl = makeCluster(nNode)		
registerDoSNOW(cl)
clusterEvalQ(cl, library(RespChamberProc))		# functions need to be loaded on remote hosts

#dsi <- subset( dsChunksClean, iChunk=iChunk[1] )
system.time(res <- ddply( 	dsChunksClean, .(iChunk), function(dsi){
			collar <- dsi$Collar[1] 
			iChunk = dsi$iChunk[1]
			print( paste(iChunk, dsi$TIMESTAMP[1]," Collar: ",collar) )
			#plot( CO2_dry ~ timeSec, dsi )
			res <- calcClosedChamberFlux(dsi, colTemp = "AirTemp", colPressure = "Pa", fRegress = c(regressFluxLinear, 
							regressFluxTanh), volume = 0.6*0.6*0.6, isEstimateLeverage = TRUE)
			# get additional environmental variables at the initial time
			dsiInitial <- dsi[ 1, ,drop=FALSE]
			cbind( data.frame( time=dsiInitial[,"TIMESTAMP"], collar=collar, CO2_flux=res[1], CO2_flux_sd=res[2] )
				, dsiInitial[,c("Chamber","AirTemp","AirPres","PAR","BodyTemp","SurTemp","SoilTemp","SoilMoist")] )
}
, .parallel=TRUE
))

#stopCluster(cl)

res




