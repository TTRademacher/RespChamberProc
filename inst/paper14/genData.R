projectDir = "inst/paper14"
library(stringr)
library(ggplot2)
library(plyr)


genSMANIEControl <- function(){
	#data(collarCodes)   # merge on conc takes too long
	rawDataPaths="M:/data/SMANIE/Fluxes/raw"
	campaignDirs = list.files( rawDataPaths, pattern="SMANIP_\\d", include.dirs=TRUE)
	campaigns <- as.integer(str_extract(campaignDirs, "\\d"))
	#iCamp <- 1
	dsl <- data.frame()
	#iCamp <- 2
	# logger data has in part records from previous measurement campaigns, but take only the recent ones
	campaignStarts <- as.POSIXctUTC(c("2014-03-17", "2014-04-14","2014-05-06","2014-05-27","2014-06-23"))
	# some logger data also appears in different files (chamber 3 in other chamber files
	# make sure to take only the recent ones of the correct logger (from file name)
	for( iCamp in seq_along(campaignDirs) ){
	#dsl <-  lapply( seq_along(campaignDirs), function(iCamp){
				cat(",",iCamp)
				fNames <- dir( file.path(rawDataPaths,campaignDirs[iCamp]), pattern="\\.dat$", full.names=TRUE )
				#fName <- fNames[1]
				for( fName in fNames ){
				#dsd <- do.call(rbind, lapply( fNames, function(fName){
							dsa <- readDat(fName)
							#gsub("Chamber(\\d+).*","\\1",c("Chamber1","Chamber12_"))
							iChamber <- gsub("(.*Chamber)(\\d+)(.*)","\\2",fName)
							dsp <- subset(dsa, TIMESTAMP >= campaignStarts[iCamp] )
							ds <- subset(dsp, Chamber == iChamber)
							if( nrow(ds)==0  ) stop("encountered zero data frame. Check file names and patterns.")
							#trace(subsetContiguous, recover)	#untrace(subsetContiguous)
							ds$CO2_dry <- corrConcDilution(ds, colConc = "CO2_LI840", colVapour = "H2O_LI840")
							dsc <- subsetContiguous(ds )		# CO2_dry must be defined
							dscc <- cbind( campaign=iCamp, dsc )
							dsl <- rbind( dsl, dscc)
							if( any(table(interaction(dsl$RECORD,dsl$Chamber, drop=TRUE)) != 1) ) stop("encountered double record")
							dsl
						}
						#))
				#if( any(table(dsd$RECORD) != 1) ) stop("encountered double record")
			}
			#)
	#ds <- do.call( rbind, dsl )
	ds <- dsl
	if( any(table(interaction(ds$RECORD,ds$Chamber, drop=TRUE)) != 1) ) stop("encountered double record")
	ds <- ds[order(ds$campaign, ds$Chamber, ds$iChunk), ]
	ds$id <- interaction( ds$campaign, ds$Chamber, ds$iChunk, drop=TRUE, lex.order = TRUE)
	ds$Pa <- ds$AirPres * 1000		# convert kPa to Pa
	SMANIEconc <- ds 
	# TAKE care overwrites
	save( SMANIEconc, file=file.path(projectDir,"SMANIEconc.Rd"))
}

.plotCampaignConcSeries <- function(varName="CO2_dry"){
	load(file.path(projectDir,"SMANIEconc.Rd"))
	ds <- SMANIEconc
	#iCamp <- 1
	#dss <- subset(ds, campaign==1 & Chamber==1)
	N <- length(unique(ds$id))
	ds <- cbind(id64= factor(floor(as.numeric(ds$id)/64)), ds)
	print(paste("Number of pages (each 64 plots):", length(unique(ds$id64))) )
	dss <- subset(ds, id64==1)
	#ds$H2O_dry <- corrConcDilution(ds, colConc = "H2O_LI840", colVapour = "H2O_LI840"); varName <- "H2O_dry"
	pdf(width=8.3,height=11.7,file=file.path(projectDir,paste0("SMANIEconc_",varName,".pdf")))
	d_ply(ds, .(id64), function(dss){
				print(paste(unique(dss$id), collapse=","))
				#dss <- dsc
				p1 <- ggplot( dss, aes(x=TIMESTAMP, y=get(varName)) ) + geom_point() + 
						facet_wrap( ~ id, scales="free")
				print(p1)
				#p1
			})
	dev.off()
}

.inspectTimeSeries <- function(){
	goodPeriods <- c(
	"1.2.19" = c()		# double series
	,"1.2.82" = c()		# nonsensical
	,"1.2.83" = 1:180		# long lag time, check doubles 
	,"1.2.87" = 1:180		# good, small gradient, also 88
	,"1.3.18" = c()			# not ventilated, nonsensical
	,"2.2.5" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	,"1.2.19" = 10:200
	)
	data(collarCodes)
	load(file.path(projectDir,"SMANIEconc.Rd"))
	ds <- SMANIEconc
	idi <- "2.2.5"
	dss <- merge( subset(ds, id==idi), collarCodes )
	tmp <- diff(dss$RECORD)
	dss <- dss[1:min(which(tmp < 1)), ]
	head(dss,5)
	plotResp(dss)
	res2 <- res <- calcClosedChamberFlux( dss, colTemp="AirTemp", volume = 0.6*0.6*0.6, debugInfo = list(useOscarsLagDectect = TRUE))
	#res1 <- res <- calcClosedChamberFlux( dss, colTemp="AirTemp", volume = 0.6*0.6*0.6)
	plotResp(dss, res)
	
}


