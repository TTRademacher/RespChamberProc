plotCampaignConcSeries <- function(
		### plot all the time series into a series of pdf-pages
		ds					##<< data frame to plot, with collumns "id", "TIMESTAMP" and \code{varName}
		,resL=NULL			##<< list with results of \code{\link{calcClosedChamberFlux}} for each id-subset in ds
		,varName="CO2_dry"	##<< variable to plot
		,idCol="iChunk"		##<< collumn name of identifier of one time series
		,timeCol="TIMESTAMP"##<< collumn name of the time collumn
		,plotsPerPage=64	##<< number of plots per page
		,fileName=paste0(varName,".pdf")	##<< fileName of the output pdf
		,colIds = c()
){
	#iCamp <- 1
	#dss <- subset(ds, campaign==1 & Chamber==1)
	N <- length(unique(ds[[idCol]]))
	ds$id <- factor(ds[[idCol]])	# drop unused factor levels
	if( length(resL) && is.null(names(resL)) ) names(resL) <- levels(ds$id)
	#(as.numeric(unique(ds$id))-1)%/%plotsPerPage+1
	dsp <- cbind(iPage= factor((as.numeric(ds$id)-1)%/%plotsPerPage+1), ds)
	print(paste("Number of pages (each ",plotsPerPage," plots): ", length(unique(dsp$iPage))), sep="" )
	dss <- subset(dsp, iPage==1)
	#ds$H2O_dry <- corrConcDilution(ds, colConc = "H2O_LI840", colVapour = "H2O_LI840"); varName <- "H2O_dry"
	pdf(width=11.7,height=8.3,file=fileName)
	d_ply(dsp, .(iPage), function(dss){
				idsPage <- unique(dss$id)
				print(paste(idsPage, collapse=","))
				# calculate times0
				dss <- dss[order(dss$id, dss[[timeCol]]),]
				dssi <- subset(dss, id==dss$id[1] )
				dss$times0  <- do.call(c, dlply( dss, .(id), function(dssi){
									times <- as.numeric(dssi[[timeCol]])
									times0 <- times - times[1]
								}))
				#dss <- dsc
				p1 <- ggplot( dss, aes_string(x="times0", y=varName) ) + geom_point() + 
						facet_wrap( ~ id, scales="free") +
						theme_gray(base_size=9) + 
						theme()
				if( length(resL) ){
					resLi <- resL[ names(resL) %in% idsPage ]
					# resCalc <- resLi[[1]]
					dfLag <- data.frame( id= names(resLi), tLag=sapply( resLi, function(resCalc){
										if( !inherits(resCalc,"try-error") && length(resCalc$stat) )	resCalc$stat["tLag"] else NA_real_	
									}))
					#idi <- names(resLi)[1]
					dfFitted <- do.call( rbind, lapply( names(resLi), function(idi){
										if( !inherits(resLi[[idi]],"try-error") && length(resLi[[idi]]$model) ){
											dsr <- data.frame( id=idi, fitted=fitted(resLi[[idi]]$model))
											times0 <- subset(dss, id==idi)$times0
											dsr$times0 <- times0[times0 >= dfLag$tLag[dfLag$id==idi]][1:nrow(dsr)]
											dsr
										} else NULL
									}))
					p1 <- p1 + 
							geom_vline( aes(xintercept=tLag), col="darkgrey", linetype="dashed", data=dfLag ) +
							geom_line( aes(y=fitted), col="red", data=dfFitted ) +
							theme()
				}
				print(p1)
				#p1
			})
	dev.off()
}
attr( plotCampaignConcSeries, "ex") <- function(){
	load(file.path(projectDir,"SMANIEconc.Rd"))
	ds0 <- ds <- SMANIEconc
	#tmp <- subset(ds0, campaign==5 )
	#ds0 <- ds <- subset(SMANIEconc, id %in% levels(SMANIEconc$id)[64+1:8] ) testing
	#ds0 <- ds <- subset(SMANIEconc, id %in% c(levels(SMANIEconc$id)[64+1:3],"1.2.19","1.3.2")  ) #testing for fit errors
	.tmp.f <- function(){
		ds <- subset(ds0,  !id %in% c(
						"1.2.19"	# nonunique times
						,"1.3.2"	# exactly zero concentrations all
				) )
	}
	#dss <- subset(ds, id == ds$id[1] )
	#dss <- subset(ds, id == "1.3.2" )
	# ds <- ds[-(1:33550), ]
	.createresCalcClosedChamberFlux <- function(){
		resCalcClosedChamberFlux <- resL <- dlply( ds, .(id), function(dss){
					cat( ",", levels(dss$id)[dss$id[1]] )
					res2 <- res <- try( 
							res <- calcClosedChamberFlux( dss, colTemp="AirTemp", volume = 0.6*0.6*0.6
									, debugInfo = list(useOscarsLagDectect = TRUE, omitAutoCorrFit=TRUE, omitEstimateLeverage = TRUE)
							)
							,silent=TRUE)
				})
		save( resCalcClosedChamberFlux, file=file.path(projectDir,"resCalcClosedChamberFlux_SMANIEconc.Rd"))
	}
	load( file.path(projectDir,"resCalcClosedChamberFlux_SMANIEconc.Rd"))
	resL <- resCalcClosedChamberFlux
	names(resL)[ sapply( resL, inherits,  "try-error") ]	# "1.2.19"  "1.3.2"   "1.3.23"  "2.2.270" "2.2.277"
	#
	fileName=file.path(projectDir,paste0("SMANIEconc_CO2_dry.pdf"))	##<< fileName
	plotCampaignConcSeries(ds, resL=resL, fileName=fileName)
}

