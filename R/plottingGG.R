plotCampaignConcSeries <- function(
		### get a series of ggplots of the time series and its fits
		ds					##<< data frame to plot, with collumns \code{idCol}, \code{timeCol} and \code{varName}
		,resL=NULL			##<< list with results of \code{\link{calcClosedChamberFlux}} for each id-subset in ds
		,varName="CO2_dry"	##<< variable to plot
		,idCol="iChunk"		##<< collumn name of identifier of one time series
		,timeCol="TIMESTAMP"##<< collumn name of the time collumn
		,plotsPerPage=64	##<< number of plots per page
		,fileName=""		##<< if non-zero length string, the fileName where all plots are printed to  #paste0(varName,".pdf")
		,colIds = c()
		,ggplotList=c()		##<< list added to each ggplot. 
){
	#iCamp <- 1
	#dss <- subset(ds, campaign==1 & Chamber==1)
	N <- length(unique(ds[[idCol]]))
	ds$id <- factor(ds[[idCol]])	# drop unused factor levels
	if( length(resL) && is.null(names(resL)) ) names(resL) <- levels(ds$id)
	#(as.numeric(unique(ds$id))-1)%/%plotsPerPage+1
	dsp <- cbind(iPage= factor((as.numeric(ds$id)-1)%/%plotsPerPage+1), ds)
	print(paste("Number of pages (each ",plotsPerPage," plots): ", length(unique(dsp$iPage))), sep="" )
	dss <- dsp[dsp$iPage==1,]
	#ds$H2O_dry <- corrConcDilution(ds, colConc = "H2O_LI840", colVapour = "H2O_LI840"); varName <- "H2O_dry"
	plotList <- dlply(dsp, "iPage", function(dss){
				idsPage <- unique(dss$id)
				print(paste(idsPage, collapse=","))
				# calculate times0
				dss <- dss[order(dss$id, dss[[timeCol]]),]
				dssi <- dss[dss$id==dss$id[1],]
				dss$times0  <- do.call(c, dlply( dss, "id", function(dssi){
									times <- as.numeric(dssi[[timeCol]])
									times0 <- times - times[1]
								}))
				#dss <- dsc
				p1 <- ggplot( dss, aes_string(x="times0", y=varName) ) + geom_point(shape=1) + 
						facet_wrap( ~id, scales="free") +
						theme_bw(base_size=9) 
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
											times0 <- dss[dss$id==idi,"times0"]
											dsr$times0 <- times0[times0 >= dfLag$tLag[dfLag$id==idi]][1:nrow(dsr)]
											dsr
										} else NULL
									}))
					p1 <- p1 + 
							geom_vline( aes_string(xintercept="tLag"), col="darkgrey", linetype="dashed", data=dfLag ) +
							geom_line( aes_string(y="fitted"), col="red", data=dfFitted ) +
							theme()
				}
				p1
			})
	if( nzchar(fileName) ){
		pdf(width=11.7,height=8.3,file=fileName)
		for( i in seq_along(plotList)){
			cat(",",i)
			print(plotList[[i]])
		}
		dev.off()
		if( isVerbose ) message("printed plots to file ",fileName)
	}
	plotList
}


