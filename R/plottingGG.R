plotCampaignConcSeries <- function(
		### get a series of ggplots of the time series and its fits
		ds					##<< data frame to plot, with collumns \code{idCol}, \code{timeCol} and \code{varName}
		,resL=NULL			##<< list with results of \code{\link{calcClosedChamberFlux}} for each id-subset in ds
		,varName="CO2_dry"	##<< variable to plot
		,idCol="iChunk"		##<< collumn name of identifier of one time series
		,timeCol="TIMESTAMP"##<< collumn name of the time collumn
		,qualityFlag=0		##<< vector of length nrow(ds) of a quality flag. For chunks where
			## this flag is not 0, subplots are dimmed.
		,fTextBR=NULL 		##<< function(resFit) to add some text to the bottom right of the plot, by default the rSquared from fitting object
		,fTextTL=NULL 		##<< function(resFit) to add some text to the top left of the plot, by default the autocorrelation from fitting object
		,fTextTR=NULL 		##<< function(resFit) to add a second text to the top right of the plot, by default the provided quality flag, where it is not zero
		,plotsPerPage=64	##<< number of plots per page
		,fileName=""		##<< if non-zero length string, the fileName where all plots are printed to  #paste0(varName,".pdf")
		,colIds = c()
		,ggplotList=c()		##<< list added to each ggplot.
		,isVerbose=TRUE
){
	# do not clutter the function declaration with these long defaults, better assing in body: 
	if( !length(fTextBR) ) fTextBR <- function(resFit){ if( is.finite(resFit$stat["r2"])) format(resFit$stat["r2"],digits=2) else ""} 
	if( !length(fTextTL) ) fTextTL <- function(resFit){ if( is.finite(resFit$stat["autoCorr"])) format(resFit$stat["autoCorr"],digits=2) else ""} 
	if( !length(fTextTR) ) fTextTR <- function(resFit){ if( length(resFit$qf) && resFit$qf != 0 ) as.character(resFit$qf) else ""}  
	#iCamp <- 1
	#dss <- subset(ds, campaign==1 & Chamber==1)
	N <- length(unique(ds[[idCol]]))
	ds$id <- factor(ds[[idCol]])	# drop unused factor levels
	if( length(resL) && is.null(names(resL)) ) names(resL) <- levels(ds$id)
	if( length(qualityFlag)==1L ) qualityFlag <- rep(qualityFlag, N)
	if( length(qualityFlag) != N ) warning("provided quality flag with different length than than number of measurement cycles in data.")
	if( length(resL) ){
		if( length(resL) != N ) warning("provided fitting results with different length than number of measurement cycles in data.")
		dsQf <- data.frame(id=names(resL), qf=factor(qualityFlag))
		ds <- merge(ds, dsQf)
	} else {
		ds$qf <- factor(0)
	}
	uniqueQf <- unique(qualityFlag)
	colCodes <- rep("lightgray", length(uniqueQf))
	colCodes[1] <- "black"
	colCodes[uniqueQf==10] <- "darkgrey"
	#(as.numeric(unique(ds$id))-1)%/%plotsPerPage+1
	dsp <- cbind(iPage= factor((as.numeric(ds$id)-1)%/%plotsPerPage+1), ds)
	message(paste("Number of pages (each ",plotsPerPage," plots): ", length(unique(dsp$iPage))), sep="" )
	dss <- dsp[dsp$iPage==1,]
	#ds$H2O_dry <- corrConcDilution(ds, colConc = "H2O_LI840", colVapour = "H2O_LI840"); varName <- "H2O_dry"
	if( length(resL) ) resL <- structure(lapply( seq_along(resL), function(i){ resLi <- resL[[i]]; resLi$qf <- qualityFlag[i]; resLi }), names=names(resL))
	plotList <- dlply(dsp, "iPage", function(dss){
				idsPage <- unique(dss$id)
				if(isVerbose) message(paste(idsPage, collapse=","))
				# calculate times0
				dss <- dss[order(dss$id, dss[[timeCol]]),]
				dssi <- dss[dss$id==dss$id[1],]
				dss$times0  <- do.call(c, dlply( dss, "id", function(dssi){
									times <- as.numeric(dssi[[timeCol]])
									times0 <- times - times[1]
								}))
				#dss <- dsc
				p1 <- ggplot( dss, aes_string(x="times0", y=varName) ) + 
						geom_point(shape=1, aes_string(col="qf"), na.rm=TRUE) + 
						facet_wrap( ~id, scales="free") +
						scale_color_manual(values=colCodes, guide = FALSE) +
						theme_bw(base_size=9) + 
						theme(panel.grid=element_blank())
						#theme(panel.grid.minor=element_blank())
				if( length(resL) ){
					iiChunk <- which(names(resL) %in% idsPage)
					resLi <- resL[iiChunk]
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
							#resFit <- resLi[[54]]
							tmp <- sapply(resLi, fTextBR )
							dfTextBR <- data.frame(id=names(tmp), text=tmp, row.names = NULL)
							tmp <- sapply(resLi, fTextTL )
							dfTextTL <- data.frame(id=names(tmp), text=tmp, row.names = NULL)
							tmp <- sapply(resLi, fTextTR )
							dfTextTR <- data.frame(id=names(tmp), text=tmp, row.names = NULL)
							p1 <- p1 + 
							geom_vline( data=dfLag, aes_string(xintercept="tLag"), col="darkgrey", linetype="dashed", na.rm=TRUE ) +
							geom_line( data=dfFitted, aes_string(y="fitted"), col="red", na.rm=TRUE  ) +
							geom_text( data=dfTextBR, aes_string(label="text"), x=+Inf, y=-Inf, hjust=1.05, vjust=0, na.rm=TRUE) +
							geom_text( data=dfTextTL, aes_string(label="text"), x=-Inf, y=+Inf, hjust=0, vjust=1, na.rm=TRUE) +
							geom_text( data=dfTextTR, aes_string(label="text"), x=+Inf, y=+Inf, hjust=1, vjust=1, na.rm=TRUE) +
							theme()
				}
				p1
			})
	if( nzchar(fileName) ){
		pdf(width=11.7,height=8.3,file=fileName)
		on.exit(dev.off())
		for( i in seq_along(plotList)){
			if( isVerbose ) message(",",i, appendLF=FALSE)
			print(plotList[[i]])
		}
		if( isVerbose ) message("\nprinted plots to file ",fileName) 
	}
	plotList
}


