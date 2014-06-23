readDat <- function(
	### Read a data loggger .dat file into a data-frame		
	fName					##<< scalar string: file name
	,nRowsFileInfo = 1		##<< integer scalar: number of lines before column information
	,nRowsColInfo = 2		##<< integer vector: number of lines with column description
	,sep=","				##<< column separator
	,...					##<< further arguments to
	,colClasses = rep(NA, ncol(colInfo))	##<< see \code{link{read.table}}
	,colsTimeStamp=1		##<< integer vector: colums with time stamp column (will be set to POSIXct
	,formatTS="%Y-%m-%d %H:%M:%S"	##<< format of the timestamp columns, see \code{\link{strptime}}, e.g. 
	,tz="UTC"				##<< specify a time zone when converting to POSIXct, default: current local e.g CET, UTC
	,na.strings=c('','NA','NAN','"NAN"')
){
	##details<< 
	## Assumes that there
	setClass("myDate", where=globalenv())
	setAs("character","myDate", function(from) as.POSIXct(from, format=formatTS, tz=tz), where=globalenv() )
	fileInfo <- readLines(fName, n=nRowsFileInfo )
	colInfo <- read.table(fName, header=TRUE, skip=nRowsFileInfo, nrows=max(1,nRowsColInfo), sep=sep, na.strings=na.strings)
	colClasses[colsTimeStamp] <- "myDate"
	rawData <- read.table(fName, header=FALSE, skip=nRowsFileInfo+1+nRowsColInfo, sep=sep, na.strings=na.strings, ...
					,colClasses=colClasses)
	colnames(rawData) <- colnames(colInfo)	
	#plot( CO2_Avg ~ TIMESTAMP, data=rawData )
	attr(rawData,"fileInfo") <- fileInfo
	attr(rawData,"colInfo") <- colInfo
	rawData
}
attr(readDat,"ex") <- function(){
	fName <- system.file("genData/chamberLoggerEx1_short.dat", package = "RespChamberProc")
	if( nzchar(fName) ){
		ds <- readDat(fName)
	}
}


subsetContiguous <- function(
	### Get contiguous subsets 
	ds						##<< data.frame of measurements 
	,colTime="TIMESTAMP"	##<< column name that of time (POSIXct)
	,colIndex="Collar"		##<< column name of index variable (factor or integer)
	,gapLength=20			##<< minimal length of a gap between subsets (seconds)
	,minNRec=20				##<< minimum number of records within one contiguous subset
	,minTime=60				##<< minimum length of time that a contiguous subsets covers
	,indexNA=0				##<< value of the index column, that signifies records not to use
	,fIsBadChunk=function(dsi) FALSE	
	### additional function taking and subset and returning a boolean value whether its a chunk to be omitted
){
	##details<< 
	## The time series in logger data consists of several chunks of concentration measurments.
	## In order to determine these chunks, either a change in an index variable (input by between the
	## measurements) or a gap in time is used.
	timeDiff <- diff(as.numeric(ds[,colTime]))
	iGaps <-  which( timeDiff > gapLength)   
	iCollarChanges <- which( diff(as.numeric(ds[,colIndex])) != 0 )
	iChunks <- c( 1, sort(union(iCollarChanges, iGaps)), nrow(ds) ) 
	dsChunksL <- lapply( 2:length(iChunks), function(i){
				cbind( iChunk=i-1, ds[(iChunks[i-1]+1):(iChunks[i]), ,drop=FALSE])
			}) 
	##details<<
	## Between the actural series of measurements, the logger may record sparse data.
	## These chunks are indicated by value \code{indexNA} in the index column or
	## by shortness of the series. Only chunks with at least \code{minNRec} records and at least 
	## \code{minTime} seconds are reported. Others are neglected.
	#
	#sapply( dsChunksL, nrow )
	#sapply( dsChunkLs, function(dsi){ dsi$Collar[1] })
	#dsia <- dsChunksL[[4]]	
	#dsia <- dsChunksL[[47]]	
	dsChunksLClean <- do.call( rbind, lapply(dsChunksL, function(dsia){
						collar <- dsia$Collar[1] 
						dsi <- dsia[is.finite(dsia$CO2_dry),]
						timeSec <- as.numeric(dsi$TIMESTAMP) - as.numeric( dsi$TIMESTAMP[1] )
						#if( collar == indexNA || nrow(dsi) < minNRec || max(timeSec) < minTime || var(dsi$CO2_dry)==0  ){
						if( collar == indexNA || nrow(dsi) < minNRec || max(timeSec) < minTime || fIsBadChunk(dsi)){
							return( NULL )
						}else return(dsi)
					}))
	### Argument \code{ds} with between-Chunk rows omitted and an additional integer column \code{iChunk} that designates the chunk number.
}
