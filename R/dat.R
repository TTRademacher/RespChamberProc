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
){
	##details<< 
	## Assumes that there
	setClass("myDate")
	setAs("character","myDate", function(from) as.POSIXct(from, format=formatTS, tz=tz) )
	fileInfo <- readLines(fName, n=nRowsFileInfo )
	colInfo <- read.table(fName, header=TRUE, skip=nRowsFileInfo, nrows=max(1,nRowsColInfo), sep=sep)
	colClasses[colsTimeStamp] <- "myDate"
	rawData <- read.table(fName, header=FALSE, skip=nRowsFileInfo+1+nRowsColInfo, sep=sep, ...
					,colClasses=colClasses)
	colnames(rawData) <- colnames(colInfo)	
	#plot( CO2_Avg ~ TIMESTAMP, data=rawData )
	attr(rawData,"fileInfo") <- fileInfo
	attr(rawData,"colInfo") <- colInfo
	rawData
}
attr(readDat,"ex") <- function(){
	ds <- readDat("tmp/chamberLoggerEx1_short.dat")
}