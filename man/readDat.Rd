\name{readDat}
\alias{readDat}
\title{readDat}
\description{Read a data loggger .dat file into a data-frame		}
\usage{readDat(fName, nRowsFileInfo = 1, nRowsColInfo = 2, sep = ",", 
    ..., colClasses = rep(NA, ncol(colInfo)), colsTimeStamp = 1, 
    formatTS = "\%Y-\%m-\%d \%H:\%M:\%S", tz = "UTC", na.strings = c("", 
        "NA", "NAN", "\\"NAN\\""))}
\arguments{
  \item{fName}{scalar string: file name}
  \item{nRowsFileInfo}{integer scalar: number of lines before column information}
  \item{nRowsColInfo}{integer vector: number of lines with column description}
  \item{sep}{column separator}
  \item{\dots}{further arguments to}
  \item{colClasses}{see \code{link{read.table}}}
  \item{colsTimeStamp}{integer vector: colums with time stamp column (will be set to POSIXct}
  \item{formatTS}{format of the timestamp columns, see \code{\link{strptime}}, e.g.}
  \item{tz}{specify a time zone when converting to POSIXct, default: current local e.g CET, UTC}
  \item{na.strings}{
}
}
\details{Assumes that there}


\author{Thomas Wutzler, Oscar Perez Priego}




\examples{
fName <- system.file("genData/chamberLoggerEx1_short.dat", package = "RespChamberProc")
if( nzchar(fName) ){
	ds <- readDat(fName)
}
}
