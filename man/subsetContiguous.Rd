\name{subsetContiguous}
\alias{subsetContiguous}
\title{subsetContiguous}
\description{Get contiguous subsets }
\usage{subsetContiguous(ds, colTime = "TIMESTAMP", colIndex = "Collar", 
    gapLength = 12, minNRec = 20, minTime = 60, indexNA = 0, 
    fIsBadChunk = function(dsi) FALSE)}
\arguments{
  \item{ds}{data.frame of measurements}
  \item{colTime}{column name that of time (POSIXct)}
  \item{colIndex}{column name of index variable (factor or integer)}
  \item{gapLength}{minimal length of a gap between subsets (seconds)}
  \item{minNRec}{minimum number of records within one contiguous subset}
  \item{minTime}{minimum length of time that a contiguous subsets covers}
  \item{indexNA}{value of the index column, that signifies records not to use}
  \item{fIsBadChunk}{additional function taking and subset and returning a boolean value whether its a chunk to be omitted}
}
\details{The time series in logger data consists of several chunks of concentration measurments.
In order to determine these chunks, either a change in an index variable (input by between the
measurements) or a gap in time is used.

Between the actural series of measurements, the logger may record sparse data.
These chunks are indicated by value \code{indexNA} in the index column or
by shortness of the series. Only chunks with at least \code{minNRec} records and at least 
\code{minTime} seconds are reported. Others are neglected.}
\value{Argument \code{ds} with between-Chunk rows omitted and an additional integer column \code{iChunk} that designates the chunk number.}

\author{Thomas Wutzler, Oscar Perez Priego}





