\name{calcClosedChamberFluxForChunks}
\alias{calcClosedChamberFluxForChunks}
\title{calcClosedChamberFluxForChunks}
\description{apply \code{\link{calcClosedChamberFlux}} for each chunk in data}
\usage{calcClosedChamberFluxForChunks(ds, colChunk = "iChunk", 
    ..., volume, volumesByChunk = volume, isVerbose = TRUE)}
\arguments{
  \item{ds}{tibble or data.frame}
  \item{colChunk}{column name (scalar string) holding a factor unique to each chunk}
  \item{\dots}{
}
  \item{volume}{volume inside the chamber im [m3]}
  \item{volumesByChunk}{Allows for adjusting the chamber volume across different chunks (argument \code{volumne} in \code{\link{calcClosedChamberFlux}})
numeric vector of length of contiguous chunks in ds}
  \item{isVerbose}{set to FALSE to avoid messages}
}



\author{Thomas Wutzler, Oscar Perez Priego}





