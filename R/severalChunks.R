calcClosedChamberFluxForChunks <- function(
	### apply \code{\link{calcClosedChamberFlux}} for each chunk in data
	ds						##<< tibble or data.frame
	,colChunk = "iChunk"	##<< column name (scalar string) holding a factor unique to each chunk
	,...
	,volume					##<< volume inside the chamber im [m3]
	,volumesByChunk=volume	##<< numeric vector of length of contiguous chunks in ds
		### Allows for adjusting the chamber volume across different chunks (argument \code{volumne} in \code{\link{calcClosedChamberFlux}})
	,isVerbose=TRUE			##<< set to FALSE to avoid messages
){
	uniqueChunks <- unique(dsChunked[[colChunk]])
	nChunk <- length(uniqueChunks)
	if( length(volumesByChunk) == 1){
		message("Using the same chamber volume for all chunks.")
		volumesByChunk <- rep(volumesByChunk, nChunk)
	} 
	if( length(volumesByChunk) != nChunk ) stop("argument volumesByChunk must be of same length as chunks in the data (",nChunk,") but was ", length(volumesByChunk))
	names(volumesByChunk) <- uniqueChunks	
	ans <- ds %>% group_by_(.dots=as.symbol(colChunk) ) %>% #summarize(nRec=n())
			do_(~{
						iChunk <- as.character(.[[colChunk]][1])
						vol = volumesByChunk[iChunk]
						calcClosedChamberFlux(.,...,volume=vol)
					})
	ans
}


