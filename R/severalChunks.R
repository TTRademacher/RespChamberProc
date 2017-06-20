calcClosedChamberFluxForChunks <- function(
	### apply \code{\link{calcClosedChamberFlux}} for each chunk in data
	ds						##<< tibble or data.frame
	,colChunk = "iChunk"	##<< column name (scalar string) holding a factor unique to each chunk
	,...
){
	ans <- ds %>% group_by_(.dots=as.symbol(colChunk) ) %>% #summarize(nRec=n())
			do_(~calcClosedChamberFlux(.,...))
	ans
}


