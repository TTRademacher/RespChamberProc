projectDir = "inst/paper14"
library(stringr)
library(ggplot2)
library(plyr)


genSMANIEControl <- function(){
	#data(collarCodes)   # merge on conc takes too long
	rawDataPaths="M:/data/SMANIE/Fluxes/raw"
	campaignDirs = list.files( rawDataPaths, pattern="SMANIP_\\d", include.dirs=TRUE)
	campaigns <- as.integer(str_extract(campaignDirs, "\\d"))
	#iCamp <- 1
	SMANIEConc <- do.call( rbind, lapply( seq_along(campaignDirs), function(iCamp){
				cat(",",iCamp)
				fNames <- dir( file.path(rawDataPaths,campaignDirs[iCamp]), pattern="dat$", full.names=TRUE )
				#fName <- fNames[1]
				dsd <- do.call(rbind, lapply( fNames, function(fName){
							ds <- readDat(fName)
							ds$CO2_dry <- corrConcDilution(ds, "CO2_LI840","H2O_LI840")
							dsc <- subsetContiguous(ds, )
							dsc$iChunk <- dsc$iChunk + 100*dsc$Chamber
							dsc
						}))
				cbind(campaign=iCamp, dsd )
			}))
	# TAKE care overwrites
	save( SMANIEconc, file=file.path(projectDir,"SMANIEconc.Rd"))
}

.plotCampaignConcSeries <- function(){
	load(file.path(projectDir,"SMANIEconc.Rd"))
	ds <- SMANIEconc
	#iCamp <- 1
	#dss <- subset(ds, campaign==1 & Chamber==1)
	.tmp.f <- function(){
		d_ply(ds, .(campaign, Chamber), function(dss){
					id <- paste(dss$campaign[1],dss$Chamber[1],sep="_")
					cat(",",id)
					#dss <- dsc
					p1 <- ggplot( dss, aes(x=TIMESTAMP, y=CO2_dry) ) + geom_point() + facet_wrap( ~ iChunk, scales="free")
					#p1
					pdf(width=8.3,height=11.7,file=file.path(projectDir,paste0("SMANIEconc_",id,".pdf")))
					print(p1)
					dev.off()
				})
	}
	id <- factor( c(ds$campaign, ds$iChunk) )
	N <- length(levels(id))
	levels(id) <- 1:N
	ds <- cbind(id64= factor(floor(as.numeric(id)/64)), ds)
	pdf(width=8.3,height=11.7,file=file.path(projectDir,paste0("SMANIEconc.pdf")))
	d_ply(ds, .(id64), function(dss){
				id <- paste(dss$campaign[1],dss$Chamber[1],sep="_")
				cat(",",id)
				#dss <- dsc
				p1 <- ggplot( dss, aes(x=TIMESTAMP, y=CO2_dry) ) + geom_point() + facet_wrap( ~ iChunk, scales="free")
				print(p1)
				#p1
			})
	dev.off()
	
}


