ds <- readDat("tmp/chamberLoggerEx1.dat")
ds$Pa <- 101		# 101 kPa = 1010 hPa = 1.01 bar

chamberLoggerEx1s <- head(ds,30)
save( chamberLoggerEx1s, file=file.path("data","chamberLoggerEx1s.RData") )

