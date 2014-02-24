ds <- readDat("tmp/chamberLoggerEx1.dat")
ds$Pa <- 101		# 101 kPa = 1010 hPa = 1.01 bar
#plot( CO2_Avg ~ TIMESTAMP, ds)
min(ds$TIMESTAMP)
ds1 <- subset( ds, TIMESTAMP <= as.POSIXct("2013-03-12"))
#ds2 <- subset(ds1, as.integer(format(ds1$TIMESTAMP,"%H")) < 10)
ds2 <- ds1[1:100,]
plot( CO2_Avg ~ TIMESTAMP, ds2)

chamberLoggerEx1s <- ds2
save( chamberLoggerEx1s, file=file.path("data","chamberLoggerEx1s.RData") )

