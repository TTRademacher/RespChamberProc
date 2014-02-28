ds <- readDat("tmp/MANIP_Ch1_0.dat", tz="CET")

summary(ds)
plot( CO2_LI840 ~ TIMESTAMP, ds)


