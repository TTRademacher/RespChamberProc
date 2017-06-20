#require(testthat)
context("readDat")

test_that("x81 example",{
			fName <- system.file("genData/Flux2_140929_1700.81x", package = "RespChamberProc")
			if( nzchar(fName) ){
				ds <- read81x(fName)
				expect_true( max(ds$iChunk) > 1 )
				expect_true( table(ds$iChunk)[1] > 1 )
				#plot( CO2 ~ Date, ds )
				#plot( CO2 ~ Date, ds[ds$iChunk==2,] )
			}
		})

test_that("x81 annotation",{
			fName <- system.file("genData/testVaryingColNumber.81x", package = "RespChamberProc")
			if( nzchar(fName) ){
				ds <- read81x(fName)
				expect_true( max(ds$iChunk) > 1 )
				expect_true( table(ds$iChunk)[1] > 1 )
				expect_true( "some Comment" %in% ds$Annotation )
				#plot( CO2 ~ Date, ds )
				#plot( CO2 ~ Date, ds[ds$iChunk==9,] )
			}
		})
