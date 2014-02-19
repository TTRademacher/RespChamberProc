# For creating the package, you may consider using the twDev package.
# It is available from 
# https://www.bgc-jena.mpg.de/bgi/index.php/Intra/ComputingCodeListPackages
# The directory must be names according to the package id
# The workspace should be the directory of the DESCRIPTION file.

library(twDev)
?twDev
loadPkg()

genRd()
runRCheck()
hgCommit("descriptive comment of your changes")
hgPush()
