# this is a template to start up an package
# replace all occurrences of myPackageId with the package-Identifier
# * DESCRIPTION
# * tests/doRUnit.R
# * inst/genData/myPackageId-package.Rd
# The name of directory must also be the package-Identifier

# For creating the package, you may consider using the twDev package.
# It is available from 
# https://www.bgc-jena.mpg.de/bgc-mdi/index.php/Intra/ComputingCodeListPackages
# The workspace should be the directory of the DESCRIPTION file.

library(twDev)
?twDev
loadPkg()

genRd()
runRCheck()
svnCommit("descriptive comment of your changes")