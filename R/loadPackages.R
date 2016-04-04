## This function loads all the require packages for functions in this R Package
## Registration. 

## If the package does not exist, it installs them from CRAN Mirror IA (USA). 
## It returns a lit of all packages loaded and a list of all packages to be 
## loaded in case of parallel computation

loadPackages <- function(CRANMirror=83){
  ## Choose USA (IA) as the CRAN mirror
  Mirrors <- getCRANmirrors(all = FALSE, local.only = FALSE)
  chooseCRANmirror(graphics = F, ind = which(Mirrors$Name == 'USA (IA)'))
  
  Packages <- c(
    'Biostrings',  ## For sequence comparison
    'boot',
    'car',
    'cluster',
    'clusterSim',
    'clValid',
    'doParallel',  ## For parallel execution with foreach
    'doSNOW',
    'fda',
    'fdakma',
    'fdasrvf',
    'fda.usc',     ## For functional data depth
    'foreach',
    'ggplot2',
    'gridBase',
    'gridExtra',
    'gtools',
    'lattice',
    'MASS',
    'matrixStats', ## For weighted row means
    'mclust',
    'png',
    'plyr',
    'reshape',
    'reshape2',
    'robustX',   ## For multivariate median
    'rpart',
    'seqinr'
  )

  ## Requiring packages and installing them if something doesnt exist
  for(Package in Packages){
    if(require(package=Package, character.only=T) == F){
      print(paste('Installing', Package))
      #try(install.packages(Package, dependencies = TRUE))
    } else{
      print(paste('Loading', Package))
      require(package=Package, character.only=T)
    }
  }

  ## For parallel processing, when passing the list of packages to load
  ## in all the cores. Could be different from Packages
  Packages_Par <- Packages
  return(list(Packages=Packages, Packages_Par=Packages_Par))
}
