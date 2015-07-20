## This function loads all the require packages for functions in this R Package
## Registration. 

## If the package does not exist, it installs them from CRAN Mirror IA (USA). 
## It returns a lit of all packages loaded and a list of all packages to be 

## loaded in case of parallel computation

loadPackages <- function(CRANMirror=83){
  ## Choose USA (IA) as the CRAN mirror
  chooseCRANmirror(graphics=F, ind=83)
  
  Packages <- c(
    'boot',
    'car',
    'cluster',
    'clusterSim',
    'clValid',
    'doParallel',  ## For parallel execution with foreach
    'fda',
    'fdakma',
    'fdasrvf',
    'foreach',
    'ggplot2',
    'grid',
    'gridExtra',
    'gtools',
    'lattice',
    'mclust',
    'plyr',
    'reshape',
    'reshape2'
  )

  ## Requiring packages and installing them if something doesnt exist
  for(Package in Packages){
    if(require(package=Package, character.only=T) == F){
      print(paste('Installing', Package))
      try(install.packages(Package, dependencies = TRUE))
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
