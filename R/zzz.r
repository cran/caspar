.packageName <- 'caspar'

.First.lib <- function(lib, pkg)
{
  library.dynam("caspar", pkg, lib)
  require(mvtnorm)
}
