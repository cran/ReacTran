\name{ReacTran-package}
\alias{ReacTran-package}
\alias{ReacTran}
\docType{package}
\title{Reactive transport modelling in 1D, 2D and 3D}
\description{
  
  R-package ReacTran contains routines that enable the development of
  reactive transport models in aquatic systems (rivers, lakes, oceans), 
  porous media (floc aggregates, sediments,...) and
  even idealized organisms (spherical cells, cylindrical worms,...).

  The geometry of the model domain is either one-dimensional, two-dimensional
  or three-dimensional.
  
  The package contains:
  \itemize{
    \item Functions to setup a finite-difference grid (1D or 2D)
    \item Functions to attach parameters and properties to this grid (1D or 2D)
    \item Functions to calculate the advective-diffusive transport term over
      the grid (1D, 2D, 3D)
    \item Utility functions
  }
  
  
}
\details{
  \tabular{ll}{
    Package: \tab ReacTran\cr
    Type: \tab Package\cr
    Version: \tab 1.3\cr
    Date: \tab 2010-11-01\cr
    License: \tab GNU Public License 2 or above\cr
  }
}
\author{
  Karline Soetaert (Maintainer)
  
  Filip Meysman
}
\examples{
\dontrun{

## show examples (see respective help pages for details)
## 1-dimensional transport
example(tran.1D)
example(tran.volume.1D)

## 2-dimensional transport
example(tran.2D)
example(tran.polar)

## 3-dimensional transport
example(tran.3D)
example(tran.cylindrical)
example(tran.spherical)

## open the directory with documents
browseURL(paste(system.file(package="ReacTran"), "/doc", sep=""))

## open the directory with fortran codes of the transport functions
browseURL(paste(system.file(package="ReacTran"), "/doc/fortran", sep=""))

## show package vignette with how to use ReacTran and how to solve PDEs
## + source code of the vignettes
vignette("ReacTran")
vignette("PDE")
edit(vignette("ReacTran"))

## a directory with fortran implementations of the transport
browseURL(paste(system.file(package="ReacTran"), "/doc/fortran", sep=""))
}
}
\keyword{ package }
\seealso{
  Functions \link{ode.1D}, \link{ode.2D}, \link{ode.3D} from package \code{deSolve} 
  to integrate the reactive-transport model
  
  Functions \link{steady.1D}, \link{steady.2D}, \link{steady.3D} from 
  package \code{rootSolve} to find the steady-state solution of the 
  reactive-transport model

  \code{\link{tran.1D}},   \code{\link{tran.2D}},  \code{\link{tran.3D}} for 
  a discretisation of the general transport equations

  \code{\link{tran.volume.1D}} for discretisation of the 1-D transport equations using finite volumes
  
  \code{\link{tran.cylindrical}}, \code{\link{tran.spherical}} 
   for a discretisation of 3-D transport equations in cylindrical and 
   spherical coordinates

  \code{\link{tran.polar}}, 
   for a discretisation of 2-D transport equations in polar coordinates  

  \code{\link{setup.grid.1D}},   \code{\link{setup.grid.2D}} for the creation of grids in 1-D and in 2-D

  \code{\link{setup.prop.1D}},   \code{\link{setup.prop.2D}} for defining properties on these grids.
}
