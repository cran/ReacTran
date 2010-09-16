\name{setup.prop.2D}
\alias{setup.prop.2D}
\alias{contour.prop.2D}
\title{
  Attaches a Property to a Two-Dimensional Grid
}

\description{
  Calculates the value of a given property at the middle of grid cells
  (\code{mid}) and at the interfaces of the grid cells (\code{int}).

  Two possibilities are available: either specifying a mathematical function
  (\code{func}) that describes the spatial dependency of the property, or
  asssuming a constant value (\code{value}). To allow for anisotropy, the
  spatial dependency can be different in the x and y direction.

  For example, in a sediment model, the routine can be used to specify the
  porosity, the mixing intensity or other parameters over the grid of the
  reactangular sediment domain.
}

\usage{
setup.prop.2D(func = NULL, value = NULL, grid, y.func = func, 
  y.value = value, ...)

\method{contour}{prop.2D}(x, grid, xyswap = FALSE, filled = FALSE, \dots)
}

\arguments{
  \item{func }{function that describes the spatial dependency in the 
    x-direction; defined as \code{func <- function (x,y, ...)}
  }
  \item{value }{constant value given to the property in the x-direction
  }
  \item{grid }{list specifying the 2D grid characteristics, see
    \code{\link{setup.grid.2D}} for details on the structure of this list
  }
  \item{y.func }{function that describes the spatial dependency in the
    y-direction. By default the same as in the x-direction.
    \code{y.func <- function (x, y, ...)}
  }
  \item{y.value }{constant value given to the property in the y-direction.
    By default the same as in the x-direction.
  }
  \item{x }{the object of class \code{prop.2D} that needs plotting
  }
  \item{filled }{if \code{TRUE}, uses \code{filled.contour}, else \code{contour}
  }
  \item{xyswap }{if \code{TRUE}, then x- and y-values are swapped and
    the y-axis is oriented from top to bottom. Useful for drawing
    vertical depth profiles
  }
  \item{... }{additional arguments that are passed on to \code{func} or
    to the method
  }
}

\value{
  A list of type \code{prop.2D} containing:
  \item{x.mid }{property value in the x-direction defined at the middle of the
    grid cells; Nx * Ny matrix (where Nx and Ny = number of cells in x, y
    direction)
  }
  \item{y.mid }{property value in the y-direction at the middle of the grid
    cells; Nx * Ny matrix
  }
  \item{x.int }{property value in the x-direction defined at the
    x-interfaces of the grid cells;
    (Nx+1)*Ny matrix
  }
  \item{y.int }{property value in the y-direction at the y-interfaces of the
    grid cells; Nx*(Ny+1) matrix
  }
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\note{
  For some properties, it does not make sense to use \code{y.func} different 
  to \code{func}. For instance, for volume fractions, AFDW.
  
  For other properties, it may be usefull to have \code{y.func} or 
  \code{y.value} different from \code{func} or \code{value}, for instance
  for velocities, surface areas, ...
}

\examples{
# Inverse quadratic function 
inv.quad <- function(x,y,a=NULL,b=NULL)
return(1/((x-a)^2+(y-b)^2))


# Construction of the 2D grid 
x.grid <- setup.grid.1D(x.up=0,L=10,N=10)
y.grid <- setup.grid.1D(x.up=0,L=10,N=10)
grid2D <- setup.grid.2D(x.grid,y.grid)

# Attaching the inverse quadratic function to the 2D grid 
(twoD<-setup.prop.2D(func=inv.quad,grid=grid2D,a=5,b=5))
contour(log(twoD$x.int))
}

\details{
  \itemize{
    \item {When the property is isotropic, the \code{x.mid} and
      \code{y.mid} values are identical. This is for example the case for
      sediment porosity.}
    \item {When the property is anisotropic, the \code{x.mid} and
      \code{y.mid} values can differ. This can be for example the case for
      the velocity, where in general, the value will differ between the x and
      y direction.}
  }
}
\keyword{utilities}

