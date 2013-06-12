\name{setup.prop.1D}
\alias{setup.prop.1D}
\alias{plot.prop.1D}
\title{
  Attaches a Property to a One-Dimensional Grid
}

\description{
  This routine calculates the value of a given property at the middle of
  the grid cells (\code{mid}) and at the interfaces of the grid cells
  (\code{int}).

  Two possibilities are available: either specifying a mathematical function
  (\code{func}) that describes the spatial dependency of the property,
  or obtaining the property from interpolation of a data series
  (via the input of the data matrix \code{xy}).

  For example, in a sediment model, \code{setup.prop.1D} can be used to
  specify the porosity, the mixing intensity or some other parameter over
  the one-dimensional grid. Similarly, in a vertical water column model, \code{setup.prop.1D} can be
  used to specify the sinking velocity of particles or other model parameters
  changing with water depth.
}

\usage{
setup.prop.1D(func = NULL, value = NULL, xy = NULL,
              interpolate = "spline", grid, ...)

\method{plot}{prop.1D}(x, grid, xyswap = FALSE, \dots)

}

\arguments{
  \item{func }{function that describes the spatial dependency. For example, one can use the functions provided in \code{\link{p.exp}} 
  }
  \item{value }{constant value given to the property (no spatial dependency)
  }
  \item{xy }{a two-column data matrix where the first column (\code{x}) provides
    the position, and the second column (\code{y}) provides the values that
    needs interpolation over the grid
  }
  \item{interpolate }{specifies how the interpolation should be done, one
    of "spline" or "linear"; only used when \code{xy} is present
  }
  \item{grid }{list specifying the 1D grid characteristics, see
    \code{\link{setup.grid.1D}} for details on the structure of this list
  }
  \item{x }{the object of class \code{prop.1D} that needs plotting
  }
  \item{xyswap }{if \code{TRUE}, then x- and y-values are swapped and
    the y-axis is oriented from top to bottom. Useful for drawing
    vertical depth profiles
  }
  \item{... }{additional arguments that are passed on to \code{func} or
    to the S3 method
  }
}

\value{
  A list of type \code{prop.1D} containing:
  \item{mid }{property value in the middle of the grid cells; vector of
    length N (where N is the number of grid cells)
  }
  \item{int }{property value at the interface of the grid cells;
    vector of length N+1
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>,
  Filip Meysman <filip.meysman@nioz.nl>
}

\examples{

# Construction of the 1D grid 

grid <- setup.grid.1D(x.up = 0, L = 10, N = 10)

# Porosity profile via function specification

P.prof <- setup.prop.1D(func = p.exp, grid = grid, y.0 = 0.9,
                        y.inf = 0.5, x.att = 3)

# Porosity profile via data series interpolation

P.data <- matrix(ncol = 2, data = c(0,3,6,10,0.9,0.65,0.55,0.5))
P.spline <- setup.prop.1D(xy = P.data, grid = grid)
P.linear <- setup.prop.1D(xy = P.data, grid = grid, interpolate = "linear")

# Plot different profiles 

plot(P.prof, grid = grid, type = "l",
     main = "setup.prop, function evaluation")
points(P.data, cex = 1.5, pch = 16) 
lines(grid$x.int, P.spline$int, lty = "dashed")
lines(grid$x.int, P.linear$int, lty = "dotdash")
}

\details{
  There are two options to carry out the data interpolation:
  \itemize{
     \item "spline" gives a smooth profile, but sometimes generates
       strange profiles - always check the result!
     \item "linear" gives a segmented profile
  }
}

\keyword{utilities}

\seealso{
  \code{\link{tran.1D}}, for a discretisation of the general transport equation in 1-D

  \code{\link{setup.grid.1D}}, the creation of grids in 1-D

  \code{\link{setup.prop.2D}} for defining properties on 2-D grids.
}
