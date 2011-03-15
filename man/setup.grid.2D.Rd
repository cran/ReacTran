\name{setup.grid.2D}
\alias{setup.grid.2D}
\title{
  Creates a Finite Difference Grid over a Two-Dimensional Rectangular Domain
}
\description{
  Creates a finite difference grid over a rectangular two-dimensional model
  domain starting from two separate one-dimensional grids (as created by
  \code{\link{setup.grid.1D}}).
}

\usage{
setup.grid.2D(x.grid = NULL, y.grid = NULL)
}

\arguments{
  \item{x.grid }{list containing the one-dimensional grid in the vertical
    direction - see \code{\link{setup.grid.1D}} for the structure of the list
  }
  \item{y.grid }{list containing the one-dimensional grid in the horizontal
    direction - see \code{\link{setup.grid.1D}} for the structure of the list
  }
}

\value{
  a list of type \code{grid.2D} containing:
  \item{x.up }{position of the upstream interface in x-direction
    (i.e. if x is vertical, the upper boundary); one value
  }
  \item{x.down }{position of the downstream interface in x-direction
    (i.e. if x is vertical, the lower boundary); one value
  }
  \item{x.mid }{position of the middle of the grid cells in x-direction;
    vector of length \code{x.N}
  }
  \item{x.int }{position of the interfaces of the grid cells in x-direction;
    vector of length \code{x.N+1}
  }
  \item{dx }{distance between adjacent cell interfaces in x-direction
    (thickness of grid cells); vector of length \code{x.N}
  }
  \item{dx.aux }{auxiliary vector containing the distance between adjacent
    cell centers; at the upstream and downstream boundary calculated as
    (\code{x[1]-x.up}) and (\code{x.down-x[x.N]}) respectively;
    vector of length \code{x.N+1}
  }
  \item{x.N }{total number of grid cells in the x direction; one value
  }
  \item{y.left }{position of the upstream interface in y-direction
    (i.e. if y us the horizontal, the left boundary); one value
  }
  \item{y.right }{position of the downstream interface in y-direction
    (i.e. if y us the horizontal, the right boundary); one value
  }
  \item{y.mid }{position of the middle of the grid cells in y-direction;
    vector of length \code{y.N}
  }
  \item{y.int }{position of the interfaces of the grid cells in y-direction;
    vector of length \code{y.N+1}
  }
  \item{dy }{distance between adjacent cell interfaces in y-direction
    (thickness of grid cells); vector of length \code{y.N}
  }
  \item{dy.aux }{auxiliary vector containing the distance between adjacent
    cell centers; at the upstream and downstream boundary calculated as
    (\code{y[1]-y.up}) and (\code{y.down-y[y.N]}) respectively;
    vector of length \code{y.N+1}
  }
  \item{y.N }{total number of grid cells in the y direction;
    one value
    }
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
# test of the setup.grid.2D functionality
x.grid  <- setup.grid.1D(x.up = 0, L = 10, N = 5)
y.grid  <- setup.grid.1D(x.up = 0, L = 20, N = 10)
(grid2D <- setup.grid.2D(x.grid, y.grid))

}
\keyword{utilities}

\seealso{
  \code{\link{tran.2D}},  for a discretisation of the general transport equation in 2-D

  \code{\link{setup.grid.1D}}, for the creation of grids in 1-D

  \code{\link{setup.prop.2D}} for defining properties on the 2-D grid.
}
