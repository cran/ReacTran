\name{setup.grid.1D}
\alias{setup.grid.1D}
\alias{plot.grid.1D}
\title{
  Creates a One-Dimensional Finite Difference Grid
}
\description{
  Subdivides the one-dimensional model domain into one or more zones that
  are each sub-divided into grid cells. The resulting grid structure can be
  used in the other \code{ReacTran} functions.

  The grid structure is characterized by the position of the middle of
  the grid cells (\code{x.mid}) and the position of the interfaces between
  grid cells (\code{x.int}).

  Distances are calculated between the interfaces (\code{dx}), i.e. the
  thickness of the grid cells. An auxiliary set of distances (\code{dx.aux})
  is calculated between the points where the concentrations are specified
  (at the center of each grid cell and the two external interfaces).

  A more complex grid consisting of multiple zones can be constructed when
  specifying the endpoints of ech zone (\code{x.down}), the interval length
  (\code{L}), and the number of layers in each zone (\code{N}) as vectors.
  In each zone, one can control the grid resolution near the upstream and
  downstream boundary.
  
  The grid resolution at the upstream interface changes according to the
  power law relation \code{dx[i+1] = min(max.dx.1,p.dx.1*dx[i])},
  where \code{p.dx.1} determines the rate of increase and \code{max.dx.1}
  puts an upper limit on the grid cell size.

  A similar formula controls the resolution at the downstream interface.
  This allows refinement of the grid near the interfaces.
  
  If only \code{x.up, N} and \code{dx.1} are specified, then the grid size
  is taken constant = \code{dx.1} (and \code{L=N*dx.1})
}

\usage{
setup.grid.1D(x.up = 0, x.down = NULL, L = NULL, N = NULL, dx.1 = NULL,
              p.dx.1 = rep(1, length(L)), max.dx.1 = L, dx.N = NULL,
              p.dx.N = rep(1, length(L)), max.dx.N = L)

\method{plot}{grid.1D}(x, \dots)

}

\arguments{
  \item{x.up }{position of the upstream interface; one value [L]
  }
  \item{x.down }{position of the endpoint of each zone; one value when the
    model domain covers only one zone (\code{x.down} = position of downstream
    interface), or a vector of length M when the model domain is divided into
    M zones (\code{x.down[M]} = position of downstream interface) [L]
  }
  \item{L }{thickness of zones; one value (model domain = one zone) or a
    vector of length M (model domain = M zones) [L]
  }
  \item{N }{number of grid cells within a zone; one value or a vector of
    length M [-]
  }
  \item{dx.1 }{size of the first grid cell in a zone; one value or a vector
    of length M [L]
  }
  \item{p.dx.1 }{power factor controlling the increase in grid cell size near
    the upstream boundary; one value or a vector of length M. The default
    value is 1 (constant grid cell size) [-]
  }
  \item{max.dx.1 }{maximum grid cell size in the upstream half of the zone;
    one value or a vector of length M [L]
  }
  \item{dx.N }{size of the last grid cell in a zone; one value or a vector
    of length M [L]
  }
  \item{p.dx.N }{power factor controlling the increase in grid cell size near the
    downstream boundary; one value or a vector of length M. The default value
    is 1 (constant grid cell size) [-]
  }
  \item{max.dx.N }{maximum grid cell size in the downstream half of the zone;
    one value or a vector of length M [L]
  }
  
  \item{x }{the object of class \code{grid.1D} that needs plotting
  }
  \item{...}{additional arguments passed to the function \link{plot}
  }
}

\value{
  a list of type \code{grid.1D} containing:

  \item{N }{the total number of grid cells [-]
  }
  \item{x.up }{position of the upstream interface; one value [L]
  }
  \item{x.down }{position of the downstream interface; one value [L]
  }
  \item{x.mid }{position of the middle of the grid cells;
    vector of length \code{N} [L]
  }
  \item{x.int }{position of the interfaces of the grid cells;
    vector of length \code{N+1} [L]
  }
  \item{dx }{distance between adjacent cell interfaces (thickness of grid
    cells); vector of length \code{N} [L]
  }
  \item{dx.aux }{auxiliary vector containing the distance between adjacent
    cell centers; at the upper and lower boundary calculated as
    (\code{x[1]-x.up}) and (\code{x.down-x[N]}) respectively; vector of
    length \code{N+1} [L]
  }
}

\author{
  Filip Meysman <filip.meysman@nioz.nl>,
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
# one zone, constant resolution
(GR <- setup.grid.1D(x.up = 0, L = 10, N = 10))
(GR <- setup.grid.1D(x.up = 0, L = 10, dx.1 = 1))
(GR <- setup.grid.1D(x.up = 0, L = 10, dx.N = 1))
plot(GR)

# one zone, constant resolution, origin not zero
(GR <- setup.grid.1D(x.up = 5, x.down = 10, N = 10))
plot(GR)

# one zone, variable resolution
(GR <- setup.grid.1D(x.up = 0, L = 10, dx.1 = 1, p.dx.1 = 1.2))
(GR <- setup.grid.1D(x.up = 0, L = 10, dx.N = 1, p.dx.N = 1.2))
plot(GR)

# one zone, variable resolution, imposed number of layers
(GR <- setup.grid.1D(x.up = 0, L = 10, N = 6, dx.1 = 1, p.dx.1 = 1.2))
(GR <- setup.grid.1D(x.up = 0, L = 10, N = 6, dx.N = 1, p.dx.N = 1.2))
plot(GR)

# one zone, higher resolution near upstream and downstream interfaces
(GR <- setup.grid.1D(x.up = 0, x.down = 10, dx.1 = 0.1,
                     p.dx.1 = 1.2, dx.N = 0.1, p.dx.N = 1.2))
plot(GR)

# one zone, higher resolution near upstream and downstream interfaces
# imposed number of layers 
(GR <- setup.grid.1D(x.up = 0, x.down = 10, N = 20, 
                     dx.1 = 0.1, p.dx.1 = 1.2, 
                     dx.N = 0.1, p.dx.N = 1.2))
plot(GR)

# two zones, higher resolution near the upstream
# and downstream interface
(GR<-setup.grid.1D(x.up = 0, L = c(5, 5), 
         dx.1 = c(0.2, 0.2), p.dx.1 = c(1.1, 1.1), 
         dx.N = c(0.2, 0.2), p.dx.N = c(1.1, 1.1)))
plot(GR)

# two zones, higher resolution near the upstream
# and downstream interface
# the number of grid cells in each zone is imposed via N
(GR <- setup.grid.1D(x.up = 0, L = c(5, 5), N = c(20, 10),
         dx.1 = c(0.2, 0.2), p.dx.1 = c(1.1, 1.1), 
         dx.N = c(0.2, 0.2), p.dx.N = c(1.1, 1.1)))
plot(GR)
}
\keyword{utilities}

\seealso{
  \code{\link{tran.1D}}, for a discretisation of the general transport equation in 1-D

  \code{\link{setup.grid.2D}} for the creation of grids in 2-D

  \code{\link{setup.prop.1D}}, for defining properties on the 1-D grid.
}
