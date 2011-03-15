\name{tran.polar}
\alias{tran.polar}
\alias{polar2cart}

\title{
  Diffusive Transport in polar (r, theta) coordinates.
}

\description{
  Estimates the transport term (i.e. the rate of change of a concentration
  due to diffusion) in a polar (r, theta)  coordinate system
}

\usage{
tran.polar (C, C.r.up = NULL, C.r.down = NULL, 
            C.theta.up = NULL, C.theta.down = NULL, 
            flux.r.up = NULL, flux.r.down = NULL, 
            flux.theta.up = NULL, flux.theta.down = NULL, 
            cyclicBnd = NULL, D.r = 1, D.theta = D.r, 
            r = NULL, theta = NULL, full.output = FALSE)

polar2cart (out, r, theta, x = NULL, y = NULL)

}

\arguments{
  \item{C }{concentration, expressed per unit volume, defined at the centre
    of each grid cell; Nr*Nteta matrix [M/L3].
  }
  \item{C.r.up }{concentration at upstream boundary in r(x)-direction;
    vector of length Nteta [M/L3].
  }
  \item{C.r.down }{concentration at downstream boundary in r(x)-direction;
    vector of length Nteta [M/L3].
  }
  \item{C.theta.up }{concentration at upstream boundary in theta-direction;
    vector of length Nr [M/L3].
  }
  \item{C.theta.down }{concentration at downstream boundary in theta-direction;
    vector of length Nr [M/L3].
  }
  \item{flux.r.up }{flux across the upstream boundary in r-direction,
    positive = INTO model domain; vector of length Ntheta [M/L2/T].
  }
  \item{flux.r.down }{flux across the downstream boundary in r-direction,
    positive = OUT of model domain; vector of length Ntheta [M/L2/T].
  }
  \item{flux.theta.up }{flux across the upstream boundary in theta-direction,
    positive = INTO model domain; vector of length Nr [M/L2/T].
  }
  \item{flux.theta.down }{flux across the downstream boundary in theta-direction,
    positive = OUT of model domain; vector of length Nr [M/L2/T].
  }
  \item{cyclicBnd }{If not \code{NULL}, the direction in which a cyclic 
    boundary is defined, i.e. \code{cyclicBnd = 1} for the \code{r} direction, 
    \code{cyclicBnd = 2} for the \code{theta} direction and 
    \code{cyclicBnd = c(1,2)} for both the \code{r} and \code{theta} direction.
  }
  \item{D.r }{diffusion coefficient in r-direction, defined on grid cell
    interfaces. One value, a vector of length (Nr+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nr+1)* Nteta matrix [L2/T].
  }
  \item{D.theta }{diffusion coefficient in theta-direction, defined on grid cell
    interfaces. One value, a vector of length (Ntheta+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nr*(Ntheta+1) matrix [L2/T].
  }
  \item{r }{position of adjacent cell interfaces in the r-direction. 
    A vector of length Nr+1 [L].
  }
  \item{theta }{position of adjacent cell interfaces in the theta-direction.
    A vector of length Ntheta+1 [L].  Theta should be within [0,2 pi]
  }
  \item{full.output }{logical flag enabling a full return of the output
    (default = \code{FALSE}; \code{TRUE} slows down execution by 20 percent).
  }
  \item{out }{output as returned by \code{tran.polar}, and which is to be
    mapped from polar to cartesian coordinates
  }
  \item{x }{The cartesian x-coordinates to whicht the polar coordinates are
    to be mapped
  }
  \item{y }{The cartesian y-coordinates to whicht the polar coordinates are
    to be mapped
  }
}
\value{
  a list containing:
  \item{dC }{the rate of change of the concentration C due to transport,
    defined in the centre of each grid cell, a Nr*Nteta matrix. [M/L3/T].
  }
  \item{C.r.up }{concentration at the upstream interface in r-direction.
    A vector of length Nteta [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{C.r.down }{concentration at the downstream interface in r-direction.
    A vector of length Nteta [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{C.theta.up }{concentration at the the upstream interface in theta-direction.
    A vector of length Nr [M/L3].  Only when \code{full.output = TRUE}.
  }
  \item{C.theta.down }{concentration at the downstream interface in theta-direction.
    A vector of length Nr [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{r.flux }{flux across the interfaces in x-direction of the grid cells.
    A (Nr+1)*Nteta matrix  [M/L2/T]. Only when \code{full.output = TRUE}.
  }
  \item{theta.flux }{flux across the interfaces in y-direction of the grid cells.
    A Nr*(Nteta+1) matrix [M/L2/T].  Only when \code{full.output = TRUE}.
  }
  \item{flux.r.up }{flux across the upstream boundary in r-direction,
    positive = INTO model domain. A vector of length Nteta [M/L2/T].
  }
  \item{flux.r.down }{flux across the downstream boundary in r-direction,
    positive = OUT of model domain. A vector of length Nteta [M/L2/T].
  }
  \item{flux.theta.up }{flux across the upstream boundary in theta-direction,
    positive = INTO model domain. A vector of length Nr [M/L2/T].
  }
  \item{flux.theta.down }{flux across the downstream boundary in theta-direction,
   positive = OUT of model domain. A vector of length Nr [M/L2/T].
  }

}

\examples{
## =============================================================================
## Testing the functions
## =============================================================================
# Parameters
F        <- 100             # input flux [micromol cm-2 yr-1]
D        <- 400             # mixing coefficient [cm2 yr-1]

# Grid definition
r.N   <- 4     # number of cells in r-direction
theta.N <- 6   # number of cells in theta-direction
r.L <- 8       # domain size r-direction [cm]
r      <- seq(0, r.L,len = r.N+1)      # cell size r-direction [cm]
theta  <- seq(0, 2*pi,len = theta.N+1) # theta-direction - theta: from 0, 2pi
 
# Intial conditions 
C <- matrix(nrow = r.N, ncol = theta.N, data = 0)

# Boundary conditions: fixed concentration  
C.r.up       <- rep(1, times = theta.N)
C.r.down     <- rep(0, times = theta.N)
C.theta.up   <- rep(1, times = r.N)
C.theta.down <- rep(0, times = r.N)

# Concentration boundary conditions
tran.polar(C = C, D.r = D, D.theta = D, 
  r = r, theta = theta,
  C.r.up = C.r.up, C.r.down = C.r.down, 
  C.theta.up = C.theta.up, C.theta.down = C.theta.down)

# Flux boundary conditions
flux.r.up <- rep(200, times = theta.N)
flux.r.down <- rep(-200, times = theta.N)
flux.theta.up <- rep(200, times = r.N)
flux.theta.down <- rep(-200, times = r.N)

tran.polar(C = C, D.r = D, r = r, theta = theta,
  flux.r.up = flux.r.up, flux.r.down = flux.r.down,
  flux.theta.up = flux.theta.up, flux.theta.down = flux.theta.down,
  full.output = TRUE)


## =============================================================================
## A model with diffusion and first-order consumption
## =============================================================================
N     <- 50          # number of grid cells
XX    <- 4           # total size
rr    <- 0.005       # consumption rate
ini   <- 1           # initial value at x=0
D     <- 400

r     <- seq (2, 4, len = N+1)
theta   <- seq(0, 2*pi, len = N+1)
theta.m <- 0.5*(theta[-1]+theta[-(N+1)])

# The model equations

Diffpolar <- function (t, y, parms)  {
  CONC  <- matrix(nrow = N, ncol = N, data = y)
  tran  <- tran.polar(CONC, D.r = D, D.theta = D, r = r, theta = theta,
        C.r.up = 0, C.r.down = 1*sin(5*theta.m), 
        cyclicBnd = 2, full.output=TRUE )
  dCONC <- tran$dC  - rr * CONC
  return (list(dCONC))
}

# solve to steady-state; cyclicBnd = 2, because of C.theta.up, C.theta.down
out <- steady.2D (y = rep(0, N*N), func = Diffpolar, parms = NULL,
                  dim = c(N, N), lrw = 1e6, cyclicBnd = 2)

image(out)

cart <- polar2cart(out, r = r, theta = theta, 
                        x = seq(-4, 4, len = 100), 
                        y = seq(-4, 4, len = 100))
image(cart)
          
}

\references{
  Soetaert and Herman, 2009. a practical guide to ecological modelling -
  using R as a simulation platform. Springer
}

\details{
  \code{tran.polar} performs (simplified) transport in polar coordinates

  The \bold{boundary conditions} are either
  \itemize{
    \item (1) zero gradient
    \item (2) fixed concentration
    \item (3) fixed flux
    \item (4) cyclic boundary
  }

  This is also the order of priority. The cyclic boundary overrules the other.
  If fixed concentration, fixed flux, and cyclicBnd are \code{NULL} then
  the boundary is zero-gradient

  A cyclic boundary condition has concentration and flux at upstream and 
  downstream boundary the same.
  
  \code{polar2cart} maps the polar coordinates to cartesian coordinates
  
  If \code{x} and \code{y} is not provided, then it will create an (x,y)
  grid based on \code{r} : \code{seq(-maxr, maxr, length.out=Nr)}, where
  \code{maxr} is the maximum value of \code{r}, and \code{Nr} is the number
  of elements in \code{r}.
}  

\keyword{utilities}

\seealso{
  \code{\link{tran.cylindrical}}, \code{\link{tran.spherical}} 
   for a discretisation of 3-D transport equations in cylindrical and 
   spherical coordinates
  
  \code{\link{tran.1D}},   \code{\link{tran.2D}},  \code{\link{tran.3D}}
}
