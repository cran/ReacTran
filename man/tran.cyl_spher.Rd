\name{tran.cylindrical}
\alias{tran.cylindrical}
\alias{tran.spherical}
\title{
  Diffusive Transport in cylindrical (r, theta, z) and spherical (r, theta, phi)
  coordinates.
}

\description{
  Estimates the transport term (i.e. the rate of change of a concentration
  due to diffusion) in a cylindrical (r, theta, z) or spherical (r, theta, phi) 
  coordinate system.
}

\usage{
tran.cylindrical (C, C.r.up = NULL, C.r.down = NULL, 
        C.theta.up = NULL, C.theta.down = NULL, 
        C.z.up = NULL, C.z.down = NULL, 
        flux.r.up = NULL, flux.r.down = NULL, 
        flux.theta.up = NULL, flux.theta.down = NULL,          
        flux.z.up = NULL, flux.z.down = NULL, 
        cyclicBnd = NULL,
        D.r = NULL, D.theta = D.r, D.z = D.r, 
        r = NULL,  theta = NULL, z = NULL)

tran.spherical (C, C.r.up = NULL, C.r.down = NULL, 
        C.theta.up = NULL, C.theta.down = NULL, 
        C.phi.up = NULL, C.phi.down = NULL, 
        flux.r.up = NULL, flux.r.down = NULL, 
        flux.theta.up = NULL, flux.theta.down = NULL,          
        flux.phi.up = NULL, flux.phi.down = NULL, 
        cyclicBnd = NULL,
        D.r = NULL, D.theta = D.r, D.phi = D.r, 
        r = NULL,  theta = NULL, phi = NULL)
}

\arguments{
  \item{C }{concentration, expressed per unit volume, defined at the centre
    of each grid cell; Nr*Nteta*Nz (cylindrica) or Nr*Ntheta*Nphi (spherical
    coordinates) array [M/L3].
  }
  \item{C.r.up }{concentration at upstream boundary in r(x)-direction; 
    one value [M/L3].
  }
  \item{C.r.down }{concentration at downstream boundary in r(x)-direction;
    one value [M/L3].
  }
  \item{C.theta.up }{concentration at upstream boundary in theta-direction;
    one value [M/L3].
  }
  \item{C.theta.down }{concentration at downstream boundary in theta-direction;
    one value [M/L3].
  }
  \item{C.z.up }{concentration at upstream boundary in z-direction (cylindrical
   coordinates); one value [M/L3].
  }
  \item{C.z.down }{concentration at downstream boundary in z-direction(cylindrical
   coordinates); one value [M/L3].
  }
  \item{C.phi.up }{concentration at upstream boundary in phi-direction (spherical
   coordinates); one value [M/L3].
  }
  \item{C.phi.down }{concentration at downstream boundary in phi-direction(spherical
   coordinates); one value [M/L3].
  }
  \item{flux.r.up }{flux across the upstream boundary in r-direction,
    positive = INTO model domain; one value [M/L2/T].
  }
  \item{flux.r.down }{flux across the downstream boundary in r-direction,
    positive = OUT of model domain; one value [M/L2/T].
  }
  \item{flux.theta.up }{flux across the upstream boundary in theta-direction,
    positive = INTO model domain; one value [M/L2/T].
  }
  \item{flux.theta.down }{flux across the downstream boundary in theta-direction,
    positive = OUT of model domain; one value [M/L2/T].
  }
  \item{flux.z.up }{flux across the upstream boundary in z-direction(cylindrical
   coordinates); positive = INTO model domain; one value [M/L2/T].
  }
  \item{flux.z.down }{flux across the downstream boundary in z-direction,
    (cylindrical coordinates); positive = OUT of model domain; one value [M/L2/T].
  }
  \item{flux.phi.up }{flux across the upstream boundary in phi-direction(spherical
   coordinates); positive = INTO model domain; one value [M/L2/T].
  }
  \item{flux.phi.down }{flux across the downstream boundary in phi-direction,
    (spherical coordinates); positive = OUT of model domain; one value [M/L2/T].
  }
  \item{cyclicBnd }{If not \code{NULL}, the direction in which a cyclic 
    boundary is defined, i.e. \code{cyclicBnd = 1} for the \code{r} direction, 
    \code{cyclicBnd = 2} for the \code{theta} direction and 
    \code{cyclicBnd = c(1,2)} for both the \code{r} and \code{theta} direction.
  }
  \item{D.r }{diffusion coefficient in r-direction, defined on grid cell
    interfaces. One value or a vector of length (Nr+1), [L2/T].
  }
  \item{D.theta }{diffusion coefficient in theta-direction, defined on grid cell
    interfaces. One value or or a vector of length (Ntheta+1), [L2/T].
  }
  \item{D.z }{diffusion coefficient in z-direction, defined on grid cell
    interfaces for cylindrical coordinates. One value or a vector of length 
    (Nz+1) [L2/T].
  }
  \item{D.phi }{diffusion coefficient in phi-direction, defined on grid cell
    interfaces for cylindrical coordinates. One value or a vector of length 
    (Nphi+1) [L2/T].
  }
  \item{r }{position of adjacent cell interfaces in the r-direction. 
    A vector of length Nr+1 [L].
  }
  \item{theta }{position of adjacent cell interfaces in the theta-direction.
    A vector of length Ntheta+1 [L].  Theta should be within [0,2 pi]
  }
  \item{z }{position of adjacent cell interfaces in the z-direction (cylindrical
   coordinates). A vector of length Nz+1 [L].
  }
  \item{phi }{position of adjacent cell interfaces in the phi-direction (spherical
   coordinates). A vector of length Nphi+1 [L]. Phi should be within [0,2 pi]
  }
}
\value{
  a list containing:
  \item{dC }{the rate of change of the concentration C due to transport,
    defined in the centre of each grid cell, a Nr*Nteta*Nz (cylindrical) or 
    Nr*Ntheta*Nphi (spherical coordinates) array. [M/L3/T].
  }
  \item{flux.r.up }{flux across the upstream boundary in r-direction,
    positive = INTO model domain. A matrix of dimension Nteta*Nz (cylindrical) 
    or Ntheta*Nphi (spherical) [M/L2/T].
  }
  \item{flux.r.down }{flux across the downstream boundary in r-direction,
    positive = OUT of model domain. A matrix of dimension Nteta*Nz (cylindrical)
    or Ntheta*Nphi (spherical) [M/L2/T].
  }
  \item{flux.theta.up }{flux across the upstream boundary in theta-direction,
    positive = INTO model domain. A matrix of dimension Nr*Nz (cylindrical) or 
    or Nr*Nphi (spherical) [M/L2/T].
  }
  \item{flux.theta.down }{flux across the downstream boundary in theta-direction,
   positive = OUT of model domain. A matrix of dimension Nr*Nz (cylindrical) or 
   Nr*Nphi (spherical) [M/L2/T].
  }
  \item{flux.z.up }{flux across the upstream boundary in z-direction,
   for cylindrical coordinates; 
   positive = OUT of model domain. A matrix of dimension Nr*Nteta [M/L2/T].
  }
  \item{flux.z.down }{flux across the downstream boundary in z-direction 
   for cylindrical coordinates; 
   positive = OUT of model domain. A matrix of dimension Nr*Nteta [M/L2/T].
  }
  \item{flux.phi.up }{flux across the upstream boundary in phi-direction,
   for spherical coordinates;
   positive = OUT of model domain. A matrix of dimension Nr*Nteta [M/L2/T].
  }
  \item{flux.phi.down }{flux across the downstream boundary in phi-direction,
   for spherical coordinates;
   positive = OUT of model domain. A matrix of dimension Nr*Nteta [M/L2/T].
  }

}

\examples{
## =============================================================================
## Testing the functions
## =============================================================================
# Grid definition
r.N     <- 4   # number of cells in r-direction
theta.N <- 6   # number of cells in theta-direction
z.N     <- 3   # number of cells in z-direction

D       <- 100 # diffusion coefficient
 
r      <- seq(0,   8,len=r.N+1)         # cell size r-direction [cm]
theta  <- seq(0,2*pi,len=theta.N+1)     # theta-direction - theta: from 0, 2pi
phi    <- seq(0,2*pi,len= z.N+1)        # phi-direction (0,2pi)
z      <- seq(0,5,len=z.N+1)            # cell size z-direction [cm]
 
# Intial conditions 
C <- array(dim = c(r.N, theta.N, z.N), data=0)

# Concentration boundary conditions
tran.cylindrical (C = C, D.r = D, D.theta = D, 
  C.r.up = 1, C.r.down = 1,
  C.theta.up = 1, C.theta.down = 1, 
  C.z.up = 1, C.z.down = 1,
  r = r, theta = theta, z = z )

tran.spherical (C = C, D.r = D, D.theta = D, 
  C.r.up = 1, C.r.down = 1, C.theta.up = 1, C.theta.down = 1, 
  C.phi.up = 1, C.phi.down = 1,
  r = r, theta = theta, phi = phi)

# Flux boundary conditions
tran.cylindrical(C = C, D.r = D, r = r, theta = theta, z = z,
  flux.r.up = 10, flux.r.down = 10,
  flux.theta.up = 10, flux.theta.down = 10,
  flux.z.up = 10, flux.z.down = 10)

tran.spherical(C = C, D.r = D, r = r, theta = theta, phi = phi,
  flux.r.up = 10, flux.r.down = 10,
  flux.theta.up = 10, flux.theta.down = 10,
  flux.phi.up = 10, flux.phi.down = 10)

# cyclic boundary conditions
tran.cylindrical(C = C, D.r = D, r = r, theta = theta, z = z,
  cyclicBnd = 1:3)
tran.spherical(C = C, D.r = D, r = r, theta = theta, phi = phi,
  cyclicBnd = 1:3)

# zero-gradient boundary conditions
tran.cylindrical(C = C, D.r = D, r = r, theta = theta, z = z)
tran.spherical(C = C, D.r = D, r = r, theta = theta, phi = phi)

## =============================================================================
## A model with diffusion and first-order consumption
## =============================================================================


N     <- 10          # number of grid cells
rr    <- 0.005       # consumption rate
D     <- 400

r       <- seq (2, 4, len=N+1)
theta   <- seq (0, 2*pi, len = N+1)
z       <- seq (0, 3, len = N+1)
phi     <- seq (0, 2*pi, len = N+1)

# The model equations
Diffcylin <- function (t, y, parms)  {
  CONC  <- array(dim = c(N, N, N), data = y)
  tran  <- tran.cylindrical(CONC, 
        D.r = D, D.theta = D, D.z = D,
        r = r, theta = theta, z = z,
        C.r.up = 0,  C.r.down = 1,
        cyclicBnd = 2)
  dCONC <- tran$dC  - rr * CONC
  return (list(dCONC))
}

Diffspher <- function (t, y, parms)  {
  CONC  <- array(dim = c(N, N, N), data = y)
  tran  <- tran.spherical (CONC, 
        D.r = D, D.theta = D, D.phi = D,
        r = r, theta = theta, phi = phi,
        C.r.up = 0,  C.r.down = 1,
        cyclicBnd = 2:3)
  dCONC <- tran$dC  - rr * CONC
  return (list(dCONC))
}

# initial condition: 0 everywhere, except in central point
y   <- array(dim = c(N, N, N), data = 0)
N2  <- ceiling(N/2)

y[N2, N2, N2] <- 100  # initial concentration in the central point...

# solve to steady-state; cyclicBnd = 2, 
outcyl <- steady.3D (y = y, func = Diffcylin, parms = NULL,
                  dim = c(N, N, N), lrw = 1e6, cyclicBnd = 2)

STDcyl <- array(dim = c(N, N, N), data = outcyl$y)
image(STDcyl[,,1])

# For spherical coordinates, cyclic Bnd = 2, 3
outspher <- steady.3D (y = y, func = Diffspher, parms = NULL, pos=TRUE,
                  dim = c(N, N, N), lrw = 1e6, cyclicBnd = 2:3)

#STDspher <- array(dim = c(N, N, N), data = outspher$y)
#image(STDspher[,,1])

\dontrun{
  image(outspher)
}          
}

%%\references{
%%}

\details{
  \code{tran.cylindrical} performs (diffusive) transport in cylindrical coordinates

  \code{tran.spherical} performs (diffusive) transport in spherical coordinates

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
  downstream boundary the same. It is useful mainly for the \code{theta} and 
  \code{phi} direction.

  ** Do  not expect too much of this equation: do not try to use it with 
  many boxes ** 
}  
\seealso{
  \code{\link{tran.polar}} 
   for a discretisation of 2-D transport equations in polar coordinates
  
  \code{\link{tran.1D}},   \code{\link{tran.2D}},  \code{\link{tran.3D}}
}

\keyword{utilities}

