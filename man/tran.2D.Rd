\name{tran.2D}
\alias{tran.2D}
\title{
  General Two-Dimensional Advective-Diffusive Transport
}

\description{
  Estimates the transport term (i.e. the rate of change of a concentration
  due to diffusion and advection) in a two-dimensional model domain.
}

\usage{
tran.2D (C, C.x.up = C[1,], C.x.down = C[nrow(C),],
         C.y.up = C[,1], C.y.down = C[ ,ncol(C)],
         flux.x.up = NULL, flux.x.down = NULL, 
         flux.y.up = NULL, flux.y.down = NULL,
         a.bl.x.up = NULL, a.bl.x.down = NULL, 
         a.bl.y.up = NULL, a.bl.y.down = NULL, 
         D.grid = NULL, D.x = NULL, D.y = D.x,
         v.grid = NULL, v.x = 0, v.y = 0,
         AFDW.grid = NULL, AFDW.x = 1, AFDW.y = AFDW.x,
         VF.grid = NULL, VF.x = 1, VF.y = VF.x,
         A.grid = NULL, A.x = 1, A.y = 1,
         grid = NULL, dx = NULL, dy = NULL,
         full.check = FALSE, full.output = FALSE)
}

\arguments{
  \item{C }{concentration, expressed per unit volume, defined at the centre
    of each grid cell; Nx*Ny matrix [M/L3].
  }
  \item{C.x.up }{concentration at upstream boundary in x-direction;
    vector of length Ny [M/L3].
  }
  \item{C.x.down }{concentration at downstream boundary in x-direction;
    vector of length Ny [M/L3].
  }
  \item{C.y.up }{concentration at upstream boundary in y-direction;
    vector of length Nx [M/L3].
  }
  \item{C.y.down }{concentration at downstream boundary in y-direction;
    vector of length Nx [M/L3].
  }
  \item{flux.x.up }{flux across the upstream boundary in x-direction,
    positive = INTO model domain; vector of length Ny [M/L2/T].
  }
  \item{flux.x.down }{flux across the downstream boundary in x-direction,
    positive = OUT of model domain; vector of length Ny [M/L2/T].
  }
  \item{flux.y.up }{flux across the upstream boundary in y-direction,
    positive = INTO model domain; vector of length Nx [M/L2/T].
  }
  \item{flux.y.down }{flux across the downstream boundary in y-direction,
    positive = OUT of model domain; vector of length Nx [M/L2/T].
  }
  \item{a.bl.x.up }{transfer coefficient across the upstream boundary layer.
     in x-direction;

     \code{Flux=a.bl.x.up*(C.x.up-C[1,])}. One value [L/T].
  }
  \item{a.bl.x.down }{transfer coefficient across the downstream boundary
    layer in x-direction;

    \code{Flux=a.bl.x.down*(C[Nx,]-C.x.down)}.
    One value [L/T].
  }
  \item{a.bl.y.up }{transfer coefficient across the upstream boundary layer.
     in y-direction;

     \code{Flux=a.bl.y.up*(C.y.up-C[,1])}. One value [L/T].
  }
  \item{a.bl.y.down }{transfer coefficient across the downstream boundary
    layer in y-direction;

    \code{Flux=a.bl.y.down*(C[,Ny]-C.y.down)}.
    One value [L/T].
  }
  \item{D.grid }{diffusion coefficient defined on all grid cell
    interfaces. A \code{prop.2D} list created by \code{\link{setup.prop.2D}} [L2/T].
    See last example for creating spatially-varying diffusion coefficients.
  }
  \item{D.x }{diffusion coefficient in x-direction, defined on grid cell
    interfaces. One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nx+1)* Ny matrix [L2/T].
  }
  \item{D.y }{diffusion coefficient in y-direction, defined on grid cell
    interfaces. One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*(Ny+1) matrix [L2/T].
  }
  \item{v.grid }{advective velocity defined on all grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    A \code{prop.2D} list created by \code{\link{setup.prop.2D}} [L/T].
  }
  \item{v.x }{advective velocity in the x-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nx+1)*Ny matrix [L/T].
  }
  \item{v.y }{advective velocity in the y-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*(Ny+1) matrix [L/T].
  }
  \item{AFDW.grid }{weight used in the finite difference scheme for advection
    in the x- and y- direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    A \code{prop.2D} list created by \code{\link{setup.prop.2D}} [-].
  }
  \item{AFDW.x }{weight used in the finite difference scheme for advection
    in the x-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nx+1)*Ny matrix [-].
  }
  \item{AFDW.y }{weight used in the finite difference scheme for advection
    in the y-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*(Ny+1) matrix [-].
  }
  \item{VF.grid }{Volume fraction. A \code{prop.2D} list created by
    \code{\link{setup.prop.2D}} [-].
  }
  \item{VF.x }{Volume fraction at the grid cell interfaces in the x-direction.
    One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nx+1)*Ny matrix [-].
  }
  \item{VF.y }{Volume fraction at the grid cell interfaces in the y-direction.
    One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*(Ny+1) matrix [-].
  }
  \item{A.grid }{Interface area. A \code{prop.2D} list created by
    \code{\link{setup.prop.2D}} [L2].
  }
  \item{A.x }{Interface area defined at the grid cell interfaces in
    the x-direction. One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nx+1)*Ny matrix [L2].
  }
  \item{A.y }{Interface area defined at the grid cell interfaces in
    the y-direction. One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*(Ny+1) matrix [L2].
  }
  \item{dx }{distance between adjacent cell interfaces in the x-direction
    (thickness of grid cells). One value or vector of length Nx [L].
  }
  \item{dy }{distance between adjacent cell interfaces in the y-direction
    (thickness of grid cells). One value or vector of length Ny [L].
  }
  \item{grid }{discretization grid, a list containing at least elements
    \code{dx}, \code{dx.aux}, \code{dy}, \code{dy.aux}
    (see \code{\link{setup.grid.2D}}) [L].
  }
  \item{full.check }{logical flag enabling a full check of the consistency
    of the arguments (default = \code{FALSE}; \code{TRUE} slows down
    execution by 50 percent).
  }
  \item{full.output }{logical flag enabling a full return of the output
    (default = \code{FALSE}; \code{TRUE} slows down execution by 20 percent).
  }
}
\value{
  a list containing:
  \item{dC }{the rate of change of the concentration C due to transport,
    defined in the centre of each grid cell, a Nx*Ny matrix. [M/L3/T].
  }
  \item{C.x.up }{concentration at the upstream interface in x-direction.
    A vector of length Ny [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{C.x.down }{concentration at the downstream interface in x-direction.
    A vector of length Ny [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{C.y.up }{concentration at the the upstream interface in y-direction.
    A vector of length Nx [M/L3].  Only when \code{full.output = TRUE}.
  }
  \item{C.y.down }{concentration at the downstream interface in y-direction.
    A vector of length Nx [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{x.flux }{flux across the interfaces in x-direction of the grid cells.
    A (Nx+1)*Ny matrix  [M/L2/T]. Only when \code{full.output = TRUE}.
  }
  \item{y.flux }{flux across the interfaces in y-direction of the grid cells.
    A Nx*(Ny+1) matrix [M/L2/T].  Only when \code{full.output = TRUE}.
  }
  \item{flux.x.up }{flux across the upstream boundary in x-direction,
    positive = INTO model domain. A vector of length Ny [M/L2/T].
  }
  \item{flux.x.down }{flux across the downstream boundary in x-direction,
    positive = OUT of model domain. A vector of length Ny [M/L2/T].
  }
  \item{flux.y.up }{flux across the upstream boundary in y-direction,
    positive = INTO model domain. A vector of length Nx [M/L2/T].
  }
  \item{flux.y.down }{flux across the downstream boundary in y-direction,
    positive = OUT of model domain. A vector of length Nx [M/L2/T].
  }

}

\author{
  Filip Meysman <filip.meysman@nioz.nl>,
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\note{
  It is much more efficient to use the \emph{grid} input rather than 
  vectors or single numbers. 
  
  Thus: to optimise the code, use \link{setup.grid.2D} to create the 
  \code{grid}, and use \link{setup.prop.2D} to create \code{D.grid},
  \code{v.grid}, \code{AFDW.grid}, \code{VF.grid}, and \code{A.grid},
  even if the values are 1 or remain constant.
  
  There is no provision (yet) to deal with \emph{cross-diffusion}. 
  Set \code{D.x} and \code{D.y} different only if cross-diffusion effects
  are unimportant. 
  
}

\examples{
## =============================================================================
## Testing the functions
## =============================================================================
# Parameters
F        <- 100             # input flux [micromol cm-2 yr-1]
por      <- 0.8             # constant porosity
D        <- 400             # mixing coefficient [cm2 yr-1]
v        <- 1               # advective velocity [cm yr-1]

# Grid definition
x.N <- 4   # number of cells in x-direction
y.N <- 6   # number of cells in y-direction
x.L <- 8   # domain size x-direction [cm]
y.L <- 24  # domain size y-direction [cm]
dx  <- x.L/x.N             # cell size x-direction [cm]
dy  <- y.L/y.N             # cell size y-direction [cm]
 
# Intial conditions 
C <- matrix(nrow = x.N, ncol = y.N, data = 0, byrow = FALSE)

# Boundary conditions: fixed concentration  
C.x.up   <- rep(1, times = y.N)
C.x.down <- rep(0, times = y.N)
C.y.up   <- rep(1, times = x.N)
C.y.down <- rep(0, times = x.N)

# Only diffusion 
tran.2D(C = C, D.x = D, D.y = D, v.x = 0, v.y = 0,
  VF.x = por, VF.y = por, dx = dx, dy = dy,
  C.x.up = C.x.up, C.x.down = C.x.down,
  C.y.up = C.y.up, C.y.down = C.y.down, full.output = TRUE)

# Strong advection, backward (default), central and forward 
#finite difference schemes 
tran.2D(C = C, D.x = D, v.x = 100*v, 
  VF.x = por, dx = dx, dy = dy,
  C.x.up = C.x.up, C.x.down = C.x.down, 
  C.y.up = C.y.up, C.y.down = C.y.down)
  
tran.2D(AFDW.x = 0.5, C = C, D.x = D, v.x = 100*v, 
  VF.x = por, dx = dx, dy = dy,
  C.x.up = C.x.up, C.x.down = C.x.down, 
  C.y.up = C.y.up, C.y.down = C.y.down)

tran.2D(AFDW.x = 0, C = C, D.x = D, v.x = 100*v, 
  VF.x = por, dx = dx, dy = dy,
  C.x.up = C.x.up, C.x.down = C.x.down, 
  C.y.up = C.y.up, C.y.down = C.y.down)

# Boundary conditions: fixed fluxes 

flux.x.up <- rep(200, times = y.N)
flux.x.down <- rep(-200, times = y.N)
flux.y.up <- rep(200, times = x.N)
flux.y.down <- rep(-200, times = x.N)
tran.2D(C = C, D.x = D, v.x = 0, 
  VF.x = por, dx = dx, dy = dy,
  flux.x.up = flux.x.up, flux.x.down = flux.x.down,
  flux.y.up = flux.y.up, flux.y.down = flux.y.down)

# Boundary conditions: convective boundary layer on all sides

a.bl <- 800   # transfer coefficient
C.x.up <- rep(1, times = (y.N)) # fixed conc at boundary layer
C.y.up <- rep(1, times = (x.N)) # fixed conc at boundary layer
tran.2D(full.output = TRUE, C = C, D.x = D, v.x = 0, 
  VF.x = por, dx = dx, dy = dy, 
  C.x.up   = C.x.up, a.bl.x.up = a.bl, 
  C.x.down = C.x.up, a.bl.x.down = a.bl, 
  C.y.up   = C.y.up, a.bl.y.up = a.bl,
  C.y.down = C.y.up, a.bl.y.down = a.bl)

# Runtime test with and without argument checking

n.iterate <-500

test1 <- function() {
  for (i in 1:n.iterate )
    ST <- tran.2D(full.check = TRUE, C = C, D.x = D, 
      v.x = 0, VF.x = por, dx = dx, dy = dy,
      C.x.up = C.x.up, a.bl.x.up = a.bl, C.x.down = C.x.down)
} 
system.time(test1())

test2 <- function() {
  for (i in 1:n.iterate )
    ST <- tran.2D(full.output = TRUE, C = C, D.x = D, 
      v.x = 0, VF.x = por, dx = dx, dy = dy,
      C.x.up = C.x.up, a.bl.x.up = a.bl, C.x.down = C.x.down)
} 
system.time(test2())

test3 <- function() {
  for (i in 1:n.iterate )
    ST <- tran.2D(full.output = TRUE, full.check = TRUE, C = C,
      D.x = D, v.x = 0, VF.x = por, dx = dx, dy = dy,
      C.x.up = C.x.up, a.bl.x.up = a.bl, C.x.down = C.x.down)
} 
system.time(test3())

## =============================================================================
## A 2-D model with diffusion in x- and y direction and first-order
## consumption - unefficient implementation
## =============================================================================

N     <- 51          # number of grid cells
XX    <- 10           # total size
dy    <- dx <- XX/N  # grid size
Dy    <- Dx <- 0.1   # diffusion coeff, X- and Y-direction
r     <- 0.005       # consumption rate
ini   <- 1           # initial value at x=0

N2  <- ceiling(N/2)
X   <- seq (dx, by = dx, len = (N2-1))
X   <- c(-rev(X), 0, X)

# The model equations

Diff2D <- function (t, y, parms)  {

 CONC  <- matrix(nrow = N, ncol = N, y)
 dCONC <- tran.2D(CONC, D.x = Dx, D.y = Dy, dx = dx, dy = dy)$dC + r * CONC

 return (list(dCONC))

}

# initial condition: 0 everywhere, except in central point
y <- matrix(nrow = N, ncol = N, data = 0)
y[N2, N2] <- ini  # initial concentration in the central point...

# solve for 10 time units
times <- 0:10
out <- ode.2D (y = y, func = Diff2D, t = times, parms = NULL,
                dim = c(N,N), lrw = 160000)

pm <- par (mfrow = c(2, 2))

# Compare solution with analytical solution...
for (i in seq(2, 11, by = 3))  {
  tt   <- times[i]
  mat  <-  matrix(nrow = N, ncol = N, 
                  data = subset(out, time == tt))
  plot(X, mat[N2,], type = "l", main = paste("time=", times[i]),
       ylab = "Conc", col = "red")
  ana <- ini*dx^2/(4*pi*Dx*tt)*exp(r*tt-X^2/(4*Dx*tt))
  points(X, ana, pch = "+")
}

legend ("bottom", col = c("red","black"), lty = c(1, NA), 
  pch = c(NA, "+"), c("tran.2D", "exact"))
par("mfrow" = pm )



## =============================================================================
## A 2-D model with diffusion in x- and y direction and first-order
## consumption - more efficient implementation, specifying ALL 2-D grids
## =============================================================================

N     <- 51          # number of grid cells
Dy    <- Dx <- 0.1   # diffusion coeff, X- and Y-direction
r     <- 0.005       # consumption rate
ini   <- 1           # initial value at x=0

x.grid    <- setup.grid.1D(x.up = -5, x.down = 5, N = N)
y.grid    <- setup.grid.1D(x.up = -5, x.down = 5, N = N)
grid2D    <- setup.grid.2D(x.grid, y.grid)

D.grid    <- setup.prop.2D(value = Dx, y.value = Dy, grid = grid2D)
v.grid    <- setup.prop.2D(value = 0, grid = grid2D)
A.grid    <- setup.prop.2D(value = 1, grid = grid2D)
AFDW.grid <- setup.prop.2D(value = 1, grid = grid2D)
VF.grid   <- setup.prop.2D(value = 1, grid = grid2D)

# The model equations - using the grids

Diff2Db <- function (t, y, parms)  {

   CONC  <- matrix(nrow = N, ncol = N, data = y)

   dCONC <- tran.2D(CONC, grid = grid2D, D.grid = D.grid, 
      A.grid = A.grid, VF.grid = VF.grid, AFDW.grid = AFDW.grid, 
      v.grid = v.grid)$dC + r * CONC
  
  return (list(dCONC))
}

# initial condition: 0 everywhere, except in central point
y <- matrix(nrow = N, ncol = N, data = 0)
y[N2,N2] <- ini  # initial concentration in the central point...

# solve for 8 time units
times <- 0:8
outb <- ode.2D (y = y, func = Diff2Db, t = times, parms = NULL,
                dim = c(N, N), lrw = 160000)

image(outb, ask = FALSE, mfrow = c(3, 3), main = paste("time", times))

## =============================================================================
## Same 2-D model, but now with spatially-variable diffusion coefficients
## =============================================================================

N     <- 51          # number of grid cells
r     <- 0.005       # consumption rate
ini   <- 1           # initial value at x=0
N2    <- ceiling(N/2)

D.grid <- list()

# Diffusion on x-interfaces
D.grid$x.int <- matrix(nrow = N+1, ncol = N, data = runif(N*(N+1)))

# Diffusion on y-interfaces
D.grid$y.int <- matrix(nrow = N, ncol = N+1, data = runif(N*(N+1)))

dx <- 10/N
dy <- 10/N

# The model equations

Diff2Dc <- function (t, y, parms)  {

   CONC  <- matrix(nrow = N, ncol = N, data = y)

   dCONC <- tran.2D(CONC, dx = dx, dy = dy, D.grid = D.grid)$dC + r * CONC

  return (list(dCONC))
}

# initial condition: 0 everywhere, except in central point
y <- matrix(nrow = N, ncol = N, data = 0)
y[N2, N2] <- ini  # initial concentration in the central point...

# solve for 8 time units
times <- 0:8
outc <- ode.2D (y = y, func = Diff2Dc, t = times, parms = NULL,
                dim = c(N, N), lrw = 160000)

outtimes <- c(1, 3, 5, 7)
image(outc, ask = FALSE, mfrow = c(2, 2), main = paste("time", outtimes),
      legend = TRUE, add.contour = TRUE, subset = time \%in\% outtimes)

}
\references{
  Soetaert and Herman, 2009. a practical guide to ecological modelling -
  using R as a simulation platform. Springer
}

\details{
  The \bold{boundary conditions} are either
  \itemize{
    \item (1) zero-gradient
    \item (2) fixed concentration
    \item (3) convective boundary layer
    \item (4) fixed flux
  }
  This is also the order of priority. The zero gradient is the default,
  the fixed flux overrules all other.
}  


\seealso{
  \code{\link{tran.polar}} for a discretisation of 2-D transport equations 
  in polar coordinates
  
  \code{\link{tran.1D}},   \code{\link{tran.3D}}
}
\keyword{utilities}
