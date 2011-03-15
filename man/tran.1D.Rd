\name{tran.1D}
\alias{tran.1D}

\title{
  General One-Dimensional Advective-Diffusive Transport
}

\description{
  Estimates the transport term (i.e. the rate of change of a concentration
  due to diffusion and advection) in a one-dimensional model of a
  liquid (volume fraction constant and equal to one) or in a porous medium
  (volume fraction variable and lower than one).

  The interfaces between grid cells can have a variable cross-sectional area,
  e.g. when modelling spherical or cylindrical geometries (see example).
}

\usage{
tran.1D(C, C.up = C[1], C.down = C[length(C)],
        flux.up = NULL, flux.down = NULL, 
        a.bl.up = NULL, a.bl.down = NULL, 
        D = 0, v = 0, AFDW = 1, VF = 1, A = 1, dx,
        full.check = FALSE, full.output = FALSE)
}

\arguments{
  \item{C }{concentration, expressed per unit of phase volume, defined at the
    centre of each grid cell. A vector of length N [M/L3]
  }
  \item{C.up }{concentration at upstream boundary. One value [M/L3]
  }
  \item{C.down }{concentration at downstream boundary. One value [M/L3]
  }
  \item{flux.up }{flux across the upstream boundary, positive = INTO model
    domain. One value, expressed per unit of total surface [M/L2/T]. 
    If \code{NULL}, the boundary is prescribed as
    a concentration or a convective transfer boundary.
  }
  \item{flux.down }{flux across the downstream boundary, positive = OUT
    of model domain. One value, expressed per unit of total surface [M/L2/T].
    If \code{NULL}, the boundary is prescribed as
    a concentration or a convective transfer boundary.
  }
  \item{a.bl.up }{convective transfer coefficient across the upstream
    boundary layer. \code{Flux = a.bl.up*(C.up-C0)}. One value [L/T]
  }
  \item{a.bl.down }{convective transfer coefficient across the downstream
    boundary layer (L). \code{Flux = a.bl.down*(CL-C.down)}.
    One value [L/T]
  }
  \item{D }{diffusion coefficient, defined on grid cell interfaces.
    One value, a vector of length N+1 [L2/T], or a \code{1D property} list; the list
    contains at least the element \code{int} (see \code{\link{setup.prop.1D}})
    [L2/T]
  }
  \item{v }{advective velocity, defined on the grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length N+1 [L/T], or a \code{1D property} list; the list
    contains at least the element \code{int} (see \code{\link{setup.prop.1D}})
    [L/T]
  }
  \item{AFDW }{weight used in the finite difference scheme for advection,
    defined on grid cell interfaces; backward = 1, centred = 0.5, forward = 0;
    default is backward. One value, a vector of length N+1, or a
    \code{1D property} list; the list contains at least the element \code{int}
    (see \code{\link{setup.prop.1D}}) [-]
  }
  \item{VF }{Volume fraction defined at the grid cell interfaces. One value,
    a vector of length N+1, or a \code{1D property} list; the list
    contains at least the elements \code{int} and \code{mid}
    (see \code{\link{setup.prop.1D}}) [-]
  }
  \item{A }{Interface area defined at the grid cell interfaces. One value,
    a vector of length N+1, or a \code{1D grid property} list; the list
    contains at least the elements \code{int} and \code{mid}
    (see \code{\link{setup.prop.1D}}) [L2]
  }
  \item{dx }{distance between adjacent cell interfaces (thickness of grid
    cells). One value, a vector of length N, or a \code{1D grid} list containing
    at least the elements
    \code{dx} and \code{dx.aux} (see \code{\link{setup.grid.1D}}) [L]
  }
  \item{full.check }{logical flag enabling a full check of the consistency
    of the arguments (default = \code{FALSE}; \code{TRUE} slows down execution
    by 50 percent)
  }
  \item{full.output }{logical flag enabling a full return of the output
    (default = \code{FALSE}; \code{TRUE} slows down execution by 20 percent)
  }
}

\value{
  \item{dC }{the rate of change of the concentration C due to transport,
    defined in the centre of each grid cell. The rate of change is expressed
    per unit of phase volume [M/L3/T]
  }
  \item{C.up }{concentration at the upstream interface. One value [M/L3]
    only when (\code{full.output} = \code{TRUE})
  }
  \item{C.down }{concentration at the downstream interface. One value [M/L3]
    only when (\code{full.output} = \code{TRUE})
  }
  \item{dif.flux }{diffusive flux across at the interface of each grid cell.
    A vector of length N+1 [M/L2/T]
    only when (\code{full.output} = \code{TRUE})
  }
  \item{adv.flux }{advective flux across at the interface of each grid cell.
    A vector of length N+1 [M/L2/T]
    only when (\code{full.output} = \code{TRUE})
  }
  \item{flux }{total flux across at the interface of each grid cell. A vector
    of length N+1 [M/L2/T].
    only when (\code{full.output} = \code{TRUE})
  }
  \item{flux.up }{flux across the upstream boundary, positive = INTO model
    domain. One value [M/L2/T]
  }
  \item{flux.down }{flux across the downstream boundary, positive = OUT of
    model domain. One value [M/L2/T]
  }
}
\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
## =============================================================================
## EXAMPLE 1: O2 and OC consumption in sediments
## =============================================================================

# this example uses only the volume fractions 
# in the reactive transport term

#====================#
# Model formulation  #
#====================#

# Monod consumption of oxygen (O2)

O2.model <- function (t = 0, O2, pars = NULL) {

  tran <- tran.1D(C = O2, C.up = C.ow.O2, D = D.grid, 
                  v = v.grid, VF = por.grid, dx = grid)$dC
  reac <- - R.O2*(O2/(Ks+O2))
  return(list(dCdt = tran + reac))
}

# First order consumption of organic carbon (OC)

OC.model <- function (t = 0, OC, pars = NULL) {

  tran <- tran.1D(C = OC, flux.up = F.OC, D = Db.grid,
                  v = v.grid, VF = svf.grid, dx = grid)$dC
  reac <- - k*OC
  return(list(dCdt = tran + reac))
}

#======================#
# Parameter definition #
#======================#

# Parameter values

F.OC    <- 25    # input flux organic carbon [micromol cm-2 yr-1]
C.ow.O2 <- 0.25  # concentration O2 in overlying water [micromol cm-3]
por     <- 0.8   # porosity
D       <- 400   # diffusion coefficient O2 [cm2 yr-1]
Db      <- 10    # mixing coefficient sediment [cm2 yr-1]
v       <- 1     # advective velocity [cm yr-1]
k       <- 1     # decay constant organic carbon [yr-1]
R.O2    <- 10    # O2 consumption rate [micromol cm-3 yr-1]
Ks      <- 0.005 # O2 consumption saturation constant 

# Grid definition

L <- 10   # depth of sediment domain [cm]
N <- 100  # number of grid layers
grid <- setup.grid.1D(x.up = 0, L = L, N = N)

# Volume fractions 

por.grid <- setup.prop.1D(value = por, grid = grid)
svf.grid <- setup.prop.1D(value = (1-por), grid = grid)
D.grid   <- setup.prop.1D(value = D, grid = grid)
Db.grid  <- setup.prop.1D(value = Db, grid = grid)
v.grid   <- setup.prop.1D(value = v, grid = grid)

#====================#
# Model solution     #
#====================#

# Initial conditions + simulation O2

yini <- rep(0, length.out = N) 
O2   <- steady.1D(y = yini, func = O2.model, nspec = 1)

# Initial conditions + simulation OC

yini <- rep(0, length.out = N) 
OC   <- steady.1D(y = yini, func = OC.model, nspec = 1)

# Plotting output, using S3 plot method of package rootSolve"

plot(O2, grid = grid$x.mid, xyswap = TRUE, main = "O2 concentration", 
     ylab = "depth [cm]", xlab = "", mfrow = c(1,2), type = "p", pch = 16)

plot(OC, grid = grid$x.mid, xyswap = TRUE, main = "C concentration", 
     ylab = "depth [cm]", xlab = "", mfrow = NULL)

## =============================================================================
## EXAMPLE 2: O2 in a cylindrical and spherical organism
## =============================================================================

# This example uses only the surface areas 
# in the reactive transport term

#====================#
# Model formulation  #
#====================#

# the numerical model - rate of change = transport-consumption
Cylinder.Model <- function(time, O2, pars)
  return (list(
    tran.1D(C = O2, C.down = BW, D = Da, A = A.cyl, dx = dx)$dC - Q
    ))

Sphere.Model <- function(time, O2, pars)
  return (list(
    tran.1D(C = O2, C.down = BW, D = Da, A = A.sphere, dx = dx)$dC - Q
    ))

#======================#
# Parameter definition #
#======================#

# parameter values

BW     <- 2      # mmol/m3,  oxygen conc in surrounding water
Da     <- 0.5    # cm2/d     effective diffusion coeff in organism
R      <- 0.0025 # cm        radius of organism
Q      <- 250000 # nM/cm3/d  oxygen consumption rate/ volume / day
L      <- 0.05   # cm        length of organism (if a cylinder)

# the numerical model

N  <- 40                              # layers in the body
dx <- R/N                             # thickness of each layer
x.mid <- seq(dx/2, by = dx, length.out = N) # distance of center to mid-layer
x.int <- seq(0, by = dx, length.out = N+1)  # distance to layer interface

# Cylindrical surfaces
A.cyl   <- 2*pi*x.int*L  # surface at mid-layer depth

# Spherical surfaces
A.sphere <- 4*pi*x.int^2 # surface of sphere, at each mid-layer

#====================#
# Model solution     #
#====================#

# the analytical solution of cylindrical and spherical model
cylinder <- function(Da, Q, BW, R, r)  BW + Q/(4*Da)*(r^2-R^2)
sphere   <- function(Da, Q, BW, R, r)  BW + Q/(6*Da)*(r^2-R^2)

# solve the model numerically for a cylinder
O2.cyl <- steady.1D (y = runif(N), name = "O2", 
      func = Cylinder.Model, nspec = 1, atol = 1e-10)

# solve the model numerically for a sphere
O2.sphere <- steady.1D (y = runif(N), name = "O2", 
      func = Sphere.Model, nspec = 1, atol = 1e-10)

#====================#
# Plotting output    #
#====================#
# Analytical solution - "observations"
Ana.cyl   <- cbind(x.mid, O2 = cylinder(Da, Q, BW, R, x.mid))
Ana.spher <- cbind(x.mid, O2 = sphere(Da, Q, BW, R, x.mid))

plot(O2.cyl, O2.sphere, grid = x.mid, lwd = 2, lty = 1, col = 1:2, 
     xlab = "distance from centre, cm", 
     ylab = "mmol/m3", main = "tran.1D",
     sub = "diffusion-reaction in a cylinder and sphere",
     obs = list(Ana.cyl, Ana.spher), obspar = list(pch = 16, col =1:2))
     
legend ("topleft", lty = c(1, NA), pch = c(NA, 18),
        c("numerical approximation", "analytical solution"))
legend ("bottomright", pch = 16, lty = 1, col = 1:2,
        c("cylinder", "sphere"))

## =============================================================================
## EXAMPLE 3: O2 consumption in a spherical aggregate
## =============================================================================

# this example uses both the surface areas and the volume fractions
# in the reactive transport term

#====================#
# Model formulation  #
#====================#

Aggregate.Model <- function(time, O2, pars) {

  tran <- tran.1D(C = O2, C.down = C.ow.O2,
                  D = D.grid, A = A.grid,
                  VF = por.grid, dx = grid )$dC

  reac <- - R.O2*(O2/(Ks+O2))*(O2>0)
  return(list(dCdt = tran + reac, consumption = -reac))

}

#======================#
# Parameter definition #
#======================#

# Parameters

C.ow.O2 <- 0.25     # concentration O2 water [micromol cm-3]
por     <- 0.8      # porosity
D       <- 400      # diffusion coefficient O2 [cm2 yr-1]
v       <- 0        # advective velocity [cm yr-1]
R.O2    <- 1000000  # O2 consumption rate [micromol cm-3 yr-1]
Ks      <- 0.005    # O2 saturation constant [micromol cm-3]

# Grid definition
R <- 0.025           # radius of the agggregate [cm]
N <- 100             # number of grid layers
grid <- setup.grid.1D(x.up = 0, L = R, N = N)

# Volume fractions 

por.grid <- setup.prop.1D(value = por, grid = grid)
D.grid   <- setup.prop.1D(value = D, grid = grid)

# Surfaces 

A.mid <- 4*pi*grid$x.mid^2  # surface of sphere at middle of grid cells
A.int <- 4*pi*grid$x.int^2  # surface of sphere at interface
A.grid <- list(int = A.int, mid = A.mid)

#====================#
# Model solution     #
#====================#

# Numerical solution: staedy state 

O2.agg <- steady.1D (runif(N), func = Aggregate.Model, nspec = 1,
                     atol = 1e-10, names = "O2")

#====================#
# Plotting output    #
#====================#

par(mfrow = c(1,1))

plot(grid$x.mid, O2.agg$y, xlab = "distance from centre, cm",
     ylab = "mmol/m3",
     main = "Diffusion-reaction of O2 in a spherical aggregate")
legend ("bottomright", pch = c(1, 18), lty = 1, col = "black",
        c("O2 concentration"))

# Similar, using S3 plot method of package rootSolve"
plot(O2.agg, grid = grid$x.mid, which = c("O2", "consumption"),
     xlab = "distance from centre, cm", ylab = c("mmol/m3","mmol/m3/d")) 

}
\references{
  Soetaert and Herman (2009). A practical guide to ecological modelling -
  using R as a simulation platform. Springer
}

\details{
  The \bold{boundary conditions} are either
  \itemize{
    \item (1) zero-gradient.
    \item (2) fixed concentration.
    \item (3) convective boundary layer.
    \item (4) fixed flux.
  }
  The above order also shows the priority. The default condition is the
  zero gradient. The fixed concentration condition overrules the zero gradient.
  The convective boundary layer condition overrules the fixed concentration
  and zero gradient. The fixed flux overrules all other specifications.

  Ensure that the boundary conditions are well defined: for instance, it 
  does not make sense to specify an influx in a boundary cell with the advection
  velocity pointing outward.
  
  \bold{Transport properties:}
  
  The \emph{diffusion coefficient} (\code{D}),
  the \emph{advective velocity} (\code{v}),
  the \emph{volume fraction} (VF), the \emph{interface surface} (\code{A}),
  and the \emph{advective finite difference weight} (\code{AFDW})
  can either be specified as one value, a vector or a 1D property list
  as generated by \code{\link{setup.prop.1D}}.

  When a vector, this vector must be of length N+1, defined at all grid
  cell interfaces, including the upper and lower boundary.

  The \bold{finite difference grid} (\code{dx}) is specified either as
  one value, a vector or a 1D grid list, as generated by \code{\link{setup.grid.1D}}. 
}
\note{ 
The advective equation is not checked for mass conservation. Sometimes, this is 
not an issue, for instance when \code{v} represents a sinking velocity of 
particles or a swimming velocity of organisms. 
In others cases however, mass conservation needs to be accounted for. 
To ensure mass conservation, the advective velocity must obey certain 
continuity constraints: in essence the product of the volume fraction (VF), 
interface surface area (A) and advective velocity (v) should be constant. 
In sediments, one can use \code{\link{setup.compaction.1D}} to ensure that 
the advective velocities for the pore water and solid phase meet these 
constraints. 

In terms of the units of concentrations and fluxes we follow the convention 
in the geosciences. 
The concentration \code{C}, \code{C.up}, \code{C.down} as well at the rate of 
change of the concentration \code{dC} are always expressed per unit of 
phase volume (i.e. per unit volume of solid or liquid). 

Total concentrations (e.g. per unit volume of bulk sediment) can be obtained by 
multiplication with the appropriate volume fraction. In contrast, fluxes are 
always expressed per unit of total interface area (so here the volume fraction 
is accounted for).     
}
\seealso{
  \code{\link{tran.volume.1D}} for a discretisation the transport equation using finite volumes.

  \code{\link{tran.2D}},   \code{\link{tran.3D}}

  \code{\link{advection.1D}}, for more sophisticated advection schemes

}

\keyword{utilities}

