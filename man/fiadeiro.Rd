\name{fiadeiro}
\alias{fiadeiro}
\title{
  Advective Finite Difference Weights
}

\description{
  Weighing coefficients used in the finite difference scheme for advection
  calculated according to Fiadeiro and Veronis (1977).

  This particular AFDW (advective finite difference weights) scheme switches
  from backward differencing (in advection dominated conditions; large Peclet
  numbers) to central differencing (under diffusion dominated conditions;
  small Peclet numbers).

  This way it forms a compromise between stability, accuracy and reduced
  numerical dispersion.
}

\usage{
fiadeiro(v, D, dx.aux = NULL, grid = list(dx.aux = dx.aux))
}

\arguments{
  \item{v }{advective velocity; either one value or a vector of length N+1,
    with N the number of grid cells [L/T]
  }
  \item{D }{diffusion coefficient; either one value or a vector of length N+1
    [L2/T]
  }
  \item{dx.aux }{auxiliary vector containing the distances between the
    locations where the concentration is defined (i.e. the grid cell centers
    and the two outer interfaces);
    either one value or a vector of length N+1 [L]
  }
  \item{grid }{discretization grid as calculated by \code{\link{setup.grid.1D}}
  }
}

\value{
  the Advective Finite Difference Weighing (AFDW) coefficients as used in 
  the transport routines \code{\link{tran.1D}} and \code{\link{tran.volume.1D}};
  either one value or a vector of length N+1 [-]
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
#============================================
# Model formulation (differential equations)
#============================================

# This is a test model to evaluate the different finite difference schemes 
# and evaluate their effect on munerical diffusion. The model describes the
# decay of organic carbon (OC) as it settles through the ocean water column.

model <- function (time, OC, pars, AFDW = 1) {
  dOC <- tran.1D(OC, flux.up = F_OC, D = D.eddy, 
                 v = v.sink, AFDW = AFDW, dx = dx)$dC - k*OC
  return(list(dOC))
}
#============================================
# Parameter set
#============================================

L <- 1000         # water depth model domain [m]
x.att <- 200      # attenuation depth of the sinking velocity [m]
v.sink.0 <- 10    # sinking velocity at the surface [m d-1]
D.eddy <- 10      # eddy diffusion coefficient [m2 d-1]
F_OC <- 10        # particle flux [mol m-2 d-1]
k <- 0.1          # decay coefficient [d-1]

## =============================================================================
## Model solution for a coarse grid (10 grid cells)
## =============================================================================

# Setting up the grid
N <- 10                               # number of grid layers 
dx <- L/N                             # thickness of boxes [m]
dx.aux <- rep(dx, N+1)                # auxilliary grid vector
x.int <- seq(from = 0, to = L, by = dx)    # water depth at box interfaces [m]
x.mid <- seq(from = dx/2, to = L, by = dx) # water depth at box centres [m]

# Exponentially declining sink velocity
v.sink <- v.sink.0 * exp(-x.int/x.att) # sink velocity [m d-1]
Pe <- v.sink * dx/D.eddy               # Peclet number

# Calculate the weighing coefficients
AFDW <- fiadeiro(v = v.sink, D = D.eddy, dx.aux = dx.aux)

par(mfrow = c(2, 1), cex.main = 1.2, cex.lab = 1.2)

# Plot the Peclet number over the grid 

plot(Pe, x.int, log = "x", pch = 19, ylim = c(L,0), xlim = c(0.1, 1000), 
     xlab = "", ylab = "depth [m]", 
     main = "Peclet number", axes = FALSE)
abline(h = 0)
axis(pos = NA, side = 2)
axis(pos = 0, side = 3)

# Plot the AFDW coefficients over the grid 

plot(AFDW, x.int, pch = 19, ylim = c(L, 0), xlim = c(0.5, 1), 
     xlab = "", ylab = "depth [m]", main = "AFDW coefficient", axes = FALSE)
abline(h = 0)
axis(pos = NA, side = 2)
axis(pos = 0, side = 3)

# Three steady-state solutions for a coarse grid based on:
# (1) backward differences (BD)
# (2) central differences (CD)
# (3) Fiadeiro & Veronis scheme (FV)

BD <- steady.1D(y = runif(N), func = model, AFDW = 1.0, nspec = 1)
CD <- steady.1D(y = runif(N), func = model, AFDW = 0.5, nspec = 1)
FV <- steady.1D(y = runif(N), func = model, AFDW = AFDW, nspec = 1)

# Plotting output - use rootSolve's plot method
plot(BD, CD, FV, grid = x.mid, xyswap = TRUE, mfrow = c(1,2), 
     xlab = "", ylab = "depth [m]", main = "conc (Low resolution grid)")

legend("bottomright", col = 1:3, lty = 1:3,
       legend = c("backward diff", "centred diff", "Fiadeiro&Veronis"))


## =============================================================================
## Model solution for a fine grid (1000 grid cells)
## =============================================================================

# Setting up the grid
N <- 1000                            # number of grid layers 
dx <- L/N                            # thickness of boxes[m]
dx.aux <- rep(dx, N+1)              # auxilliary grid vector
x.int <- seq(from = 0, to = L, by = dx)      # water depth at box interfaces [m]
x.mid <- seq(from = dx/2, to = L, by = dx)   # water depth at box centres [m]

# Exponetially declining sink velocity
v.sink <- v.sink.0 * exp(-x.int/x.att) # sink velocity [m d-1]
Pe <- v.sink * dx/D.eddy               # Peclet number

# Calculate the weighing coefficients
AFDW <- fiadeiro(v = v.sink, D = D.eddy, dx.aux = dx.aux)

# Three steady-state solutions for a coarse grid based on:
# (1) backward differences (BD)
# (2) centered differences (CD)
# (3) Fiadeiro & Veronis scheme (FV)

BD <- steady.1D(y = runif(N), func = model, AFDW = 1.0, nspec = 1)
CD <- steady.1D(y = runif(N), func = model, AFDW = 0.5, nspec = 1)
FV <- steady.1D(y = runif(N), func = model, AFDW = AFDW, nspec = 1)

# Plotting output
plot(BD, CD, FV, grid = x.mid, xyswap = TRUE, mfrow = NULL, 
     xlab = "", ylab = "depth [m]", main = "conc (High resolution grid)")

legend("bottomright", col = 1:3, lty = 1:3,
       legend = c("backward diff", "centred diff", "Fiadeiro&Veronis"))

# Results and conclusions:
# - For the fine grid, all three solutions are identical
# - For the coarse grid, the BD and FV solutions show numerical dispersion
#   while the CD provides more accurate results
}

\references{
  \itemize{
    \item Fiadeiro ME and Veronis G (1977) Weighted-mean schemes for
      finite-difference approximation to advection-diffusion equation.
      Tellus 29, 512-522.
    \item Boudreau (1997) Diagnetic models and their implementation.
      Chapter 8: Numerical Methods. Springer.
  }
}

\details{
  The Fiadeiro and Veronis (1977) scheme adapts the differencing method to
  the local situation (checks for advection or diffusion dominance).

  Finite difference schemes are based on following rationale:
  \itemize{
    \item When using forward differences (AFDW = 0), the scheme is first
      order accurate, creates a low level of (artificial) numerical dispersion,
      but is highly unstable (state variables may become negative).
    \item When using backward differences (AFDW = 1), the scheme is first
      order accurate, is universally stable (state variables always remain
      positive), but the scheme creates high levels of numerical dispersion.
    \item When using central differences (AFDW = 0.5), the scheme is second
      order accurate, is not universally stable, and has a moderate level of
      numerical dispersion, but state variables may become negative.
  }

  Because of the instability issue, forward schemes should be avoided.
  Because of the higher accuracy, the central scheme is preferred over the
  backward scheme.

  The central scheme is stable when sufficient physical dispersion is present,
  it may become unstable when advection is the only transport process.

  The Fiadeiro and Veronis (1977) scheme takes this into account: it uses
  central differencing when possible (when physical dispersion is high enough),
  and switches to backward differing when needed (when advection dominates).
  The switching is determined by the Peclet number
  
  \code{Pe = abs(v)*dx.aux/D}
  
  \itemize{
    \item the higher the diffusion \code{D} (\code{Pe > 1}), the closer the
      AFDW coefficients are to 0.5 (central differencing)
    \item the higher the advection \code{v} (\code{Pe < 1}), the closer the
      AFDW coefficients are to 1 (backward differencing)
  }
}

\note{
  \itemize{
    \item If the state variables (concentrations) decline in the direction
      of the 1D axis, then the central difference scheme will be stable.
      If this is known a prioiri, then central differencing is
      preferred over the fiadeiro scheme.
    \item Each scheme will always create some numerical diffusion. This
      principally depends on the resolution of the grid (i.e. larger 
      \code{dx.aux} values create higher numerical diffusion). In order to reduce numerical dispersion, one should
      increase the grid resolution (i.e. decrease \code{dx.aux}).
  }
}
\keyword{utilities}

