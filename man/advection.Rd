\name{advection.1D}
\alias{advection.1D}
\alias{advection.volume.1D}
\title{
  One-Dimensional Advection Equation
}

\description{
  Estimates the advection term in a one-dimensional model of a
  liquid (volume fraction constant and equal to one) or in a porous medium
  (volume fraction variable and lower than one).

  The interfaces between grid cells can have a variable cross-sectional area,
  e.g. when modelling spherical or cylindrical geometries (see example).
  
  TVD (total variation diminishing) slope limiters ensure monotonic and 
  positive schemes in the presence of strong gradients. 
  
  \code{advection.1-D}: uses finite differences.
   
  This implies the use of velocity (length per time) and fluxes 
  (mass per unit of area per unit of time). 

 
  \code{advection.volume.1D}
  Estimates the volumetric advection term in a one-dimensional model 
  of an aquatic system (river, estuary). This routine is particularly 
  suited for modelling channels (like rivers, estuaries) where the 
  cross-sectional area changes, and hence the velocity changes. 
  
  Volumetric transport implies the use of flows (mass per unit of time). 
  
  When solved dynamically, the euler method should be used, unless the
  first-order upstream method is used.  
}

\usage{
advection.1D(C, C.up = NULL, C.down = NULL,
  flux.up = NULL, flux.down = NULL, v, VF = 1, A = 1, dx, 
  dt.default = 1, adv.method = c("muscl", "super", "quick", "p3", "up"),
  full.check = FALSE)

advection.volume.1D(C, C.up = C[1], C.down = C[length(C)],
  F.up = NULL, F.down = NULL, flow, V, 
  dt.default = 1, adv.method = c("muscl", "super", "quick", "p3", "up"),
  full.check = FALSE)
}

\arguments{
  \item{C }{concentration, expressed per unit of phase volume, defined at the
    centre of each grid cell. A vector of length N [M/L3].
  }
  \item{C.up }{concentration at upstream boundary. One value [M/L3]. If \code{NULL},
    and \code{flux.up} is also \code{NULL}, then a zero-gradient boundary 
    is assumed, i.e. \code{C.up = C[1]}.
  }
  \item{C.down }{concentration at downstream boundary. One value [M/L3]. If \code{NULL},
    and \code{flux.down} is also \code{NULL}, then a zero-gradient boundary 
    is assumed, i.e. \code{C.down = C[length(C)]}.
  }
  \item{flux.up }{flux across the upstream boundary, positive = INTO model
    domain. One value, expressed per unit of total surface [M/L2/T]. 
    If \code{NULL}, the boundary is prescribed as
    a concentration boundary.
  }
  \item{flux.down }{flux across the downstream boundary, positive = OUT
    of model domain. One value, expressed per unit of total surface [M/L2/T].
    If \code{NULL}, the boundary is prescribed as
    a concentration boundary.
  }
  \item{F.up }{total input across the upstream boundary, positive = INTO model
    domain; used with \code{advection.volume.1D}. 
    One value, expressed in [M/T]. 
    If \code{NULL}, the boundary is prescribed as
    a concentration boundary.
  }
  \item{F.down }{total input across the downstream boundary, positive = OUT
    of model domain; used with \code{advection.volume.1D}. 
    One value, expressed in [M/T].
    If \code{NULL}, the boundary is prescribed as
    a concentration boundary.
  }
  \item{v }{advective velocity, defined on the grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length N+1 [L/T], or a \code{1D property} list; the list
    contains at least the element \code{int} (see \code{\link{setup.prop.1D}})
    [L/T]. Used with \code{advection.1D}.
    
  }
  \item{flow }{water flow rate, defined on grid cell interfaces. 
    One value, a vector of length N+1, or a list as defined by 
    \code{setup.prop.1D} [L^3/T]. 
    Used with \code{advection.volume.1D}.
  }
  \item{VF }{Volume fraction defined at the grid cell interfaces. One value,
    a vector of length N+1, or a \code{1D property} list; the list
    contains at least the elements \code{int} and \code{mid}
    (see \code{\link{setup.prop.1D}}) [-].
  }
  \item{A }{Interface area defined at the grid cell interfaces. One value,
    a vector of length N+1, or a \code{1D grid property} list; the list
    contains at least the elements \code{int} and \code{mid}
    (see \code{\link{setup.prop.1D}}) [L^2].
  }
  \item{dx }{distance between adjacent cell interfaces (thickness of grid
    cells). One value, a vector of length N, or a \code{1D grid} list containing
    at least the elements
    \code{dx} and \code{dx.aux} (see \code{\link{setup.grid.1D}}) [L].
  }
  \item{dt.default }{timestep to be used, if it cannot be estimated (e.g.
    when calculating steady-state conditions.
  }
  \item{V }{volume of cells. One value, or a vector of length N [L^3].
  }
  \item{adv.method }{the advection method, slope limiter used to reduce the 
    numerical dispersion. One of "quick","muscl","super","p3","up" - see details.
  }
  \item{full.check }{logical flag enabling a full check of the consistency
    of the arguments (default = \code{FALSE}; \code{TRUE} slows down execution
    by 50 percent).
  }
}

\value{
  \item{dC }{the rate of change of the concentration C due to advective 
    transport, defined in the centre of each grid cell. 
    The rate of change is expressed per unit of (phase) volume [M/L^3/T].
  }
  \item{adv.flux }{advective flux across at the interface of each grid cell.
    A vector of length N+1 [M/L2/T] - only for \code{advection.1D}.
  }
  \item{flux.up }{flux across the upstream boundary, positive = INTO model
    domain. One value [M/L2/T] - only for \code{advection.1D}.
  }
  \item{flux.down }{flux across the downstream boundary, positive = OUT of
    model domain. One value [M/L2/T] - only for \code{advection.1D}.
  }
  \item{adv.F }{advective mass flow across at the interface of each grid cell.
    A vector of length N+1 [M/T] - only for \code{advection.volume.1D}.
  }
  \item{F.up }{mass flow across the upstream boundary, positive = INTO model
    domain. One value [M/T] - only for \code{advection.volume.1D}.
  }
  \item{F.down }{flux across the downstream boundary, positive = OUT of
    model domain. One value [M/T] - only for \code{advection.volume.1D}.
  }
  \item{it }{number of split time iterations that were necessary.
  }  
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{

## =============================================================================
## EXAMPLE 1: Testing the various methods - moving a square pulse  
## use of advection.1D
## The tests as in Pietrzak 
## =============================================================================

#--------------------#
# Model formulation  #
#--------------------#
model <- function (t, y, parms,...) {

  adv <- advection.1D(y, v = v, dx = dx, 
     C.up = y[n], C.down = y[1], ...)  # out on one side -> in at other
  return(list(adv$dC))

}

#--------------------#
# Parameters         #
#--------------------#

n     <- 100
dx    <- 100/n
y     <- c(rep(1, 5), rep(2, 20), rep(1, n-25))
v     <- 2 
times <- 0:300   # 3 times out and in

#--------------------#
# model solution     #
#--------------------#

## a plotting function
plotfun <- function (Out, ...) {
  plot(Out[1, -1], type = "l", col = "red", ylab = "y", xlab = "x", ...)
  lines(Out[nrow(Out), 2:(1+n)])
}

# courant number = 2
pm   <- par(mfrow = c(2, 2))

## third order TVD, quickest limiter
out <- ode.1D(y = y, times = times, func = model, parms = 0, hini = 1,
              method = "euler", nspec = 1, adv.method = "quick")

plotfun(out, main = "quickest, euler") 

## third-order ustream-biased polynomial
out2 <- ode.1D(y = y, times = times, func = model, parms = 0, hini = 1,
              method = "euler", nspec = 1, adv.method = "p3")

plotfun(out2, main = "p3, euler")

## third order TVD, superbee limiter
out3 <- ode.1D(y = y, times = times, func = model, parms = 0, hini = 1,
              method = "euler", nspec = 1, adv.method = "super")

plotfun(out3, main = "superbee, euler")

## third order TVD, muscl limiter
out4 <- ode.1D(y = y, times = times, func = model, parms = 0, hini = 1,
              method = "euler", nspec = 1, adv.method = "muscl")

plotfun(out4, main = "muscl, euler")

## =============================================================================
## upstream, different time-steps , i.e. different courant number
## =============================================================================
out <- ode.1D(y = y, times = times, func = model, parms = 0, 
              method = "euler", nspec = 1, adv.method = "up")

plotfun(out, main = "upstream, courant number = 2")

out2 <- ode.1D(y = y, times = times, func = model, parms = 0, hini = 0.5,  
               method = "euler", nspec = 1, adv.method = "up")
plotfun(out2, main = "upstream, courant number = 1")


## Now muscl scheme, velocity against x-axis
y     <- rev(c(rep(0, 5), rep(1, 20), rep(0, n-25)))
v     <- -2.0
out6 <- ode.1D(y = y, times = times, func = model, parms = 0, hini = 1,
               method = "euler", nspec = 1, adv.method = "muscl")

plotfun(out6, main = "muscl, reversed velocity, , courant number = 1")

image(out6, mfrow = NULL)
par(mfrow = pm)


## =============================================================================
## EXAMPLE 2: moving a square pulse in a widening river  
## use of advection.volume.1D
## =============================================================================

#--------------------#
# Model formulation  #
#--------------------#

river.model <- function (t=0, C, pars = NULL, ...)  {

 tran <- advection.volume.1D(C = C, C.up = 0,
                     flow = flow, V = Volume,...)
 return(list(dCdt = tran$dC, F.down = tran$F.down, F.up = tran$F.up))
}

#--------------------#
# Parameters         #
#--------------------#

# Initialising morphology river:

nbox          <- 100                # number of grid cells
lengthRiver   <- 100000             # [m]
BoxLength     <- lengthRiver / nbox # [m]

Distance      <- seq(BoxLength/2, by = BoxLength, len = nbox)   # [m]

# Cross sectional area: sigmoid function of distance [m2]
CrossArea <- 4000 + 72000 * Distance^5 /(Distance^5+50000^5)

# Volume of boxes                          (m3)
Volume    <- CrossArea*BoxLength

# Transport coefficients
flow      <- 1000*24*3600   # m3/d, main river upstream inflow

#--------------------#
# Model solution     #
#--------------------#

pm   <- par(mfrow=c(2,2))

# a square pulse
yini <-  c(rep(10, 10), rep(0, nbox-10))

## third order TVD, muscl limiter
Conc <- ode.1D(y = yini, fun = river.model,  method = "euler", hini = 1,
              parms = NULL, nspec = 1, times = 0:40, adv.method = "muscl")

image(Conc, main = "muscl", mfrow = NULL)
plot(Conc[30, 2:(1+nbox)], type = "l", lwd = 2, xlab = "x", ylab = "C", 
     main = "muscl after 30 days")

## simple upstream differencing
Conc2<- ode.1D(y = yini, fun = river.model, method = "euler", hini = 1,
              parms = NULL, nspec = 1, times = 0:40, adv.method = "up")

image(Conc2, main = "upstream", mfrow = NULL)
plot(Conc2[30, 2:(1+nbox)], type = "l", lwd = 2, xlab = "x", ylab = "C", 
     main = "upstream after 30 days")

par(mfrow = pm)


# Note: the more sophisticated the scheme, the more mass lost/created
# increase tolerances to improve this.

CC <- Conc[ , 2:(1+nbox)]
MASS <- t(CC)*Volume     
colSums(MASS)

## =============================================================================
## EXAMPLE 3: A steady-state solution 
## use of advection.volume.1D
## =============================================================================

Sink <- function (t, y, parms, ...) {
  C1 <- y[1:N]
  C2 <- y[(N+1):(2*N)]
  C3 <- y[(2*N+1):(3*N)]
  # Rate of change= Flux gradient and first-order consumption                  
  
  # upstream can be implemented in two ways:
 
  dC1  <- advection.1D(C1, v = sink, dx = dx, 
           C.up = 100, adv.method = "up", ...)$dC - decay*C1

# same, using tran.1D  
#  dC1  <- tran.1D(C1, v = sink, dx = dx, 
#        C.up = 100, D = 0)$dC -
#           decay*C1

  dC2  <- advection.1D(C2, v = sink, dx = dx, 
        C.up = 100, adv.method = "p3", ...)$dC -
           decay*C2

  dC3  <- advection.1D(C3, v = sink, dx = dx, 
        C.up = 100, adv.method = "muscl", ...)$dC -
           decay*C3

  list(c(dC1, dC2, dC3))
}

N     <- 10
L     <- 1000
dx    <- L/N                          # thickness of boxes
sink  <- 10
decay <- 0.1
out <- steady.1D(runif(3*N), func = Sink, names = c("C1", "C2", "C3"),
        parms = NULL, nspec = 3, bandwidth = 2)
matplot(out$y, 1:N, type = "l", ylim = c(10, 0), lwd = 2, 
  main = "Steady-state")          
legend("bottomright", col = 1:3, lty = 1:3, 
  c("upstream", "p3", "muscl"))


}
\references{
Pietrzak J (1998) The use of TVD limiters for forward-in-time
  upstream-biased advection schemes in ocean modeling. Monthly
  Weather Review 126: 812 .. 830

Hundsdorfer W and Verwer JG (2003)
 Numerical Solution of Time-Dependent Advection-Diffusion-Reaction Equations. 
 Springer Series in Computational Mathematics, Springer-Verlag, Berlin, 471 pages
  
Burchard H, Bolding K, Villarreal MR (1999) GOTM, a general
ocean turbulence model. Theory, applications and test cases.
Tech Rep EUR 18745 EN. European Commission  

Leonard BP (1988) Simple high accuracy resolution program for convective modeling
       of discontinuities. Int. J. Numer. Meth.Fluids 8: 1291--1318.

Roe PL (1985) Some contributions to the modeling of discontinuous flows.
       Lect. Notes Appl. Math. 22: 163-193.

van Leer B. (1979) Towards the ultimate conservative difference scheme V. A second 
               order sequel to Godunov's method. J. Comput. Phys. 32: 101-136
}

\details{
  This implementation is based on the GOTM code

  The \bold{boundary conditions} are either
  \itemize{
    \item zero-gradient.
    \item fixed concentration.
    \item fixed flux.
  }
  The above order also shows the priority. The default condition is the
  zero gradient. The fixed concentration condition overrules the zero gradient.
  The fixed flux overrules the other specifications.

  Ensure that the boundary conditions are well defined: for instance, it 
  does not make sense to specify an influx in a boundary cell with the advection
  velocity pointing outward.
    
  \bold{Transport properties:}
  
  The \emph{advective velocity} (\code{v}),
  the \emph{volume fraction} (VF), and the \emph{interface surface} (\code{A}),
  can either be specified as one value, a vector, or a 1D property list
  as generated by \code{\link{setup.prop.1D}}.

  When a vector, this vector must be of length N+1, defined at all grid
  cell interfaces, including the upper and lower boundary.

  The \bold{finite difference grid} (\code{dx}) is specified either as
  one value, a vector or a 1D grid list, as generated by \code{\link{setup.grid.1D}}. 
  
  Several slope limiters are implemented to obtain monotonic and positive 
  schemes also in the presence of strong gradients, i.e. to reduce the effect
  of numerical dispersion. The methods are (Pietrzak, 1989, Hundsdorfer and 
  Verwer, 2003):
  \itemize{
    \item "quick": third-order scheme (TVD) with ULTIMATE QUICKEST limiter 
     (quadratic upstream interpolation for convective kinematics with 
     estimated stream terms) (Leonard, 1988)
    \item "muscl": third-order scheme (TVD) with MUSCL limiter (monotonic upstream 
      centered schemes for conservation laws) (van Leer, 1979).
    \item "super": third-order scheme (TVD) with Superbee limiter (method=Superbee)
      (Roe, 1985)
    \item "p3": third-order upstream-biased polynomial scheme (method=P3)
    \item "up": first-order upstream ( method=UPSTREAM) - this is the same
      method as implemented in \link{tran.1D} or \link{tran.volume.1D}
  }
  where "TVD" means a total variation diminishing scheme

  Some schemes may produce artificial steepening. Scheme "p3" is not necessarily
  monotone (may produce negative concentrations!).

  If during a certain time step the maximum Courant number is larger
  than one, a split iteration will be carried out which guarantees that the
  split step Courant numbers are just smaller than 1. The maximal number of such
  iterations is set to 100.
  
  These limiters are supposed to work with explicit methods (\link{euler}). However, 
  they will also work with implicit methods, although less effectively.
  Integrate \link{ode.1D} only if the model is stiff (see first example).
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
\keyword{utilities}

\seealso{
  \code{\link{tran.1D}}, for a discretisation of the general transport equations
}
