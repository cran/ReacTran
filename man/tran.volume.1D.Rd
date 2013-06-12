\name{tran.volume.1D}
\alias{tran.volume.1D}
\alias{tran.volume.2D}
\alias{tran.volume.3D}
\title{
  1-D, 2-D and 3-D Volumetric Advective-Diffusive Transport in an Aquatic System
}

\description{
  Estimates the volumetric transport term (i.e. the rate of change of the
  concentration due to diffusion and advection) in a 1-D, 2-D or 3-D model of
  an aquatic system (river, estuary).

  Volumetric transport implies the use of flows (mass per unit of time) rather
  than fluxes (mass per unit of area per unit of time) as is done in
  \code{\link{tran.1D}}, \code{\link{tran.2D}} or \code{\link{tran.3D}}.

  The \code{tran.volume.xD} routines are particularly suited for modelling
  channels (like rivers, estuaries) where the cross-sectional area changes,
  but where this area change needs not to be explicitly modelled as such.

  Another difference with \code{tran.1D} is that the \code{tran.volume.1D} 
  routine also allows lateral water or lateral mass input (as from side rivers 
  or diffusive lateral ground water inflow).
  
  The \code{tran.volume.2D} routine can check for water balance and assume an
  in- or efflux in case the net flows in and out of a box are not = 0
}

\usage{
tran.volume.1D(C, C.up = C[1], C.down = C[length(C)],
               C.lat = C, F.up = NULL, F.down = NULL, F.lat = NULL,
               Disp,	flow = 0, flow.lat = NULL, AFDW = 1,
               V = NULL, full.check = FALSE, full.output = FALSE)

tran.volume.2D(C, C.x.up = C[1, ], C.x.down = C[nrow(C), ], 
               C.y.up = C[, 1], C.y.down = C[, ncol(C)], 
               C.z = C, masscons = TRUE, 
               F.x.up = NULL, F.x.down = NULL, 
               F.y.up = NULL, F.y.down = NULL, 
               Disp.grid = NULL, Disp.x = NULL, Disp.y = Disp.x, 
               flow.grid = NULL, flow.x = NULL, flow.y = NULL, 
               AFDW.grid = NULL, AFDW.x = 1, AFDW.y = AFDW.x,         
               V = NULL, full.check = FALSE, full.output = FALSE) 

tran.volume.3D(C,  C.x.up = C[1, , ], C.x.down = C[dim(C)[1], , ],                
               C.y.up = C[, 1, ], C.y.down = C[, dim(C)[2], ], 
               C.z.up = C[, , 1], C.z.down = C[, , dim(C)[3]], 
               F.x.up = NULL, F.x.down = NULL,  
               F.y.up = NULL, F.y.down = NULL, 
               F.z.up = NULL, F.z.down = NULL,         
               Disp.grid = NULL, 
               Disp.x = NULL, Disp.y = Disp.x, Disp.z = Disp.x,      
               flow.grid = NULL, flow.x = 0, flow.y = 0, flow.z = 0, 
               AFDW.grid = NULL, AFDW.x = 1, AFDW.y = AFDW.x, 
               AFDW.z = AFDW.x,                 
               V = NULL, full.check = FALSE, full.output = FALSE) 
}

\arguments{
  \item{C }{tracer concentration, defined at the centre of the grid cells.
    A vector of length N [M/L3] (tran.volume.1D), 
    a matrix of dimension Nr*Nc (tran.volume.2D) or 
    an Nx*Ny*Nz array (tran.volume.3D) [M/L3].
  }
  \item{C.up }{tracer concentration at the upstream interface. 
   One value [M/L3].
  }
  \item{C.down }{tracer concentration at downstream interface. One value [M/L3].
  }
  \item{C.lat }{tracer concentration in the lateral input, defined at
    grid cell centres. One value, a vector of length N, or a
    list as defined by \code{\link{setup.prop.1D}} [M/L3].
    The default is  \code{C.lat = C}, (a zero-gradient condition).
    Setting \code{C.lat=0}, together with a positive \code{F.lat} will
    lead to dilution of the tracer concentration in the
    grid cells.
  }
  \item{C.x.up }{concentration at upstream boundary in x-direction;
    vector of length Ny (2D) or matrix of dimensions Ny*Nz (3D) [M/L3].
  }
  \item{C.x.down }{concentration at downstream boundary in x-direction;
    vector of length Ny (2D) or matrix of dimensions Ny*Nz (3D) [M/L3].
  }
  \item{C.y.up }{concentration at upstream boundary in y-direction;
    vector of length Nx (2D) or matrix of dimensions Nx*Nz (3D) [M/L3].
  }
  \item{C.y.down }{concentration at downstream boundary in y-direction;
   vector of length Nx (2D) or matrix of dimensions Nx*Nz (3D) [M/L3].
  }
  \item{C.z.up }{concentration at upstream boundary in z-direction;
    matrix of dimensions Nx*Ny [M/L3].
  }
  \item{C.z.down }{concentration at downstream boundary in z-direction;
    matrix of dimensions Nx*Ny [M/L3].
  }
  \item{C.z }{concentration at boundary in z-direction for 2-D models where
      \code{masscons} = TRUE. Matrix of dimensions Nx*Ny [M/L3].
  }
  \item{masscons }{When \code{TRUE}, will check flow balance in 2D model. 
    The flow in the third direction will then be estimated.
  }
  \item{F.up }{total tracer input at the upstream interface. One value [M/T].
  }
  \item{F.down }{total tracer input at downstream interface. One value [M/T].
  }
  \item{F.lat }{total lateral tracer input, defined at grid cell centres.
    One value, a vector of length N, or a 1D list property as defined by \code{\link{setup.prop.1D}},[M/T].
  }
  \item{F.x.up }{total tracer input at the upstream interface in x-direction.
  positive = INTO model domain. A vector of length Ny (2D) or a matrix of dimensions Ny*Nz (3D)  [M/T].
  }
  \item{F.x.down }{total tracer input at downstream interface in x-direction. positive = INTO model domain. A vector of length Ny (2D) or a matrix of dimensions Ny*Nz (3D)  [M/T].
  }
  \item{F.y.up }{total tracer input at the upstream interface in y-direction.
  positive = INTO model domain. A vector of length Nx (2D) or a matrix of dimensions Nx*Nz (3D)  [M/T].
  }
  \item{F.y.down }{total tracer input at downstream interface in y-direction. positive = INTO model domain. A vector of length Nx (2D) or a matrix of dimensions Nx*Nz (3D)  [M/T].
  }
  \item{F.z.up }{total tracer input at the upstream interface in z-direction.
  positive = INTO model domain. A matrix of dimensions Nx*Ny    [M/T].
  }
  \item{F.z.down }{total tracer input at downstream interface in z-direction. positive = INTO model domain. A matrix of dimensions Nx*Ny [M/T].
  }
  \item{Disp.grid }{BULK dispersion coefficients defined on all grid cell
    interfaces. For \code{tran.volume.2D}, should contain two matrices, x.int (dimension (Nx+1)*Ny)  and y.int (dimension Nx * (Ny+1)). 
    For tran.volume.3D should contain three arrays x.int (dim = (Nx+1)*Ny*Nz), y.int (dim = Nx*(Ny+1)*Nz), and z.int (dim = Nx*Ny*(Nz+1))
  }
  \item{Disp }{BULK dispersion coefficient, defined on grid cell interfaces.
    One value, a vector of length N+1, or a 1D list property as defined by \code{\link{setup.prop.1D}} [L3/T].
  }
  \item{Disp.x }{BULK dispersion coefficient in x-direction, defined on grid cell interfaces. One value, a vector of length (Nx+1), a \code{prop.1D} list created by \code{\link{setup.prop.1D}}, a (Nx+1)* Ny matrix (2D) or a  Nx*(Ny+1)*Nz array (3D)   [L3/T].
  }
  \item{Disp.y }{BULK dispersion coefficient in y-direction, defined on grid cell
    interfaces. One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*(Ny+1) matrix (2D) or a Nx*(Ny+1)*Nz array  (3D)[L3/T].
  }
  \item{Disp.z }{BULK dispersion coefficient in z-direction, defined on grid cell
    interfaces. One value, a vector of length (Nz+1), or a Nx*Ny*(Nz+1) array   [L3/T].
  }
  \item{flow }{water flow rate, defined on grid cell interfaces. One value, a vector of length N+1, or a list as defined by \code{\link{setup.prop.1D}} [L3/T]. 
   If \code{flow.lat} is not \code{NULL} the \code{flow} should be one value containing the flow rate at the upstream boundary. 
   If \code{flow.lat} is not \code{NULL} then \code{flow} must be a vector or a list.
  }
  \item{flow.lat }{lateral water flow rate [L3/T] into each volume box, defined at grid cell centres. One value, a vector of
    length N, or a list as defined by \code{\link{setup.prop.1D}}. If \code{flow.lat} has a value, then
    \code{flow} should be the flow rate at the upstream interface (one value).
    For each grid cell, the \code{flow} at the downstream side of a grid cell is 
    then estimated by water balance (adding \code{flow.lat} in the cell to 
    flow rate at the upstream side of the grid cell). If \code{flow.lat} is \code{NULL}, then it is determined by water balance 
    from \code{flow}.
  }
  \item{flow.grid }{flow rates defined on all grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    Should contain elements x.int, y.int, z.int (3-D), arrays with the values on the
    interfaces in x, y and z-direction  [L3/T].
  }
  \item{flow.x }{flow rates in the x-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}} (2D),
    a (Nx+1)*Ny matrix (2D) or a (Nx+1)*Ny*Nz array (3D) [L3/T].
  }
  \item{flow.y }{flow rates in the y-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}} (2D),
    a Nx*(Ny+1) matrix (2D) or a Nx*(Ny+1)*Nz array [L3/T].
  }
  \item{flow.z }{flow rates in the z-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length (Nz+1),
    or a Nx*Ny*(Nz+1) array [L3/T].
  }
  \item{AFDW }{weight used in the finite difference scheme for advection,
    defined on grid cell interfaces; backward = 1, centred = 0.5, forward = 0;
    default is backward. One value, a vector of length N+1, or a
    list as defined by \code{\link{setup.prop.1D}} [-]. 
  }
  \item{AFDW.grid }{weight used in the finite difference scheme for advection
    in the x-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    For \code{tran.volume.3D} should contain elements x.int, y.int, z.int (3D),  for \code{tran.volume.2D} should contain elements x.int and y.int. [-].
  }
  \item{AFDW.x }{weight used in the finite difference scheme for advection
    in the x-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
     a (Nx+1)*Ny matrix (2D) or a (Nx+1)*Ny*Nz array (3D) [-].
  }
  \item{AFDW.y }{weight used in the finite difference scheme for advection
    in the y-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    a Nx*(Ny+1) matrix (2D) or a Nx*(Ny+1)*Nz array [-].
  }
  \item{AFDW.z }{weight used in the finite difference scheme for advection
    in the z-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    One value, a vector of length (Nz+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*Ny*(Nz+1) array [-].
  }
  \item{V }{grid cell volume, defined at grid cell centres [L3]. One value, a
    vector of length N, or a list as defined by \code{\link{setup.prop.1D}}.
  }
  \item{full.check }{logical flag enabling a full check of the consistency
    of the arguments (default = \code{FALSE}; \code{TRUE} slows down execution
    by 50 percent).
  }
  \item{full.output }{logical flag enabling a full return of the output
    (default = \code{FALSE}; \code{TRUE} slows down execution by 20 percent).
  }
}

\value{
  \item{dC }{the rate of change of the concentration C due to transport,
    defined in the centre of each grid cell [M/L3/T].
  }
  \item{F.up }{mass flow across the upstream boundary, positive = INTO
    model domain. One value [M/T].
  }
  \item{F.down }{mass flow across the downstream boundary, positive = OUT
    of model domain. One value [M/T].
  }
  \item{F.lat }{lateral mass input per volume box, positive = INTO model
    domain. A vector of length N [M/T].
  }
  \item{flow }{water flow across the interface of each grid cell. A vector
    of length N+1 [L3/T]. Only provided when (\code{full.output} = \code{TRUE}
  }
  \item{flow.up }{water flow across the upstream (external) boundary, positive = INTO
    model domain. One value [L3/T]. Only provided when (\code{full.output} = \code{TRUE})
  }
  \item{flow.down }{water flow across the downstream (external) boundary, positive = OUT
    of model domain. One value [L3/T]. Only provided when
    (\code{full.output} = \code{TRUE})
  }
  \item{flow.lat }{lateral water input on each volume box, positive = INTO model
    domain. A vector of length N [L3/T]. Only provided when
    (\code{full.output} = \code{TRUE})
  }
  \item{F }{mass flow across at the interface of each grid cell. A vector
    of length N+1 [M/T]. Only provided when (\code{full.output} = \code{TRUE})
  }
}

\author{
  Filip Meysman <filip.meysman@nioz.nl>,
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{

## =============================================================================
##  EXAMPLE : organic carbon (OC) decay in a widening estuary
## =============================================================================

# Two scenarios are simulated: the baseline includes only input 
# of organic matter upstream. The second scenario simulates the 
# input of an important side river half way the estuary.  

#====================#
# Model formulation  #
#====================#

river.model <- function (t = 0, OC, pars = NULL) {

  tran <- tran.volume.1D(C = OC, F.up = F.OC, F.lat = F.lat,
          Disp = Disp, flow = flow.up, flow.lat = flow.lat, 
          V = Volume, full.output = TRUE) 

  reac <- - k*OC
  return(list(dCdt = tran$dC + reac, Flow = tran$flow))
}

#======================#
# Parameter definition #
#======================#

# Initialising morphology estuary: 

nbox          <- 500     # number of grid cells
lengthEstuary <- 100000  # length of estuary [m]
BoxLength     <- lengthEstuary/nbox # [m]
Distance      <- seq(BoxLength/2, by = BoxLength, len =nbox) # [m]
Int.Distance  <- seq(0, by = BoxLength, len = (nbox+1))      # [m]

# Cross sectional area: sigmoid function of estuarine distance [m2]
CrossArea <- 4000 + 72000 * Distance^5 /(Distance^5+50000^5)

# Volume of boxes                          (m3)
Volume  <- CrossArea*BoxLength

# Transport coefficients
Disp    <- 1000   # m3/s, bulk dispersion coefficient
flow.up  <- 180    # m3/s, main river upstream inflow
flow.lat.0  <- 180    # m3/s, side river inflow

F.OC    <- 180               # input organic carbon [mol s-1]
F.lat.0 <- 180              # lateral input organic carbon [mol s-1]

k       <- 10/(365*24*3600)  # decay constant organic carbon [s-1]


#====================#
# Model solution     #
#====================#
#scenario 1: without lateral input
F.lat    <- rep(0, length.out = nbox)
flow.lat <- rep(0, length.out = nbox)

Conc1 <- steady.1D(runif(nbox), fun = river.model, nspec = 1, name = "OC")   

#scenario 2: with lateral input
F.lat <- F.lat.0 * dnorm(x =Distance/lengthEstuary,
                         mean = Distance[nbox/2]/lengthEstuary, 
                         sd = 1/20, log = FALSE)/nbox 
flow.lat <- flow.lat.0 * dnorm(x = Distance/lengthEstuary,
                               mean = Distance[nbox/2]/lengthEstuary, 
                               sd = 1/20, log = FALSE)/nbox 

Conc2 <- steady.1D(runif(nbox), fun = river.model, nspec = 1, name = "OC")   

#====================#
# Plotting output    #
#====================#
# use S3 plot method
plot(Conc1, Conc2, grid = Distance/1000, which = "OC", 
     mfrow = c(2, 1), lwd = 2, xlab = "distance [km]", 
     main = "Organic carbon decay in the estuary",
     ylab = "OC Concentration [mM]")
       
plot(Conc1, Conc2, grid = Int.Distance/1000, which = "Flow", 
     mfrow = NULL, lwd = 2, xlab = "distance [km]", 
     main = "Longitudinal change in the water flow rate",
     ylab = "Flow rate [m3 s-1]")  

legend ("topright", lty = 1:2, col = 1:2, lwd = 2,
        c("baseline", "+ side river input"))

}
\references{
  Soetaert and Herman (2009) A practical guide to ecological modelling -
  using R as a simulation platform. Springer.
}

\details{
  The \bold{boundary conditions} are of type
  \itemize{
    \item 1. zero-gradient (default)
    \item 2. fixed concentration
    \item 3. fixed input
  }
  The \emph{bulk dispersion coefficient} (\code{Disp}) and the \emph{flow rate}
  (\code{flow}) can be either one value or a vector of length N+1, defined at
  all grid cell interfaces, including upstream and downstream boundary.

  The spatial discretisation is given by the volume of each box (\code{V}),
  which can be one value or a vector of length N+1, defined at the centre of
  each grid cell.

  The water flow is mass conservative. Over each volume box, the routine
  calculates internally the downstream outflow of water in terms of the
  upstream inflow and the lateral inflow.
}

\seealso{
  \code{\link{tran.1D}}

  \code{\link{advection.volume.1D}}, for more sophisticated advection schemes

}

\keyword{utilities}

