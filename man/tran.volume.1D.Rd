\name{tran.volume.1D}
\alias{tran.volume.1D}
\title{
  1-D Volumetric Advective-Diffusive Transport in an Aquatic System
}

\description{
  Estimates the volumetric transport term (i.e. the rate of change of the
  concentration due to diffusion and advection) in a one-dimensional model of
  an aquatic system (river, estuary).

  Volumetric transport implies the use of flows (mass per unit of time) rather
  than fluxes (mass per unit of area per unit of time) as is done in
  \code{\link{tran.1D}}.

  The \code{tran.volume.1D} routine is particularly suited for modelling
  channels (like rivers, estuaries) where the cross-sectional area changes,
  but where this area change needs not to be explicitly modelled as such.

  Another difference with \code{tran.1D} is that the present routine also
  allows lateral water or lateral mass input (as from side rivers or diffusive lateral ground water inflow).
}

\usage{
tran.volume.1D(C, C.up=C[1], C.down=C[length(C)],
   C.lat=C, F.up=NULL, F.down=NULL, F.lat=NULL,
   Disp,	flow = 0, flow.lat=NULL, AFDW = 1,
   V=NULL, full.check = FALSE, full.output = FALSE)

}

\arguments{
  \item{C }{tracer concentration, defined at the centre of the grid cells.
    A vector of length N [M/L3].
  }
  \item{C.up }{tracer concentration at the upstream interface. One value [M/L3].
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
  \item{F.up }{total tracer input at the upstream interface. One value [M/T].
  }
  \item{F.down }{total tracer input at downstream interface. One value [M/T].
  }
  \item{F.lat }{total lateral tracer input, defined at grid cell centres.
    One value, a vector of length N, or a 1D list property as defined by \code{\link{setup.prop.1D}},[M/T].
  }
  \item{Disp }{BULK dispersion coefficient, defined on grid cell interfaces.
    One value, a vector of length N+1, or a 1D list property as defined by \code{\link{setup.prop.1D}} [L3/T].
  }
  \item{flow }{water flow rate, defined on grid cell interfaces. One value, a vector of length N+1, or a list as defined by \code{\link{setup.prop.1D}} [L3/T]. If \code{flow.lat} is not \code{NULL} the \code{flow} should be one value containing the flow rate at the upstream boundary. If \code{flow.lat} is not \code{NULL} then \code{flow} must be a vector or a list.
  }
  \item{flow.lat }{lateral water flow rate [L3/T] into each volume box, defined at grid cell centres. One value, a vector of
    length N, or a list as defined by \code{\link{setup.prop.1D}}. If \code{flow.lat} has a value, then
    \code{flow} should be the flow rate at the upstream interface (one value).
    For each grid cell, the \code{flow} at the downstream side of a grid cell is 
    then estimated by water balance (adding \code{flow.lat} in the cell to 
    flow rate at the upstream side of the grid cell). If \code{flow.lat} is \code{NULL}, then it is determined by water balance 
    from \code{flow}.
  }
  \item{AFDW }{weight used in the finite difference scheme for advection,
    defined on grid cell interfaces; backward = 1, centred = 0.5, forward = 0;
    default is backward. One value, a vector of length N+1, or a
    list as defined by \code{\link{setup.prop.1D}} [-]. 
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
  \item{F }{mass flow across at the interface of each grid cell. A vector
    of length N+1 [M/T]. Only provided when (\code{full.output} = \code{TRUE}
  }
  \item{flow.up }{water flow across the upstream boundary, positive = INTO
    model domain. One value [L3/T]. Only provided when
    (\code{full.output} = \code{TRUE})
  }
  \item{flow.down }{water flow across the downstream boundary, positive = OUT
    of model domain. One value [L3/T]. Only provided when
    (\code{full.output} = \code{TRUE})
  }
  \item{F.lat }{lateral water input on each volume box, positive = INTO model
    domain. A vector of length N [L3/T]. Only provided when
    (\code{full.output} = \code{TRUE})
  }
  \item{F }{water flow across at the interface of each grid cell. A vector
    of length N+1 [M/T]. Only provided when (\code{full.output} = \code{TRUE})
  }
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
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

river.model <- function (t=0,OC,pars=NULL)
{
tran <- tran.volume.1D(C=OC,F.up=F.OC,F.lat=F.lat,Disp=Disp,
flow=flow.up,flow.lat=flow.lat,V=Volume)$dC 
reac <- - k*OC
return(list(dCdt = tran + reac))
}

#======================#
# Parameter definition #
#======================#

# Initialising morphology estuary: 

nbox          <- 500     # number of grid cells
lengthEstuary <- 100000  # length of estuary [m]
BoxLength     <- lengthEstuary/nbox # [m]
Distance      <- seq(BoxLength/2, by=BoxLength, len=nbox) # [m]
Int.Distance  <- seq(0, by=BoxLength, len=(nbox+1)) # [m]

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
F.lat <- rep(0,length.out=nbox)
flow.lat <- rep(0,length.out=nbox)
Conc1 <- steady.band(runif(nbox),fun=river.model,nspec=1)$y   
tran1 <- tran.volume.1D(C=Conc1,F.up=F.OC,F.lat=F.lat,Disp=Disp,
flow=flow.up,flow.lat=flow.lat,V=Volume,full.output=TRUE)

#scenario 1: with lateral input
F.lat <- F.lat.0*dnorm(x=Distance/lengthEstuary,
mean = Distance[nbox/2]/lengthEstuary, sd = 1/20, log = FALSE)/nbox 
flow.lat <- flow.lat.0*dnorm(x=Distance/lengthEstuary,
mean = Distance[nbox/2]/lengthEstuary, sd = 1/20, log = FALSE)/nbox 
Conc2 <- steady.band(runif(nbox),fun=river.model,nspec=1)$y   
tran2 <- tran.volume.1D(C=Conc2,F.up=F.OC,F.lat=F.lat,Disp=Disp,
flow=flow.up,flow.lat=flow.lat,V=Volume,full.output=TRUE)

#====================#
# Plotting output    #
#====================#
par(mfrow=c(2,1))

matplot(Distance/1000,cbind(Conc1,Conc2),lwd=2,
main="Organic carbon decay in the estuary",xlab="distance [km]",
ylab="OC Concentration [mM]",
type="l")  
legend ("topright",lty=1,col=c("black","red"),
        c("baseline","with side river input"))

matplot(Int.Distance/1000,cbind(tran1$flow,tran2$flow),lwd=2,
main="Longitudinal change in the water flow rate",xlab="distance [km]",
ylab="Flow rate [m3 s-1]", ylim=c(0,400),
type="l")  
legend ("bottomright",lty=1,col=c("black","red"),
        c("baseline","with side river input"))
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

