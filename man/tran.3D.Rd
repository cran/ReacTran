\name{tran.3D}
\alias{tran.3D}
\title{
  General Three-Dimensional Advective-Diffusive Transport
}

\description{
  Estimates the transport term (i.e. the rate of change of a concentration
  due to diffusion and advection) in a three-dimensional rectangular
  model domain.
  
  Do not use with too many boxes!
}

\usage{
tran.3D (C, C.x.up=C[1,,], C.x.down=C[dim(C)[1],,],
  C.y.up=C[,1,],  C.y.down=C[,dim(C)[2],],
  C.z.up=C[,,1],  C.z.down=C[,,dim(C)[3]],
  flux.x.up=NULL, flux.x.down=NULL,
  flux.y.up=NULL, flux.y.down=NULL,
  flux.z.up=NULL, flux.z.down=NULL,
  a.bl.x.up=NULL, a.bl.x.down=NULL, 
  a.bl.y.up=NULL, a.bl.y.down=NULL, 
  a.bl.z.up=NULL, a.bl.z.down=NULL, 
  D.grid=NULL, D.x=NULL, D.y=D.x, D.z=D.x,
  v.grid=NULL, v.x=0, v.y=0, v.z=0,
  AFDW.grid=NULL, AFDW.x=1, AFDW.y=AFDW.x, AFDW.z=AFDW.x,
  VF.grid=NULL, VF.x=1, VF.y=VF.x, VF.z=VF.x,
  A.grid=NULL, A.x=1, A.y=1, A.z=1,
  grid=NULL, dx=NULL, dy=NULL, dz=NULL,
  full.check = FALSE, full.output = FALSE)
}

\arguments{
  \item{C }{concentration, expressed per unit volume, defined at the centre
    of each grid cell; Nx*Ny*Nz matrix [M/L3].
  }
  \item{C.x.up }{concentration at upstream boundary in x-direction;
    matrix of dimensions Ny*Nz [M/L3].
  }
  \item{C.x.down }{concentration at downstream boundary in x-direction;
    matrix of dimensions Ny*Nz [M/L3].
  }
  \item{C.y.up }{concentration at upstream boundary in y-direction;
    matrix of dimensions Nx*Nz [M/L3].
  }
  \item{C.y.down }{concentration at downstream boundary in y-direction;
   matrix of dimensions Nx*Nz [M/L3].
  }
  \item{C.z.up }{concentration at upstream boundary in z-direction;
    matrix of dimensions Nx*Ny [M/L3].
  }
  \item{C.z.down }{concentration at downstream boundary in z-direction;
    matrix of dimensions Nx*Ny [M/L3].
  }
  \item{flux.x.up }{flux across the upstream boundary in x-direction,
    positive = INTO model domain; matrix of dimensions Ny*Nz [M/L2/T].
  }
  \item{flux.x.down }{flux across the downstream boundary in x-direction,
    positive = OUT of model domain; matrix of dimensions Ny*Nz [M/L2/T].
  }
  \item{flux.y.up }{flux across the upstream boundary in y-direction,
    positive = INTO model domain; matrix of dimensions Nx*Nz [M/L2/T].
  }
  \item{flux.y.down }{flux across the downstream boundary in y-direction,
    positive = OUT of model domain; matrix of dimensions Nx*Nz [M/L2/T].
  }
  \item{flux.z.up }{flux across the upstream boundary in z-direction,
    positive = INTO model domain; matrix of dimensions Nx*Ny [M/L2/T].
  }
  \item{flux.z.down }{flux across the downstream boundary in z-direction,
    positive = OUT of model domain; matrix of dimensions Nx*Ny [M/L2/T].
  }
  \item{a.bl.x.up }{transfer coefficient across the upstream boundary layer.
     in x-direction

     \code{Flux=a.bl.x.up*(C.x.up-C[1,,])}. One value [L/T].
  }
  \item{a.bl.x.down }{transfer coefficient across the downstream boundary
    layer in x-direction;

    \code{Flux=a.bl.x.down*(C[Nx,,]-C.x.down)}.
    One value [L/T].
  }
  \item{a.bl.y.up }{transfer coefficient across the upstream boundary layer.
     in y-direction

     \code{Flux=a.bl.y.up*(C.y.up-C[,1,])}. One value [L/T].
  }
  \item{a.bl.y.down }{transfer coefficient across the downstream boundary
    layer in y-direction;

    \code{Flux=a.bl.y.down*(C[,Ny,]-C.y.down)}.
    One value [L/T].
  }
  \item{a.bl.z.up }{transfer coefficient across the upstream boundary layer.
     in y-direction

     \code{Flux=a.bl.y.up*(C.y.up-C[,,1])}. One value [L/T].
  }
  \item{a.bl.z.down }{transfer coefficient across the downstream boundary
    layer in z-direction;

    \code{Flux=a.bl.z.down*(C[,,Nz]-C.z.down)}.
    One value [L/T].
  }
  \item{D.grid }{diffusion coefficient defined on all grid cell
    interfaces.
    Should contain elements x.int, y.int, z.int, arrays with the values on the
    interfaces in x, y and z-direction, and with dimensions
    (Nx+1)*Ny*Nz, Nx*(Ny+1)*Nz and Nx*Ny*(Nz+1) respectively. [L2/T].
  }
  \item{D.x }{diffusion coefficient in x-direction, defined on grid cell
    interfaces. One value, a vector of length (Nx+1),
    or a (Nx+1)* Ny *Nz array [L2/T].
  }
  \item{D.y }{diffusion coefficient in y-direction, defined on grid cell
    interfaces. One value, a vector of length (Ny+1),
    or a Nx*(Ny+1)*Nz array [L2/T].
  }
  \item{D.z }{diffusion coefficient in z-direction, defined on grid cell
    interfaces. One value, a vector of length (Nz+1),
    or a Nx*Ny*(Nz+1) array [L2/T].
  }
  \item{v.grid }{advective velocity defined on all grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    Should contain elements x.int, y.int, z.int, arrays with the values on the
    interfaces in x, y and z-direction, and with dimensions
    (Nx+1)*Ny*Nz, Nx*(Ny+1)*Nz and Nx*Ny*(Nz+1) respectively.  [L/T].
  }
  \item{v.x }{advective velocity in the x-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length (Nx+1),
    or a (Nx+1)*Ny*Nz array [L/T].
  }
  \item{v.y }{advective velocity in the y-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length (Ny+1),
    or a Nx*(Ny+1)*Nz array [L/T].
  }
  \item{v.z }{advective velocity in the z-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value, a vector of length (Nz+1),
    or a Nx*Ny*(Nz+1) array [L/T].
  }
  \item{AFDW.grid }{weight used in the finite difference scheme for advection
    in the x-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    Should contain elements x.int, y.int, z.int, arrays with the values on the
    interfaces in x, y and z-direction, and with dimensions
    (Nx+1)*Ny*Nz, Nx*(Ny+1)*Nz and Nx*Ny*(Nz+1) respectively. [-].
  }
  \item{AFDW.x }{weight used in the finite difference scheme for advection
    in the x-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nx+1)*Ny*Nz array [-].
  }
  \item{AFDW.y }{weight used in the finite difference scheme for advection
    in the y-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*(Ny+1)*Nz array [-].
  }
  \item{AFDW.z }{weight used in the finite difference scheme for advection
    in the z-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward.
    One value, a vector of length (Nz+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*Ny*(Nz+1) array [-].
  }
  \item{VF.grid }{Volume fraction. A list.
    Should contain elements x.int, y.int, z.int, arrays with the values on the
    interfaces in x, y and z-direction, and with dimensions
    (Nx+1)*Ny*Nz, Nx*(Ny+1)*Nz and Nx*Ny*(Nz+1) respectively. [-].
  }
  \item{VF.x }{Volume fraction at the grid cell interfaces in the x-direction.
    One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nx+1)*Ny*Nz array [-].
  }
  \item{VF.y }{Volume fraction at the grid cell interfaces in the y-direction.
    One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*(Ny+1)*Nz array [-].
  }
  \item{VF.z }{Volume fraction at the grid cell interfaces in the z-direction.
    One value, a vector of length (Nz+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*Ny*(Nz+1) array [-].
  }
  \item{A.grid }{Interface area, a list.
    Should contain elements x.int, y.int, z.int, arrays with the values on the
    interfaces in x, y and z-direction, and with dimensions
    (Nx+1)*Ny*Nz, Nx*(Ny+1)*Nz and Nx*Ny*(Nz+1) respectively. [L2].
  }
  \item{A.x }{Interface area defined at the grid cell interfaces in
    the x-direction. One value, a vector of length (Nx+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nx+1)*Ny*Nz array [L2].
  }
  \item{A.y }{Interface area defined at the grid cell interfaces in
    the y-direction. One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*(Ny+1)*Nz array [L2].
  }
  \item{A.z }{Interface area defined at the grid cell interfaces in
    the z-direction. One value, a vector of length (Nz+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nx*Ny*(Nz+1) array [L2].
  }
  \item{dx }{distance between adjacent cell interfaces in the x-direction
    (thickness of grid cells). One value or vector of length Nx [L].
  }
  \item{dy }{distance between adjacent cell interfaces in the y-direction
    (thickness of grid cells). One value or vector of length Ny [L].
  }
  \item{dz }{distance between adjacent cell interfaces in the z-direction
    (thickness of grid cells). One value or vector of length Nz [L].
  }
  \item{grid }{discretization grid, a list containing at least elements
    \code{dx}, \code{dx.aux}, \code{dy}, \code{dy.aux}, \code{dz}, \code{dz.aux}
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
    defined in the centre of each grid cell, an array with dimension
    Nx*Ny*Nz [M/L3/T].
  }
  \item{C.x.up }{concentration at the upstream interface in x-direction.
    A matrix of dimension Ny*Nz [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{C.x.down }{concentration at the downstream interface in x-direction.
   A matrix of dimension Ny*Nz [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{C.y.up }{concentration at the upstream interface in y-direction.
    A matrix of dimension Nx*Nz [M/L3].  Only when \code{full.output = TRUE}.
  }
  \item{C.y.down }{concentration at the downstream interface in y-direction.
    A matrix of dimension Nx*Nz [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{C.z.up }{concentration at the upstream interface in z-direction.
    A matrix of dimension Nx*Ny [M/L3].  Only when \code{full.output = TRUE}.
  }
  \item{C.z.down }{concentration at the downstream interface in z-direction.
    A matrix of dimension Nx*Ny [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{x.flux }{flux across the interfaces in x-direction of the grid cells.
    A (Nx+1)*Ny*Nz array [M/L2/T]. Only when \code{full.output = TRUE}.
  }
  \item{y.flux }{flux across the interfaces in y-direction of the grid cells.
    A Nx*(Ny+1)*Nz array [M/L2/T].  Only when \code{full.output = TRUE}.
  }
  \item{z.flux }{flux across the interfaces in z-direction of the grid cells.
    A Nx*Ny*(Nz+1) array [M/L2/T].  Only when \code{full.output = TRUE}.
  }
  \item{flux.x.up }{flux across the upstream boundary in x-direction,
    positive = INTO model domain. A matrix of dimension Ny*Nz [M/L2/T].
  }
  \item{flux.x.down }{flux across the downstream boundary in x-direction,
    positive = OUT of model domain. A matrix of dimension Ny*Nz [M/L2/T].
  }
  \item{flux.y.up }{flux across the upstream boundary in y-direction,
    positive = INTO model domain. A matrix of dimension Nx*Nz [M/L2/T].
  }
  \item{flux.y.down }{flux across the downstream boundary in y-direction,
    positive = OUT of model domain. A matrix of dimension Nx*Nz [M/L2/T].
  }
  \item{flux.z.up }{flux across the upstream boundary in z-direction,
    positive = INTO model domain. A matrix of dimension Nx*Ny [M/L2/T].
  }
  \item{flux.z.down }{flux across the downstream boundary in z-direction,
    positive = OUT of model domain. A matrix of dimension Nx*Ny [M/L2/T].
  }

}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
## =============================================================================
## Diffusion in 3-D; imposed boundary conditions
## =============================================================================
diffusion3D <- function(t,Y,par) {
  yy    <- array(dim=c(n,n,n),data=Y)  # vector to 3-D array
  dY   <- -r*yy        # consumption
  BND   <- matrix(nr=n,nc=n,1)   # boundary concentration

  dY <- dY + tran.3D(C=yy,
      C.x.up=BND, C.y.up=BND, C.z.up=BND,
      C.x.down=BND,C.y.down=BND,C.z.down=BND,
      D.x=Dx, D.y=Dy, D.z=Dz,
      dx = dx, dy=dy, dz=dz,full.check=TRUE)$dC
  return(list(as.vector(dY)))
}

# parameters
dy    <- dx <- dz <-1   # grid size
Dy    <- Dx <- Dz <-1   # diffusion coeff, X- and Y-direction
r     <- 0.025     # consumption rate

n  <- 10
y  <- array(dim=c(n,n,n),data=10.)

print(system.time(
  ST3 <- steady.3D(y,func=diffusion3D,parms=NULL,pos=TRUE,dimens=c(n,n,n),
                 lrw=2000000, verbose=TRUE)
))

y <- array(dim=c(n,n,n),data=ST3$y)
filled.contour(y[,,n/2],color.palette=terrain.colors)


}

\references{
  Soetaert and Herman, a practical guide to ecological modelling - using R as
  a simulation platform, 2009. Springer
}

\details{
  Do not use this with too large grid.

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
  \code{\link{tran.cylindrical}}, \code{\link{tran.spherical}} 
   for a discretisation of 3-D transport equations in cylindrical and 
   spherical coordinates
  
  \code{\link{tran.1D}},   \code{\link{tran.2D}}
}

\keyword{utilities}

