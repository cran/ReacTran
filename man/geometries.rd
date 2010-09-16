\name{g.sphere}
\alias{g.sphere}
\alias{g.spheroid}
\alias{g.cylinder}

\title{
  Surface Area and Volume of Geometrical Objects
}

\description{
  \itemize{
    \item \code{g.sphere} the surface and volume of a sphere
    \item \code{g.spheroid} the surface and volume of a spheroid
    \item \code{g.cylinder} the surface and volume of a cylinder;
      note that the surface area calculation ignores the top and bottom.
  }
}

\usage{

g.sphere(x)
g.spheroid (x, b=1)
g.cylinder (x, L=1)

}

\arguments{
  \item{x }{the radius 
  }
  \item{b }{the ratio of long/short radius of the spheroid;
    if b<1: the spheroid is oblate.
  }
  \item{L }{the length of the cylinder
  }
}

\value{
  A list containing:
  \item{surf }{the surface area 
  }
  \item{vol }{the volume 
  }
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{

 mf <- par(mfrow=c(3,2))
 x <- seq(from=0,to=1,length=10)
 plot(x, g.sphere(x)$surf,main="sphere surface")
 plot(x, g.sphere(x)$vol,main="sphere volume")
 plot(x, g.spheroid(x,b=0.5)$surf,main="spheroid surface")
 plot(x, g.spheroid(x,b=0.5)$vol,main="spheroid volume")
 plot(x, g.cylinder(x,L=1)$surf,main="cylinder surface")
 plot(x, g.cylinder(x,L=1)$vol,main="cylinder volume")
 par("mfrow"=mf)
}
\keyword{utilities}


