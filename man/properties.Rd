\name{p.exp}
\alias{p.exp}
\alias{p.lin}
\alias{p.sig}

\title{
  Common Properties with Distance, to be used with setup.prop.1D
}

\description{
  Functions that define an y-property as a function of the one-dimensional
  x-coordinate. These routines can be used to specify properties and parameters
  as a function of distance, e.g. depth in the water column or the sediment.

  They make a transition from an upper (or upstream) zone, with value
  \code{y.0} to a lower zone with a value \code{y.inf}.

  Particularly useful in combination with \link{setup.prop.1D}
  
  \itemize{
    \item \code{p.exp}: exponentially decreasing transition 
      \deqn{
         y = y_{\inf} + (y_0-y_{\inf}) \exp(-\max(0,x-x_0)/x_a)
      }{
      y=y0+(y0-yinf)*exp(-max(0,(x-x0))/xa}

    \item \code{p.lin}: linearly decreasing transition 
      \deqn{
         y = y_0;  y = y_0 - (y_0-y_{inf})*(x-x_L)/x_{att}) ;
         y = y_{inf}
        }{y=y0 ; y=y0-(y0-yinf)(x-xl)/xatt ; y = yinf}
        
    for \eqn{0 \leq x \leq x_L}, \eqn{x_L \leq x \leq x_L + x_{att}} and
     \eqn{(x \geq x_L + x.att )} respectively.
     
    \item \code{p.sig}: sigmoidal decreasing transition
      \deqn{
         y = y_{inf} + (y_0-y_{inf})\frac{\exp(-(x-x_L)/
         (0.25 x_{att}))}{(1+\exp(-(x-x_L))/(0.25 x_{att}))})
        }{y=yinf+(y0-yinf)exp(-(x-xL)/(0.25xatt)) /
        (1+exp(-(x-xL)/(0.25xatt)))}
  }
}

\usage{
p.exp(x, y.0 = 1, y.inf = 0.5, x.L = 0, x.att = 1)
p.lin(x, y.0 = 1, y.inf = 0.5, x.L = 0, x.att = 1)
p.sig(x, y.0 = 1, y.inf = 0.5, x.L = 0, x.att = 1)
}

\arguments{
  \item{x }{the x-values for which the property has to be calculated.
  }
  \item{y.0 }{the y-value at the origin 
  }
  \item{y.inf }{the y-value at infinity 
  }
  \item{x.L }{the x-coordinate where the transition zone starts;
    for \code{x <= x.0}, the value will be equal
    to \code{y.0}. For \code{x >> x.L + x.att} the value will
    tend to \code{y.inf}
  }
  \item{x.att }{attenuation coefficient in exponential decrease,
    or the size of the transition zone in the linear and sigmoid decrease
  }
}

\value{
  the property value, estimated for each x-value.
}

\details{
  For \code{p.lin}, the width of the transition zone equals \code{x.att} and
  the depth where the transition zone starts is \code{x.L}.
  
  For \code{p.sig}, \code{x.L} is located the middle of the smooth transition zone of approaximate width \code{x.att}.  
  
  For \code{p.exp}, there is no clearly demarcated transition zone;
  there is an abrupt change at \code{x.L} after which the property
  exponentially changes from \code{y.0} towards \code{y.L} with attenuation
  coefficient \code{x.att}; the larger \code{x.att} the less steep the change.

}
\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{

x <- seq(0, 5, len = 100)
plot(x, p.exp(x, x.L = 2), 
     xlab = "x.coordinate", ylab = "y value", ylim = c(0, 1))
lines(x, p.lin(x, x.L = 2), col = "blue")
lines(x, p.sig(x, x.L = 2), col = "red")
}
 
\keyword{utilities}


