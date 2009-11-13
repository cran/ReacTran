
##==============================================================================
## Attaches a property to a 2D grid
##==============================================================================

setup.prop.2D <- function(func = NULL, y.func = func, value = NULL, y.value = value, grid, ...) {

  ## check input
  gn <- names(grid)
  if (! "x.mid"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.mid")
  if (! "x.int"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.int")
  if (! "y.mid"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.mid")
  if (! "y.int"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.int")

  if (is.null(func) && is.null(value))
    stop("error in setup.prop: function and value should not be both NULL")
  if (is.null(y.func) && is.null(y.value))
    stop("error in setup.prop: y.function and y.value should not be both NULL")

  if (!is.null(value)) { # profile specification via value
    x.int <- matrix(nrow=length(grid$x.int),ncol=length(grid$y.mid),data=value)
    y.int <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.int),data=y.value)
    x.mid <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.mid),data=value)
    y.mid <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.mid),data=y.value)
  }

  if (!is.null(func)) { # profile specification via function
    x.int <- matrix(nrow=length(grid$x.int),ncol=length(grid$y.mid))
    x.mid <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.mid))

    for (i in 1:length(grid$x.int))
    {
      for (j in 1:length(grid$y.mid))
          x.int[i,j] <- func(grid$x.int[i],grid$y.mid[j],...)
    }
    for (i in 1:length(grid$x.mid))
    {
      for (j in 1:length(grid$y.mid))
          x.mid[i,j] <- func(grid$x.mid[i],grid$y.mid[j],...)
    }
  }

  if (!is.null(y.func)) { # profile specification via function
    y.int <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.int))
    y.mid <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.mid))

    for (i in 1:length(grid$x.mid))
    {
      for (j in 1:length(grid$y.int))
          y.int[i,j] <- y.func(grid$x.mid[i],grid$y.int[j],...)
    }
    for (i in 1:length(grid$x.mid))
    {
      for (j in 1:length(grid$y.mid))
          y.mid[i,j] <- y.func(grid$x.mid[i],grid$y.mid[j],...)
    }
  }
  
  Res <- list(x.mid = x.mid,
              y.mid = y.mid, 
              x.int = x.int,
              y.int = y.int)
  class(Res) <- "prop.2D"
  return(Res)
}

##==============================================================================
## S3 method: Plotting of a two-dimensional grid property
##==============================================================================

contour.prop.2D <- function(x, grid, xyswap = FALSE, filled = FALSE, ...) {
  if (! filled) {
  if (xyswap)
    contour(x=grid$y.mid,y=rev(-grid$x.mid),z=t(x$x.mid),...)
  else
    contour(x=grid$x.mid,y=grid$y.mid,z=x$x.mid,...)
  } else {
    if (xyswap)
      filled.contour(x=grid$y.mid,y=rev(-grid$x.mid),z=t(x$x.mid),...)
    else
      filled.contour(x=grid$x.mid,y=grid$y.mid,z=x$x.mid,...)
  }
}
