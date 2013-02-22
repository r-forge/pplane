\name{phaseArrows}
\alias{phaseArrows}
\title{phaseArrows}
\description{plotting phase space vectors}
\usage{phaseArrows(fDeriv, xlims = c(-3, 3), ylims = xlims, resol = 20, 
    isJitter = FALSE, add = F, arrowHeads = 0.04, arrowLength = 0.5, 
    col = rev(heat.colors(150))[-(1:50)], logLength = FALSE, 
    dimnames = NULL, ...)}
\arguments{
  \item{fDeriv}{derivative function \code{function(x,y,...)}}
  \item{xlims}{numeric vector (2): range of the x values}
  \item{ylims}{numeric vector (2): range of the y values}
  \item{resol}{scalar integer: number of points in x and y range}
  \item{isJitter}{set to TRUE to jitter vector starting positions to avoid overplotting 
(however calculating contours does not work then anymore) }
  \item{add}{set to TRUE to add arrows to an existing plot}
  \item{arrowHeads}{size of the arrow heads (in inches), set t0 0 for heads, 0.04 gives nice results}
  \item{arrowLength}{length of the vectors, set to 0 to scale vectors with magnitude of flow}
  \item{col}{color scale}
  \item{logLength}{set to TRUE to calculate colors or length for log of the flow strength}
  \item{dimnames}{character vector(2): labels of the x and y axis, to overwrite defaults}
  \item{\dots}{further arguments to \code{fun}, such as \code{parms}}
}

\value{result of \code{\link{evalDerivGrid}}: a list with entries \itemize{
\item z: numeric matrix vector (resol x resol,2): calculated flow at each of the grid (labels before before jittering)
\item xy: numeric matrix (resol,resol,2): x and y values of the grid
}    }

\author{Thomas Wutzler}




\examples{
fDeriv <- predatorprey(lambda=3, epsilon=2, delta=3, eta=2)
#
# set up a plotting window
# windows(width = 4.6, height = 3.2, pointsize = 10); par( las = 1, mar = c(2, 3.3, 0, 0) + 0.3, tck = 0.02, mgp = c(1.1, 0.2, 0))
#
# default: strength is colored, no arrow heads, same length
tmp <- phaseArrows( fDeriv, c(-2,5),c(-2,5) );
#phaseContours(tmp$z, tmp$xy, add=TRUE)
phaseNullclines(tmp$z, tmp$xy, add=TRUE)
drawTrajectories( fDeriv, tEnd=3, x0=list(x=1,y=2) )    # initial starting point by script
#drawTrajectories( fDeriv, tEnd=3, x0=locator(2,"p") )   # set the starting point in the graph: need to click several times
#
# add arrow heads and use colors for log scale, scale vector length
phaseArrows( fDeriv, c(-2,5),c(-2,5), logLength=TRUE, arrowHeads=0.04, arrowLength=0 );
#
# for background, sometimes a decent color is useful
tmp <- phaseArrows( fDeriv, c(-2,5),c(-2,5), col="grey" );
#
# may use parallel calculation of flow field and trajectories
if( FALSE ){    # do not run on R CMC check
    require(snowfall)
    tmp <- phaseArrows( fDeriv, c(-2,5),c(-2,5), useSnowfall=TRUE );       # using parallel calculation
    drawTrajectories( fDeriv, tEnd=3, x0=list(x=1,y=c(2:5)), fLapply=sfLapply )    # initial starting point by script
}
}
