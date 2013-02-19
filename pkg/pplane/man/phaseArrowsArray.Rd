\name{phaseArrowsArray}
\alias{phaseArrowsArray}
\title{phaseArrowsArray}
\description{plotting phase space vectors}
\usage{phaseArrowsArray(z, xy, add = F, arrowHeads = 0, arrowLength = 0.5, 
    col = rev(heat.colors(150))[-(1:50)], logLength = FALSE, 
    dimnames)}
\arguments{
  \item{z}{numeric matrix (dimx,dimy,2): calculated flow at a grid, see}
  \item{xy}{numeric matrix (dimx,dimy,2): (x,y) of each gridpoint}
  \item{add}{set to TRUE to add arrows to an existing plot}
  \item{arrowHeads}{size of the arrow heads (in inches), defaults to no heads, 0.04 gives nice results}
  \item{arrowLength}{length of the vectors, set to 0 to scale vectors with magnitude of flow}
  \item{col}{color scale}
  \item{logLength}{set to TRUE to calculate colors or length for log of the flow strength}
  \item{dimnames}{character vector(2): labels of the x and y axis}
}

\value{invisible numeric matrix (2,resol,resol): calculated flow at each of the grid points before jittering}

\author{Thomas Wutzler}





