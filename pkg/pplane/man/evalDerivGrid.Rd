\name{evalDerivGrid}
\alias{evalDerivGrid}
\title{evalDerivGrid}
\description{evaluating derivative function at grid}
\usage{evalDerivGrid(fDeriv, xlims, ylims, resol = 11, isJitter = TRUE, 
    ..., useSnowfall = FALSE)}
\arguments{
  \item{fDeriv}{derivative function \code{function(x,y,...)}}
  \item{xlims}{numeric vector (2): range of the x values}
  \item{ylims}{numeric vector (2): range of the y values}
  \item{resol}{scalar integer: number of points in x and y range}
  \item{isJitter}{set to FALSE to avoid jittering vector starting positions}
  \item{\dots}{further arguments to \code{fun}, such as \code{parms}}
  \item{useSnowfall}{set to TRUE to use parallel execution using snowfall}
}

\value{list with entries
\item{z}{numeric matrix vector (resol x resol,2): calculated flow at each of the grid (labels before before jittering)}
\item{xy}{numeric matrix (resol,resol,2): x and y values of the grid
}}

\author{Thomas Wutzler}





