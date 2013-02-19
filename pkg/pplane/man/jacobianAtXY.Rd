\name{jacobianAtXY}
\alias{jacobianAtXY}
\title{jacobianAtXY}
\description{Numerical approximation of the Jacobian at a point.}
\usage{jacobianAtXY(fDeriv, x = NULL, y = NULL, h = 1e-06)}
\arguments{
  \item{fDeriv}{derivative function \code{function(x,y,...)}}
  \item{x}{numeric scalar x coordinate}
  \item{y}{numeric scalar y corrdinate}
  \item{h}{step to approximate derivative}
}

\value{numeric matrix (2x2) of numerical approximation of Jacobian }

\author{Thomas Wutzler}





