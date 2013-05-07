\name{fhn}
\alias{fhn}
\title{fhn}
\description{fhn}
\usage{fhn(gamma = 2.5, epsilon = 1, a = 0.3)}
\arguments{
  \item{gamma}{coefficient, see details}
  \item{epsilon}{coefficient, see details}
  \item{a}{coefficient, see details}
}
\details{derivate function for system \cr
dx = -x*(x-a)*(x-1) - y \cr
dy = epsilon*(x - gamma*y)}
\value{closure calculating the derivative based on vector arguments (x,y).}

\author{Thomas Wutzler}





