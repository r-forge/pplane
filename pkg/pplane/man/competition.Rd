\name{competition}
\alias{competition}
\title{competition}
\description{A competition system in which two populations of animals compete for the same limited resource}
\usage{competition(mu = 2, lambda = 2, Kx = 1000, Ky = 500)}
\arguments{
  \item{mu}{coefficient, see details}
  \item{lambda}{coefficient, see details}
  \item{Kx}{coefficient, see details}
  \item{Ky}{coefficient, see details}
}
\details{derivate function for system
\cr dx = mu*(1-(x+y)/Kx)*x
\cr dy = lambda*(1-(x+y)/Ky)*y}
\value{closure calculating the derivative based on vector arguments (x,y).}

\author{Thomas Wutzler}





