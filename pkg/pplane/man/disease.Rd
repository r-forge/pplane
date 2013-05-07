\name{disease}
\alias{disease}
\title{disease}
\description{disease}
\usage{disease(b = 1, mu = 1, C = 2)}
\arguments{
  \item{b}{coefficient, see details}
  \item{mu}{coefficient, see details}
  \item{C}{coefficient, see details}
}
\details{derivate function for system \cr
dx = -x*(x-a)*(x-1) - y \cr
dy = epsilon*(x - gamma*y)}
\value{closure calculating the derivative based on vector arguments (x,y).}

\author{Thomas Wutzler}





