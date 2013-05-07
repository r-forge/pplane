\name{predatorprey}
\alias{predatorprey}
\title{predatorprey}
\description{Lotka-Volterra system}
\usage{predatorprey(lambda = 1, epsilon = 0.001, delta = 1, eta = 0.001)}
\arguments{
  \item{lambda}{coefficient for dx = f(x)}
  \item{epsilon}{coefficient for dx = f(xy)}
  \item{delta}{coefficient for dy = g(y)}
  \item{eta}{coefficient for dy = g(xy)}
}
\details{derivate function for system \cr
dx = f(x,y) = (lambda - epsilon*y)*x \cr
dy = g(x,y) = (-delta + eta*x)*y}
\value{closure calculating the derivative based on vector arguments (x,y).}

\author{Thomas Wutzler}




\examples{
fDeriv <- predatorprey(lambda=3, epsilon=2, delta=3, eta=2)
fDeriv(1,2)
if( FALSE ){    # do not run on R CMC check
    require(snowfall)
    # test if closure works also on remote process
    sfInit(TRUE,2)
    argsL <- list( p1=c(1,2), p2=c(2,2) )
    sfLapply( argsL, fDeriv )   # works :)
}
}
