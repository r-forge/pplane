\name{linear}
\alias{linear}
\title{linear}
\description{damped harmonic oscillator, e.g. The dynamics of a spring-mass system with air friction. }
\usage{linear(x0, y0, a, b, c, d, e = 0, f = 0)}
\arguments{
  \item{x0}{initial positions (x0,y0)}
  \item{y0}{initial positions}
  \item{a}{entry (1,1) of the transformation matrix}
  \item{b}{entry (1,2) of the transformation matrix}
  \item{c}{entry (2,1) of the transformation matrix}
  \item{d}{entry (2,2) of the transformation matrix}
  \item{e}{x-offset for affine systems}
  \item{f}{y-offset for affine systems}
}
\details{In physics, this sort of system is called a damped
harmonic oscillator. The x variable is the spring position,
the y variable is the spring velocity.

derivate function for system
\cr dx = a*(x-x0) + b*(y-y0) + e
\cr dy = c*(x-x0) + d*(y-y0) + f}
\value{closure calculating the derivative based on vector arguments (x,y).}

\author{Thomas Wutzler}





