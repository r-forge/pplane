\name{drawTrajectories}
\alias{drawTrajectories}
\title{drawTrajectories}
\description{evaluate and draw trajectories in the phase plane}
\usage{drawTrajectories(fDeriv, tEnd = 1, tCut = 60, ..., tStart = 0, 
    x0 = list(x = c(1), y = c(1)), color = "black", arrowHeads = 0.06, 
    qArrow = 0.1, fLapply = lapply, fOde = lsoda)}
\arguments{
  \item{fDeriv}{derivative function \code{function(x,y,...)}}
  \item{tEnd}{numeric scalar: end time}
  \item{tCut}{numeric scalar: number of points of the trajectory}
  \item{\dots}{further arguments to \code{fDeriv}}
  \item{tStart}{start time}
  \item{x0}{list of starting positions, such as returned by \code{\link{locator}(2,"p")}, with entries\describe{
\item{x}{numeric vector of x positions}
\item{y}{numeric vector of y positions}
}}
  \item{color}{colour of the trajectory}
  \item{arrowHeads}{size of the arrow heads (in inches), set t0 0 to avoid}
  \item{qArrow}{quantile of the timepoints at which an arrow is drawn}
  \item{fLapply}{apply function (use \code{sfLapply} for parallel calculation of trajectories)}
  \item{fOde}{function to solve the forward problem}
}

\value{invisible list of \code{\link{lsoda}} output for each trajectory}

\author{Thomas Wutzler}



\seealso{\code{\link{phaseArrows}}, \code{\link{pplane}} }

