\name{pplane-package}
\alias{pplane-package}\alias{pplane}
\title{Phase plane plots of 2D nonlinear systems.}
\description{Plotting the the flow vector field of a 2D ODE system. With additions of trajectories and null clines.}
\details{The main function is \code{\link{phaseArrows}}, which plots vectors of the flow field of the phase plane.Further functionality of this package are \itemize{\item \code{\link{phaseNullclines}}: draw nullclines, estimated by contour line at level 0\item \code{\link{drawTrajectories}}: evaluate and draw trajectories in the phase plane\item \code{\link{phaseContours}}: draw contour lines for the y and x component of the derivative\item \code{\link{jacobianAtXY}}: Numerical approximation of the Jacobian at a point.}%itemizeSee also functionality of other packages:\itemize{\item \code{\link{eigen}}: Computes eigenvalues and eigenvectors of real (double, integer, logical) or complex matrices. }%itemizeExample systems provided with the package are: \itemize{\item \code{\link{predatorprey}}: Lotka-Volterra system \item \code{\link{competition}}: A competition system in which two populations of animals compete for the same limited resource  \item \code{\link{newtoncooling}}:   \item \code{\link{fhn}}:   \item \code{\link{disease}}:   \item \code{\link{linear}}: damped harmonic oscillator, e.g. The dynamics of a spring-mass system with air friction.   }%itemize%Further misc functionality: \itemize{%}%itemizeThis package draws from ideas and code of http://www.macalester.edu/~kaplan/math194/pplane.pdf of Prof. Daniel Kaplan. The original code is provided in origR/pplane.R in the installation directory.}%details\author{Thomas Wutzler}

\keyword{ package }


