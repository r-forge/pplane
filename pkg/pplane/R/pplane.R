jacobianAtXY <- function(
    ### Numerical approximation of the Jacobian at a point.
    fDeriv      ##<< derivative function \code{function(x,y,...)}
    ,x=NULL     ##<< numeric scalar x coordinate
    ,y=NULL     ##<< numeric scalar y corrdinate
    ,h=.000001  ##<< step to approximate derivative  
){
    if (is.null(x) | is.null(y)) {
        x0 <- locator(n=1);
        x <- x0$x; 
        y <- x0$y;  
    }
    foo <- fDeriv(x,y);
    foox <- fDeriv(x+h,y);
    fooy <- fDeriv(x,y+h);
    A <- (foox[1] - foo[1])/h;
    B <- (fooy[1] - foo[1])/h;
    C <- (foox[2] - foo[2])/h;
    D <- (fooy[2] - foo[2])/h;
    ##value<< numeric matrix (2x2) of numerical approximation of Jacobian 
    return(matrix( c(A,B,C,D ), 2,2, byrow=T))
}


evalDerivGrid <- function(
    ### evaluating derivative function at grid
    fDeriv      ##<< derivative function \code{function(x,y,...)}
    ,xlims      ##<< numeric vector (2): range of the x values
    ,ylims      ##<< numeric vector (2): range of the y values
    ,resol=11   ##<< scalar integer: number of points in x and y range
    ,isJitter=TRUE      ##<< set to FALSE to avoid jittering vector starting positions 
    ,...                ##<< further arguments to \code{fun}, such as \code{parms}
    ,useSnowfall=FALSE  ##<< set to TRUE to use parallel execution using snowfall
){
    xGrid <- seq(xlims[1],xlims[2], length=resol)
    yGrid <- seq(ylims[1],ylims[2], length=resol)
    x <- matrix(xGrid, byrow=T, resol,resol);
    y <- matrix(yGrid, byrow=F, resol, resol);
    if( isTRUE(isJitter) ){
        # jitter 
        xspace <- abs(diff(xlims))/(resol*5)
        yspace <- abs(diff(ylims))/(resol*5)
        npts <- resol*resol;
        x <- x + matrix(runif(npts, -xspace, xspace),resol,resol)
        y <- y + matrix(runif(npts, -yspace, yspace),resol,resol)
    }
    xy <- abind(x,y,rev.along=0)
    #z <- apply( xy, c(1,2), function(X){ fun(0, X, ...)[[1]] })
    zv1 <- fDeriv( x[1], y[1] )     # in order to get proper variable names
    zv <- if( isTRUE(useSnowfall)){
        ncpu = sfCpus()
        nChunks <- length(c(x)) %/% ncpu
        iChunk <- 1+(seq_along(c(x))-1) %/% nChunks
        fRemote <- function(iCpu,iChunk, fDeriv, xv, yv, ...){
            ii <- which(iChunk==iCpu)
            matrix( fDeriv(xv[ii], yv[ii]), ncol=2 )
        }
        res <- sfLapply( x=1:ncpu, fun=fRemote, iChunk=iChunk, fDeriv=fDeriv, xv=c(x), yv=c(y), ...)
        abind( res, along=1)
    }else fDeriv( c(x), c(y) )
    z <- array( zv, dim=c(resol,resol,2) )
    dimnames(z) <- list( x=NULL, y=NULL, names(zv1))
    ##value<< list with entries
    list(
       z=z      ##<< numeric matrix vector (resol x resol,2): calculated flow at each of the grid (labels before before jittering)
       ,xy=xy   ##<< numeric matrix (resol,resol,2): x and y values of the grid
       )
}

phaseArrowsArray <- function(
    ### plotting phase space vectors
    z                   ##<< numeric matrix (dimx,dimy,2): calculated flow at a grid, see 
    ,xy                 ##<< numeric matrix (dimx,dimy,2): (x,y) of each gridpoint
    ,add=F              ##<< set to TRUE to add arrows to an existing plot
    ,arrowHeads=0       ##<< size of the arrow heads (in inches), defaults to no heads, 0.04 gives nice results 
    ,arrowLength=0.5    ##<< length of the vectors, set to 0 to scale vectors with magnitude of flow
    ,col=rev(heat.colors(150))[-(1:50)]  ##<< color scale  
    ,logLength=FALSE    ##<< set to TRUE to calculate colors or length for log of the flow strength
    ,dimnames           ##<< character vector(2): labels of the x and y axis
){
    if( length(dimnames)==0 )
        dimnames = if( length(rownames(z)) ) rownames(z) else c("dx","dy") 
    xlims <- range( xy[,,1])
    ylims <- range( xy[,,2])
    resol <- nrow(xy)
    z1 <- z[,,1]
    z2 <- z[,,2]
    lens <- sqrt(z1^2 + z2^2)
    lens2 <- if( arrowLength==0){
        maxx <- max(abs(z1))
        maxy <- max(abs(z2))
        dt2 <- min( abs(diff(xlims))/maxx, abs(diff(ylims))/maxy)/resol
        lens2 <- (lens/max(lens) +0.1)/(dt2)
    }else{
        dt <- min( abs(diff(xlims)), abs(diff(ylims)))/resol
        lens2 <- lens / dt / arrowLength
    }
    nCol <- length(col)
    colVal <- if(nCol==1 ) col else {
        lens3 <- if(isTRUE(logLength)) log(lens) else lens
        col[ round( .twRescale( abs(lens3), to=c(1, nCol) ) ) ]
    }
    if (add==F) {
        plot(1,xlim=xlims, ylim=ylims, type='n', xlab=dimnames[1], ylab=dimnames[2] )
    }
    arrows(c(xy[,,1]), c(xy[,,2]), c(xy[,,1]+z1/lens2), c(xy[,,2]+z2/lens2), length=arrowHeads, col=colVal)
    ##value<< invisible numeric matrix (2,resol,resol): calculated flow at each of the grid points before jittering
    invisible(NULL)
}

phaseArrows <- function(
        ### plotting phase space vectors
        fDeriv      ##<< derivative function \code{function(x,y,...)}
        ,xlims      ##<< numeric vector (2): range of the x values
        ,ylims      ##<< numeric vector (2): range of the y values
        ,resol=20   ##<< scalar integer: number of points in x and y range
        ,isJitter=TRUE      ##<< set to FALSE to avoid jittering vector starting positions 
        ,add=F              ##<< set to TRUE to add arrows to an existing plot
        ,arrowHeads=0       ##<< size of the arrow heads (in inches), defaults to no heads, 0.04 gives nice results 
        ,arrowLength=0.5    ##<< length of the vectors, set to 0 to scale vectors with magnitude of flow
        ,col=rev(heat.colors(150))[-(1:50)]  ##<< color scale  
        ,logLength=FALSE    ##<< set to TRUE to calculate colors or length for log of the flow strength
        ,dimnames=NULL      ##<< character vector(2): labels of the x and y axis, to overwrite defaults
        ,...                ##<< further arguments to \code{fun}, such as \code{parms}
){
    res <- evalDerivGrid(fDeriv=fDeriv, xlims=xlims, ylims=ylims, resol=resol, isJitter=isJitter,...)
    phaseArrowsArray(z=res$z, xy=res$xy, add=add, arrowHeads=arrowHeads, arrowLength=arrowLength, col=col, logLength=logLength, dimnames=dimnames)
    ##value<< result of \code{\link{evalDerivGrid}}: a list with entries \itemize{
    ## \item z: numeric matrix vector (resol x resol,2): calculated flow at each of the grid (labels before before jittering)
    ## \item xy: numeric matrix (resol,resol,2): x and y values of the grid
    ## }    
    invisible(res)
}
attr( phaseArrows,"ex") <- function(){
    fDeriv <- predatorprey(lambda=3, epsilon=2, delta=3, eta=2)
    #
    # set up a plotting window
    windows(width = 4.6, height = 3.2, pointsize = 10)
    par( las = 1, mar = c(2, 3.3, 0, 0) + 0.3, tck = 0.02, mgp = c(1.1, 0.2, 0))
    #
    # default: strength is colored, no arrow heads, same length
    tmp <- phaseArrows( fDeriv, c(-2,5),c(-2,5) );
    #phaseContours(tmp$z, tmp$xy, add=TRUE)
    phaseNullclines(tmp$z, tmp$xy, add=TRUE)
    drawTrajectories( fDeriv, tEnd=3, x0=list(x=1,y=2) )    # initial starting point by script
    #drawTrajectories( fDeriv, tEnd=3, x0=locator(2,"p") )   # set the starting point in the graph: need to click several times
    #
    # add arrow heads and use colors for log scale, scale vector length
    phaseArrows( fDeriv, c(-2,5),c(-2,5), logLength=TRUE, arrowHeads=0.04, arrowLength=0 );
    #
    # for background, sometimes a decent color is useful
    tmp <- phaseArrows( fDeriv, c(-2,5),c(-2,5), col="grey" );
    #
    # may use parallel calculation of flow field and trajectories
    if( FALSE ){    # do not run on R CMC check
        require(snowfall)
        tmp <- phaseArrows( fDeriv, c(-2,5),c(-2,5), useSnowfall=TRUE );       # using parallel calculation
        drawTrajectories( fDeriv, tEnd=3, x0=list(x=1,y=c(2:5)), fLapply=sfLapply )    # initial starting point by script
    }
}

phaseContours <- function(
    ### draw contour lines for the y and x component of the derivative
    z       ##<< result component of \code{\link{evalDerivGrid}}
    ,xy     ##<< result component of \code{\link{evalDerivGrid}}
    ,add=F  ##<< set to TRUE to add arrows to an existing plot 
    ,colors=c('red', 'blue')    ##<< colors of the contourlines for x and y direction of the derivative
    , ...       ##<< further arguments to \code{\link{contour}}, such as \code{levels}
) {
    z1 <- z[,,1]
    z2 <- z[,,2]
    x <- xy[,,1]
    y <- xy[,,2]
    contour(x[1,],y[,1],z1, add=add, col=colors[1], ...);
    contour(x[1,],y[,1],z2, add=T, col=colors[2],...); 
}

phaseNullclines <- function(
    ### draw nullclines, estimated by contour line at level 0
    ...     ##<< arguments to \code{\link{phaseContours}}
){
    phaseContours(..., levels=c(0))
}


drawTrajectories <- function(
    ### evaluate and draw trajectories in the phase plane
    fDeriv      ##<< derivative function \code{function(x,y,...)}
    , tEnd=1    ##<< numeric scalar: end time
    , tCut=60   ##<< numeric scalar: number of points of the trajectory
    , ...       ##<< further arguments to \code{fDeriv}
    , tStart=0  ##<< start time
    , x0=list( ##<< list of starting positions, such as returned by \code{\link{locator}(2,"p")}, with entries 
            ##describe<<
            x=c(1)      ##<< numeric vector of x positions
            ,y=c(1)     ##<< numeric vector of y positions
        ##end<<
    )  
    , color = "black"	##<< colour of the trajectory
    , arrowHeads=0.06   ##<< size of the arrow heads (in inches), set t0 0 to avoid
    , qArrow = 0.1      ##<< quantile of the timepoints at which an arrow is drawn
    , fLapply=lapply    ##<< apply function (use \code{sfLapply} for parallel calculation of trajectories)
    , fOde=lsoda        ##<< function to solve the forward problem
) {
    ##seealso<< \code{\link{phaseArrows}}, \code{\link{pplane}} 
    #print(paste("Click", loc.num, "initial values"))
    #x0 <- locator(loc.num, "p")
    func <- function(time, y, parms){
        deriv <- parms$fDeriv( y[1], y[2], ... )
        list( deriv )
    }
    #i <- 1
    # make vectors the correct length
    nTr <- max(sapply( x0, length ))
    x <- numeric(nTr);    x[] <- x0[[1]]
    y <- numeric(nTr);    y[] <- x0[[2]]
    c(fOde, fDeriv)        # evaluate in current context
    traj <- fLapply( 1:nTr, function(i){        
        y <-c(x=x[i], y=y[i]) 
        out <- as.data.frame(fOde(func=func, y=y, parms=list(fDeriv=fDeriv), times = seq(tStart, tEnd, length = tCut)))
    })
    sapply( traj, function(out){
                lines(out$x, out$y, col = color)
                iArrow <- round(qArrow*nrow(out))
                outn <- out[iArrow+c(0,1),][,-1]
                arrows( outn[1,1], outn[1,2], outn[2,1], outn[2,2], length=arrowHeads, col=color )
            })
    ##value<< invisible list of \code{\link{lsoda}} output for each trajectory
    return(invisible(traj))
}




