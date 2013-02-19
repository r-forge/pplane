
# start at c(200,30)
predatorprey <- function(
    ### Lotka-Volterra system
    lambda=1            ##<< coefficient for dx = f(x) 
    , epsilon=.001      ##<< coefficient for dx = f(xy)
    ,delta=1            ##<< coefficient for dy = g(y)
    , eta=.001          ##<< coefficient for dy = g(xy)
){
    ##details<< derivate function for system 
    ## \cr dx = f(x,y) = (lambda - epsilon*y)*x
    ## \cr dy = g(x,y) = (-delta + eta*x)*y
    function(x,y=NULL){
        if (is.null(y)) {
            y <- x[2]; x <- x[1];
        } 
        dx <- (lambda - epsilon*y)*x;
        dy <- (-delta + eta*x)*y;
        return( c(dx=dx, dy=dy) )
    }
    ##value<< closure calculating the derivative based on vector arguments (x,y).
}
attr(predatorprey,"ex") <- function(){
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

competition <- function(
    ### A competition system in which two populations of animals compete for the same limited resource
    mu=2        ##<< coefficient, see details 
    , lambda=2  ##<< coefficient, see details
    , Kx=1000   ##<< coefficient, see details
    , Ky=500    ##<< coefficient, see details
) {
    ##details<< 
    ## derivate function for system
    ## \cr dx = mu*(1-(x+y)/Kx)*x
    ## \cr dy = lambda*(1-(x+y)/Ky)*y
    function(x,y=NULL){
        if (is.null(y)) {
            y<- x[2]; x <- x[1];
        }
        dx <- mu*(1-(x+y)/Kx)*x;
        dy <- lambda*(1-(x+y)/Ky)*y;
        return( c(dx, dy) );
    }
    ##value<< closure calculating the derivative based on vector arguments (x,y).
}

newtoncooling <- function(
    ### NewtonColling
    a1  ##<< coefficient, see details
    ,a2 ##<< coefficient, see details
) {
    ##details<< 
    ## derivate function for system
    ## \cr dx = a1*(y-x)
    ## \cr dy = a2*(x-y)
    function(x,y=NULL){
        if (is.null(y)) {
            y<- x[2]; x <- x[1];
        }
        dx <- a1*(y-x);
        dy <- a2*(x-y);
        return( c(dx, dy) );
    }
    ##value<< closure calculating the derivative based on vector arguments (x,y).
}

fhn <- function(
    ### fhn
    gamma=2.5   ##<< coefficient, see details
    , epsilon=1 ##<< coefficient, see details
    , a=0.3     ##<< coefficient, see details
) {
    ##details<< 
    ## derivate function for system
    ## \cr dx = -x*(x-a)*(x-1) - y
    ## \cr dy = epsilon*(x - gamma*y)
    function(x,y=NULL){
        if (is.null(y)) {
            y<- x[2]; x <- x[1];
        }
#  gamma = 2.5;
#  epsilon=1;
#  a = 0.3;
        dx <- -x*(x-a)*(x-1) - y;
        dy <- epsilon*(x - gamma*y);
        return( c(dx, dy) );
    }
    ##value<< closure calculating the derivative based on vector arguments (x,y).
}

disease <- function(
    ### disease
    b=1     ##<< coefficient, see details
    ,mu=1   ##<< coefficient, see details
    ,C=2    ##<< coefficient, see details
){
    ##details<< 
    ## derivate function for system
    ## \cr dx = -x*(x-a)*(x-1) - y
    ## \cr dy = epsilon*(x - gamma*y)
    function(x,y=NULL){
        if (is.null(y)) {
            y<- x[2]; x <- x[1];
        }
        dx <- b-C*x*y;
        dy <- C*x*y - mu*y;
        return( c(dx, dy) );
    }
    ##value<< closure calculating the derivative based on vector arguments (x,y).
}

linear <- function(
        ### damped harmonic oscillator, e.g. The dynamics of a spring-mass system with air friction. 
        x0      ##<< initial positions (x0,y0)
        ,y0     ##<< initial positions
        ,a      ##<< entry (1,1) of the transformation matrix
        ,b      ##<< entry (1,2) of the transformation matrix
        ,c      ##<< entry (2,1) of the transformation matrix
        ,d      ##<< entry (2,2) of the transformation matrix
        ,e=0    ##<< x-offset for affine systems
        ,f=0    ##<< y-offset for affine systems
) {
    ##details<<
    ## In physics, this sort of system is called a damped
    ## harmonic oscillator. The x variable is the spring position,
    ## the y variable is the spring velocity.
    ##
    ## derivate function for system
    ## \cr dx = a*(x-x0) + b*(y-y0) + e
    ## \cr dy = c*(x-x0) + d*(y-y0) + f
    foo <- function(x,y=NULL) {
        if (is.null(y)) {
            y<- x[2]; x <- x[1];
        }
        dx <- a*(x-x0) + b*(y-y0) + e;
        dy <- c*(x-x0) + d*(y-y0) + f;
        return( c(dx, dy) );
    }
    ##value<< closure calculating the derivative based on vector arguments (x,y).
}


