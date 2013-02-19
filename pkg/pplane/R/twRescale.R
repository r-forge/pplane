.inside <- function (x, interval){  x >= interval[1] & x <= interval[2] }


.twRescale <- function (
        ### Rescale numeric vector to have specified minimum and maximum. 
        x               ##<< data to rescale
        ,to = c(0, 1)   ##<< range to scale to
        ,from =         ##<< range to scale from, defaults to range of data
                range(x[is.finite(x)], na.rm = TRUE)    
        ,clip = TRUE    ##<< should values be clipped to specified range?
){
    ##details<<
    ## adapted from package ggplot2 to avoid package redundancies
    
    ##author<< Hadley Wickham <h.wickham@gmail.com>, Thomas Wutzler
    
    ##details<< 
    ## If from[1] == from[2] then the mean of interval to is returned.
    if( length(from) != 2 || length(to) != 2)
        stop("twRescale: arguments to and from must be ranges.")
    if ( from[1] == from[2] ) 
        return( rep( mean(to), length(x) ) )
    if (is.factor(x)) {
        warning("twRescale: Categorical variable automatically converted to continuous", 
                call. = FALSE)
        x <- as.numeric(x)
    }
    scaled <- (x - from[1])/diff(from) * diff(to) + to[1]
    if (clip) {
        ifelse(!is.finite(scaled) | .inside(scaled,to), scaled, 
                NA)
    }
    else {
        scaled
    }
}

