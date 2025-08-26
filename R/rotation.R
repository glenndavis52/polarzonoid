

rotationfromarcs <- function( arcmat, tol=5.e-9 )
    {
    #pac = 'onion'
    #if( ! requireNamespace( pac, quietly=TRUE ) )
    #    {
    #    log_level( ERROR, "required package '%s' could not be loaded.", pac )
    #    return(NULL)
    #    }

    arcmat   = prepareNxM( arcmat, 2 )
    if( is.null(arcmat) )    return(NULL)

    #   check for 1 or 2 arcs
    ok  = nrow(arcmat) %in% 1:2

    if( ! ok )
        {
        log_level( ERROR, "nrow(arcmat) = %d, which is invalid.", nrow(arcmat) )
        return(NULL)
        }

    #   check the total arc length
    L   = sum( arcmat[ ,2] )

    if( tol < abs(L - pi) )
        {
        log_level( ERROR, "Total length L=%g is invalid.  tol=%g < %g = abs(L - pi).", L, tol, abs(L - pi) )
        return(NULL)
        }

    #   set n=2 to force u to have length 5
    u   = spherefromarcs( arcmat, n=2 )

    if( is.null(u) )    return(NULL)

    # cat( "u = ", spherefromboundary( p ), '\n' )
    
    # u   = rev( u[1:4] )   #just for testing

    out = rotationfromquat( u[ 1:4 ]  ) # ignore the last one which is near 0 because L is near pi

    # out = onion::as.orthogonal( onion::as.quaternion( u, single=TRUE ) )

    # out = rotmatrix( quaternion(u) )@x [ , , 1]

    return( out )
    }


arcsfromrotation <- function( rotation, tol=1.e-7 )   # tol=5.e-9 is too small for rotations around the x axis [1,0,0]
    {
    #pac = 'onion'
    #if( ! requireNamespace( pac, quietly=TRUE ) )
    #    {
    #    log_level( ERROR, "required package '%s' could not be loaded.", pac )
    #    return(NULL)
    #    }
    
    ok  = is.numeric(rotation)  &&  all( dim(rotation)==c(3,3) )
    
    if( ! ok )
        {
        log_level( ERROR, "argument rotation is invalid; it is not a 3x3 numeric matrix." )
        return(NULL)
        }

    q   = quatfromrotation( rotation )
    
    # q   = rev( q )   #just for testing
    
    p   = boundaryfromsphere( q )
    if( is.null(p) )    return(NULL)

    arcmat  = arcsfromboundary( p, tol=tol )

    return( arcmat )
    }


rotationfromquat    <- function( quat )
    {
    #   copied from onion/R/orthogonal.R  by Robin Hankin
    
    #   quat    = unitize( quat )
    
    r   = quat[1]
    i   = quat[2]
    j   = quat[3]
    k   = quat[4]
  
    out =   c(
            1-2*(j^2+k^2) ,   2*(i*j-k*r) ,   2*(i*k+j*r),
            2*(i*j+k*r)   , 1-2*(i^2+k^2) ,   2*(j*k-i*r),
            2*(i*k-j*r)   ,   2*(j*k+i*r) ,  1-2*(i^2+j^2)
            )
            
    matrix( out, 3, 3, byrow=TRUE )
    }
    
rotationaroundaxis <- function( axis, theta )
    {
    #   validate axis
    if( length(axis) != 3 )
        {
        log_level( ERROR, "argument axis is invalid.  length(axis) = %d, but it must be 3.", length(axis) )
        return(NULL)
        }
        
    axis    = unitize(axis)
    if( is.na(axis[1]) )
        {
        log_level( ERROR, "argument axis is invalid.  It is 0." )
        return(NULL)
        }

    rotationfromquat( c( cos(theta/2), sin(theta/2)*axis ) )
    }
    

quatfromrotation    <- function( M )
    {
    #   copied from onion/R/orthogonal.R  by Robin Hankin
    #   omitted test for orthogonality
    
    trace   = M[1,1]+M[2,2]+M[3,3]
    
    if ( 0 < trace ) { 
        S <- 2*sqrt(1 + M[1,1] + M[2,2] + M[3,3])
        qw <- S/4
        qx <- (M[3,2] - M[2,3]) / S
        qy <- (M[1,3] - M[3,1]) / S 
        qz <- (M[2,1] - M[1,2]) / S 
    } else if ((M[1,1] > M[2,2]) && (M[1,1] > M[3,3])) { 
        S <- 2*sqrt(1 + M[1,1] - M[2,2] - M[3,3])
        qw <- (M[3,2] - M[2,3]) / S
        qx <- S/4
        qy <- (M[1,2] + M[2,1]) / S 
        qz <- (M[1,3] + M[3,1]) / S
    } else if (M[2,2] > M[3,3]) { 
        S <- 2*sqrt(1 + M[2,2] - M[1,1] - M[3,3])
        qw <- (M[1,3] - M[3,1]) / S
        qx <- (M[1,2] + M[2,1]) / S
        qy <- S/4
        qz <- (M[2,3] + M[3,2]) / S
    } else { 
        S <- 2*sqrt(1 + M[3,3] - M[1,1] - M[2,2])
        qw <- (M[2,1] - M[1,2]) / S
        qx <- (M[1,3] + M[3,1]) / S
        qy <- (M[2,3] + M[3,2]) / S
        qz <- S/4
    }
    
    return( c(qw,qx,qy,qz) )
    }



#########################       deadwood below      #######################################

rotationfromquat_old    <- function( quat )
    {
    #onion::as.orthogonal( onion::as.quaternion( quat, single=TRUE ) )
    }

quatfromrotation_old    <- function( rotation )
    {
    #as.numeric( onion::matrix2quaternion( rotation ) )
    }
