

###########     argument processing     ##############
#
#   A   a non-empty numeric NxM matrix, or something that can be converted to be one
#
#   Nmin    the minimum allowed number of rows
#
#   returns such a matrix, or NULL in case of error

prepareNxM  <-  function( A, M, Nmin=1 )
    {
    ok  = is.numeric(A) &&  M*Nmin<=length(A)  &&  (length(dim(A))<=2)  # &&  (0<M)

    ok  = ok  &&  ifelse( is.matrix(A), ncol(A)==M, ((length(A) %% M)==0)  )

    if( ! ok )
        {
        #print( "prepareNx3" )
        #print( sys.frames() )
        mess    = substr( paste0(as.character(A),collapse=','), 1, 10 )
        #arglist = list( ERROR, "A must be a non-empty numeric Nx3 matrix (with N>=%d). A='%s...'", mess )
        #do.call( log.string, arglist, envir=parent.frame(n=3) )
        #myfun   = log.string
        #environment(myfun) = parent.frame(3)

        Aname = deparse(substitute(A))

        #   notice hack with 2L to make log.string() print name of parent function
        #log.string( c(ERROR,2L), "Argument '%s' must be a non-empty numeric Nx%d matrix (with N>=%d). %s='%s...'",
        #                            Aname, M, Nmin, Aname, mess )

        log_level( ERROR, "Argument '%s' must be a non-empty numeric Nx%d matrix (with N>=%d). %s='%s...'",
                                    Aname, M, Nmin, Aname, mess, .topcall=sys.call(-2L) )
        return(NULL)
        }

    if( ! is.matrix(A) )
        A = matrix( A, ncol=M, byrow=TRUE )

    return( A )
    }



#   u, v    unit vectors of dimension n
#
#   returns nxn rotation matrix that takes u to v, and fixes all vectors ortho to u and v
#
#   Characteristic Classes, Milnor and Stasheff, p. 77

rotationshortest <- function( u, v )
    {
    n   = length(u)
    if( length(v) != n )    return(NULL)

    uv  = u + v
    if( all(uv == 0) )  return(NULL)

    out = diag(n)  -  (uv %o% uv)/(1 + sum(u*v))  +  2*(v %o% u)

    return(out)
    }



#   u   a unit vector
#
#   returns a rotation matrix that rotates u to the closest pole:
#       either (0,0,...,+1) or (0,0,...,-1)   in case of a tie, then +1
#   the rotation is the shortest one

rotation2pole <- function( u )
    {
    m   = length(u)
    s   = sign(u[m])
    if( s == 0 )    s = 1
    pole    = c( numeric(m-1), s )

    return( rotationshortest(u,pole) )
    }

unitize <- function( x, tol=5.e-14 )
    {
    r2  = sum( x^2 )
    if( r2 == 0 )    return( rep(NA_real_,length(x)) )

    r   = sqrt(r2)

    if( r <= tol )
        {
        log_level( WARN, "length of x is %g <= %g", r, tol )
        }

    return( x / r )
    }


#   x       a real vector, typically of even length
#           if odd length, then the last one is ignored
#
#   returns a complex vector, the real and imag parts are taken alternatingly from vec

complexfromrealvector <- function( x )
    {
    n   = as.integer( length(x)/2 )

    if( n == 0 )    return( complex(0) )

    even    = 2L * ( 1:n )
    odd     = even - 1L

    return( complex( real=x[odd], imaginary=x[even] ) )
    }


#   z   complex real vector
#
#   returns real,imaginary parts alternating
realfromcomplexvector <- function( z )
    {
    if( ! is.complex(z) )
        {
        log_level( ERROR, "z is invalid; it is not complex." )
        return( NA_real_ )
        }
        
    mat = rbind( Re(z), Im(z) )

    return( as.numeric(mat) )
    }

#   z   complex real vector, of length N
#
#   return 2 x 2N matrix

real2x2Nfromcomplex <- function( z )
    {
    if( ! is.complex(z) )
        {
        log_level( ERROR, "z is invalid; it is not complex." )
        return( NA_real_ )
        }
        
    x   = Re(z)
    y   = Im(z)
    
    mat1    = rbind( x, -y )
    mat2    = rbind( y, x )
    
    out = rbind( as.numeric(mat1), as.numeric(mat2) )
    
    return( out )
    }


#   z   complex vector
#   returns absolute value squared   = a real number
abs2 <- function( z )
    {
    Re(z)^2 + Im(z)^2
    }

#   returns "inner product" of 2 complex numbers = a real number
#
#   this is equal to  (z1*Conj(z2) + Conj(z1)*z2)
#   and is 2 * the standard inner product

iprod <- function( z1, z2 )
    {
    2 * ( Re(z1)*Re(z2) + Im(z1)*Im(z2) )
    }

    
#   u0, u1  unit vectors of the same dimension.  The unit property is not verified
#   tau     vector of interpolating numbers, usually in [0,1]
#
#   returns a length(tau) x length(u0) matrix, with interpolated unit vectors in the rows

slerp <- function( u0, u1, thetamax=pi/36, tau=NULL )
    {
    theta   = anglebetween( u0, u1, unitized=TRUE )
    
    if( is.na(theta) )  return(NULL)
    
    if( theta == pi )
        {
        #   antipodal points
        log_level( ERROR, "u0 and u1 are antipodal." )
        return(NULL)
        }
        
        
    if( is.null(tau) )
        {
        #   compute tau from thetamax
        if( thetamax <= 0 )
            {
            log_level( ERROR, "thetamax = %g is invalid.", thetamax )
            return(NULL)
            }
        
        if( 0 < theta )
            {
            k   = ceiling( theta / thetamax )
            tau = (0:k) / k
            }
        else
            tau = 0
        }
        
    ok  = is.numeric(tau)  &&  0 < length(tau)
    if( ! ok )
        {
        log_level( ERROR, "tau is invalid. It must be a non-empty numeric vector." )
        return(NULL)
        }
    
    m       = length( u0 )
    count   = length(tau)
    
    if( theta == 0 )
        {
        #   u0 and u1 are equal
        out = matrix( u0, nrow=count, ncol=m, byrow=TRUE )
        return( out )
        }

    #   weight is count x 2
    weight  = cbind( sin( (1-tau)*theta ), sin( tau*theta ) ) / sin(theta)
        
    #   u0u1 is  2 x m
    u0u1    = rbind( u0, u1 )

    #   out is count x m
    out = weight %*% u0u1
    
    rownames(out)   = as.character( tau )

    return( out )
    }
    
    
    
#   vec1 and vec2     non-zero vectors of the same dimension
#
anglebetween  <-  function( vec1, vec2, unitized=FALSE, eps=5.e-14 )
    {
    if( length(vec1) != length(vec2) )
        {
        log_level( WARN, "length(vec1)=%d  !=  %d=length(vec2)", length(vec1), length(vec2) )
        return(NA_real_)
        }

    q   = sum( vec1*vec2 )

    if( ! unitized )
        {
        len1    = sqrt( sum(vec1^2) )
        len2    = sqrt( sum(vec2^2) )     #;    print( denom )

        denom   = len1 * len2

        if( abs(denom) < eps )    return(NA_real_)

        q   = q / denom  #; print(q)
        }

    if( abs(q) < 0.99 )
        {
        #   the usual case uses acos
        out = acos(q)
        }
    else
        {
        #   use asin instead
        if( ! unitized )
            {
            vec1   = vec1 / len1
            vec2   = vec2 / len2
            }

        if( q < 0 ) vec2 = -vec2

        d   = vec1 - vec2
        d   = sqrt( sum(d*d) )

        out = 2 * asin( d/2 )

        if( q < 0 ) out = pi - out
        }

    return(out)
    }
    
    
    
findRunsTRUE <- function( mask, circular=TRUE )
    {
    #   put sentinels on either end, to make things far simpler
    dif = diff( c(FALSE,mask,FALSE) )
    
    start   = which( dif ==  1 )
    stop    = which( dif == -1 )
    
    if( length(start) != length(stop) )
        {
        log_level( FATAL, "Internal error.  length(start)=%d != %d=length(stop)",
                                    length(start), length(stop) )
        return(NULL)
        }
        
    stop    = stop - 1L
    
    if( circular  &&  2<=length(start) )
        {
        m   = length(start)
        
        if( start[1]==1  &&  stop[m]==length(mask) )
            {
            #   merge first and last
            start[1]    = start[m]
            start   = start[ 1:(m-1) ]
            stop    = stop[ 1:(m-1) ]
            }
        }
        
    return( cbind( start=start, stop=stop ) )
    }
        
    