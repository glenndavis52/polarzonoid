

p.arcmatcolnames    = c( 'center', 'length' )


#   arcmat  an Nx2 matrix, with data on N disjoint arcs
#           in each row, the first number is the center of the arc
#           the second number is the length of the arc
#
#   m       dimension of space containing zonoid, it must be odd and >= 3
#
#   gapmin  for verifying that arcmat[] has disjoint arcs
#           the default gapmin=0 allows abutting arcs
#           if gapmin<0 some overlap is allowed.  If gapmin == -Inf, then test is skipped
#
#   returns a list with these items:
#       *) arcmat   the given arcs
#       *) point    the corresponding point on the polar zonoid, defined by trig polynomials
#       *) normal   outward normal of a supporting plane at point


boundaryfromarcs <- function( arcmat, n=NULL, gapmin=0 )
    {
    arcmat   = prepareNxM( arcmat, 2 )
    if( is.null(arcmat) )    return(NULL)

    #   check individual lengths
    neg = (arcmat[ ,2] < 0)
    if( any(neg) )
        {
        log_level( ERROR, "%d of %d arcs have negative length.", sum(neg), length(neg) )
        return(NULL)
        }

    if( gapmin != -Inf )
        {
        #   validation check
        gm  = gapminimum( arcmat )

        if( gm < gapmin )
            {
            log_level( ERROR, "Arcs are not disjoint. minimum gap = %g < %g.", gm, gapmin )
            return(NULL)
            }
        }

    ltotal  = sum( arcmat[ ,2] )

    improper    = ltotal==0  ||  ltotal==2*pi

    if( improper )
        nmin    = 0L
    else
        nmin    = nrow(arcmat)

    if( is.null(n) )    n   = nmin

    #   check that n is large enough
    if( n < nmin )
        {
        log_level( ERROR, "Argument n=%d is invalid, because it is < %d.", n, nmin )
        return(NULL)
        }


    m   = 2L*n + 1L

    if( improper )
        {
        #   return the 'south pole' or 'north pole'

        out     = numeric(m)
        out[m]  = ltotal

        #out$normal  = numeric(m)
        #out$normal[m]   = ltotal/pi - 1

        return(out)
        }

    #   since improper is FALSE,  n > 0

    point   = numeric(0)

    theta   = endpointsfromarcs( arcmat )

    for( j in 1:n )
        {
        cossin  = numeric(2)

        for( k in 1:nrow(arcmat) )
            {
            thetalo = j*( theta[k,1] )
            thetahi = j*( theta[k,2] )

            xypair  =  c( sinmy(thetahi) - sinmy(thetalo), -cosmy(thetahi) + cosmy(thetalo) ) / j

            cossin  = cossin + xypair
            }

        point = c( point, cossin )
        }

    #   append the total length
    point[m]  = ltotal

    names(point) = boundarynames(m)

    out = point

    if( FALSE )
    {
    if( m <= 7 )
        {
        #   use the known gradient of the known boundary function
        #   this is more reliable
        out$normal  = boundarynormal( point )
        }
    else if( 2*length(theta)+1 == length(point) )
        {
        #   use the jacobian at theta, and then SVD
        #   this can fail if the rcond is too small
        out$normal  = supportingnormal( theta, point )
        }
    else
        #   cannot be determined
        out$normal  = NA_real_
    }


    #   out$bvalue  = boundaryfunction(point)

    return( out )
    }



#   theta   n x 2 matrix with thetamin and thetamax in the rows
#           these are the endpoints of the arcs that map to the boundary point
#
#   point   point on the boundary of the zonoid
#           this is only used to ensure that the normal is outward-pointing
#
#   tol     tolerance for smallest singular value
#
#   returns 2*n+1 normal to the zonoid at point, unitized
#
#   if the boundary is not differentiable at the corresponding point on the boundary
#   the function returns NA_real_

supportingnormal <- function( theta, point, tol=5.e-14 )
    {
    n   = nrow(theta)

    #   the mapping is from 2*n parameters to (2*n + 1) parameters

    #   make jacobian matrix that is 2*n+1  x  2*n
    m   = 2L*n + 1L

    if( length(point) != m )
        {
        log_level( ERROR, "Internal error.  length(point) = %d  !=  %d = 2*nrow(theta)+1.", length(point), m )
        return(NULL)
        }

    jac = array( 0, dim=c(m,2*n) )

    for( i in 1:n )
        {
        for( k in 1:n )
            {
            thetalo = i * theta[k,1]
            thetahi = i * theta[k,2]

            jac[ 2*i-1, 2*k-1 ] = -cosmy(thetalo)
            jac[ 2*i-1, 2*k ]   =  cosmy(thetahi)
            jac[ 2*i,   2*k-1 ] = -sinmy(thetalo)
            jac[ 2*i,   2*k ]   =  sinmy(thetahi)
            }
        }

    #   now the last row
    for( k in 1:n )
        {
        jac[ m, 2*k-1]  = -1
        jac[ m, 2*k  ]  = +1
        }

    #print( jac )

    #   use svd() to find a unit vector that is orthogonal to all the columns
    res = base::svd( jac, nu=m, nv=0 )  #; print( res )

    rcond   = res$d[m-1] / res$d[1]

    if( rcond <= tol )
        {
        log_level( WARN, "rcond = %g <= %g.  Supporting normal is set to NA.", rcond, tol )
        return(NA_real_)
        }

    normal  = res$u[ ,m]

    testdot = sum( normal * point ) - normal[m]*pi
    if( testdot < 0 )
        #   change from inward to outward
        normal  = -normal

    names(normal)   =   names(point)

    return( normal )
    }



#   m   an odd number

#   returns a character vector of length m

boundarynames <- function( m )
    {
    n   = (m-1)/2   # of arcs

    out = character(0)

    for( i in seq_len(n) )
        out = c( out, paste( c('x','y'), i, sep='' ) )

    out = c( out, 'L' )

    return( out )
    }



#   arcmat  an Nx2 matrix, with data on N arcs
#           in each row, the first number is the center of the arc
#           the second number is the length of the arc
#
#   n       the given set of arcs is taken to be in \eqn{A_n}
#
#   returns corresponding point on the polar zonoid, then projected
#   onto unit sphere centered at the center of the zonoid = (0,0,0,0, ... , pi)

spherefromarcs <- function( arcmat, n=NULL, gapmin=0 )
    {
    p   = boundaryfromarcs( arcmat, n=n, gapmin=gapmin )

    if( is.null(p) )  return(NULL)

    p   = as.numeric(p)

    return( spherefromboundary(p) )
    }



#   arcmat  an Nx2 matrix, with data on n0 arcs.  n0 = nrow(arcmat)
#           in each row, the first number is the center of the arc
#           the second number is the length of the arc
#
#   n       the given set of arcs is taken to be in \eqn{A_n}
#
#   returns a list with these items:
#       u       the point on the sphere that arcmat maps to - a unit vector
#       tangent 2n+1  x  2n0 matrix.  The jacobian of the map at arcmat.
#               the 2*n0 columns are tangent to the stratum A_n0
#       normal  2n+1 x 2*(n-n0) matrix.  The columns are an orthonormal basis
#               for the normal space to the stratum A_n0 inside A_n
#
#   if necessary, one of the normals is flipped so that
#       det( cbind(u,tangent,normal) ) > 0

spherefromarcs_plus <- function( arcmat, n=NULL, gapmin=5.e-10 )
    {
    arcmat   = prepareNxM( arcmat, 2 )
    if( is.null(arcmat) )    return(NULL)

    if( gapmin <= 0 )
        {
        log_level( WARN, "gapmin = %g <= 0 is invalid.", gapmin )
        return(NULL)
        }

    p   = boundaryfromarcs( arcmat, n=n, gapmin=gapmin )

    if( is.null(p) )  return(NULL)

    p   = as.numeric(p)

    u   = spherefromboundary(p)

    n0    = nrow(arcmat)

    if( is.null(n) )    n = n0

    #   make jacobian matrix that is 2*n+1  x  2*n0
    theta   = endpointsfromarcs( arcmat )

    m   = 2L*n + 1L

    jac = array( 0, dim=c(m,2*n0) )

    for( i in 1:n )
        {
        for( k in 1:n0 )
            {
            thetalo = i * theta[k,1]
            thetahi = i * theta[k,2]

            jac[ 2*i-1, 2*k-1 ] = -cosmy(thetalo)
            jac[ 2*i-1, 2*k ]   =  cosmy(thetahi)
            jac[ 2*i,   2*k-1 ] = -sinmy(thetalo)
            jac[ 2*i,   2*k ]   =  sinmy(thetahi)
            }
        }

    #   now the last row
    for( k in 1:n0 )
        {
        jac[ m, 2*k-1]  = -1
        jac[ m, 2*k  ]  = +1
        }

    #   the columns of jac are tangent to the boundary of the zonoid
    #   to get the tangent vectors to the sphere,
    #   project columns of jac onto the plane normal to u

    ujac    = as.numeric( u %*% jac )

    tangent = jac - u %o% ujac

    if( n == n0 )
        {
        #   set normal to a matrix with 0 columns
        out = list( u=u, tangent=tangent, normal=matrix(0,nrow=m,ncol=0 ) )
        return( out )
        }

    mat = cbind( u, tangent )   # mat has 2*n0+1 columns

    #   use svd() to find 2*(n-n0) unit vectors that are orthogonal to all the columns of mat
    res = base::svd( mat, nu=m, nv=0 ) # ; print( res )

    rcond   = res$d[ length(res$d) ] / res$d[1]

    tol = 1.e-14

    if( tol < rcond )
        {
        normal  = res$u[ , (m - 2*(n-n0) +1):m ]

        #   check the orientation and make it + if necessary
        mat = cbind( mat, normal )

        #   mat[,] is now square, m x m
        if( base::determinant( mat )$sign < 0 )
            #   reverse the sign of last normal vector
            normal[ ,ncol(normal)]  = -normal[ ,ncol(normal)]
        }
    else
        {
        normal  = NA_real_
        log_level( WARN, "rcond = %g <= %g.  normal[,] is set to NA.", rcond, tol )
        }

    out = list( u=u, tangent=tangent, normal=normal )

    return( out )
    }




#   p   a point known to be on the boundary of the zonoid
#       p should come from boundaryfromsphere()  OR   boundaryfromarcs()
#       if p has even dimension, it is interpreted as being on the equator of the zonoid of one higher dimension

spherefromboundary <- function( p )
    {
    m   = length(p)

    if( (m %% 2) == 1 )
        {
        #   odd is the usual cases
        x       = p
        x[m]    = x[m] - pi
        }
    else
        {
        #   if even, interpret p as on the equator, just append 0
        x   = c( p, 0 )
        }

    x   = as.numeric( x )

    #   x now has odd length
    return( unitize( x ) )
    }



#   x   a non-zero vector.  If even-dimensional, then 0 is appended to make it odd-dimensional
#       x is usually on the sphere, but if not it is re-unitized to make it so

boundaryfromsphere <- function( x )
    {
    u   = unitize(x)

    if( any( is.na(u) ) )
        {
        log_level( ERROR, "argument x is 0, which is invalid." )
        return(NULL)
        }

    if( (length(u) %% 2L) == 0L )
        #   x is even-dimensional, so interpret as being on equator.  Append a 0.
        u   = c( u, 0 )

    #   u is now odd-dimensional
    m   = length(u)

    n   = (m-1L) / 2L       # n is the number of arcs

    if( 3 < n )
        {
        log_level( ERROR, "The number of arcs is n=%d, but in this package version it must be 1, 2, or 3.", n )
        return(NULL)
        }

    center  = c( numeric(m-1), pi )

    if( u[m] == 0  &&  n < 3 )
        {
        #   in this special case, there is a closed-form expression for the root
        if( n == 1 )
            # only 1 arc, n=1 and m=3
            root   = 2
        else
            {
            #   2 arcs, n=2, and m=5
            beta    = min( sum(u[1:2]^2), 1 )   # do not allow even a tiny bit more than 1

            root   = 4 / (sqrt(1 + 3*beta) + sqrt(1 - beta))
            }

        out = center + root*u

        names(out)  = boundarynames(m)
        
        log_level( INFO, "computed root=%g using closed form expression.  f(root) = %g.",
                                root, boundaryfunction(out) )

        return( out )
        }


    #   the usual case, use uniroot()
    myfun <- function( s )  { boundaryfunction( center + s*u ) }

    #   the root must be in the interval [rootmin,rootmax]
    if( n == 1 )
        rootmin = 2 - 0.1   # expand a little
    else
        rootmin = 1

    rootmax = 2.2 # pi

    #   expand the interval just a bit
    #  rootmin = rootmin - 0.1

    # cat( "myfun(endpoints) = ", myfun(rootmin), "  ", myfun(rootmax), '\n' )


    #   in the next call, tol is much smaller than the default
    res = try( stats::uniroot( myfun, c(rootmin,rootmax), tol=.Machine$double.eps^0.75, extendInt="upX" ), silent=FALSE )   #; print( str(res) )

    if( ! inherits(res,"try-error" ) )    # class(res) == "try-error" )
        {
        #   success
        root   = res$root

        log_level( INFO, "computed transverse root = %g using stats::uniroot().  f(root) = %g.",
                                root, res$f.root )
        }
    else
        {
        #   failure
        root    = NA_real_

        log_level( WARN, "stats::uniroot() failed. myfun() at endpoints.  myfun(%g)=%g   myfun(%g)=%g.",
                            rootmin, myfun(rootmin), rootmax, myfun(rootmax) )

        if( FALSE )
            {
            # print( dumpalong( seq( rootmin, rootmax, len=101 ), u ) )
            }

        #  cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
        return( NULL )
        }

    if( is.finite(root)  &&  root <= pi )
        {
        #   the root is probably valid
        out = center + root*u

        names(out)  = boundarynames(m)

        return( out )
        }


    #   the computed transverse root is probably bogus
    #   the true root is probably non-transverse,  where the x-axis is tangent to the curve
    #   find the maximum on the interval [rootmin,pi]
    res = stats::optimize( myfun, interval=c(rootmin,pi), maximum=TRUE, tol=.Machine$double.eps^0.50 )

    tol = 5.e-12
    
    if( -tol < res$objective )
        {
        #  success
        root    = res$maximum

        log_level( INFO, "computed non-transverse root = %g using stats::optimize().    f(root) = %g > %g.",
                                    root, res$objective, -tol )

        out = center + root*u

        names(out)  = boundarynames(m)
        
        return( out )
        }

        
    #   failure
    log_level( ERROR, "stats::optimize() failed to compute root.  argmax = %g and f(argmax) = %g >= %g.",
                                    res$maximum, res$objective, -tol )

    return( NULL )
    }


dumpalong <- function( svec, u )
    {
    m   = length(u)

    center  = c( numeric(m-1), pi )

    count   = length(svec)

    lhs     = complex(count)
    rhs     = numeric(count)
    value   = numeric(count)

    for( k in 1:count )
        {
        x   = center + svec[k]*u
        z   = complexfromrealvector( x )
        L   = x[ length(x) ]

        res = boundaryparts( z, L )
        lhs[k]  = res$lhs
        rhs[k]  = res$rhs
        value[k] = abs(res$lhs) - res$rhs
        }

    out = data.frame( svec=svec )
    out$lhs = lhs
    out$rhs = rhs
    out$value   = value

    return(out)
    }



#   p   a non-zero vector, on the boundary of the zonoid, with length 2 to 7
#       p should come from boundaryfromarcs() or boundaryfromsphere().
#       verification of this boundary property is minimal.
#       if the length is even, a pi is appended to  put the point on the "equator"
#       so the number of arcs must be 0, 1, 2, or 3
#
#   tol if the distance from p to a sub-stratum of the boundary is less than tol
#       then take p to actually *be* in this substratum and compute the arcs from there
#
#   returns a matrix arcmat.  Or NULL in case of error.

arcsfromboundary <- function( p, tol=5.e-9 )
    {
    if( (length(p) %% 2) == 0 )   p = c( p, pi )    #  put the point on the "equator"

    m   = length(p)

    n   = (m-1L) / 2L       # n is the number of arcs

    ok  = n %in% (0:3)
    if( ! ok )
        {
        log_level( ERROR, "The number of arcs is %d, which is invalid in this version of the package.", n )
        return(NULL)
        }

    L   = p[m]

    if( FALSE  &&  m == 5 )
        {
        #   test whether the 2 arcs compute arcs will be abutting
        test    = sqrt( sum( p[1:2]^2 ) ) - crd(p[5])

        if( abs(test) <= tol )
            {
            #   just return a single arc, ignore the double angle parts 3 & 4
            p   = p[ c(1,2,5) ]
            m   = length(p)
            n   = (m-1L) / 2L
            }
        }

    normal  = boundarynormal( p, tol=tol )

    if( is.null(normal)  ||  is.na( normal[1] ) )
        {
        log_level( WARN, "The boundary normal cannot be computed.\n" )
        return(NULL)
        }

    #  convert normal to coefficients of a trig polynomial
    ablist  = trigpolyfromvector( normal )

    #   the length of theta is always even, and possibly 0
    theta   = trigpolyroot( ablist ) #; print( theta )

    out = arcsfromendpointsandpoly( theta, ablist )

    #  cat( "normal = ", normal, "   theta = ", theta, '\n' )

    return( out )
    }



#   u   a non-zero vector, with length 2 or more
#   if the length is even, a 0 is appended to make the length odd
#   typically it is a point in the unit sphere, but if not it is unitized automatically
#
#   returns a matrix arcmat, to invert spherefromarcs()
#   the number of arcs is floor(length(u)/2) or less.
#   if u is the "north pole", a single improper arc of length 2*pi is returned


arcsfromsphere <- function( u )
    {
    p   = boundaryfromsphere( u )

    if( is.null(p) )    return(NULL)

    arcmat  = arcsfromboundary( p )

    return( arcmat )
    }



    
    
#   z   complex vector of length n = # of arcs
#   L   total length of arcs
#
#   return list with:
#       lhs     complex number
#       rhs     real number

boundaryparts <- function( z, L )
    {
    n   = length(z)

    if( n == 1 )
        {
        lhs = z[1]
        rhs = 2*sinmy(L/2)
        }
    else if( n == 2 )
        {
        #   lhs is complex
        lhs = 2*sinmy(L/2)*z[2] - cosmy(L/2)*z[1]^2

        # rhs is real
        # if x is *on* the boundary, rhs>=0 and is 0 iff 2 arcs degenerate to 1
        rhs = 4*sinmy(L/2)^2 - abs2(z[1])

        #x1  = x[1]
        #y1  = x[2]
        #x2  = x[3]
        #y2  = x[4]

        #lhs_old = sqrt( (2*sinmy(L/2)*x2 - cosmy(L/2) * (x1^2 - y1^2))^2  +   (2*sinmy(L/2)*y2 - cosmy(L/2) * 2*x1*y1)^2  )
        #rhs_old = 4*sinmy(L/2)^2 - (x1^2 + y1^2)

        #old = lhs_old - max( rhs_old, 0 )
        }
    else if( n == 3 )
        {
        z1  = z[1]
        z2  = z[2]
        z3  = z[3]

        sine    = sinmy(L/2)
        cosine  = cosmy(L/2)
        z1abs2  = abs2(z1)

        #   lhs and rhs were computed with help from Macaulay2

        lhs = 12*(4*sine^2 - z1abs2)*z3  +  (8*cosine^2 - z1abs2 + 4)*z1^3  -  48*sine*cosine*z1*z2  +  12*z2^2*Conj(z1)    # Conj(z2) missing

        #   the rhs is actually real
        #   if x is *on* the boundary, rhs>=0 and is 0 iff 3 arcs degenerate to 2
        rhs = 6*sine*( z1abs2^2  -  8*z1abs2  -  4*abs2(z2) )  +  12*cosine*iprod( z1^2, z2 )  +  96*sine^3         # z3 is missing
        }

    out = list( lhs=lhs, rhs=rhs )

    return( out )
    }



#   x   a point in R^(2*n+1)
#       m=length(x) must be positive and odd

boundaryfunction <- function( x )
    {
    m   = length(x)

    ok  = (1 <= m)  &&  ((m %% 2L)==1L)

    if( ! ok )
        {
        log_level( ERROR, "length(x) = %d is invalid.", m )
        return( NA_real_ )
        }

    #   n is the number of arcs
    n   = as.integer( (m-1)/2 )

    Lgiven  = x[m]

    #   clamp Lgiven to [0,2*pi], and use L in calculations
    L   = min( max( Lgiven, 0 ), 2*pi )


    if( FALSE  &&  L == pi )
        {
        #   a big simplification
        even    = 2L *( 1:n )
        odd     = even - 1L

        zmod2   = as.numeric( x[odd]^2 + x[even]^2 )

        if( n == 3 )
            {
            zmod = sqrt( zmod2 )

            out = (1/4)*zmod[1]^3  +  (1/2)*zmod[2]^2  +  zmod[3] - 2

            return( out )
            }

        out = sqrt(zmod2[n]) - 2

        if( 2 <= n )
            {
            for( k in (n-1):1 )
                out = zmod2[k]^((n+1-k)/2) + 2*out
            }

        return( out )
        }

    if( n == 0 )
        {
        #   then m=length(x) must be 1
        #   the zonoid is the interval [0,2*pi]
        out = x*(x - 2*pi)/(2*pi)
        }
    else if( n == 1 )
        {
        L   = Lgiven

        out = as.numeric( sqrt( sum( x[1:2]^2 ) ) )

        if( L <= 0 )
            out = out - L
        else if( L <= 2*pi )
            out = out - crd(L)
        else
            out = out + (L - 2*pi)
        }
    else if( n %in% 2:3 )
        {
        z   = complexfromrealvector( x )

        res = boundaryparts( z, L )

        out = abs(res$lhs) - res$rhs    # max(rhs,0)

        if( L != Lgiven )
            #   Lgiven is outside the prescribed interval, so add penalty
            out = out + (Lgiven - L)^2

        names(out)  = NULL
        }
    else
        {
        log_level( WARN, "The number of arcs n=%d is too large. n must be 0,1,2, or 3.  Returning NA.", n )
        return( NA_real_ )
        }

    return( out )
    }

#   x   a point in R^(2n+1) where n is 0, 1 or 2 or 3
#       x is also assumed to be on the boundary of the polar zonoid, but this is not verified.
#       x should really come from boundaryfromarcs() or boundaryfromsphere()
#
#   tol if the distance from x to a sub-stratum of the boundary is less than tol
#       then take x to actually *be* in this substratum and return the dimension of that substratum
#
#   returns the dimension of the stratum that contains x.
#   this is 2*|arcs|

stratumdimension <- function( x, tol=5.e-9 )
    {
    m   = length(x)

    ok  = (1 <= m) && (m <= 7)  &&  ((m %% 2L)==1L)

    if( ! ok )
        {
        log_level( ERROR, "length(x) = %d is invalid.", m )
        return( NA_integer_ )
        }

    if( m == 1 )    return( 0L )    # trivial case, x must be 0 or 2*pi.  x is one of the "poles".

    z   = complexfromrealvector( x )

    L   = x[m]      # sum of all arc lengths

    #   n is the number of arcs
    n   = length(z)

    for( k in 1:n )
        {
        res = boundaryparts( z[1:k], L )

        #cat( "k=", k, "  rhs=", res$rhs, '\n' )

        #   the rhs is the "distance" in the stratification
        if( res$rhs <= tol )    return( 2L * (k - 1L) )
        }

    return( 2L * n )
    }





#   x   a point in R^(2n+1) where n is 0, 1 or 2 or 3
#       x is also assumed to be on the boundary of the polar zonoid, but this is not verified.
#       x should really come from boundaryfromarcs() or boundaryfromsphere()
#       A fundamental property of convex bodies is that every boundary point
#       has a supporting normal
#
#   tol if the distance from x to a sub-stratum of the boundary is less than tol
#       then take x to actually *be* in this substratum and compute the normal there
#
#   returns a supporting normal to boundary at x, unitized.
#   in the generic case, this normal is the gradient of the boundaryfunction() at x, and then unitized

boundarynormal <- function( x, tol=5.e-9 )
    {
    out = boundarygradient( x, tol )

    if( is.na(out[1]) ) return( out )

    return( unitize(out) )
    }


boundarygradient <- function( x, tol=5.e-9 )
    {
    m   = length(x)

    ok  = (1 <= m) && (m <= 7)  &&  ((m %% 2L)==1L)

    if( ! ok )
        {
        #   log_level( ERROR, "length(x) = %d is invalid.", m )
        return( rep(NA_real_,m) )
        }

    L   = x[m]      # sum of all arc lengths

    #n   = as.integer( (m-1)/2 )

    #   test for being close to a substratum
    #   n is the number of arcs
    n   = stratumdimension( x, tol=tol ) / 2L

    if( is.na(n) )  return( rep(NA_real_,m) )

    out = numeric(m)

    if( n == 0 )
        {
        #   L is either 0 or 2*pi
        out[m]  = L/pi - 1  # -1 or +1
        }
    else if( n == 1 )
        {
        #   the gradient here is simple
        out[1:2] = unitize( x[1:2] )

        out[m]  = -cosmy( L/2 )
        }
    else if( n == 2 )
        {
        #---  using complex numbers
        z       = complexfromrealvector( x[1:4] )       # z[] has length 2
        sine    = sinmy(L/2)
        cosine  = cosmy(L/2)

        #   gradient of the LHS
        lhs     = cosine*z[1]^2  -  2*sine*z[2]
        abslhs  = abs(lhs)
        if( abslhs == 0 )
            {
            # this means that 2 arcs are abutting.  There is a crease on the boundary
            # this should not happen if tol > 0,
            # since 2 abutting arcs should have forced n=1 in the substrata test
            log_level( WARN, "LHS is 0, so given point x = %s is at a crease.", paste(x[1:4],collapse=',') )            
            return( rep(NA_real_,m) )
            }

        zL  = -sine*z[1]^2/2 - cosine*z[2]

        jac2x5  = cbind( real2x2Nfromcomplex( c( 2*cosine*z[1], -2*sine ) ), realfromcomplexvector(zL) )  #; print( jac2x5 )

        jac1x2  = realfromcomplexvector(lhs)/abslhs     #; print( jac1x2 )

        gradlhs = as.numeric(  jac1x2  %*%  jac2x5  )

        #   gradient of the RHS
        gradrhs = c( realfromcomplexvector( c( -2*z[1], 0 ) ), 2*sinmy(L) )     # length is 5

        grad    = gradlhs - gradrhs

        out = grad

        if( FALSE )
            {
            #   using real numbers
            x1  = x[1]
            y1  = x[2]
            x2  = x[3]
            y2  = x[4]

            C   = 2*sinmy(L/2)*x2 - cosmy(L/2)*(x1^2 - y1^2)
            D   = 2*sinmy(L/2)*y2 - cosmy(L/2)*(2*x1*y1)

            G   = sqrt( C^2 + D^2 )
            if( G == 0 )
                {
                # this means that 2 arcs are abutting.  There is a crease on the boundary
                # this should not happen if tol > 0,
                # since 2 abutting arcs should have forced n=1 in the substrata test
                log_level( WARN, "LHS is 0, so given point x = %s is at a crease.", paste(x,collapse=',') )
                return( rep(NA_real_,m) )
                }

            out[1]  = (2*C*2*x1  +  2*D*2*y1) * (-cosmy(L/2))

            out[2]  = (-2*C*2*y1  +  2*D*2*x1) * (-cosmy(L/2))

            out[3]  = 2*C*2*sinmy(L/2)

            out[4]  = 2*D*2*sinmy(L/2)

            out[5]  = (2*x2^2 + 2*y2^2 - 0.5*(x1^2 + y1^2)^2)*sinmy(L)  -  2*(x2*(x1^2 - y1^2) + y2*2*x1*y1)*cosmy(L)

            #   add correction from chain rule
            out = (0.5 * 1/G) * out

            #   and subtract the gradient of the second term
            out[1]  = out[1] + 2*x1
            out[2]  = out[2] + 2*y1
            out[5]  = out[5] - 2*sinmy(L)

            #   compare out and grad
            #cat( "out = ", out, '\n' )
            #cat( "grad = ", grad, '\n' )
            #cat( "delta = ", max( abs(out-grad) ), '\n' )
            }
        }
    else if( n == 3 )
        {
        #---  using complex numbers
        z       = complexfromrealvector( x[1:6] )       # z[] has length 3

        lhs     = boundaryparts( z, L )$lhs
        abslhs  = abs(lhs)
        if( abslhs == 0 )
            {
            # this means that 2 arcs are abutting.  There is a crease on the boundary
            # this should not happen if tol > 0,
            # since 2 abutting arcs should have forced n=1 in the substrata test
            log_level( WARN, "LHS is 0, so given point x = %s is at a crease.", paste(x[1:6],collapse=',') )
            return( rep(NA_real_,m) )
            }
            
        sine    = sinmy(L/2)
        cosine  = cosmy(L/2)
        abs2z1  = abs2( z[1] )
        xy1     = realfromcomplexvector( z[1] )

        #   gradient of the LHS
        pz1 =  8*cosine^2*3*z[1]^2  -  abs2z1*3*z[1]^2  +  12*z[1]^2  - 48*sine*cosine*z[2]   # more to be added later

        pz2 = -48*sine*cosine*z[1]  +  24*Conj(z[1])*z[2]

        pz3 = 48*sine^2 - 12*abs2z1

        pL  = 48*sine*cosine*z[3]  -  8*sine*cosine*z[1]^3  - 24*(cosine^2 - sine^2)*z[1]*z[2]

        #   start with just 2x2
        jac2x7  = real2x2Nfromcomplex( pz1 )

        jac2x7  = jac2x7  -  real2x2Nfromcomplex( 12*z[3] + z[1]^3 )  %*%  rbind( 2*xy1, 0 )

        jac2x7  = jac2x7  +  real2x2Nfromcomplex( 12*z[2]^2 )  %*%  matrix( c(1,0,0,-1), 2, 2 )

        #   now add 5 more columns
        jac2x7  = cbind( jac2x7, real2x2Nfromcomplex( c(pz2,pz3) ), realfromcomplexvector(pL) )

        jac1x2  = realfromcomplexvector(lhs)/abslhs     #; print( jac1x2 )

        gradlhs = as.numeric(  jac1x2  %*%  jac2x7  )           # length is 7


        #   gradient of the RHS
        zgrad   = c( -48*sine*2*z[1] + 6*sine*4*abs2z1*z[1] + 12*cosine*4*Conj(z[1])*z[2], -24*sine*2*z[2] + 12*cosine*2*z[1]^2, 0 )

        partL   = 144*sine^2*cosine  -  24*cosine*abs2z1  +  3*cosine*abs2z1^2  -  12*cosine*abs2(z[2])  -  6*sine*iprod(z[1]^2,z[2])

        gradrhs = c( realfromcomplexvector(zgrad), partL )      # length is 7

        grad    = gradlhs - gradrhs

        out = grad
        }

    return(out)
    }





#
#   direction   RxM matrix, with the R directions in the rows, direction (0,0) is invalid
#
#   returns a data.frame with R rows and these columns:
#       direction   the given matrix of directions
#       value       the value of the support function of the zonoid, in the given direction
#       argmax      the point on the boundary of the zonoid where the max is taken
#       arcs        number of arcs for argmax

support <- function( direction )
    {
    ok  = is.numeric(direction)  &&  0<length(direction)

    if( ! ok )
        {
        log_level( ERROR, "argument direction is invalid." )
        return(NULL)
        }

    if( ! is.matrix(direction) )
        dim(direction)  = c(1,length(direction))

    if( (ncol(direction) %% 2) == 0 )
        #   the number of columns is even, so append a column of 0s
        direction   = cbind( direction, 0 )

    m       = ncol(direction)   # m is odd

    count   = nrow(direction)

    value   = rep( NA_real_, count )
    argmax  = matrix( NA_real_, count, m )
    arcs    = rep( NA_integer_, count )

    #  convert normal to coefficients of a trig polynomial
    n       = (m-1L)/2L
    even    = 2L * ( 1:n )
    odd     = even - 1L

    for( i in 1:count )
        {
        normal  = direction[i, ]

        if( all(normal==0) )    next

        ablist  = list( a0=normal[m], a=normal[odd], b=normal[even] )   #; print( ablist )

        #   length(theta) is guaranteed to be even
        theta   = trigpolyroot( ablist )   #; print( theta )

        arcmat = arcsfromendpointsandpoly( theta, ablist )  #; print( arcmat )

        p   = boundaryfromarcs( arcmat, n=n )       #;  print(p)

        #z   = res$point

        value[i]    = sum( normal * p )
        argmax[i, ] = p
        arcs[i]     = length(theta)/2L      #; nrow( arcmat )
        }

    rnames  = rownames(direction)
    if( is.null(rnames)  ||  anyDuplicated(rnames)!=0 ) rnames = 1:count

    out = data.frame( row.names=rnames )
    out$direction   = direction
    out$value       = value
    out$argmax      = argmax
    out$arcs        = arcs

    return(out)
    }

