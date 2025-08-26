
library( polarzonoid )

#library( logger )
#log_threshold( TRACE, "polarzonoid" )


options( width=160 )

#   n       number of arcs
#   count   # of random points on the unit sphere, inside R^(2*n+1)
#   tol     round-trip tolerance
#
#   maps from sphere to boundary of zonoid, and back again
#   return TRUE or FALSE

testinversionsphere  <- function( n, count=1000, tol=5.e-11 )
    {
    set.seed(0)

    m   = 2L*n + 1L

    delta   = numeric(count)

    bvalue  = numeric(count)

    #   make count random points on the sphere
    for( i in 1:count )
        {
        u   = polarzonoid:::unitize( rnorm( m ) )

        x   = boundaryfromsphere( u )

        if( is.null(x) )
            {
            return(FALSE)
            }

        bvalue[i]   = polarzonoid:::boundaryfunction( x )

        up  = spherefromboundary( x )

        delta[i]    = max( abs(u-up) )
        }

    maxdelta = max(delta)

    cat( "testinversionsphere().  n =", n,  "   delta max = ", maxdelta , '\n', file=stderr() )

    ran = range(bvalue)
    cat( "testinversionsphere().  n =", n,  "   bvalue range = ", ran, '\n', file=stderr() )


    if( tol < maxdelta )
        {
        mask    = tol < delta

        mess    = sprintf( "ERR.  testinversionsphere().  n=%d. %g of %d sphere inversions are invalid.  tol=%g\n",
                        n, sum(mask), length(mask), tol )
        cat( mess, file=stderr() )

        return(FALSE)
        }

    if( tol < max( abs(ran) ) )
        {
        mask    = tol < abs(bvalue)

        mess    = sprintf( "ERR.  testinversionsphere().  n=%d. %g of %d boundary values are invalid.  tol=%g\n",
                        n, sum(mask), length(mask), tol )
        cat( mess, file=stderr() )

        return(FALSE)
        }

    return(TRUE)
    }



#   nrad    radius of loop around torus "core"
#   count   # of points in loop
#   tol     round-trip tolerance
#
#   maps from points on loop in sphere to arcs, and back again
#   return TRUE or FALSE

testinversionsphereloop  <- function( rad2=0.1, count=360, tol=5.e-12 )
    {
    n   = 2L

    delta   = numeric(count)

    rad1    = sqrt( 1 - rad2^2 )

    #   make count random points on the sphere
    for( i in 0:(count-1) )
        {
        theta   = (i*360)/count     # theta is in degrees

        u   = c( rad1, 0, rad2 * polarzonoid:::cosdeg(theta), rad2 * polarzonoid:::sindeg(theta), 0 )

        arcmat  = arcsfromsphere( u )

        up  = spherefromarcs( arcmat )

        delta[i]    = max( abs(u-up) )

        # cat( "u=", u,  "   up=", up, '\n' )
        }

    maxdelta = max(delta)

    cat( "testinversionsphereloop().  n =", n,  "   delta max = ", maxdelta , '\n', file=stderr() )

    if( tol < maxdelta )
        {
        mask    = tol < delta

        mess    = sprintf( "ERR.  testinversionsphereloop().  n=%d. %g of %d sphere loop inversions are invalid.  tol=%g\n",
                        n, sum(mask), length(mask), tol )
        cat( mess, file=stderr() )

        warnings()

        return(FALSE)
        }

    return(TRUE)
    }


#   n       # of arcs in each set
#   count   # of random sets of arcs with n arcs
#   tol     # for the round-trip discrepancy

#   maps from arcs to boundary of zonoid, and back again
#   return TRUE or FALSE

testinversionarcs  <- function( n, count=1000, tol=5.e-9 )
    {
    set.seed(0)

    delta  = numeric(count)

    for( i in 1:count )
        {
        # cat( "---------------------\n" )

        arcmat1 = polarzonoid:::randomarcs( n )  #; print( arcmat1 )

        p   = boundaryfromarcs( arcmat1, gapmin=-Inf )
        if( is.null(p) )
            return(FALSE)

        #p   = res$point

        arcmat2 = arcsfromboundary( p ) #; print( arcmat2 )
        if( is.null(arcmat2) )
            {
            print( arcmat1 )
            return(FALSE)
            }

        d   = arcsdistance( arcmat1, arcmat2 )
        
        if( is.null(d) )
            {
            print( arcmat1 )
            print( arcmat2 )
            return(FALSE)
            }
            
        delta[i]    = d

        # cat( "delta = ", delta[i], '\n' )
        }

    maxdelta = max(delta)

    cat( "testinversionarcs().  n =", n,  "   delta max = ", maxdelta , '\n', file=stderr() )

    if( tol < maxdelta )
        {
        mask    = tol < delta

        mess    = sprintf( "ERR.  testinversionarcs().  n=%d. %g of %d arc inversions are invalid.  tol=%g\n",
                        n, sum(mask), length(mask), tol )
        cat( mess, file=stderr() )

        return(FALSE)
        }

    return( TRUE )
    }

#   return TRUE or FALSE
#
#   single arc   -->   boundary   -->   arcs
#
#   test whether round trip is accurate, within tol

testinversion_singlearc  <- function( count=1000, tol=5.e-11 )
    {
    set.seed(0)

    delta  = numeric(count)

    for( i in 1:count )
        {
        arcmat1 = polarzonoid:::randomarcs( 1 )

        #   set n=2 to force a point in 5D, which would normally come from 2 arcs and not 1
        #   but now comes from a single arc
        p = boundaryfromarcs( arcmat1, n=2L, gapmin=-Inf )
        if( is.null(p) )
            return(FALSE)

        # p   = res$point

        #   but p is on a non-differentiable point on the boundary
        #   and arcmat1p should have only one arc in it
        arcmat1p = arcsfromboundary( p )
        if( is.null(arcmat1p) )
            {
            print( arcmat1 )
            return(FALSE)
            }

        if( nrow(arcmat1p) != 1 )
            {
            cat( "row(arcmat1p) = ",nrow(arcmat1p), " but expected 1.\n" )
            return(FALSE)
            }

        delta[i]    = arcsdistance( arcmat1, arcmat1p )
        }

    maxdelta = max(delta)

    cat( "testinversion_singlearc().  ", "   delta max = ", maxdelta , '\n', file=stderr() )

    if( tol < maxdelta )
        {
        mask    = tol < delta

        mess    = sprintf( "ERR.  testinversion_singlearc().   %g of %d arc inversions are invalid.  tol=%g\n",
                         sum(mask), length(mask), tol )
        cat( mess, file=stderr() )

        return(FALSE)
        }

    return( TRUE )
    }





if( ! testinversionsphere(0,count=10) )  stop( "testinversionsphere(0) failed !" )

if( ! testinversionsphere(1) )  stop( "testinversionsphere(1) failed !" )

if( ! testinversionsphere(2) )  stop( "testinversionsphere(2) failed !" )

if( ! testinversionsphere(3,count=5000,tol=5.e-10) )  stop( "testinversionsphere(3) failed !" )



if( ! testinversionarcs(0,count=10) )  stop( "testinversionarcs(0) failed !" )

if( ! testinversionarcs(1,tol=5.e-12) )  stop( "testinversionarcs(1) failed !" )

if( ! testinversionarcs(2,tol=5.e-8) )  stop( "testinversionarcs(2) failed !" )

if( ! testinversionarcs(3,tol=5.e-7) )  stop( "testinversionarcs(3) failed !" )



if( ! testinversion_singlearc() )  stop( "testinversion_singlearc() failed !" )

if( ! testinversionsphereloop() )  stop( "testinversionsphereloop() failed !" )


cat( "Passed all inversion tests !\n", file=stderr() )
