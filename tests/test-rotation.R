
library( polarzonoid )

#library( logger )
#log_threshold( TRACE, "polarzonoid" )


options( width=160 )

#   count   # of random rotation 3x3 matrices  
#   tol     round-trip tolerance
#
#   maps from rotation to arcs, and back again
#   return TRUE or FALSE

testinversionrotation  <- function( count=1000, tol=5.e-11 )
    {
    set.seed(0)

    delta   = numeric(count)

    #   make count random points on the sphere
    for( i in 1:count )
        {
        quat    = polarzonoid:::unitize( rnorm( 4 ) )
        
        M   = polarzonoid:::rotationfromquat( quat )

        arcmat  = arcsfromrotation( M )

        if( is.null(arcmat) )
            {
            return(FALSE)
            }

        Mp  = rotationfromarcs( arcmat )

        delta[i]    = max( abs(M-Mp) )
        }

    maxdelta = max(delta)

    cat( "testinversionsphere().    maxdelta = ", maxdelta , '\n', file=stderr() )


    if( tol < maxdelta )
        {
        mask    = tol < delta

        mess    = sprintf( "ERR.  testinversionrotation().   %g of %d rotation inversions are invalid.  tol=%g. maxdelta=%g\n",
                        sum(mask), length(mask), tol, maxdelta )
        cat( mess, file=stderr() )

        return(FALSE)
        }

    return(TRUE)
    }





#   return TRUE or FALSE
testinversion_singlearc  <- function( count=1000, tol=5.e-11 )
    {
    set.seed(0)

    delta  = numeric(count)

    for( i in 1:count )
        {
        arcmat1 = polarzonoid:::randomarcs( 1, L=pi )


        M = rotationfromarcs( arcmat1  )
        if( is.null(M) )
            return(FALSE)

        # p   = res$point

        #   in the next call, tol is tricky
        #   tol is passed to arcsfromboundary()
        #   tol=5.e-9 sometimes yields 2 arcs, not the right substratum
        #   tol=5.e-8 leads to "boundary gradient cannot be computed"
        #   tol=1.e-7 and 1.e-6  are OK
        arcmat1p = arcsfromrotation( M , tol=1.e-7 )    
        if( is.null(arcmat1p) )
            {
            print( arcmat1 )
            return(FALSE)
            }
            
        if( nrow(arcmat1p) != 1 )
            {
            cat( "For test arcs i=", i, "   row(arcmat1p) = ",nrow(arcmat1p), " but expected 1.\n" )
            print( arcmat1 )
            print( arcmat1p )
            return(FALSE)
            }            
 
        delta[i]    = arcsdistance( arcmat1, arcmat1p )        

        if( 6 < delta[i] )
            {
            #   probably got the complementary arc, so take complement
            arcmat1p    = complementaryarcs( arcmat1p )
            delta[i]    = arcsdistance( arcmat1, arcmat1p ) 
            }
            
        # cat( sprintf( "delta[%d] = %g\n", i, delta[i] ) )                
        }

    maxdelta = max(delta)

    cat( "testinversion_singlearc().  ", "   maxdelta = ", maxdelta , '\n', file=stderr() )

    if( tol < maxdelta )
        {
        mask    = tol < delta

        mess    = sprintf( "ERR.  testinversion_singlearc().   %g of %d single arc inversions are invalid.  tol=%g\n",
                         sum(mask), length(mask), tol )
        cat( mess, file=stderr() )

        return(FALSE)
        }

    return( TRUE )
    }



if( ! testinversion_singlearc(count=1000) )  stop( "testinversion_singlearc() failed !" )



if( ! testinversionrotation(count=1000) )  stop( "testinversionrotation() failed !" )





cat( "Passed all rotation inversion tests !\n", file=stderr() )
