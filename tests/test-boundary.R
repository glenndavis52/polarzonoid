
library( polarzonoid )

options( width=160 )   


#   n       # of arcs
#   count   # of reps
#   tol     # for the discrepancy

#   return TRUE or FALSE    
testboundaryvalue  <- function( n, count=1000, tol=5.e-13 )
    {
    set.seed(0)
    
    bvalue  = numeric(count)
    
    for( i in 1:count )
        {
        arcmat  = polarzonoid:::randomarcs( n )

        p   = boundaryfromarcs( arcmat, gapmin=-Inf )
        
        bvalue[i]   = polarzonoid:::boundaryfunction( p )
        }
        
    ran = range(bvalue)
    cat( "n =", n,  "   bvalue range = ", ran, '\n', file=stderr() )
    
    if( tol < max( abs(ran) ) )
        {
        mask    = tol < abs(bvalue)
        
        mess    = sprintf( "ERR.  testboundaryvalue().  n=%d. %d of %d boundary values are invalid.\n",
                        n, sum(mask), length(mask) )
        cat( mess, file=stderr() )       
        
        return(FALSE)        
        }
    
    return(TRUE)
    }
    
    
#   n   number of arcs, this does not work if n==0
#    
#   return TRUE or FALSE    

testboundarynormal  <- function( n, count=1000, tol=5.e-12 )
    {
    set.seed(0)
    
    delta   = numeric(count)
    
    for( i in 1:count )
        {
        arcmat  = polarzonoid:::randomarcs( n )

        p   = boundaryfromarcs( arcmat, gapmin=-Inf )
        
        #   find normal using implicitization
        normal_imp  = polarzonoid:::boundarynormal( p )
        
        #   find normal using SVD        
        theta   = polarzonoid:::endpointsfromarcs( arcmat )
        
        normal_svd  = polarzonoid:::supportingnormal( theta, p )
        
        delta[i]    = max( abs(normal_imp - normal_svd) )
        }
        
    maxdelta    = max(delta)
    
    cat( "n =", n,  "   boundary normal max(delta) = ", maxdelta, '\n', file=stderr() )
            

    if( tol < maxdelta )
        {
        mask    = tol < delta
        
        mess    = sprintf( "ERR.  testboundarynormal().  n=%d. %d of %d boundary normals are invalid.  max(delta)=%g\n",
                        n, sum(mask), length(mask), maxdelta )
        cat( mess, file=stderr() )       
        
        return(FALSE)
        }
    
    return(TRUE)
    }
    
    
if( ! testboundaryvalue(0,count=10) )  stop( "testboundaryvalue(0) failed !" )
    
if( ! testboundaryvalue(1) )  stop( "testboundaryvalue(1) failed !" )
    
if( ! testboundaryvalue(2) )  stop( "testboundaryvalue(2) failed !" )

if( ! testboundaryvalue(3,tol=5.e-12) )  stop( "testboundaryvalue(3) failed !" )



if( ! testboundarynormal(1) )  stop( "testboundarynormal(1) failed !" )
   
if( ! testboundarynormal(2) )  stop( "testboundarynormal(2) failed !" )
   
if( ! testboundarynormal(3,tol=5.e-10) )  stop( "testboundarynormal(3,tol=5.e-10) failed !" )
   
    
cat( "Passed all boundary tests !\n", file=stderr() )       
