
#   theta   an increasing vector of angles in [0,2*pi)
#           it must have even length
#           for a sequence of length 2*n, there are n arcs,
#           Exception: when n=0, there is one empty or full arc
#   comp    take complementary arcs.
#           When TRUE, 0 is in an arc.  When FALSE, 0 is not in an arc.
#
#   returns an arcmat[,]

arcsfromendpoints <- function( theta, comp=FALSE )
    {
    if( (length(theta) %% 2) == 1L )
        {
        log_level( ERROR, "length(theta) = %d is invalid.  It must be even.", length(theta) )
        return(NULL)
        }

    n   = length(theta) / 2L    # the number of arcs

    if( n == 0 )
        {
        #   in this special case, make a single improper arc, of length 0 or 2*pi
        out = matrix( 0, nrow=1, ncol=2 )
        colnames( out ) = p.arcmatcolnames

        #out[1,1]    = 0     # center
        if( comp )  out[1,2]    = 2*pi

        return( out )
        }

    out = matrix( 0, nrow=n, ncol=2 )
    colnames( out ) = p.arcmatcolnames

    even    = seq( 2L, length(theta), by=2L )
    odd     = even - 1L

    if( comp )
        #   shift last to the beginning
        theta   = c( theta[2*n] - 2*pi, theta[1:(2*n-1)] )

    out[ ,1]    = (theta[odd] + theta[even]) / 2
    out[ ,2]    = theta[even] - theta[odd]

    return( out )
    }



arcsfromendpointsandpoly <- function( theta, ablist )
    {
    out = arcsfromendpoints( theta, TRUE )

    #   the trig polynomial must be positive on all arcs
    #   complement arcs if necessary

    #   find the longest arc
    imax    = which.max( out[ ,2] )
    center  = out[imax,1]  # center of longest arc

    if( evaltrigpoly(ablist,center) < 0 )
        #  cat( "complementing arcs\n" )
        out = arcsfromendpoints( theta, FALSE )

    return( out )
    }




#   arcmat  Nx2 matrix with center and length in the rows
#
#   returns Nx2 matrix with thetamin and thetamax in the rows
#   the arcs are sorted in counterclockwise order

endpointsfromarcs <- function( arcmat )
    {
    #   ensure these are in proper order
    center  = arcmat[ ,1]  %%  (2*pi)
    perm    = order( center )
    arcmat  = arcmat[ perm,  , drop=FALSE]

    n   = nrow(arcmat)

    out = matrix( 0, n, 2 )

    colnames(out)   = c('thetamin','thetamax')

    out[ ,1]    = arcmat[ ,1] - arcmat[ ,2]/2
    out[ ,2]    = arcmat[ ,1] + arcmat[ ,2]/2

    return( out )
    }



#   arcmat  mx2 matrix with center,length in the rows
#
#   returns mx2 matrix with the complementary arcs in the same format
#
#   this function does NOT use arcsfromendpoints() or endpointsfromarcs(), though it could :-)
#
#   for the improper arcs (empty arc or full circle) the center does not matter
#   and the center is set to 0
#
#   if arcs are not strictly disjoint, returns NULL

complementaryarcs <- function( arcmat )
    {
    arcmat   = prepareNxM( arcmat, 2 )

    n   = nrow( arcmat )

    if( n == 1 )
        {
        #   only 1 arc is easy
        #   this section also handles the improper arcs
        out     = matrix( 0, nrow=1, ncol=2 )

        out[1,1]    = (arcmat[1,1] + pi)  %%  (2*pi)    # antipodal
        out[1,2]    = 2*pi - arcmat[1,2]

        if( out[1,2] %in% c(0,2*pi) )   out[1,1] = 0    # improper arc

        colnames(out)   = p.arcmatcolnames

        return(out)
        }


    theta   = endpointsfromarcs( arcmat )

    #   compute thetamin and thetamax for the complementary arcs
    #   we use the known fact that for elements in the nx2 matrix theta,
    #   the row index changes fastest
    thetamin    = theta[ (n+1L):(2*n) ]
    thetamax    = theta[ c(2:n,1L) ]

    len     = (thetamax - thetamin) %% (2*pi)
    center  = (thetamin + len/2) %% (2*pi)

    out     = cbind( center, len )

    colnames(out)   = p.arcmatcolnames

    return( out )
    }





#   arcmat  Nx2 matrix defining arcs in center,length form

#   tests whether the arcs are *strictly* disjoint. abutting arcs are not strictly disjoint
#   returns TRUE or FALSE

disjointarcs <- function( arcmat )
    {
    if( nrow(arcmat) == 1 ) return(TRUE)

    gapmin  = gapminimum( arcmat )

    return( 0 < gapmin )
    }





#   arcmat  Nx2 matrix defining arcs in center,length form

#   returns smallest gap between consecutive arcs
#   if arcs abutt, this is 0
#   if arcs overlap, this is negative
#   if only 1 arc, returns the length of the complementary arc

gapminimum <- function( arcmat )
    {
    m   = nrow( arcmat )

    if( m == 1 )
        {
        #   only 1 arc is easy
        return( 2*pi - arcmat[1,2] )
        }

    center  = arcmat[ ,1]
    len     = arcmat[ ,2]

    center  = center %% (2*pi)
    perm    = order( center )

    center  = center[perm]
    len     = len[perm]

    #cat( "center = ", center, '\n' )
    #cat( "len = ", len, '\n' )

    out = 2*pi

    for( k in 1:m )
        {
        if( k == m )
            {
            end0    = center[k]-2*pi + len[k]/2
            end1    = center[1] - len[1]/2
            }
        else
            {
            end0    = center[k] + len[k]/2
            end1    = center[k+1] - len[k+1]/2
            }

        #cat( "k=", k, "  end0=", end0,  "  end1=", end1, '\n' )

        #   the gap is end1 - end0
        out = min( end1 - end0, out )
        }

    return( out )
    }



#   n   the # of arcs.  0 means return one of the 2 improper arcs.
#   L   total length of the arcs, we must have L < 2*pi
#       if NA_real_ then total length is random
#
#   returns arcmat[] with arcs strictly disjoint

randomarcs  <- function( n, L=NA_real_ )
    {
    if( n == 0 )
        {
        #   an improper arc, ignore L in this case
        L   = ifelse( runif(1) < 0.5, 0, 2*pi )

        arcmat  = cbind( center=0, length=L )

        return( arcmat )
        }

    if( is.na(L) )
        {
        #   no total length given; this is easier

        theta   = sort( runif( 2L*n, min=0, max=2*pi ) )

        comp    = ( runif(1) < 0.5 )

        arcmat     = arcsfromendpoints( theta, comp )

        return( arcmat )
        }


    #   L is given, so this is harder, with lots of scaling
    ok  = (0 < L) && (L < 2*pi)
    if( ! ok )
        {
        log_level( ERROR, "L = %g is invalid.", L )
        return(NULL)
        }

    theta   = runif( 1, min=0, max=2*pi )

    if( n == 1 )
        {
        arcmat  = cbind( center=theta, length=L )
        return( arcmat )
        }

    #   create unscaled lengths and gaps
    len = runif( n )
    gap = runif( n )

    G   = 2*pi - L      # total gap length

    #   now scale lengths and gaps to have the desired sums
    len = (L * len) / sum(len)
    gap = (G * gap) / sum(gap)      #; cat( "gap=", gap, '\n' )

    delta   = sum(len) - L
    if( delta != 0 )
        {
        #cat( "delta=", delta, '\n' )
        imin    = which.min(len)
        len[imin]   = len[imin] - delta

        #delta2  = sum(len) - L
        #if( delta2 != 0 )
        #    cat( "delta2=", delta2, '\n' )
        }

    center  = len/2 + c( 0, cumsum(len + gap)[1:(n-1)] )

    #   now rotate a random amount
    center  = (center + theta) %% (2*pi)

    #   now sort ascending
    perm    = order(center)

    arcmat  = cbind( center=center[perm], length=len[perm] )

    return( arcmat )
    }



#   zmat    Nx2 complex matrix
#           each row has starting endpoint z1 and ending z2, on unit circle

centerlengthmat <- function( zmat )
    {
    n   = length(zmat)/2

    if( ! is.matrix(zmat) )     zmat = matrix( zmat, nrow=n, ncol=2, byrow=TRUE )

    #n   = nrow(zmat)

    out = matrix( 0, nrow=n, ncol=2 )
    colnames( out ) = p.arcmatcolnames

    for( k in 1:n )
        {
        if( any( is.na(zmat[k, ]) ) )   next

        theta   = Arg( zmat[k, ] )

        L = (theta[2] - theta[1]) %% (2*pi)

        theta[2]    = theta[1] + L

        out[k,1]    = sum(theta)/2
        out[k,2]    = L
        }

    return( out )
    }



#   dat an Nx2 matrix, either
#       a real matrix with center angle, and arc length
#           or
#       a complex matrix with starting endpoint z1 and ending z2, on unit circle

plotarcs <- function( arcmat, labels=FALSE, main=NULL, margintext=NA, rad=1, lwd=3, pch=20, cex=1.5, add=FALSE, ...  )
    {
    if( is.complex(arcmat) )
        arcmat  = centerlengthmat(arcmat)

    arcmat   = prepareNxM( arcmat, 2 )
    if( is.null(arcmat) )    return( invisible(FALSE) )
    
    vararg  = list(...)
    
   

    if( ! add )
        {
        #   initialize the plot
        plot.default( c(-1,1), c(-1,1), type='n', xlab='', ylab='', asp=1, las=1, tcl=0, mgp=c(3,0.25,0) )
        title( xlab='x', line=1.5 )
        title( ylab='y', line=1.5 )

        grid( lty=1 )
        abline( h=0, v=0 )

        points( 0, 0, pch=20 )

        theta   = seq( 0, 2*pi, len=361 )

        lines( cos(theta), sin(theta) )
        }


    dolines = is.numeric(lwd)  &&  length(lwd)==1  &&  0<=lwd

    for( k in 1:nrow(arcmat) )
        {
        L   = arcmat[k,2]

        if( L == 0 )    next    # the empty arc

        count   = ceiling( 361 * L/(2*pi) )

        center = arcmat[k,1]

        theta   = seq( center-L/2, center+L/2, len=count )

        if( dolines )
            lines( rad*cos(theta), rad*sin(theta), lwd=lwd, ... )

        theta   = theta[ c(1,count) ]

        if( L < 2*pi )
            #   not the full circle
            points( rad*cos(theta), rad*sin(theta), pch=pch, cex=cex, ... )

        if( labels )
            {
            label   = sprintf( "L%d=%g", k, L )
            text( rad*cos(center), rad*sin(center), label, cex=0.75, adj=c(1,0) )
            }
        }

    L   = sum( arcmat[ ,2] )

    if( ! add )
        {
        #   chordsum    = sum( chordsfromarcs(arcmat) )

        if( 0 < L ) n   = nrow(arcmat)
        else        n   = 0L

        if( is.null(main) )
            main    = sprintf( "|arcs|=%d  total L=%g  disjoint=%s",
                                n, L, disjointarcs(arcmat) )
        if( ! is.na(main) )
            title( main=main )

        if( ! is.na(margintext) )
            title( xlab=margintext, line=3 )
        }

    return( invisible(TRUE) )
    }





arcsintersection <- function( arcmat1, arcmat2 )
    {
    arcscombo( arcmat1, arcmat2, sumind=2L )
    }


arcsunion <- function( arcmat1, arcmat2 )
    {
    arcscombo( arcmat1, arcmat2, sumind=1:2 )
    }


arcssymmdiff <- function( arcmat1, arcmat2 )
    {
    arcscombo( arcmat1, arcmat2, sumind=1L )
    }


#   returns the total length of the arcs in the symmetric difference
#
#   this is the L1 metric

arcsdistance <- function( arcmat1, arcmat2 )
    {
    arcmat  = arcssymmdiff( arcmat1, arcmat2 )
    if( is.null(arcmat) )   return(NULL)

    return( sum( arcmat[ ,2] ) )
    }




#   arcmat1 1st set of arcs
#   arcmat2 2nd set of arcs
#   sumind  desired sum of the indicator functions
#               2L      =>  intersection
#               1:2     =>  union
#               1L      =>  symmetric-difference
#
#   computes the jumps in the sum of the indicator functions
#   returns the set of arcs where this sum is in sumind

arcscombo <- function( arcmat1, arcmat2, sumind=2L )
    {
    theta1  = endpointsfromarcs( prepareNxM(arcmat1,2) )
    theta2  = endpointsfromarcs( prepareNxM(arcmat2,2) )

    if( is.null(theta1)  ||  is.null(theta2) )  return(NULL)

    
    df  = stepdata( theta1, theta2 )
    if( is.null(df) )   return(NULL)
    
    # print( df )
    
    
    #   start the scan at the midpoint of the arc from the largest theta to 2*pi
    
    thetastart  = (df$theta[ nrow(df) ] + 2*pi) / 2  # this is guaranteed to be < 2*pi

    #   in the next line TRUE converts to 1L, and FALSE converts to 0L
    sumstart    = pointsinarcs(thetastart,arcmat1) + pointsinarcs(thetastart,arcmat2)

    log_level( TRACE, "Starting the scan at theta=%g, where sum function = %d.", thetastart, sumstart )
    
    
    df  = collapse_stepdata( df )
    if( is.null(df) )   return(NULL)

    #  print( df )
    
    m   = nrow(df)

    log_level( TRACE, "The number of jumps = %d.", m )
    
    #   the number of jumps in df must be even    
    if( FALSE  &&  m %% 2 == 1L )
        {
        log_level( FATAL, "Internal Error. The number of jumps = %d, which is invalid because it is odd.", m )
        return(NULL)
        }

    if( m == 0 )
        {
        #   no jumps, it means that the sum function is constant
        if( arcmat1[1,2]==0  &&  arcmat2[1,2]==0 )
            sumexp  = 0L
        else if( arcmat1[1,2]==2*pi  &&  arcmat2[1,2]==2*pi )
            sumexp  = 2L
        else
            sumexp  = 1L
            
        if( sumstart != sumexp )
            {
            log_level( FATAL, "Internal Error. The sum function is constant, and expected sum is %d, but it is %d.",
                                sumexp, sumstart )
            return(NULL)
            }

        #   the combo arcmat is improper
        #   the length is either 0 (empty arc) or 2*pi (full circle)
        #   is thetastart in the set ?
        inset   = sumstart %in% sumind
        L   = ifelse( inset, 2*pi, 0 )

        arcmat  = matrix( c(0,L), nrow=1 )

        return( arcmat )
        }

    #   there are jumps in the sum of the indicator functions,
    #   so now we have to do real work !
    #   the columns of df are theta and transition


    thesum  = sumstart + cumsum( df$transition ) #;  cat( "thesum=", thesum, '\n' )

    #   check that all values are in 0:2
    ran = range(thesum)

    if( ran[1] < 0  ||  2 < ran[2] )
        {
        log_level( FATAL, "Internal Error.  range(thesum) = %d,%d,  which is outside 0:2.",
                        ran[1], ran[2] )
        return(NULL)
        }

    inset   = thesum %in% sumind #;  cat( "inset=", inset, '\n' )

    if( all(inset == inset[1]) )
        {
        #   all T   or   all F, means an improper arc
        L   = ifelse( inset[1], 2*pi, 0 )

        arcmat  = matrix( c(0,L), nrow=1 )

        return( arcmat )
        }


    #   now find the runs of TRUE in inset
    ss  = findRunsTRUE( inset ) #; print( ss )

    #   for each run of TRUEs, there is an output arc
    #arcs    = nrow(ss)

    istart  = ss[ ,1]
    istop   = (ss[ ,2]  %%  length(inset) ) + 1L

    thetamin    = df$theta[istart]
    thetamax    = df$theta[istop]

    len     = (thetamax - thetamin)  %%  (2*pi)
    center  = (thetamin + len/2)  %%  (2*pi)

    out     = cbind( center, len )

    colnames(out)   = p.arcmatcolnames

    return(out)
    }



#   theta   vector of points on the circle
#   arcmat  center and length matrix
#
#   returns distance from the point to the union of the arcs
#   if the point is in an arc, the distance is 0
#
#   if arcmat is the full circle, returns 0
#   if arcmat is the empty arc, returns Inf

distancetoarcs <- function( theta, arcmat )
    {
    m   = length(theta)

    out = rep( Inf, m )

    n   = nrow(arcmat)

    for( k in 1:m )
        {
        for( i in 1:n )
            {
            len = arcmat[i,2]
            if( len == 0 )
                {
                #   the empty arc, leave it Inf
                break
                }

            d   = distancecircle( theta[k], arcmat[i,1] ) - len/2

            if( d <= 0 )
                {
                #   theta is inside an arc, so distance is 0 and we're done
                out[k]  = 0
                break
                }

            out[k]  = min( d, out[k] )
            }
        }

    return( out )
    }


#   theta   a vector of angles, which are points on the unit circle
#   arcmat  closed arcs in center,length format
#
#   returns TRUE iff the point is in one of the arcs

pointsinarcs <- function( theta, arcmat )
    {
    return( distancetoarcs( theta, arcmat ) == 0 )
    }



#   theta1  a point on the unit circle
#   theta2  another point
#
#   returns the distance between theta1 and theta2, along the shorter arc
#   the distance is always in [0,pi]

distancecircle  <- function( theta1, theta2 )
    {
    out = abs( theta1 - theta2 ) %% (2*pi)

    if( pi < out )  out = 2*pi - out

    return( out )
    }




stepdata    <- function( theta1, theta2 )
    {
    min1    = data.frame( theta=theta1[ ,1], transition=+1L )  #, idx=1L )
    max1    = data.frame( theta=theta1[ ,2], transition=-1L )  #, idx=1L )

    min2    = data.frame( theta=theta2[ ,1], transition=+1L )  #, idx=2L )
    max2    = data.frame( theta=theta2[ ,2], transition=-1L )  #, idx=2L )

    out = rbind( min1, max1, min2, max2 )

    rownames(out)   = 1:nrow(out)

    out$theta   = out$theta %% (2*pi)   # wrap to standard range [0,2*pi)

    perm    = order( out$theta )

    out = out[ perm, , drop=FALSE ]

    return( out )
    }



#   df  a data.frame returned by stepdata()

collapse_stepdata <- function( df )
    {
    #   df$theta is already in increasing order
    datarle = rle( df$theta )

    if( all(datarle$lengths==1 ) )   return(df)  # theta is strictly increasing, so no change needed

    out     = df

    #   find the non-trivial groups
    idxdup  = which( 1 < datarle$lengths )

    csum    = cumsum( datarle$lengths )

    #   make logical mask for which rows to delete
    mask    = logical( nrow(df) )

    for( i in idxdup )
        {
        #   find the consecutive indexes for this group
        kstart  = csum[i] - datarle$lengths[i] + 1L
        kstop   = csum[i]

        idx = kstart:kstop

        transition  = sum( out$transition[idx] )

        if( transition == 0 )
            #   mark all the rows for deletion
            mask[ idx ] = TRUE
        else
            {
            #   set the last row to the sum
            out$transition[kstop]   = transition

            #   mark the other rows for deletion
            mask[ kstart:(kstop-1L) ] = TRUE
            }
        }

    #   now delete the marked rows
    out = out[ ! mask, ]

    #   check that the sum of the transitions is 0
    sumtrans    = sum( out$transition )
    if( sumtrans != 0 )
        {
        log_level( FATAL, "Internal Error. The sum of the transitions = %d, but it should be 0.", sumtrans )
        return(NULL)
        }


    return( out )
    }





#############################   deadwood below  ############################################


#   arcmat  mx2 matrix with center,length in the rows
#
#   returns mx2 matrix with the complementary arcs in the same format
#
#   this function does NOT use arcsfromendpoints() or endpointsfromarcs(), though it could :-)
#
#   for the improper arcs (empty arc or full circle) the center does not matter
#   and the center is set to 0
#
#   if arcs are not strictly disjoint, returns NULL

complementaryarcs_old <- function( arcmat )
    {
    m   = nrow( arcmat )

    center  = arcmat[ ,1]
    len     = arcmat[ ,2]

    out = matrix( 0, nrow=m, ncol=2 )
    colnames(out)   = colnames(arcmat)

    if( m == 1 )
        {
        #   only 1 arc is easy
        out[1,1]    = (arcmat[1,1] + pi)  %%  (2*pi)    # antipodal
        out[1,2]    = 2*pi - arcmat[1,2]

        if( out[1,2] %in% c(0,2*pi) )   out[1,1] = 0    # improper arc

        return(out)
        }

    center  = sort( center %% (2*pi) )

    for( k in 1:m )
        {
        end0    = center[k] + len[k]/2

        if( k == m )
            {
            end1    = center[1] - len[1]/2
            if( 0 < end1 )  end0 = end0 - 2*pi
            }
        else
            end1    = center[k+1] - len[k+1]/2

        out[k,1]    = (end0 + end1) / 2
        out[k,2]    = end1 - end0

        # cat( "k=", k, "  end0=", end0,  "  end1=", end1, '\n' )

        if( out[k,2] <= 0 )
            {
            log_level( ERROR, "Arcs are not disjoint. end1=%g <= %g=end0", end1, end0 )
            return( NULL )
            }
        }

    return( out )
    }




#   theta   a vector of angles, which are points on the unit circle
#   arcmat  closed arcs in center,length format
#
#   returns TRUE iff the point is in the interior of one of the arcs

pointinarcs_old <- function( theta, arcmat )
    {
    m   = length(theta)

    out = logical(m)    # initially all FALSE

    for( k in 1:m )
        {
        for( i in 1:nrow(arcmat) )
            {
            center  = arcmat[i,1]
            radius  = arcmat[i,2] / 2   # radius is length/2

            #   rotate theta[k] by -center
            #   to get thetarot in (-pi,pi]
            thetarot    = (theta[k] - center) %% (2*pi)
            if( pi <= thetarot ) thetarot = thetarot - (2*pi)

            if( -radius <= thetarot  &&  thetarot < radius )
                {
                #   theta[k] is inside the i'th half-open arc
                out[k]  = TRUE
                break
                }
            }
        }

    return( out )
    }

    
    
#   returns distance between arcs in the 'quotient' topology

distancearcsQ <- function( arcmat1, arcmat2 )
    {
    if( nrow(arcmat1) != nrow(arcmat2) )    return( Inf )   #  TODO: fix later

    #   'normalize' both sets of centers
    perm1   = order( arcmat1[ ,1] %% (2*pi) )
    perm2   = order( arcmat2[ ,1] %% (2*pi) )

    #   compute distance between centers, on the circle
    distcen = (arcmat1[perm1,1] - arcmat2[perm2,1]) %% (2*pi)
    mask    = pi <= distcen
    distcen[mask]  = 2*pi - distcen[mask]

    #   compute distance between lengths
    distlen = abs(arcmat1[perm1,2] - arcmat2[perm2,2])

    return( max( distcen, distlen, na.rm=TRUE ) )
    }
