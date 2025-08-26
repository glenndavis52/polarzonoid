

#   coeff   complex coefficients of a trig polynomial, returned from trigpolyfromroots()
#           length is 2*n+1
#           exponents from +n to -n , OR   -n to +n   (does not matter)
#
#   returns a non-zero point on the line of symmetry

linegenerator <- function( coeff )
    {
    if( length(coeff) %% 2 != 1L )
        {
        log_level( WARN, "length = %d is invalid.  It must be even.", length(coeff) )
        return(NA_complex_)
        }

    n   = floor( length(coeff) / 2 )

if( TRUE )
    {
    pospart = coeff[ n+1 + 1:n ]            # increasing,  1 to n
    negpart = coeff[ n+1 + ((-1):(-n)) ]    # decreasing, -1 to -n

    #   compute the unrotated complex a's, they all lie on the line of symmetry
    a       = negpart + pospart
    amod    = abs(a)
    ia      = which.max( amod )

    #   compute the unrotated complex b's, they all lie on the line that is orthogonal to the line of symmetry
    b       = negpart - pospart
    bmod    = abs(b)
    ib      = which.max( bmod )

    if( bmod[ib] <= amod[ia] )
        out     = a[ia]
    else
        out     = b[ib] * 1i

    #cat( "a =", a, '\n' )
    #cat( "b =", b, '\n' )
    #cat( "a0 =", coeff[n+1], '\n' )
    #cat( "out =", out, '\n' )
    }
else
    {
    #   the middle coefficient, exponent 0, is at n+1
    out = coeff[n+1]

    if( abs(out) == 0 )
        {
        #   TODO:  handle this case too; it should be possible
        log_level( WARN, "out=0 is invalid.\n" )
        return(NA_complex_)
        }
    }

    return( out )
    }


#   coeff   complex coefficients from -n to +n
#           the line of symmetry must be the real axis
#           this symmetry property is enforced by trigpolyfromroots()
#
#   returns list with 3 items
#       a0  the constant coeff
#       a   coeffs of cos(i*t), from 1 to n
#       b   coeffs of sin(i*t), from 1 to n

realfromcomplexcoeff <- function( coeff )
    {
    n   = floor( length(coeff) / 2 )

    out = list( a0=Re(coeff[n+1]) )


if( TRUE )
    {
    pospart = coeff[ n+1 + 1:n ]            # increasing,  1 to n
    negpart = coeff[ n+1 + ((-1):(-n)) ]    # decreasing, -1 to -n

    out$a   = Re( negpart + pospart )
    out$b   = Im( negpart - pospart )
    }
else
    {
    out$a   = numeric(n)
    out$b   = numeric(n)

    for( k in 1:n )
        {
        out$a[k] = Re( coeff[n+1-k]  +  coeff[n+1+k] )

        out$b[k] = Im( coeff[n+1-k]  -  coeff[n+1+k] )
        }
    }

    return( out )
    }


#   ablist  a list with 3 items
#           a0  the constant coeff
#           a   coeffs of cos(i*t), i = 1 to n
#           b   coeffs of sin(i*t), i = 1 to n
#   all these coeffs are real
#
#   returns vector of 2*n+1 complex numbers, from c_{-n} to c_n

complexfromrealcoeff <- function( ablist )
    {
    negpart = complex( real=ablist$a, imaginary=ablist$b ) / 2
    pospart = complex( real=ablist$a, imaginary=-ablist$b ) / 2

    out = c( rev(negpart), ablist$a0, pospart )

    return( out )
    }

#   vec     a vector of odd length
#           the coefficients are in order a1,b1,  a2,b2, .... ,  an,bn, a0
    
trigpolyfromvector <- function( vec )
    {
    m       = length(vec)   # m is odd
    
    n       = (m-1L)/2L
    
    if( 0 < n )
        {
        #   the usual case
        even    = 2 * (1:n)
        odd     = even - 1L

        ablist  = list( a0=vec[m], a=vec[odd], b=vec[even] )
        }
    else
        {
        #   trivial case, a polynomial of degree 0
        ablist  = list( a0=vec[m], a=numeric(0), b=numeric(0) )
        }
        
    return( ablist )
    }
    
    
    

#   coeff   complex coefficients, starting with leading term and ending in constant term
#   zvec    vector of z's at which to evaluate, usually on the unit circle

evalpoly    <- function( coeff, zvec )
    {
    out = 0     # ; rep(0,length(zvec))

    for( k in 1:length(coeff) )
        {
        out  = coeff[k] + out*zvec
        }

    return( out )
    }


#   ablist  a0, a, and b.  all real
#   theta   vector of angles, usually in [0,2*pi]

evaltrigpoly    <- function( ablist, theta )
    {   
    n   = length(ablist$a)
    
    if( n == 0 )    return( ablist$a0 )
    
    out =  ablist$a0 + sum( ablist$a * cosmy( (1:n)*theta ) )  +  sum( ablist$b * sinmy( (1:n)*theta ) )

    return( out )
    }






#   ablist  a list with 3 items
#           a0  the constant coeff
#           a   coeffs of cos(i*t), i = 1 to n
#           b   coeffs of sin(i*t), i = 1 to n
#   all these coeffs must be real
#
#   tol     tolerance for the imaginary part, to extract pure reals
#
#   returns vector of real roots in [0,2*pi) in increasing order
#   the length is always even, and possibly 0
#   roots with multiplicities are repeated
#
#   if ablist is invalid, or if all coefficients are 0, it returns NULL

trigpolyroot <- function( ablist, tol=.Machine$double.eps^0.5 )
    {
    #   must be a list of length 3
    ok  = is.list(ablist) && length(ablist)==3 && all( sapply(ablist,is.numeric) )   # &&  all( 0<lengths(ablist) )
    if( ! ok )
        {
        log_level( ERROR, "Argument ablist is invalid. It must be a list of 3 numeric vectors." )
        return(NULL)
        }

    if( ! all( names(ablist)==c('a0','a','b') ) )
        {
        log_level( ERROR, "Argument ablist is invalid. The names are not 'a0','a','b'." )
        return(NULL)
        }

    if( length(ablist$a0) != 1 )
        {
        log_level( ERROR, "Argument ablist is invalid. a0 has length %d, but it must be 1.", length(ablist$a0) )
        return(NULL)
        }

    n   = length( ablist$a )

    if( length(ablist$b) != n )
        {
        log_level( ERROR, "The lengths of a and b are not equal." )
        return(NULL)
        }

    coeff   = complexfromrealcoeff( ablist )

    if( all( coeff==0 ) )
        {
        log_level( WARN, "All coefficients of the polynomial are 0." )
        return( NULL )
        }
        
    if( n == 0 )    return( complex(0) )

    #   get the complex roots
    root    = base::polyroot( coeff )  #;  print( root )

    #   convert to real
    out = base::log(root) / (1i)    #;    print( out )

    absim       = abs( Im(out) )

    maskreal    = absim <= tol

    if( ! any( maskreal ) )
        # no real roots
        return( numeric(0) )

    out = out[maskreal]

    if( length(out) %% 2 == 1L )
        {
        #   an odd number of roots, possibly only 1 root
        #   delete the one with the largest imaginary part
        absim   = absim[maskreal]   # make absim[] match out
        imax    = which.max( absim )

        out = out[-imax]
        }

    out = sort( Re(out) %% (2*pi) )

    return( out )
    }


#   ablist  list of real coeffs, as returned from realfromcomplexcoeff()

plottrigpoly <- function( ablist, rootvec=NULL )
    {
    theta   = seq( 0, 2*pi, length.out=361 )

    n   = length( ablist$a )

    y   = rep( ablist$a0, length(theta) )

    for( k in seq_len(n) )
        {
        y   = y + ablist$a[k] * cos( k * theta )    +    ablist$b[k] * sin( k * theta )
        }

    plot( theta, y, type='l' )

    abline( h=0 )

    if( ! is.null(rootvec) )
        points( rootvec, numeric(length(rootvec)) )
    }



#   z1, r1      center and radius of circle #1 in complex plane
#   z2, r2      center and radius of circle #1 in complex plane
#
#   0 < d(z1,z2) <= r1 + r2
#
#   returns vector 2 complex numbers, or NA_complex_

circleintersection  <- function( z1, r1, z2, r2 )
    {
    ok  = (0 <= r1) && (0 <= r2)
    if( ! ok )
        {
        return( NA_complex_ )
        }

    delta   = z2 - z1

    d   = abs( delta )

    if( r1 == 0 )
        {
        if( abs(d - r2) <= 5.e-16 )
            #   z1 is within tolerance of the other circle
            return( c(z1,z1) )
        else
            return( NA_complex_ )
        }

    if( r2 == 0 )
        {
        if( abs(d - r1) <= 5.e-16 )
            #   z2 is within tolerance of the other circle
            return( c(z2,z2) )
        else
            return( NA_complex_ )
        }


    ok  = (0 < d) && (d <= r1+r2)
    if( ! ok )
        {
        return( NA_complex_ )
        }

    lambda  = (d^2 - r2^2 + r1^2)/ (2*d^2)

    #   the "radical center"
    zmid    = (1-lambda)*z1 + lambda*z2

    mu  = sqrt( max( r1^2 - (lambda*d)^2, 0 ) ) / d

    #   rotate delta pi/2, and scale by mu
    offset  = mu * (delta * 1i)

    out = c( zmid - offset, zmid + offset )

    return( out )
    }