
#   x   an arc length <= 2*pi
#
#   returns the length of the chord, which is = 2*sin(x/2)
#
crd <- function( x )
    {
    2*sinmy( x/2 )
    }

#   len     length of a chord <= 2
acrd <- function( len )
    {
    2 * base::asin( len/2 )
    }


cosmy <- function( x )
    {
    base::cospi( x / pi )
    }
    
sinmy <- function( x )
    {
    base::sinpi( x / pi )
    }    
    
tanmy <- function( x )
    {
    (1 - cosmy(2*x)) / sinmy(2*x)
    }    
    
cotmy <- function( x )
    {
    sinmy(2*x) / (1 - cosmy(2*x))
    }        
    
    
cosdeg <- function( x )
    {
    base::cospi( x/180 )
    }
    
sindeg <- function( x )
    {
    base::sinpi( x/180 )
    }
    
    