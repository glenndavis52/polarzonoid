


.onLoad <- function( libname, pkgname )
    {
    if( requireNamespace( 'logger', quietly=FALSE ) )
        {
        #   log_formatter( formatter_mine )
        log_formatter( logger::formatter_sprintf, namespace="polarzonoid" )     # force sprintf(), even if glue is installed
        log_layout( layout_mine, namespace="polarzonoid" )                      # put fn() between timestamp and the msg    
        log_threshold( WARN, namespace="polarzonoid" )                          # default is INFO
        }
    }



.onAttach <- function( libname, pkgname )
    {

    }
