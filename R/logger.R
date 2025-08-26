
#   put fn() between timestamp and the msg    
layout_mine <- structure(
    function(level, msg, namespace="polarzonoid",
                                    .logcall = sys.call(), .topcall = sys.call(-1), .topenv = parent.frame())
        {
        # cat( "obj_addr()=", obj_addr( .topcall[[1L]] ), '\n' )
        # cat( "deparse1 =", deparse1( .topcall[[1L]] ), '\n' )
        
        fn  = deparse1( .topcall[[1L]] )
        
        paste0( attr(level, 'level'), ' [', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '] ', namespace, "::", fn, '(). ', msg )
        },
    generator = quote(layout_mine())
)
