\name{arc operations}


\alias{complementaryarcs}
\alias{arcsintersection}
\alias{arcsunion}
\alias{arcssymmdiff}
\alias{arcsdistance}


\title{
Boolean arc operations, and also the distance between two collections of arcs.
}


\description{
Perform Boolean operations on collections of arcs in the circle.
Arcs are always considered to be \emph{closed},
i.e. to contain their endpoints.
}

\usage{
complementaryarcs( arcmat )
arcsintersection( arcmat1, arcmat2 )
arcsunion( arcmat1, arcmat2 )
arcssymmdiff( arcmat1, arcmat2 )

arcsdistance( arcmat1, arcmat2 )
}
\arguments{
\item{arcmat}{
an Nx2 matrix with an arc definition in each row;
so the total number of arcs is N.
The 1st number in the row is the \emph{center} of the arc,
and the 2nd number is the \emph{length} of the arc;
both in radians.\cr
The length can be 0 or \eqn{2\pi},
defining the empty arc and the full circle respectively.
For these improper arcs, there must be only 1 row,
and the center is ignored.

The given arcs must be strictly disjoint.
}

\item{arcmat1, arcmat2}{
2 matrices that define 2 collections of arcs, as in the previous argument.
The number of arcs in each collection are not required to be equal.
}
}



\value{
\code{complementaryarcs()} returns
a matrix of the same size, which represents
the closure of the complement of the union of the given arcs,
as a subset of the circle.
The given arcs must be strictly disjoint;
if not the the function logs and error and returns \code{NULL}.

\code{arcsintersection()}, \code{arcsunion()}, and \code{arcssymmdiff()}  return
the intersection, union, and symmetric difference of the 2 given arc collections,
respectively.

\code{arcsdistance()} returns the distance between two collections of arcs.
This distance is the sum of the arc lengths of the symmetric difference.
}



\seealso{
\code{\link{plotarcs}()},
\code{\link{disjointarcs}()}
}





% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

