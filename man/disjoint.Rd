\name{arc disjointness}

\alias{disjointarcs}


\title{
arc disjointness
}


\description{
Test whether a set of arcs are pairwise strictly disjoint.
}

\usage{
disjointarcs( arcmat )
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
}
}

\value{
\code{disjointarcs()} returns a logical - whether
the given arcs are pairwise strictly disjoint.
If 2 arcs overlap, or are abutting, the function returns \code{FALSE}.
}

\seealso{
\code{\link{plotarcs}()},
\code{\link{complementaryarcs}()},
\code{\link{arcsintersection}()},
\code{\link{arcsunion}()},
\code{\link{arcssymmdiff}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

