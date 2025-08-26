\name{boundary and sphere}

\alias{boundaryfromsphere}
\alias{spherefromboundary}

\title{
The Homeomorphism between the Boundary of the Polar Zonoid
in \eqn{\mathbb{R}^{2n+1}}
and the sphere \eqn{\mathbb{S}^{2n}}
}

\description{
This section calculates the homeomorphism from
the space of \eqn{n} or fewer arcs on the circle
to the boundary of the polar zonoid, and its inverse.
In this version of the package, \eqn{n} must be
0, 1, 2, or 3.
}

\usage{
spherefromboundary( p )
boundaryfromsphere( x )
}

\arguments{
\item{p}{
an vector of length \eqn{2n{+}1} on the boundary of
the polar zonoid in that dimension.
\code{p} can also be even-dimensional, and in that case \eqn{\pi}
is appended to make \code{p} odd-dimensional.
}


\item{x}{
a non-zero vector of length \eqn{2n{+}1},
which is then unitized to put it on the sphere.
It can also be even-dimensional, and in that case 0
is appended to make \code{x} odd-dimensional.
}

}



\value{
\code{spherefromboundary()} returns
a unit vector in \eqn{\mathbb{S}^{2n}}.
In case of error, the function returns \code{NULL}.

\code{boundaryfromsphere()} returns the computed point on the boundary of the zonoid.
Names are assigned indicating the corresponding term
in the trigonometric polynomial.
In case of error, the function returns \code{NULL}.
}


\details{
\code{spherefromboundary()} is simply
central projection of the given point onto the unit sphere
centered at (0,0,...,0,0,\eqn{\pi}),
which is the center of the zonoid.
In this direction there is no restriction on \eqn{n}.


\code{boundaryfromsphere()} is much harder.
One must find the intersection of two objects:
1) the ray based at the center of the zonoid in the direction \code{x},
and 2) the boundary of the zonoid.
To do this, an implicit formula for the boundary has been
programmed, but only when 
\eqn{n} is 0, 1, 2, or 3.

These two functions are inverses of each other.
}






\seealso{
\code{\link{spherefromarcs}()},
\code{\link{arcsfromsphere}()},
\code{\link{boundaryfromarcs}()},
\code{\link{arcsfromboundary}()}
}





% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

