\name{arcs and boundary}

\alias{boundaryfromarcs}
\alias{arcsfromboundary}

\title{
The Homeomorphism between the Space of \eqn{n} or fewer arcs,
and the boundary of the Polar Zonoid
}

\description{
This section calculates the homeomorphism from
the space of \eqn{n} or fewer arcs on the circle,
which is denoted by \eqn{A_n},
to the boundary of the polar zonoid, and its inverse.
}

\usage{
boundaryfromarcs( arcmat, n=NULL, gapmin=0 )
arcsfromboundary( p, tol=5.e-9 )
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
For these 2 improper arcs, there must be only 1 row,
and the center is ignored.
}

\item{n}{
the given set of arcs is taken to be in \eqn{A_n}.
If \code{n} is \code{NULL}, then \code{n} is set to
the number of rows in \code{arcmat},
unless \code{arcmat} defines an improper arc,
when \code{n} is set to 0.
If \code{n} is not \code{NULL}, then we must have
\code{n} \eqn{\ge} the default value.

The returned vector is in \eqn{\mathbb{R}^{2n+1}}.
}

\item{gapmin}{
the minimum gap between arcs in \code{arcmat} that is valid.
The default \code{gapmin=0} allows abutting arcs,
but not overlapping arcs.
If one is sure that the arcs are not overlapping,
then set \code{gapmin=-Inf} and this validation check
is skipped, which saves a tiny bit of time.
}

\item{p}{
a vector in \eqn{\mathbb{R}^{2N+1}} or \eqn{\mathbb{R}^{2N}}.
In the latter case (an even-dimensional space),
a \eqn{\pi} is appended to it to make p odd-dimensional.
The vector must be on the boundary of the polar zonoid in
\eqn{\mathbb{R}^{2N+1}}, which usually means that it is
computed by either
\code{boundaryfromarcs()} or \code{boundaryfromsphere()}.
In this version of the package,
valid values for \eqn{N} are 0, 1, 2, and 3.
}

\item{tol}{
if the approximate distance from \code{p} to a
sub-stratum of the boundary is less than \code{tol},
then take \code{p} to actually \emph{be} in this substratum.
This step is necessary because the boundary is not smooth
at the sub-strata.
}

}



\value{
\code{boundaryfromarcs()} maps from \eqn{A_n}
to \eqn{\mathbb{R}^{2n+1}}.
It returns the computed point on the boundary of the zonoid.
Names are assigned indicating the corresponding term
in the trigonometric polynomial.
\cr
In case of error, the function returns \code{NULL}.

If \eqn{m} is the length of \code{p},
then \code{arcsfromboundary()}
maps from \eqn{\mathbb{R}^{m}} to \eqn{A_N}.
It returns an Nx2 matrix defining N arcs as above.
Because \code{p} might be in a substratum of the boundary,
N might be less than expected, which is \eqn{(m-1)/2}.
\cr
In case of error, the function returns \code{NULL}.
}


\details{
In \code{boundaryfromarcs()}, the calculation of the returned point
is a straightforward integration of
trigonometric functions over the arcs.
The last component of the vector is simply the sum of
the arc lengths.
\cr
Let \eqn{h} denote \code{boundaryfromarcs()}.
If \eqn{s} denotes the full circle,
then \eqn{h(s)} = (0,...,0,\eqn{2\pi}).
If \eqn{\phi} denotes the empty arc, then \eqn{h(\phi)} = (0,...,0,0).


\code{arcsfromboundary()} first calculates the outward pointing
normal at \code{p}.
In this version of the package, the implicitization of the boundary,
and thus the normal, is only available when \eqn{N} is 0, 1, 2, or 3.
This normal determines a trigonometric polynomial
whose roots are calculated with \code{trigpolyroot()}.
These roots are the endpoints of the arcs.

These two functions are inverses of each other.
}


\examples{
# make two disjoint arcs
arcmat = matrix( c(pi/4,pi/2,  pi,pi/4), 2, 2, byrow=TRUE ) ; arcmat

##           [,1]      [,2]
## [1,] 0.7853982 1.5707963
## [2,] 3.1415927 0.7853982

plotarcs( arcmat )

#  map to boundary of the zonoid
b = boundaryfromarcs( arcmat ) ;  b

##        x1        y1        x2        y2         L 
## 0.2346331 1.0000000 0.7071068 1.0000000 2.3561945 

#  map b back to arcs, and compare with original
arcsdistance( arcmat, arcsfromboundary(b) )

## [1] 2.220446e-16

# so the round trip returns to original pair of arcs, up to numerical precision
}



\seealso{
\code{\link{spherefromarcs}()},
\code{\link{arcsfromsphere}()},
\code{\link{boundaryfromsphere}()},
\code{\link{spherefromboundary}()},
\code{\link{complementaryarcs}()},
\code{\link{trigpolyroot}()}
}





% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

