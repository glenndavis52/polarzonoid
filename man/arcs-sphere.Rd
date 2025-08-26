\name{arcs and sphere}

\alias{arcsfromsphere}
\alias{spherefromarcs}
\alias{spherefromarcs_plus}


\title{
The Homeomorphism between the Space of \eqn{n} or fewer arcs,
and the sphere \eqn{\mathbb{S}^{2n}}
}

\description{
This section calculates the natural homeomorphism from
the space of \eqn{n} or fewer arcs on the circle,
denoted by \eqn{A_n}, 
to the sphere \eqn{\mathbb{S}^{2n}}, and its inverse.
}

\usage{
spherefromarcs( arcmat, n=NULL, gapmin=0 )
spherefromarcs_plus( arcmat, n=NULL, gapmin=5e-10 )
arcsfromsphere( u )
}
\arguments{
\item{arcmat}{
an \code{n0 x 2} matrix with an arc definition in each row;
so the total number of arcs is \code{n0}.
The 1st number in the row is the \emph{center} of the arc,
and the 2nd number is the \emph{length} of the arc;
both in radians.\cr
The length can be 0 or \eqn{2\pi},
defining the empty arc and the full circle respectively.
For these improper arcs, there must be only 1 row,
and the center is ignored.
}

\item{n}{
the given set of arcs is taken to be in \eqn{A_n}.
If \code{n} is \code{NULL}, then \code{n} is set to
\code{n0 = nrow(arcmat)},
unless \code{arcmat} defines an improper arc,
when \code{n} is set to 0.
If \code{n} is not \code{NULL}, then we must have
\code{n} \eqn{\ge} \code{n0}.

The returned vector is in \eqn{\mathbb{R}^{2n+1}}.
}


\item{gapmin}{
the minimum gap between arcs in \code{arcmat} that is valid.
The value \code{gapmin=0} allows abutting arcs,
but not overlapping arcs.
If one is sure that the arcs are not overlapping,
then set \code{gapmin=-Inf} and this validation check
is skipped, which saves a tiny bit of time.

In \code{spherefromarcs_plus()} abutting arcs are invalid,
so \code{gapmin} must be positive.
}


\item{u}{
a unit vector in \eqn{\mathbb{R}^{2N+1}} or \eqn{\mathbb{R}^{2N}}.
In the latter case (an even-dimensional space),
a 0 is appended to make \code{u} odd-dimensional.
}
}



\value{
\code{spherefromarcs()} maps from \eqn{A_n} to \eqn{\mathbb{S}^{2n}};
it returns a unit vector in \eqn{\mathbb{R}^{2n+1}}.
It is simply the composition of \code{boundaryfromarcs()}
and \code{spherefromboundary()},
with a little error checking.
\cr
In case of error, it returns \code{NULL}.

\code{spherefromarcs_plus()} is the same as \code{spherefromarcs()},
except it returns additional data.
It returns a list with these items:
\item{u}{
the same unit vector returned by \code{spherefromarcs()}.
Its length is \code{2*n+1}.
}
\item{tangent}{
a \code{2*n+1 x 2*n0} matrix whose columns are tangent, at \code{u}, to the substratum \eqn{\mathbb{S}^{2 n_0} \hookrightarrow \mathbb{S}^{2n}}.
This is the jacobian of the map \eqn{A_{n_0} \to \mathbb{S}^{2 n_0}},
when the space of arcs is parameterized by the endpoints of the arcs.
}
\item{normal}{
a \code{2*n+1 x 2*(n-n0)} matrix whose columns are normal, at \code{u},
to the substratum \eqn{\mathbb{S}^{2 n_0} \hookrightarrow \mathbb{S}^{2n}}.
These columns are an orthonormal basis for the subspace.
If \code{n==n0} there are 0 columns.
}
The matrix \code{cbind(u,tangent,normal)} is square.
When \code{n>n0}, the last normal vector is flipped if necessary,
so that the determinant of the square matrix is positive.

If \eqn{m} is the length of \code{u},
then \code{arcsfromsphere()}
maps from \eqn{\mathbb{S}^{m-1}} to \eqn{A_N}.
It returns an Nx2 matrix defining N arcs as above.
Because the space of arcs is stratified,
N might be less than expected, which is \eqn{(m-1)/2}.
It is simply the composition of \code{boundaryfromsphere()}
and \code{arcsfromboundary()},
with a little error checking.
In this version of the package,
valid values for \eqn{m} are 1,3,5, and 7.
\cr
In case of error, it returns \code{NULL}.
}


\details{
These first and last functions are inverses of each other.

Let \eqn{a} be a set of \eqn{n} strictly disjoint arcs in \eqn{A_n},
and denote the \eqn{n} complementary arcs by \eqn{\bar{a}}.
Let \eqn{\alpha : \mathbb{S}^{2n} \to \mathbb{S}^{2n}}
denote the antipodal map.
Let \eqn{h} denote \code{spherefromarcs()}.
Then \eqn{h( \bar{a} ) = \alpha( h(a) )}.
A fancier way to say this:
the antipodal map \eqn{\alpha} on \eqn{\mathbb{S}^{2n}}
and the complementary map \eqn{ a \mapsto \bar{a} } on \eqn{A_n}
are \emph{conjugate}.

If \eqn{s} is the full circle, then \eqn{h(s)} = (0,...,0,1).
If \eqn{\phi} is the empty arc, then \eqn{h(\phi)} = (0,...,0,-1).
}



\seealso{
\code{\link{boundaryfromarcs}()},
\code{\link{spherefromboundary}()},
\code{\link{boundaryfromsphere}()},
\code{\link{arcsfromboundary}()},
\code{\link{complementaryarcs}()}

}





% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

