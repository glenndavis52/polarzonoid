\name{rotation}

\alias{arcsfromrotation}
\alias{rotationfromarcs}
\alias{rotationaroundaxis}


\title{
The homeomorphism between balanced configurations and rotation matrices
}

\description{
In the
\href{../doc/rotation.html}{\bold{Real Projective Spaces and 3x3 Rotation Matrices}}
vignette,
it is shown that there is a homeomorphism between \eqn{P_4^\pi},
the space of balanced configurations of 4 or 2 points on the circle,
and \eqn{\text{SO(3)}}, the space of 3x3 rotation matrices.
These functions implement the homeomorphim and its inverse.
}

\usage{
rotationfromarcs( arcmat, tol=5.e-9 )
arcsfromrotation( rotation, tol=1.e-7 )
rotationaroundaxis( axis, theta )
}

\arguments{
\item{arcmat}{
a \code{2x2} or \code{1x2} matrix with an arc definition in each row;
so the number of arcs is 2 or 1.
The 1st number in the row is the \emph{center} of the arc,
and the 2nd number is the \emph{length} of the arc;
both in radians.\cr
The total length must be \eqn{\pi}, and this is tested.
If there is only one arc, it must be a semicircle.
Only the endpoints of \code{arcmat} matter;
the complement of the arcs maps to the same rotation.
}

\item{tol}{
In \code{rotationfromarcs()} \code{tol} is the tolerance for the test that the total
arc length is \eqn{\pi}.\cr
In \code{arcsfromrotation()} \code{tol} is the tolerance for testing that the mapped point
in \eqn{\partial Z_4} is close to the substratum  \eqn{\partial Z_2}.
This is to avoid generating 2 arcs, one of which is very tiny.
}

\item{rotation}{
a 3x3 rotation matrix.
The numeric 3x3 matrix property is verified, but, in the interest of speed,
orthogonality and orientation-preservation is not verified.
}

\item{axis}{
a non-zero 3-vector defining an oriented axis in \eqn{\mathbb{R}^3}.
\code{axis} is not required to be a unit vector.
}

\item{theta}{
an angle in radians.  
The computed rotation matrix rotates by this angle.
}
}



\value{
\code{rotationfromarcs()} maps from \eqn{P_4^\pi} to \eqn{\text{SO(3)}};
it returns a 3x3 rotation matrix.

\code{arcsfromrotation()} maps from \eqn{\text{SO(3)}} to \eqn{P_4^\pi};
the returned set of arcs (either 1 or 2 of them) is only defined up to complementation.
Since all we want is the endpoints, the ambiguity does not matter.
The total length of the arcs is \eqn{\pi}.

\code{rotationaroundaxis()} returns a 3x3 rotation matrix
which rotates \code{theta} radians around the given \code{axis}.

In case of error, all these functions return \code{NULL}.
}


\details{
As currently implemented, a 2-point configuration in \eqn{P_2^\pi},
whose points must be antipodal, maps to a rotation around the x-axis.
The configuration of the 2 points \eqn{(0,\pm 1)} maps to the 3x3 identity matrix.
See the
\href{../doc/rotation.html}{\bold{Real Projective Spaces and 3x3 Rotation Matrices}}
vignette,
for an animation of this.
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

