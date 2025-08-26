\name{slerp}

\alias{slerp}

\title{
Spherical Linear Interpolation (slerp) on the Unit Sphere
}

\description{
Let \eqn{u_0} and \eqn{u_1} be two unit vectors on the unit sphere \eqn{\mathbb{S}^{m-1}}, that are not antipodal.

The primary purpose of \code{slerp()} is to compute uniformly spaced
points along the great circle arc from \eqn{u_0} to \eqn{u_1}.
The angle step \eqn{\Delta \theta} is computed to be
as large as possible and also satisfy \eqn{\Delta \theta \le \theta_{max}},
where \eqn{\theta_{max}} is a positive user-given parameter.
If the distance between \eqn{u_0} and \eqn{u_1} is \eqn{\theta}, then
\eqn{\Delta \theta = \theta / \text{ceiling}( \theta / \theta_{max})}.

The secondary purpose is to compute points on the great circle arc
from user supplied interpolation parameters \eqn{\tau} = a numeric
vector with values typically in [0,1].
The value 0 generates \eqn{u_0} and 1 generates \eqn{u_1}.
For example to generate \code{n} uniformly space points,
set \eqn{\tau} = \code{ (0:(n-1))/(n-1) }.
}

\usage{
slerp( u0, u1, thetamax=pi/36, tau=NULL )
}

\arguments{
\item{u0, u1}{
numeric vectors with same length (verified) and unit length (not verified).
If they are antipodal, the great circle arc is undefined and the
function returns \code{NULL}.
}

\item{thetamax}{
the maximum angular step between successive points on the arc.
This value must be positive.
The default corresponds to 5\eqn{^\circ}.
}

\item{tau}{
a numeric vector of interpolation parameters.
The default is \code{NULL} which means to ignore it and use \code{thetamax}.
If \code{tau} is not \code{NULL}, then \code{thetamax} is ignored.
}
}



\value{
\code{slerp()} returns an \code{n x m} matrix with computed points
in the rows.
\code{n} is the number of points, and
\code{m} is the length of \code{u0} and \code{u1}.


When using \code{thetamax},
\code{n = theta / ceiling(theta/thetamax) + 1},
where \code{theta} is the angle between \code{u0} and \code{u1}.
Buf if \code{u0} and \code{u1} are equal, then \code{n=1}.

When using \code{tau}, \code{n=length(tau)}.
But if \code{tau} is empty, the function returns \code{NULL}.

In case of error, e.g. \code{u0} and \code{u1} are antipodal,
the function returns \code{NULL}.
}


\details{
If \bold{both} \code{thetamax} and \code{tau} are not \code{NULL},
\code{tau} has priority and \code{thetamax} is ignored.

The interpolation formula uses the \code{sin()} function.
For details see \bold{Shoemake}.
}

\references{
Shoemake, Ken.
\emph{Animating Rotation with Quaternion Curves}.
pp. 245-254.
SIGGRAPH 1985.
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

