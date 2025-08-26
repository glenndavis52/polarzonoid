\name{support}

\alias{support}

\title{
The Support Function for the Polar Zonoid in \eqn{\mathbb{R}^{2n+1}}
}

\description{
This section calculates the classical support function
for convex bodies, for the specific case of the polar zonoid.
}

\usage{
support( direction )
}

\arguments{
\item{direction}{
an \code{Rx(2n+1)} matrix with non-zero direction vectors in the
\code{R} rows.
If the number of columns is even, a column of 0s is appended
to make it odd.
\code{direction} can also be a vector of length \code{M},
which is then converted to a \code{1xM} matrix.
}

}



\value{
\code{support()} returns a data frame with \code{R} rows and these columns:
  \item{direction}{the given matrix of directions}
  \item{value}{the value of the support function of the zonoid, in the given direction}
  \item{argmax}{the point on the boundary of the zonoid where the maximum is taken}
  \item{arcs}{the number of arcs for \code{argmax}. This number is always \eqn{\le} \code{n}.}
\cr
In case of error, the function returns \code{NULL}.
}


\details{
Each direction is converted to a non-zero trignometric polynomial.
The roots are computed using \code{trigpolyroot()} and the the roots
are converted to arcs, and thus \code{arcs}.
The function \code{boundaryfromarcs()} is then used to compute
\code{argmax}, and then \code{value}.
}


\seealso{
\code{\link{boundaryfromarcs}()},
\code{\link{trigpolyroot}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

