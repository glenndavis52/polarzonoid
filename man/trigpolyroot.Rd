\name{trigpolyroot}

\alias{trigpolyroot}


\title{
compute roots of a trigonometric polynomial
}

\description{
Given a trigonometric polynomial of degree \eqn{n},
defined on the real interval \eqn{[0,2\pi]}.
The function transforms it to a polynomial
of degree \eqn{2n} in the standard way,
and then solves that polynomial.
Only the real roots of the original trigonometric polynomial are returned.
}

\usage{
trigpolyroot( ablist, tol=.Machine$double.eps^0.5 )
}
\arguments{
\item{ablist}{
a list with 3 items
\tabular{ll}{
\code{a0} \tab the constant term of the trigonometric polynomial \cr
\code{a} \tab real coefficients of \eqn{cos(jt), ~ j = 1..n }\cr
\code{b} \tab real coefficients of \eqn{sin(jt), ~ j = 1..n }
}
}

\item{tol}{
tolerance for the imaginary part; if the absolute value of the
imaginary part of the root is less than or equal to \code{tol},
then the root is considered to be real, and is returned
}
}


\details{
The given trigonometric polynomial is converted,
using a change of variable and Euler's formula,
to a standard polynomial of degree \eqn{2n}, see \bold{Wikipedia}.
This polynomial is solved using \code{base::polyroot()}.
After the inverse change of variable,
the imaginary roots are removed using \code{tol}.
}


\value{
\code{trigpolyroot()} returns the real roots in the interval
\eqn{[0,2\pi)} and in increasing order (as a numeric vector).
There may be duplicates, which correspond to multiple roots.
It is up to the user to handle multiple roots.

The number of roots returned is always even, and in the set
\eqn{\{ 0,2,4,..,2n \}}.

If \eqn{n} is 0, or all the coefficients are 0, it is an error.
In case of error, the function returns \code{NULL}.
}


\source{
Wikipedia.
\bold{Trigonometric polynomial}.
\url{https://en.wikipedia.org/wiki/Trigonometric_polynomial}
}


\seealso{
\code{\link[base]{polyroot}()}
}





% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

