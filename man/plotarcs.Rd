\name{plotarcs}

\alias{plotarcs}


\title{
plot a collection of arcs
}

\description{
The function \code{plotarcs()} plots a circle
in the plane, with the given arcs overlayed with a thicker line.
}


\usage{
plotarcs( arcmat, labels=FALSE, main=NULL, margintext=NA, rad=1, 
                        lwd=3, pch=20, cex=1.5, add=FALSE, ...  )
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
For improper arcs the endpoints are not drawn.
}

\item{labels}{
if \code{TRUE} then add a text label near the center of
each arc
}

\item{main}{
Text to put above the plot, passed to \code{title()}.
If \code{main=NULL} then a suitable descriptive title is constructed.
}

\item{margintext}{
if not \code{NA} then add this text under the x-axis
}

\item{rad}{
radius of the arcs, which are centered at (0,0).
A smaller radius can be used with \code{add=TRUE} to avoid overlap
with the first set of arcs.
}


\item{lwd}{
line width of the arcs.  The default \code{lwd=3} makes the arcs
bolder to distinguish them from the underlying circle,
which is drawn with \code{lwd=1}.
If\code{lwd=NA} or \code{lwd<0} then drawing is skipped and only the endpoints are drawn.
}

\item{pch}{
symbol for the endpoints, see \code{points()}
}

\item{cex}{
expansion factor for the endpoints, see \code{points()}
}

\item{add}{
when \code{add=TRUE} the arcs are added to an existing plot
}

\item{...}{
extra arguments for the both the arcs and their endpoints,
passed to both \code{lines()} and \code{points()}.
For example, the color \code{col}.
}

}



\value{
\code{plotarcs()} returns \code{TRUE} or \code{FALSE}.
}





