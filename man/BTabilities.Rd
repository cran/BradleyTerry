\name{BTabilities}
\alias{BTabilities}
\title{ estimated abilities in a Bradley-Terry model }
\description{
   Extracts the \code{abilities} component from a model object of class \code{\link{BTm}}.
}
\usage{
BTabilities(model)
}
\arguments{
  \item{model}{a model object for which \code{inherits(model, "BTm")} is
    \code{TRUE}}
}
\value{
  A two-column numeric matrix, with columns named "ability" and "se";
  one row for each player.
}
\author{ David Firth }
\references{
Firth, D. (2005)  Bradley-Terry models in R.  \emph{Journal of Statistical
Software}, to appear.
}
\seealso{ \code{\link{BTm}}, \code{\link{BTresiduals}} }
\examples{
data(citations)
origin <- factor(c("UK", "USA", "USA", "UK"))
citeModel2 <- BTm(citations ~ origin)
BTabilities(citeModel2)
}
\keyword{ models }
