\name{BTresiduals}
\alias{BTresiduals}
\title{ player-specific residuals from a Bradley-Terry model }
\description{
  Computes player-specific residuals from a model object of class
  \code{\link{BTm}}, suitable for diagnostic checking of a predictor involving
  player-level covariates.
}
\usage{
BTresiduals(model)
}
\arguments{
  \item{model}{a model object for which \code{inherits(model, "BTm")} is
    \code{TRUE}}
}
\details{
  The residuals returned by \code{BTresiduals} are weighted means of
  working residuals, with weights equal to the binomial denominators in
  the fitted model.  These are suitable for diagnostic model checking,
  for example plotting against candidate predictors.
  }
\value{
  A numeric vector of length equal to the number of players, with a
  "weights" attribute.
}
\author{ David Firth }
\references{
Firth, D. (2005)  Bradley-Terry models in R.  \emph{Journal of Statistical
Software}, to appear.
}
\seealso{ \code{\link{BTm}}, \code{\link{BTabilities}} }
\examples{
data(citations)
origin <- factor(c("UK", "USA", "USA", "UK"))
citeModel2 <- BTm(citations ~ origin)
BTresiduals(citeModel2)
}
\keyword{ models }
