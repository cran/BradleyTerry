\name{baseball}
\alias{baseball}
\docType{data}
\title{ Baseball data from Agresti (2002) }
\description{
  Baseball results for games in the 1987 season between 7 teams
  in the Eastern Division of the American League.
}
\usage{data(baseball)}
\format{
  A data frame with 98 observations on the following 4 variables.
  \describe{
    \item{\code{winner}}{a factor with levels \code{Milwaukee} \code{Detroit} \code{Toronto} \code{New York} \code{Boston} \code{Cleveland} \code{Baltimore}}
    \item{\code{loser}}{a factor with levels \code{Milwaukee} \code{Detroit} \code{Toronto} \code{New York} \code{Boston} \code{Cleveland} \code{Baltimore}}
    \item{\code{Freq}}{a numeric vector}
    \item{\code{home.adv}}{a numeric vector}
  }
}
\details{
  The \code{home.adv} variable is \code{1} for games where the home team won, 
  and \code{-1} for games where the away team won.
}
\source{
page 438 of Agresti, A (2002) \emph{Categorcal Data Analysis} (2nd Edn.).  
Wiley, New York.
}
\references{
Firth, D. (2005)  Bradley-Terry models in R.  \emph{Journal of Statistical
Software}, to appear.
}
\seealso{\code{\link{BTm}}}
\examples{
data(baseball)
##  The data in collapsed tabular form as on p438 of Agresti
xtabs(Freq ~ winner + loser, baseball)
##  Simple Bradley-Terry model as in Agresti p437
print(baseballModel <- BTm(baseball ~ ..))
##  Introduce order effect as in Agresti p438
update(baseballModel, order.effect = baseball$home.adv) 
}
\keyword{datasets}
