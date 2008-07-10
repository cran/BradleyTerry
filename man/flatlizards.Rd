\name{flatlizards}
\alias{flatlizards}
\docType{data}
\title{ Augrabies Male Flat Lizards: Contest Results and Predictor
  Variables }
\description{
  Data collected at Augrabies Falls National Park (South Africa) in
  September-October 2002, on the contest performance and background
  attributes of 77 male flat lizards (\emph{Platysaurus broadleyi}).  The
  results of exactly 100 contests were recorded, along with various
  measurements made on each lizard.  Full details of the study are in
  Whiting et al. (2006).
}
\usage{data(flatlizards)}
\format{
  This dataset is a list containing two data frames: \code{flatlizards$contests}
  and \code{flatlizards$predictors}.

  The \code{flatlizards$contests} data frame has 100 observations on the following 2 variables:
  \describe{
    \item{\code{winner}}{a factor with 77 levels \code{lizard003} ... \code{lizard189}}
    \item{\code{loser}}{a factor with the same 77 levels \code{lizard003} ... \code{lizard189}}
  }

  The \code{flatlizards$predictors} data frame has 77 observations (one
  for each of the 77 lizards) on the following 18 variables:
  \describe{
    \item{\code{id}}{an integer identifier for each lizard}
    \item{\code{throat.PC1}}{numeric, the first principal component of
  the throat spectrum}
    \item{\code{throat.PC2}}{numeric, the second principal component of
  the throat spectrum}
    \item{\code{throat.PC3}}{numeric, the third principal component of
  the throat spectrum}
    \item{\code{frontleg.PC1}}{numeric, the first principal component of
  the front-leg spectrum}
    \item{\code{frontleg.PC2}}{numeric, the second principal component of
  the front-leg spectrum}
    \item{\code{frontleg.PC3}}{numeric, the third principal component of
  the front-leg spectrum}
    \item{\code{badge.PC1}}{numeric, the first principal component of
  the ventral colour patch spectrum}
    \item{\code{badge.PC2}}{numeric, the second principal component of
  the ventral colour patch spectrum}
    \item{\code{badge.PC3}}{numeric, the third principal component of
  the ventral colour patch spectrum}
    \item{\code{badge.size}}{numeric, a measure of the area of the
  ventral colour patch}
    \item{\code{testosterone}}{numeric, a measure of blood testosterone concentration}
    \item{\code{SVL}}{numeric, the snout-vent length of the lizard}
    \item{\code{head.length}}{numeric, head length}
    \item{\code{head.width}}{numeric, head width}
    \item{\code{head.height}}{numeric, head height}
    \item{\code{condition}}{numeric, a measure of body condition}
    \item{\code{repro.tactic}}{a factor indicating reproductive tactic; 
  levels are \code{resident} and \code{floater}}
  }
}
\details{
  There were no duplicate contests (no pair of lizards was seen fighting
  more than once), and there were no tied contests (the result of each
  contest was clear).

  The variables \code{head.length}, \code{head.width},
  \code{head.height} and \code{condition} were all computed as
  residuals (of directly measured head length, head width, head height
  and body mass index, respectively) from simple least-squares regressions
  on \code{SVL}.

  Values of some predictors are missing (\code{NA}) for some lizards,
  `at random', because of instrument problems unconnected with the value of the
  measurement being made.
}
\source{
  The data were collected by Dr Martin Whiting, University of the Witwatersrand,
  \url{http://web.wits.ac.za/Academic/Science/APES/Research/MWLab/People/},
  and they appear here with his kind permission.
}
\references{
 Whiting, M.J., Stuart-Fox, D.M., O'Connor, D., Firth, D., Bennett,
 N.C. and Blomberg, S.P. (2006). Ultraviolet signals ultra-aggression in
  a lizard. \emph{Animal Behaviour} \bold{72}, 353-363. (\url{http://dx.doi.org/10.1016/j.anbehav.2005.10.018})
}
\examples{
library(BradleyTerry)
data(flatlizards)
##
##  Fit the standard Bradley-Terry model, using the bias-reduced
##  maximum likelihood method:
##
BTmodel <- BTm(flatlizards$contests ~ .., br = TRUE)
summary(BTmodel)
##
##  That's fairly useless, though, because of the rather small
##  amount of data on each lizard.  And really the scientific
##  interest is not in the abilities of these particular 77
##  lizards, but in the relationship between ability and the
##  measured predictor variables.
##
##  So next fit (by maximum likelihood) a "structured" B-T model in
##  which abilities are determined by a linear predictor.  This 
##  reproduces results reported in Table 1 of Whiting et al. (2006).
##
BTmodel2 <- BTm(flatlizards$contests ~ throat.PC1 + throat.PC3 +
  head.length + SVL, data = flatlizards$predictors)
summary(BTmodel2)
}
\keyword{datasets}
