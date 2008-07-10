\name{BTm}
\alias{BTm}
\alias{drop1.BTm}
\alias{add1.BTm}
\alias{terms.BTm}
\alias{formula.BTm}
\title{ Bradley-Terry model and extensions }
\description{
  Fits Bradley-Terry models for pair comparison data, including
  models with structured scores, order effect and missing covariate data.  
  Fits
  by either maximum likelihood or maximum penalized likelihood
  (with Jeffreys-prior penalty).
}
\usage{
BTm(formula, refcat = NULL, offset = NULL, contrasts = NULL,
data = NULL, subset = NULL, br = FALSE, order.effect =
NULL, ...)
}
\arguments{
  \item{formula}{ A model formula for the Bradley-Terry \dQuote{ability}
    parameters.  Response variable should be a data frame containing
    contest results: columns include factors named "winner" (by default
    column 1) and "loser" (by default column 2), and possibly also
    "Freq", a numeric vector recording the
    frequency of each contest result (taken to be 1 if omitted).
    Variables on
    RHS have length equal to the number of players.  The special RHS
    formula \code{..} specifies the standard Bradley-Terry model
    with unstructured abilities.}
  \item{refcat}{ Character.  Which is the \dQuote{reference} player?  Only used with
    \code{..} on the RHS of the formula (otherwise ignored if supplied).
    Default is the name of the alphanumerically-first player.}
  \item{offset}{ An optional offset term in the model. A vector of
    length equal to the number of players.}
  \item{contrasts}{ As for \code{\link{glm}}. }
  \item{data}{ A data frame, in which RHS variables can be found. }
  \item{subset}{  An optional logical or numeric vector specifying a 
  subset of 
  observations (ie, a subset of rows of the response dataframe) 
  to be used in the fitting process. }
  \item{br}{ Logical.  If \code{TRUE}, fitting will be by penalized
    maximum likelihood as in Firth (1992, 1993), using
    \code{\link[brglm]{brglm}}. Default is fitting by maximum likelihood.}
  \item{order.effect}{An optional vector, numeric, indicating an
    order effect to be estimated in the model (for example, a home
    advantage effect).  Values should be 1 where contest winner has
    the advantage, -1 where loser has the advantage, and 0 where neither
    player is advantaged.}
  \item{\dots}{Other arguments for \code{glm} or \code{brglm} }
}
\details{
  No allowance is made for tied contests.

  Aside from the possibility of an order effect, contest-specific
  predictors are not catered for by \code{BTm}.  However, the availability of
  the \code{model} and \code{x} components in objects of class
  \code{BTm} allows a model fitted with only player-specific predictors to be
  manipulated subsequently to include further terms involving
  contest-specific predictors.

}
\value{
  An object of class \code{c("BTm", "glm", "lm")}, or of class
  \code{c("BTm", "brglm", "glm", "lm")} if \code{br = TRUE}.  Components are as
  for \code{glm} or \code{brglm}, with additionally
  \item{x0}{Model matrix for the formula as supplied (rather than
    the model matrix actually used in the subsidiary call to \code{glm}
    or \code{brglm}, which is included as component \code{x} if the call
    includes \code{x = TRUE}).  One row for each player.}
  \item{offset0}{The supplied offset vector, if any was specified.  
    One element for each player.  (The offset vector actually used
    in the subsidiary call to \code{glm}
    or \code{brglm} is included as component \code{offset}.)}
  \item{y0}{The data frame of contest winners and losers, containing 
  only those rows actually used in fitting the model.}
  \item{order.effect}{The values of order.effect, if specified.}
  \item{abilities}{A two-column matrix of estimated abilities, with the
    ability for \code{refcat} set to zero if \code{refcat} is specified
    (otherwise the first player has zero ability).  First column is
    estimated ability, second column is the standard error for that
    estimate.  One row for each player.}
}
\seealso{
  \code{\link{BTresiduals}}, \code{\link{BTabilities}}
  }
\references{

Agresti, A (2002)  \emph{Categorical Data Analysis} (2nd ed).  New York: Wiley.

Firth, D. (1992)  Bias reduction, the Jeffreys prior and GLIM. In 
\emph{Advances in GLIM and Statistical Modelling}, Eds. L Fahrmeir,
B J Francis, R Gilchrist and G Tutz, pp91--100.  New York: Springer.

Firth, D. (1993)  Bias reduction of maximum likelihood estimates.
  \emph{Biometrika} \bold{80}, 27--38.
  
Firth, D. (2005)  Bradley-Terry models in R.  \emph{Journal of Statistical
Software}, to appear.
  
Stigler, S. (1994)  Citation patterns in the journals of statistics 
and probability.  \emph{Statistical Science} \bold{9}, 94--108.
  
}
\author{ David Firth }
\note{
  Methods specific to the \code{BTm} class of models are
  \itemize{
    \item \code{add1.BTm}
    \item \code{drop1.BTm}
    \item \code{terms.BTm}
    \item \code{formula.BTm}
  }
  Others are inherited from \code{glm} or \code{lm}.
}
\examples{
##
##  Statistics journal citation data from Stigler (1994)
##  -- see also Agresti (2002, p448)
data(citations)
                 
##  First fit the "standard" Bradley-Terry model
print(citeModel <- BTm(citations ~ ..))

##  Now the same thing with a different "reference" journal
update(citeModel, . ~ ., refcat = "JASA")

##  Is the "citeability" of a journal predicted by its country of origin?
origin <- factor(c("UK", "USA", "USA", "UK"))
print(citeModel2 <- BTm(citations ~ origin))

##  Hmm... not so sure about the origin of "Comm Statist" ...
is.na(origin[2]) <- TRUE
citeModel2 <- update(citeModel2, . ~ .)

##  Now an example with an order effect -- see Agresti (2002) p438
data(baseball)

##  Simple Bradley-Terry model as in Agresti p437
print(baseballModel <- BTm(baseball ~ ..))

##  Introduce order effect as in Agresti p438
update(baseballModel, order.effect = baseball$home.adv)   

}
\keyword{ models }
