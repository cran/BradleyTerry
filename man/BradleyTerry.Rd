\name{BradleyTerry}
\alias{BTm}
\alias{update.BTm}
\alias{drop1.BTm}
\alias{add1.BTm}
\alias{terms.BTm}
\alias{formula.BTm}
\alias{BTresiduals}
\alias{BTabilities}
\title{ Bradley-Terry model and extensions }
\description{
  Fits Bradley-Terry models for pair comparison data, including
  models with structured scores and missing covariate data.  Fits
  by either maximum likelihood or maximum penalized likelihood
  (with Jeffreys-prior penalty).
}
\usage{
BTm(formula, refcat = NULL, offset = NULL, contrasts = NULL,
data = NULL, weights = NULL, subset = NULL, br = FALSE, order.effect =
NULL, ...)

BTresiduals(model)

BTabilities(model)
}
\arguments{
  \item{formula}{ A model formula for the Bradley-Terry \dQuote{ability}
    parameters.  Response variable should be a data frame containing
    contest results: columns include factors named "winner" (by default
    column 1) and     "loser" (by default column 2), and possibly also
    "Freq", a numeric vector recording the
    frequency of each contest result (taken to be 1 if omitted).
    Variables on
    RHS have length equal to the number of players.  The special RHS
    formula \code{..} specifies the standard Bradley-Terry model
    with unstructured abilities.}
  \item{refcat}{ Which is the \dQuote{reference} player?  Only used with
    \code{..} on the RHS of the formula.  Default is player 1.}
  \item{offset}{ An optional offset term in the model. A vector of
    length equal to the number of players.}
  \item{contrasts}{ As for \code{\link{glm}}. }
  \item{data}{ A data frame, in which RHS variables can be found. }
  \item{weights}{ An optional vector of weights.  Has no effect on the fit. }
  \item{subset}{  An optional logical or numeric vector specifying a subset of observations to be
          used in the fitting process. }
  \item{br}{ Logical.  If \code{TRUE}, fitting will be by penalized
    maximum likelihood as in Firth (1992, 1993), using
    \code{\link[brlr]{brlr}} from
    package brlr. Default is fitting by maximum likelihood.}
  \item{order.effect}{An optional vector, numeric, indicating an
    order effect to be estimated in the model (for example, a home
    advantage effect).  Values should be 1 where contest winner has
    the advantage, -1 where loser has the advantage, and 0 where neither
    player is advantaged.}
  \item{\dots}{other arguments for \code{glm} or \code{brlr} }
  \item{model}{a model object for which \code{inherits(model, "BTm")} is
    \code{TRUE}}
}
\details{
  No allowance is made for tied contests.

  Aside from the possibility of an order effect, contest-specific
  predictors are not catered for by \code{BTm}.  However, the availability of
  the \code{model} and \code{x} components in objects of class
  \code{BTm} allows a model fitted with only player-specific predictors to be
  manipulated subsequently to include further terms involving
  contest-specific predictors.

  The residuals returned by \code{BTresiduals} are weighted means of
  working residuals, with weights equal to the binomial denominators in
  the fitted model.  These are suitable for diagnostic model checking,
  for example plotting against candidate predictors.
}
\value{
  An object of class \code{c("BTm", "glm", "lm")}, or of class
  \code{c("BTm", "brlr", "glm", "lm")} if \code{br = TRUE}.  Components are as
  for \code{glm} or \code{brlr}, with additionally
  \item{x0}{Model matrix for the formula as supplied (rather than
    the model matrix actually used in the subsidiary call to \code{glm}
    or \code{brlr}, which is included as component \code{x} if the call
    includes \code{x = TRUE}).}
  \item{abilities}{A two-column matrix of estimated abilities, with the
    ability for \code{refcat} set to zero if \code{refcat} is specified
    (otherwise the first player has zero ability).  First column is
    estimated ability, second column is the standard error for that
    estmate.  One row for each player.}

  Function \code{BTresiduals} returns a vector of residuals, one for each
  \dQuote{player}.

  Function \code{BTabilities} simply returns the \code{abilities}
  component of the specified model.
  
}
\references{

Agresti, A (2002)  \emph{Categorical Data Analysis} (2nd ed).  New York: Wiley.

Firth, D. (1992)  Bias reduction, the Jeffreys prior and GLIM. In 
\emph{Advances in GLIM and Statistical Modelling}, Eds. L Fahrmeir,
B J Francis, R Gilchrist and G Tutz, pp91--100.  New York: Springer.

Firth, D. (1993)  Bias reduction of maximum likelihood estimates.
  \emph{Biometrika} \bold{80}, 27--38.
  
Stigler, S. (1994)  Citation patterns in the journals of statistics 
and probability.  \emph{Statistical Science} \bold{9}, 94--108.
  
}
\author{ David Firth }
\note{
  Methods specific to the \code{BTm} class of models are
  \itemize{
    \item \code{add1.BTm}
    \item \code{drop1.BTm}
    \item \code{update.BTm}
  }
  Others are inherited from \code{glm} or \code{lm}.
}
\examples{
##
##  Statistics journal citation data from Stigler (1994)
##  -- see also Agresti (2002, p448)
journals <- c("Biometrika", "Comm Statist", "JASA", "JRSS-B")
citedata <- matrix(c(NA,33,320,284,
                     730,NA,813,276,
                     498,68,NA,325,
                     221,17,142,NA),
                   4,4,
                   dimnames = list(winner = journals, loser = journals))
citedata <- as.data.frame.table(citedata)
                 
##  First fit the "standard" Bradley-Terry model
print(mymodel <- BTm(citedata ~ ..))

##  Now the same thing with a different "reference" journal
print(mymodel2 <- update(mymodel, . ~ ., refcat = "JASA"))

##  Is the "citeability" of a journal predicted by its country of origin?
origin <- factor(c("UK", "USA", "USA", "UK"))
print(mynewmodel <- BTm(citedata ~ origin))

##  Hmm... not so sure about the origin of "Comm Statist" ...
is.na(origin[2]) <- TRUE
mynewmodel2 <- update(mynewmodel, . ~ .)
BTabilities(mynewmodel2)

##  Now an example with an order effect -- see Agresti (2002) p438
teams <- c("Milwaukee", "Detroit", "Toronto",
           "New York", "Boston", "Cleveland", "Baltimore")
baseball <- data.frame(
  winner = factor(rep(rep(1:7, rep(7,7)), 2), labels =  teams),
  loser = factor(rep(rep(1:7, 7), 2), labels = teams),
  Freq = c(NA,4,4,4,6,4,6,  # home wins
           3,NA,4,4,6,6,4,
           2,4,NA,2,4,4,6,
           3,5,2,NA,4,4,6,
           5,2,3,4,NA,5,6,
           2,3,3,4,4,NA,2,
           2,1,1,2,1,3,NA,
           NA,3,5,3,1,5,5,  # away wins
           3,NA,3,1,5,3,5,
           2,2,NA,5,3,4,6,
           3,3,4,NA,2,3,4,
           1,0,3,3,NA,2,6,
           2,1,2,2,2,NA,4,
           0,3,0,1,0,4,NA),
   home.adv = rep(c(1, -1), c(49,49)))

##  Simple Bradley-Terry model as in Agresti p437
print(baseball.model <- BTm(baseball ~ ..))
BTresiduals(baseball.model)

##  Introduce order effect as in Agresti p438
update(baseball.model, order.effect = baseball$home.adv)   

}
\keyword{ models }
