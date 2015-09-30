\name{wald.strata}

\alias{wald.strata}

\title{ Wald Test for Stratification Factors }

\description{
  It performs a Wald test to test the hypothesis of equal spline coefficients among strata 
  in the approximation of baseline risk function.
}

\usage{
wald.strata(fit)
} 

\arguments{
  \item{fit}{an object of class \code{jointModel} with \code{method = "spline-PH-GH"}
  and with a strata specification in the survival part.}
}

\value{
  an object of class \code{wald.strata} with components:
  \item{alternative}{a character string naming the alternative.}
  \item{Result}{a numeric matrix with the results of the Wald test.}
}

\note{
This test is valid when the same knots have been used across strata.
}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\examples{
\dontrun{
fitLME <- lme(log(serBilir) ~ drug * year - drug, random = ~ year | id, 
    data = pbc2)
fitSURV <- coxph(Surv(years, status2) ~ drug + strata(hepatomegaly), 
    data = pbc2.id, x = TRUE)
fit.pbc <- jointModel(fitLME, fitSURV, timeVar = "year", method = "spline-PH-aGH")
wald.strata(fit.pbc)
}
}

\keyword{multivariate}
\keyword{regression}

