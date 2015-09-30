\name{summary.weibull.frailty}

\alias{summary.weibull.frailty}

\title{ Summary Method for weibull.frailty Objects}

\description{
  Summarizes the fit of a Weibull model with Gamma frailties
}

\usage{
\method{summary}{weibull.frailty}(object, sand.se = FALSE, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{weibull.frailty}.}
  \item{sand.se}{logical; if \code{TRUE}, sandwich standard errors are also produced.}
  \item{\dots}{ additional arguments; currently none is used.}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{weibull.frailty}}}

\examples{
fit <- weibull.frailty(Surv(time, status) ~ age + sex, kidney)
summary(fit)
summary(fit, TRUE)
}

\keyword{methods}

