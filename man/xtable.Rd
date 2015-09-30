\name{xtable}

\alias{xtable.jointModel}

\title{xtable Method from Joint Models.}

\description{
  produces a LaTeX table with the results of a joint model using package xtable.
}

\usage{
\method{xtable}{jointModel}(x, caption = NULL, label = NULL, align = NULL,
    digits = NULL, display = NULL, which = c("all", "Longitudinal", "Event"),
    varNames.Long = NULL, varNames.Event = NULL, p.values = TRUE, 
    digits.pval = 4, \dots)
}

\arguments{
  \item{x}{an object inheriting from class \code{jointModel}.}
  \item{caption}{the \code{caption} argument of \code{xtable()}.}
  \item{label}{the \code{label} argument of \code{xtable()}.}
  \item{align}{the \code{align} argument of \code{xtable()}.}
  \item{digits}{the \code{digits} argument of \code{xtable()}.}
  \item{display}{the \code{display} argument of \code{xtable()}.}
  \item{which}{a character string indicating which results to include in the
      LaTeX table. Options are all results, the results of longitudinal
      submodel or the results of the survival submodel.}
  \item{varNames.Long}{a character vector of the variable names for the
      longitudinal submodel.}
  \item{varNames.Event}{a character vector of the variable names for the
      survival submodel.}
  \item{p.values}{logical; should p-values be included in the table.}
  \item{digits.pval}{a numeric scalare denoting the number of significance
      digits in the \eqn{p}-value.}
  \item{\dots}{additional arguments; currently none is used.}
}

\value{A LaTeX code chunk with the results of the joint modeling analysis.}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{jointModel}}}

\examples{
\dontrun{
require(xtable)
prothro$t0 <- as.numeric(prothro$time == 0)
lmeFit <- lme(pro ~ treat * (time + t0), random = ~ time | id, data = prothro)
survFit <- coxph(Surv(Time, death) ~ treat, data = prothros, x = TRUE)
jointFit <- jointModel(lmeFit, survFit, timeVar = "time", 
    method = "weibull-PH-aGH")

xtable(jointFit, math.style.negative = TRUE)
}
}

\keyword{ methods }
