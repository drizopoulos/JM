\name{ranef}

\alias{ranef.jointModel}

\title{Random Effects Estimates for Joint Models}

\description{
Extracts the random effects estimates from a fitted joint model.
}

\usage{
\method{ranef}{jointModel}(object, type = c("mean", "mode"), postVar = FALSE, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{type}{what type of empirical Bayes estimates to compute, the mean of the posterior distribution or
	the mode of the posterior distribution.}
  \item{postVar}{logical; if \code{TRUE} the variance of the posterior distribution is also returned. When
	\code{type == "mode"}, then this equals \eqn{\{- \partial^2 \log p(b_i | T_i, 
	\delta_i, y_i) / \partial b_i^\top \partial b_i \}^{-1}}.}
  \item{\dots}{additional arguments; currently none is used.}
}

\value{
  a numeric matrix with rows denoting the individuals and columns the random effects (e.g., intercepts, slopes, etc.). 
  If \code{postVar = TRUE}, the numeric matrix has an extra attribute ``postVar".
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.
}

\seealso{\code{\link{coef.jointModel}}, \code{\link{fixef.jointModel}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ 1 | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")
ranef(fitJOINT)
}
}

\keyword{methods}
