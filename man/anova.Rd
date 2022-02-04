\name{anova}

\alias{anova.jointModel}

\title{Anova Method for Fitted Joint Models}

\description{
  Produces marginal Wald tests or Performs a likelihood 
	ratio test between two nested joint models.
}

\usage{
\method{anova}{jointModel}(object, object2, test = TRUE, 
   process = c("both", "Longitudinal", "Event"), L = NULL, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}, nested in \code{object2}.}
  \item{object2}{an object inheriting from class \code{jointModel}.}
  \item{test}{logical; if \code{TRUE} the likelihood ratio test is performed.}
  \item{process}{for which of the two submodels to produce the marginal Wald tests table.}
  \item{L}{a numeric matrix of appropriate dimensions defining the contrasts of interest.}
  \item{\dots}{additional arguments; currently none is used.}
}

\value{
  An object of class \code{aov.jointModel} with components,
    \item{nam0}{the name of \code{object}.}
    \item{L0}{the log-likelihood under the null hypothesis (\code{object}).}
    \item{aic0}{the AIC value for the model given by \code{object}.}
    \item{bic0}{the BIC value for the model given by \code{object}. }
    \item{nam1}{the name of \code{object2}.}
    \item{L1}{the log-likelihood under the alternative hypothesis (\code{object2}).}
    \item{aic1}{the AIC value for the model given by \code{object2}.}
    \item{bic1}{the BIC value for the model given by \code{object2}.}
    \item{df}{the degrees of freedom for the test (i.e., the difference in the number of parameters).}
    \item{LRT}{the value of the Likelihood Ratio Test statistic (returned if \code{test = TRUE}).}
    \item{p.value}{the \eqn{p}-value of the test (returned if \code{test = TRUE}).}
    \item{aovTab.Y}{a data.frame with the marginal Wald tests for the longitudinal process; 
	produced only when \code{object2} is missing.}
    \item{aovTab.T}{a data.frame with the marginal Wald tests for the event process; 
	produced only when \code{object2} is missing.}
    \item{aovTab.L}{a data.frame with the marginal Wald tests for the user-defined contrasts matrix; 
	produced only when \code{object2} is missing and \code{L} is not \code{NULL}.}
}

\section{Warning}{
    The code minimally checks whether the models are nested! The user is responsible to supply nested models in
    order the LRT to be valid.
}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2010) JM: An R Package for the Joint Modelling of Longitudinal and Time-to-Event Data. 
\emph{Journal of Statistical Software} \bold{35} (9), 1--33. \doi{10.18637/jss.v035.i09}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{jointModel}}}

\examples{
\dontrun{
# linear mixed model fit without treatment effect
fitLME.null <- lme(sqrt(CD4) ~ obstime, 
    random = ~ 1 | patient, data = aids)
# cox model fit without treatment effect
fitCOX.null <- coxph(Surv(Time, death) ~ 1, 
    data = aids.id, x = TRUE)
# joint model fit without treatment effect
fitJOINT.null <- jointModel(fitLME.null, fitCOX.null, 
    timeVar = "obstime", method = "weibull-PH-aGH")

# linear mixed model fit with treatment effect
fitLME.alt <- lme(sqrt(CD4) ~ obstime * drug - drug, 
    random = ~ 1 | patient, data = aids)
# cox model fit with treatment effect
fitCOX.alt <- coxph(Surv(Time, death) ~ drug, 
    data = aids.id, x = TRUE)
# joint model fit with treatment effect
fitJOINT.alt <- jointModel(fitLME.alt, fitCOX.alt, timeVar = "obstime", 
    method = "weibull-PH-aGH")

# likelihood ratio test for treatment effect
anova(fitJOINT.null, fitJOINT.alt)
}
}

\keyword{methods}
