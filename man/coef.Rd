\name{coef}

\alias{coef.jointModel}
\alias{fixef.jointModel}

\title{Estimated Coefficients for Joint Models}

\description{
  Extracts estimated coefficients from fitted joint models.
}

\usage{
\method{coef}{jointModel}(object, process = c("Longitudinal", "Event"), 
    include.splineCoefs = FALSE, \dots)
\method{fixef}{jointModel}(object, process = c("Longitudinal", "Event"), 
    include.splineCoefs = FALSE, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{process}{for which model (i.e., linear mixed model or survival model) to extract the estimated 
    coefficients.}
  \item{include.splineCoefs}{logical; if \code{TRUE} and the method argument in \code{jointModel()} is 
    \code{"ch-Laplace"}, the estimated B-spline coefficients are included as well.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  When \code{process = "Event"} both methods return the same output. However, for \code{process = "Longitudinal"},
  the \code{coef()} method returns the subject-specific coefficients, whereas \code{fixef()} only the fixed effects.
}

\value{A numeric vector or a matrix of the estimated parameters for the fitted model.}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{ranef.jointModel}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(sqrt(CD4) ~ obstime * drug - drug, 
    random = ~ 1 | patient, data = aids)
# cox model fit
fitCOX <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# joint model fit
fitJOINT <- jointModel(fitLME, fitCOX, 
    timeVar = "obstime")

# fixed effects for the longitudinal process
fixef(fitJOINT)

# fixed effects + random effects estimates for the longitudinal 
# process
coef(fitJOINT)

# fixed effects for the event process
fixef(fitJOINT, process = "Event")
coef(fitJOINT, process = "Event")
}
}

\keyword{ methods }
