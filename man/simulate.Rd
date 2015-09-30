\name{simulate}

\alias{simulateJM}
\alias{simulate.jointModel}

\title{Simulate from Joint Models.}

\description{
  simulate longitudinal responses and event times from joint models
}

\usage{
simulateJM(nsim, nsub, thetas, times, formulas, Data = NULL, 
    method = c("weibull-PH", "weibull-AFT", "piecewise-PH", "spline-PH"), 
    lag = 0, censoring = "uniform", max.FUtime = NULL, seed = NULL, 
    return.ranef = FALSE)
\method{simulate}{jointModel}(object, nsim, seed = NULL, times = NULL, 
    Data = NULL, \dots)
}

\arguments{
  \item{nsim}{number of data sets to be simulated.}
  \item{nsub}{the number of subjects in each data set.}
  \item{thetas}{a list with the parameter values. This should be of the same structure as 
    the \code{coefficients} component returned by \code{jointModel()}.}
  \item{times}{a numeric vector denoting the time points at which longitudinal measurements
    are planned to be taken.}
  \item{formulas}{a list with components: \code{Yfixed} a formula for the fixed-effects part
    of the linear mixed model, \code{Yrandom} a formula for the random-effects part of the
    linear mixed model, \code{Tfixed} a formula for the baseline covariates part of the 
    survival submodel, \code{timeVar} a character string indicating the name of the time 
    variable in the linear mixed model.}
  \item{Data}{a data frame containing any covariates used in the formulas defined in the
    \code{formulas} argument.}
  \item{method}{a character string indicating from what type of survival submodel to simulate.
    There are the same options as the ones provided by \code{\link{jointModel}}.}
  \item{lag}{a numeric value denoting a lagged effect; the same as the \code{lag} 
    argument of \code{\link{jointModel}}.}
  \item{censoring}{a character string or a numeric vector.}
  \item{max.FUtime}{a numeric value denoting the maximum follow-up time for the study. The default
    is \code{max(times) + 2 * IQR(times)}.}
  \item{seed}{an object specifying if and how the random number generator should 
    be initialized ('seeded'). It could be either \code{NULL} or an integer that 
    will be used in a call to \code{set.seed()} before simulating the response 
    vectors. If set, the value is saved as the "seed" attribute of the returned value.}
  \item{return.ranef}{logical; if \code{TRUE}, each component of the returned list has the attributed
    \code{"ranef"} that contains the random-effects values used in the simulation.}
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{\dots}{additional arguments; currently none is used.}
}

\value{A list of length \code{nsim} of data frames that contains the simulated responses 
for the longitudinal process "y", the simulated event times "Time", the event indicator 
"Event", and the subject identification number "id". If extra covariates were assumed, 
these are also included.}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{jointModel}}}

\examples{
\dontrun{
prothro$t0 <- as.numeric(prothro$time == 0)
lmeFit <- lme(pro ~ treat * (time + t0), random = ~ time | id, data = prothro)
survFit <- coxph(Surv(Time, death) ~ treat, data = prothros, x = TRUE)
jointFit <- jointModel(lmeFit, survFit, timeVar = "time", 
    method = "weibull-PH-aGH")

newData <- simulate(jointFit, nsim = 1, times = seq(0, 11, len = 15))
newData
}
}

\keyword{ methods }
