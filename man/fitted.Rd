\name{fitted}

\alias{fitted.jointModel}

\title{Fitted Values for Joint Models}

\description{
  Calculates fitted values for joint models.
}

\usage{
\method{fitted}{jointModel}(object, process = c("Longitudinal", "Event"), 
type = c("Marginal", "Subject", "EventTime", "Slope"), scale = c("survival", 
"cumulative-Hazard", "log-cumulative-Hazard"), M = 200, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{process}{for which model (i.e., linear mixed model or survival model) to calculate the fitted values.}
  \item{type}{what type of fitted values to calculate for the survival outcome. See \bold{Details}.}
  \item{scale}{in which scale to calculate; relevant only when \code{process = "Event"}.}
  \item{M}{how many times to simulate random effects; see \bold{Details} for more info.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  For \code{process = "Longitudinal"}, let \eqn{X} denote the design matrix for the fixed effects \eqn{\beta}, and 
  \eqn{Z} the design matrix for the random effects \eqn{b}. Then for \code{type = "Marginal"} the fitted values are 
  \eqn{X \hat{\beta},} whereas for \code{type = "Subject"} they are \eqn{X \hat{\beta} + Z \hat{b}}. For \code{type = "EventTime"} 
  is the same as \code{type = "Subject"} but evaluated at the observed event times. Finally, \code{type == "Slope"}
  returns \eqn{Xs \hat{\beta} + Zs \hat{b}} where \eqn{Xs} and \eqn{Zs} denote the fixed- and random-effects design
  matrices corresponding to the slope term which is specified in the \code{derivForm} argument of \code{\link{jointModel}}.
  
  For \code{process = "Event"} and \code{type = "Subject"} the linear predictor conditional on the random effects 
  estimates is calculated for each sample unit. Depending on the value of the \code{scale} argument the fitted survival 
  function, cumulative hazard function or log cumulative hazard function is returned. For \code{type = "Marginal"}, 
  random effects values for each sample unit are simulated \code{M} times from a normal distribution with zero mean and 
  covariance matrix the estimated covariance matrix for the random effects. The marginal survival function for the 
  \eqn{i}th sample unit is approximated by \deqn{S_i(t) = \int S_i(t | b_i) p(b_i) db_i \approx (1/M) \sum_{m = 1}^M 
  S_i(t | b_{im}),} where \eqn{p(b_i)} denotes the normal probability density function, and \eqn{b_{im}} the \eqn{m}th 
  simulated value for the random effect of the \eqn{i}th sample unit. The cumulative hazard and log cumulative hazard 
  functions are calculated as \eqn{H_i(t) = - \log S_i(t)} and \eqn{\log H_i(t) = \log \{ - \log S_i(t)\},} respectively.
}

\value{
  a numeric vector of fitted values.
}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2010) JM: An R Package for the Joint Modelling of Longitudinal and Time-to-Event Data. 
\emph{Journal of Statistical Software} \bold{35} (9), 1--33. \url{http://www.jstatsoft.org/v35/i09/}
}


\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{residuals.jointModel}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * year, 
    random = ~ 1 | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, 
    data = pbc2.id, x = TRUE)
# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")

# fitted for the longitudinal process
head(cbind(
    "Marg" = fitted(fitJOINT), 
    "Subj" = fitted(fitJOINT, type = "Subject")
))

# fitted for the event process - survival function
head(cbind(
    "Marg" = fitted(fitJOINT, process = "Ev"), 
    "Subj" = fitted(fitJOINT, process = "Ev", type = "Subject")
))

# fitted for the event process - cumulative hazard function
head(cbind(
    "Marg" = fitted(fitJOINT, process = "Ev", 
        scale = "cumulative-Hazard"), 
    "Subj" = fitted(fitJOINT, process = "Ev", type = "Subject", 
        scale = "cumulative-Hazard")
))
}
}

\keyword{methods}
