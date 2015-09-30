\name{survfitJM}

\alias{survfitJM}
\alias{survfitJM.jointModel}

\title{Prediction in Joint Models}

\description{ 
  This function computes the conditional probability of surviving later times than the last observed time for which a
  longitudinal measurement was available.
}

\usage{
survfitJM(object, newdata, \dots)

\method{survfitJM}{jointModel}(object, newdata, idVar = "id", simulate = TRUE, survTimes = NULL, 
	last.time = NULL, M = 200, CI.levels = c(0.025, 0.975), scale = 1.6, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{newdata}{a data frame that contains the longitudinal and covariate information for the subjects for which prediction 
    of survival probabilities is required. The names of the variables in this data frame must be the same as in the data frames that
    were used to fit the linear mixed effects model (using \code{lme()}) and the survival model (using \code{coxph()} or \code{survreg()})
    that were supplied as the two first argument of \code{\link{jointModel}}. In addition, this data frame should contain a variable 
    that identifies the different subjects (see also argument \code{idVar}).}
  \item{idVar}{the name of the variable in \code{newdata} that identifies the different subjects.}
  \item{simulate}{logical; if \code{TRUE}, a Monte Carlo approach is used to estimate survival probabilities. If \code{FALSE},
    a first order estimator is used instead. (see \bold{Details})}
  \item{survTimes}{a numeric vector of times for which prediction survival probabilities are to be computed.}
  \item{last.time}{a numeric vector or character string. This specifies the known time at which each of the subjects in \code{newdat}
    was known to be alive. If \code{NULL}, then this is automatically taken as the last time each subject provided a longitudinal
    measurement. If a numeric vector, then it is assumed to contain this last time point for each subject. If a character string, then
    it should be a variable in the data frame \code{newdata}.}
  \item{M}{integer denoting how many Monte Carlo samples to use -- see \bold{Details}.}
  \item{CI.levels}{a numeric vector of length two that specifies which quantiles to use for the calculation of confidence interval for the 
    predicted probabilities -- see \bold{Details}.}
  \item{scale}{a numeric scalar that controls the acceptance rate of the Metropolis-Hastings algorithm -- see \bold{Details}.}
   \item{\dots}{additional arguments; currently none is used.}
}

\details{
  Based on a fitted joint model (represented by \code{object}), and a history of longitudinal responses 
  \eqn{\tilde{y}_i(t) = \{y_i(s), 0 \leq s \leq t\}}{tilde{y_i}(t) = {y_i(s), 0 \leq s \leq t}} and a covariates vector \eqn{x_i} (stored in
  \code{newdata}), this function provides estimates of \eqn{Pr(T_i > u | T_i > t, \tilde{y}_i(t), x_i)}{Pr(T_i > u | T_i > t, 
  tilde{y}_i(t), x_i)}, i.e., the conditional probability of surviving time \eqn{u} given that subject \eqn{i}, with covariate information
  \eqn{x_i}, has survived up to time \eqn{t} and has provided longitudinal the measurements \eqn{\tilde{y}_i(t)}{tilde{y}_i(t)}.
  
  To estimate \eqn{Pr(T_i > u | T_i > t, \tilde{y}_i(t), x_i)}{Pr(T_i > u | T_i > t, tilde{y}_i(t), x_i)} and if \code{simulate = TRUE}, a 
  Monte Carlo procedure is followed with the following steps:
  \describe{
    \item{Step 1:}{Simulate new parameter values, say \eqn{\theta^*}, from \eqn{N(\hat{\theta}, C(\hat{\theta}))},
        where \eqn{\hat{\theta}} are the MLEs and \eqn{C(\hat{\theta})} their large sample covariance matrix, which are extracted from
        \code{object}.}
    \item{Step 2:}{Simulate random effects values, say \eqn{b_i^*}, from their posterior distribution given survival up to time \eqn{t},
        the vector of longitudinal responses \eqn{\tilde{y}_i(t)} and \eqn{\theta^*}. This is achieved using a Metropolis-Hastings algorithm with
        independent proposals from a properly centered and scaled multivariate \eqn{t} distribution. The \code{scale} argument controls the 
        acceptance rate for this algorithm.}
    \item{Step 3}{Using \eqn{\theta^*} and \eqn{b_i^*}, compute \eqn{Pr(T_i > u | T_i > t, b_i^*, x_i; \theta^*)}{Pr(T_i > 
        u | T_i > t, b_i^*, x_i; \theta^*)}.}
    \item{Step 4:}{Repeat Steps 1-3 \code{M} times.}
  }
  Based on the \code{M} estimates of the conditional probabilities, we compute useful summary statistics, such as their mean, median, and
  quantiles (to produce a confidence interval).

  If \code{simulate = FALSE}, then survival probabilities are estimated using the formula \deqn{Pr(T_i > u | T_i > t, \hat{b}_i, x_i; 
    \hat{\theta}),}{Pr(T_i > u | T_i > t, hat{b}_i, x_i; hat{\theta}),} where \eqn{\hat{\theta}} denotes the MLEs as above, and \eqn{\hat{b}_i}
    denotes the empirical Bayes estimates.
}

\note{
Predicted probabilities are not computed for joint models with \code{method = "ch-Laplace"} and \code{method = "Cox-PH-GH"}.
}

\value{
  A list of class \code{survfitJM} with components:
  \item{summaries}{a list with elements numeric matrices with numeric summaries of the predicted probabilities for each subject.}
  \item{survTimes}{a copy of the \code{survTimes} argument.}
  \item{last.time}{a numeric vector with the time of the last available longitudinal measurement of each subject.}
  \item{obs.times}{a list with elements numeric vectors denoting the timings of the longitudinal measurements for each subject.}
  \item{y}{a list with elements numeric vectors denoting the longitudinal responses for each subject.}
  \item{full.results}{a list with elements numeric matrices with predicted probabilities for each subject in each replication of the Monte Carlo
    scheme described above.}
  \item{success.rate}{a numeric vector with the success rates of the Metropolis-Hastings algorithm described above for each subject.}
  \item{scale}{a copy of the \code{scale} argument.}
}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for longitudinal and time-to-event data. 
\emph{Biometrics} \bold{67}, 819--829.

Rizopoulos, D. (2010) JM: An R Package for the Joint Modelling of Longitudinal and Time-to-Event Data. 
\emph{Journal of Statistical Software} \bold{35} (9), 1--33. \url{http://www.jstatsoft.org/v35/i09/}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{jointModel}}, \code{\link{plot.survfitJM}}}

\examples{
# linear mixed model fit
fitLME <- lme(sqrt(CD4) ~ obstime + obstime:drug, 
    random = ~ 1 | patient, data = aids)
# cox model fit
fitCOX <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# joint model fit
fitJOINT <- jointModel(fitLME, fitCOX, 
    timeVar = "obstime", method = "weibull-PH-aGH")

# sample of the patients who are still alive
ND <- aids[aids$patient == "141", ]
ss <- survfitJM(fitJOINT, newdata = ND, idVar = "patient", M = 50)
ss
}

\keyword{methods}
