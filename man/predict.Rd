\name{predict}

\alias{predict.jointModel}

\title{Predictions for Joint Models}

\description{
  Calculates predicted values for the longitudinal part of a joint model.
}

\usage{
\method{predict}{jointModel}(object, newdata, type = c("Marginal", "Subject"),
    interval = c("none", "confidence", "prediction"), level = 0.95, idVar = "id", 
    FtTimes = NULL, M = 300, returnData = FALSE, scale = 1.6, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{newdata}{a data frame in which to look for variables with which to predict.}
  \item{type}{a character string indicating the type of predictions to compute,
        marginal or subject-specific. See \bold{Details}.}
  \item{interval}{a character string indicating what type of intervals should be computed.}
  \item{level}{a numeric scalar denoting the tolerance/confidence level.}
  \item{idVar}{a character string indicating the name of the variable in 
        \code{newdata} that corresponds to the subject identifier; required 
        when \code{type = "Subject"}.}
  \item{FtTimes}{a list with components numeric vectors denoting the time points 
        for which we wish to compute subject-specific predictions after the last
        available measurement provided in \code{newdata}. For each subject in
        \code{newdata} the default is a sequence of 25 equally spaced time points 
        from the last available measurement to the maximum follow-up time of all 
        subjects (plus a small quantity). This argument is only used when 
        \code{type = "Subject"}.}
  \item{M}{numeric scalar denoting the number of Monte Carlo samples. 
        See \bold{Details}.}
  \item{returnData}{logical; if \code{TRUE} the data frame supplied in
        \code{newdata} is returned augmented with the outputs of the function.}
  \item{scale}{a numeric value setting the scaling of the covariance matrix
        of the empirical Bayes estimates in the Metropolis step during the 
        Monte Carlo sampling.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  When \code{type = "Marginal"}, this function computes predicted values for the 
  fixed-effects part of the longitudinal submodel. In particular,
  let \eqn{X} denote the fixed-effects design matrix calculated using
  \code{newdata}. The \code{predict()} calculates \eqn{\hat{y} = X \hat{\beta}},
  and if \code{interval = "confidence"}, \eqn{var(\hat{y}) = X V X^t}, with \eqn{V}
  denoting the covariance matrix of \eqn{\hat{\beta}}. The confidence interval
  is based on a normal approximation. Confidence intervals are constructed under 
  the normal approximation.
  
  When \code{type = "Subject"}, this functions computes subject-specific 
  predictions for the longitudinal outcome based on the joint model. 
  This accomplished with a Monte Carlo simulation scheme, similar to the one
  described in \code{\link{survfitJM}}. The only difference is in Step 3, where
  for \code{interval = "confidence"} \eqn{y_i^* = X_i \beta^* + Z_i b_i^*}, whereas
  for \code{interval = "prediction"} \eqn{y_i^*} is a random vector from a normal 
  distribution with mean \eqn{X_i \beta^* + Z_i b_i^*} and standard deviation 
  \eqn{\sigma^*}. Based on this Monte Carlo simulation scheme we take as 
  estimate of \eqn{\hat{y}_i} the average of the \code{M} estimates \eqn{y_i^*} 
  from each Monte Carlo sample. Confidence intervals are constructed using the 
  percentiles of \eqn{y_i^*} from the Monte Carlo samples.
}

\value{
  If \code{se.fit = FALSE} a numeric vector of predicted values, otherwise a
  list with components \code{pred} the predicted values, \code{se.fit} the
  standard error for the fitted values, and \code{low} and \code{upp} the lower
  and upper limits of the confidence interval. If \code{returnData = TRUE}, it
  returns the data frame \code{newdata} with the previously mentioned components
  added.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.
}

\seealso{\code{\link{jointModel}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * year, 
    random = ~ year | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, 
    data = pbc2.id, x = TRUE)
# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")

DF <- with(pbc2, expand.grid(drug = levels(drug),
    year = seq(min(year), max(year), len = 100)))
Ps <- predict(fitJOINT, DF, interval = "confidence", return = TRUE)
require(lattice)
xyplot(pred + low + upp ~ year | drug, data = Ps,
    type = "l", col = c(2,1,1), lty = c(1,2,2), lwd = 2,
    ylab = "Average log serum Bilirubin")


# Subject-specific predictions
ND <- pbc2[pbc2$id == 2, ]
Ps.ss <- predict(fitJOINT, ND, type = "Subject",
  interval = "confidence", return = TRUE)
require(lattice)
xyplot(pred + low + upp ~ year | id, data = Ps.ss,
    type = "l", col = c(2,1,1), lty = c(1,2,2), lwd = 2,
    ylab = "Average log serum Bilirubin")
}
}

\keyword{methods}
