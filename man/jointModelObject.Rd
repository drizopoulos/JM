\name{jointModelObject}

\alias{jointModelObject}

\title{Fitted jointModel Object}

\description{
An object returned by the \code{jointModel} function, inheriting from class \code{jointModel} and representing a fitted
joint model for longitudinal and time-to-event data. Objects of this class have methods for the generic functions 
\code{anova}, \code{coef}, \code{fitted}, \code{fixed.effects}, \code{logLik}, \code{plot}, \code{print}, 
\code{random.effects}, \code{residuals}, \code{summary}, and \code{vcov}.
}

\value{
  The following components must be included in a legitimate \code{jointModel} object.
  \item{coefficients}{a list with the estimated coefficients. The components of this list are:
    \describe{
        \item{betas}{the vector of fixed effects for the linear mixed effects model.}
        \item{sigma}{the measurement error standard deviation for the linear mixed effects model.}
        \item{gammas}{the vector of baseline covariates for the survival model.}
        \item{alpha}{the association parameter(s).}
        \item{Dalpha}{the association parameter(s) corresponding to the slope of the true trajectory.}
        \item{sigma.t}{the scale parameter for the Weibull survival model; returned only when
            \code{method = "weibull-PH-GH"} or \code{method = "weibull-AFT-GH"}.}
        \item{xi}{the parameter of the piecewise constant baseline hazard; returned only when
            \code{method = "piecewise-PH-GH"}.}
        \item{gamma.bs}{the coefficients of the B-splines use to approximate the baseline hazard;
            returned only when \code{method = "spline-PH-GH"}.}
        \item{lambda0}{a two-column numeric matrix with the first column containing the estimated baseline hazard values,
            and the second the unique sorted event times; returned only when \code{method = "Cox-PH-GH"}.}
        \item{D}{the variance-covariance matrix of the random effects.}
    }
  }
  \item{Hessian}{the Hessian matrix evaluated at the estimated parameter values.}
  \item{logLik}{the log-likelihood value.}
  \item{EB}{a list with components:
    \describe{
        \item{post.b}{the estimated random effects values.}
        \item{post.vb}{the estimated variance for the random effects estimates.}
        \item{Zb}{the estimated random effects part of the linear predictor for the longitudinal outcome (i.e., \eqn{Z}
            is the design matrix for the random effects \eqn{b}).}
        \item{Ztimeb}{the estimated random effects part of the linear predictor for the survival outcome (i.e.,
            evaluated at the observed event times).}
        \item{Ztime2b}{the estimated random effects part of the linear predictor for the survival outcome (i.e.,
            for the \eqn{i}th sample unit is evaluated at all event times that are less or equal to the \eqn{i}th
            observed event time); returned only when \code{method = "Cox-PH-GH"}.}
    }
  }
  \item{knots}{the numeric vector of the knots positions; returned only when \code{method = "spline-PH-GH"},
    \code{method = "piecewise-PH-GH"} or \code{method = "ch-Laplace"}.}
  \item{iters}{the number of iterations in the optimization algorithm.}
  \item{convergence}{convergence identifier: 0 corresponds to successful convergence, whereas 1 to a problem
    (i.e., when 1, usually more iterations are required).}
  \item{n}{the number of sample units.}
  \item{N}{the total number of repeated measurements for the longitudinal outcome.}
  \item{ni}{a vector with the number of repeated measurements for each sample unit.}
  \item{d}{a numeric vector with 0 denoting censored observation and 1 events.}
  \item{id}{the grouping vector for the longitudinal responses.}
  \item{x}{a list with the design matrices for the longitudinal and event processes.}
  \item{y}{a list with the response vectors for the longitudinal and event processes.}
  \item{data.id}{a \code{data.frame} containing the variables for the linear mixed effects model at the time of the
    event.}
  \item{method}{the value of the \code{method} argument.}
  \item{termsY}{the \code{terms} component of the \code{lmeObject}.}
  \item{termsT}{the \code{terms} component of the \code{survObject}.}
  \item{formYx}{the formula for the fixed effects part of the longitudinal model.}
  \item{formYz}{the formula for the random effects part of the longitudinal model.}
  \item{formT}{the formula for the survival model.}
  \item{timeVar}{the value of the \code{timeVar} argument}
  \item{control}{the value of the \code{control} argument.}
  \item{parameterization}{the value of the \code{parameterization} argument.}
  \item{interFact}{the value of the \code{interFact} argument}
  \item{derivForm}{the value of the \code{derivForm} argument.}
  \item{lag}{the value of the \code{lag} argument.}
  \item{call}{the matched call.}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{
    \code{\link{jointModel}}
}

\keyword{multivariate}
\keyword{regression}
