\name{weibull.frailty}

\alias{weibull.frailty}

\title{ Weibull Model with Gamma Frailties }

\description{
  Fits a Weibull model with Gamma frailties for multivariate survival data under maximum likelihood
}

\usage{
weibull.frailty(formula = formula(data), data = parent.frame(), 
    id = "id", subset, na.action, init, control = list())
} 

\arguments{
  \item{formula}{ an object of class \code{formula}: a symbolic description of the model to be fitted. The response must 
    be a survival object as returned by function \code{Surv()}.}
  \item{data}{ an optional data frame containing the variables specified in the model. }
  \item{id}{ either a character string denoting a variable name in \code{data} or a numeric vector specifying which event times belong to
    the same cluster (e.g., hospital, patient, etc.).}
  \item{subset}{ an optional vector specifying a subset of observations to be used in the fitting process. }
  \item{na.action}{ what to do with missing values. }
  \item{init}{ a numeric vector of length \eqn{p + 3} of initial values. The first \eqn{p} elements should correspond to the regression coefficients 
    for the covariates, and the last \eqn{3} to log-scale, log-shape, and log-frailty-variance, respectively. See \bold{Details}.}
  \item{control}{a list of control values with components:
    \describe{
        \item{optimizer}{a character string indicating which optimizer to use; options are "optim" (default) and
            "nlminb".}
        \item{parscale}{the \code{parscale} control argument for \code{optim()}, or the \code{scale} argument for 
            \code{nlminb()}. It should be a numeric vector of length equal to the number of parameters. Default is 0.01
            for all parameters.}
        \item{maxit}{the maximum number of iterations. Default is 500.}
        \item{numeriDeriv}{a character string indicating which type of numerical derivative to use to compute the
            Hessian matrix; options are "fd" denoting the forward difference approximation, and "cd" (default)
            denoting the central difference approximation.}
        \item{eps.Hes}{tolerance value used in the numerical derivative method. Default is 1e-03.}
    }
  }
}

\details{
  The fitted model is defined as follows: \deqn{\lambda(t_i | \omega_i) = \lambda_0(t_i) \omega_i \exp(x_i^T \beta),}{
  \lambda(t_i | \omega_i) = \lambda_0(t_i) \omega_i \exp(x_i^T \beta),} where \eqn{i} denotes the subject, \eqn{\lambda(\cdot)}{\lambda(.)} 
  denotes the hazard function, conditionally on the frailty \eqn{\omega_i}, \eqn{x_i} is a vector of covariates with corresponding regression 
  coefficients \eqn{\beta}, and \eqn{\lambda_0(\cdot)}{\lambda_0(.)} is the Weibull baseline hazard defined as \eqn{\lambda_0(t) = shape * 
  scale * t^{shape -1}}. Finally, for the frailties we assume \eqn{\omega_i \sim Gamma(\eta, \eta)}{\omega_i ~ Gamma(\eta, \eta)}, with 
  \eqn{\eta^{-1}} denoting the unknown variance of \eqn{\omega_i}'s.
}

\note{
\code{weibull.frailty()} currently supports only right-censored data.
}

\value{
  an object of class \code{weibull.frailty} with components:
  \item{coefficients}{a list with the estimated coefficients values. The components of this list are: \code{betas}, \code{scale}, \code{shape}, 
    and \code{var.frailty}, and correspond to the coefficients with the same name.}
  \item{hessian}{the hessian matrix at convergence. For the shape, scale, and var-frailty parameters the Hessian is computed on the log scale.}
  \item{logLik}{the log-likelihood value.}
  \item{control}{a copy of the \code{control} argument.}
  \item{y}{an object of class \code{Surv} containing the observed event times and the censoring indicator.}
  \item{x}{the design matrix of the model.}
  \item{id}{a numeric vector specifying which event times belong to the same cluster.}
  \item{nam.id}{the value of argument \code{id}, if that was a character string.}
  \item{terms}{the term component of the fitted model.}
  \item{data}{a copy of \code{data} or the created \code{model.frame}.}
  \item{call}{the matched call.}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\examples{
weibull.frailty(Surv(time, status) ~ age + sex, kidney)
}

\keyword{multivariate}
\keyword{regression}

