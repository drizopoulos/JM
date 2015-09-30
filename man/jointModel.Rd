\name{jointModel}

\alias{jointModel}

\title{Joint Models for Longitudinal and Survival Data}

\description{
This function fits shared parameter models for the joint modelling of normal longitudinal responses and time-to-event
data under a maximum likelihood approach. Various options for the survival model are available.
}

\usage{
jointModel(lmeObject, survObject, timeVar, 
  parameterization = c("value", "slope", "both"),
  method = c("weibull-PH-aGH", "weibull-PH-GH", "weibull-AFT-aGH", 
    "weibull-AFT-GH", "piecewise-PH-aGH", "piecewise-PH-GH", 
    "Cox-PH-aGH", "Cox-PH-GH", "spline-PH-aGH", "spline-PH-GH", 
    "ch-Laplace"),
  interFact = NULL, derivForm = NULL, lag = 0, scaleWB = NULL,
  CompRisk = FALSE, init = NULL, control = list(), \dots)
}

\arguments{
  \item{lmeObject}{an object inheriting from class \code{lme} (see also \bold{Note}).}
  \item{survObject}{an object inheriting from class \code{coxph} or class \code{survreg}. In the call to \code{coxph()}
    or \code{survreg()}, you need to specify the argument \code{x = TRUE} such that the design matrix is contained in 
    the object fit. See \bold{Examples}.}
  \item{timeVar}{a character string indicating the time variable in the linear mixed effects model.}
  \item{parameterization}{a character string indicating the type of parameterization. See \bold{Details}}
  \item{method}{a character string specifying the type of joint model to fit. See \bold{Details}.}
  \item{interFact}{a list with components \code{value} a formula for the interaction terms corresponding to the
    \code{value} parameterization, \code{slope} a formula for the interaction terms corresponding to the
    \code{slope} parameterization, \code{data} a data frame containing these variables (this should have the same
    number of rows and ordering of subjects, as the one in \code{survObject}).}
  \item{derivForm}{a list with components \code{fixed} a formula representing the derivative of the fixed-effects part of the 
    liner mixed model with respect to time, \code{indFixed} a numeric vector indicating which fixed effects of \code{lmeObject} 
    correspond to the derivative, \code{random} a formula representing the derivative of the random-effects part of the 
    liner mixed model with respect to time, and \code{indRamdom} a numeric vector indicating which random effects of \code{lmeObject} 
    correspond to the derivative. When a random intercepts linear mixed model is assumed, then \code{random = ~ 1} and 
    \code{indRandom = FALSE}. Required only when \code{parameterization == "slope"} or \code{parameterization == "both"}. 
    See \bold{Examples}.}
  \item{lag}{a numeric scalar denoting a lag effect in the time-dependent covariate represented by the mixed model; default is 0.}
  \item{scaleWB}{a numeric scalar denoting a fixed value for the scale parameter of the Weibull hazard; used only when 
    \code{method = "weibull-AFT-GH"} or \code{method = "weibull-PH-GH"}. The default \code{NULL} means that the scale 
    parameter is estimated.}
  \item{CompRisk}{logical; should a competing risks joint model be fitted.}
  \item{init}{a named list of user-specified initial values:
    \describe{
        \item{betas}{the vector of fixed effects for the linear mixed effects model.}
        \item{sigma}{the measurement error standard deviation for the linear mixed effects model.}
        \item{D}{the variance-covariance matrix of the random effects.}
        \item{gammas}{the vector of baseline covariates for the survival model. For \code{method = "ch-Laplace"} this vector should 
            first contain initial values for the sorted B-spline coefficients used to model the log cumulative baseline hazard.}
        \item{alpha}{the association parameters.}
        \item{Dalpha}{the association parameters for the true slopes parameterization}
        \item{xi}{the vector of baseline risk function values within the intervals specified by the knots; specified only 
            when \code{method = "piecewise-PH-GH"}.}
        \item{gammas.bs}{the vector of spline coefficients; specified only when \code{method = "spline-PH-GH"}.}
        \item{sigma.t}{the scale parameter for the Weibull baseline risk function; specified only when 
            \code{method = "weibull-AFT-GH"} or \code{method = "weibull-PH-GH"}.}
        \item{lambda0}{a vector of the baseline hazard values at the sorted unique event times; specified only when 
            \code{method = "Cox-PH-GH"}.}
    }
    When this list of initial values does not contain some of these components or contains components 
    not of the appropriate length, then the default initial values are used instead.
  }
  \item{control}{a list of control values with components:
    \describe{
        \item{only.EM}{logical; if \code{TRUE} only the EM algorithm is used in the optimization, otherwise if 
            convergence has not been achieved a quasi-Newton algorithm is initiated. Default is \code{FALSE} except for 
            \code{method = "Cox-PH-GH"} for which only the EM algorithm is available.}
        \item{iter.EM}{the number of EM iterations. Default is 50 except for \code{method = "Cox-PH-GH"} for which
            the default is 200.}
        \item{iter.qN}{the number of quasi-Newton iterations. Default is 150.}
        \item{optimizer}{a character string indicating which optimizer to use; options are "optim" (default) and
            "nlminb".}
        \item{tol1}{tolerance value for convergence in the parameters; see \bold{Details}. Default is 1e-03.}
        \item{tol2}{tolerance value for convergence in the parameters; see \bold{Details}. Default is 1e-04.}
        \item{tol3}{tolerance value for convergence in the log-likelihood; see \bold{Details}. Default is 
            \code{sqrt(.Machine$double.eps)}.}
        \item{numeriDeriv}{a character string indicating which type of numerical derivative to use to compute the
            Hessian matrix; options are "fd" (default) denoting the forward difference approximation, and "cd" 
            denoting the central difference approximation.}
        \item{eps.Hes}{tolerance value used in the numerical derivative method. Default is 1e-06; if you choose 
            \code{numeriDeriv = "cd"} a larger value (e.g., 1e-04) is suggested.}
        \item{parscale}{the \code{parscale} control argument for \code{optim()}, or the \code{scale} argument for 
            \code{nlminb()}. It should be a numeric vector of length equal to the number of parameters. Default is 0.01
            for all parameters.}
        \item{step.max}{tolerance value for the maximum step size in the Newton-Raphson algorithm used to update the 
            parameters of the survival submodel for \code{method = "ch-Laplace"}. Default is 0.1.}
        \item{backtrackSteps}{the number of backtrack steps to use when updating the parameters of the survival submodel
            under \code{method = "ch-Laplace"}.}
        \item{knots}{a numeric vector of the knots positions for the piecewise constant baseline risk function of for 
            the log times used in the B-splines approximation of the log cumulative baseline hazard; therefore, this argument 
            is relevant only when \code{method = "piecewise-PH-GH"}, \code{method = "spline-PH-GH"} or \code{method = "ch-Laplace"}. 
            The default is to place equally-spaced \code{lng.in.kn} knots in the quantiles of the observed event times. For stratified models
            fitted with \code{method = "spline-PH-GH"} this should be a list with elements numeric vectors of knots positions for each strata.}
	\item{ObsTimes.knots}{logical; if \code{TRUE} (default) the positions of the knots are specified based in the observed event times, 
	     otherwise the positions of the knots are specified using only the true event times.}
	\item{lng.in.kn}{the number of internal knots; relevant only when when \code{method = "piecewise-PH-GH"} where it 
            denotes the number of internal knots for the piecewise constant baseline risk function or when \code{method = "spline-PH-GH"} 
            or \code{method = "ch-Laplace"} where it denotes the number of internal knots for B-splines approximation of the log 
            baseline hazard. Default is 6 when \code{method = "piecewise-PH-GH"} and 5 otherwise.}
 	\item{equal.strata.knots}{logical; if \code{TRUE} (the default), then the same knots are used in the approximation of the 
            baseline risk function in different strata when \code{method = "spline-PH-GH"}.}
        \item{ord}{a positive integer denoting the order of the B-splines used to approximate the log cumulative hazard
            (default is 4); relevant only when \code{method = "spline-PH-GH"} or \code{method = "ch-Laplace"}.}
        \item{typeGH}{a character string indicating the type of Gauss-Hermite rule to be used. Options are "simple"
            and "adaptive". The default is "simple" but it is turned to adaptive when the user specifies in the
            \code{method} argument an option that contains \code{aGH}.}
        \item{GHk}{the number of Gauss-Hermite quadrature points used to approximate the integrals over the random 
            effects. The default is 15 for one- or two-dimensional integration and for \eqn{N < 2000}, and 9 otherwise for the 
            simple Gauss-Hermite rule, and 5 for one-, two-dimensional or three-dimensional integration and for \eqn{N < 2000}, 
            and 3 otherwise for the pseudo adaptive Gauss-Hermite rule, where \eqn{N} denotes the total number of longitudinal
            measurements.}
        \item{GKk}{the number of Gauss-Kronrod points used to approximate the integral involved in the calculation of the survival function. Two 
            options are available, namely 7 or 15. For \code{method = "weibull-PH-GH"}, \code{method = "weibull-AFT-GH"} and
            \code{method = "spline-PH-GH"} 15 are used, whereas for \code{method = "piecewise-PH-GH"} 7.}
        \item{verbose}{logical; if \code{TRUE}, the parameter estimates and the log-likelihood value are printed during
            the optimization procedure. Default is \code{FALSE}.}
    }
  }
  \item{\dots}{options passed to the \code{control} argument.}
}

\details{
Function \code{jointModel} fits joint models for longitudinal and survival data (more detailed information about the formulation of these
models can be found in Rizopoulos (2010)). For the longitudinal responses the linear mixed effects model represented by the \code{lmeObject} is 
assumed. For the survival times let \eqn{w_i} denote the vector of baseline covariates in \code{survObject}, with associated parameter vector 
\eqn{\gamma}, \eqn{m_i(t)} the value of the longitudinal outcome at time point \eqn{t} as approximated by the linear mixed model 
(i.e., \eqn{m_i(t)} equals the fixed-effects part \code{+} random-effects part of the linear mixed effects model for sample unit \eqn{i}), 
\eqn{\alpha} the association parameter for \eqn{m_i(t)}, \eqn{m_i'(t)} the derivative of \eqn{m_i(t)} with respect to \eqn{t}, and 
\eqn{\alpha_d} the association parameter for \eqn{m_i'(t)}. Then, for \code{method = "weibull-AFT-GH"} a time-dependent Weibull model under 
the accelerated failure time formulation is assumed. For \code{method = "weibull-PH-GH"} a time-dependent relative risk model is postulated 
with a Weibull baseline risk function. For \code{method = "piecewise-PH-GH"} a time-dependent relative risk model is postulated with a 
piecewise constant baseline risk function. For \code{method = "spline-PH-GH"} a time-dependent relative risk model is assumed in which the 
log baseline risk function is approximated using B-splines. For \code{method = "ch-Laplace"} an additive model on the log cumulative hazard 
scale is assumed (see Rizopoulos et al., 2009 for more info). Finally, for \code{method = "Cox-PH-GH"} a time-dependent relative risk model 
is assumed where the baseline risk function is left unspecified (Wulfsohn and Tsiatis, 1997). For all these options the linear predictor for the
survival submodel is written as \deqn{\eta = \gamma^\top w_i + \alpha m_i\{max(t-k, 0)\},}{eta = gamma' w_i + alpha m_i{max(t-k,0)},} when 
\code{parameterization = "value"}, \deqn{\eta = \gamma^\top w_i + \alpha_s m_i'\{max(t-k, 0)\},}{eta = gamma' w_i + alpha m_i'{max(t-k, 0)}} 
when \code{parameterization = "slope"}, and \deqn{\eta = \gamma^\top w_i + \alpha m_i\{max(t-k, 0)\} + \alpha_s m_i'\{max(t-k, 0)\},}{eta = 
gamma' w_i + alpha m_i{max(t-k, 0)} + \alpha_s m_i'{max(t-k, 0)},} when \code{parameterization = "both"}, where in all the above the value 
of \eqn{k} is specified by the \code{lag} argument and \eqn{m_i'(t) = d m_i(t) / dt}. If \code{interFact} is specified, then
\eqn{m_i\{max(t-k, 0)\}} and/or \eqn{m_i'\{max(t-k, 0)\}} are multiplied with the design matrices derived from the formulas
supplied as the first two arguments of \code{interFact}, respectively. In this case \eqn{\alpha} and/or \eqn{\alpha_s} become vectors of
association parameters.

For \code{method = "spline-PH-GH"} it is also allowed to include stratification factors. These should be included in the specification of
the \code{survObject} using function \code{strata()}. Note that in this case \code{survObject} must only be a 'coxph' object.

For all survival models except for the time-dependent proportional hazards model, the optimization algorithm starts 
with EM iterations, and if convergence is not achieved, it switches to quasi-Newton iterations (i.e., BFGS in 
\code{optim()} or \code{nlminb()}, depending on the value of the \code{optimizer} control argument). For \code{method = "Cox-PH-GH"} only the 
EM algorithm is used. During the EM iterations, convergence is declared if either of the following two conditions is satisfied: (i) 
\eqn{ L(\theta^{it}) - L(\theta^{it - 1}) < tol_3 \{ | L(\theta^{it - 1}) | + tol_3 \} }, or (ii) 
\eqn{ \max \{ | \theta^{it} - \theta^{it - 1} | / ( | \theta^{it - 1} | + tol_1) \} < tol_2}, where \eqn{\theta^{it}} and 
\eqn{\theta^{it - 1}} is the vector of parameter values at the current and previous iterations, respectively, and \eqn{L(.)} is the 
log-likelihood function. The values for \eqn{tol_1}, \eqn{tol_2} and \eqn{tol_3} are specified via the \code{control} argument. During the 
quasi-Newton iterations, the default convergence criteria of either \code{optim()} or \code{nlminb()} are used.

The required integrals are approximated using the standard Gauss-Hermite quadrature rule when the chosen option for the \code{method}
argument contains the string "GH", and the (pseudo) adaptive Gauss-Hermite rule when the chosen option for the \code{method}
argument contains the string "aGH". For \code{method = "ch-Laplace"} the fully exponential Laplace approximation described in 
Rizopoulos et al. (2009) is used.  The (pseudo) adaptive Gauss-Hermite and the Laplace approximation are particularly useful when 
high-dimensional random effects vectors are considered (e.g., when modelling nonlinear subject-specific trajectories with splines 
or high-order polynomials).
}

\value{
See \code{\link{jointModelObject}} for the components of the fit.
}

\references{
Henderson, R., Diggle, P. and Dobson, A. (2000) Joint modelling of longitudinal measurements and event time data.
\emph{Biostatistics} \bold{1}, 465--480.

Hsieh, F., Tseng, Y.-K. and Wang, J.-L. (2006) Joint modeling of survival and longitudinal data: Likelihood 
approach revisited. \emph{Biometrics} \bold{62}, 1037--1043.

Rizopoulos, D. (2012a) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2012b) Fast fitting of joint models for longitudinal and event time data using a 
pseudo-adaptive Gaussian quadrature rule. \emph{Computational Statistics and Data Analysis} \bold{56}, 491--501.

Rizopoulos, D. (2011) Dynamic predictions and prospective accuracy in joint models for longitudinal 
and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.

Rizopoulos, D. (2010) JM: An R package for the joint modelling of longitudinal and time-to-event data. 
\emph{Journal of Statistical Software} \bold{35} (9), 1--33. \url{http://www.jstatsoft.org/v35/i09/}

Rizopoulos, D., Verbeke, G. and Lesaffre, E. (2009) Fully exponential Laplace approximations for the joint
modelling of survival and longitudinal data. \emph{Journal of the Royal Statistical Society, Series B} \bold{71}, 
637--654.

Rizopoulos, D., Verbeke, G. and Molenberghs, G. (2010) Multiple-imputation-based residuals and diagnostic plots
for joint models of longitudinal and survival outcomes. \emph{Biometrics} \bold{66}, 20--29.

Tsiatis, A. and Davidian, M. (2004) Joint modeling of longitudinal and time-to-event data: an overview. 
\emph{Statistica Sinica} \bold{14}, 809--834.

Wulfsohn, M. and Tsiatis, A. (1997) A joint model for survival and longitudinal data measured with error. 
\emph{Biometrics} \bold{53}, 330--339.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\note{
1. The \code{lmeObject} argument should represent a linear mixed model object with a simple random-effects 
structure, i.e., only the \code{pdDiag()} class is currently allowed.

2. The \code{lmeObject} object should not contain any within-group correlation structure (i.e., \code{correlation} 
argument of \code{lme()}) or within-group heteroscedasticity structure (i.e., \code{weights} argument of \code{lme()}).

3. It is assumed that the linear mixed effects model \code{lmeObject} and the survival model \code{survObject} have been 
fitted to the same subjects. Moreover, it is assumed that the ordering of the subjects is the same for both 
\code{lmeObject} and \code{survObject}, i.e., that the first line in the data frame containing the event times 
corresponds to the first set of lines identified by the grouping variable in the data frame containing the repeated 
measurements, and so on.

4. In the \code{print} and \code{summary} generic functions for class \code{jointModel}, the estimated coefficients (and
standard errors for the \code{summary} generic) for the event process are augmented with the element "Assoct" that 
corresponds to the association parameter \eqn{\alpha} and the element "Assoct.s" that corresponds to the parameter 
\eqn{\alpha_s} when \code{parameterization} is \code{"slope"} or \code{"both"} (see \bold{Details}).

5. The standard errors returned by the \code{summary} generic function for class \code{jointModel} when 
\code{method = "Cox-PH-GH"} are based on the profile score vector (i.e., given the NPMLE for the unspecified baseline 
hazard). Hsieh et al. (2006) have noted that these standard errors are underestimated.

6. As it is the case for all types of mixed models that require numerical integration, it is advisable (especially in
difficult datasets) to check the stability of the maximum likelihood estimates with an increasing number of 
Gauss-Hermite quadrature points.

7. It is assumed that the scale of the time variable (e.g., days, months years) is the same in both \code{lmeObject} and \code{survObject}.
}

\seealso{
    \code{\link{jointModelObject}},
    \code{\link{anova.jointModel}},
    \code{\link{coef.jointModel}},
    \code{\link{fixef.jointModel}},
    \code{\link{ranef.jointModel}},
    \code{\link{fitted.jointModel}},
    \code{\link{residuals.jointModel}},
    \code{\link{plot.jointModel}},
    \code{\link{survfitJM}},
    \code{\link{rocJM}},
    \code{\link{dynCJM}},
    \code{\link{aucJM}},
    \code{\link{prederrJM}}
}

\examples{
\dontrun{
# linear mixed model fit (random intercepts)
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ 1 | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)
# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")
fitJOINT
summary(fitJOINT)

# linear mixed model fit (random intercepts + random slopes)
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ year | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)
# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")
fitJOINT
summary(fitJOINT)

# we also include an interaction term of log(serBilir) with drug
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year",
    interFact = list(value = ~ drug, data = pbc2.id))
fitJOINT
summary(fitJOINT)


# a joint model in which the risk for and event depends both on the true value of
# marker and the true value of the slope of the longitudinal trajectory
lmeFit <- lme(sqrt(CD4) ~ obstime * drug, random = ~ obstime | patient, data = aids)
coxFit <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# to fit this model we need to specify the 'derivForm' argument, which is a list
# with first component the derivative of the fixed-effects formula of 'lmeFit' with
# respect to 'obstime', second component the indicator of which fixed-effects 
# coefficients correspond to the previous defined formula, third component the 
# derivative of the random-effects formula of 'lmeFit' with respect to 'obstime', 
# and fourth component the indicator of which random-effects correspond to the 
# previous defined formula
dForm <- list(fixed = ~ 1 + drug, indFixed = c(2, 4), random = ~ 1, indRandom = 2)
jointModel(lmeFit, coxFit, timeVar = "obstime", method = "spline-PH-aGH",
  parameterization = "both", derivForm = dForm)


# Competing Risks joint model
# we first expand the PBC dataset in the competing risks long format
# with two competing risks being death and transplantation
pbc2.idCR <- crLong(pbc2.id, "status", "alive")

# we fit the linear mixed model as before
lmeFit.pbc <- lme(log(serBilir) ~ drug * ns(year, 3), 
    random = list(id = pdDiag(form = ~ ns(year, 3))), data = pbc2)

# however, for the survival model we need to use the data in the long
# format, and include the competing risks indicator as a stratification
# factor. We also take interactions of the baseline covariates with the
# stratification factor in order to allow the effect of these covariates
# to be different for each competing risk
coxCRFit.pbc <- coxph(Surv(years, status2) ~ (drug + sex)*strata + strata(strata), 
    data = pbc2.idCR, x = TRUE)

# the corresponding joint model is fitted simply by including the above
# two submodels as main arguments, setting argument CompRisk to TRUE, 
# and choosing as method = "spline-PH-aGH". Similarly as above, we also 
# include strata as an interaction factor to allow serum bilirubin to 
# have a different effect for each of the two competing risks
jmCRFit.pbc <- jointModel(lmeFit.pbc, coxCRFit.pbc, timeVar = "year", 
    method = "spline-PH-aGH", 
    interFact = list(value = ~ strata, data = pbc2.idCR), 
    CompRisk = TRUE)
summary(jmCRFit.pbc)

# linear mixed model fit
fitLME <- lme(sqrt(CD4) ~ obstime * drug - drug, 
    random = ~ 1 | patient, data = aids)
# cox model fit
fitCOX <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)
# joint model fit with a spline-approximated baseline hazard function
fitJOINT <- jointModel(fitLME, fitCOX, 
    timeVar = "obstime", method = "spline-PH-aGH")
fitJOINT
summary(fitJOINT)
}
}

\keyword{multivariate}
\keyword{regression}
