\name{residuals}

\alias{residuals.jointModel}

\title{Residuals for Joint Models}

\description{
  Calculates residuals for joint models.
}

\usage{
\method{residuals}{jointModel}(object, process = c("Longitudinal", "Event"), 
  type = c("Marginal", "Subject", "stand-Marginal", 
  "stand-Subject", "Martingale", "nullMartingale", "CoxSnell", "AFT"), 
  MI = FALSE, M = 50, time.points = NULL, return.data = FALSE, 
  \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{process}{for which model (i.e., linear mixed model or survival model) to calculate residuals.}
  \item{type}{what type of residuals to calculate. See \bold{Details}.}
  \item{MI}{logical; if \code{TRUE} multiple-imputation-based residuals are calculated.}
  \item{M}{integer denoting how many imputations to use for the MI residuals.}
  \item{time.points}{for fixed visit times, this should be a numeric vector with the unique times points at which 
    longitudinal measurements are supposed to be taken; if \code{NULL}, then the code attempts to extract these unique time points
    using the design matrix for the fixed effects of the longitudinal model and the value of the \code{timeVar} argument of
    \code{\link{jointModel}}. For random visit times, this should be an object of class \code{\link{weibull.frailty}} that represents the fit of
    Weibull model with Gamma frailties for the visiting process. The user may also augment the object \code{weibull.frailty} with the following 
    two attributes: \code{"prev.y"} denoting the variable name for the previous longitudinal responses, and \code{"tmax"} denoting the end of 
    the study.}
  \item{return.data}{logical; if \code{TRUE} and \code{MI = TRUE} and fixed visit times are considered, then the multiply imputed data sets 
    are returned.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  When \code{process = "Longitudinal"}, residuals are calculated for the longitudinal outcomes. In particular, if 
  \code{type = "Marginal"} these are \eqn{e_{ij} = y_{ij} - x_{ij}^T \hat{\beta}}, whereas for \code{type = "Subject"}, 
  \eqn{e_{ij} = y_{ij} - x_{ij}^T \hat{\beta} - z_{ij}^T b_i}, where \eqn{i} denotes the subject and \eqn{j} the 
  measurement, \eqn{y_{ij}} the longitudinal responses, \eqn{x_{ij}^T} and \eqn{z_{ij}^T} the corresponding rows of the
  fixed and random effects design matrices, respectively, and \eqn{\beta} and \eqn{b_i} denote the fixed effects 
  and random effects components. If \code{type = "stand-Marginal"} or \code{type = "stand-Subject"}, the above defined 
  residuals are divided by the estimated standard deviation of the corresponding error term. If \code{MI = TRUE}, multiple-imputation-based 
  residuals are calculated for the longitudinal process; for more information regarding these residuals, check Rizopoulos et al. (2009).
  
  When \code{process = "Event"}, residuals are calculated for the survival outcome. Martingale residuals are available 
  for all options for the survival submodel (for the different options of survival submodel, check the \code{method} 
  argument of \code{\link{jointModel}}). when option \code{type = "nullMartingale"} is invoked, the martingale residuals
  are calculated with the coefficient(s) that correspond to the marker set to zero. Cox-Snell residuals (Cox and Snell, 
  1968) are available for the Weibull model and the additive log cumulative hazard model. AFT residuals are only 
  available for the Weibull model.
}

\note{
The multiple-imputation-based residuals are not available for joint models with \code{method = "Cox-PH-GH"}.
}

\value{
  If \code{MI = FALSE}, a numeric vector of residual values. Otherwise a list with components:
  \item{fitted.values}{the fitted values for the observed data.}
  \item{residuals}{the residuals for the observed data.}
  \item{fitted.valsM}{the fitted values for the missing data.}
  \item{resid.valsM}{the multiply imputed residuals for the missing longitudinal responses.}
  \item{mean.resid.valsM}{the average of the multiply imputed residuals for the missing longitudinal responses; returned only if 
    fixed visit times are considered.}
  \item{dataM}{if \code{return.data = TRUE} and fixed visit times are considered, then it returns the data set with the simulated response 
    values for the longitudinal outcome, for each of the multiple imputations.}    
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\references{
Cox, D. and Snell, E. (1968) A general definition of residuals. \emph{Journal of the Royal Statistical Society, 
Series B} \bold{30}, 248--275.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D., Verbeke, G. and Molenberghs, G. (2010) Multiple-imputation-based residuals and diagnostic plots 
for joint models of longitudinal and survival outcomes. \emph{Biometrics} \bold{66}, 20--29.

Rizopoulos, D. (2010) JM: An R Package for the Joint Modelling of Longitudinal and Time-to-Event Data. 
\emph{Journal of Statistical Software} \bold{35} (9), 1--33. \url{http://www.jstatsoft.org/v35/i09/}
}

\seealso{\code{\link{fitted.jointModel}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(sqrt(CD4) ~ obstime * drug - drug, 
    random = ~ 1 | patient, data = aids)
# cox model fit
fitCOX <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# joint model fit, under the additive log cumulative hazard model
fitJOINT <- jointModel(fitLME, fitCOX, 
    timeVar = "obstime")

# residuals for the longitudinal outcome
head(cbind(
    "Marginal" = residuals(fitJOINT),
    "std-Marginal" = residuals(fitJOINT, type = "stand-Marginal"),
    "Subject" = residuals(fitJOINT, type = "Subject"),
    "std-Subject" = residuals(fitJOINT, type = "stand-Subject")
))

# residuals for the survival outcome
head(cbind(
    "Martingale" = residuals(fitJOINT, process = "Event", type = "Martingale"),
    "CoxSnell" = residuals(fitJOINT, process = "Event", type = "CoxSnell")
))
}
}

\keyword{methods}
