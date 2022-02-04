\name{rocJM}
\alias{rocJM}
\alias{rocJM.jointModel}

\title{ Predictive Accuracy Measures for Longitudinal Markers under a Joint Modelling Framework }

\description{
  It computes sensitivity, specificity, ROC and AUC measures for joint models.
}

\usage{
rocJM(object, dt, data, idVar = "id", directionSmaller = NULL, cc = NULL, min.cc = NULL,
  max.cc = NULL, optThr = c("sens*spec", "youden"), 
  diffType = c("absolute", "relative"), abs.diff = 0, rel.diff = 1, 
  M = 300, burn.in = 100, scale = 1.6)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{dt}{a numeric vector indicating the lengths of the time intervals of primary interest within which we want to 
        distinguish between subjects who died within the intervals from subjects who survived longer than that.}
  \item{data}{a data frame that contains the baseline covariates for the longitudinal and survival submodels,
        including a case identifier (i.e., the variable denoted by the argument \code{idVar}), the time points on
        which longitudinal measurements are assumed to be taken (this should have the same name as in the argument
        \code{timeVar} of \code{\link{jointModel}}).}
  \item{idVar}{the name of the variable in \code{data} that identifies the different generic subjects to be considered.}
  \item{directionSmaller}{logical; if \code{TRUE}, then smaller values for the longitudinal outcome are associated 
		with higher risk for an event.}
  \item{cc}{a numeric vector of threshold values for the longitudinal marker; if \code{NULL}, this is computed using 
        a regular sequence based on percentiles of the observed marker values.}
  \item{min.cc}{the start of the regular sequence for the threshold values for the longitudinal marker; 
        see argument \code{cc} above.}
  \item{max.cc}{the end of the regular sequence for the threshold values for the longitudinal marker; 
        see argument \code{cc} above.}
  \item{optThr}{character string defining how the optimal threshold is to be computed. The default chooses the
	cut-point for the marker that maximizes the product of sensitivity and specificity. Option \code{"youden"}
	chooses the cut-point that maximizes Youden's index that equals sensitivity + specificity - 1.}
  \item{diffType}{character string defining the type of prediction rule. See \bold{Details}.}
  \item{abs.diff}{a numeric vector of absolute differences in the definition of composite prediction rules.}
  \item{rel.diff}{a numeric vector of relative differences in the definition of composite prediction rules.}
  \item{M}{a numeric scalar denoting the number of Monte Carlo samples.}
  \item{burn.in}{a numeric scalar denoting the iterations to discard.}
  \item{scale}{a numeric scalar that controls the acceptance rate of the Metropolis-Hastings algorithm. See \bold{Details}.}
}

\details{
(\bold{Note:} the following contain some math formulas, which are better viewed in the pdf version of the manual
accessible at \url{https://cran.r-project.org/package=JM}.)

Assume that we have collected longitudinal measurements \eqn{Y_i(t) = \{y_i(s); 0 \leq s \leq t\}} up to time point \eqn{t} for
subject \eqn{i}. We are interested in events occurring in the medically relevant time frame \eqn{(t, t + \Delta t]} within which the
physician can take an action to improve the survival chance of the patient. Using an appropriate function of the marker history
\eqn{Y_i(t)}, we can define a prediction rule to discriminate between patients of high and low risk for an event.  For instance, 
for in HIV infected patients, we could consider values of CD4 cell count smaller than a specific threshold as predictive for death. 
Since we are in a longitudinal context, we have the flexibility of determining which values of the longitudinal history of the 
patient will contribute to the specification of the prediction rule. That is, we could define a composite prediction rule that is not 
based only on the last available measurement but rather on the last two or last three measurements of a patient. Furthermore, it 
could be of relevance to consider different threshold values for each of these measurements, for instance, we could define as success 
the event that the pre-last CD4 cell count is \eqn{c} and the last one \eqn{0.5c}, indicating that a 50\% decrease is strongly 
indicative for death. Under this setting we define sensitivity and specificity as,
\deqn{Pr \bigl \{ {\cal S}_i(t, k, c) \mid T_i^* > t, T_i^* \in (t, t + \Delta t] \bigr \},} 
and \deqn{Pr \bigl \{ {\cal F}_i(t, k, c) \mid T_i^* > t, T_i^* > t +
\Delta t \bigr \},} respectively, where we term \eqn{{\cal S}_i(t, k, c) = \{y_i(s) \leq c_s; k \leq s \leq t\}} as success 
(i.e., occurrence of the event of interest), and \eqn{{\cal F}_i(t, k, c) = \{y_i(s) > c_s; k \leq s \leq t\}} as a failure,
\eqn{T_i^*} denotes the time-to-event, and \eqn{\Delta t} the length of the medically relevant time window (specified by argument 
\code{dt}). The cut values for the marker \eqn{c} are specified by the \code{cc}, \code{min.cc} and \code{max.cc} arguments. Two types of
composite prediction rules can be defined depending on the value of the \code{diffType} argument. Absolute prediction rules in which, between
successive measurements there is an absolute difference of between the cut values, and relative prediction rules in which there is a 
relative difference between successive measurements of the marker. The lag values for these differences are defined by the \code{abs.diff}
and \code{rel.diff} arguments. Some illustrative examples:
\describe{
  \item{Ex1:}{keeping the defaults we define a simple rule that is only based on the last available marker measurement.}
  \item{Ex2:}{to define a prediction rule that is based on the last two available measurements using the same cut values (e.g.,
        if a patient had two successive measurements below a medically relevant threshold), we need to set \code{abs.diff = c(0, 0)}.}
  \item{Ex3:}{to define a prediction rule that is based on the last two available measurements using a drop of 5 units between the cut 
        values (e.g., the pre-last measurement is \eqn{c} and the last \eqn{c-5}), we need to set \code{abs.diff = c(0, -5)}.}
  \item{Ex4:}{to define a prediction rule that is based on the last two available measurements using a drop of 20\% units between the cut 
        values (e.g., the pre-last measurement is \eqn{c} and the last \eqn{0.8c}), we need to set \code{diffType = "relative"} and 
        \code{rel.diff = c(0, 0.8)}.}
}

The estimation of the above defined probabilities is achieved with a Monte Carlo scheme similar to the one described in
\code{\link{survfitJM}}. The number of Monte Carlo samples is defined by the \code{M} argument, and the burn-in iterations for
the Metropolis-Hastings algorithm using the \code{burn.in} argument.

More details can be found in Rizopoulos (2011).
}

\value{
  An object of class \code{rocJM} is a list with components,
  \item{MCresults}{a list of length the number of distinct cases in \code{data}. Each component of this
        list is again a list with four components the estimated Sensitivity \code{Sens} and its standard 
        error \code{seSens}, and the estimated Specificity \code{Spec} and its standard error \code{seSpec}.
        All these four components are matrices with rows corresponding to the different \code{dt} values and
        columns corresponding to the different \code{cc} values.}
  \item{AUCs}{a numeric vector of estimated areas under the ROC curves for the different values of \code{dt}.}
  \item{optThr}{a numeric vector with the optimal threshold values for the markers for the different
        \code{dt} under the choice made in argument \code{optThr}.}
  \item{times}{a list of length the number of distinct cases in \code{data} with components numeric vectors
        of the time points at which longitudinal measurements are supposed to be taken.}
  \item{dt}{a copy of the \code{dt} argument.}
  \item{M}{a copy of the \code{M} argument.}
  \item{diffType}{a copy of the \code{diffType} argument.}
  \item{abs.diff}{a copy of the \code{abs.diff} argument.}
  \item{rel.diff}{a copy of the \code{rel.diff} argument.}
  \item{cc}{a copy of the \code{cc} argument.}
  \item{min.cc}{a copy of the \code{min.cc} argument.}
  \item{max.cc}{a copy of the \code{max.cc} argument.}
  \item{success.rate}{a numeric matrix with the success rates of the Metropolis-Hastings algorithm described above.}
}

\references{
Heagerty, P. and Zheng, Y. (2005). Survival model predictive accuracy and ROC curves.
\emph{Biometrics} \bold{61}, 92--105.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for longitudinal 
and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.

Rizopoulos, D. (2010) JM: An R package for the joint modelling of longitudinal and time-to-event data. 
\emph{Journal of Statistical Software} \bold{35} (9), 1--33. \doi{10.18637/jss.v035.i09}

Zheng, Y. and Heagerty, P. (2007). Prospective accuracy for longitudinal markers.
\emph{Biometrics} \bold{63}, 332--341.
}

\author{
    Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
    \code{\link{plot.rocJM}},
    \code{\link{survfitJM}},
    \code{\link{dynCJM}},
    \code{\link{aucJM}},
    \code{\link{prederrJM}},
    \code{\link{jointModel}}
    
}

\examples{
\dontrun{
fitLME <- lme(sqrt(CD4) ~ obstime * (drug + AZT + prevOI + gender), 
    random = ~ obstime | patient, data = aids)
fitSURV <- coxph(Surv(Time, death) ~ drug + AZT + prevOI + gender, 
    data = aids.id, x = TRUE)
fit.aids <- jointModel(fitLME, fitSURV, timeVar = "obstime", 
    method = "piecewise-PH-aGH")

# the following will take some time to execute...
ND <- aids[aids$patient == "7", ]
roc <- rocJM(fit.aids, dt = c(2, 4, 8), ND, idVar = "patient")
roc
}
}

\keyword{methods}
