\name{plot.rocJM}

\alias{plot.rocJM}

\title{Plot Method for rocJM Objects}

\description{
  Produces plots of ROC curves and the corresponding areas under the curve.
}

\usage{
\method{plot}{rocJM}(x, which = NULL, type = c("ROC", "AUC"), 
  ndt = "all", main = NULL, caption = NULL, xlab = NULL, 
  ylab = NULL, ask = NULL, legend = FALSE, lx = NULL, ly = NULL, 
  lty = NULL, col = NULL, cex.caption = 0.8, cex.axis = NULL, 
  cex.lab = NULL, cex.main = NULL, \dots)
}

\arguments{
  \item{x}{an object inheriting from class \code{rocJM}.}
  \item{which}{a numeric vector specifying for which generic subjects to produce the plots. 
    This refers to the different cases identified by the \code{idVar} argument in \code{\link{rocJM}}.}
  \item{type}{a character string specifying which plot to produce the ROC curves or the areas under
    the ROC curves.}
  \item{ndt}{the character string \code{"all"} or a numeric scalar specifying for which time windows 
    (\code{dt} argument of \code{\link{rocJM}}) to produce the plots.}
  \item{main}{a character string specifying the title in the plot.}
  \item{caption}{a character string specifying a caption in the plot.}
  \item{xlab}{a character string specifying the x-axis label in the plot.}
  \item{ylab}{a character string specifying the y-axis label in the plot.}
  \item{ask}{logical; if \code{TRUE}, the user is asked before each plot, see \code{par()}.}
  \item{legend}{logical; if \code{TRUE}, a legend is included in the plot.}
  \item{lx,ly}{the \code{x} and \code{y} arguments of \code{legend()}.}
  \item{lty}{what types of lines to use.}
  \item{col}{which colors to use.}
  \item{cex.caption}{font size for the caption.}
  \item{cex.axis, cex.lab, cex.main}{graphical parameters; see \code{par} for more info.}
  \item{\dots}{extra graphical parameters passed to \code{plot()}.}
}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for 
longitudinal and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.
}


\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{rocJM}}}

\examples{
\dontrun{
fitLME <- lme(sqrt(CD4) ~ obstime + obstime:(drug + AZT + prevOI + gender), 
    random = ~ obstime | patient, data = aids)
fitSURV <- coxph(Surv(Time, death) ~ drug + AZT + prevOI + gender, 
    data = aids.id, x = TRUE)
fit.aids <- jointModel(fitLME, fitSURV, timeVar = "obstime", 
    method = "piecewise-PH-aGH")

ND <- aids[aids$patient == "7", ]
roc <- rocJM(fit.aids, dt = c(2, 4, 8), ND, idVar = "patient")
plot(roc, lwd = 2, legend = TRUE)
plot(roc, type = "AUC")
}
}

\keyword{methods}
