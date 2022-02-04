\name{plot.survfitJM}

\alias{plot.survfitJM}

\title{Plot Method for survfitJM Objects}

\description{
  Produces plots of conditional probabilities of survival.
}

\usage{
\method{plot}{survfitJM}(x, estimator = c("both", "mean", "median"), 
    which = NULL, fun = NULL, conf.int = FALSE, 
    fill.area = FALSE, col.area = "grey", col.abline = "black", col.points = "black",
    add.last.time.axis.tick = FALSE, include.y = FALSE, main = NULL, 
    xlab = NULL, ylab = NULL, ylab2 = NULL, lty = NULL, col = NULL, 
    lwd = NULL, pch = NULL, ask = NULL, legend = FALSE, \dots,
    cex.axis.z = 1, cex.lab.z = 1)
}

\arguments{
  \item{x}{an object inheriting from class \code{survfitJM}.}
  \item{estimator}{character string specifying, whether to include in the plot the mean of the conditional probabilities of survival,
    the median or both. The mean and median are taken as estimates of these conditional probabilities over the M replications of the
    Monte Carlo scheme described in \code{\link{survfitJM}}.}
  \item{which}{a numeric or character vector specifying for which subjects to produce the plot. If a character vector, then is
    should contain a subset of the values of the \code{idVar} variable of the \code{newdata} argument of \code{\link{survfitJM}}.}
  \item{fun}{a vectorized function defining a transformation of the survival curve. For example with \code{fun=log} the log-survival curve 
    is drawn.}
  \item{conf.int}{logical; if \code{TRUE}, then a pointwise confidence interval is included in the plot.}
  \item{fill.area}{logical; if \code{TRUE} the area defined by the confidence interval of the survival function is put in color.}
  \item{col.area}{the color of the area defined by the confidence interval of the survival function.}
 \item{col.abline,col.points}{the color for the vertical line and the points when \code{include.y} is \code{TRUE}.} 
\item{add.last.time.axis.tick}{logical; if \code{TRUE}, a tick is added in the x-axis for the last available time point for which a
    longitudinal measurement was available.}
  \item{include.y}{logical; if \code{TRUE}, two plots are produced per subject, i.e., the plot of conditional probabilities of survival 
    and a scatterplot of his longitudinal measurements.}
  \item{main}{a character string specifying the title in the plot.}
  \item{xlab}{a character string specifying the x-axis label in the plot.}
  \item{ylab}{a character string specifying the y-axis label in the plot.}
  \item{ylab2}{a character string specifying the y-axis label in the plotm when \code{include.y = TRUE}.}
  \item{lty}{what types of lines to use.}
  \item{col}{which colors to use.}
  \item{lwd}{the thickness of the lines.}
  \item{pch}{the type of points to use.}  
  \item{ask}{logical; if \code{TRUE}, the user is asked before each plot, see \code{par()}.}
  \item{legend}{logical; if \code{TRUE}, a legend is included in the plot.}
  \item{cex.axis.z, cex.lab.z}{the par \code{cex} argument for the axis at side 4, when \code{include.y = TRUE}.}
  \item{\dots}{extra graphical parameters passed to \code{plot()}.}
}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for longitudinal and time-to-event data. 
\emph{Biometrics} \bold{67}, 819--829.

Rizopoulos, D. (2010) JM: An R Package for the Joint Modelling of Longitudinal and Time-to-Event Data. 
\emph{Journal of Statistical Software} \bold{35} (9), 1--33. \doi{10.18637/jss.v035.i09}
}


\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{survfitJM}}}

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
plot(ss)
plot(ss, include.y = TRUE, add.last.time.axis.tick = TRUE, legend = TRUE)
}

\keyword{methods}
