% Help page of a function that infercalcR implements and inferguidR re-exports
% (R/reexports.R). Maintained by hand: roxygen2 has no source block for it here.
\name{binormal_rates}
\alias{binormal_rates}
\title{Analytic rates of the binormal model at given cutoffs}
\usage{
binormal_rates(cutoff, auc, sigma = 1)
}
\arguments{
\item{cutoff}{Numeric vector of real-valued cutoffs (alert if score >=
cutoff) on the score scale standardised to the non-event distribution.}

\item{auc}{Area under the ROC curve in (0, 1).}

\item{sigma}{Ratio of the event to the non-event standard deviation (> 0).}
}
\value{
A data frame with \code{cutoff}, \code{tpr} and \code{fpr}, computed from the
normal upper tails (no discretisation).
}
\description{
Analytic rates of the binormal model at given cutoffs
}
\examples{
binormal_rates(c(0, 1, 2), auc = 0.8, sigma = 1)
}
\seealso{
\itemize{
\item \code{\link[=classification_spec]{classification_spec()}}
\item \code{\link[=evaluate_cutoff]{evaluate_cutoff()}}
\item \code{\link[=plot_classification]{plot_classification()}}
}
}
\concept{classification}
