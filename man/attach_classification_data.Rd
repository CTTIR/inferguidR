% Help page of a function that infercalcR implements and inferguidR re-exports
% (R/reexports.R). Maintained by hand: roxygen2 has no source block for it here.
\name{attach_classification_data}
\alias{attach_classification_data}
\title{Reattach scored data to a classification specification}
\usage{
attach_classification_data(spec, data)
}
\arguments{
\item{spec}{A \code{classification_spec} with input kind \code{scored_data}.}

\item{data}{The original data frame.}
}
\value{
The specification with the data attached.
}
\description{
Specifications serialised to JSON carry a SHA-256 fingerprint of the
analysed data, never the rows. To evaluate such a specification again,
reattach the original data: the columns named in the model records are read
and the fingerprint must match exactly.
}
\examples{
d <- data.frame(y = c(1, 1, 0, 0, 0), s = c(3, 2, 2, 1, 0))
spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 2)
back <- ig_from_json(ig_to_json(analysis_envelope(spec)))$spec
spec2 <- attach_classification_data(back, d)
identical(spec2$spec_hash, spec$spec_hash)
}
\seealso{
\itemize{
\item \code{\link[=classification_spec]{classification_spec()}}
\item \code{\link[=evaluate_cutoff]{evaluate_cutoff()}}
\item \code{\link[=plot_classification]{plot_classification()}}
}
}
\concept{classification}
