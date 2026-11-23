% Help page of a function that infercalcR implements and inferguidR re-exports
% (R/reexports.R). Maintained by hand: roxygen2 has no source block for it here.
\name{ig_status}
\alias{ig_status}
\title{Create a typed status}
\usage{
ig_status(code, reason = "", details = list())
}
\arguments{
\item{code}{One of \code{\link[=ig_status_codes]{ig_status_codes()}}.}

\item{reason}{A human-readable explanation. Required (non-empty) for every
code except \code{"ok"}.}

\item{details}{A named list of structured details (for example the list of
problems found during validation). Values must be serialisable data.}
}
\value{
An object of class \code{ig_status}.
}
\description{
Create a typed status
}
\examples{
ig_status("unsupported", "engine for t_generic not implemented yet")
}
\seealso{
Other status helpers:
\code{\link[=ig_status_codes]{ig_status_codes()}},
\code{\link[=is_ig_status]{is_ig_status()}}
}
\concept{status helpers}
