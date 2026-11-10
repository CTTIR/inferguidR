% Compatibility help for S3 methods provided by infercalcR.
\name{ig_immutable}
\alias{ig_immutable}
\alias{$<-.ig_frozen}
\alias{[[<-.ig_frozen}
\alias{[<-.ig_frozen}
\title{Immutable inferguidR objects}

\arguments{
\item{x}{An immutable inferguidR object.}

\item{name, i, ...}{Ignored.}

\item{value}{Ignored.}
}
\value{
Never returns; always signals an error.
}
\description{
Specifications, results, envelopes, report snapshots and recommendations
cannot be modified in place. Assignment with \verb{$<-}, \verb{[[<-} or \verb{[<-}
signals an \code{invalid_input} error; use \code{\link[=update.power_spec]{update.power_spec()}} or the
constructors to create a new object.
}
\keyword{internal}
\note{The assignment methods for the legacy \code{ig_frozen} class are provided and registered by \pkg{infercalcR}, which is loaded with this package. The historical class and error behavior are preserved.}
