% Help page of a function that infercalcR implements and inferguidR re-exports
% (R/reexports.R). Maintained by hand: roxygen2 has no source block for it here.
\name{calculator_functions}
\alias{calculator_functions}
\title{Functions of the distribution calculator}
\usage{
calculator_functions()
}
\value{
A data frame with columns \code{name}, \code{arguments}, \code{n_args},
\code{category} (\code{function} or \code{distribution}), \code{description} and
\code{gpower_note} (\code{NA} when there is nothing to note).
}
\description{
Lists the functions that \code{\link[=distribution_calculator]{distribution_calculator()}} accepts, with their
arguments in the order of the 'G\*Power' 3.1.9.7 calculator, a description
and notes on reference behaviour where 'G\*Power' differs.
}
\examples{
fns <- calculator_functions()
fns[fns$category == "distribution", c("name", "arguments")][1:6, ]
}
\seealso{
\code{\link[=distribution_calculator]{distribution_calculator()}}
}
