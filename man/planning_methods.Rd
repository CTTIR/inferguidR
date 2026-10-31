\name{planning_methods}
\alias{planning_methods}
\title{Explicit relations between observed procedures and planning engines}
\description{List reviewed procedure identities and distinguish matched planning assumptions, approximations, different-procedure comparators and unavailable mappings.}
\usage{planning_methods(method_id = NULL)}
\arguments{\item{method_id}{Optional character vector of observed-method IDs. NULL returns all inventoried identities.}}
\details{This is a relation catalogue, not automatic conversion of observed data to prospective assumptions. The existing advisor_power_spec bridge remains specifically scoped to Boschloo. A comparator or approximation is never silently selected. Welch and pooled tests, exact and asymptotic McNemar, and finite-sample rank inference retain distinct identities.}
\value{A data frame with mapping version, observed method, planning engine, relation, required prospective inputs, assumptions and supported conversion.}
\examples{planning_methods(c("t_means_two_independent_welch", "t_means_two_independent_student"))}
