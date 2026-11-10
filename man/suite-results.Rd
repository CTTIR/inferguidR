\name{suite-results}
\alias{as_inference.power_result}
\alias{result_data.power_result}
\alias{result_metadata.power_result}
\title{Access planning results through the shared inference contract}
\description{View existing power results as finalized planning rows without changing historical objects or their hashes.}
\usage{
\method{as_inference}{power_result}(x, ...)
\method{result_data}{power_result}(x, ...)
\method{result_metadata}{power_result}(x, ...)
}
\arguments{\item{x}{A power_result.}\item{...}{Reserved for future methods; unused.}}
\details{Solved, achieved and numerical planning quantities become explicitly named planning rows. No artificial p-values or confidence intervals are introduced. The complete original result is retained as payload. Old JSON must first be read with inferguidR::ig_from_json(), which loads the registered methods. Convert once with infercalcR::as_inference() to use the common result in a session without guidance installed. Historical JSON round trips are separate from new result provenance.}
\value{A common planning_result, data.frame or metadata list.}
