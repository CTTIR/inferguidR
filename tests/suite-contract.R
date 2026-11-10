library(inferguidR)
s <- analysis_spec("t_means_two_independent", "a_priori",
  inputs = list(alpha = .05, power = .8, d = .5), options = list(tails = "two"))
x <- solve_design(s)
stopifnot(x$status$code == "ok")
before <- ig_to_json(x)
common <- infercalcR::as_inference(x)
stopifnot(infercalcR::result_metadata(common)$kind == "planning",
  identical(ig_to_json(x), before), all(is.na(infercalcR::result_data(common)$p.value)),
  identical(infercalcR::result_data(x), infercalcR::result_data(common)))
back <- ig_from_json(before)
stopifnot(identical(infercalcR::result_data(back), infercalcR::result_data(x)))
mapping <- planning_methods(c("t_means_two_independent_welch", "t_means_two_independent_student"))
stopifnot(identical(mapping$relation, c("comparator", "same_procedure_under_assumptions")))
rank <- planning_methods("rank_sum_mann_whitney")
stopifnot(rank$relation == "approximation")
mc <- planning_methods(c("mcnemar_exact_conditional", "mcnemar_asymptotic_cc"))
stopifnot(!identical(mc$relation[1], mc$relation[2]))
cat("Shared planning access and explicit method relations checked\n")
