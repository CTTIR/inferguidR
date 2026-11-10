# Synthetic teaching example: fix the hypothesis and method before analysis.
# Independent evaluation units are required; no survey or fractional weights.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative), length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Estimand: probability of a positive difference among nonzero differences.
# Interpreting this as a median test needs the corresponding continuity/design assumptions.
analyse <- function(x, mu0, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("sign_test_one_sample", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, mu0 = mu0)$payload$native
}
example_x <- c(5,6,7,4,3,8,7,2,5,9)
example_reference <- 4
example_result <- analyse(example_x,example_reference,alternative,confidence_level)
print(example_result)
