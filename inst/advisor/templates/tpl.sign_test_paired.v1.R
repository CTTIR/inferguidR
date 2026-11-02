# Synthetic teaching example: fix the hypothesis and method before analysis.
# Independent evaluation units are required; no survey or fractional weights.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative), length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Estimand: probability of a positive difference among nonzero differences.
# Interpreting this as a median test needs the corresponding continuity/design assumptions.
analyse <- function(data, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("sign_test_paired", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data)$payload$native
}
example_data <- data.frame(id=paste0("unit",1:10),x=c(5,6,7,5,3,8,7,4,5,9),y=c(4,4,6,5,4,6,6,5,4,7))
example_result <- analyse(example_data,alternative,confidence_level)
print(example_result)
