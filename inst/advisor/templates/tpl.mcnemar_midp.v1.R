# Synthetic teaching example: fix the hypothesis and method before analysis.
# Independent evaluation units are required; no survey or fractional weights.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative), length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# n10: x positive, y negative; n01: x negative, y positive, in the same units.
# The first and second measurements and their positive labels must be fixed in advance.
analyse <- function(counts, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("mcnemar_midp", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, counts = counts)$payload$native
}
example_counts <- c(n00=20,n01=3,n10=9,n11=15)
example_result <- analyse(example_counts,alternative,confidence_level)
print(example_result)
