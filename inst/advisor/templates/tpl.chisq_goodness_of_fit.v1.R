# Synthetic teaching example for prespecified count-data inference.
# Counts must represent independent units, except for explicitly paired McNemar tables.
# No missing observations, survey weights, fitted null probabilities or automatic test switching.
check_settings <- function(alternative, conf.level) {
 stopifnot(identical(alternative,"two.sided"),is.numeric(conf.level),length(conf.level)==1L,
           is.finite(conf.level),conf.level>0,conf.level<1)
}
# Null: the prespecified category probabilities p; no parameters estimated from these counts.
# Pearson's omnibus statistic uses a right-tail chi-squared reference distribution.
analyse <- function(counts, p, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("chisq_goodness_of_fit", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, counts = counts, p = p)$payload$native
}
example_result<-analyse(c(A=20,B=15,C=25),c(.3,.3,.4),alternative,confidence_level)
print(example_result)
