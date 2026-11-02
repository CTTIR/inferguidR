# Synthetic teaching example. The selected method and hypothesis are fixed before analysis.
# Reject missing and infinite inputs; document an explicit missingness strategy first.
check_settings <- function(alternative, conf.level) {
 stopifnot(length(alternative)==1L, alternative %in% c("two.sided","greater","less"),
           is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Independent binary trials with one fixed event definition and reference probability.
# The interval estimates the event probability, not a difference from the reference.
analyse <- function(events, n, p0, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("binomial_exact_one_sample", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, events = events, n = n, p0 = p0)$payload$native
}
example_result <- analyse(7,20,.2,alternative,confidence_level)
print(example_result)
graphics::barplot(c(event=7,non_event=13),ylab="Number of independent trials")
