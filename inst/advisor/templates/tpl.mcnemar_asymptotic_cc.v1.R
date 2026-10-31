# Synthetic teaching example for prespecified count-data inference.
# Counts must represent independent units, except for explicitly paired McNemar tables.
# No missing observations, survey weights, fitted null probabilities or automatic test switching.
check_settings <- function(alternative, conf.level) {
 stopifnot(identical(alternative,"two.sided"),is.numeric(conf.level),length(conf.level)==1L,
           is.finite(conf.level),conf.level>0,conf.level<1)
}
# Rows are time1 and columns time2, with event first at both times.
# Each participant contributes exactly one cell. The null is equal discordance probabilities.
# The descriptive marginal difference is time1 event probability minus time2 event probability.
analyse <- function(counts, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("mcnemar_asymptotic_cc", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, counts = counts)$payload$native
}
example_counts<-rbind(Event=c(Event=45,Non_event=20),Non_event=c(Event=10,Non_event=25))
example_result<-analyse(example_counts,alternative,confidence_level)
print(example_result)
