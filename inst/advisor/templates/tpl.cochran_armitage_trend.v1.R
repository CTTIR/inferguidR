# Independent binomial groups with an order and scores fixed before observing outcomes.
# Score spacing is part of the hypothesis; arbitrary monotone recoding can change the test.
analyse <- function(events, totals, scores, scores_prespecified = FALSE, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("cochran_armitage_trend", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, events = events, totals = totals, scores = scores, scores_prespecified = scores_prespecified)$payload$native
}
example_result<-analyse(c(low=3,medium=8,high=13),c(low=20,medium=25,high=30),c(low=0,medium=1,high=2),scores_prespecified=TRUE,alternative=alternative,conf.level=confidence_level)
print(example_result)
