# Independent ordinal/continuous samples: theta = P(X > Y) + 0.5 P(X = Y).
# This tests stochastic equality, not equality of medians or identical distributions.
# A t approximation is used explicitly; no automatic exact/permutation or normal fallback.
analyse <- function(x, y, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("brunner_munzel", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, y = y)$payload$native
}
example_x<-(seq_len(14)*7+seq_len(14)^2)%%19/3
example_y<-(seq_len(11)*11+2*seq_len(11)^2)%%23/3
example_result<-analyse(example_x,example_y,alternative=alternative,conf.level=confidence_level)
print(example_result)
