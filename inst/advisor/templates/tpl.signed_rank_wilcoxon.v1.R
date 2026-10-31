# Synthetic teaching example: select the hypothesis before inspecting results.
# Independent units; finite interval-level observations, not ordinal category codes.
check_settings <- function(alternative,conf.level,exact) {
 stopifnot(is.character(alternative),length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  is.logical(exact),length(exact)==1L,!is.na(exact))
}
# Symmetry about mu, not a mean test or an unrestricted median test.
# exact=TRUE is a deliberate small-sample, no-zero/no-tie policy. Request
# exact=FALSE explicitly for the continuity-corrected normal approximation.
signed_analysis <- function(x, mu, alternative, conf.level, exact) {
    spec <- infercalcR::method_spec("signed_rank_wilcoxon_one_sample", alternative, conf.level)
    infercalcR::run_inference(spec, x = x, mu = mu, exact = exact)$payload$native
}
# One wide row per confirmed unique participant, difference x minus y.
analyse <- function(data, alternative = "two.sided", conf.level = 0.95, exact = TRUE) {
    spec <- infercalcR::method_spec("signed_rank_wilcoxon", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, exact = exact)$payload$native
}
example_data <- data.frame(id=paste0("unit",1:10),x=c(-1,2,-3,4,5,-6,7,8,9,10),y=0)
example_result <- analyse(example_data,alternative,confidence_level)
print(example_result)
