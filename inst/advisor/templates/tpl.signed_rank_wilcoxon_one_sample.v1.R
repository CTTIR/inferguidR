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
analyse <- function(x, mu, alternative = "two.sided", conf.level = 0.95, exact = TRUE) {
    spec <- infercalcR::method_spec("signed_rank_wilcoxon_one_sample", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, mu = mu, exact = exact)$payload$native
}
example_x <- c(-1,2,-3,4,5,-6,7,8,9,10)
example_result <- analyse(example_x,0,alternative,confidence_level)
print(example_result)
