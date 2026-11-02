# Synthetic teaching example: pre-specify the hypothesis and two independent groups.
# Rank inference tests identical distributions, not P(superiority)=.5 under
# unequal shapes/spreads. A median/mean interpretation needs stronger assumptions.
check_settings <- function(alternative,conf.level,exact,location_shift) {
 stopifnot(is.character(alternative),length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1,
  is.logical(exact),length(exact)==1L,!is.na(exact),
  is.logical(location_shift),length(location_shift)==1L,!is.na(location_shift))
}
# Ordered numeric codes may represent ordinal levels when location_shift=FALSE.
# location_shift=TRUE is an explicit interval-scale, common-shape shift assumption;
# never enable it for ordinal codes. It requests a separate shift interval, not
# an interval for the probability of superiority.
analyse <- function(x, y, alternative = "two.sided", conf.level = 0.95, exact = TRUE, location_shift = FALSE) {
    spec <- infercalcR::method_spec("rank_sum_mann_whitney", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, y = y, exact = exact, location_shift = location_shift)$payload$native
}
example_x <- c(2,4,6,8,10)
example_y <- c(1,3,5,7,9)
example_result <- analyse(example_x,example_y,alternative,confidence_level,location_shift=FALSE)
print(example_result)
