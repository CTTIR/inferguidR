# Synthetic teaching example for an independently prespecified spread question.
# Do not use a spread test as a preliminary rule for selecting a mean test.
check_settings <- function(alternative, conf.level, omnibus=FALSE) {
 stopifnot(is.character(alternative),length(alternative)==1L,!is.na(alternative),
  alternative%in%c("two.sided","greater","less"),is.numeric(conf.level),length(conf.level)==1L,
  is.finite(conf.level),conf.level>0,conf.level<1)
 if(omnibus) stopifnot(identical(alternative,"two.sided"))
}
# Estimand: population variance of x divided by population variance of y.
# Exact F inference requires independent normal samples; it is sensitive to non-normality.
analyse <- function(x, y, ratio0 = 1, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("variance_f_test_two", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, y = y, ratio0 = ratio0)$payload$native
}
example_x<-c(1,2,4,6,7,9);example_y<-c(2,3,4,4,5,6)
example_result<-analyse(example_x,example_y,alternative=alternative,conf.level=confidence_level)
print(example_result)
