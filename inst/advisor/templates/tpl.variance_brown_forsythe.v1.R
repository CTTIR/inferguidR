# Synthetic teaching example for an independently prespecified spread question.
# Do not use a spread test as a preliminary rule for selecting a mean test.
check_settings <- function(alternative, conf.level, omnibus=FALSE) {
 stopifnot(is.character(alternative),length(alternative)==1L,!is.na(alternative),
  alternative%in%c("two.sided","greater","less"),is.numeric(conf.level),length(conf.level)==1L,
  is.finite(conf.level),conf.level>0,conf.level<1)
 if(omnibus) stopifnot(identical(alternative,"two.sided"))
}
# Brown-Forsythe: ordinary one-way ANOVA on absolute deviations from group medians.
# This is the median-centred Levene test, not the Welch/Brown-Forsythe test of means.
# No sampling interval for a common variance ratio is supplied by this omnibus test.
analyse <- function(x, group, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("variance_brown_forsythe", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, group = group)$payload$native
}
example_x<-c(1,2,3,4,5, 1,2,4,7,9, 2,3,3,4,4)
example_group<-rep(c("A","B","C"),each=5)
example_result<-analyse(example_x,example_group,alternative,confidence_level)
print(example_result)
