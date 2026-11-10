# Synthetic teaching example for prespecified count-data inference.
# Counts must represent independent units, except for explicitly paired McNemar tables.
# No missing observations, survey weights, fitted null probabilities or automatic test switching.
check_settings <- function(alternative, conf.level) {
 stopifnot(identical(alternative,"two.sided"),is.numeric(conf.level),length(conf.level)==1L,
           is.finite(conf.level),conf.level>0,conf.level<1)
}
# Rows are independent groups and columns are mutually exclusive outcome categories.
# Null: row and column independence. Pearson's test is uncorrected and omnibus.
analyse <- function(counts, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("chisq_independence_rxc", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, counts = counts)$payload$native
}
example_counts<-rbind(Group1=c(A=30,B=40,C=30),Group2=c(A=45,B=35,C=20))
example_result<-analyse(example_counts,alternative,confidence_level)
print(example_result)
