# Synthetic teaching example: Fisher's fixed-margin probability-ordering test.
# Independent units, mutually exclusive categories; both margins conditioned on.
# The alternative is omnibus. No odds-ratio interval is claimed for a general table.
# This teaching script limits total count to 100 and each dimension to 5.
# Larger exact problems need an explicit resource plan, not automatic approximation.
analyse <- function(counts, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("fisher_exact_rxc", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, counts = counts)$payload$native
}
example_counts<-rbind(Group1=c(A=3,B=1,C=4),Group2=c(A=1,B=5,C=2))
example_result<-analyse(example_counts,alternative,confidence_level)
print(example_result)
