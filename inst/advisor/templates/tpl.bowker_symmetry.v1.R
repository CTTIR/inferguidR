# Two matched occasions, same named categories in the same order on both axes.
# Cells count participants; each participant contributes once to a confirmed pair.
# Independent participants, complete observed pairs, chi-squared asymptotic inference.
# Missing pairs must be handled and reported before constructing this aggregate table.
analyse <- function(counts, alternative = "two.sided", conf.level = 0.95, n_excluded = 0L) {
    spec <- infercalcR::method_spec("bowker_symmetry", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, counts = counts, n_excluded = n_excluded)$payload$native
}
example_counts<-matrix(c(20,4,6,8,25,3,5,7,18),3,byrow=TRUE,
 dimnames=list(c("A","B","C"),c("A","B","C")))
example_result<-analyse(example_counts,alternative,confidence_level)
print(example_result)
