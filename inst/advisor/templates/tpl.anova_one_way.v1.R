# Synthetic independent-group omnibus example; diagnostics never choose a test.
# Estimand: equality of population means, not a signed pairwise contrast.
# Independent normal errors are required for exact pooled F inference;
# Welch uses approximate degrees of freedom and permits unequal variances.
# Equal variances for the pooled test must be justified before inspecting data.
analyse <- function(x, group, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("anova_one_way", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, group = group)$payload$native
}
example_x<-c(1,2,3,4,5, 2,3,5,8,9, 4,5,6,6,7)
example_group<-rep(c("A","B","C"),each=5)
example_result<-analyse(example_x,example_group,alternative,confidence_level)
print(example_result)
