# Synthetic independent-group omnibus example; diagnostics never choose a test.
# Rank omnibus test under identical group distributions, with average tied ranks.
# It is not a test of medians unless a common-shape location-shift model is justified.
# A rejection with unequal shapes/spreads need not establish stochastic ordering.
analyse <- function(x, group, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("kruskal_wallis", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, group = group)$payload$native
}
example_x<-c(1,2,3,4,5, 2,3,5,8,9, 4,5,6,6,7)
example_group<-rep(c("A","B","C"),each=5)
example_result<-analyse(example_x,example_group,alternative,confidence_level)
print(example_result)
