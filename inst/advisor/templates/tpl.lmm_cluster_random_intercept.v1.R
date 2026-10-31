# Prespecified Gaussian random-intercept model; covariance choice is explicit.
# REML estimates, marginal fixed-effect F test, nlme containment degrees of freedom.
# No comparison of REML likelihoods across different fixed-effect models.
analyse <- function(data, independent_clusters = FALSE, random_intercept_justified = FALSE, alternative = "two.sided",
    conf.level = 0.95, missing_outcomes = "error") {
    spec <- infercalcR::method_spec("lmm_cluster_random_intercept", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, independent_clusters = independent_clusters, random_intercept_justified = random_intercept_justified,
        missing_outcomes = missing_outcomes)$payload$native
}
example_data<-expand.grid(replicate=1:3,cluster=1:12);example_data$id<-seq_len(nrow(example_data))
example_data$factor<-factor(ifelse(example_data$cluster<=6,'control','treatment'),levels=c('control','treatment'))
example_data$y<-rep(c(-3,-2,-1,1,2,3,0,1,2,4,5,6),each=3)+rep(c(-1,0,1),12)
example_result<-analyse(example_data,TRUE,TRUE,alternative,confidence_level)
print(example_result)
