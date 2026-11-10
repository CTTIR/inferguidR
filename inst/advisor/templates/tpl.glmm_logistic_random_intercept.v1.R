# Prespecified two-group clustered binary example, with explicit estimand.
# A unique id identifies each observation; cluster identifies dependent units.
analyse <- function(data, independent_clusters = FALSE, model_assumptions_justified = FALSE, alternative = "two.sided",
    conf.level = 0.95, missing_outcomes = "error") {
    spec <- infercalcR::method_spec("glmm_logistic_random_intercept", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, independent_clusters = independent_clusters, model_assumptions_justified = model_assumptions_justified,
        missing_outcomes = missing_outcomes)$payload$native
}
counts<-c(0,1,1,2,2,3,4,4,5,6,0,2,3,5,6,1,2,3,3,4,4,5,5,6,6,2,3,4,5,6)
example_data<-data.frame(id=1:180,cluster=rep(1:30,each=6),group=factor(rep(c('control','treatment'),each=90),levels=c('control','treatment')),y=unlist(lapply(counts,function(s)c(rep(1,s),rep(0,6-s)))))
example_result<-analyse(example_data,TRUE,TRUE,alternative,confidence_level)
print(example_result)
