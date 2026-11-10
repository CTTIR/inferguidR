# Prespecified numeric predictors, independent observations and explicit event/exposure coding.
# Asymptotic likelihood-ratio omnibus test and marginal Wald coefficient intervals.
# No automatic selection, robust/penalized fallback, or causal interpretation.
analyse <- function(data, predictors, independent_units = FALSE, alternative = "two.sided", conf.level = 0.95,
    complete_cases = FALSE) {
    spec <- infercalcR::method_spec("regression_negative_binomial", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, predictors = predictors, independent_units = independent_units,
        complete_cases = complete_cases)$payload$native
}
example_data<-data.frame(id=1:18,x1=rep(c(-1,0,1),6),x2=rep(c(0,1),each=9),exposure=rep(c(1,2,3),6),y=c(0,0,18,1,12,0,0,2,16,0,21,2,9,0,34,1,3,0))
example_result<-analyse(example_data,c('x1','x2'),TRUE,alternative,confidence_level)
print(example_result)
