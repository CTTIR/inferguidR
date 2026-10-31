# Prespecified numeric predictors, independent observations and explicit event/exposure coding.
# Asymptotic likelihood-ratio omnibus test and marginal profile-likelihood slope intervals.
# Marginal Wald coefficient intervals are retained separately on the link scale.
# No automatic selection, robust/penalized fallback, or causal interpretation.
analyse <- function(data, predictors, independent_units = FALSE, alternative = "two.sided", conf.level = 0.95,
    complete_cases = FALSE) {
    spec <- infercalcR::method_spec("regression_logistic_binary", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, predictors = predictors, independent_units = independent_units,
        complete_cases = complete_cases)$payload$native
}
example_data<-data.frame(id=1:24,x1=rep(c(-1,0,1),8),x2=rep(c(0,1),each=12),y=c(0,0,1,0,1,1,0,1,0,1,0,1,0,1,1,1,0,1,0,1,0,1,1,1))
example_result<-analyse(example_data,c('x1','x2'),TRUE,alternative,confidence_level)
print(example_result)
