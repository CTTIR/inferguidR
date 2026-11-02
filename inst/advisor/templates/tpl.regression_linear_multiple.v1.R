# Prespecified numeric predictors and an intercept; independent units, no automatic selection.
# Coefficients describe conditional-mean associations, not causal effects without design support.
# Classical t/F inference assumes linear mean, constant residual variance and normal errors.
# Residual diagnostics do not automatically switch models or standard errors.
analyse <- function(data, predictors = c("x1", "x2"), tested_predictors = NULL, alternative = "two.sided",
    conf.level = 0.95, missing = "error") {
    spec <- infercalcR::method_spec("regression_linear_multiple", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, predictors = predictors, tested_predictors = tested_predictors,
        missing = missing)$payload$native
}
example_data<-data.frame(id=1:12,x1=1:12,x2=c(1,0,2,1,3,0,1,2,0,3,2,1),y=c(3,4,4,7,6,9,8,12,10,13,12,16))
example_result<-analyse(example_data,tested_predictors="x2",alternative=alternative,conf.level=confidence_level)
print(example_result)
