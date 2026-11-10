# Prespecified numeric predictors and explicit category coding; no selection.
# Independent observations; asymptotic likelihood-ratio and marginal Wald inference.
analyse <- function(data, predictors, independent_units = FALSE, reference_category = NULL, alternative = "two.sided",
    conf.level = 0.95, complete_cases = FALSE) {
    spec <- infercalcR::method_spec("regression_multinomial", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, predictors = predictors, independent_units = independent_units,
        reference_category = reference_category, complete_cases = complete_cases)$payload$native
}
example_data<-data.frame(id=1:60,x1=rep(c(-1,0,1),each=20),y=factor(c(rep(1,10),rep(2,6),rep(3,4),rep(1,6),rep(2,8),rep(3,6),rep(1,3),rep(2,6),rep(3,11)),levels=1:3,ordered=FALSE))
example_result<-analyse(example_data,'x1',TRUE,reference_category='1',alternative=alternative,conf.level=confidence_level)
print(example_result)
