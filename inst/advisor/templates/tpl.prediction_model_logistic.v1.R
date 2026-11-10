# Prediction development with prespecified predictors and bootstrap optimism correction.
# No variable selection, cutoff optimization, or external-validation claim.
analyse <- function(data, predictors, independent_units = FALSE, alternative = NULL, conf.level = 0.95,
    B = 200L, seed = 104729L, complete_cases = FALSE, resamples = NULL) {
    if (!is.null(alternative))
        stop("This procedure has no hypothesis alternative")
    spec <- infercalcR::method_spec("prediction_model_logistic", alternative = "two.sided", conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, predictors = predictors, independent_units = independent_units,
        B = B, seed = seed, complete_cases = complete_cases, resamples = resamples)$payload$native
}
predict_new <- function(result, newdata, interval = c("none", "risk"), conf.level = result$requested_confidence_level) {
    interval <- match.arg(interval)
    common <- infercalcR::new_inference_result(data.frame(), list(method_id = "prediction_model_logistic"),
        payload = list(native = result))
    value <- infercalcR::predict_inference(common, newdata, interval = interval, conf.level = conf.level)$payload$native
    if (identical(interval, "none"))
        value$estimate
    else value
}
example_data<-data.frame(id=1:200,x1=rep(-4:5,20),x2=rep(c(0,1),each=100));example_data$y<-as.numeric(((seq_len(200)*37)%%101)/101 < stats::plogis(-.4+.18*example_data$x1+.35*example_data$x2))
example_result<-analyse(example_data,c('x1','x2'),TRUE,conf.level=confidence_level,B=200L)
print(example_result[c('status','method','apparent','corrected','monte_carlo_se','B','n_used')])
