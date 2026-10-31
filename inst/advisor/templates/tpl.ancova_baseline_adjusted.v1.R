# Prespecified numeric predictors and an intercept; independent units, no automatic selection.
# Coefficients describe conditional-mean associations, not causal effects without design support.
# Classical t/F inference assumes linear mean, constant residual variance and normal errors.
# Residual diagnostics do not automatically switch models or standard errors.
# The adjustment set and common slopes must be prespecified, not selected by diagnostics.
analyse <- function(data, covariates, group_levels, common_slopes, conf.level = 0.95, missing = "error",
    alternative = "two.sided") {
    spec <- infercalcR::method_spec("ancova_baseline_adjusted", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, covariates = covariates, group_levels = group_levels,
        common_slopes = common_slopes, missing = missing)$payload$native
}
example_data<-data.frame(id=1:18,group=rep(c("control","intervention","comparison"),6),baseline=c(1,3,2,4,5,2,6,3,5,4,7,8,5,9,6,8,7,10),y=c(4,7,5,6,9,6,10,8,9,7,11,14,9,13,11,12,14,16))
example_groups<-c("control","intervention","comparison")[seq_len(example_group_count)]
example_data<-example_data[example_data$group %in% example_groups,]
example_result<-analyse(example_data,"baseline",example_groups,common_slopes=TRUE,conf.level=confidence_level,alternative=alternative)
print(example_result)
