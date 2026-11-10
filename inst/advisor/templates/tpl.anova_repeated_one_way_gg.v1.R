# Synthetic complete repeated-measures example, one row per independent subject.
# Prespecified one-way occasion effect, with routine Greenhouse-Geisser correction.
# No Mauchly pretest, covariance-model selection or automatic post-hoc comparisons.
analyse <- function(data, occasion_columns, independent_subjects = FALSE, alternative = "two.sided",
    conf.level = 0.95, complete_cases = FALSE) {
    spec <- infercalcR::method_spec("anova_repeated_one_way_gg", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, occasion_columns = occasion_columns, independent_subjects = independent_subjects,
        complete_cases = complete_cases)$payload$native
}
example_data<-data.frame(id=1:6,before=c(5,7,6,8,4,9),during=c(7,8,6,11,6,10),after=c(8,8,9,12,7,13))
example_result<-analyse(example_data,c('before','during','after'),TRUE,alternative,confidence_level)
print(example_result)
