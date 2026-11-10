# Synthetic fixed-factor example: prespecify the factorial model and hypotheses.
# Two or three crossed categorical factors; full interactions, no empty cells.
# Type-III hypotheses use sum-to-zero contrasts and unweighted marginal means.
# All term tests are reported, with Holm adjustment across this declared family.
check_settings <- function(alternative,conf.level) {
 stopifnot(identical(alternative,"two.sided"),is.numeric(conf.level),length(conf.level)==1L,
  is.finite(conf.level),conf.level>0,conf.level<1)
}
analyse <- function(data, factor_columns, independent_units = FALSE, equal_variance_justified = FALSE,
    alternative = "two.sided", conf.level = 0.95, complete_cases = FALSE) {
    spec <- infercalcR::method_spec("anova_factorial", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, factor_columns = factor_columns, independent_units = independent_units,
        equal_variance_justified = equal_variance_justified, complete_cases = complete_cases)$payload$native
}
example_data<-expand.grid(A=factor(c("a1","a2")),B=factor(c("b1","b2","b3")),replicate=1:3)
example_data$id<-seq_len(nrow(example_data))
example_data$y<-2*as.integer(example_data$A)+as.integer(example_data$B)+
 as.integer(example_data$A)*as.integer(example_data$B)+c(-1,0,1)[example_data$replicate]
example_result<-analyse(example_data,c("A","B"),TRUE,TRUE,alternative,confidence_level)
print(example_result)
