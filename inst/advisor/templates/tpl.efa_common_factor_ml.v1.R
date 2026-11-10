# Exploratory common-factor parallel analysis with squared multiple correlations.
# Independent column permutations preserve marginal values and break association.
# The result supports a factor-count decision; it never fits or selects an ML model.
parallel_analysis <- function(data, variables, independent_units = FALSE, continuous_measurements = FALSE,
    B = 199L, seed = 104729L, quantile_probability = 0.95, complete_cases = FALSE, permutations = NULL) {
    infercalcR::parallel_analysis(data = data, variables = variables, independent_units = independent_units,
        continuous_measurements = continuous_measurements, B = B, seed = seed, quantile_probability = quantile_probability,
        complete_cases = complete_cases, permutations = permutations)$payload$native
}

# Prespecified number of common factors; maximum likelihood and oblique promax rotation.
# Numeric continuous measurements; no PCA substitution or Kaiser-only selection.
analyse <- function(data, variables, factors, independent_units = FALSE, continuous_measurements = FALSE,
    alternative = "two.sided", conf.level = 0.95, complete_cases = FALSE) {
    spec <- infercalcR::method_spec("efa_common_factor_ml", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, variables = variables, factors = factors, independent_units = independent_units,
        continuous_measurements = continuous_measurements, complete_cases = complete_cases)$payload$native
}
t<-seq_len(120);f1<-sin(t*.13);f2<-.3*f1+cos(t*.21)
example_data<-data.frame(id=t)
for(j in 1:6)example_data[[paste0('v',j)]]<-(if(j<=3).8*f1+.2*f2 else .2*f1+.8*f2)+.7*sin(t*(.61+j*.1))
example_result<-analyse(example_data,paste0('v',1:6),2,TRUE,TRUE)
print(example_result)

example_parallel<-parallel_analysis(example_data,paste0('v',1:6),TRUE,TRUE,B=100L)
print(example_parallel[c('status','observed_eigenvalues','null_thresholds','suggested_leading_factors','limitations')])
