# Exploratory composite scores, not latent factors or confirmatory inference.
# Choose variables and retained component count explicitly; no eigenvalue cutoff rule.
analyse <- function(data, variables = c("x1", "x2", "x3"), components = 2L, standardise = TRUE, missing = "error") {
    spec <- infercalcR::method_spec("pca_components", alternative = "two.sided", conf.level = 0.95)
    infercalcR::run_inference(spec, data = data, variables = variables, components = components, standardise = standardise,
        missing = missing)$payload$native
}
example_data<-data.frame(id=1:8,x1=c(1,2,3,5,4,7,6,8),x2=c(2,1,4,3,6,5,8,7),x3=c(3,5,2,7,4,9,6,8))
example_result<-analyse(example_data,components=2L)
print(example_result)
