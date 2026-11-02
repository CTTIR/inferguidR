# Exploratory clustering with prespecified scaling, k and subsampling stability.
# Case identifiers establish deterministic tie ordering; labels have no intrinsic meaning.
analyse <- function(data, variables, k, standardize, independent_units = FALSE, alternative = "two.sided",
    conf.level = 0.95, complete_cases = FALSE, B = 50L, seed = 104729L) {
    spec <- infercalcR::method_spec("cluster_kmeans", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, variables = variables, k = k, standardize = standardize,
        independent_units = independent_units, complete_cases = complete_cases, B = B, seed = seed)$payload$native
}
example_data<-data.frame(id=sprintf('case%02d',1:18),x=c(-5,-4.6,-5.2,-4.8,-5.1,-4.9,0,.4,-.2,.2,-.1,.1,5,5.4,4.8,5.2,4.9,5.1),y=c(-4,-4.2,-3.8,-4.1,-3.9,-4,1,.8,1.2,.9,1.1,1,5,4.8,5.2,4.9,5.1,5))
example_result<-analyse(example_data,c('x','y'),3,TRUE,TRUE)
print(example_result[c('status','method','k','mean_silhouette','stability')])
