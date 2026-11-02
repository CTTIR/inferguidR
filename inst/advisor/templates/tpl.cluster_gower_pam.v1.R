# Explicit mixed measurement types, weights and k; no silent type downgrading.
analyse <- function(data, variables, types, weights, k, independent_units = FALSE, alternative = "two.sided",
    conf.level = 0.95, missing = "error", B = 50L, seed = 104729L) {
    spec <- infercalcR::method_spec("cluster_gower_pam", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, variables = variables, types = types, weights = weights,
        k = k, independent_units = independent_units, missing = missing, B = B, seed = seed)$payload$native
}
example_data<-data.frame(id=sprintf('case%02d',1:12),amount=c(1,2,1.5,2.5,1.8,2.2,8,9,8.5,9.5,8.8,9.2),kind=factor(rep(c('a','b'),each=6)),grade=ordered(rep(c('low','high'),each=6),levels=c('low','middle','high')),flag=rep(c(0,1),each=6))
example_types<-c(amount='continuous',kind='nominal',grade='ordinal',flag='asymmetric_binary')
example_result<-analyse(example_data,names(example_types),example_types,stats::setNames(rep(1,4),names(example_types)),2,TRUE)
print(example_result[c('status','method','medoid_ids','mean_silhouette','stability')])
