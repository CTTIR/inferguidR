# Ordered-group rank trend under exchangeability of independent observations.
# Choose the group order and alternative before inspecting outcomes.
analyse <- function(data, group_order, order_prespecified, independent_units, seed, alternative = "two.sided",
    B = 1999L, conf.level = 0.95, mc.level = 0.99, complete_cases = FALSE) {
    spec <- infercalcR::method_spec("jonckheere_terpstra", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, group_order = group_order, order_prespecified = order_prespecified,
        independent_units = independent_units, seed = seed, B = B, mc.level = mc.level, complete_cases = complete_cases)$payload$native
}
example_data<-data.frame(id=paste0("unit",1:9),group=rep(c("low","middle","high"),each=3),y=c(1,3,2,2,4,3,5,4,6))
example_result<-analyse(example_data,c("low","middle","high"),TRUE,TRUE,seed=104729L,
 alternative=alternative,B=1999L,conf.level=confidence_level)
print(example_result)
