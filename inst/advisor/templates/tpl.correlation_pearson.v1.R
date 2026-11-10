# Synthetic teaching example. The selected method and hypothesis are fixed before analysis.
# Reject missing and infinite inputs; document an explicit missingness strategy first.
check_settings <- function(alternative, conf.level) {
 stopifnot(length(alternative)==1L, alternative %in% c("two.sided","greater","less"),
           is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Pearson correlation assumes a linear association; its t test assumes bivariate normal independent pairs.
# x and y must already be matched by a confirmed unique unit ID; row matching is never guessed.
# Greater means positive association; less means negative association.
analyse <- function(data, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("correlation_pearson", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data)$payload$native
}
example_data <- data.frame(id=paste0("unit_",1:8),x=c(1,2,3,4,5,6,7,8),y=c(2,1,4,3,7,5,8,6))
example_result <- analyse(example_data,alternative,confidence_level)
print(example_result)
graphics::plot(example_data$x,example_data$y,xlab="x",ylab="y",main="Matched observations")
