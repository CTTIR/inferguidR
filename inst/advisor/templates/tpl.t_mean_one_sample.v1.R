# Synthetic teaching example. The selected method and hypothesis are fixed before analysis.
# Reject missing and infinite inputs; document an explicit missingness strategy first.
check_settings <- function(alternative, conf.level) {
 stopifnot(length(alternative)==1L, alternative %in% c("two.sided","greater","less"),
           is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Estimand: population mean minus the prespecified reference mean.
# Independence is a design assumption; a QQ plot cannot establish it.
analyse <- function(x, mu0, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("t_mean_one_sample", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, mu0 = mu0)$payload$native
}
example_x <- c(5.2,6.1,4.9,7.3,5.8,6.5,7.1,4.7)
example_reference <- 5
example_result <- analyse(example_x,example_reference,alternative,confidence_level)
print(example_result)
stats::qqnorm(example_x,main="Measurements: normal-model diagnostic")
stats::qqline(example_x)
