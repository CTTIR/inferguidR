# Synthetic teaching example. The selected method and hypothesis are fixed before analysis.
# Reject missing and infinite inputs; document an explicit missingness strategy first.
check_settings <- function(alternative, conf.level) {
 stopifnot(length(alternative)==1L, alternative %in% c("two.sided","greater","less"),
           is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Estimand: population mean minus the prespecified reference mean.
# Independence is a design assumption; a QQ plot cannot establish it.
analyse <- function(x, mu0, alternative="two.sided", conf.level=.95) {
 check_settings(alternative,conf.level)
 stopifnot(is.numeric(x),length(x)>=2L,all(is.finite(x)),is.numeric(mu0),length(mu0)==1L,is.finite(mu0))
 fit <- stats::t.test(x,mu=mu0,alternative=alternative,conf.level=conf.level)
 list(method=fit$method,n=length(x),mean=mean(x),reference=mu0,mean_difference=mean(x)-mu0,
      confidence_interval=fit$conf.int-mu0,statistic=unname(fit$statistic),df=unname(fit$parameter),p_value=fit$p.value)
}
example_x <- c(5.2,6.1,4.9,7.3,5.8,6.5,7.1,4.7)
example_reference <- 5
example_result <- analyse(example_x,example_reference,alternative,confidence_level)
print(example_result)
stats::qqnorm(example_x,main="Measurements: normal-model diagnostic")
stats::qqline(example_x)
