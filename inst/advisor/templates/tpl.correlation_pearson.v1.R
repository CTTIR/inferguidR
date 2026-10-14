# Synthetic teaching example. The selected method and hypothesis are fixed before analysis.
# Reject missing and infinite inputs; document an explicit missingness strategy first.
check_settings <- function(alternative, conf.level) {
 stopifnot(length(alternative)==1L, alternative %in% c("two.sided","greater","less"),
           is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# Pearson correlation assumes a linear association; its t test assumes bivariate normal independent pairs.
# x and y must already be matched by a confirmed unique unit ID; row matching is never guessed.
# Greater means positive association; less means negative association.
analyse <- function(data,alternative="two.sided",conf.level=.95) {
 check_settings(alternative,conf.level)
 stopifnot(is.data.frame(data),all(c("id","x","y") %in% names(data)),nrow(data)>=4,
           !anyNA(data$id),!anyDuplicated(data$id),is.numeric(data$x),is.numeric(data$y),
           all(is.finite(data$x)),all(is.finite(data$y)),stats::sd(data$x)>0,stats::sd(data$y)>0)
 fit <- stats::cor.test(data$x,data$y,method="pearson",alternative=alternative,conf.level=conf.level)
 list(method=fit$method,n=nrow(data),estimate=unname(fit$estimate),statistic=unname(fit$statistic),
      confidence_interval=if(is.null(fit$conf.int)) NULL else fit$conf.int,
      interval_status=if(is.null(fit$conf.int)) "not_provided_by_method" else "available",
      p_value=fit$p.value,alternative=alternative)
}
example_data <- data.frame(id=paste0("unit_",1:8),x=c(1,2,3,4,5,6,7,8),y=c(2,1,4,3,7,5,8,6))
example_result <- analyse(example_data,alternative,confidence_level)
print(example_result)
graphics::plot(example_data$x,example_data$y,xlab="x",ylab="y",main="Matched observations")
