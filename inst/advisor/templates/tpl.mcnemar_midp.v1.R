# Synthetic teaching example: fix the hypothesis and method before analysis.
# Independent evaluation units are required; no survey or fractional weights.
check_settings <- function(alternative, conf.level) {
 stopifnot(is.character(alternative), length(alternative)==1L,
  alternative %in% c("two.sided","greater","less"),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
}
# n10: x positive, y negative; n01: x negative, y positive, in the same units.
# The first and second measurements and their positive labels must be fixed in advance.
check_counts <- function(counts) {
 stopifnot(is.numeric(counts),is.null(dim(counts)),length(counts)==4L,!anyNA(counts),all(is.finite(counts)),
  all(counts>=0),all(counts==floor(counts)),!anyDuplicated(names(counts)),
  setequal(names(counts),c("n00","n01","n10","n11")),sum(counts)>0)
 counts[c("n00","n01","n10","n11")]
}
analyse <- function(counts,alternative="two.sided",conf.level=.95) {
 check_settings(alternative,conf.level);counts <- check_counts(counts)
 b <- unname(counts["n10"]);c <- unname(counts["n01"]);n <- b+c
 if(n==0) return(list(status="not_applicable",reason="No discordant pairs; no conditional test information.",
  n_total=sum(counts),n_discordant=0,p_value=NULL,confidence_interval=NULL))
 fit <- stats::binom.test(b,n,p=.5,alternative=alternative,conf.level=conf.level)
 mass <- stats::dbinom(b,n,.5)
 lower <- stats::pbinom(b-1,n,.5)+.5*mass
 upper <- stats::pbinom(b,n,.5,lower.tail=FALSE)+.5*mass
 p <- switch(alternative,two.sided=min(1,2*min(lower,upper)),greater=upper,less=lower)
 list(status="ok",method="Conditional McNemar mid-p test",n_total=sum(counts),n_discordant=n,
  n_x_only=b,n_y_only=c,conditional_positive_probability=b/n,
  marginal_difference=unname((b-c)/sum(counts)),p_value=p,alternative=alternative,
  confidence_interval=as.numeric(fit$conf.int),confidence_level=conf.level,
  interval_method="Clopper-Pearson for P(x-only | discordant); not an inversion of the mid-p test",
  limitations="Mid-p does not guarantee exact type-I error control. Interval concerns discordance probability, not marginal difference.")
}
example_counts <- c(n00=20,n01=3,n10=9,n11=15)
example_result <- analyse(example_counts,alternative,confidence_level)
print(example_result)
