# Synthetic teaching example; replace only after checking the recorded design.
# No normality or variance pre-test changes the selected method.
# Missing or infinite measurements are rejected; choose and document a missingness strategy first.
# The mean difference is x minus y; a one-sided 'greater' hypothesis means x > y.
# alpha is not inferred from the observed result. The confidence level is explicit.
# Welch two-sample t-test: group variances are estimated separately.
analyse <- function(x, y, alternative = "two.sided", conf.level = 0.95) {
  stopifnot(is.numeric(x), is.numeric(y), length(x) >= 2L, length(y) >= 2L,
            all(is.finite(x)), all(is.finite(y)),
            length(alternative) == 1L, alternative %in% c("two.sided", "greater", "less"),
            is.numeric(conf.level), length(conf.level) == 1L,
            is.finite(conf.level), conf.level > 0, conf.level < 1)
  fit <- stats::t.test(x, y, paired = FALSE, var.equal = FALSE,
                       alternative = alternative, mu = 0, conf.level = conf.level)
  list(method = fit$method, n = c(x = length(x), y = length(y)),
       mean_difference = mean(x) - mean(y), confidence_interval = fit$conf.int,
       statistic = unname(fit$statistic), df = unname(fit$parameter),
       p_value = fit$p.value, alternative = fit$alternative)
}
example_x <- c(5.2, 6.1, 4.9, 7.3, 5.8, 6.5, 7.1, 4.7)
example_y <- c(3.8, 5.9, 4.2, 5.1, 6.7, 4.8, 3.9, 5.3, 6.0)
example_result <- analyse(example_x, example_y, alternative = alternative,
                          conf.level = confidence_level)
print(example_result)
plot_example <- function() {
  old <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(old))
  stats::qqnorm(example_x - mean(example_x), main = "Group x: centered measurements")
  stats::qqline(example_x - mean(example_x))
  stats::qqnorm(example_y - mean(example_y), main = "Group y: centered measurements")
  stats::qqline(example_y - mean(example_y))
}
# Inspect tail behaviour and outliers in context; these plots do not verify independence.
plot_example()
