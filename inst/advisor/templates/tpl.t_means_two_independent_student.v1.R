# Synthetic teaching example; replace only after checking the recorded design.
# No normality or variance pre-test changes the selected method.
# Missing or infinite measurements are rejected; choose and document a missingness strategy first.
# The mean difference is x minus y; a one-sided 'greater' hypothesis means x > y.
# alpha is not inferred from the observed result. The confidence level is explicit.
# Pooled two-sample t-test: equal population variances must be justified in advance.
analyse <- function(x, y, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("t_means_two_independent_student", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, y = y)$payload$native
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
