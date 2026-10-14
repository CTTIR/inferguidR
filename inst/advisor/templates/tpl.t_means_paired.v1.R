# Synthetic teaching example; replace only after checking the recorded design.
# No normality or variance pre-test changes the selected method.
# Missing or infinite measurements are rejected; choose and document a missingness strategy first.
# The mean difference is x minus y; a one-sided 'greater' hypothesis means x > y.
# alpha is not inferred from the observed result. The confidence level is explicit.
# Paired t-test: one wide row per confirmed ID; x and y are matched occasions.
# Pair IDs are supplied explicitly, never guessed from row order in two separate tables.
analyse <- function(data, alternative = "two.sided", conf.level = 0.95) {
  stopifnot(is.data.frame(data), all(c("id", "x", "y") %in% names(data)),
            nrow(data) >= 2L, !anyNA(data$id), !anyDuplicated(data$id),
            is.numeric(data$x), is.numeric(data$y),
            all(is.finite(data$x)), all(is.finite(data$y)),
            length(alternative) == 1L, alternative %in% c("two.sided", "greater", "less"),
            is.numeric(conf.level), length(conf.level) == 1L,
            is.finite(conf.level), conf.level > 0, conf.level < 1)
  fit <- stats::t.test(data$x, data$y, paired = TRUE, alternative = alternative,
                       mu = 0, conf.level = conf.level)
  list(method = fit$method, n_pairs = nrow(data),
       mean_difference = mean(data$x - data$y), confidence_interval = fit$conf.int,
       statistic = unname(fit$statistic), df = unname(fit$parameter),
       p_value = fit$p.value, alternative = fit$alternative)
}
example_data <- data.frame(id = paste0("participant_", seq_len(8)),
 x = c(5.2, 6.1, 4.9, 7.3, 5.8, 6.5, 7.1, 4.7),
 y = c(3.8, 5.9, 4.2, 5.1, 6.7, 4.8, 3.9, 5.3))
example_result <- analyse(example_data, alternative = alternative,
                          conf.level = confidence_level)
print(example_result)
plot_example <- function() {
  differences <- example_data$x - example_data$y
  stats::qqnorm(differences, main = "Within-pair differences (x minus y)")
  stats::qqline(differences)
}
# Normality concerns within-pair differences, not the two raw occasion distributions.
# This plot does not verify independence between participants.
plot_example()
