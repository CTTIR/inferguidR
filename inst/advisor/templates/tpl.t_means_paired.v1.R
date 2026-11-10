# Synthetic teaching example; replace only after checking the recorded design.
# No normality or variance pre-test changes the selected method.
# Missing or infinite measurements are rejected; choose and document a missingness strategy first.
# The mean difference is x minus y; a one-sided 'greater' hypothesis means x > y.
# alpha is not inferred from the observed result. The confidence level is explicit.
# Paired t-test: one wide row per confirmed ID; x and y are matched occasions.
# Pair IDs are supplied explicitly, never guessed from row order in two separate tables.
analyse <- function(data, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("t_means_paired", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data)$payload$native
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
