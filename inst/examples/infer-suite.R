# Synthetic, explicit-design examples for matching installed suite candidates.

spec <- infercalcR::method_spec("t_means_two_independent_welch")
result <- infercalcR::run_inference(spec,
  x = c(5.2, 6.1, 4.9, 7.3, 5.8, 6.5, 7.1, 4.7),
  y = c(3.8, 5.9, 4.2, 5.1, 6.7, 4.8, 3.9, 5.3, 6.0))
infercalcR::result_data(result)
infereportR::report_table(result)
infereportR::report_text(result)

plot_data <- data.frame(group = rep(c("x", "y"), c(8, 9)),
  value = c(5.2, 6.1, 4.9, 7.3, 5.8, 6.5, 7.1, 4.7,
            3.8, 5.9, 4.2, 5.1, 6.7, 4.8, 3.9, 5.3, 6.0))
positions <- data.frame(row_id = infercalcR::result_data(result)$row_id[1],
  xmin = 1, xmax = 2, y_position = 8)
ggplot2::ggplot(plot_data, ggplot2::aes(group, value)) +
  ggplot2::geom_boxplot() +
  infereportR::stat_inference(result, positions, rows = positions$row_id)

control <- infercalcR::boschloo_control()
planned <- inferguidR::boschloo_power(15, 15, .2, .6,
  alpha = .05, variant = "minlike", control = control)
observed <- infercalcR::boschloo_test(3, 15, 9, 15,
  alpha = .05, variant = "minlike", groups = c("control", "treatment"),
  event = "response", control = control)
infereportR::report_table(infercalcR::as_inference(observed))

d <- data.frame(id = paste0("unit", 1:6),
  truth = c("yes", "no", "yes", "no", "yes", "no"),
  score = c(.9, .2, .7, .4, .6, .3))
spec <- infercalcR::classification_spec(d, truth = "truth",
  positive_label = "yes", scores = "score", id = "id", cutoff = .5,
  score_scale = "probability", evaluation_split = "held_out_test")
result <- infercalcR::evaluate_cutoff(spec)
infereportR::report_table(infercalcR::as_inference(result))
