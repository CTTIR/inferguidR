test_that("metric families preserve coordinates and disclose omitted values", {
 skip_if_not_installed("ggplot2")
 r <- evaluate_cutoff(classification_examples("scored_synthetic")$spec)
 for (fam in names(.ig_cls_families)) {
  p <- plot_classification(r, "metric_vs_cutoff", metric_family = fam)
  pd <- classification_plot_data(r, "metric_vs_cutoff")
  expect_identical(attr(p, "inferguidR_plot_data"), pd)
  expect_true(all(p$data$status == "ok" & is.finite(p$data$value)))
  expect_silent(ggplot2::ggplot_build(p))
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomPath"), logical(1))))
 }
 p <- plot_classification(r,"metric_vs_cutoff",metric_family="summary")
 expect_equal(p$coordinates$limits$y,c(-1,1))
 p <- plot_classification(r,"metric_vs_cutoff",metric_family="likelihood_ratios")
 expect_true(all(p$data$value>0))
 expect_match(p$labels$caption,"data table")
 expect_error(plot_classification(r,"metric_vs_cutoff",metric_family="absent"),"supported")
})

test_that("log metric pairs omit rather than clamp undefined and infinite values", {
 skip_if_not_installed("ggplot2")
 r <- evaluate_cutoff(classification_examples("scored_synthetic")$spec)
 p <- plot_classification(r,"metric_pair",metric_x="lr_pos",metric_y="ppv",log_x=TRUE,log_y=TRUE)
 pd <- attr(p,"inferguidR_plot_data")
 keep <- pd$data$points$x_status=="ok" & pd$data$points$y_status=="ok" &
   is.finite(pd$data$points$x) & is.finite(pd$data$points$y) & pd$data$points$x>0 & pd$data$points$y>0
 expect_true(any(!keep))
 expect_identical(p$data$set_index,pd$data$points$set_index[keep])
 expect_equal(p$data$px,pd$data$points$x[keep])
 expect_equal(p$data$py,pd$data$points$y[keep])
 expect_match(p$labels$caption,as.character(sum(!keep)),fixed=TRUE)
 expect_silent(ggplot2::ggplot_build(p))
 expect_error(plot_classification(r,"metric_pair",log_x=NA),"TRUE or FALSE")
 expect_error(plot_classification(r,"metric_pair",metric_x="invented"),"metric ids")
})

test_that("dependency view uses the dictionary and preserves statuses", {
 skip_if_not_installed("ggplot2")
 r <- evaluate_cutoff(classification_examples("counts_fixed")$spec)
 p <- plot_classification(r,"metric_network",metric_focus="ppv")
 pd <- classification_plot_data(r,"metric_network")
 expect_identical(attr(p,"inferguidR_plot_data"),pd)
 parents <- unique(pd$data$edges$from[pd$data$edges$to=="ppv"])
 expect_setequal(p$data$id,c(parents,"ppv"))
 expect_true(all(p$data$value==pd$data$nodes$value[match(p$data$id,pd$data$nodes$id)]))
 expect_silent(ggplot2::ggplot_build(p))
 expect_match(p$labels$caption,"not[[:space:]]+causal")
 expect_error(plot_classification(r,"metric_network",metric_focus="made_up"),"operating-point metric")
 leaf <- plot_classification(r,"metric_network",metric_focus="tp")
 expect_equal(nrow(leaf$data),1)
 expect_silent(ggplot2::ggplot_build(leaf))
})


test_that("every metric dependency view has valid dictionary nodes", {
 skip_if_not_installed("ggplot2")
 r <- evaluate_cutoff(classification_examples("counts_fixed")$spec)
 for (id in .ig_cls_operating_ids()) {
  p <- plot_classification(r,"metric_network",metric_focus=id)
  expect_false(anyNA(p$data$id),info=id)
  expect_false(anyNA(p$data$status),info=id)
  expect_silent(ggplot2::ggplot_build(p))
 }
})

test_that("entirely nonpositive metric pairs do not invent log coordinates", {
 skip_if_not_installed("ggplot2")
 d <- data.frame(y=c(1,1,0,0),s=c(0,0,1,1))
 r <- evaluate_cutoff(classification_spec(d,truth="y",positive_label=1,scores="s",cutoff=.5,
   direction="higher_is_positive",comparator=">="))
 expect_error(plot_classification(r,"metric_pair",metric_x="mcc",metric_y="youden_j",log_x=TRUE),
              "No values are drawable")
 expect_silent(ggplot2::ggplot_build(plot_classification(r,"metric_pair",metric_x="mcc",metric_y="youden_j")))
})
