.report_rehash <- function(x) {
 x <- unclass(x);x$snapshot_hash <- NULL;x$snapshot_hash <- .ig_hash_value(x);x
}

test_that("rehashed reports reject malformed and inconsistent nested content", {
 r <- boschloo_test(1,1,0,1,"greater")
 x <- create_report_snapshot(analysis_envelope(r$spec,r))
 mutations <- list(
  function(y) {y$results[[1]]$p_value <- -.1;y},
  function(y) {y$results[[1]]$p_upper <- 2;y},
  function(y) {y$specifications[[1]]$analysis_hash <- paste(rep("0",64),collapse="");y},
  function(y) {y$metrics[[1]]$p_value <- .8;y},
  function(y) {y$specifications <- rep(y$specifications,2);y$results <- rep(y$results,2);y},
  function(y) {y$provenance <- "not a record";y},
  function(y) {y$results[[1]]$status <- NULL;y},
  function(y) {y$results[[1]]$unrecognised_field <- TRUE;y},
  function(y) {y$visualisations[[1]]$data$event <- list(-1,1);y},
  function(y) {y$comparisons[[1]] <- list(p_value=.5);y},
  function(y) {y$completeness <- "incomplete";y},
  function(y) {y$report_kind <- "model_cutoff";y},
  function(y) {y$privacy$score_curve_tables_included <- TRUE;y},
  function(y) {y$visualisations[[1]]$analysis <- .5;y},
  function(y) {y$display_timezone <- "invented";y})
 for(mutate in mutations) {
  y <- .report_rehash(mutate(unclass(x)))
  expect_error(.ig_report_document_from_list(y),class="inferguidR_invalid_input")
  expect_error(ig_from_json(.ig_canonical_json(y)),class="inferguidR_invalid_input")
 }
})

test_that("domain checks remain effective if analysis hashes are also recomputed", {
 r <- boschloo_test(1,1,0,1,"greater")
 x <- create_report_snapshot(analysis_envelope(r$spec,r))
 y <- unclass(x);y$results[[1]]$p_value <- -.1;y$metrics[[1]]$p_value <- -.1
 reconstructed <- y$results[[1]];reconstructed$spec <- r$spec;class(reconstructed) <- class(r)
 y$specifications[[1]]$analysis_hash <- analysis_envelope(r$spec,reconstructed)$analysis_hash
 y$analysis_hash <- .ig_hash_value(list(report_kind=y$report_kind,
   envelopes=I(vapply(y$specifications,`[[`,"","analysis_hash")),recommendation=NULL))
 expect_error(.ig_report_document_from_list(.report_rehash(y)),"p-value")
})

test_that("valid aggregate classifier reports and missing settings remain portable", {
 for(name in c("scored_synthetic","counts_fixed","binormal_equal_variance","binormal_unequal_variance")) {
  s <- classification_examples(name)$spec;r <- evaluate_cutoff(s)
  x <- create_report_snapshot(analysis_envelope(s,r),report_kind="model_cutoff")
  expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
 }
 s <- analysis_spec("t_means_one_sample","post_hoc",list(d=.5,alpha=.05,n_total=30),list(tails="two"))
 x <- create_report_snapshot(analysis_envelope(s),report_kind="incomplete_settings")
 expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
})

test_that("portable report schema checks nested records and probability domains", {
 skip_if_not_installed("jsonvalidate")
 r <- boschloo_test(1,1,0,1,"greater")
 observed <- analysis_envelope(r$spec,r)
 s <- analysis_spec("t_means_one_sample","post_hoc",list(d=.5,alpha=.05,n_total=30),list(tails="two"))
 power <- analysis_envelope(s,calculate_power(s))
 for(e in list(observed,power)) {
  x <- create_report_snapshot(e)
  expect_true(.ig_validate_schema(ig_to_json(x),"report_document")$valid)
  y <- unclass(x);y$results[[1]]$unexpected <- TRUE
  expect_false(.ig_validate_schema(.ig_canonical_json(y),"report_document")$valid)
 }
 for(name in c("scored_synthetic","counts_fixed","binormal_equal_variance")) {
  s <- classification_examples(name)$spec;x <- create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
  expect_true(.ig_validate_schema(ig_to_json(x),"report_document")$valid,info=name)
 }
 x <- create_report_snapshot(analysis_envelope(power$spec),report_kind="incomplete_settings")
 expect_true(.ig_validate_schema(ig_to_json(x),"report_document")$valid)
 x <- unclass(create_report_snapshot(observed));x$results[[1]]$p_value <- -.1
 expect_false(.ig_validate_schema(.ig_canonical_json(x),"report_document")$valid)
})
