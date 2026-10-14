test_that("mixed model scales validate separately and enable only applicable metrics", {
 d <- data.frame(y=c(1,0,1,0),raw=c(2,-1,4,-3),risk=c(.8,.2,.7,.1))
 s <- classification_spec(d,truth="y",positive_label=1,scores=c(raw_model="raw",risk_model="risk"),cutoff=.5,
   score_scale=c(risk_model="probability",raw_model="raw_score"))
 expect_identical(vapply(s$models,`[[`,"","score_scale"),c("raw_score","probability"))
 expect_match(s$models[[1]]$score_meaning,"raw scale",fixed=TRUE)
 expect_match(s$models[[2]]$score_meaning,"predicted probability",fixed=TRUE)
 r <- evaluate_cutoff(s)
 p <- r$models[[2]]$summaries
 expect_equal(p$value[p$metric_id=="brier_score"],mean((d$risk-d$y)^2),tolerance=1e-14)
 raw <- r$models[[1]]$summaries
 expect_identical(raw$status[raw$metric_id=="brier_score"],"not_applicable")
 expect_equal(r$models[[1]]$operating_point$tp,2)
 expect_equal(r$models[[1]]$operating_point$tn,2)
 back <- ig_from_json(ig_to_json(analysis_envelope(s)))$spec
 attached <- attach_classification_data(back,d)
 expect_identical(attached$spec_hash,s$spec_hash)
 expect_identical(evaluate_cutoff(attached)$models,r$models)
 d$risk[1] <- 2
 expect_error(attach_classification_data(back,d),"outside",class="inferguidR_invalid_input")
})

test_that("score scale maps cannot silently omit, reorder or duplicate models", {
 d <- data.frame(y=c(1,0),a=c(.8,.2),b=c(.7,.3))
 make <- function(scale) classification_spec(d,truth="y",positive_label=1,scores=c(A="a",B="b"),cutoff=.5,score_scale=scale)
 for(scale in list(c("raw_score","probability"),c(A="probability"),c(A="raw_score",A="probability"),
   c(A="probability",C="raw_score"),c(A=NA_character_,B="probability"),c(A="unknown",B="probability")))
  expect_error(make(scale),class="inferguidR_invalid_input")
 expect_identical(make("probability")$spec_hash,make(c(B="probability",A="probability"))$spec_hash)
 s <- classification_spec(d,truth="y",positive_label=1,scores=c(A="a",B="b"),cutoff=.5)
 expect_identical(s$spec_hash,make("raw_score")$spec_hash)
 d$a <- matrix(c(.8,.2,.9,.1),2)
 expect_error(classification_spec(d,truth="y",positive_label=1,scores="a",cutoff=.5),"numeric vector",class="inferguidR_invalid_input")
})

test_that("reattachment preserves the declared negative outcome label", {
 d <- data.frame(y=c("event","control","event","control"),score=c(.8,.2,.7,.1))
 s <- classification_spec(d,truth="y",positive_label="event",negative_label="control",scores="score",cutoff=.5)
 back <- ig_from_json(ig_to_json(analysis_envelope(s)))$spec
 d$y[d$y=="control"] <- "other"
 expect_error(attach_classification_data(back,d),"declared labels",class="inferguidR_invalid_input")
})

test_that("per-model app scale choices override only the selected column", {
 d <- data.frame(y=c(1,0),raw=c(2,-1),risk=c(.8,.2))
 i <- list(example="uploaded",selection="manual",cutoff="0.5",truth_column="y",positive="1",
   score_columns=c("raw","risk"),score_scale="raw_score")
 i[[.ig_classifier_model_key("risk","scale")]] <- "probability"
 s <- .ig_classifier_custom_spec(i,"en",d)
 expect_identical(s$models[[1]]$score_scale,"raw_score")
 expect_identical(s$models[[2]]$score_scale,"probability")
 i[[.ig_classifier_model_key("raw","scale")]] <- "probability"
 expect_error(.ig_classifier_custom_spec(i,"en",d),"outside",class="inferguidR_invalid_input")
})
