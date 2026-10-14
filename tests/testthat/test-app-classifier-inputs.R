test_that("CSV readers preserve names and reject malformed structures", {
 p <- tempfile(fileext=".csv");on.exit(unlink(p),add=TRUE)
 writeLines(c("outcome;score;id","0;0,2;a","1;0,8;b"),p)
 d <- .ig_classifier_read_csv(p,";",",")
 expect_equal(d$score,c(.2,.8))
 expect_error(.ig_classifier_read_csv(p,",",","),class="inferguidR_invalid_input")
 writeLines(c("score,score","1,2"),p)
 expect_error(.ig_classifier_read_csv(p),class="inferguidR_invalid_input")
})

test_that("uploaded models use explicit outcome, directions and constraints", {
 d <- data.frame(outcome=c(0,0,1,1),a=c(.1,.4,.6,.9),b=c(.3,.7,.4,.8))
 i <- list(example="uploaded",cutoff="0.5",selection="manual",truth_column="outcome",score_columns=c("a","b"),positive="1")
 s <- .ig_classifier_custom_spec(i,"en",d)
 r <- evaluate_cutoff(s)
 expect_length(r$models,2)
 expect_identical(vapply(r$models, `[[`, "", "model_id"), c("model_1", "model_2"))
 expect_equal(r$models[[1]]$operating_point$tp,2)
 expect_equal(r$models[[2]]$operating_point$tp,1)
 i$selection <- "constraint";i$tpr_min <- "1";i$fpr_max <- "0.5";i$optimise <- "ppv"
 r <- evaluate_cutoff(.ig_classifier_custom_spec(i,"en",d))
 for (m in r$models) {
   expect_identical(m$selection$status, "ok")
   expect_equal(m$operating_point$tp, 2)
   expect_lte(m$operating_point$fp, 1)
 }
 i$score_columns <- "outcome"
 expect_error(.ig_classifier_custom_spec(i,"en",d),class="inferguidR_invalid_input")
 i$score_columns <- "a";i$selection_data <- "predeclared_rule"
 expect_error(.ig_classifier_custom_spec(i,"en",d),class="inferguidR_invalid_input")
})

test_that("counts and theory retain their different estimation bases", {
 i <- list(example="custom_counts",tp="40",fn="10",fp="20",tn="130")
 s <- .ig_classifier_custom_spec(i,"en")
 expect_equal(s$counts$tp,40)
 i$selection <- "max_f1"
 expect_error(.ig_classifier_custom_spec(i,"en"),class="inferguidR_invalid_input")
 i <- list(example="custom_theory",cutoff="1",selection="manual",auc="0,8",sigma="1",theory_n="1000",theory_prevalence="0,1")
 s <- .ig_classifier_custom_spec(i,"de")
 expect_equal(s$theory$auc,.8)
 expect_identical(s$operating_rule$selection_data,"predeclared_rule")
 expect_identical(evaluate_cutoff(s)$models[[1]]$count_kind,"expected")
})
