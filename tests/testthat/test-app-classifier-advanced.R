test_that("all declared constraints and capacity units reach the operating rule", {
 i <- list(example="uploaded",selection="constraint",optimise="tpr",tie_break="most_alerts",
   tpr_min="0.2",tnr_min="0.3",ppv_min="0.4",npv_min="0.5",fpr_max="0.6",fnr_max="0.7",fdr_max="0.8",alert_rate_max="0.9")
 r <- .ig_classifier_input_rule(i,"en",NULL)
 expect_setequal(names(r$constraints),c("tpr_min","tnr_min","ppv_min","npv_min","fpr_max","fnr_max","fdr_max","alert_rate_max"))
 expect_equal(r$constraints$npv_min,.5);expect_equal(r$constraints$fdr_max,.8)
 expect_identical(r$tie_break,"most_alerts")
 i$selection <- "capacity";i$capacity_kind <- "max_alert_rate";i$max_alert_rate <- "0,25"
 for(nm in .ig_cls_constraint_names) i[[nm]] <- ""
 r <- .ig_classifier_input_rule(i,"de",NULL)
 expect_equal(r$capacity,list(max_alert_rate=.25))
 i$max_alert_rate <- "1,1"
 expect_error(.ig_classifier_input_rule(i,"de",NULL),class="inferguidR_invalid_input")
 i$selection <- "cost_utility";i$cost_fn <- "5";i$cost_fp <- "2";i$cost_tp <- "-1";i$cost_tn <- "-0,5"
 r <- .ig_classifier_input_rule(i,"de",NULL)
 expect_equal(r$costs,list(fn=5,fp=2,tn=-.5,tp=-1))
})

test_that("per-model direction, comparator and metadata stay with the score column", {
 d <- data.frame(y=c(1,1,0,0),a=c(.5,.8,.5,.1),b=c(.5,.2,.5,.9))
 i <- list(example="uploaded",selection="manual",cutoff="0.5",truth_column="y",score_columns=c("a","b"),positive="1")
 set <- function(column,field,value) i[[.ig_classifier_model_key(column,field)]] <<- value
 set("a","name","Model A");set("a","version","v2");set("a","meaning","Raw test score")
 set("a","population","Synthetic");set("a","unit","unit");set("a","comparator",">")
 set("b","direction","lower_is_positive");set("b","comparator","<=")
 s <- .ig_classifier_custom_spec(i,"en",d);r <- evaluate_cutoff(s)
 expect_identical(s$models[[1]]$name,"Model A");expect_identical(s$models[[1]]$version,"v2")
 expect_equal(unlist(r$models[[1]]$operating_point[c("tp","fn","fp","tn")]),c(tp=1,fn=1,fp=0,tn=2))
 expect_equal(unlist(r$models[[2]]$operating_point[c("tp","fn","fp","tn")]),c(tp=2,fn=0,fp=1,tn=1))
 i$score_columns <- rev(i$score_columns)
 reordered <- .ig_classifier_custom_spec(i,"en",d)
 expect_identical(reordered$models[[2]]$name,"Model A")
 expect_identical(reordered$models[[1]]$direction,"lower_is_positive")
 set("a","direction","lower_is_positive")
 expect_error(.ig_classifier_custom_spec(i,"en",d),class="inferguidR_invalid_input")
})

test_that("weight and cluster semantics remain explicit", {
 d <- data.frame(y=c(1,1,0,0),score=c(.9,.8,.2,.1),w=c(2,1,3,1),cluster=c("a","a","b","b"))
 i <- list(example="uploaded",selection="manual",cutoff="0.5",truth_column="y",score_columns="score",positive="1",
   weight_column="w",cluster_column="cluster",sampling="case_control")
 s <- .ig_classifier_custom_spec(i,"en",d);r <- evaluate_cutoff(s)
 expect_equal(r$models[[1]]$operating_point$tp,3)
 expect_equal(r$models[[1]]$operating_point$tn,4)
 expect_identical(s$models[[1]]$weight_semantics,"frequency")
 expect_identical(s$prevalence$sampling,"case_control")
 d$w[1] <- .5
 expect_error(.ig_classifier_custom_spec(i,"en",d),class="inferguidR_unsupported")
 i$weight_column <- "score"
 expect_error(.ig_classifier_custom_spec(i,"en",d),class="inferguidR_invalid_input")
})

test_that("uncertainty controls retain explicit methods, seeds and locale", {
 i <- list(confidence_level="0,9",proportion_methods="wilson",auc_interval="none",bootstrap_replicates="100",
           bootstrap_seed="42",bootstrap_strata="none",bootstrap_metrics=c("tpr","roc_auc"))
 u <- .ig_classifier_uncertainty_input(i,"de")
 expect_equal(u$level,.9);expect_identical(u$proportion_methods,"wilson")
 expect_identical(u$auc_method,"none");expect_equal(u$bootstrap_replicates,100)
 expect_equal(u$seed,42);expect_identical(u$bootstrap_strata,"none")
 expect_identical(u$bootstrap_metrics,c("tpr","roc_auc"))
 i$bootstrap_seed <- ""
 expect_error(.ig_classifier_uncertainty_input(i,"de"),class="inferguidR_invalid_input")
 i$bootstrap_replicates <- "0";i$proportion_methods <- character()
 u <- .ig_classifier_uncertainty_input(i,"de")
 expect_null(u$proportion_methods);expect_null(u$bootstrap_metrics)
 i$confidence_level <- "95e-2"
 expect_equal(.ig_classifier_uncertainty_input(i,"de")$level,.95)
 expect_equal(.ig_classifier_uncertainty_input(i,"en")$level,.95)
})

test_that("uncertainty changes invalidate prior results and seeded reevaluation is reproducible", {
 shiny::testServer(.ig_classifier_workspace_server,args=list(language=shiny::reactive("en")), {
  session$setInputs(example="scored_synthetic",selection="manual",cutoff="0.3",scenario=FALSE,
    confidence_level="0.9",proportion_methods="wilson",auc_interval="none",bootstrap_replicates="100",
    bootstrap_seed="42",bootstrap_metrics="tpr",bootstrap_strata="class")
  session$setInputs(evaluate=1)
  old <- result()
  expect_s3_class(old,"classification_result")
  expect_equal(candidate()$uncertainty$level,.9)
  expect_identical(candidate()$spec_hash,old$spec_hash)
  session$setInputs(confidence_level="0.8")
  expect_false(identical(candidate()$spec_hash,result()$spec_hash))
  session$setInputs(evaluate=2)
  first <- result()
  session$setInputs(evaluate=3)
  expect_identical(result()$models[[1]]$intervals,first$models[[1]]$intervals)
  expect_identical(candidate()$spec_hash,result()$spec_hash)
 })
})

test_that("fixed counts check constraints without claiming a measured score cutoff", {
 i <- list(example="custom_counts",selection="manual",tp="4",fn="1",fp="2",tn="3",tpr_min="0.9")
 s <- .ig_classifier_custom_spec(i,"en")
 expect_identical(s$operating_rule$cutoff_source,"binary_test_convention")
 expect_false(evaluate_cutoff(s)$models[[1]]$selection$constraints_met)
 i$tpr_min <- "0.8"
 expect_true(evaluate_cutoff(.ig_classifier_custom_spec(i,"en"))$models[[1]]$selection$constraints_met)
})
