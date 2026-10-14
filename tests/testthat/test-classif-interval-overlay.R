interval_example <- function(cluster=FALSE) {
 d <- data.frame(y=c(rep(1,5),rep(0,5)),score=c(.9,.8,.7,.6,.2,.8,.7,.3,.2,.1),
                 cluster=rep(letters[1:5],2))
 evaluate_cutoff(classification_spec(d,truth="y",positive_label=1,scores="score",cutoff=.5,
   cluster=if(cluster) "cluster" else NULL))
}

test_that("empty uncertainty and comparison tables state why they are unavailable", {
 s <- classification_spec(counts=c(tp=4,fn=1,fp=2,tn=3),
   uncertainty=classification_uncertainty(proportion_methods=character(),auc_method="none"))
 r <- evaluate_cutoff(s)
 d <- .ig_classifier_interval_table(r,"test","observed","en")
 expect_identical(d$status,"not_computed")
 expect_match(d$reason,"No intervals were computed")
 d <- .ig_classifier_comparison_table(r,"de")
 expect_identical(d$inference_status,"not_applicable")
 expect_match(d$reason,"Mindestens zwei")
})

test_that("ROC and PR point intervals use the correct marginal denominators", {
 r <- interval_example()
 for(panel in c("roc","pr")) for(method in c("wilson","clopper_pearson")) {
  p <- plot_classification(r,panel,interval_method=method)
  d <- attr(p,"inferguidR_plot_data")$data$intervals
  expect_identical(d$metric_id,if(panel=="roc") c("fpr","tpr") else c("tpr","ppv"))
  k <- if(panel=="roc") c(2,4) else c(4,4)
  n <- if(panel=="roc") c(5,5) else c(5,6)
  if(method=="wilson") {
   z<-qnorm(.975);ph<-k/n;centre<-(ph+z^2/(2*n))/(1+z^2/n)
   half<-z*sqrt(ph*(1-ph)/n+z^2/(4*n^2))/(1+z^2/n)
   lo<-centre-half;hi<-centre+half
  } else {lo<-qbeta(.025,k,n-k+1);hi<-qbeta(.975,k+1,n-k)}
  expect_equal(d$lower,lo,tolerance=1e-12);expect_equal(d$upper,hi,tolerance=1e-12)
  expect_identical(d$status,c("ok","ok"))
  expect_equal(d$x[1],lo[1]);expect_equal(d$xend[1],hi[1])
  expect_equal(d$y[2],lo[2]);expect_equal(d$yend[2],hi[2])
  built<-ggplot2::ggplot_build(p)
  expect_equal(nrow(tail(built$data,1)[[1]]),2)
  expect_match(gsub("[[:space:]]+"," ",p$labels$caption),"neither a joint confidence")
 }
})

test_that("interval display uses German decimals without changing numeric records", {
 d <- data.frame(metric_id="tpr",lower=.123456,upper=.987654,level=.95)
 shown <- .ig_classifier_display_table(d,"de")
 expect_identical(shown$Untergrenze,"0,123456")
 expect_identical(shown$Konfidenzniveau,"0,95")
 expect_equal(d$lower,.123456)
})

test_that("scenario, clustered and uncomputed intervals are never invented", {
 r <- interval_example()
 s <- update(classification_evaluated_spec(r),prevalence=.1)
 scenario <- evaluate_cutoff(s)
 p <- plot_classification(scenario,"pr",basis="scenario",interval_method="wilson")
 d <- attr(p,"inferguidR_plot_data")$data$intervals
 expect_true(all(d$status=="not_applicable"));expect_true(all(is.na(d$lower)))
 tab <- .ig_classifier_interval_table(scenario,"score","scenario","de")
 expect_true(all(tab$status=="not_applicable"));expect_true(all(grepl("Szenariowerte",tab$reason)))
 p <- plot_classification(interval_example(TRUE),"roc",interval_method="wilson")
 d <- attr(p,"inferguidR_plot_data")$data$intervals
 expect_true(all(d$status=="unsupported"));expect_true(all(is.na(d$lower)))
 p <- plot_classification(r,"roc",interval_method="bootstrap_percentile")
 expect_true(all(attr(p,"inferguidR_plot_data")$data$intervals$status=="not_computed"))
 expect_error(plot_classification(r,"roc",interval_method="unknown"),class="inferguidR_invalid_input")
 expect_error(plot_classification(r,"mosaic",interval_method="wilson"),class="inferguidR_invalid_input")
})

test_that("undefined current PR points retain reasons without interval segments", {
 r <- interval_example();r<-evaluate_cutoff(update(classification_evaluated_spec(r),cutoff=10))
 p <- plot_classification(r,"pr",interval_method="wilson")
 d <- attr(p,"inferguidR_plot_data")$data$intervals
 expect_true(all(d$status=="not_applicable"));expect_true(all(is.na(d$lower)))
 expect_true(all(grepl("No drawable",d$reason)))
})

test_that("the workbench exposes uncertainty and paired-comparison tables", {
 shiny::testServer(.ig_classifier_workspace_server,args=list(language=shiny::reactive("en")), {
  session$setInputs(example="scored_synthetic",selection="manual",cutoff="0.3",scenario=FALSE,
    curve_panel="roc",interval_overlay="wilson")
  session$setInputs(evaluate=1)
  expect_match(output$intervals,"wilson")
  expect_match(output$comparisons,"delong_paired")
  expect_match(output$comparisons,"model_a")
  expect_true(all(curve_data()$data$intervals$status=="ok"))
  session$setInputs(interval_overlay="bootstrap_percentile")
  expect_true(all(curve_data()$data$intervals$status=="not_computed"))
  expect_true(all(grepl("not computed",curve_data()$data$intervals$reason)))
  session$setInputs(scenario=TRUE,prevalence="0.1");session$setInputs(evaluate=2,basis="scenario")
  expect_match(output$intervals,"Scenario values are expected")
 })
})
