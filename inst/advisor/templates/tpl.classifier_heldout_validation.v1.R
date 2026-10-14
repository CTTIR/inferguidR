# Synthetic independent-unit example. Explicit score direction is never chosen to improve AUC.
# Identifiers must have the same meaning across splits. No training or resampling is performed.
# A selected symbolic alert-none rule has no transferable numeric cutoff and is refused.
validate_data <- function(d) {
 stopifnot(is.data.frame(d),all(c("id","truth","score") %in% names(d)),nrow(d)>=2L,
   is.character(d$id),!anyNA(d$id),all(nzchar(d$id)),!anyDuplicated(d$id),
   is.numeric(d$truth),!anyNA(d$truth),all(d$truth %in% c(0,1)),
   is.numeric(d$score),all(is.finite(d$score)))
 invisible(d)
}
make_spec <- function(d,cutoff,direction,comparator,split,origin,conf.level,intervals) {
 validate_data(d)
 stopifnot(is.character(direction),length(direction)==1L,direction %in% c("higher_is_positive","lower_is_positive"),
  is.character(comparator),length(comparator)==1L,comparator %in% if(direction=="higher_is_positive") c(">=",">") else c("<=","<"))
 inferguidR::classification_spec(d,truth="truth",positive_label=1,negative_label=0,scores="score",id="id",
  models=inferguidR::classifier_model("score",score_meaning="declared numeric classifier score"),
  direction=direction,comparator=comparator,score_scale="raw_score",sampling="cross_sectional",evaluation_split=split,
  rule=inferguidR::cutoff_rule(cutoff=cutoff,selection_data=origin),
  uncertainty=inferguidR::classification_uncertainty(level=conf.level,
   proportion_methods=if(intervals) "clopper_pearson" else character(),auc_method="none"))
}
example_data<-data.frame(id=paste0("evaluation_",1:8),truth=c(0,1,0,1,0,1,0,1),score=c(.1,.2,.2,.6,.7,.8,.8,.9))
analyse <- function(tuning,evaluation,direction,comparator,independent_splits,conf.level=.95) {
 stopifnot(isTRUE(independent_splits));validate_data(tuning);validate_data(evaluation)
 stopifnot(length(unique(tuning$truth))==2L)
 train<-make_spec(tuning,0,direction,comparator,"held_out_tuning","held_out_tuning",conf.level,FALSE)
 test<-make_spec(evaluation,0,direction,comparator,"held_out_test","held_out_tuning",conf.level,TRUE)
 selected<-inferguidR::select_cutoff(train,"max_youden",tie_break="fewest_alerts")
 result<-inferguidR::evaluate_cutoff(test,cutoff=selected,include_curves=TRUE)
 list(method="Declared held-out evaluation of a tuning-selected Youden rule",status=result$status$code,
  selected_rule=selected,operating_point=result$models[[1]]$operating_point,
  metrics=inferguidR::classification_table(result,"metrics"),intervals=inferguidR::classification_table(result,"intervals"),
  summaries=inferguidR::classification_table(result,"summaries"),specification=inferguidR::classification_evaluated_spec(result),result=result,
  limitations="The prespecified example maximises Youden J on tuning data only, with fewest-alerts tie break. Evaluation outcomes never select the rule. Known shared unit IDs are rejected; disjoint IDs and a declaration cannot prove independence or exclude model-training leakage. Identifiers must share a namespace. Conditional marginal Clopper-Pearson intervals omit training and tuning uncertainty; no resampling validation is implemented here.")
}
example_tuning<-data.frame(id=paste0("tuning_",1:8),truth=c(0,0,1,0,1,0,1,1),score=c(.1,.2,.3,.4,.6,.7,.8,.9))
example_result<-analyse(example_tuning,example_data,"higher_is_positive",">=",independent_splits=TRUE,conf.level=confidence_level)
print(example_result)
