# Synthetic independent-unit example. Explicit score direction is never chosen to improve AUC.
# Identifiers must have the same meaning across splits. No training or resampling is performed.
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
analyse <- function(data,direction,comparator) {
 spec<-make_spec(data,0,direction,comparator,"same_data","same_data",.95,FALSE)
 result<-inferguidR::evaluate_cutoff(spec,include_curves=TRUE)
 list(method="Descriptive cutoff exploration",status=result$status$code,
  decision_sets=inferguidR::classification_table(result,"curve"),summaries=inferguidR::classification_table(result,"summaries"),
  specification=spec,result=result,
  limitations="All attainable decision sets are enumerated on these data. No optimal cutoff is certified, no held-out validation or sampling interval is claimed. Scores must belong to independent units; no weights, clusters, missing-data deletion or probability calibration is implied.")
}
example_result<-analyse(example_data,"higher_is_positive",">=")
print(example_result)
