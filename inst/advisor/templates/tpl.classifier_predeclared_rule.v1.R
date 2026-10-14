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
analyse <- function(data,cutoff,direction,comparator,predeclared,conf.level=.95) {
 stopifnot(isTRUE(predeclared))
 spec<-make_spec(data,cutoff,direction,comparator,"same_data","predeclared_rule",conf.level,TRUE)
 result<-inferguidR::evaluate_cutoff(spec,include_curves=TRUE)
 list(method="Declared pre-fixed cutoff evaluation",status=result$status$code,
  operating_point=result$models[[1]]$operating_point,metrics=inferguidR::classification_table(result,"metrics"),
  intervals=inferguidR::classification_table(result,"intervals"),summaries=inferguidR::classification_table(result,"summaries"),
  specification=spec,result=result,
  limitations="The analyst declares that this cutoff preceded inspection of these outcomes. This declaration is not proof. Marginal Clopper-Pearson intervals condition on the fixed rule and independent units; they do not include training or selection uncertainty. No validation of the fitted score model is established by fixing a cutoff.")
}
example_result<-analyse(example_data,.5,"higher_is_positive",">=",predeclared=TRUE,conf.level=confidence_level)
print(example_result)
