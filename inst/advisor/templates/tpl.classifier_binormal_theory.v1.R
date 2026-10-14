# Hypothetical score distributions: S0 ~ N(0,1), S1 ~ N(d,sigma^2).
# Higher scores are positive; cutoff is on the scale standardised to non-events.
# Expected N scales counts only and is not a validated sample-size recommendation.
analyse <- function(auc,sigma,cutoff,prevalence,n_total) {
 stopifnot(is.numeric(cutoff),length(cutoff)==1L,is.finite(cutoff))
 spec<-inferguidR::classification_spec(theory=list(auc=auc,sigma=sigma,n_total=n_total),
   cutoff=cutoff,prevalence=prevalence,prevalence_source="assumed",
   models=inferguidR::classifier_model("binormal",name="Hypothetical binormal model",score_meaning="score standardised to the non-event distribution"),
   direction="higher_is_positive",comparator=">=",
   uncertainty=inferguidR::classification_uncertainty(proportion_methods=character(),auc_method="none"))
 result<-inferguidR::evaluate_cutoff(spec,include_curves=FALSE)
 m<-result$models[[1]]
 list(method="Hypothetical binormal cutoff scenario",status=result$status$code,
   model=inferguidR::binormal_model(auc,sigma),cutoff=cutoff,
   analytic_rates=inferguidR::binormal_rates(cutoff,auc,sigma),
   prevalence=prevalence,n_total=n_total,count_kind=m$count_kind,
   operating_point=m$operating_point,metrics=m$metrics,summaries=m$summaries,
   confidence_interval=NULL,interval_status="not_applicable",specification=spec,result=result,
   limitations="Assumed score laws and prevalence, not observed performance or external validation. No sampling interval is inferred from expected counts. Changing N scales expected counts without changing analytic rates. Class-conditional rates remain defined at prevalence endpoints even when count-based conditional metrics have zero denominators. No classifier training, cutoff optimisation or sample-size planning is performed.")
}
example_result<-analyse(auc=.8,sigma=1.5,cutoff=1,prevalence=.1,n_total=1000)
print(example_result)
