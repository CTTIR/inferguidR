# Guidance consumes canonical numerical and identity providers from infercalcR.
# Historical schema dispatch remains local; no numerical engines are duplicated.
.bl_control_record <- function(control) infercalcR::boschloo_control_record(control)
.bl_control_from_record <- function(rec) infercalcR::boschloo_control_from_record(rec)
.bl_check_control <- function(control) infercalcR::boschloo_validate_control(control)
.bl_check_design <- function(n1,n2) infercalcR::boschloo_validate_design(n1,n2)
.bl_variant <- function(variant) {
  if(is.null(variant)) .ig_abort_invalid("variant must be supplied explicitly")
  infercalcR::boschloo_variant_details(variant)
}
.bl_variant_info <- infercalcR::boschloo_variant_details()
.bl_evidence <- infercalcR::boschloo_engine_info()$evidence
.bl_max_tables <- infercalcR::boschloo_engine_info()$max_tables
.bl_provenance <- function() infercalcR::boschloo_engine_info()$provenance
.bl_tables <- function(n1,n2) infercalcR::binomial_table_space(n1,n2)
.bl_region <- function(n1,n2,alpha,info,control) infercalcR::boschloo_design_kernel(n1,n2,alpha,info$variant,control)
.bl_power_idx <- function(idx,n1,n2,p1,p2) infercalcR::binomial_region_probability(idx,n1,n2,p1,p2)
.bl_region_object <- function(reg,control) infercalcR::boschloo_region_record(reg,control)
.bl_region_status <- function(reg) infercalcR::boschloo_kernel_status(reg)
.bl_fisher_p <- function(n1,n2,info,control,convention="same") infercalcR::fisher_ordering_probabilities(n1,n2,info$variant,control,convention)
.bl_variant_record <- function(info) infercalcR::boschloo_variant_record(info)
.bl_control_call <- function(control) infercalcR::boschloo_control_call(control)
.bl_plotdata <- function(plot_id,title,data,fixed,scope,reference=list(),notes=character(),status=ig_status("ok")) infercalcR::boschloo_plot_data_record(plot_id,title,data,fixed,scope,reference,notes,status)
.ig_cls_dictionary <- function() infercalcR::classification_dictionary()
.ig_cls_operating_ids <- function() infercalcR::classification_metric_ids()
.ig_cls_reason_text <- function(code,lang="en") infercalcR::classification_reason(code,lang)
.ig_cls_short_label <- function(ids,lang) infercalcR::classification_short_label(ids,lang)
.ig_cls_check_thresholds <- function(threshold_probability) infercalcR::classification_check_thresholds(threshold_probability)
.ig_cls_check_bins <- function(bins) infercalcR::classification_check_bins(bins)
.ig_cls_check_partial <- function(partial_auc) infercalcR::classification_check_partial(partial_auc)
.ig_cls_metric_statuses <- infercalcR::classification_contract()$metric_statuses
.ig_cls_interval_statuses <- infercalcR::classification_contract()$interval_statuses
.ig_cls_families <- infercalcR::classification_contract()$plot_families
.ig_cls_table_columns <- infercalcR::classification_contract()$table_columns
.ig_cls_bootstrap_curve_ids <- infercalcR::classification_contract()$bootstrap_curve_ids
.ig_cls_constraint_names <- infercalcR::classification_contract()$constraint_names
.ig_sha256 <- function(txt) infercalcR::legacy_sha256(txt)
.ig_json_parse <- function(txt) infercalcR::legacy_parse_json(txt)
# Guidance converts its own power/advisor/report objects through retained dispatch.
.ig_canonical_json <- function(x) infercalcR::legacy_canonical_json(.ig_to_plain(x))
.ig_hash_value <- function(x) infercalcR::legacy_hash_value(.ig_to_plain(x))
.ig_json_max_bytes <- infercalcR::legacy_json_limits()$max_bytes
.ig_format_double <- function(x) infercalcR::legacy_format_double(x)
.ig_parse_decimal <- function(txt) infercalcR::legacy_parse_decimal(txt)
.ig_doubles <- function(x) infercalcR::legacy_double_values(x)
.ig_finite_paths <- function(x,path) infercalcR::legacy_nonfinite_paths(x,path)
.ig_r_literal <- function(v) infercalcR::legacy_r_literal(v)
.ig_unclass_frozen <- function(x) infercalcR::legacy_unclass_frozen(x)
.ig_p_norm <- function(q,mean=0,sd=1,lower_tail=TRUE,log_p=FALSE) infercalcR::normal_probability(q,mean,sd,lower_tail,log_p)
.ig_q_norm <- function(p,mean=0,sd=1,lower_tail=TRUE,log_p=FALSE) infercalcR::normal_quantile(p,mean,sd,lower_tail,log_p)
.ig_p_t <- function(q,df,ncp=0,lower_tail=TRUE,log_p=FALSE) infercalcR::t_probability(q,df,ncp,lower_tail,log_p)
.ig_q_t <- function(p,df,ncp=0,lower_tail=TRUE,log_p=FALSE) infercalcR::t_quantile(p,df,ncp,lower_tail,log_p)
.ig_p_chisq <- function(q,df,ncp=0,lower_tail=TRUE,log_p=FALSE) infercalcR::chisq_probability(q,df,ncp,lower_tail,log_p)
.ig_q_chisq <- function(p,df,ncp=0,lower_tail=TRUE,log_p=FALSE) infercalcR::chisq_quantile(p,df,ncp,lower_tail,log_p)
.ig_p_f <- function(q,df1,df2,ncp=0,lower_tail=TRUE,log_p=FALSE) infercalcR::f_probability(q,df1,df2,ncp,lower_tail,log_p)
.ig_q_f <- function(p,df1,df2,ncp=0,lower_tail=TRUE,log_p=FALSE) infercalcR::f_quantile(p,df1,df2,ncp,lower_tail,log_p)
.ig_p_binom <- function(q,size,prob,lower_tail=TRUE,log_p=FALSE) infercalcR::binomial_probability(q,size,prob,lower_tail,log_p)
.ig_logspace_add <- function(a,b) infercalcR::log_probability_add(a,b)
.ig_logspace_sum <- function(x) infercalcR::log_probability_sum(x)
.ig_invert_cdf <- function(logcdf,p,lower_tail,log_p,support,start) infercalcR::invert_log_cdf(logcdf,p,lower_tail,log_p,support,start)
.ig_f_lower_log_bound <- function(crit,df1,df2,ncp) infercalcR::f_lower_log_bound(crit,df1,df2,ncp)
.ig_calc_nct_pdf <- function(x,df,nc) infercalcR::noncentral_t_density(x,df,nc)
