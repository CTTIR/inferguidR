# Validate retained classifier aggregates without executing stored analysis code.
.ig_report_validate_classifier <- function(r, s) {
 bad <- function(message) .ig_abort("invalid_input",paste("Invalid classifier report:",message))
 scalar <- function(x) is.numeric(x)&&length(x)==1L&&is.finite(x)
 equal <- function(a,b) identical(.ig_canonical_json(a),.ig_canonical_json(b))
 near <- function(a,b) scalar(a)&&scalar(b)&&abs(a-b)<=1e-10*max(1,abs(a),abs(b))
 fields <- function(x,required,optional=character()) {
  if(!is.list(x)||anyDuplicated(names(x))||!all(required %in% names(x))||any(!names(x)%in%c(required,optional))) bad("missing or unexpected nested fields")
 }
 dict <- .ig_cls_dictionary()$metrics
 cell <- function(x,j) if(is.null(x)) NULL else if(is.list(x)) x[[j]] else x[j]
 table <- function(x,kind) {
  if(is.null(x)) return(invisible(NULL))
  columns <- .ig_cls_table_columns[[kind]];fields(x,names(columns))
  sizes <- vapply(x,function(v) max(1L,length(v)),0L)
  if(length(unique(sizes))!=1L) bad("table columns have inconsistent lengths")
  n <- sizes[1L]
  for(j in seq_len(n)) {
   for(k in names(columns)) {
    v <- cell(x[[k]],j)
    if(!is.null(v) && !switch(columns[[k]],character=.ig_is_string(v),double=scalar(v),logical=is.logical(v)&&length(v)==1L&&!is.na(v))) bad("invalid table cell type")
   }
   id <- cell(x$metric_id,j);st <- cell(x$status,j)
   statuses<-if(kind=="intervals") .ig_cls_interval_statuses else .ig_cls_metric_statuses
   if(!.ig_is_string(id)||!id%in%names(dict)||!.ig_is_string(st)||!st%in%statuses) bad("unknown metric or status")
   vals <- if(kind=="intervals") c("lower","upper") else "value"
   for(k in vals) {
    v <- cell(x[[k]],j)
    if((st=="ok"&&!scalar(v))||(st!="ok"&&!is.null(v))) bad("metric value disagrees with status")
    if(!is.null(v)) {
     domain<-dict[[id]]$domain
     if((!is.null(domain$lower)&&v<domain$lower-1e-10)||(!is.null(domain$upper)&&v>domain$upper+1e-10)) bad("metric is outside its domain")
    }
   }
   if(st!="ok" && !.ig_is_string(cell(x$reason_code,j))) bad("unavailable metric needs a reason")
   if(kind=="intervals") {
    level<-cell(x$level,j)
    if(!scalar(level)||level<=0||level>=1) bad("invalid interval confidence level")
    if(st=="ok"&&cell(x$lower,j)>cell(x$upper,j)) bad("reversed interval")
   }
  }
  invisible(NULL)
 }
 counts <- function(x,kind) {
  if(!is.list(x)||!all(c("tp","fn","fp","tn")%in%names(x))) bad("missing confusion cells")
  cts<-unlist(x[c("tp","fn","fp","tn")],use.names=TRUE)
  if(length(cts)!=4L||!is.numeric(cts)||any(!is.finite(cts))||any(cts<0)||sum(cts)<=0) bad("invalid confusion cells")
  if(kind=="observed"&&any(cts!=floor(cts))) bad("observed cells must be integers")
  totals<-list(n_events=cts[['tp']]+cts[['fn']],n_non_events=cts[['fp']]+cts[['tn']],
    n_alerts=cts[['tp']]+cts[['fp']],n_non_alerts=cts[['fn']]+cts[['tn']],n_total=sum(cts))
  for(k in intersect(names(totals),names(x))) if(!near(x[[k]],totals[[k]])) bad("confusion margins disagree with cells")
  c(as.list(cts),totals)
 }
 common_metrics <- function(x,cts) {
  if(is.null(x)) return(invisible(NULL))
  ratios<-list(prevalence=c(cts$n_events,cts$n_total),alert_rate=c(cts$n_alerts,cts$n_total),
   tpr=c(cts$tp,cts$n_events),fnr=c(cts$fn,cts$n_events),tnr=c(cts$tn,cts$n_non_events),
   fpr=c(cts$fp,cts$n_non_events),ppv=c(cts$tp,cts$n_alerts),fdr=c(cts$fp,cts$n_alerts),
   npv=c(cts$tn,cts$n_non_alerts),for_rate=c(cts$fn,cts$n_non_alerts),accuracy=c(cts$tp+cts$tn,cts$n_total))
  ids<-unlist(x$metric_id,use.names=FALSE)
  if(anyDuplicated(ids)) bad("duplicate operating metric")
  for(id in intersect(ids,c(names(cts),names(ratios)))) {
   j<-match(id,ids);st<-cell(x$status,j);v<-cell(x$value,j)
   if(id%in%names(cts)) {expected<-cts[[id]];defined<-TRUE} else {
    pair<-ratios[[id]];defined<-pair[2]>0;expected<-if(defined) pair[1]/pair[2] else NULL
   }
   if(defined && (st!="ok"||!near(v,expected))) bad("metric disagrees with confusion cells")
   if(!defined && (st!="undefined"||!is.null(v))) bad("zero denominator must remain undefined")
  }
 }
 if(!identical(r$input_kind,s$input_kind)) bad("input kind disagrees with specification")
 if(!is.list(r$models)||length(r$models)!=length(s$models)) bad("model inventory disagrees with specification")
 ids<-vapply(r$models,function(m) if(.ig_is_string(m$model_id)) m$model_id else "","")
 expected_ids<-vapply(s$models,`[[`,"","model_id")
 if(anyDuplicated(ids)||!setequal(ids,expected_ids)) bad("model identifiers disagree with specification")
 for(m in r$models) {
  fields(m,c("model_id","status"),c("name","version","score_meaning","score_scale","direction","comparator",
    "tie_rule","evaluation_split","count_kind","operating_point","metrics","intervals","summaries","selection","scenario"))
  sm<-s$models[[match(m$model_id,expected_ids)]]
  for(k in c("name","version","score_meaning","score_scale","direction","comparator","tie_rule","evaluation_split"))
   if(!equal(m[[k]],sm[[k]])) bad("model metadata disagrees with specification")
  fields(m$status,c("code","reason","details"))
  if(!.ig_is_string(m$status$code)||!m$status$code%in%ig_status_codes()) bad("unknown model status")
  if(!.ig_is_string(m$count_kind)||!m$count_kind%in%c("observed","frequency_weighted","expected")) bad("unknown count kind")
  for(kind in c("metrics","summaries","intervals")) table(m[[kind]],kind)
  if(!is.null(m$operating_point)) {
   if(m$status$code!="ok") bad("unavailable model has operating counts")
   op<-m$operating_point
   fields(op,c("tp","fn","fp","tn","n_events","n_non_events","n_alerts","n_non_alerts","n_total","count_kind","comparator"),
     c("cutoff","cutoff_kind","rule_text","cutoff_scale"))
   if(!identical(op$count_kind,m$count_kind)||!identical(op$comparator,m$comparator)) bad("operating-point metadata mismatch")
   cts<-counts(op,m$count_kind);common_metrics(m$metrics,cts)
   if(s$input_kind=="counts_2x2") for(k in c("tp","fn","fp","tn"))
    if(!near(op[[k]],s$counts[[k]])) bad("operating cells disagree with supplied counts")
  }
  if(!is.null(m$scenario)) {
   sc<-m$scenario
   fields(sc,c("status","reason_code","prevalence","n_total","count_kind","source","assumption"),c("counts","metrics"))
   if(!.ig_is_string(sc$status)||!sc$status%in%.ig_cls_metric_statuses) bad("unknown scenario status")
   if(sc$status=="ok") {
    if(!identical(sc$count_kind,"expected")||!scalar(sc$prevalence)||sc$prevalence<0||sc$prevalence>1) bad("invalid scenario declaration")
    cts<-counts(sc$counts,"expected");if(!near(sc$n_total,cts$n_total)||!near(cts$n_events/cts$n_total,sc$prevalence)) bad("scenario margins disagree")
    table(sc$metrics,"metrics");common_metrics(sc$metrics,cts)
   }
  }
 }
 .ig_report_validate_comparisons(r,s)
 invisible(TRUE)
}
