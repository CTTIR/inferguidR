# Validate portable content without evaluating stored calls or recomputing analyses.
.ig_report_validate_content <- function(x) {
  bad <- function(message) .ig_abort("invalid_input", paste("Invalid report:", message))
  fields <- function(value, required, optional = character()) {
    if (!is.list(value) || anyDuplicated(names(value)) ||
        !all(required %in% names(value)) || any(!names(value) %in% c(required, optional)))
      bad("missing or unexpected nested fields")
  }
  equal <- function(a,b) identical(.ig_canonical_json(a), .ig_canonical_json(b))
  scalar <- function(v) is.numeric(v) && length(v)==1L && is.finite(v)
  probability <- function(v) scalar(v) && v>=0 && v<=1
  array <- function(v) is.list(v) && (length(v)==0L || is.null(names(v)))
  for (name in c("specifications","results","metrics","uncertainty","comparisons","sweeps",
                 "visualisations","validation","reproducibility"))
    if (!array(x[[name]])) bad(paste(name,"must be an array"))
  if (length(x$sweeps) && identical(x$schema_version,"inferguidR.report_document/0.1.0")) bad("sweep content is not supported by this report schema version")
  .ig_report_validate_sweeps(x)
  .ig_report_validate_resources(x)
  n <- length(x$specifications)
  if (!n || any(vapply(x[c("results","metrics","uncertainty","comparisons","validation","reproducibility")],length,0L)!=n))
    bad("per-analysis sections have inconsistent lengths")
  if (!.ig_is_string(x$title) || !.ig_is_string(x$snapshot_id) ||
      !.ig_is_string(x$display_timezone) || !x$display_timezone %in% c("UTC",OlsonNames()) ||
      !.ig_is_string(x$report_kind) || !x$report_kind %in% .ig_report_kinds ||
      !.ig_is_string(x$completeness) || !x$completeness %in% c("complete","incomplete")) bad("invalid metadata")
  fields(x$provenance,c("package","package_version","r_version","platform"))
  if (!all(vapply(x$provenance,.ig_is_string,TRUE))) bad("invalid runtime provenance")
  probability_schema <- identical(x$schema_version,"inferguidR.report_document/0.4.0")
  fields(x$privacy,c("raw_rows_included","score_curve_tables_included","note",
    if(probability_schema) "probability_aggregates_included"))
  if (probability_schema && (!is.logical(x$privacy$probability_aggregates_included)||
      length(x$privacy$probability_aggregates_included)!=1L||is.na(x$privacy$probability_aggregates_included))) bad("probability aggregates require an explicit inclusion choice")
  if (!isFALSE(x$privacy$raw_rows_included) || !isFALSE(x$privacy$score_curve_tables_included) ||
      !.ig_is_string(x$privacy$note)) bad("invalid privacy declaration")
  if (x$report_kind=="single_analysis" && n!=1L)
    bad("analysis count does not match report scope")
  kinds <- vapply(x$specifications,function(entry) {
    fields(entry,c("analysis_kind","analysis_hash","specification"))
    if (!.ig_is_string(entry$analysis_kind) || !entry$analysis_kind %in% .ig_analysis_kinds ||
        !.ig_is_string(entry$analysis_hash) || !grepl("^[0-9a-f]{64}$",entry$analysis_hash)) bad("invalid analysis reference")
    entry$analysis_kind
  },"")
  multiple_models <- n==1L && kinds[1]=="classification" &&
    is.list(x$specifications[[1]]$specification) && length(x$specifications[[1]]$specification$models)>=2L
  if(x$report_kind=="comparison" && n<2L && !multiple_models) bad("comparison needs multiple analyses or classifier models")
  if(x$report_kind=="model_cutoff" && !"classification" %in% kinds) bad("model scope needs classifier content")
  codes <- character(n)
  for(i in seq_len(n)) {
    entry <- x$specifications[[i]];r <- x$results[[i]];kind <- kinds[i]
    s <- if(kind=="power") .ig_spec_from_list(entry$specification) else .ig_aux_spec_from_list(entry$specification,kind)
    if(!is.list(r) || !is.list(r$status) || !.ig_is_string(r$status$code)) bad("missing result status")
    codes[i] <- r$status$code
    missing <- identical(codes[i],"not_computed")
    if(missing) {
      fields(r,"status");fields(r$status,c("code","reason"))
      if(!.ig_is_string(r$status$reason)) bad("missing result needs a reason")
    } else if(!codes[i] %in% ig_status_codes()) bad("unknown result status")
    if(kind=="classification" && !missing)
      fields(r,c("status","input_kind","models"),c("prevalence","evaluation","theory","settings",
        "diagnostics","limitations","comparisons","validation","data_fingerprint"))
    if(kind=="classification" && !missing) .ig_report_validate_classifier(r,s)
    original <- NULL
    if(kind=="power" && !missing) original <- .ig_result_from_list(r)
    if(kind=="observed_test") {
      if(!s$variant_id %in% c("boschloo_minlike","boschloo_central","boschloo_lower","boschloo_upper")) bad("unsupported observed method")
      d <- s$inputs
      for(name in c("x1","n1","x2","n2")) if(!scalar(d[[name]]) || d[[name]]<0 || d[[name]]!=floor(d[[name]])) bad("invalid observed counts")
      if(d$n1<1 || d$n2<1 || d$x1>d$n1 || d$x2>d$n2) bad("invalid observed group sizes")
      if(!missing) {
        fields(r,c("schema_version","object_type","method","status","spec_hash","variant","data",
          "p_value","p_upper","alpha","decision","smaller_one_sided","fisher_p_value","fisher_note",
          "sides","undefined","settings","scope","provenance","r_call"))
        if(!identical(r$method,"boschloo") || !identical(r$spec_hash,s$spec_hash)) bad("observed method/specification mismatch")
        if(codes[i]=="ok" && (!probability(r$p_value) || !probability(r$p_upper) || r$p_upper<r$p_value ||
            !probability(r$fisher_p_value))) bad("invalid p-value or numerical bound")
        if(!is.null(r$alpha) && (!probability(r$alpha) || r$alpha==0 || r$alpha==1)) bad("invalid significance level")
        if(!equal(r$data[c("x1","n1","x2","n2")],d[c("x1","n1","x2","n2")])) bad("observed counts disagree with specification")
        original <- r;original$spec <- s;class(original) <- c("observed_test_result","ig_frozen")
      }
    }
    if(kind!="classification" || missing) {
      rebuilt <- analysis_envelope(s,original)
      if(!identical(rebuilt$analysis_hash,entry$analysis_hash)) bad("analysis hash does not match retained content")
    }
    expected <- if(kind=="power") list(analysis=i,achieved=r$achieved,solved=r$solved,quantities=r$quantities) else
      if(kind=="observed_test") list(analysis=i,p_value=r$p_value,p_upper=r$p_upper,alpha=r$alpha,decision=r$decision,variant=s$options$variant) else
        list(analysis=i,models=lapply(r$models,function(m) list(model_id=m$model_id,metrics=m$metrics,summaries=m$summaries,scenario=m$scenario)))
    if(!equal(x$metrics[[i]],expected)) bad("metric summary disagrees with result")
    if(!equal(x$comparisons[[i]],r$comparisons)) bad("comparison summary disagrees with result")
    retained <- if(missing) NULL else r
    if(kind=="classification") {
      expected_uncertainty<-list(analysis=i,settings=s$uncertainty,
        models=lapply(r$models,function(m) list(model_id=m$model_id,intervals=m$intervals)))
      expected_validation<-list(analysis=i,evidence=retained$validation)
      reproduction<-x$reproducibility[[i]]
      fields(reproduction,c("analysis","r_code","private_data_required","note"))
      if((missing && !is.null(reproduction$r_code)) ||
         (!missing && !.ig_is_string(reproduction$r_code))) bad("invalid classifier reproduction text")
      expected_reproduction<-list(analysis=i,r_code=reproduction$r_code,
        private_data_required=identical(s$input_kind,"scored_data"),
        note="Reattach the original private input using its fingerprint before recomputing scored-data results.")
    } else if(kind=="power") {
      expected_uncertainty<-list(analysis=i,status="not_applicable",
        reason="Prospective power calculation; no sampling interval for the design inputs.")
      expected_validation<-list(analysis=i,status="research_preview",
        reason="Complete numerical release certification remains pending.",evidence=retained$validation)
      expected_reproduction<-list(analysis=i,r_code=if(missing) .ig_spec_call(s) else retained$r_call,
        private_data_required=FALSE)
    } else {
      expected_uncertainty<-list(analysis=i,status="not_computed",
        reason="No sampling confidence interval is computed. p_value and p_upper bound numerical optimisation error, not sampling uncertainty.",
        numerical_bounds=list(lower=retained$p_value,upper=retained$p_upper),evidence=retained$sides)
      expected_validation<-list(analysis=i,status="research_preview",
        reason="Complete numerical release certification remains pending.",
        computation_status=retained$status,settings=s$options)
      expected_reproduction<-list(analysis=i,r_code=retained$r_call,private_data_required=FALSE,
        note="Independent binomial groups; counts and event orientation are retained.")
    }
    normalise <- function(value) .ig_json_parse(.ig_canonical_json(.ig_report_plain(value)))
    if(!equal(x$uncertainty[[i]],normalise(expected_uncertainty))) bad("uncertainty summary disagrees with result")
    if(!equal(x$validation[[i]],normalise(expected_validation))) bad("validation summary disagrees with retained evidence")
    if(!equal(x$reproducibility[[i]],normalise(expected_reproduction))) bad("reproduction metadata disagrees with retained analysis")
  }
  complete <- all(codes=="ok") && all(vapply(x$sweeps,function(sw) all(vapply(sw$rows,function(row) identical(row$status,"ok"),TRUE)),TRUE))
  if(!identical(x$completeness,if(complete) "complete" else "incomplete") ||
      (!complete && x$report_kind!="incomplete_settings")) bad("invalid completeness claim")
  rec <- if(is.null(x$recommendation)) NULL else .ig_recommendation_from_list(x$recommendation)
  expected_hash <- .ig_hash_value(list(report_kind=x$report_kind,
    envelopes=I(vapply(x$specifications,`[[`,"","analysis_hash")),
    recommendation=if(is.null(rec)) NULL else rec$recommendation_hash))
  if(!identical(x$analysis_hash,expected_hash)) bad("report analysis hash mismatch")
  for(v in x$visualisations) {
    fields(v,c("analysis","plot_type","status","axes","data","caption"),c("reason","fixed_parameters","critical","model_id","basis","sweep_index","target"))
    if(!is.list(v) || !scalar(v$analysis) || v$analysis!=floor(v$analysis) || v$analysis<1 || v$analysis>n ||
        !.ig_is_string(v$plot_type) || !.ig_is_string(v$status)) bad("invalid visualisation reference")
    if(!v$plot_type %in% c("distribution","confusion_counts","two_group_counts","classifier_mosaic","classifier_prevalence","classifier_comparison","classifier_roc_point","classifier_pr_point","classifier_calibration","classifier_net_benefit","power_sweep","power_lattice")) bad("unknown plot type")
    if(v$plot_type=="power_sweep") {
      expected <- .ig_report_sweep_figures(x$sweeps,x$locale)
      matching <- Filter(function(z) identical(as.double(z$sweep_index),as.double(v$sweep_index)) && identical(z$target,v$target),expected)
      if(length(matching)!=1L || !equal(v,.ig_report_plain(matching[[1]]))) bad("sweep figure disagrees with retained rows")
    }
    if(v$plot_type=="distribution") {
      if(kinds[v$analysis]!="power" || !identical(x$results[[v$analysis]]$status$code,"ok")) bad("distribution figure needs a completed power result")
      result<-.ig_result_from_list(x$results[[v$analysis]])
      pd<-power_plot_data(result)
      expected<-list(analysis=v$analysis,plot_type="distribution",status=pd$status,reason=pd$reason,
        axes=list(x="test statistic",y=pd$ordinate),fixed_parameters=result$spec$inputs,
        critical=pd$critical,data=pd$data,
        caption=paste(c(pd$notes,"Displayed densities or probability masses describe the fixed design."),collapse=" "))
      # Older snapshots may retain a failed display even when a newer renderer succeeds.
      # Preserve the explicit absence of coordinates; never substitute a new figure.
      if (identical(x$schema_version,"inferguidR.report_document/0.1.0") &&
          identical(v$status,"unavailable") && array(v$data) && !length(v$data) &&
          .ig_is_string(v$reason) && startsWith(v$reason,"Distribution display could not be computed:")) {
        expected$status <- "unavailable"
        expected$reason <- v$reason
        expected$axes$y <- NULL
        expected$axes <- list(x="test statistic", y=NULL)
        expected$data <- list()
        expected$caption <- "Displayed densities or probability masses describe the fixed design."
      }
      expected<-.ig_json_parse(.ig_canonical_json(.ig_report_plain(expected)))
      # Display-coordinate portability tolerance only; no engine/reference certification.
      coordinate_equal<-function(a,b) scalar(a)&&scalar(b)&&abs(a-b)<=1e-12+1e-10*max(abs(a),abs(b))
      matches <- function(candidate) {
        if(!equal(v[setdiff(names(v),"data")],candidate[setdiff(names(candidate),"data")]) ||
           length(v$data)!=length(candidate$data)) return(FALSE)
        for(j in seq_along(candidate$data)) {
          row<-v$data[[j]];ref<-candidate$data[[j]]
          fields(row,c("hypothesis","x","y","reject"))
          if(!identical(row$hypothesis,ref$hypothesis)||!identical(row$reject,ref$reject)||
             !coordinate_equal(row$x,ref$x)||!coordinate_equal(row$y,ref$y)||row$y<0) return(FALSE)
        }
        TRUE
      }
      matched <- matches(expected)
      if (!matched && identical(x$schema_version,"inferguidR.report_document/0.1.0") &&
          identical(v$status,"ok")) {
        legacy <- .ig_power_plot_data(result, legacy_t = TRUE)
        if (identical(legacy$status,"ok")) {
          previous <- list(analysis=v$analysis,plot_type="distribution",status=legacy$status,reason=legacy$reason,
            axes=list(x="test statistic",y=legacy$ordinate),fixed_parameters=result$spec$inputs,
            critical=legacy$critical,data=legacy$data,
            caption=paste(c(legacy$notes,"Displayed densities or probability masses describe the fixed design."),collapse=" "))
          previous <- .ig_json_parse(.ig_canonical_json(.ig_report_plain(previous)))
          matched <- matches(previous)
        }
      }
      if (!matched) bad("distribution metadata or coordinates disagree with retained result")
    }
    if(v$plot_type=="classifier_comparison") {
      if(kinds[v$analysis]!="classification" || !is.list(v$data) || !scalar(v$data$panel_index)) bad("comparison figure needs a classifier panel")
      expected<-.ig_report_classifier_comparison_figures(x$results[[v$analysis]],v$analysis,x$locale)
      j<-v$data$panel_index
      if(j!=floor(j) || j<1 || j>length(expected) || !equal(v,.ig_report_plain(expected[[j]]))) bad("comparison figure disagrees with retained pairs")
    }
    if(v$plot_type %in% c("classifier_roc_point","classifier_pr_point")) {
      if(kinds[v$analysis]!="classification" || !.ig_is_string(v$model_id)) bad("operating figure needs a classifier model")
      models<-x$results[[v$analysis]]$models;ids<-vapply(models,`[[`,"","model_id");j<-match(v$model_id,ids)
      if(is.na(j)) bad("unknown operating figure model")
      expected<-.ig_report_classifier_operating_figures(models[[j]],v$analysis,x$locale)
      same<-vapply(expected,function(z) identical(z$plot_type,v$plot_type),FALSE)
      if(sum(same)!=1L || !equal(v,.ig_report_plain(expected[[which(same)]]))) bad("operating figure disagrees with retained cells")
    }
    if(v$plot_type %in% c("classifier_mosaic","classifier_prevalence")) {
      if(kinds[v$analysis]!="classification" || !.ig_is_string(v$model_id)) bad("aggregate figure needs a classifier model")
      models<-x$results[[v$analysis]]$models;ids<-vapply(models,`[[`,"","model_id");j<-match(v$model_id,ids)
      if(is.na(j)) bad("unknown aggregate figure model")
      expected<-.ig_report_classifier_aggregate_figures(models[[j]],v$analysis,x$locale)
      same<-vapply(expected,function(z) identical(z$plot_type,v$plot_type)&&identical(z$basis,v$basis),FALSE)
      if(sum(same)!=1L || !equal(v,.ig_report_plain(expected[[which(same)]]))) bad("aggregate figure disagrees with retained cells")
    }
    if(v$plot_type=="confusion_counts") {
      if(kinds[v$analysis]!="classification" || !.ig_is_string(v$model_id)) bad("classifier figure needs a model reference")
      models<-x$results[[v$analysis]]$models
      ids<-vapply(models,`[[`,"","model_id");j<-match(v$model_id,ids)
      if(is.na(j)||is.null(models[[j]]$operating_point)) bad("unknown plotted model or operating point")
      if(!equal(v$data,models[[j]]$operating_point[c("tp","fn","fp","tn")])) bad("plotted classifier cells disagree with result")
    }
    if(v$plot_type=="two_group_counts") {
      if(kinds[v$analysis]!="observed_test") bad("group count figure needs observed-test input")
      d <- x$specifications[[v$analysis]]$specification$inputs
      expected <- list(group_labels=d$group_labels,event_label=d$event_label,
        event=c(d$x1,d$x2),non_event=c(d$n1-d$x1,d$n2-d$x2))
      if(!equal(v$data,expected)) bad("plotted counts disagree with specification")
    }
  }
  for(i in seq_len(n)) {
    distributions<-Filter(function(v) v$analysis==i && v$plot_type=="distribution",x$visualisations)
    expected_distribution<-as.integer(kinds[i]=="power" && identical(codes[i],"ok") &&
      !(identical(x$schema_version,"inferguidR.report_document/0.4.0")&&identical(x$specifications[[i]]$specification$method_id,.bl_method_id)))
    if(length(distributions)!=expected_distribution) bad("distribution figure is missing or duplicated")
    points<-Filter(function(v) v$analysis==i && v$plot_type %in% c("classifier_roc_point","classifier_pr_point"),x$visualisations)
    if(length(points)) {
      expected<-unlist(lapply(x$results[[i]]$models,function(m) .ig_report_classifier_operating_figures(m,i,x$locale)),recursive=FALSE)
      counts<-vapply(expected,function(z) sum(vapply(points,function(v) identical(v$model_id,z$model_id)&&identical(v$plot_type,z$plot_type),FALSE)),0L)
      if(length(points)!=length(expected) || any(counts!=1L)) bad("operating point figures are incomplete or duplicated")
    }
    shown<-Filter(function(v) v$analysis==i && v$plot_type=="classifier_comparison",x$visualisations)
    if(length(shown)) {
      expected<-.ig_report_classifier_comparison_figures(x$results[[i]],i,x$locale)
      ids<-vapply(shown,function(v) v$data$panel_index,0)
      if(length(shown)!=length(expected) || !setequal(ids,seq_along(expected))) bad("comparison panels are incomplete or duplicated")
    }
  }
  expected_sweeps <- .ig_report_sweep_figures(x$sweeps,x$locale)
  actual_sweeps <- Filter(function(v) v$plot_type=="power_sweep",x$visualisations)
  if(!equal(actual_sweeps,.ig_report_plain(expected_sweeps))) bad("sweep figures are missing or duplicated")
  .ig_report_validate_probability(x)
  .ig_report_validate_lattice(x)
  invisible(TRUE)
}
