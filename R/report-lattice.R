.ig_report_lattice_metadata <- function(result) {
  list(n1=result$group_sizes[["n1"]],n2=result$group_sizes[["n2"]],
    alpha=result$spec$inputs$alpha,p1=result$spec$inputs$p1,p2=result$spec$inputs$p2,
    variant=result$spec$options$variant,settings=result$method_outputs$settings,
    region_hash=result$method_outputs$region_hash,
    null_probability_at=result$method_outputs$attained_size_argmax %||% .5)
}

.ig_report_lattice_figure <- function(result,analysis,basis,data,status,reason,locale) {
  de <- locale=="de"
  list(analysis=analysis,plot_type="power_lattice",basis=basis,status=status,reason=reason,
    axes=list(x=if(de) "Ereignisse in Gruppe 1" else "Events in group 1",
      y=if(de) "Ereignisse in Gruppe 2" else "Events in group 2"),
    fixed_parameters=.ig_report_lattice_metadata(result),data=data,
    caption=if(de) paste("Alle moeglichen Tabellen; Kreuze markieren Ablehnung von H0, keine beobachteten Daten.",
      "H0 verwendet das gespeicherte numerische Maximum der Null-Ablehnungswahrscheinlichkeit; H1 die vorgegebenen Ereigniswahrscheinlichkeiten.",
      "Getrennte bedingte Wahrscheinlichkeiten, keine Vorhersagen oder Unsicherheitsintervalle. Kreisgroessen entsprechen Tabellenwahrscheinlichkeiten.") else
      paste("Full possible-table support; crosses mark rejection of H0, not observed data.",
        "H0 uses the retained numerical maximum of null rejection probability; H1 uses the specified event probabilities.",
        "Separate conditional probabilities, not predictions or uncertainty intervals. Circle areas represent table probabilities."))
}

.ig_report_lattice_figures <- function(result,analysis,locale) {
  pd <- power_plot_data(result,"lattice")
  lapply(c("H0","H1"),function(basis) {
    data <- if(pd$status=="ok") pd$data[pd$data$hypothesis==basis,c("x1","x2","probability","reject"),drop=FALSE] else list()
    .ig_report_lattice_figure(result,analysis,basis,data,pd$status,pd$reason,locale)
  })
}

# Readback uses frozen coordinates, binomial arithmetic and the retained region
# hash. It does not rerun a Boschloo region solver or replace the recorded grid.
.ig_report_validate_lattice <- function(x) {
  bad <- function(message) .ig_abort("invalid_input",paste("Invalid report lattice:",message))
  equal <- function(a,b) identical(.ig_canonical_json(a),.ig_canonical_json(b))
  scalar <- function(v) is.numeric(v)&&length(v)==1L&&is.finite(v)
  near <- function(a,b) scalar(a)&&scalar(b)&&abs(a-b)<=1e-12+1e-10*abs(b)
  for(i in seq_along(x$specifications)) {
    figures <- Filter(function(v) identical(v$plot_type,"power_lattice")&&v$analysis==i,x$visualisations)
    entry <- x$specifications[[i]]
    expected <- identical(x$schema_version,"inferguidR.report_document/0.4.0")&&
      identical(entry$analysis_kind,"power")&&identical(entry$specification$method_id,.bl_method_id)&&
      identical(x$results[[i]]$status$code,"ok")
    if(!expected) {if(length(figures)) bad("lattice needs a completed Boschloo result and schema0.4");next}
    if(length(figures)!=2L || !identical(vapply(figures,`[[`,"","basis"),c("H0","H1"))) bad("missing or duplicated conditional distributions")
    r <- .ig_result_from_list(x$results[[i]]);meta <- .ig_report_lattice_metadata(r)
    n1 <- meta$n1;n2 <- meta$n2;count <- (n1+1)*(n2+1)
    rejects <- NULL
    for(v in figures) {
      if(!v$status %in% c("ok","unavailable")) bad("invalid display status")
      header <- .ig_report_lattice_figure(r,i,v$basis,v$data,v$status,v$reason,x$locale)
      if(!equal(v,.ig_report_plain(header))) bad("altered lattice metadata")
      if(v$status=="unavailable") {
        if(length(v$data)||!.ig_is_string(v$reason)||!nzchar(v$reason)) bad("unavailable grid must retain an explicit reason and no coordinates")
        next
      }
      if(count>5000||!is.null(v$reason)||!is.list(v$data)||length(v$data)!=count) bad("invalid lattice support")
      wanted_x1 <- rep(0:n1,n2+1);wanted_x2 <- rep(0:n2,each=n1+1)
      probs <- numeric(count);rejected <- logical(count)
      p1 <- if(v$basis=="H0") meta$null_probability_at else meta$p1
      p2 <- if(v$basis=="H0") meta$null_probability_at else meta$p2
      for(j in seq_along(v$data)) {
        row <- v$data[[j]]
        if(!is.list(row)||anyDuplicated(names(row))||!setequal(names(row),c("x1","x2","probability","reject"))||
          !scalar(row$x1)||!scalar(row$x2)||row$x1!=wanted_x1[j]||row$x2!=wanted_x2[j]||
          !is.logical(row$reject)||length(row$reject)!=1L||is.na(row$reject)||
          !near(row$probability,stats::dbinom(row$x1,n1,p1)*stats::dbinom(row$x2,n2,p2))||row$probability<0)
          bad("invalid coordinate or table probability")
        probs[j]<-row$probability;rejected[j]<-row$reject
      }
      if(!near(sum(probs),1)) bad("conditional probabilities do not sum to one")
      target <- if(v$basis=="H0") r$achieved$alpha_attained else r$achieved$power
      if(!near(sum(probs[rejected]),target)) bad("rejection mass disagrees with result")
      if(!is.null(rejects)&&!identical(rejects,rejected)) bad("conditional panels use different rejection sets")
      rejects<-rejected
      region_hash <- .ig_hash_value(list(design=list(n1=n1,n2=n2,n_total=n1+n2,n_tables=count),
        alpha=meta$alpha,variant_id=r$method_outputs$variant_id,settings=meta$settings,
        rejected=list(x1=I(wanted_x1[rejected]),x2=I(wanted_x2[rejected]))))
      if(!identical(region_hash,meta$region_hash)) bad("rejection set does not match retained region hash")
    }
    if(!identical(figures[[1]]$status,figures[[2]]$status)) bad("conditional panels disagree about availability")
  }
  invisible(TRUE)
}

.ig_report_draw_lattice <- function(v,locale) {
  d<-v$data;meta<-v$fixed_parameters;de<-locale=="de"
  x<-vapply(d,`[[`,0,"x1");y<-vapply(d,`[[`,0,"x2");p<-vapply(d,`[[`,0,"probability")
  rejected<-vapply(d,`[[`,TRUE,"reject")
  graphics::plot(c(-.5,meta$n1+.5),c(-.5,meta$n2+.5),type="n",xlab=v$axes$x,ylab=v$axes$y,
    main=paste(if(de) "Boschloo-Tabellengitter" else "Boschloo table lattice",v$basis),asp=1,xaxt="n",yaxt="n")
  ticks<-function(n) {
    at<-if(n<=10) 0:n else unique(round(pretty(c(0,n),n=5)))
    at[at>=0&at<=n]
  }
  graphics::axis(1,at=ticks(meta$n1));graphics::axis(2,at=ticks(meta$n2))
  graphics::symbols(x,y,circles=sqrt(p/pi),inches=FALSE,add=TRUE,bg="#167C80",fg="#167C80")
  graphics::points(x[rejected],y[rejected],pch=4,cex=.8)
  graphics::mtext(paste(if(de) "Variante:" else "Variant:",meta$variant,
    if(v$basis=="H0") paste("| p =",format(meta$null_probability_at)) else paste("| p1 =",meta$p1,"p2 =",meta$p2)),side=3,line=.6,cex=.8)
  invisible(NULL)
}
