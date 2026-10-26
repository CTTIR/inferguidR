rl_report <- function(locale="en") {
 s<-analysis_spec("boschloo_proportions_two_independent","post_hoc",
   list(p1=.2,p2=.7,alpha=.13,n1=2,n2=2),list(variant="central"))
 create_report_snapshot(analysis_envelope(s,calculate_power(s)),locale=locale)
}
rl_rehash <- function(x) {
 x<-unclass(x);x$snapshot_hash<-NULL;x$snapshot_hash<-ig_hash(x);ig_to_json(x)
}
test_that("unequal-size and solved Boschloo variants retain numerical null maxima", {
 for(variant in c("central","minlike","less","greater")) {
  s<-analysis_spec("boschloo_proportions_two_independent","post_hoc",
    list(p1=.23,p2=.64,alpha=.05,n1=3,n2=5),list(variant=variant))
  r<-calculate_power(s);expect_identical(r$status$code,"ok")
  doc<-create_report_snapshot(analysis_envelope(s,r))
  expect_identical(ig_from_json(ig_to_json(doc)),doc)
  for(v in doc$visualisations) {
   expect_length(v$data,24)
   expect_equal(sum(vapply(Filter(function(row)row$reject,v$data),`[[`,0,"probability")),
     if(v$basis=="H0") r$achieved$alpha_attained else r$achieved$power,tolerance=1e-10)
   expect_equal(v$fixed_parameters$null_probability_at,r$method_outputs$attained_size_argmax)
  }
 }
 s<-analysis_spec("boschloo_proportions_two_independent","a_priori",
   list(p1=.1,p2=.8,alpha=.1,power=.5,allocation_ratio=1.5,n_max=20),
   list(variant="minlike",allocation="n2_ceiling"))
 r<-calculate_power(s);expect_identical(r$status$code,"ok")
 doc<-create_report_snapshot(analysis_envelope(s,r))
 expect_identical(ig_from_json(ig_to_json(doc)),doc)
 expect_equal(doc$visualisations[[1]]$fixed_parameters$n1,unname(r$group_sizes['n1']))
 expect_equal(doc$visualisations[[1]]$fixed_parameters$n2,unname(r$group_sizes['n2']))
})

test_that("an unavailable bounded display is frozen without manufacturing a lattice", {
 # Exercise the genuine display-cap guard with a synthetic oversized design;
 # use its unavailable result on an otherwise valid report fixture to avoid
 # computing an expensive rejection region merely to test display handling.
 pd<-.ig_power_boschloo_lattice(list(method_id=.bl_method_id,
   group_sizes=c(n1=100L,n2=100L)),list(status="unavailable",data=data.frame()))
 expect_match(pd$reason,"5000-table display limit",fixed=TRUE)
 with_mocked_bindings({
  doc<-rl_report()
  expect_identical(vapply(doc$visualisations,`[[`,"","status"),rep("unavailable",2))
  expect_true(all(vapply(doc$visualisations,function(v)length(v$data)==0,logical(1))))
  expect_identical(ig_from_json(ig_to_json(doc)),doc)
 },power_plot_data=function(...)pd,.package="inferguidR")
})
test_that("Boschloo reports freeze both exact conditional lattices without solver replay", {
 for(lang in c("en","de")) {
  doc<-rl_report(lang)
  expect_identical(doc$schema_version,"inferguidR.report_document/0.4.0")
  expect_false(doc$privacy$probability_aggregates_included)
  expect_length(doc$probability_summaries,0)
  lattice<-Filter(function(v)v$plot_type=="power_lattice",doc$visualisations)
  expect_length(lattice,2)
  expect_identical(vapply(lattice,`[[`,"","basis"),c("H0","H1"))
  for(v in lattice) {
   expect_length(v$data,9)
   expect_equal(sum(vapply(v$data,`[[`,0,"probability")),1,tolerance=1e-12)
   rejected<-Filter(function(r)r$reject,v$data)
   expect_length(rejected,2)
   expect_setequal(vapply(rejected,function(r)paste(r$x1,r$x2),""),c("0 2","2 0"))
  }
  alt<-lattice[[2]]$data
  expected<-.8^2*.7^2+.2^2*.3^2
  expect_equal(sum(vapply(Filter(function(r)r$reject,alt),`[[`,0,"probability")),expected)
  with_mocked_bindings({
   expect_identical(ig_from_json(ig_to_json(doc)),doc)
   if(report_pdf_capabilities()$available) {
    f<-tempfile(fileext=".pdf");on.exit(unlink(f),add=TRUE)
    export_report(doc,f,"pdf");expect_gt(file.size(f),1000)
   }
  },power_plot_data=function(...) stop("must not recompute display"),.package="inferguidR")
 }
})
test_that("resealed lattice tampering cannot change support, mass, region or metadata", {
 doc<-rl_report()
 for(mutate in list(
  function(x){x$visualisations[[1]]$data[[1]]$x1<-1;x},
  function(x){x$visualisations[[1]]$data[[1]]$probability<-.9;x},
  function(x){x$visualisations[[1]]$data[[1]]$reject<-TRUE;x},
  function(x){x$visualisations[[1]]$fixed_parameters$null_probability_at<-.3;x},
  function(x){x$visualisations[[1]]$data<-x$visualisations[[1]]$data[-1];x},
  function(x){x$visualisations<-x$visualisations[1];x},
  function(x){x$visualisations[[1]]$data[[1]]$raw_rows<-list(1);x},
  function(x){x$schema_version<-"inferguidR.report_document/0.3.0";x}
 )) expect_error(ig_from_json(rl_rehash(mutate(unclass(doc)))),class="inferguidR_invalid_input")
})

test_that("general XY reports retain different modes and failed input rows", {
 for(mode in c("post_hoc","a_priori")) {
  s<-analysis_spec("t_means_one_sample",mode,
   if(mode=="post_hoc") list(d=.5,n_total=30,alpha=.05) else list(d=.5,power=.8,alpha=.05),list(tails="two"))
  r<-calculate_power(s);sw<-sweep_design(s,"alpha",c(.01,0,.05))
  doc<-create_report_snapshot(analysis_envelope(s,r),report_kind="incomplete_settings",sweeps=list(list(analysis=1,sweep=sw)))
  plot<-Filter(function(v)v$plot_type=="power_sweep",doc$visualisations)[[1]]
  expect_identical(vapply(plot$data,`[[`,"","status"),c("ok","invalid_input","ok"))
  expect_null(plot$data[[2]]$y)
  expect_equal(vapply(plot$data,`[[`,0,"x"),c(.01,0,.05))
  expect_identical(ig_from_json(ig_to_json(doc)),doc)
 }
})
