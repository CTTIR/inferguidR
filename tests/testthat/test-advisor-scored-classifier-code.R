.scored_template <- function(variant) {
 e<-new.env();e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.",variant,".v1.R"),package="inferguidR"),e)));e
}
.scored_data <- function(i=0L) data.frame(id=paste0("unit_",1:12),
 truth=rep(c(0,1),6),score=((1:12)*(i+3L) %% 11)/10)
.manual_cells <- function(d,cut,op) {
 positive<-switch(op,">="=d$score>=cut,">"=d$score>cut,"<="=d$score<=cut,"<"=d$score<cut)
 c(tp=sum(positive & d$truth==1),fn=sum(!positive & d$truth==1),fp=sum(positive & d$truth==0),tn=sum(!positive & d$truth==0))
}
test_that("scored templates preserve manual decisions and empirical pairwise AUC", {
 x<-.scored_template("classifier_cutoff_exploration");f<-.scored_template("classifier_predeclared_rule")
 for(i in 0:3) for(op in c(">=",">","<=","<")) {
  d<-.scored_data(i);direction<-if(substr(op,1,1)==">") "higher_is_positive" else "lower_is_positive"
  explored<-x$analyse(d,direction,op);rows<-explored$decision_sets
  expected<-t(vapply(seq_len(nrow(rows)),function(j) {
   if(rows$cutoff_kind[j]=="alert_none") return(c(tp=0,fn=6,fp=0,tn=6))
   if(rows$cutoff_kind[j]=="alert_all") return(c(tp=6,fn=0,fp=6,tn=0))
   .manual_cells(d,rows$cutoff[j],op)
  },numeric(4)))
  expect_equal(unname(as.matrix(rows[c("tp","fn","fp","tn")])),unname(expected))
  expect_equal(nrow(rows),length(unique(d$score))+1L)
  positive<-d$score[d$truth==1];negative<-d$score[d$truth==0]
  differences<-outer(positive,negative,"-");if(direction=="lower_is_positive") differences<--differences
  auc<-mean((differences>0)+.5*(differences==0))
  summary<-explored$summaries
  expect_equal(summary$value[summary$metric_id=="roc_auc"],auc)
  expect_null(explored$result$models[[1]]$intervals)
  for(cut in c(-1,.2,.5,2)) {
   result<-f$analyse(d,cut,direction,op,TRUE)
   expect_equal(unlist(result$operating_point[c("tp","fn","fp","tn")]),.manual_cells(d,cut,op))
   expect_identical(result$result$evaluation$label,"predeclared_rule")
  }
 }
})
test_that("held-out rule selection depends only on tuning outcomes", {
 h<-.scored_template("classifier_heldout_validation")
 for(op in c(">=",">","<=","<")) {
  tuning<-.scored_data();evaluation<-.scored_data(1);evaluation$id<-paste0("test_",1:12)
  direction<-if(substr(op,1,1)==">") "higher_is_positive" else "lower_is_positive"
  if(direction=="lower_is_positive") {tuning$score<--tuning$score;evaluation$score<--evaluation$score}
  r<-h$analyse(tuning,evaluation,direction,op,TRUE)
  grid<-classification_table(r$selected_rule,"curve")
  j<-grid$tp/6-grid$fp/6;indices<-which(abs(j-max(j))<1e-14)
  best<-indices[which.min(grid$alerts[indices])]
  selected<-r$selected_rule$models[[1]]$operating_point
  expect_equal(unlist(selected[c("tp","fn","fp","tn")]),unlist(grid[best,c("tp","fn","fp","tn")]),ignore_attr=TRUE)
  changed<-evaluation;changed$truth<-1-changed$truth
  r2<-h$analyse(tuning,changed,direction,op,TRUE)
  expect_identical(r$selected_rule,r2$selected_rule)
  expect_identical(r$operating_point$cutoff,r2$operating_point$cutoff)
  expect_identical(r$result$evaluation$transfer$overlap_check$unit_identity,"no_overlap_detected")
  expect_true("split_independence_unverified" %in% r$result$limitations$code)
  expect_error(h$analyse(tuning,tuning,direction,op,TRUE),"same data|share unit identifiers")
  overlap<-evaluation;overlap$id[1]<-tuning$id[1]
  expect_error(h$analyse(tuning,overlap,direction,op,TRUE),"share unit identifiers")
  expect_error(h$analyse(tuning,evaluation,direction,op,FALSE))
  flat<-tuning;flat$score<-0
  if(op %in% c(">=","<=")) expect_error(h$analyse(flat,evaluation,direction,op,TRUE),"no numeric cutoff") else expect_equal(h$analyse(flat,evaluation,direction,op,TRUE)$operating_point$cutoff,0)
 }
})
test_that("scored examples reject malformed data and keep endpoint metrics typed", {
 f<-.scored_template("classifier_predeclared_rule");x<-.scored_template("classifier_cutoff_exploration");d<-.scored_data()
 expect_error(f$analyse(d,.5,"higher_is_positive",">=",FALSE))
 for(kind in c("missing","infinite","duplicate","truth","emptyid")) {
  bad<-d
  if(kind=="missing") bad$score[1]<-NA_real_
  if(kind=="infinite") bad$score[1]<-Inf
  if(kind=="duplicate") bad$id[2]<-bad$id[1]
  if(kind=="truth") bad$truth[1]<-2
  if(kind=="emptyid") bad$id[1]<-""
  expect_error(x$analyse(bad,"higher_is_positive",">="))
  expect_error(f$analyse(bad,.5,"higher_is_positive",">=",TRUE))
 }
 expect_error(x$analyse(d,"lower_is_positive",">="))
 expect_error(f$analyse(d,NA_real_,"higher_is_positive",">=",TRUE))
 expect_error(f$analyse(d,.5,"higher_is_positive",">=",TRUE,conf.level=1))
 d$truth<-0;r<-f$analyse(d,.5,"higher_is_positive",">=",TRUE)
 expect_identical(r$metrics$status[r$metrics$metric_id=="tpr"],"undefined")
})
test_that("public ENDE scored classifier exports state their evaluation design", {
 variants<-c(explore="classifier_cutoff_exploration",predeclared="classifier_predeclared_rule",heldout="classifier_heldout_validation")
 for(purpose in names(variants)) for(lang in c("en","de")) {
  s<-advisor_replay(c(q.aim="performance",q.stage="analysis",q.perf_inputs="truth_and_scores",q.perf_purpose=purpose),language=lang)
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,unname(variants[purpose]))
  code<-advisor_r_code(s);e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
  expect_identical(e$example_result$status,"ok")
  expect_false(exists("alternative",e,inherits=FALSE))
  expect_identical(exists("confidence_level",e,inherits=FALSE),purpose!="explore")
  expect_error(advisor_r_code(s,alternative="two.sided"),class="inferguidR_invalid_input")
 }
})
