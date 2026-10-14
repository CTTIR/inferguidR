test_that("report summaries separate observed and expected counts without changing snapshots", {
 s<-classification_spec(counts=c(tp=40,fn=10,fp=20,tn=130),cutoff=.5,prevalence=.1)
 x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)),report_kind="model_cutoff")
 before<-ig_to_json(x);text<-paste(.ig_report_summary_lines(x),collapse="\n")
 expect_match(text,"TP: 40",fixed=TRUE)
 expect_match(text,"PREVALENCE SCENARIO (expected counts)",fixed=TRUE)
 expect_match(text,"Expected TP: 16",fixed=TRUE)
 expect_match(text,"Expected FN: 4",fixed=TRUE)
 expect_match(text,"Scenario PPV: 0.4",fixed=TRUE)
 expect_match(text,"no sampling interval",fixed=TRUE)
 expect_identical(ig_to_json(x),before)
 all_text<-paste(.ig_report_text(x),collapse="\n")
 expect_match(all_text,"RESULT SUMMARY",fixed=TRUE)
 expect_match(all_text,"TECHNICAL APPENDIX",fixed=TRUE)
 expect_match(all_text,x$snapshot_id,fixed=TRUE)
})

test_that("observed summaries identify decision, event and numerical uncertainty", {
 r<-boschloo_test(1,1,0,1,"greater",alpha=.05,groups=c("Treatment","Control"),event="Response")
 x<-create_report_snapshot(analysis_envelope(r$spec,r),locale="en")
 text<-paste(.ig_report_summary_lines(x),collapse="\n")
 expect_match(text,"Treatment: 1 / 1 events",fixed=TRUE)
 expect_match(text,"Control: 0 / 1 events",fixed=TRUE)
 expect_match(text,"Boschloo p-value: 0.25",fixed=TRUE)
 expect_match(text,"Do not reject the null hypothesis",fixed=TRUE)
 expect_match(text,"not sampling uncertainty",fixed=TRUE)
 x<-create_report_snapshot(analysis_envelope(r$spec,r),locale="de")
 text<-paste(.ig_report_summary_lines(x),collapse="\n")
 expect_match(text,"ERGEBNIS\u00dcBERSICHT",fixed=TRUE)
 expect_match(text,"Nullhypothese nicht verwerfen",fixed=TRUE)
 expect_match(text,"Stichprobenunsicherheit",fixed=TRUE)
})

test_that("PDF text width is bounded for wide and Unicode literals", {
 skip_if_not(report_pdf_capabilities()$available)
 file<-tempfile(fileext=".pdf");grDevices::cairo_pdf(file,width=8.27,height=11.69)
 on.exit({grDevices::dev.off();unlink(file)},add=TRUE)
 for(family in c("sans","mono")) {
  for(original in c(paste(rep("W",240),collapse=""),paste(rep("\u00dcberpr\u00fcfung",60),collapse=" "))) {
   lines<-.ig_report_wrap_pdf_lines(original,width_in=7.1,size=11,family=family)
   expect_gt(length(lines),1L)
   widths<-vapply(lines,function(line) grid::convertWidth(grid::grobWidth(grid::textGrob(line,gp=grid::gpar(fontsize=11,fontfamily=family))),"inches",valueOnly=TRUE),0)
   expect_true(all(widths<=7.1))
   expect_identical(gsub(" ","",paste(lines,collapse=""),fixed=TRUE),gsub(" ","",original,fixed=TRUE))
  }
 }
})

test_that("incomplete summaries and comparison scope make no invented test claim", {
 s<-analysis_spec("t_means_one_sample","a_priori",list(d=.5,alpha=.05,power=.8),list(tails="two"))
 x<-create_report_snapshot(analysis_envelope(s),report_kind="incomplete_settings")
 expect_match(paste(.ig_report_summary_lines(x),collapse="\n"),"not_computed",fixed=TRUE)
 e<-analysis_envelope(s,calculate_power(s))
 x<-create_report_snapshot(list(e,e),report_kind="comparison")
 expect_match(paste(.ig_report_summary_lines(x),collapse="\n"),"no automatic paired test",fixed=TRUE)
})


test_that("classifier summary includes actual AUC and PR summary identifiers", {
 d<-data.frame(y=c(1,1,0,0),score=c(.9,.8,.2,.1))
 s<-classification_spec(d,truth="y",positive_label=1,scores="score",cutoff=.5,score_scale="probability")
 x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
 text<-paste(.ig_report_summary_lines(x),collapse="\n")
 expect_match(text,"ROC_AUC: 1 [ok]",fixed=TRUE)
 expect_match(text,"AVERAGE_PRECISION: 1 [ok]",fixed=TRUE)
 expect_match(text,"PR_AREA_TRAPEZOID: 1 [ok]",fixed=TRUE)
 expect_match(text,"BRIER_SCORE: 0.025 [ok]",fixed=TRUE)
})
