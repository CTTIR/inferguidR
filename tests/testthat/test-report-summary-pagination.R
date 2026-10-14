test_that("text pagination preserves all content without nearly empty trailing pages", {
 for(capacity in c(50L,54L)) for(n in c(0L,1L,49L,50L,51L,54L,55L,99L,100L,101L,107L,108L,109L)) {
  original<-paste("line",seq_len(n));if(n==0L) original<-character()
  chunks<-.ig_report_page_chunks(original,capacity)
  expect_identical(unlist(chunks,use.names=FALSE),if(n==0L) NULL else original)
  expect_true(all(lengths(chunks)<=capacity))
  if(length(chunks)>1L) expect_gte(length(tail(chunks,1L)[[1]]),floor(capacity/2))
 }
})
test_that("held-out summaries disclose checked and unavailable identities in ENDE", {
 tuning<-classification_examples("scored_synthetic",part="tuning",models="model_a")$spec
 target<-classification_examples("scored_synthetic",part="test",models="model_a")$spec
 selected<-select_cutoff(tuning,"max_youden")
 imported<-ig_from_json(ig_to_json(analysis_envelope(classification_evaluated_spec(selected),selected)))$result
 for(source in list(selected,imported)) for(lang in c("en","de")) {
  r<-evaluate_cutoff(target,cutoff=source)
  report<-create_report_snapshot(analysis_envelope(classification_evaluated_spec(r),r),locale=lang)
  before<-ig_to_json(report);lines<-paste(.ig_report_summary_lines(report),collapse="\n")
  expect_match(lines,if(lang=="en") "Declared held-out evaluation" else "Deklarierte Validierung",fixed=TRUE)
  expect_match(lines,if(lang=="en") "not established" else "belegt weder",fixed=TRUE)
  checked<-identical(r$evaluation$transfer$overlap_check$unit_identity,"no_overlap_detected")
  expect_match(lines,if(lang=="en") {if(checked) "no overlap detected" else "not checked"} else {if(checked) "keine Ueberschneidung" else "nicht geprueft"},fixed=TRUE)
  expect_identical(ig_to_json(report),before)
 }
})
test_that("declaration and descriptive labels survive summary rendering without invention", {
 s<-classification_examples("scored_synthetic",models="model_a")$spec
 for(lang in c("en","de")) {
  fixed<-update(s,rule=cutoff_rule(cutoff=.5,selection_data="predeclared_rule"))
  r<-evaluate_cutoff(fixed)
  report<-create_report_snapshot(analysis_envelope(fixed,r),locale=lang)
  expect_match(paste(.ig_report_summary_lines(report),collapse="\n"),if(lang=="en") "not independently verified" else "nicht unabhaengig geprueft",fixed=TRUE)
  r<-select_cutoff(s,"max_youden")
  report<-create_report_snapshot(analysis_envelope(classification_evaluated_spec(r),r),locale=lang)
  expect_match(paste(.ig_report_summary_lines(report),collapse="\n"),if(lang=="en") "Descriptive evaluation" else "Deskriptive Auswertung",fixed=TRUE)
  report<-create_report_snapshot(analysis_envelope(s),report_kind="incomplete_settings",locale=lang)
  text<-paste(.ig_report_summary_lines(report),collapse="\n")
  expect_match(text,"not_computed",fixed=TRUE)
  expect_false(grepl("Declared held-out|Deklarierte Validierung",text))
 }
})
