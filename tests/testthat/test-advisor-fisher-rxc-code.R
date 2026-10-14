.fisher_rxc_env <- function() {
 e<-new.env(parent=baseenv());e$alternative<-"two.sided";e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates","tpl.fisher_exact_rxc.v1.R",package="inferguidR"),e)));e
}

test_that("general-table Fisher examples match exact rational enumeration", {
 e<-.fisher_rxc_env()
 cases<-jsonlite::fromJSON(test_path("fixtures","fisher-rxc-cases.json"),simplifyVector=FALSE)
 refs<-read.csv(test_path("fixtures","fisher-rxc-independent.csv"))
 for(i in seq_along(cases)) {
  a<-do.call(rbind,lapply(cases[[i]],unlist));r<-e$analyse(a)
  expect_identical(r$status,"ok")
  expect_lte(abs(r$p_value-refs$p_value[i]),1e-12)
  expect_null(r$confidence_interval)
 }
})

test_that("Fisher examples reject malformed counts and directional substitutions", {
 e<-.fisher_rxc_env();valid<-rbind(c(3,1,4),c(1,5,2))
 bad<-list(as.vector(valid),matrix(0,2,3),matrix(1,1,3),matrix(1,6,2),valid+.5,
   valid*100,rbind(c(-1,1,2),c(2,3,4)),rbind(c(0,0,0),c(1,1,1)),cbind(valid,0))
 for(value in c(NA,Inf)) {a<-valid;a[1]<-value;bad[[length(bad)+1L]]<-a}
 for(a in bad) expect_error(e$analyse(a))
 for(level in c(0,1,NA,Inf)) expect_error(e$analyse(valid,conf.level=level))
 expect_error(e$analyse(valid,alternative="greater"))
 expect_error(e$analyse(valid,alternative="less"))
})

test_that("Fisher RxC exports through sparse-table advisor routes in both languages", {
 a<-c(q.aim="compare",q.stage="analysis",q.outcome_type="nominal",q.unit="individuals",
   q.structure="independent",q.independence_check="no",q.n_factors="one",q.n_groups="more_than_two",
   q.covariates="no",q.small_expected="yes",q.missingness="none")
 for(lang in c("en","de")) {
  s<-advisor_replay(a,language=lang);r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,"fisher_exact_rxc")
  code<-advisor_r_code(s);expect_match(code,r$recommendation_hash,fixed=TRUE)
  e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
  expect_identical(e$example_result$status,"ok")
  expect_null(e$example_result$confidence_interval)
  expect_error(advisor_r_code(s,alternative="greater"),class="inferguidR_invalid_input")
 }
})
