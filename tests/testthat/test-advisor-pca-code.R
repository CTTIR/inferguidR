.pca_env <- function() {
 e<-new.env(parent=baseenv())
 invisible(capture.output(sys.source(system.file("advisor","templates","tpl.pca_components.v1.R",package="inferguidR"),e)))
 e
}
test_that("PCA agrees with independent high-precision correlation and covariance eigensystems", {
 e<-.pca_env();cases<-jsonlite::fromJSON(test_path("fixtures","pca-cases.json"),simplifyVector=FALSE)
 refs<-jsonlite::fromJSON(test_path("fixtures","pca-independent.json"),simplifyVector=FALSE)
 for(i in seq_along(cases)) {
  z<-cases[[i]];x<-do.call(rbind,lapply(z$x,unlist));colnames(x)<-paste0("x",seq_len(z$p));d<-data.frame(id=seq_len(z$n),x)
  a<-e$analyse(d,variables=colnames(x),components=z$components,standardise=z$standardise);ref<-refs[[i]]
  expect_equal(unname(a$eigenvalues),as.numeric(unlist(ref$eigenvalues)),tolerance=2e-11)
  expect_equal(unname(a$variance_proportion),as.numeric(unlist(ref$proportion)),tolerance=2e-11)
  expect_equal(unname(a$weights %*% t(a$weights)),do.call(rbind,lapply(ref$projector,function(v) as.numeric(unlist(v)))),tolerance=2e-11)
  expect_equal(a$reconstruction_residual_sum_squares,as.numeric(ref$residual_ss),tolerance=2e-11)
  expect_equal(unname(crossprod(a$weights)),diag(z$components),tolerance=2e-12)
  expect_equal(sum(a$variance_proportion),1)
  expect_equal(a$n,z$n)
  expect_identical(a$standardised,z$standardise)
  expect_null(a$p_value)
 }
})
test_that("PCA identifies rows and preserves scaling and degenerate subspaces", {
 e<-.pca_env();d<-e$example_data;r<-e$analyse(d)
 a<-d;a$x1<-a$x1*100+17;changed<-e$analyse(a)
 expect_equal(changed$eigenvalues,r$eigenvalues,tolerance=1e-12)
 expect_equal(changed$weights%*%t(changed$weights),r$weights%*%t(r$weights),tolerance=1e-12)
 expect_equal(e$analyse(d[nrow(d):1,])$eigenvalues,r$eigenvalues,tolerance=1e-12)
 a<-d;a$x2[2]<-NA;expect_error(e$analyse(a));b<-e$analyse(a,missing="complete_cases")
 expect_identical(b$n_excluded,1L);expect_identical(b$row_ids,d$id[-2]);expect_equal(b$eigenvalues,e$analyse(d[-2,])$eigenvalues)
 a<-d;a$id[2]<-a$id[1];expect_error(e$analyse(a))
 a<-d;a$x2[2]<-Inf;expect_error(e$analyse(a,missing="complete_cases"))
 a<-d;a$x2<-1;expect_error(e$analyse(a));expect_error(e$analyse(a,standardise=FALSE))
 expect_error(e$analyse(d,components=0));expect_error(e$analyse(d,components=4));expect_error(e$analyse(d,components=1.5))
 expect_error(e$analyse(d,variables=c("x1","x1")));expect_error(e$analyse(d,variables=c("id","x1")))
 expect_error(e$analyse(d,missing="impute"));expect_error(e$analyse(d,standardise=NA));expect_error(e$analyse(d[1:2,]))
 a<-d;a$x3<-a$x1+a$x2;b<-e$analyse(a,components=2L)
 expect_identical(b$numerical_rank,2L);expect_lt(b$reconstruction_residual_sum_squares,1e-24)
 # Equal eigenvalues: compare the full plane, never the arbitrary individual axes.
 a<-data.frame(id=1:4,x1=c(1,-1,0,0),x2=c(0,0,1,-1));b<-e$analyse(a,variables=c("x1","x2"),components=2L)
 expect_equal(b$eigenvalues,c(1,1));expect_equal(unname(b$weights%*%t(b$weights)),diag(2),tolerance=1e-12)
 expect_equal(b$retained_variance,1);expect_lt(b$reconstruction_residual_sum_squares,1e-24)
})
test_that("PCA exports the exact exploratory primary or override route", {
 for(lang in c("en","de")) for(purpose in c("composite","latent")) {
  a<-c(q.aim="structure",q.stage="analysis",q.structure_goal="variables",q.reduction_purpose=purpose,q.missingness="none")
  s<-advisor_replay(a,language=lang);if(purpose=="latent") s<-advisor_override(s,"pca_components")
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,"pca_components")
  code<-advisor_r_code(s);e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
  expect_identical(e$example_result$status,"ok");expect_null(e$example_result$p_value)
  expect_false(exists("alternative",e,inherits=FALSE));expect_false(exists("confidence_level",e,inherits=FALSE))
  expect_error(advisor_r_code(s,alternative="two.sided"),class="inferguidR_invalid_input")
 }
})
