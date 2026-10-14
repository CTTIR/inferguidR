test_that("exact advisor examples preserve method, direction and interval estimand", {
 ids <- c("sign_test_one_sample", "sign_test_paired", "mcnemar_exact_conditional",
          "mcnemar_midp", "fisher_exact_2x2")
 for(id in ids) for(alt in c("two.sided", "greater", "less")) {
  e <- new.env(); e$alternative <- alt; e$confidence_level <- .95
  path <- system.file("advisor", "templates", paste0("tpl.",id,".v1.R"), package="inferguidR")
  expect_true(nzchar(path)); invisible(capture.output(sys.source(path,e)))
  r <- e$example_result
  if(grepl("sign_test",id)) {
   k <- r$n_positive; n <- r$n_nonzero
  } else if(grepl("mcnemar",id)) {
   k <- r$n_x_only; n <- r$n_discordant
  } else {
   tab <- e$example_table; k <- sum(tab[1,]); m <- sum(tab[,1]); n <- sum(tab[,2]); x <- tab[1,1]
   j <- max(0,k-n):min(k,m); w <- choose(m,j)*choose(n,k-j); w <- w/sum(w)
   expected <- switch(alt,less=sum(w[j<=x]),greater=sum(w[j>=x]),two.sided=sum(w[w<=w[j==x]*(1+1e-7)]))
  }
  if(id!="fisher_exact_2x2") {
   j <- 0:n; w <- choose(n,j)/2^n; weight <- if(id=="mcnemar_midp") .5 else 1
   lo <- sum(w[j<k])+weight*w[k+1]; hi <- sum(w[j>k])+weight*w[k+1]
   expected <- switch(alt,less=lo,greater=hi,two.sided=min(1,2*min(lo,hi)))
   a <- if(alt=="two.sided") .025 else .05
   ci <- c(if(k==0 || alt=="less") 0 else qbeta(a,k,n-k+1),
           if(k==n || alt=="greater") 1 else qbeta(1-a,k+1,n-k))
   expect_equal(r$confidence_interval,ci,tolerance=1e-12)
  }
  expect_equal(r$p_value,expected,tolerance=1e-12)
  expect_identical(r$alternative,alt)
  expect_error(e$check_settings(alt,1))
  if(id=="sign_test_one_sample") {
   expect_identical(e$analyse(c(0,0),0)$status,"not_applicable")
   expect_error(e$analyse(c(1,NA),0))
  } else if(id=="sign_test_paired") {
   expect_error(e$analyse(data.frame(id=c(1,1),x=1:2,y=0)))
   expect_identical(e$analyse(data.frame(id=1:2,x=0,y=0))$status,"not_applicable")
  } else if(grepl("mcnemar",id)) {
   expect_identical(e$analyse(c(n00=2,n01=0,n10=0,n11=3))$status,"not_applicable")
   expect_error(e$analyse(c(n00=2,n01=.5,n10=2,n11=3)))
   expect_error(e$analyse(c(n00=2,n01=-1,n10=2,n11=3)))
   if(id=="mcnemar_midp") expect_match(r$interval_method,"not an inversion")
  } else {
   expect_error(e$analyse(matrix(0,2,2)))
   expect_error(e$analyse(matrix(c(1,2,3,.5),2,2)))
  }
 }
})

test_that("exact advisor examples export through actual routes and allowed overrides", {
 one <- c(q.aim="compare",q.stage="analysis",q.outcome_type="ordinal",q.unit="individuals",
          q.structure="one_sample",q.reference_value="yes",q.direction="two_sided",q.missingness="none")
 paired <- c(q.aim="compare",q.stage="analysis",q.outcome_type="ordinal",q.unit="individuals",
             q.structure="paired",q.pairing_id="yes",q.direction="two_sided",q.missingness="none")
 binary <- paired; binary["q.outcome_type"] <- "binary"
 fisher <- c(q.aim="compare",q.stage="analysis",q.outcome_type="binary",q.unit="individuals",
             q.structure="independent",q.independence_check="no",q.n_factors="one",q.n_groups="two",
             q.covariates="no",q.binary_margins="both_margins_fixed",q.direction="two_sided",q.missingness="none")
 cases <- list(sign_test_one_sample=one,sign_test_paired=paired,mcnemar_exact_conditional=binary,
               mcnemar_midp=binary,fisher_exact_2x2=fisher)
 for(id in names(cases)) for(direction in c("two_sided","one_sided_prespecified")) {
  answers <- cases[[id]]; answers["q.direction"] <- direction
  s <- advisor_replay(answers)
  if(id=="mcnemar_midp") s <- advisor_override(s,id)
  r <- advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,id)
  expect_identical(r$r_implementation$template_status,"available")
  alts <- if(direction=="two_sided") "two.sided" else c("greater","less")
  if(direction!="two_sided") expect_error(advisor_r_code(s),class="inferguidR_invalid_input")
  for(alt in alts) {
   code <- advisor_r_code(s,alternative=alt)
   expect_match(code,r$recommendation_hash,fixed=TRUE)
   e <- new.env();invisible(capture.output(eval(parse(text=code),e)))
   expect_identical(e$example_result$alternative,alt)
   expect_true(is.finite(e$example_result$p_value))
  }
 }
})
