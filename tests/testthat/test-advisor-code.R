code_session <- function(direction = "two_sided") {
 advisor_replay(c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="independent",
  q.independence_check="no",q.n_factors="one",q.n_groups="two",q.covariates="no",q.estimand_continuous="mean",
  q.variance="not_justified",q.direction=direction,q.missingness="none"))
}

test_that("code export retains exact method and rejects unresolved direction", {
 s <- code_session()
 before <- advisor_save(s)
 code <- advisor_r_code(s)
 expect_match(code,"var.equal = FALSE",fixed=TRUE)
 expect_match(code,advisor_recommendation(s)$recommendation_hash,fixed=TRUE)
 expect_match(code,"q.structure = independent",fixed=TRUE)
 expect_match(code,"Template SHA-256:",fixed=TRUE)
 expect_identical(advisor_save(s),before)
 expect_error(advisor_r_code(s,"less"),class="inferguidR_invalid_input")
 expect_error(advisor_r_code(code_session("one_sided_prespecified")),class="inferguidR_invalid_input")
 expect_error(advisor_r_code(code_session("one_sided_prespecified"),"two.sided"),class="inferguidR_invalid_input")
 expect_match(advisor_r_code(code_session("one_sided_prespecified"),"greater"),'alternative <- "greater"',fixed=TRUE)
 for(level in list(0,1,NA_real_,Inf,c(.8,.9),".95"))
  expect_error(advisor_r_code(s,conf.level=level),class="inferguidR_invalid_input")
 expect_error(advisor_r_code(advisor_session()),class="inferguidR_unsupported")
 expect_match(advisor_r_code(advisor_override(s,"rank_sum_mann_whitney")),"probability_superiority",fixed=TRUE)
 student <- advisor_r_code(advisor_override(s,"t_means_two_independent_student"))
 expect_match(student,"var.equal = TRUE",fixed=TRUE)
 expect_match(student,"Status: conditional",fixed=TRUE)
 expect_match(student,"override_variant",fixed=TRUE)
})

test_that("bundled t examples agree with manual calculations and reject malformed data", {
 variants <- c("t_means_two_independent_welch","t_means_two_independent_student","t_means_paired")
 for(variant in variants) for(alt in c("two.sided","greater","less")) {
  env <- new.env(parent=globalenv());env$alternative <- alt;env$confidence_level <- .95
  file <- system.file("advisor","templates",paste0("tpl.",variant,".v1.R"),package="inferguidR")
  expect_true(nzchar(file))
  plotfile <- tempfile(fileext=".pdf")
  grDevices::pdf(plotfile)
  tryCatch(suppressMessages(capture.output(sys.source(file,env))),finally=grDevices::dev.off())
  expect_gt(file.info(plotfile)$size,1000)
  unlink(plotfile)
  paired <- identical(variant,"t_means_paired")
  if(paired) {
   x <- env$example_data$x-env$example_data$y
   delta <- mean(x);se <- stats::sd(x)/sqrt(length(x));df <- length(x)-1
   d <- env$example_data;d$id[2] <- d$id[1];expect_error(env$analyse(d))
   d <- env$example_data;d$x[1] <- NA_real_;expect_error(env$analyse(d))
   d <- env$example_data;d$y[1] <- Inf;expect_error(env$analyse(d))
   d <- env$example_data;d$id[1] <- NA_character_;expect_error(env$analyse(d))
  } else {
   x <- env$example_x;y <- env$example_y;nx <- length(x);ny <- length(y);delta <- mean(x)-mean(y)
   if(identical(variant,"t_means_two_independent_student")) {
    df <- nx+ny-2;se <- sqrt(((nx-1)*stats::var(x)+(ny-1)*stats::var(y))/df*(1/nx+1/ny))
   } else {
    vx <- stats::var(x)/nx;vy <- stats::var(y)/ny
    se <- sqrt(vx+vy);df <- (vx+vy)^2/(vx^2/(nx-1)+vy^2/(ny-1))
   }
   expect_error(env$analyse(c(1,NA),c(3,4)))
   expect_error(env$analyse(c(1,Inf),c(3,4)))
   expect_error(env$analyse(1,c(3,4)))
  }
  p <- switch(alt,two.sided=2*stats::pt(-abs(delta/se),df),greater=stats::pt(delta/se,df,lower.tail=FALSE),less=stats::pt(delta/se,df))
  ci <- switch(alt,two.sided=delta+c(-1,1)*stats::qt(.975,df)*se,
   greater=c(delta-stats::qt(.95,df)*se,Inf),less=c(-Inf,delta+stats::qt(.95,df)*se))
  expect_equal(env$example_result$df,df)
  expect_equal(env$example_result$mean_difference,delta)
  expect_equal(env$example_result$statistic,delta/se)
  expect_equal(env$example_result$p_value,p)
  expect_equal(as.numeric(env$example_result$confidence_interval),ci)
 }
})
