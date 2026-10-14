test_that("continuous error bands separate two tails and retain exact critical edges", {
  x <- seq(-3,3,by=.5)
  data <- data.frame(hypothesis=rep(c("H0","H1"),each=length(x)),x=rep(x,2),
    y=c(dnorm(x),dnorm(x,1)),reject=rep(abs(x)>=2,2))
  b <- .ig_power_error_bands(data,c(-2,2),"density")
  a <- b[b$error=="alpha",]; beta <- b[b$error=="beta",]
  expect_setequal(unique(a$hypothesis),"H0")
  expect_setequal(unique(beta$hypothesis),"H1")
  expect_length(unique(a$group),2L)
  expect_length(unique(beta$group),1L)
  expect_true(all(abs(a$x)>=2))
  expect_equal(range(beta$x),c(-2,2))
  expect_identical(beta$y,dnorm(beta$x,1))
  for(side in c("lower","upper")) {
    d<-data;d$reject<-if(side=="lower")d$x<=-2 else d$x>=2
    z<-.ig_power_error_bands(d,if(side=="lower")-2 else 2,"density")
    expect_length(unique(z$group),2L)
    expect_true(all(z$x[z$error=="alpha"] %in% d$x[d$hypothesis=="H0" & d$reject]))
    expect_true(if(side=="lower")min(z$x[z$error=="beta"])==-2 else max(z$x[z$error=="beta"])==2)
  }
})

test_that("discrete error marks respect exact decisions without continuous boundary extension", {
  x<-0:8
  d<-data.frame(hypothesis=rep(c("H0","H1"),each=9),x=rep(x,2),
    y=c(dbinom(x,8,.3),dbinom(x,8,.7)),reject=rep(x<=0|x>=6,2))
  z<-.ig_power_error_bands(d,c(0,6),"probability_mass")
  expect_identical(z$x[z$error=="alpha"],c(0L,6L,7L,8L))
  expect_identical(z$x[z$error=="beta"],1:5)
  expect_equal(sum(z$y[z$error=="alpha"]),dbinom(0,8,.3)+pbinom(5,8,.3,lower.tail=FALSE))
  expect_equal(sum(z$y[z$error=="beta"]),pbinom(5,8,.7)-dbinom(0,8,.7))
  d$reject<-FALSE
  z<-.ig_power_error_bands(d,c(-Inf,Inf),"probability_mass")
  expect_true(all(z$error=="beta"));expect_equal(sum(z$y),1)
  d$reject<-TRUE
  z<-.ig_power_error_bands(d,c(-Inf,Inf),"probability_mass")
  expect_true(all(z$error=="alpha"));expect_equal(sum(z$y),1)
})

test_that("rendering preserves engine probabilities and immutable source rows", {
  skip_if_not_installed("ggplot2")
  cases<-list(
    analysis_spec("t_means_one_sample","post_hoc",list(d=.5,alpha=.05,n_total=30),list(tails="two")),
    analysis_spec("t_means_one_sample","post_hoc",list(d=-.5,alpha=.05,n_total=30),list(tails="one")),
    analysis_spec("exact_proportion_one_sample_binomial","post_hoc",list(g=.2,p0=.5,n_total=20,alpha=.05),list(tails="two")))
  for(s in cases) {
    r<-calculate_power(s);before<-ig_to_json(r);pd<-power_plot_data(r)
    for(lang in c("en","de")) {
      p<-plot_power(r,lang=lang)
      expect_identical(attr(p,"inferguidR_plot_data"),pd)
      expect_identical(unname(p$scales$get_scales("fill")$labels),unname(.ig_power_error_labels(r$achieved$alpha_attained,r$achieved$beta,lang)))
      expect_no_error(ggplot2::ggplot_build(p))
      standard_pdf <- tempfile(fileext=".pdf")
      expect_no_warning({
        grDevices::pdf(standard_pdf)
        tryCatch(print(p), finally=grDevices::dev.off())
      })
      expect_gt(file.size(standard_pdf),1000)
      unlink(standard_pdf)
      doc<-create_report_snapshot(analysis_envelope(s,r),locale=lang)
      hash<-doc$snapshot_hash
      if(report_pdf_capabilities()$available) {
        file<-tempfile(fileext=".pdf");on.exit(unlink(file),add=TRUE)
        expect_no_error(export_report(doc,file,"pdf"))
        expect_gt(file.size(file),1000)
      }
      expect_identical(doc$snapshot_hash,hash)
    }
    expect_identical(ig_to_json(r),before)
  }
})
