test_that("additional advisor examples match independent expected values and reject malformed inputs", {
files <- vapply(c("t_mean_one_sample","binomial_exact_one_sample","prop_one_sample_score", "correlation_pearson","correlation_spearman","correlation_kendall"), function(id) system.file("advisor","templates",paste0("tpl.",id,".v1.R"),package="inferguidR"), "")
expect_true(all(nzchar(files)))
assert <- function(x) expect_true(x)
reject <- function(expr) assert(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
for(file in files) for(alt in c('two.sided','greater','less')) {
 e <- new.env();e$alternative<-alt;e$confidence_level<-.95
 plotfile <- tempfile(fileext='.pdf')
 grDevices::pdf(plotfile)
 tryCatch(capture.output(sys.source(file,e)),finally={grDevices::dev.off();unlink(plotfile)})
 r<-e$example_result
 if(grepl('t_mean_one_sample',file)) {
  x<-e$example_x;delta<-mean(x)-5;se<-sd(x)/sqrt(length(x));df<-length(x)-1
  z<-delta/se
  p<-switch(alt,two.sided=2*pt(-abs(z),df),greater=pt(z,df,lower.tail=FALSE),less=pt(z,df))
  reject(e$analyse(c(1,NA),5));reject(e$analyse(c(1,2),Inf))
 } else if(grepl('binomial_exact',file)) {
  mass<-dbinom(0:20,20,.2)
  p<-switch(alt,two.sided=sum(mass[mass<=dbinom(7,20,.2)*(1+1e-7)]),greater=pbinom(6,20,.2,lower.tail=FALSE),less=pbinom(7,20,.2))
  reject(e$analyse(1.5,20,.2));reject(e$analyse(21,20,.2))
 } else if(grepl('prop_one_sample',file)) {
  z<-(35/100-.2)/sqrt(.2*.8/100)
  p<-switch(alt,two.sided=2*pnorm(-abs(z)),greater=pnorm(z,lower.tail=FALSE),less=pnorm(z))
  reject(e$analyse(7,0,.2));reject(e$analyse(7,20,1))
 } else {
  x<-e$example_data$x;y<-e$example_data$y;n<-length(x)
  if(grepl('kendall',file)) {
   pairs<-combn(n,2);tau<-mean(sign(x[pairs[1,]]-x[pairs[2,]])*sign(y[pairs[1,]]-y[pairs[2,]]))
   z<-tau/sqrt(2*(2*n+5)/(9*n*(n-1)))
   p<-switch(alt,two.sided=2*pnorm(-abs(z)),greater=pnorm(z,lower.tail=FALSE),less=pnorm(z))
  } else {
   if(grepl('spearman',file)){x<-rank(x);y<-rank(y)}
   rho<-sum((x-mean(x))*(y-mean(y)))/sqrt(sum((x-mean(x))^2)*sum((y-mean(y))^2))
   z<-rho*sqrt((n-2)/(1-rho^2))
   p<-switch(alt,two.sided=2*pt(-abs(z),n-2),greater=pt(z,n-2,lower.tail=FALSE),less=pt(z,n-2))
  }
  d<-e$example_data;d$id[2]<-d$id[1];reject(e$analyse(d))
  d<-e$example_data;d$x[1]<-NA_real_;reject(e$analyse(d))
 }
 ci <- NULL
 if(grepl('t_mean_one_sample',file)) ci <- switch(alt,two.sided=delta+c(-1,1)*qt(.975,df)*se,greater=c(delta-qt(.95,df)*se,Inf),less=c(-Inf,delta+qt(.95,df)*se))
 if(grepl('binomial_exact',file)) ci <- switch(alt,two.sided=c(qbeta(.025,7,14),qbeta(.975,8,13)),greater=c(qbeta(.05,7,14),1),less=c(0,qbeta(.95,8,13)))
 if(grepl('prop_one_sample',file)) {
  zc<-qnorm(if(alt=='two.sided') .975 else .95);ph<-35/100
  centre<-(ph+zc^2/200)/(1+zc^2/100);half<-zc*sqrt(ph*(1-ph)/100+zc^2/(4*100^2))/(1+zc^2/100)
  ci<-switch(alt,two.sided=centre+c(-half,half),greater=c(centre-half,1),less=c(0,centre+half))
 }
 if(grepl('correlation_pearson',file)) {
  zz<-atanh(rho);sez<-1/sqrt(n-3)
  ci<-switch(alt,two.sided=tanh(zz+c(-1,1)*qnorm(.975)*sez),greater=c(tanh(zz-qnorm(.95)*sez),1),less=c(-1,tanh(zz+qnorm(.95)*sez)))
 }
 if(!is.null(ci)) assert(isTRUE(all.equal(as.numeric(r$confidence_interval),ci,tolerance=1e-12)))
 assert(isTRUE(all.equal(r$p_value,p,tolerance=1e-12)))
 reject(e$check_settings(alt,1))
}
})

test_that("additional templates export through their actual advisor routes", {
 mean_answers <- c(q.aim='compare',q.stage='analysis',q.outcome_type='continuous',q.unit='individuals',q.structure='one_sample',q.reference_value='yes',q.direction='two_sided',q.missingness='none')
 binary_answers <- mean_answers;binary_answers['q.outcome_type'] <- 'binary'
 cor_answers <- c(q.aim='associate',q.stage='analysis',q.unit='individuals',q.records_per_individual='one',q.assoc_types='both_continuous',q.assoc_form='linear',q.direction='two_sided',q.missingness='none')
 rank_answers <- cor_answers;rank_answers['q.assoc_form'] <- 'monotonic'
 sessions <- list(advisor_replay(mean_answers),advisor_replay(binary_answers),
   advisor_override(advisor_replay(binary_answers),'prop_one_sample_score'),
   advisor_replay(cor_answers),advisor_replay(rank_answers),
   advisor_override(advisor_replay(rank_answers),'correlation_kendall'))
 for(s in sessions) {
  r <- advisor_recommendation(s)
  expect_identical(r$r_implementation$template_status,'available')
  code <- advisor_r_code(s)
  expect_match(code,r$preferred$analysis_variant_id,fixed=TRUE)
  expect_match(code,r$recommendation_hash,fixed=TRUE)
  e <- new.env(parent=globalenv());plotfile<-tempfile(fileext='.pdf');grDevices::pdf(plotfile)
  tryCatch(capture.output(eval(parse(text=code),e)),finally={grDevices::dev.off();unlink(plotfile)})
  expect_true(is.finite(e$example_result$p_value))
 }
})
